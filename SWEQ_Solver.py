import math
import sys
import os

try:
    import clr
    from System import Array, Double, DateTime
except ImportError:
    pass

# --- 1. CONFIGURATION & THERMODYNAMIC CONSTANTS ---

# Safety limits and convergence parameters
CO2_CRIT_P_PSI = 1070.6   # CO2 critical pressure for safety clamping
CONVERGENCE_TOL = 1e-9    # Target error for charge balance
MAX_ITER = 200            # Max iterations for pH bisection
MIN_VAL = 1e-16           # Numerical floor to prevent log(0)
DAVIES_LIMIT = 0.5        # Reliability limit for Davies activity model

# Physical Conversion Factors
PSI_TO_PA = 6894.76
ATM_TO_PA = 101325.0
MW = {'NH3': 17.031, 'H2S': 34.08, 'CO2': 44.01, 'H2O': 18.015}

# Edwards (1978) Equilibrium Constant Coefficients
# ln(K) = A + B/T + C/T^2 + D/T^3 + E/T^4 (T in Rankine)
EDWARDS_PARAMS = {
    'NH3':   [1.587,   11160.0,   0.0,         0.0,        0.0],
    'H2S':   [-293.88, 683858.0, -6.27125e8,  2.555e11,   -3.91757e13],
    'HS':    [-220.07, 258273.0, -1.8396e8,   6.809e10,   -1.0267e13],
    'CO2':   [-241.79, 536256.0, -4.8123e8,   1.94e11,    -2.96445e13],
    'HCO3':  [-294.74, 364385.0, -2.841e8,    1.23323e11, -2.0759e13],
    'Kw':    [39.5554, -177822.0, 1.843e8,    -0.8541e11, 1.4292e13]
}

# String matching tags for DWSIM component identification
COMP_MAP = {
    'NH3': ['ammonia', 'nh3', 'amonia'],
    'H2S': ['hydrogen sulfide', 'h2s', 'sulfeto'],
    'CO2': ['carbon dioxide', 'co2', 'dioxido'],
    'H2O': ['water', 'h2o', 'agua']
}

# --- 2. THERMODYNAMIC CORE FUNCTIONS ---

def calc_ln_k(coeffs, T_R):
    """Calculates ln(K) for a reaction at temperature T (Rankine)."""
    A, B, C, D, E = coeffs
    return A + B/T_R + C/(T_R**2) + D/(T_R**3) + E/(T_R**4)

def calc_davies_A(T_K):
    """Calculates the temperature-dependent A parameter for the Davies equation."""
    t = T_K - 273.15
    return 0.4913 + 6.08e-4 * t + 5.95e-6 * t**2

def get_activity_coefficients(I, A):
    """Returns activity coefficients (gamma) for ions with z=1 and z=2."""
    if I < 1e-9: return 1.0, 1.0
    sqI = math.sqrt(I)
    # Davies modification: f(I) = -A * (sqrt(I)/(1+sqrt(I)) - 0.3*I)
    f = -A * (sqI / (1.0 + sqI) - 0.3 * I)
    return 10**f, 10**(f * 4) # z=1 gamma, z=2 gamma (z^2 dependence)

def solve_charge_balance(ph, K, molals, A_param):
    """Calculates chemical speciation and electroneutrality error for a given pH."""
    h_ion = max(10**(-ph), MIN_VAL)
    m_NH3_total, m_H2S_total, m_CO2_total = molals
    
    # Initialize activity coefficients
    g1, g2, ionic_strength = 1.0, 1.0, 0.0
    
    # Loop to converge Activity Coefficients with Ionic Strength
    for _ in range(5):
        # Correct K values for activity
        k_kw = K['Kw'] / (g1**2)
        k_nh3 = K['NH3'] # NH3 hydrolysis k is often used without gammas in this form
        k_h2s = K['H2S'] / (g1**2)
        k_hs = K['HS'] / g2
        k_co2 = K['CO2'] / (g1**2)
        k_hco3 = K['HCO3'] / g2
        
        oh_ion = k_kw / h_ion
        
        # Ammonia system
        nh3_aq = m_NH3_total / (1 + k_nh3 * h_ion)
        nh4_ion = nh3_aq * k_nh3 * h_ion
        
        # Sulfide system (H2S -> HS- -> S--)
        den_s = 1 + (k_h2s/h_ion) + (k_h2s * k_hs / (h_ion**2))
        h2s_aq = m_H2S_total / den_s
        hs_ion = h2s_aq * (k_h2s / h_ion)
        s_ion = hs_ion * (k_hs / h_ion)
        
        # Carbonate system (CO2 -> HCO3- -> CO3--)
        den_c = 1 + (k_co2/h_ion) + (k_co2 * k_hco3 / (h_ion**2))
        co2_aq = m_CO2_total / den_c
        hco3_ion = co2_aq * (k_co2 / h_ion)
        co3_ion = hco3_ion * (k_hco3 / h_ion)
        
        # Update Ionic Strength: I = 0.5 * sum(m_i * z_i^2)
        ionic_strength = 0.5 * (h_ion + oh_ion + nh4_ion + hs_ion + hco3_ion + 4*s_ion + 4*co3_ion)
        g1, g2 = get_activity_coefficients(ionic_strength, A_param)
        
    # Sum of charges: positive - negative
    error = (h_ion + nh4_ion) - (oh_ion + hs_ion + 2*s_ion + hco3_ion + 2*co3_ion)
    
    speciation = {
        'NH3': nh3_aq, 'NH4': nh4_ion, 'H2S': h2s_aq, 'HS': hs_ion, 'S': s_ion, 
        'CO2': co2_aq, 'HCO3': hco3_ion, 'CO3': co3_ion, 'I': ionic_strength, 
        'g1': g1, 'g2': g2
    }
    return error, speciation

def calculate_equilibrium(T_K, m_in):
    """Orchestrates pH bisection and fugacity calculations."""
    if T_K < 273.15 or any(v < 0 for v in m_in.values()):
        return {'error': "Invalid Physical Inputs"}
        
    T_R = T_K * 1.8
    K_vals = {k: math.exp(calc_ln_k(v, T_R)) for k, v in EDWARDS_PARAMS.items()}
    A_davies = calc_davies_A(T_K)
    molal_vec = (m_in['NH3'], m_in['H2S'], m_in['CO2'])
    
    # pH Bisection Solver
    low, high = 0.0, 14.0
    for _ in range(MAX_ITER):
        ph_guess = (low + high) / 2.0
        err, spec = solve_charge_balance(ph_guess, K_vals, molal_vec, A_davies)
        if abs(err) < CONVERGENCE_TOL: break
        if err > 0: low = ph_guess
        else: high = ph_guess
    
    # Henry's Law Predictions
    nh3, h2s, co2 = spec['NH3'], spec['H2S'], spec['CO2']
    
    # 1. NH3 Partial Pressure (Edwards Interaction Parameters)
    ln_p_nh3 = 178.339 - 15517.91/T_R - 25.6767*math.log(T_R) + 0.01966*T_R + \
               (131.4/T_R - 0.1682)*nh3 + 0.06*(2*m_in['CO2'] + m_in['H2S'])
    p_nh3 = math.exp(ln_p_nh3) * nh3
    
    # 2. H2S Partial Pressure (Rigorous Edwards Henry Correlation)
    H_h2s_coeffs = [158.36, -9172.1, -23.01, 0.022] # Result in PSIA
    ln_H_h2s = H_h2s_coeffs[0] + H_h2s_coeffs[1]/T_K + H_h2s_coeffs[2]*math.log(T_K) + H_h2s_coeffs[3]*T_K
    p_h2s = math.exp(ln_H_h2s - 0.05*nh3) * h2s
    
    # 3. CO2 Partial Pressure (With Trace Species numerical guard)
    p_co2 = 0.0
    clamped = False
    co2_limit = m_in['CO2'] * 0.0001
    if m_in['CO2'] > MIN_VAL and co2 > max(co2_limit, 1e-10):
        ln_h_co2 = 301.68 - 34096.6/T_R + 1.2285e8/(T_R**2) - 6.4752e10/(T_R**3) + 1.1557e13/(T_R**4)
        if (m_in['NH3'] + m_in['CO2'] + m_in['H2S']) > MIN_VAL:
            ln_h_co2 += -0.09*(nh3 - m_in['CO2'] - m_in['H2S'])
        
        val = ln_h_co2 + math.log(co2)
        if val > math.log(CO2_CRIT_P_PSI): 
            p_co2 = CO2_CRIT_P_PSI
            clamped = True
        else:
            p_co2 = math.exp(val)
            
    # 4. Water Vapor Pressure (Kell Equation)
    p_water = (10**(5.20389 - 1733.926/(T_K - 39.485))) * 14.5038
    
    bubble_pressure_pa = (p_nh3 + p_h2s + p_co2 + p_water) * PSI_TO_PA
    
    return {
        'ph': ph_guess, 
        'P_bubble': bubble_pressure_pa, 
        'clamped': clamped, 
        'pp': {'NH3': p_nh3, 'H2S': p_h2s, 'CO2': p_co2, 'H2O': p_water}, 
        'liq': spec
    }

# --- 3. PHYSICAL PROPERTY ENGINES ---

def calculate_density(T_K, P_Pa, spec, kg_water, mass_total_kg):
    """Predicts liquid mixture density using Ion-Specific Partial Molar Volumes."""
    t_c = T_K - 273.15
    # Pure water density
    rho_w = 1000 * (1 - (t_c + 288.9414)/(508929.2 * (t_c + 68.12963)) * (t_c - 3.9863)**2)
    
    # Infinite dilution partial molar volumes (cm3/mol)
    v_inf = {'NH4': 18.0, 'HS': 20.0, 'OH': -4.0, 'H': 0.0, 'NH3': 24.5, 'H2S': 35.0}
    
    # Initial volume from water (Liters)
    vol_L = kg_water / (rho_w / 1000.0)
    
    # Add solute volume contributions (molality * kg_w * Vinf / 1000)
    vol_L += (spec['NH4'] * kg_water * v_inf['NH4']) / 1000.0
    vol_L += (spec['HS']  * kg_water * v_inf['HS'])  / 1000.0
    vol_L += (spec['NH3'] * kg_water * v_inf['NH3']) / 1000.0
    vol_L += (spec['H2S'] * kg_water * v_inf['H2S']) / 1000.0
    
    # Resulting density (kg/m3)
    rho_calculated = (mass_total_kg / vol_L) * 1000.0
    
    # Compressibility correction: rho = rho0 * (1 + beta*dP)
    return rho_calculated * (1 + 4.5e-10 * (P_Pa - 101325.0))

# --- 4. REPORT GENERATION ---

def generate_datasheet(res, T, P_op, flows_h, total_h, rho):
    """Builds the professional executive report string."""
    l = res['liq']; pp = res['pp']
    W = 80; HR = "=" * W; SR = "-" * W
    
    is_flashing = res['P_bubble'] > (P_op * 1.01)
    
    # Precise mg/L (ppm vol) calculation using the predicted density
    h2s_mgL = (flows_h['H2S'] / total_h) * rho
    nh3_mgL = (flows_h['NH3'] / total_h) * rho

    lines = [HR, " SWEQ - SOUR WATER EQUILIBRIUM SOLVER v7.4.0 ".center(W), " Hardened Industrial Release ".center(W), HR]
    try: lines.append(" Date: %s " % DateTime.Now.ToString("yyyy-MM-dd HH:mm").center(W))
    except: pass
    lines.append((" User: %-16s Model: Edwards (1978) " % os.environ.get('USERNAME', 'Engineer')).center(W))
    lines.append("\n")

    lines.append(" 1. EXECUTIVE SUMMARY & SAFETY CHECK ".ljust(W))
    lines.append(SR)
    if is_flashing:
        lines.append(" /" + "!"*76 + "\\")
        lines.append(" |" + "!!! UNSTABLE - FLASHING DETECTED !!!".center(76) + "|")
        lines.append(" \\" + "!"*76 + "/")
    else:
        lines.append(" STATUS: STABLE LIQUID")

    lines.append("")
    lines.append(" %-25s | %-15s | %-10s | %-20s" % ("PARAMETER", "VALUE", "UNIT", "COMMENT"))
    lines.append(" %s-+-%s-+-%s-+-%s" % ("-"*25, "-"*15, "-"*10, "-"*20))
    lines.append(" %-25s | %15.2f | %-10s | Inlet Condition" % ("Temperature", T-273.15, "C"))
    lines.append(" %-25s | %15.3f | %-10s | Operating P" % ("System Pressure", P_op/ATM_TO_PA, "atm"))
    lines.append(" %-25s | %15.3f | %-10s | Equilibrium P" % ("Bubble Pressure", res['P_bubble']/ATM_TO_PA, "atm"))
    lines.append(" %-25s | %15.4f | %-10s | Charge Balance" % ("Calculated pH", res['ph'], "-"))
    lines.append(" %-25s | %15.4f | %-10s | Activity Driver" % ("Ionic Strength (I)", l['I'], "molal"))
    lines.append("\n")

    lines.append(" 2. STREAM COMPOSITION & FLOWS ".ljust(W))
    lines.append(SR)
    lines.append(" %-25s | %15s | %-10s " % ("COMPONENT", "FLOW (kg/h)", "MASS %"))
    lines.append(" %s-+-%s-+-%s" % ("-"*25, "-"*15, "-"*10))
    lines.append(" %-25s | %15.2f | %10.2f " % ("Water (H2O)", flows_h['H2O'], (flows_h['H2O']/total_h)*100))
    for k in ['NH3', 'H2S', 'CO2']:
        lines.append(" %-25s | %15.2f | %10.2f " % (k, flows_h[k], (flows_h[k]/total_h)*100))
    lines.append(" %-25s | %15.2f | %10.2f " % ("TOTAL", total_h, 100.00))
    lines.append("\n")

    lines.append(" 3. CHEMICAL SPECIATION (LIQUID PHASE) ".ljust(W))
    lines.append(SR)
    lines.append(" Concentrations in molality (mol/kg H2O)")
    lines.append("")
    lines.append("  AMMONIA SYSTEM   ::  NH3(aq): %.4f  <==>  NH4+: %.4f" % (l['NH3'], l['NH4']))
    lines.append("  SULFIDE SYSTEM   ::  H2S(aq): %.4f  <==>  HS-:  %.4f  <==>  S--: %.2e" % (l['H2S'], l['HS'], l['S']))
    lines.append("  CARBONATE SYSTEM ::  CO2(aq): %.4f  <==>  HCO3: %.4f  <==>  CO3: %.4f" % (l['CO2'], l['HCO3'], l['CO3']))
    lines.append("\n")

    lines.append(" 4. ENVIRONMENTAL & VAPOR PREDICTIONS ".ljust(W))
    lines.append(SR)
    lines.append(" %-38s | %s" % ("ENVIRONMENTAL (Liquid Effluent)", "VAPOR PHASE (Equilibrium Partial P)"))
    lines.append(" %s-+-%s" % ("-"*38, "-"*39))
    lines.append(" H2S Total: %15.1f mg/L      | p(H2S): %10.4f psia" % (h2s_mgL, pp['H2S']))
    lines.append(" NH3 Total: %15.1f mg/L      | p(NH3): %10.4f psia" % (nh3_mgL, pp['NH3']))
    lines.append(" Density:   %15.2f kg/m3     | p(CO2): %10.4f psia" % (rho, pp['CO2']))
    lines.append("                                        | p(H2O): %10.4f psia" % (pp['H2O']))
    
    if l['I'] > DAVIES_LIMIT:
        lines.append("\n > Activity Driver exceeds Davies limit (0.5m). Use caution.")
    if res['clamped']:
        lines.append(" > Warning: CO2 Fugacity Clamped at Critical Pressure.")
    lines.append(HR)
    lines.append(" End of Report")
    return "\n".join(lines)

# --- 5. MAIN EXECUTION BRIDGE ---

def Main():
    """DWSIM Main Entry Point."""
    if 'ims1' not in globals(): return
    T = ims1.GetTemperature()
    P_op = ims1.GetPressure()
    F_molar = ims1.GetMolarFlow()
    if F_molar < MIN_VAL: return
    
    ids = ims1.ComponentIds
    comp = ims1.GetOverallComposition()
    m_mol_s = {'NH3': 0.0, 'H2S': 0.0, 'CO2': 0.0, 'H2O': 0.0}
    indices = {}
    
    # Map components from DWSIM stream
    for i, name in enumerate(ids):
        nc = name.lower()
        for k, tags in COMP_MAP.items():
            if any(t in nc for t in tags):
                m_mol_s[k] = F_molar * comp[i]
                indices[k] = i
                break
                
    if m_mol_s['H2O'] < MIN_VAL: return
    
    # Convert to Molality (mol/kg solvent)
    kg_water = m_mol_s['H2O'] * MW['H2O'] / 1000.0
    molalities = {k: v / kg_water for k, v in m_mol_s.items() if k != 'H2O'}
    
    # Solve Thermodynamics
    results = calculate_equilibrium(T, molalities)
    if results.get('error'): return

    # Manage DWSIM Output Streams (Vapor oms1, Liquid oms2)
    P_out = min(results['P_bubble'], 1500 * ATM_TO_PA)
    y_raw = [0.0] * len(ids)
    for k, i in indices.items():
        if k in results['pp']:
            y_raw[i] = (results['pp'][k] * PSI_TO_PA) / P_out
    
    # Normalize vapor phase composition
    total_y = sum(y_raw)
    y_norm = [v / total_y if total_y > 0 else comp[i] for i, v in enumerate(y_raw)]
    
    # Write to oms1 (Vapor Trace)
    oms1.SetTemperature(T); oms1.SetPressure(P_out); oms1.SetOverallComposition(Array[Double](y_norm)); oms1.SetMolarFlow(1e-8)
    # Write to oms2 (Liquid Bulk)
    oms2.SetTemperature(T); oms2.SetPressure(P_out); oms2.SetOverallComposition(ims1.GetOverallComposition()); oms2.SetMolarFlow(F_molar - 1e-8)
    oms1.Calculate(); oms2.Calculate()

    # Calculate final mass flows and predictive density
    m_kg_s = {k: v * MW.get(k, 28) / 1000.0 for k, v in m_mol_s.items()}
    total_mass_kg_s = sum(m_kg_s.values())
    rho_predicted = calculate_density(T, P_out, results['liq'], kg_water, total_mass_kg_s)

    flows_kg_h = {k: m * 3600 for k, m in m_kg_s.items()}
    total_kg_h = sum(flows_kg_h.values())
    
    # Report Output
    try:
        desktop_path = os.path.join(os.path.expanduser('~'), 'Desktop', 'SWEQ_Datasheet.txt')
        with open(desktop_path, 'w') as f:
            f.write(generate_datasheet(results, T, P_op, flows_kg_h, total_kg_h, rho_predicted))
        if os.name == 'nt': os.startfile(desktop_path)
    except Exception as e:
        ims1.FlowSheet.WriteMessage("SWEQ Error: " + str(e))

Main()
