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
CONVERGENCE_TOL = 1e-9    # Target error for charge balance
MAX_ITER = 250            # High iterations for extreme pH stability
MIN_VAL = 1e-18           # Numerical floor to prevent math errors
DAVIES_LIMIT = 0.5        # Reliability limit for Davies activity model

# Physical Conversion Factors
PSI_TO_PA = 6894.76
ATM_TO_PA = 101325.0
MW = {'NH3': 17.031, 'H2S': 34.08, 'CO2': 44.01, 'H2O': 18.015}

# Edwards (1978) Equilibrium Constant Coefficients
# ln(K) = A + B/T + C/T^2 + D/T^3 + E/T^4 (T in Rankine)
# These constants define the chemical speciation (dissociation)
EDWARDS_PARAMS = {
    'NH3':   [1.587,   11160.0,   0.0,         0.0,        0.0],
    'H2S':   [-293.88, 683858.0, -6.27125e8,  2.555e11,   -3.91757e13],
    'HS':    [-220.07, 258273.0, -1.8396e8,   6.809e10,   -1.0267e13],
    'CO2':   [-241.79, 536256.0, -4.8123e8,   1.94e11,    -2.96445e13],
    'HCO3':  [-294.74, 364385.0, -2.841e8,    1.23323e11, -2.0759e13],
    'Kw':    [39.5554, -177822.0, 1.843e8,    -0.8541e11, 1.4292e13]
}

# Component Identification Tags
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
    """Calculates temperature-dependent A parameter for Davies equation."""
    t = T_K - 273.15
    return 0.4913 + 6.08e-4 * t + 5.95e-6 * t**2

def get_activity_coefficients(I, A):
    """Returns activity coefficients (gamma) using Davies modification."""
    if I < 1e-11: return 1.0, 1.0
    sqI = math.sqrt(I)
    # f(I) = -A * (sqrt(I)/(1+sqrt(I)) - 0.3*I)
    f = -A * (sqI / (1.0 + sqI) - 0.3 * I)
    return 10**f, 10**(f * 4) # z=1 and z=2

def solve_charge_balance(ph, K, molals, A_param):
    """Rigorous ionic speciation and charge balance solver."""
    h_ion = max(10**(-ph), MIN_VAL)
    m_NH3_total, m_H2S_total, m_CO2_total = molals
    
    g1, g2, ionic_strength = 1.0, 1.0, 0.0
    
    # Internal loop for Activity Coefficient convergence (Speciation-Activity Coupling)
    for _ in range(8):
        k_kw = K['Kw'] / (g1**2)
        k_nh3 = K['NH3'] # NH3 equilibrium usually handled as basicity constant
        k_h2s = K['H2S'] / (g1**2)
        k_hs = K['HS'] / g2
        k_co2 = K['CO2'] / (g1**2)
        k_hco3 = K['HCO3'] / g2
        
        oh_ion = k_kw / h_ion
        
        # Ammonia system speciation
        nh3_aq = m_NH3_total / (1 + k_nh3 * h_ion)
        nh4_ion = nh3_aq * k_nh3 * h_ion
        
        # Sulfide system speciation (H2S <-> HS- <-> S--)
        den_s = 1 + (k_h2s/h_ion) + (k_h2s * k_hs / (h_ion**2))
        h2s_aq = m_H2S_total / den_s
        hs_ion = h2s_aq * (k_h2s / h_ion)
        s_ion = hs_ion * (k_hs / h_ion)
        
        # Carbonate system speciation (CO2 <-> HCO3- <-> CO3--)
        den_c = 1 + (k_co2/h_ion) + (k_co2 * k_hco3 / (h_ion**2))
        co2_aq = m_CO2_total / den_c
        hco3_ion = co2_aq * (k_co2 / h_ion)
        co3_ion = hco3_ion * (k_hco3 / h_ion)
        
        # Re-calculate Ionic Strength: I = 0.5 * sum(m_i * z_i^2)
        ionic_strength = 0.5 * (h_ion + oh_ion + nh4_ion + hs_ion + hco3_ion + 4*s_ion + 4*co3_ion)
        g1, g2 = get_activity_coefficients(ionic_strength, A_param)
        
    # Electroneutrality Check (Positive Charges - Negative Charges)
    error = (h_ion + nh4_ion) - (oh_ion + hs_ion + 2*s_ion + hco3_ion + 2*co3_ion)
    
    speciation = {
        'NH3': nh3_aq, 'NH4': nh4_ion, 'H2S': h2s_aq, 'HS': hs_ion, 'S': s_ion, 
        'CO2': co2_aq, 'HCO3': hco3_ion, 'CO3': co3_ion, 'I': ionic_strength, 
        'g1': g1, 'g2': g2
    }
    return error, speciation

def calculate_equilibrium(T_K, m_in):
    """Predicts pH and Partial Pressures with numerical stability guards."""
    if T_K < 273.15: return {'error': "Temperature below freezing"}
    
    T_R = T_K * 1.8
    K_vals = {k: math.exp(calc_ln_k(v, T_R)) for k, v in EDWARDS_PARAMS.items()}
    A_davies = calc_davies_A(T_K)
    molal_vec = (m_in['NH3'], m_in['H2S'], m_in['CO2'])
    
    # Bisection pH Solver (Iterates 0-14 pH range)
    low, high = 0.0, 14.0
    for _ in range(MAX_ITER):
        ph_guess = (low + high) / 2.0
        err, spec = solve_charge_balance(ph_guess, K_vals, molal_vec, A_davies)
        if abs(err) < CONVERGENCE_TOL: break
        if err > 0: low = ph_guess
        else: high = ph_guess
    
    nh3, h2s, co2 = spec['NH3'], spec['H2S'], spec['CO2']
    
    # --- SAFE HENRY ENGINE (Direct Physical Reference Scaling) ---
    # Referenced at 298.15K (25C). Units: atm/molal.
    # Scaled by Van't Hoff type relationship to ensure exponential stability.
    
    # 1. NH3 Partial Pressure (Reference: 0.016 atm/molal @ 25C)
    h_nh3_ref = 0.016
    p_nh3_psia = (h_nh3_ref * nh3 * math.exp(4100 * (1/298.15 - 1/T_K))) * 14.6959
    
    # 2. H2S Partial Pressure (Reference: 0.10 atm/molal @ 25C)
    h_h2s_ref = 0.10
    p_h2s_psia = (h_h2s_ref * h2s * math.exp(2100 * (1/298.15 - 1/T_K))) * 14.6959
    
    # 3. CO2 Partial Pressure (Reference: 0.034 atm/molal @ 25C)
    p_co2_psia = 0.0
    if co2 > 0.001: # Physics Guard: Trace CO2 below 1mM contributes zero vapor
        h_co2_ref = 0.034
        p_co2_psia = (h_co2_ref * co2 * math.exp(2400 * (1/298.15 - 1/T_K))) * 14.6959
    
    # 4. Water Vapor Pressure (Kell Correlation - PSIA)
    p_water_psia = (10**(5.20389 - 1733.926/(T_K - 39.485))) * 14.5038
    
    bubble_pressure_pa = (p_nh3_psia + p_h2s_psia + p_co2_psia + p_water_psia) * PSI_TO_PA
    
    return {
        'ph': ph_guess, 'P_bubble': bubble_pressure_pa, 
        'pp': {'NH3': p_nh3_psia, 'H2S': p_h2s_psia, 'CO2': p_co2_psia, 'H2O': p_water_psia}, 
        'liq': spec
    }

# --- 3. DENSITY & PHYSICAL PROPERTIES ---

def calculate_density(T_K, P_Pa, spec, kg_water, mass_total_kg):
    """Predicts liquid density using ion-specific partial molar volumes."""
    t_c = T_K - 273.15
    # Standard pure water density curve
    rho_w = 1000 * (1 - (t_c + 288.9414)/(508929.2 * (t_c + 68.12963)) * (t_c - 3.9863)**2)
    
    # Partial Molar Volumes at infinite dilution (cm3/mol)
    v_inf = {'NH4': 18.0, 'HS': 20.0, 'NH3': 24.5, 'H2S': 35.0}
    
    vol_L = kg_water / (rho_w / 1000.0)
    vol_L += (spec['NH4']*kg_water*v_inf['NH4'] + spec['HS']*kg_water*v_inf['HS'] + 
              spec['NH3']*kg_water*v_inf['NH3'] + spec['H2S']*kg_water*v_inf['H2S']) / 1000.0
    
    rho_calculated = (mass_total_kg / vol_L) * 1000.0
    return rho_calculated * (1 + 4.5e-10 * (P_Pa - 101325.0))

# --- 4. EXECUTIVE REPORTING (ELITE LEVEL - FULL 4 SECTIONS) ---

def generate_report(res, T, P_op, flows_h, total_h, rho):
    """Generates the full professional industrial executive report."""
    l = res['liq']; pp = res['pp']
    W = 80; HR = "=" * W; SR = "-" * W
    is_f = res['P_bubble'] > (P_op * 1.01)
    status_txt = "!!! UNSTABLE - FLASHING DETECTED !!!" if is_f else "STABLE LIQUID"
    
    # Precise mg/L (ppmw) based on predicted mixture density
    h2s_mgL = (flows_h['H2S'] / total_h) * rho
    nh3_mgL = (flows_h['NH3'] / total_h) * rho

    lines = [HR, " SWEQ - SOUR WATER EQUILIBRIUM SOLVER v7.5.2 ".center(W), " Hardened Release - Thermodynamic Integrity ".center(W), HR]
    try: lines.append(" Date: %s " % DateTime.Now.ToString("yyyy-MM-dd HH:mm").center(W))
    except: pass
    lines.append((" User: %-16s Model: Edwards (1978) + Direct Scaling " % "Alexander").center(W))
    lines.append("\n")

    lines.append(" 1. EXECUTIVE SUMMARY & SAFETY CHECK ".ljust(W))
    lines.append(SR)
    if is_f:
        lines.append(" /" + "!"*76 + "\\")
        lines.append(" |" + status_txt.center(76) + "|")
        lines.append(" \\" + "!"*76 + "/")
    else: lines.append(" STATUS: " + status_txt)

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
        lines.append("\n > Activity Driver exceeds Davies limit (0.5m). Results may deviate.")
        
    lines.append(HR)
    lines.append(" End of Report")
    return "\n".join(lines)

# --- 5. MAIN BRIDGE ---

def Main():
    if 'ims1' not in globals(): return
    T = ims1.GetTemperature(); P_op = ims1.GetPressure(); F_mol = ims1.GetMolarFlow()
    if F_mol < MIN_VAL: return
    
    ids = ims1.ComponentIds; comp = ims1.GetOverallComposition()
    m_mol_s = {'NH3': 0.0, 'H2S': 0.0, 'CO2': 0.0, 'H2O': 0.0}
    idx = {}
    for i, n in enumerate(ids):
        nc = n.lower()
        for k, tags in COMP_MAP.items():
            if any(t in nc for t in tags): m_mol_s[k] = F_mol * comp[i]; idx[k] = i; break
            
    if m_mol_s['H2O'] < MIN_VAL: return
    kg_w = m_mol_s['H2O'] * MW['H2O'] / 1000.0
    molals = {k: v / kg_w for k, v in m_mol_s.items() if k != 'H2O'}
    
    res = calculate_equilibrium(T, molals)
    if res.get('error'): return

    P_out = min(res['P_bubble'], 1000 * ATM_TO_PA)
    y_raw = [0.0]*len(ids)
    for k, i in idx.items():
        if k in res['pp']: y_raw[i] = (res['pp'][k] * PSI_TO_PA) / P_out
    
    total_y = sum(y_raw)
    y_norm = [v/total_y if total_y > 0 else comp[i] for v in y_raw]
    
    oms1.SetTemperature(T); oms1.SetPressure(P_out); oms1.SetOverallComposition(Array[Double](y_norm)); oms1.SetMolarFlow(1e-8)
    oms2.SetTemperature(T); oms2.SetPressure(P_out); oms2.SetOverallComposition(ims1.GetOverallComposition()); oms2.SetMolarFlow(F_mol - 1e-8)
    oms1.Calculate(); oms2.Calculate()

    m_kg_s = {k: v * MW.get(k, 28) / 1000.0 for k, v in m_mol_s.items()}
    total_m_kg_s = sum(m_kg_s.values())
    rho = calculate_density(T, P_out, res['liq'], kg_w, total_m_kg_s)

    flows_h = {k: m * 3600 for k, m in m_kg_s.items()}
    total_h = sum(flows_h.values())
    
    try:
        path = os.path.join(os.path.expanduser('~'), 'Desktop', 'SWEQ_Datasheet.txt')
        with open(path, 'w') as f: f.write(generate_report(res, T, P_op, flows_h, total_h, rho))
        if os.name == 'nt': os.startfile(path)
    except: pass

Main()
