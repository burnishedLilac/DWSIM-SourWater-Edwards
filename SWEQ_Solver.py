import math
import sys
import os

# Try importing .NET libraries for DWSIM interaction
# Fails gracefully if running in a standard Python shell for testing
try:
    import clr
    from System import Array, Double
    from System import DateTime
except ImportError:
    pass

# --- Configuration & Constants ---
CO2_CRIT_P_PSI = 1070.6
CONVERGENCE_TOL = 1e-9
MAX_ITER = 200
TRACES = 1e-16
DAVIES_LIMIT = 0.5  # Molality limit for valid activity correction

# Unit Conversions
PSI_TO_PA = 6894.76
ATM_TO_PA = 101325.0

MW = {
    'NH3': 17.031, 
    'H2S': 34.08, 
    'CO2': 44.01, 
    'H2O': 18.015
}

# Edwards et al. (1978) Correlation Coefficients (A-E)
EDWARDS_PARAMS = {
    'NH3':   [1.587,   11160.0,   0.0,         0.0,        0.0],
    'H2S':   [-293.88, 683858.0, -6.27125e8,  2.555e11,   -3.91757e13],
    'HS':    [-220.07, 258273.0, -1.8396e8,   6.809e10,   -1.0267e13],
    'CO2':   [-241.79, 536256.0, -4.8123e8,   1.94e11,    -2.96445e13],
    'HCO3':  [-294.74, 364385.0, -2.841e8,    1.23323e11, -2.0759e13],
    'Kw':    [39.5554, -177822.0, 1.843e8,    -0.8541e11, 1.4292e13]
}

# Maps DWSIM component IDs (flexible) to internal keys
COMP_MAP = {
    'NH3': ['ammonia', 'nh3', 'amonia'],
    'H2S': ['hydrogen sulfide', 'h2s', 'sulfeto'],
    'CO2': ['carbon dioxide', 'co2', 'dioxido'],
    'H2O': ['water', 'h2o', 'agua']
}

# --- Core Logic ---

def calc_ln_k(coeffs, T_Rankine):
    """Computes equilibrium constant ln(K) based on Edwards correlation."""
    A, B, C, D, E = coeffs
    return A + B/T_Rankine + C/(T_Rankine**2) + D/(T_Rankine**3) + E/(T_Rankine**4)

def calc_water_dielectric_A(T_Kelvin):
    """Returns the Davies 'A' parameter adjusted for temperature."""
    t = T_Kelvin - 273.15
    return 0.4913 + 6.08e-4 * t + 5.95e-6 * t**2

def get_gammas(ionic_strength, A_param):
    """Calculates activity coefficients using Davies Equation."""
    if ionic_strength < 1e-9: return 1.0, 1.0
    sq_I = math.sqrt(ionic_strength)
    exponent = -A_param * (sq_I / (1.0 + sq_I) - 0.3 * ionic_strength)
    return 10**(exponent), 10**(exponent * 4)

def solve_charge_balance(ph_guess, K, molals, A_param):
    """Iteratively solves speciation for a given pH."""
    h_ion = max(10**(-ph_guess), TRACES)
    m_NH3, m_H2S, m_CO2 = molals
    
    g1, g2 = 1.0, 1.0
    I = 0.0
    
    # Convergence loop for Activity Coefficients
    for _ in range(3):
        k_kw   = K['Kw'] / (g1**2)
        k_nh3  = K['NH3']
        k_h2s  = K['H2S'] / (g1**2)
        k_hs   = K['HS'] / g2
        k_co2  = K['CO2'] / (g1**2)
        k_hco3 = K['HCO3'] / g2
        
        oh_ion = k_kw / h_ion
        
        nh3_aq = m_NH3 / (1 + k_nh3 * h_ion)
        nh4_ion = nh3_aq * k_nh3 * h_ion
        
        den_s = 1 + (k_h2s/h_ion) + (k_h2s * k_hs / (h_ion**2))
        h2s_aq = m_H2S / den_s
        hs_ion = h2s_aq * (k_h2s / h_ion)
        s_ion  = hs_ion * (k_hs / h_ion)
        
        den_c = 1 + (k_co2/h_ion) + (k_co2 * k_hco3 / (h_ion**2))
        co2_aq = m_CO2 / den_c
        hco3_ion = co2_aq * (k_co2 / h_ion)
        co3_ion  = hco3_ion * (k_hco3 / h_ion)
        
        I = 0.5 * (h_ion + nh4_ion + oh_ion + hs_ion + 4*s_ion + hco3_ion + 4*co3_ion)
        g1, g2 = get_gammas(I, A_param)
        
    balance_error = (h_ion + nh4_ion) - (oh_ion + hs_ion + 2*s_ion + hco3_ion + 2*co3_ion)
    
    results = {
        'NH3': nh3_aq, 'NH4': nh4_ion,
        'H2S': h2s_aq, 'HS': hs_ion, 'S': s_ion,
        'CO2': co2_aq, 'HCO3': hco3_ion, 'CO3': co3_ion,
        'I': I, 'g1': g1, 'g2': g2
    }
    return balance_error, results

def calculate_equilibrium(T_K, molality_in):
    """Main solver driver."""
    if T_K < 273.15 or T_K > 523.15:
        return {'error': "Temperature %.2f K out of bounds" % T_K}
    
    T_R = T_K * 1.8
    K_vals = {k: math.exp(calc_ln_k(v, T_R)) for k, v in EDWARDS_PARAMS.items()}
    A_davies = calc_water_dielectric_A(T_K)
    
    m_vals = (molality_in['NH3'], molality_in['H2S'], molality_in['CO2'])
    
    ph_low, ph_high = 0.0, 14.0
    err_low, _ = solve_charge_balance(ph_low, K_vals, m_vals, A_davies)
    err_high, _ = solve_charge_balance(ph_high, K_vals, m_vals, A_davies)

    if err_low * err_high > 0:
        ph = 0.0 if abs(err_low) < abs(err_high) else 14.0
        _, spec = solve_charge_balance(ph, K_vals, m_vals, A_davies)
    else:
        ph = 7.0
        spec = {}
        for _ in range(MAX_ITER):
            err, spec = solve_charge_balance(ph, K_vals, m_vals, A_davies)
            if abs(err) < CONVERGENCE_TOL: break
            if err > 0: ph_low = ph
            else: ph_high = ph
            ph = (ph_low + ph_high) / 2.0

    nh3, h2s, co2 = spec['NH3'], spec['H2S'], spec['CO2']
    
    # Henry's Law
    ln_p_nh3 = 178.339 - 15517.91/T_R - 25.6767*math.log(T_R) + 0.01966*T_R + (131.4/T_R - 0.1682)*nh3 + 0.06*(2*molality_in['CO2'] + molality_in['H2S'])
    p_nh3 = math.exp(ln_p_nh3) * nh3
    
    ln_p_h2s = 8.8 + 0.015 * (T_R - 560)
    p_h2s = math.exp(ln_p_h2s + (-0.05*nh3 + (0.965 - 486.0/T_R)*molality_in['CO2'])) * h2s
    
    clamped = False
    p_co2 = 0.0
    if molality_in['CO2'] > TRACES:
        ln_h_co2 = 301.68 - 34096.6/T_R + 1.2285e8/(T_R**2) - 6.4752e10/(T_R**3) + 1.1557e13/(T_R**4)
        if (molality_in['NH3'] + molality_in['CO2'] + molality_in['H2S']) > TRACES:
            ln_h_co2 += -0.09 * (nh3 - molality_in['CO2'] - molality_in['H2S'])
        limit = math.log(CO2_CRIT_P_PSI)
        try:
            val = ln_h_co2 + math.log(max(co2, 1e-20))
            if val > limit:
                p_co2 = CO2_CRIT_P_PSI; clamped = True
            else:
                p_co2 = math.exp(ln_h_co2) * co2
        except: p_co2 = CO2_CRIT_P_PSI; clamped = True

    p_w = (10**(5.20389 - 1733.926/(T_K - 39.485))) * 14.5038
    p_total_pa = (p_nh3 + p_h2s + p_co2 + p_w) * PSI_TO_PA

    return {
        'ph': ph, 'P_calc': p_total_pa, 'clamped': clamped, 'A_davies': A_davies,
        'pp': {'NH3': p_nh3, 'H2S': p_h2s, 'CO2': p_co2, 'H2O': p_w},
        'liq': spec
    }

def generate_report(res, inputs):
    """Generates the Executive Datasheet (Restored v7.1 Layout)."""
    
    T, P, mass_flows, total_mass, rho = inputs
    l = res['liq']
    pp = res['pp']
    
    is_flashing = res['P_calc'] > (P * 1.01)
    status_txt = "STABLE LIQUID"
    if is_flashing: status_txt = "!!! UNSTABLE - FLASHING DETECTED !!!"
    
    rho_kg_L = rho / 1000.0
    h2s_mgL = (mass_flows['H2S'] / total_mass * 1e6) * rho_kg_L
    nh3_mgL = (mass_flows['NH3'] / total_mass * 1e6) * rho_kg_L

    W = 80
    HR = "=" * W
    SR = "-" * W
    
    lines = []
    
    # HEADER
    lines.append(HR)
    lines.append(" SWEQ - SOUR WATER EQUILIBRIUM SOLVER v7.2".center(W))
    lines.append(" Thermodynamic Datasheet".center(W))
    lines.append(HR)
    try:
        from System import DateTime
        dt = DateTime.Now.ToString("yyyy-MM-dd HH:mm")
    except:
        dt = "Current"
        
    lines.append((" Date: %-16s Model: Edwards (1978)" % dt).center(W))
    lines.append((" User: %-16s Activity: Davies (T-Dep)" % os.environ.get('USERNAME', 'Engineer')).center(W))
    lines.append("\n")

    # EXECUTIVE SUMMARY
    lines.append(" 1. EXECUTIVE SUMMARY & SAFETY CHECK".ljust(W))
    lines.append(SR)
    
    if is_flashing:
        lines.append(" /" + "!"*76 + "\\")
        lines.append(" |" + status_txt.center(76) + "|")
        lines.append(" \\" + "!"*76 + "/")
    else:
        lines.append(" STATUS: " + status_txt)
        
    lines.append("")
    lines.append(" %-25s | %-15s | %-10s | %-20s" % ("PARAMETER", "VALUE", "UNIT", "COMMENT"))
    lines.append(" %s-+-%s-+-%s-+-%s" % ("-"*25, "-"*15, "-"*10, "-"*20))
    lines.append(" %-25s | %15.2f | %-10s | Inlet Condition" % ("Temperature", T-273.15, "C"))
    lines.append(" %-25s | %15.3f | %-10s | System Pressure" % ("Operating Pressure", P/ATM_TO_PA, "atm"))
    lines.append(" %-25s | %15.3f | %-10s | Equilibrium Press" % ("Bubble Pressure", res['P_calc']/ATM_TO_PA, "atm"))
    lines.append(" %-25s | %15.4f | %-10s | Charge Balance" % ("Calculated pH", res['ph'], "-"))
    lines.append(" %-25s | %15.4f | %-10s | Activity Driver" % ("Ionic Strength (I)", l['I'], "molal"))
    lines.append("\n")

    # STREAM COMPOSITION (Corrected Water-First Logic)
    lines.append(" 2. STREAM COMPOSITION & FLOWS".ljust(W))
    lines.append(SR)
    lines.append(" %-25s | %-15s | %-10s " % ("COMPONENT", "MASS FLOW (kg/h)", "MASS %"))
    lines.append(" %s-+-%s-+-%s" % ("-"*25, "-"*15, "-"*10))
    
    # 1. Print Water First
    w_mass = mass_flows.get('H2O', 0.0)
    w_pct = w_mass / total_mass * 100 if total_mass > 0 else 0
    lines.append(" %-25s | %15.2f | %10.2f " % ("Water (H2O)", w_mass*3600, w_pct))
    
    # 2. Print Solutes
    for k, v in mass_flows.items():
        if k == 'H2O': continue # Skip water as it is already printed
        m_pct = v / total_mass * 100 if total_mass > 0 else 0
        lines.append(" %-25s | %15.2f | %10.2f " % (k, v*3600, m_pct))
        
    lines.append(" %-25s | %15.2f | %10.2f " % ("TOTAL", total_mass*3600, 100.00))
    lines.append("\n")

    # SPECIATION
    lines.append(" 3. CHEMICAL SPECIATION (LIQUID PHASE)".ljust(W))
    lines.append(SR)
    lines.append(" Concentrations in molality (mol/kg H2O)")
    lines.append("")
    lines.append("  AMMONIA SYSTEM   ::  NH3(aq): %.4f  <==>  NH4+: %.4f" % (l['NH3'], l['NH4']))
    lines.append("  SULFIDE SYSTEM   ::  H2S(aq): %.4f  <==>  HS-:  %.4f  <==>  S--: %.2e" % (l['H2S'], l['HS'], l['S']))
    lines.append("  CARBONATE SYSTEM ::  CO2(aq): %.4f  <==>  HCO3: %.4f  <==>  CO3: %.4f" % (l['CO2'], l['HCO3'], l['CO3']))
    lines.append("")
    lines.append("  > Activity Coeffs (Davies):  Gamma(1) = %.4f  |  Gamma(2) = %.4f" % (l['g1'], l['g2']))
    if l['I'] > DAVIES_LIMIT:
        lines.append("  > NOTE: Ionic strength %.2fm exceeds Davies limit (%.1fm)." % (l['I'], DAVIES_LIMIT))
    lines.append("\n")

    # ENVIRONMENTAL
    lines.append(" 4. ENVIRONMENTAL & VAPOR PREDICTIONS".ljust(W))
    lines.append(SR)
    lines.append(" %-38s | %s" % ("ENVIRONMENTAL (Liquid Effluent)", "VAPOR PHASE (Equilibrium Partial P)"))
    lines.append(" %s-+-%s" % ("-"*38, "-"*39))
    lines.append(" H2S Total: %15.1f mg/L      | p(H2S): %10.4f psia" % (h2s_mgL, pp['H2S']))
    lines.append(" NH3 Total: %15.1f mg/L      | p(NH3): %10.4f psia" % (nh3_mgL, pp['NH3']))
    lines.append(" Density:   %15.2f kg/m3     | p(CO2): %10.4f psia" % (rho, pp['CO2']))
    lines.append("                                        | p(H2O): %10.4f psia" % (pp['H2O']))
    
    if res['clamped']:
        lines.append(" %-38s | [!] CO2 Saturation Clamp Active" % "")
        
    lines.append(HR)
    lines.append(" End of Report")
    
    return "\n".join(lines)

def Main():
    """Bridge function between DWSIM and Python Logic."""
    try:
        if 'ims1' not in globals(): return 
        
        T = ims1.GetTemperature()
        P = ims1.GetPressure()
        F_mol = ims1.GetMolarFlow()
        if F_mol < TRACES: return
        
        comp = ims1.GetOverallComposition()
        ids = ims1.ComponentIds
        
        rho = 1000.0
        if hasattr(ims1, 'GetMassDensity'): rho = ims1.GetMassDensity()

        mol_in = {'NH3': 0.0, 'H2S': 0.0, 'CO2': 0.0, 'H2O': 0.0}
        idx_map = {}
        active_indices = []
        
        for i, name in enumerate(ids):
            n_clean = name.lower()
            mol_flow = F_mol * comp[i]
            for key, tags in COMP_MAP.items():
                if any(t in n_clean for t in tags):
                    mol_in[key] = mol_flow
                    idx_map[key] = i
                    active_indices.append(i)
                    break
        
        if mol_in['H2O'] < TRACES: return 

        kg_water = mol_in['H2O'] * MW['H2O'] / 1000.0
        molality_in = {k: v/kg_water for k, v in mol_in.items() if k != 'H2O'}
        
        res = calculate_equilibrium(T, molality_in)
        
        if res.get('error'):
            ims1.FlowSheet.WriteMessage("SWEQ Error: " + res['error'])
            return

        P_out = min(res['P_calc'], 1500 * ATM_TO_PA)
        pp = res['pp']
        y_out = [0.0] * len(ids)
        
        for k, idx in idx_map.items():
            if k in pp: y_out[idx] = (pp[k] * PSI_TO_PA) / P_out
        
        for i in range(len(ids)):
            if i not in active_indices: y_out[i] = comp[i]
        
        total_y = sum(y_out)
        if total_y > 0: y_out = [val/total_y for val in y_out]
        else: y_out = comp

        oms1.SetTemperature(T)
        oms1.SetPressure(P_out)
        oms1.SetOverallComposition(Array[Double](y_out))
        oms1.SetMolarFlow(1e-8)

        oms2.SetTemperature(T)
        oms2.SetPressure(P_out)
        oms2.SetOverallComposition(Array[Double](list(comp)))
        oms2.SetMolarFlow(F_mol - 1e-8)
        
        oms1.Calculate()
        oms2.Calculate()

        try:
            # FIX: Divide by 1000.0 to get kg/s
            mass_flows = {k: v * MW.get(k, 28) / 1000.0 for k, v in mol_in.items()}
            
            # FIX: Total mass was adding water twice. mass_flows ALREADY includes H2O.
            total_mass = sum(mass_flows.values()) 
            
            report_text = generate_report(res, (T, P, mass_flows, total_mass, rho))
            
            path = os.path.join(os.environ['USERPROFILE'], 'Desktop', 'SWEQ_Datasheet.txt')
            with open(path, 'w') as f: f.write(report_text)
            
            os.startfile(path)
            ims1.FlowSheet.WriteMessage("SWEQ Report generated at: " + path)
            
        except Exception as e:
            ims1.FlowSheet.WriteMessage("Report Generation Failed: " + str(e))

    except Exception as e:
        if 'ims1' in globals():
            ims1.FlowSheet.WriteMessage("SWEQ Runtime Error: " + str(e))

Main()
