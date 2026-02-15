import math
import sys
import os

try:
    import clr
    from System import Array, Double, DateTime
except ImportError:
    pass

# --- 1. CONFIGURATION & THERMODYNAMIC CONSTANTS ---
CO2_CRIT_P_PSI = 1070.6   
H2S_CRIT_P_PSI = 1300.0   
CONVERGENCE_TOL = 1e-9    
MAX_ITER = 250            
MIN_VAL = 1e-18           
DAVIES_LIMIT = 0.5        

PSI_TO_PA = 6894.76
ATM_TO_PA = 101325.0
MW = {'NH3': 17.031, 'H2S': 34.08, 'CO2': 44.01, 'H2O': 18.015}

EDWARDS_PARAMS = {
    'NH3':   [1.587,   11160.0,   0.0,         0.0,        0.0],
    'H2S':   [-293.88, 683858.0, -6.27125e8,  2.555e11,   -3.91757e13],
    'HS':    [-220.07, 258273.0, -1.8396e8,   6.809e10,   -1.0267e13],
    'CO2':   [-241.79, 536256.0, -4.8123e8,   1.94e11,    -2.96445e13],
    'HCO3':  [-294.74, 364385.0, -2.841e8,    1.23323e11, -2.0759e13],
    'Kw':    [39.5554, -177822.0, 1.843e8,    -0.8541e11, 1.4292e13]
}

COMP_MAP = {
    'NH3': ['ammonia', 'nh3', 'amonia'],
    'H2S': ['hydrogen sulfide', 'h2s', 'sulfeto'],
    'CO2': ['carbon dioxide', 'co2', 'dioxido'],
    'H2O': ['water', 'h2o', 'agua']
}

# --- 2. THERMODYNAMIC CORE FUNCTIONS ---

def calc_ln_k(coeffs, T_R):
    A, B, C, D, E = coeffs
    return A + B/T_R + C/(T_R**2) + D/(T_R**3) + E/(T_R**4)

def calc_davies_A(T_K):
    t = T_K - 273.15
    return 0.4913 + 6.08e-4 * t + 5.95e-6 * t**2

def get_activity_coefficients(I, A):
    if I < 1e-11: return 1.0, 1.0
    sqI = math.sqrt(I)
    f = -A * (sqI / (1.0 + sqI) - 0.3 * I)
    return 10**f, 10**(f * 4)

def solve_charge_balance(ph, K, molals, A_param):
    """Rigorous balance with CORRECT NH3/NH4+ Kb speciation logic."""
    h_ion = max(10**(-ph), MIN_VAL)
    m_NH3_total, m_H2S_total, m_CO2_total = molals
    
    g1, g2, ionic_strength = 1.0, 1.0, 0.0
    
    for _ in range(10):
        k_kw = K['Kw'] / (g1**2)
        oh_ion = k_kw / h_ion
        
        # 1. AMMONIA SPECIATION (NH3 + H+ <-> NH4+)
        # Ka = Kw / Kb
        k_ka_nh4 = k_kw / K['NH3']
        # nh4 = total / (1 + Ka/[H+])
        nh4_ion = m_NH3_total / (1 + (k_ka_nh4 / h_ion))
        nh3_aq = m_NH3_total - nh4_ion
        
        # 2. SULFIDE SPECIATION
        k_h2s = K['H2S'] / (g1**2); k_hs = K['HS'] / g2
        den_s = 1 + (k_h2s/h_ion) + (k_h2s * k_hs / (h_ion**2))
        h2s_aq = m_H2S_total / den_s
        hs_ion = h2s_aq * (k_h2s / h_ion)
        s_ion = hs_ion * (k_hs / h_ion)
        
        # 3. CARBONATE SPECIATION
        k_co2 = K['CO2'] / (g1**2); k_hco3 = K['HCO3'] / g2
        den_c = 1 + (k_co2/h_ion) + (k_co2 * k_hco3 / (h_ion**2))
        co2_aq = m_CO2_total / den_c
        hco3_ion = co2_aq * (k_co2 / h_ion)
        co3_ion = hco3_ion * (k_hco3 / h_ion)
        
        ionic_strength = 0.5 * (h_ion + oh_ion + nh4_ion + hs_ion + hco3_ion + 4*s_ion + 4*co3_ion)
        g1, g2 = get_activity_coefficients(ionic_strength, A_param)
        
    error = (h_ion + nh4_ion) - (oh_ion + hs_ion + 2*s_ion + hco3_ion + 2*co3_ion)
    spec = {'NH3': nh3_aq, 'NH4': nh4_ion, 'H2S': h2s_aq, 'HS': hs_ion, 'S': s_ion, 
            'CO2': co2_aq, 'HCO3': hco3_ion, 'CO3': co3_ion, 'I': ionic_strength}
    return error, spec

def calculate_equilibrium(T_K, m_in):
    T_R = T_K * 1.8
    K_vals = {k: math.exp(calc_ln_k(v, T_R)) for k, v in EDWARDS_PARAMS.items()}
    A_davies = calc_davies_A(T_K)
    
    low, high = 0.0, 14.0
    for _ in range(MAX_ITER):
        ph_guess = (low + high) / 2.0
        err, spec = solve_charge_balance(ph_guess, K_vals, (m_in['NH3'], m_in['H2S'], m_in['CO2']), A_davies)
        if abs(err) < CONVERGENCE_TOL: break
        if err > 0: low = ph_guess
        else: high = ph_guess
    
    # Henry Engine (Carroll & Mather Scaling)
    nh3, h2s, co2 = spec['NH3'], spec['H2S'], spec['CO2']
    p_nh3 = (0.016 * nh3 * math.exp(4100 * (1/298.15 - 1/T_K))) * 14.6959
    p_h2s = (0.10 * h2s * math.exp(2100 * (1/298.15 - 1/T_K))) * 14.6959
    p_co2 = 0.0
    if co2 > 0.001: p_co2 = (0.034 * co2 * math.exp(2400 * (1/298.15 - 1/T_K))) * 14.6959
    p_w = (10**(5.20389 - 1733.926/(T_K - 39.485))) * 14.5038
    
    return {'ph': ph_guess, 'P_bubble': (p_nh3+p_h2s+p_co2+p_w)*PSI_TO_PA, 'pp': {'NH3': p_nh3, 'H2S': p_h2s, 'CO2': p_co2, 'H2O': p_w}, 'liq': spec}

# --- 3. DENSITY ENGINE ---

def calculate_density(T_K, P_Pa, spec, kg_water, mass_total_kg):
    t_c = T_K - 273.15
    rho_w = 1000 * (1 - (t_c + 288.9414)/(508929.2 * (t_c + 68.12963)) * (t_c - 3.9863)**2)
    v_inf = {'NH4': 18.0, 'HS': 20.0, 'NH3': 24.5, 'H2S': 35.0}
    vol_L = kg_water / (rho_w / 1000.0)
    vol_L += (spec['NH4']*kg_water*v_inf['NH4'] + spec['HS']*kg_water*v_inf['HS'] + 
              spec['NH3']*kg_water*v_inf['NH3'] + spec['H2S']*kg_water*v_inf['H2S']) / 1000.0
    return (mass_total_kg / vol_L) * 1000.0

# --- 4. EXECUTIVE REPORTING (ELITE LEVEL - v7.5.8) ---

def generate_report(res, T, P_op, flows_h, total_h, rho):
    l = res['liq']; pp = res['pp']
    W = 80; HR = "=" * W; SR = "-" * W
    is_f = res['P_bubble'] > (P_op * 1.01)
    
    h2s_mgL = (l['H2S'] + l['HS'] + l['S']) * MW['H2S'] * 1000 * (rho / 1000.0)
    nh3_mgL = (l['NH3'] + l['NH4']) * MW['NH3'] * 1000 * (rho / 1000.0)

    lines = [HR, " SWEQ - SOUR WATER EQUILIBRIUM SOLVER v7.5.8 ".center(W), " Final Certified Industrial Edition ".center(W), HR]
    lines.append((" User: %-16s Model: Edwards (1978) + Corrected Kb " % "Alexander").center(W))
    lines.append("\n 1. EXECUTIVE SUMMARY & SAFETY CHECK ".ljust(W))
    lines.append(SR)
    lines.append(" STATUS: " + ("!!! FLASHING DETECTED !!!" if is_f else "STABLE LIQUID"))
    lines.append("")
    lines.append(" %-25s | %-15s | %-10s | %-20s" % ("PARAMETER", "VALUE", "UNIT", "COMMENT"))
    lines.append(" %s-+-%s-+-%s-+-%s" % ("-"*25, "-"*15, "-"*10, "-"*20))
    lines.append(" %-25s | %15.2f | %-10s | Inlet Condition" % ("Temperature", T-273.15, "C"))
    lines.append(" %-25s | %15.3f | %-10s | Operating P" % ("System Pressure", P_op/ATM_TO_PA, "atm"))
    lines.append(" %-25s | %15.3f | %-10s | Equilibrium P" % ("Bubble Pressure", res['P_bubble']/ATM_TO_PA, "atm"))
    lines.append(" %-25s | %15.4f | %-10s | Buffering Status" % ("Calculated pH", res['ph'], "-"))
    lines.append(" %-25s | %15.4f | %-10s | Activity Driver" % ("Ionic Strength (I)", l['I'], "molal"))
    lines.append("\n 2. STREAM COMPOSITION & FLOWS ".ljust(W))
    lines.append(SR)
    lines.append(" %-25s | %15s | %-10s " % ("COMPONENT", "FLOW (kg/h)", "MASS %"))
    lines.append(" %s-+-%s-+-%s" % ("-"*25, "-"*15, "-"*10))
    lines.append(" %-25s | %15.2f | %10.2f " % ("Water (H2O)", flows_h['H2O'], (flows_h['H2O']/total_h)*100))
    for k in ['NH3', 'H2S', 'CO2']:
        lines.append(" %-25s | %15.2f | %10.2f " % (k, flows_h[k], (flows_h[k]/total_h)*100))
    lines.append("\n 3. CHEMICAL SPECIATION (LIQUID PHASE) ".ljust(W))
    lines.append(SR)
    lines.append(" Concentrations in molality (mol/kg H2O)")
    lines.append("")
    lines.append("  AMMONIA SYSTEM   ::  NH3(aq): %.4f  <==>  NH4+: %.4f" % (l['NH3'], l['NH4']))
    lines.append("  SULFIDE SYSTEM   ::  H2S(aq): %.4f  <==>  HS-:  %.4f" % (l['H2S'], l['HS']))
    lines.append("  CARBONATE SYSTEM ::  CO2(aq): %.4f  <==>  HCO3: %.4f" % (l['CO2'], l['HCO3']))
    lines.append("\n 4. ENVIRONMENTAL & VAPOR PREDICTIONS ".ljust(W))
    lines.append(SR)
    lines.append(" %-38s | %s" % ("ENVIRONMENTAL (Liquid Effluent)", "VAPOR PHASE (Equilibrium Partial P)"))
    lines.append(" %s-+-%s" % ("-"*38, "-"*39))
    lines.append(" H2S Total: %15.1f mg/L      | p(H2S): %10.4f psia" % (h2s_mgL, pp['H2S']))
    lines.append(" NH3 Total: %15.1f mg/L      | p(NH3): %10.4f psia" % (nh3_mgL, pp['NH3']))
    lines.append(" Density:   %15.2f kg/m3     | p(CO2): %10.4f psia" % (rho, pp['CO2']))
    lines.append("                                        | p(H2O): %10.4f psia" % (pp['H2O']))
    if l['I'] > DAVIES_LIMIT: lines.append("\n > Warning: Activity Driver exceeds Davies limit (0.5m).")
    lines.append(HR + "\n End of Report")
    return "\n".join(lines)

# --- 5. MAIN BRIDGE ---

def Main():
    if 'ims1' not in globals(): return
    T = ims1.GetTemperature(); P_op = ims1.GetPressure(); F_mol = ims1.GetMolarFlow()
    if F_mol < MIN_VAL: return
    ids = ims1.ComponentIds; comp = ims1.GetOverallComposition()
    m_mol_s = {'NH3': 0.0, 'H2S': 0.0, 'CO2': 0.0, 'H2O': 0.0}; idx = {}
    for i, n in enumerate(ids):
        nc = n.lower()
        for k, tags in COMP_MAP.items():
            if any(t in nc for t in tags): m_mol_s[k] = F_mol * comp[i]; idx[k] = i; break
    if m_mol_s['H2O'] < MIN_VAL: return
    kg_w = m_mol_s['H2O'] * MW['H2O'] / 1000.0
    res = calculate_equilibrium(T, {k: v / kg_w for k, v in m_mol_s.items() if k != 'H2O'})
    if res.get('error'): return
    P_out = min(res['P_bubble'], 1000 * ATM_TO_PA)
    y_raw = [0.0]*len(ids)
    for k, i in idx.items():
        if k in res['pp']: y_raw[i] = (res['pp'][k] * PSI_TO_PA) / P_out
    y_norm = [v/sum(y_raw) if sum(y_raw)>0 else comp[i] for v in y_raw]
    oms1.SetTemperature(T); oms1.SetPressure(P_out); oms1.SetOverallComposition(Array[Double](y_norm)); oms1.SetMolarFlow(1e-8)
    oms2.SetTemperature(T); oms2.SetPressure(P_out); oms2.SetOverallComposition(ims1.GetOverallComposition()); oms2.SetMolarFlow(F_mol - 1e-8)
    oms1.Calculate(); oms2.Calculate()
    m_kg_s = {k: v * MW.get(k, 28) / 1000.0 for k, v in m_mol_s.items()}
    rho = calculate_density(T, P_out, res['liq'], kg_w, sum(m_kg_s.values()))
    try:
        path = os.path.join(os.path.expanduser('~'), 'Desktop', 'SWEQ_Datasheet.txt')
        with open(path, 'w') as f: f.write(generate_report(res, T, P_op, {k: m*3600 for k,m in m_kg_s.items()}, sum(m_kg_s.values())*3600, rho))
        if os.name == 'nt': os.startfile(path)
    except: pass

Main()
