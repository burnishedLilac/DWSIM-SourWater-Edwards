import math
import sys
import os

try:
    import clr
    from System import Array, Double, DateTime
except:
    pass

# --- Constants & Config ---
CO2_CRIT_P_PSI = 1070.6
CONVERGENCE_TOL = 1e-9
MAX_ITER = 200
MIN_VAL = 1e-16
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

# --- Core Thermodynamics ---

def calc_ln_k(coeffs, T_R):
    A, B, C, D, E = coeffs
    return A + B/T_R + C/(T_R**2) + D/(T_R**3) + E/(T_R**4)

def calc_davies_A(T_K):
    t = T_K - 273.15
    return 0.4913 + 6.08e-4 * t + 5.95e-6 * t**2

def get_gammas(I, A):
    if I < 1e-9: return 1.0, 1.0
    sqI = math.sqrt(I)
    f = -A * (sqI / (1.0 + sqI) - 0.3 * I)
    return 10**f, 10**(f * 4)

def solve_charge_balance(ph, K, molals, A_param):
    h = max(10**(-ph), MIN_VAL)
    m_NH3, m_H2S, m_CO2 = molals
    g1, g2, I = 1.0, 1.0, 0.0
    for _ in range(5):
        k_kw = K['Kw']/(g1**2); k_nh3 = K['NH3']
        k_h2s = K['H2S']/(g1**2); k_hs = K['HS']/g2
        k_co2 = K['CO2']/(g1**2); k_hco3 = K['HCO3']/g2
        oh = k_kw / h
        nh3 = m_NH3 / (1 + k_nh3 * h); nh4 = nh3 * k_nh3 * h
        den_s = 1 + (k_h2s/h) + (k_h2s * k_hs / (h**2))
        h2s = m_H2S / den_s; hs = h2s * (k_h2s / h); s = hs * (k_hs / h)
        den_c = 1 + (k_co2/h) + (k_co2 * k_hco3 / (h**2))
        co2 = m_CO2 / den_c; hco3 = co2 * (k_co2 / h); co3 = hco3 * (k_hco3 / h)
        I = 0.5 * (h + oh + nh4 + hs + hco3 + 4*s + 4*co3)
        g1, g2 = get_gammas(I, A_param)
    bal = (h + nh4) - (oh + hs + 2*s + hco3 + 2*co3)
    spec = {'NH3': nh3, 'NH4': nh4, 'H2S': h2s, 'HS': hs, 'S': s, 
            'CO2': co2, 'HCO3': hco3, 'CO3': co3, 'I': I, 'g1': g1, 'g2': g2}
    return bal, spec

def calculate_equilibrium(T_K, m_in):
    if T_K < 273.15 or any(v < 0 for v in m_in.values()):
        return {'error': "Invalid Physical Inputs"}
    T_R = T_K * 1.8
    K_v = {k: math.exp(calc_ln_k(v, T_R)) for k, v in EDWARDS_PARAMS.items()}
    A_d = calc_davies_A(T_K)
    m_v = (m_in['NH3'], m_in['H2S'], m_in['CO2'])
    low, high = 0.0, 14.0
    err_l, _ = solve_charge_balance(low, K_v, m_v, A_d)
    err_h, _ = solve_charge_balance(high, K_v, m_v, A_d)
    if err_l * err_h > 0:
        ph = 0.0 if abs(err_l) < abs(err_h) else 14.0
        _, s = solve_charge_balance(ph, K_v, m_v, A_d)
    else:
        ph = 7.0
        for _ in range(MAX_ITER):
            err, s = solve_charge_balance(ph, K_v, m_v, A_d)
            if abs(err) < CONVERGENCE_TOL: break
            if err > 0: low = ph
            else: high = ph
            ph = (low + high) / 2.0
    nh3, h2s, co2 = s['NH3'], s['H2S'], s['CO2']
    ln_nh3 = 178.339 - 15517.91/T_R - 25.6767*math.log(T_R) + 0.01966*T_R + (131.4/T_R - 0.1682)*nh3 + 0.06*(2*m_in['CO2'] + m_in['H2S'])
    p_nh3 = math.exp(ln_nh3) * nh3
    ln_h2s = 8.8 + 0.015 * (T_R - 560)
    p_h2s = math.exp(ln_h2s + (-0.05*nh3 + (0.965 - 486.0/T_R)*m_in['CO2'])) * h2s
    p_co2 = 0.0; clamped = False
    if m_in['CO2'] > MIN_VAL:
        ln_h_co2 = 301.68 - 34096.6/T_R + 1.2285e8/(T_R**2) - 6.4752e10/(T_R**3) + 1.1557e13/(T_R**4)
        if (m_in['NH3'] + m_in['CO2'] + m_in['H2S']) > MIN_VAL: ln_h_co2 += -0.09*(nh3 - m_in['CO2'] - m_in['H2S'])
        val = ln_h_co2 + math.log(max(co2, 1e-20))
        if val > math.log(CO2_CRIT_P_PSI): p_co2 = CO2_CRIT_P_PSI; clamped = True
        else: p_co2 = math.exp(val)
    p_w = (10**(5.20389 - 1733.926/(T_K - 39.485))) * 14.5038
    p_tot_pa = (p_nh3 + p_h2s + p_co2 + p_w) * PSI_TO_PA
    return {'ph': ph, 'P_calc': p_tot_pa, 'clamped': clamped, 'pp': {'NH3': p_nh3, 'H2S': p_h2s, 'CO2': p_co2, 'H2O': p_w}, 'liq': s}

# --- Predictive Density Engine ---

def calculate_density_rigorous(T_K, P_Pa, spec, kg_water, mass_total_kg):
    """
    Predicts density using ion-specific partial molar volumes.
    """
    t_c = T_K - 273.15
    # Pure water density (Kell correlation)
    rho_w = 1000 * (1 - (t_c + 288.9414)/(508929.2 * (t_c + 68.12963)) * (t_c - 3.9863)**2)
    
    # Infinite dilution partial molar volumes (cm3/mol)
    # Ref: Adjusted from Laliberte & Cooper ionic parameters
    v_inf = {
        'NH4': 18.0, 
        'HS':  20.0, 
        'OH':  -4.0, 
        'H':   0.0,
        'NH3': 24.5,
        'H2S': 35.0
    }
    
    # Volume of water in Liters
    V_L = kg_water / (rho_w / 1000.0)
    
    # Contributions from each species (n * V_inf / 1000)
    # Molality * kg_water = mols
    V_L += (spec['NH4'] * kg_water * v_inf['NH4']) / 1000.0
    V_L += (spec['HS']  * kg_water * v_inf['HS'])  / 1000.0
    V_L += (spec['NH3'] * kg_water * v_inf['NH3']) / 1000.0
    V_L += (spec['H2S'] * kg_water * v_inf['H2S']) / 1000.0
    
    # Final density calculation (kg / L) * 1000 = kg/m3
    rho_calc = (mass_total_kg / V_L) * 1000.0
    
    # Compressibility correction
    return rho_calc * (1 + 4.5e-10 * (P_Pa - 101325.0))

# --- Reporting ---

def generate_report(res, T, P_op, flows_h, total_h, rho):
    l = res['liq']; pp = res['pp']
    W = 80; HR = "=" * W; SR = "-" * W
    is_f = res['P_calc'] > (P_op * 1.01)
    status_txt = "STABLE LIQUID"
    if is_f: status_txt = "!!! UNSTABLE - FLASHING DETECTED !!!"
    
    # Rigorous environmental concentration using calculated density
    h2s_mgL = (flows_h['H2S'] / total_h) * rho
    nh3_mgL = (flows_h['NH3'] / total_h) * rho

    lines = [HR, " SWEQ - SOUR WATER EQUILIBRIUM SOLVER v7.3.6 ".center(W), " Ion-Specific Density Release ".center(W), HR]
    try: lines.append(" Date: %s " % DateTime.Now.ToString("yyyy-MM-dd HH:mm").center(W))
    except: pass
    lines.append((" User: %-16s Model: Edwards (1978) " % os.environ.get('USERNAME', 'Engineer')).center(W))
    lines.append("\n")

    lines.append(" 1. EXECUTIVE SUMMARY & SAFETY CHECK ".ljust(W))
    lines.append(SR)
    if is_f:
        lines.append(" /" + "!"*76 + "\\")
        lines.append(" |" + status_txt.center(76) + "|")
        lines.append(" \\" + "!"*76 + "/")
    else:
        lines.append(" STATUS: " + status_txt)

    lines.append("")
    lines.append(" %-25s | %-15s | %-10s | %-20s" % ("PARAMETER", "VALUE", "UNIT", "COMMENT"))
    lines.append(" %s-+-%s-+-%s-+-%s" % ("-"*25, "-"*15, "-"*10, "-"*20))
    lines.append(" %-25s | %15.2f | %-10s | Inlet Condition" % ("Temperature", T-273.15, "C"))
    lines.append(" %-25s | %15.3f | %-10s | Operating P" % ("System Pressure", P_op/ATM_TO_PA, "atm"))
    lines.append(" %-25s | %15.3f | %-10s | Equilibrium P" % ("Bubble Pressure", res['P_calc']/ATM_TO_PA, "atm"))
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
        lines.append("\n > Activity Driver exceeds Davies limit (0.5m).")
    lines.append(HR)
    lines.append(" End of Report")
    return "\n".join(lines)

# --- Entry Point ---

def Main():
    if 'ims1' not in globals(): return
    T = ims1.GetTemperature(); P = ims1.GetPressure(); F = ims1.GetMolarFlow()
    if F < MIN_VAL: return
    
    ids = ims1.ComponentIds; comp = ims1.GetOverallComposition()
    m_in = {'NH3': 0.0, 'H2S': 0.0, 'CO2': 0.0, 'H2O': 0.0}
    idx_m = {}
    for i, n in enumerate(ids):
        nc = n.lower()
        for k, tags in COMP_MAP.items():
            if any(t in nc for t in tags): m_in[k] = F * comp[i]; idx_m[k] = i; break
            
    if m_in['H2O'] < MIN_VAL: return
    kg_w = m_in['H2O'] * MW['H2O'] / 1000.0
    molals = {k: v/kg_w for k, v in m_in.items() if k != 'H2O'}
    
    res = calculate_equilibrium(T, molals)
    if res.get('error'): return

    P_o = min(res['P_calc'], 1500 * ATM_TO_PA)
    y_raw = [0.0]*len(ids)
    for k, i in idx_m.items():
        if k in res['pp']: y_raw[i] = (res['pp'][k] * PSI_TO_PA) / P_o
    
    total_y = sum(y_raw)
    y_norm = [v/total_y if total_y > 0 else comp[i] for i, v in enumerate(y_raw)]
    
    oms1.SetTemperature(T); oms1.SetPressure(P_o); oms1.SetOverallComposition(Array[Double](y_norm)); oms1.SetMolarFlow(1e-8)
    oms2.SetTemperature(T); oms2.SetPressure(P_o); oms2.SetOverallComposition(ims1.GetOverallComposition()); oms2.SetMolarFlow(F - 1e-8)
    oms1.Calculate(); oms2.Calculate()

    # RIGOROUS DENSITY ENGINE
    m_flows_kg_s = {k: v * MW.get(k, 28) / 1000.0 for k, v in m_in.items()}
    total_mass_kg_s = sum(m_flows_kg_s.values())
    rho = calculate_density_rigorous(T, P_o, res['liq'], kg_w, total_mass_kg_s)

    f_h = {k: m * 3600 for k, m in m_flows_kg_s.items()}
    t_h = sum(f_h.values())
    
    try:
        desktop = os.path.join(os.path.expanduser('~'), 'Desktop')
        path = os.path.join(desktop, 'SWEQ_Datasheet.txt')
        with open(path, 'w') as f: f.write(generate_report(res, T, P, f_h, t_h, rho))
        if os.name == 'nt': os.startfile(path)
    except: pass

Main()
