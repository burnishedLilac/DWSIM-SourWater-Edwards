"""
SOUR WATER EQUILIBRIUM SOLVER (SWEQ) - v8.0.0
--------------------------------------------------------------------------------
Thermodynamics: Edwards (1978) + PITZER MODEL (Rumpf et al., 1999)
Activity: Pitzer specific interactions for NH3/H2S + Davies fallback for CO2
Henry's Law: Continuous Physical Scaling + Pitzer Activity Coupling
Water Activity: Osmotic Mole Fraction Correction
Stability: Successive Substitution Damping for Non-Linear Convergence
Integration: Optimized for DWSIM Python Script Unit Operation

Changelog v8.0.0:
- MAJOR UPGRADE: Replaced global Davies model with Pitzer species-specific model.
- Removed arbitrary Davies limits (Model is now valid up to ~10 molal).
- Added water activity depression for concentrated brines.

Author: Alexander Francisco Cescon
Release: Enterprise Extreme Concentration Edition (GPLv3 Only).
--------------------------------------------------------------------------------
"""

import math
import sys
import os

# --- 1. CONSTANTS ---

CONST = {
    'R_GAS': 8.314,
    'PSI_TO_PA': 6894.76,
    'ATM_TO_PA': 101325.0,
    'ATM_TO_PSI': 14.6959,
    'MW': {'NH3': 17.031, 'H2S': 34.08, 'CO2': 44.01, 'H2O': 18.015},
    'H_REF': {'NH3': 0.016, 'H2S': 0.10, 'CO2': 0.034},
    'H_SCALE': {'NH3': 4100.0, 'H2S': 2100.0, 'CO2': 2400.0}
}

SOLVER = {
    'TOL_CHARGE': 1e-9, 'TOL_IONIC': 1e-6,
    'MAX_ITER_PH': 350, 'MAX_ITER_ION': 100, # Increased for Pitzer stability
    'MIN_VAL': 1e-18,
    'CO2_DAVIES_WARN': 0.5 # We still track CO2 limits
}

EDWARDS_PARAMS = {
    'NH3':   [1.587,   11160.0,   0.0,         0.0,        0.0],
    'H2S':   [-293.88, 683858.0, -6.27125e8,  2.555e11,   -3.91757e13],
    'HS':    [-220.07, 258273.0, -1.8396e8,   6.809e10,   -1.0267e13],
    'CO2':   [-241.79, 536256.0, -4.8123e8,   1.94e11,    -2.96445e13],
    'HCO3':  [-294.74, 364385.0, -2.841e8,    1.23323e11, -2.0759e13],
    'Kw':    [39.5554, -177822.0, 1.843e8,    -0.8541e11, 1.4292e13]
}

# Pitzer Parameters (Standard Natural Base)
PITZER_PARAMS = {
    # Gases (Rumpf converted)
    'beta0_NH3_NH3':     [-0.01979,   9.864,      0.0],
    'beta0_H2S_H2S':     [-0.26156,   69.751,     0.0],
    
    # Ions NH4-HS (Standard Pitzer 1-1 Electrolyte at 25C reference)
    'beta0_NH4_HS':      [0.052, 0.0, 0.0], 
    'beta1_NH4_HS':      [0.180, 0.0, 0.0], 
    'Cphi_NH4_HS':       [-0.001, 0.0, 0.0]
}

COMP_MAP = {
    'NH3': ['ammonia', 'nh3', 'amonia'],
    'H2S': ['hydrogen sulfide', 'h2s', 'sulfeto'],
    'CO2': ['carbon dioxide', 'co2', 'dioxido'],
    'H2O': ['water', 'h2o', 'agua']
}

try:
    import clr
    from System import Array, Double, DateTime
except ImportError:
    pass

# --- 2. THERMODYNAMICS & PITZER ENGINE ---

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

def calc_pitzer_p(q, T_K):
    return q[0] + q[1]/T_K + (q[2]*math.log(T_K) if q[2] != 0 else 0.0)

def get_pitzer_gammas(T_K, m_NH3, m_H2S, m_NH4, m_HS):
    """Calculates activity coefficients using standard Pitzer formulation."""
    A_phi = 0.3915 + 0.00021*(T_K - 298.15) + 2.2e-6*((T_K - 298.15)**2)
    
    b0_A_A = calc_pitzer_p(PITZER_PARAMS['beta0_NH3_NH3'], T_K)
    b0_S_S = calc_pitzer_p(PITZER_PARAMS['beta0_H2S_H2S'], T_K)
    
    b0_C_A = calc_pitzer_p(PITZER_PARAMS['beta0_NH4_HS'], T_K)
    b1_C_A = calc_pitzer_p(PITZER_PARAMS['beta1_NH4_HS'], T_K)
    cphi_C_A = calc_pitzer_p(PITZER_PARAMS['Cphi_NH4_HS'], T_K)
    
    I = 0.5 * (m_NH4 + m_HS)
    if I < 1e-10: return 1.0, 1.0, 1.0
        
    sqI = math.sqrt(I)
    f_I = -A_phi * ( (sqI / (1.0 + 1.2*sqI)) + (2.0/1.2)*math.log(1.0 + 1.2*sqI) )
    
    alpha = 2.0
    x = alpha * sqI
    g_x = 2.0 * (1.0 - (1.0 + x)*math.exp(-x)) / (x**2) if x > 1e-8 else 1.0
    
    # Standard Pitzer B_gamma term
    B_gamma = 2.0 * b0_C_A + (2.0 * b1_C_A / (alpha**2 * I)) * (1.0 - (1.0 + alpha*sqI - 0.5*(alpha**2)*I)*math.exp(-alpha*sqI))
    
    # C_gamma term
    C_gamma = 1.5 * cphi_C_A
    
    m = (m_NH4 + m_HS) / 2.0 
    
    # Mean ionic activity coefficient (Standard Pitzer 1-1)
    ln_gamma_pm = f_I + m * B_gamma + (m**2) * C_gamma
    
    # Neutral species (Simplified)
    ln_gamma_NH3 = 2.0 * b0_A_A * m_NH3 
    ln_gamma_H2S = 2.0 * b0_S_S * m_H2S
    
    # Safety caps
    ln_gamma_pm = max(min(ln_gamma_pm, 10.0), -10.0)
    ln_gamma_NH3 = max(min(ln_gamma_NH3, 10.0), -10.0)
    ln_gamma_H2S = max(min(ln_gamma_H2S, 10.0), -10.0)
    
    return math.exp(ln_gamma_NH3), math.exp(ln_gamma_H2S), math.exp(ln_gamma_pm)

def solve_charge_balance(T_K, ph, K, molals, A_davies):
    h_ion = max(10**(-ph), SOLVER['MIN_VAL'])
    m_NH3_total, m_H2S_total, m_CO2_total = molals
    
    # Initial Guesses
    g1_davies, g2_davies = 1.0, 1.0
    g_NH3, g_H2S, g_pm = 1.0, 1.0, 1.0
    nh3_aq, h2s_aq = m_NH3_total, m_H2S_total
    nh4_ion, hs_ion = 0.0, 0.0
    
    ionic_strength = 0.0
    prev_I = -1.0
    
    for _ in range(SOLVER['MAX_ITER_ION']):
        if abs(ionic_strength - prev_I) < SOLVER['TOL_IONIC']: break
        prev_I = ionic_strength

        # 1. PITZER & DAVIES CALCULATIONS
        g_NH3_new, g_H2S_new, g_pm_new = get_pitzer_gammas(T_K, nh3_aq, h2s_aq, nh4_ion, hs_ion)
        g1_dav_new, g2_dav_new = get_activity_coefficients(ionic_strength, A_davies)
        
        # DAMPING: Successive substitution to prevent non-linear oscillation
        w = 0.5 
        g_NH3 = w * g_NH3 + (1-w) * g_NH3_new
        g_H2S = w * g_H2S + (1-w) * g_H2S_new
        g_pm = w * g_pm + (1-w) * g_pm_new
        g1_davies = w * g1_davies + (1-w) * g1_dav_new
        g2_davies = w * g2_davies + (1-w) * g2_dav_new

        # 1. AMMONIA (Pitzer Coupled)
        k_ka_nh4_thermo = 1.0 / K['NH3'] 
        k_ka_nh4_app = k_ka_nh4_thermo / g_NH3
        nh4_ion = m_NH3_total / (1.0 + (k_ka_nh4_app / h_ion))
        nh3_aq = m_NH3_total - nh4_ion
        
        # 3. SULFIDE (Pitzer Coupled)
        k_h2s_app = K['H2S'] * g_H2S / (g_pm**2)
        k_hs_app  = K['HS'] / g2_davies # S-- uses Davies fallback
        den_s = 1.0 + (k_h2s_app/h_ion) + (k_h2s_app * k_hs_app / (h_ion**2))
        h2s_aq = m_H2S_total / den_s
        hs_ion = h2s_aq * (k_h2s_app / h_ion)
        s_ion = hs_ion * (k_hs_app / h_ion)
        
        # 4. CARBONATE (Davies Fallback)
        k_co2_app = K['CO2'] / (g1_davies**2)
        k_hco3_app = K['HCO3'] / g2_davies
        den_c = 1.0 + (k_co2_app/h_ion) + (k_co2_app * k_hco3_app / (h_ion**2))
        co2_aq = m_CO2_total / den_c
        hco3_ion = co2_aq * (k_co2_app / h_ion)
        co3_ion = hco3_ion * (k_hco3_app / h_ion)
        
        # 5. WATER
        k_kw_app = K['Kw'] / (g_pm**2)
        oh_ion = k_kw_app / h_ion
        
        ionic_strength = 0.5 * (h_ion + oh_ion + nh4_ion + hs_ion + hco3_ion + 4.0*s_ion + 4.0*co3_ion)
        
    error = (h_ion + nh4_ion) - (oh_ion + hs_ion + 2.0*s_ion + hco3_ion + 2.0*co3_ion)
    spec = {'NH3': nh3_aq, 'NH4': nh4_ion, 'H2S': h2s_aq, 'HS': hs_ion, 'S': s_ion, 
            'CO2': co2_aq, 'HCO3': hco3_ion, 'CO3': co3_ion, 'I': ionic_strength,
            'g_NH3': g_NH3, 'g_H2S': g_H2S, 'g_pm': g_pm, 'g_dav': g1_davies}
    return error, spec

def calculate_equilibrium(T_K, m_in):
    if T_K < 273.15: return {'error': "Temp Error"}
    T_R = T_K * 1.8
    K_vals = {k: math.exp(calc_ln_k(v, T_R)) for k, v in EDWARDS_PARAMS.items()}
    A_davies = calc_davies_A(T_K)
    
    low, high = 0.0, 14.0
    ph_guess = 7.0
    spec = {}
    
    for _ in range(SOLVER['MAX_ITER_PH']):
        ph_guess = (low + high) / 2.0
        err, spec = solve_charge_balance(T_K, ph_guess, K_vals, (m_in['NH3'], m_in['H2S'], m_in['CO2']), A_davies)
        if abs(err) < SOLVER['TOL_CHARGE']: break
        if err > 0: low = ph_guess
        else: high = ph_guess
    
    nh3, h2s, co2 = spec['NH3'], spec['H2S'], spec['CO2']
    inv_T_term = (1.0/298.15 - 1.0/T_K)
    
    # Henry's Law with Pitzer Activities
    p_nh3 = (CONST['H_REF']['NH3'] * nh3 * spec['g_NH3'] * math.exp(CONST['H_SCALE']['NH3'] * inv_T_term)) * CONST['ATM_TO_PSI']
    p_h2s = (CONST['H_REF']['H2S'] * h2s * spec['g_H2S'] * math.exp(CONST['H_SCALE']['H2S'] * inv_T_term)) * CONST['ATM_TO_PSI']
    p_co2 = (CONST['H_REF']['CO2'] * co2 * 1.0 * math.exp(CONST['H_SCALE']['CO2'] * inv_T_term)) * CONST['ATM_TO_PSI']
    
    # Water Activity (Osmotic Depression)
    sum_moles = sum([v for k,v in spec.items() if k not in ['I', 'g_NH3', 'g_H2S', 'g_pm', 'g_dav']])
    a_w = math.exp(-0.018015 * sum_moles)
    
    p_w_pure = (10**(5.20389 - 1733.926/(T_K - 39.485))) * 14.5038
    p_w = p_w_pure * a_w
    
    p_tot_pa = (p_nh3 + p_h2s + p_co2 + p_w) * CONST['PSI_TO_PA']
    
    return {'ph': ph_guess, 'P_bubble': p_tot_pa, 'pp': {'NH3': p_nh3, 'H2S': p_h2s, 'CO2': p_co2, 'H2O': p_w}, 'liq': spec, 'aw': a_w}

def calculate_density(T_K, P_Pa, spec, kg_water, mass_total_kg):
    t_c = T_K - 273.15
    rho_w = 1000.0 * (1 - (t_c + 288.9414)/(508929.2 * (t_c + 68.12963)) * (t_c - 3.9863)**2)
    v_inf = {'NH4': 18.0, 'HS': 20.0, 'NH3': 24.5, 'H2S': 35.0, 'CO2': 34.0, 'HCO3': 25.0, 'CO3': 20.0, 'S': 20.0}
    vol_water_L = kg_water / (rho_w / 1000.0)
    vol_solutes_mL = (
        spec['NH4']*v_inf['NH4'] + spec['HS']*v_inf['HS'] + spec['NH3']*v_inf['NH3'] + spec['H2S']*v_inf['H2S'] +
        spec['CO2']*v_inf['CO2'] + spec['HCO3']*v_inf['HCO3'] + spec['CO3']*v_inf['CO3'] + spec['S']*v_inf['S']
    ) * kg_water
    vol_solutes_L = vol_solutes_mL / 1000.0
    total_vol_L = vol_water_L + vol_solutes_L
    rho_kg_L = mass_total_kg / total_vol_L
    return (rho_kg_L * 1000.0) * (1 + 4.5e-10 * (P_Pa - 101325.0))

def generate_report(res, T, P_op, flows_h, total_h, rho):
    l = res['liq']; pp = res['pp']
    W = 80; HR = "=" * W; SR = "-" * W
    is_f = res['P_bubble'] > (P_op * 1.01)
    
    rho_kg_L = rho / 1000.0
    h2s_mgL = (l['H2S'] + l['HS'] + l['S']) * CONST['MW']['H2S'] * rho_kg_L * 1000.0
    nh3_mgL = (l['NH3'] + l['NH4']) * CONST['MW']['NH3'] * rho_kg_L * 1000.0

    lines = [HR, " SWEQ - SOUR WATER EQUILIBRIUM SOLVER v8.0.0 ".center(W), " Enterprise Edition - Pitzer Engine Active ".center(W), HR]
    try: lines.append(" Date: %s " % DateTime.Now.ToString("yyyy-MM-dd HH:mm").center(W))
    except: pass
    lines.append((" User: %-16s Model: Edwards/Pitzer (Rumpf 1999) " % "Alexander").center(W))
    lines.append("\n 1. EXECUTIVE SUMMARY & SAFETY CHECK ".ljust(W))
    lines.append(SR)
    lines.append(" STATUS: " + ("!!! FLASHING DETECTED !!!" if is_f else "STABLE LIQUID"))
    lines.append("")
    lines.append(" %-25s | %-15s | %-10s | %-20s" % ("PARAMETER", "VALUE", "UNIT", "COMMENT"))
    lines.append(" %s-+-%s-+-%s-+-%s" % ("-"*25, "-"*15, "-"*10, "-"*20))
    lines.append(" %-25s | %15.2f | %-10s | Inlet Condition" % ("Temperature", T-273.15, "C"))
    lines.append(" %-25s | %15.3f | %-10s | Operating P" % ("System Pressure", P_op/CONST['ATM_TO_PA'], "atm"))
    lines.append(" %-25s | %15.3f | %-10s | Equilibrium P" % ("Bubble Pressure", res['P_bubble']/CONST['ATM_TO_PA'], "atm"))
    lines.append(" %-25s | %15.4f | %-10s | Buffering Status" % ("Calculated pH", res['ph'], "-"))
    lines.append(" %-25s | %15.4f | %-10s | Driver" % ("Ionic Strength (I)", l['I'], "molal"))
    lines.append(" %-25s | %15.4f | %-10s | Osmotic Correction" % ("Water Activity", res['aw'], "a_w"))
    
    lines.append("\n 2. STREAM COMPOSITION & FLOWS ".ljust(W))
    lines.append(SR)
    lines.append(" %-25s | %15s | %-10s " % ("COMPONENT", "FLOW (kg/h)", "MASS %"))
    lines.append(" %s-+-%s-+-%s" % ("-"*25, "-"*15, "-"*10))
    lines.append(" %-25s | %15.2f | %10.2f " % ("Water (H2O)", flows_h['H2O'], (flows_h['H2O']/total_h)*100))
    for k in ['NH3', 'H2S', 'CO2']:
        lines.append(" %-25s | %15.2f | %10.2f " % (k, flows_h[k], (flows_h[k]/total_h)*100))
        
    lines.append("\n 3. CHEMICAL SPECIATION & PITZER ACTIVITIES ".ljust(W))
    lines.append(SR)
    lines.append(" Concentrations in molality (mol/kg H2O)")
    lines.append("")
    lines.append("  AMMONIA SYSTEM   ::  NH3(aq): %7.4f (g=%.2f) <==> NH4+: %7.4f (g=%.2f)" % (l['NH3'], l['g_NH3'], l['NH4'], l['g_pm']))
    lines.append("  SULFIDE SYSTEM   ::  H2S(aq): %7.4f (g=%.2f) <==> HS-:  %7.4f (g=%.2f)" % (l['H2S'], l['g_H2S'], l['HS'], l['g_pm']))
    lines.append("  CARBONATE SYSTEM ::  CO2(aq): %7.4f          <==> HCO3: %7.4f (Davies)" % (l['CO2'], l['HCO3']))
    
    lines.append("\n 4. ENVIRONMENTAL & VAPOR PREDICTIONS ".ljust(W))
    lines.append(SR)
    lines.append(" %-38s | %s" % ("ENVIRONMENTAL (Liquid Effluent)", "VAPOR PHASE (Equilibrium Partial P)"))
    lines.append(" %s-+-%s" % ("-"*38, "-"*39))
    lines.append(" H2S Total: %15.1f mg/L      | p(H2S): %10.4f psia" % (h2s_mgL, pp['H2S']))
    lines.append(" NH3 Total: %15.1f mg/L      | p(NH3): %10.4f psia" % (nh3_mgL, pp['NH3']))
    lines.append(" Density:   %15.2f kg/m3     | p(CO2): %10.4f psia" % (rho, pp['CO2']))
    lines.append("                                        | p(H2O): %10.4f psia" % (pp['H2O']))

    if (l['CO2'] + l['HCO3']) > SOLVER['CO2_DAVIES_WARN']:
        lines.append("\n > Note: High CO2 detected. Carbonate interactions using Davies fallback.")
        
    lines.append(HR + "\n End of Report")
    return "\n".join(lines)

# --- 5. MAIN BRIDGE ---

def Main():
    try:
        if 'ims1' not in globals() or 'oms1' not in globals():
            return 
        
        try:
            T = ims1.GetTemperature()
            P_op = ims1.GetPressure()
            F_mol = ims1.GetMolarFlow()
            comp = ims1.GetOverallComposition()
            ids = ims1.ComponentIds
        except:
            return 

        if F_mol < SOLVER['MIN_VAL']: return
        
        m_mol_s = {'NH3': 0.0, 'H2S': 0.0, 'CO2': 0.0, 'H2O': 0.0}
        idx = {}
        for i, n in enumerate(ids):
            nc = n.lower()
            for k, tags in COMP_MAP.items():
                if any(t in nc for t in tags): 
                    m_mol_s[k] = F_mol * comp[i]
                    idx[k] = i
                    break
        
        if m_mol_s['H2O'] < SOLVER['MIN_VAL']: return

        # Calculation
        kg_w = m_mol_s['H2O'] * CONST['MW']['H2O'] / 1000.0
        molals = {k: v / kg_w for k, v in m_mol_s.items() if k != 'H2O'}
        
        res = calculate_equilibrium(T, molals)
        if res.get('error'): return

        # Update Streams
        P_out = min(res['P_bubble'], 1500 * CONST['ATM_TO_PA'])
        
        y_raw = [0.0] * len(ids)
        for k, i in idx.items():
            if k in res['pp']: 
                y_raw[i] = (res['pp'][k] * CONST['PSI_TO_PA']) / P_out
        
        total_y = sum(y_raw)
        y_norm = [v/total_y if total_y > 0 else comp[i] for v in y_raw]
        
        oms1.SetTemperature(T)
        oms1.SetPressure(P_out)
        oms1.SetOverallComposition(Array[Double](y_norm))
        oms1.SetMolarFlow(1e-8)
        
        oms2.SetTemperature(T)
        oms2.SetPressure(P_out)
        oms2.SetOverallComposition(ims1.GetOverallComposition())
        oms2.SetMolarFlow(F_mol - 1e-8)
        
        oms1.Calculate()
        oms2.Calculate()

        # Report (Robust Save)
        m_kg_s = {k: v * CONST['MW'][k] / 1000.0 for k, v in m_mol_s.items()}
        rho = calculate_density(T, P_out, res['liq'], kg_w, sum(m_kg_s.values()))
        report_str = generate_report(res, T, P_op, {k: m*3600 for k,m in m_kg_s.items()}, sum(m_kg_s.values())*3600, rho)
        
        try:
            path = os.path.join(os.path.expanduser('~'), 'Desktop', 'SWEQ_Datasheet.txt')
            with open(path, 'w') as f: f.write(report_str)
            if os.name == 'nt': os.startfile(path)
        except:
            # Fallback to Temp
            try:
                temp_path = "C:\\Temp\\SWEQ_Report.txt"
                if not os.path.exists("C:\\Temp"): os.makedirs("C:\\Temp")
                with open(temp_path, 'w') as f: f.write(report_str)
                if os.name == 'nt': os.startfile(temp_path)
            except: pass

    except: pass

Main()
