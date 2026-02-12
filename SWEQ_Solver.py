"""
DWSIM Custom Unit Operation: Sour Water Equilibrium (SWEQ) Solver
Thermodynamic Engine: Edwards, Newman, and Prausnitz (1978)

File: SWEQ_Solver.py
Latest Update: Implementation of v6.7 logic (Convergence Check + Robust Scoping)

Features:
- Coupled Phase & Chemical Equilibrium
- Ionic Speciation (NH4+, HS-, HCO3-)
- Environmental Reporting (ppm)
- CO2 Saturation Clamp & Threshold Logic
"""

import math
import os
import System
import clr
from System import Array, Double

try:
    clr.AddReference("System.Windows.Forms")
    from System.Windows.Forms import MessageBox
except: pass

# --- Constants & Configuration ---
P_CRIT_CO2_PSIA = 1070.6
CO2_THRESHOLD = 1e-6 
MAX_ITER = 200
TOL_CONVERGENCE = 1e-9

EDWARDS_COEFFS = {
    'NH3':   [1.587,   11160.0,   0.0,         0.0,        0.0],
    'H2S':   [-293.88, 683858.0, -6.27125e8,  2.555e11,   -3.91757e13],
    'HS':    [-220.07, 258273.0, -1.8396e8,   6.809e10,   -1.0267e13],
    'CO2':   [-241.79, 536256.0, -4.8123e8,   1.94e11,    -2.96445e13],
    'HCO3':  [-294.74, 364385.0, -2.841e8,    1.23323e11, -2.0759e13],
    'Kw':    [39.5554, -177822.0, 1.843e8,    -0.8541e11, 1.4292e13]
}
SYNONYMS = {
    'NH3': ['ammonia', 'nh3', 'amonia', 'nitrogen hydride'],
    'H2S': ['hydrogen sulfide', 'h2s', 'sulfeto'],
    'CO2': ['carbon dioxide', 'co2', 'dioxido'],
    'H2O': ['water', 'h2o', 'agua']
}
MW = { 'NH3': 17.031, 'H2S': 34.08, 'CO2': 44.01, 'H2O': 18.015, 'DEFAULT': 28.0 }
PSI_TO_PA = 6894.76
ATM_TO_PA = 101325.0

def get_ln_k(coeffs, TR):
    """Calculates ln(K) using temperature-dependent polynomial."""
    A, B, C, D, E = coeffs
    return A + B/TR + C/(TR**2) + D/(TR**3) + E/(TR**4)

def run_equilibrium(TK, m_NH3, m_H2S, m_CO2):
    TR = TK * 1.8 # Convert Kelvin to Rankine
    K = {k: math.exp(get_ln_k(v, TR)) for k, v in EDWARDS_COEFFS.items()}
    
    # Bisection Solver Initialization
    ph = 9.0; ph_lo = 0.0; ph_hi = 14.0
    
    # Explicit initialization (Safety against scope leaks)
    nh3_aq = nh4_ion = h2s_aq = hs_ion = co2_aq = hco3_ion = 0.0

    for i in range(MAX_ITER):
        h_val = max(10**(-ph), 1e-16) # Prevent underflow
        oh_val = K['Kw'] / h_val
        
        # Ammonia Balance
        nh3_aq = m_NH3 / (1 + K['NH3'] * h_val)
        nh4_ion = nh3_aq * K['NH3'] * h_val
        
        # Sulfide Balance
        d_h2s = 1 + (K['H2S']/h_val) + (K['H2S']*K['HS']/(h_val**2))
        h2s_aq = m_H2S / d_h2s
        hs_ion = h2s_aq * (K['H2S']/h_val)
        
        # Carbonate Balance
        d_co2 = 1 + (K['CO2']/h_val) + (K['CO2']*K['HCO3']/(h_val**2))
        co2_aq = m_CO2 / d_co2
        hco3_ion = co2_aq * (K['CO2']/h_val)
        
        # Charge Balance (Electroneutrality)
        balance = (h_val + nh4_ion) - (oh_val + hs_ion + hco3_ion)
        
        # Convergence Check (Optimization)
        if abs(balance) < TOL_CONVERGENCE:
            break

        if balance > 0: ph_lo = ph
        else: ph_hi = ph
        ph = (ph_lo + ph_hi) / 2.0

    # --- Partial Pressure Calculations (Henry's Law + Interaction Parameters) ---
    
    # NH3 Correlation (Edwards, 1978)
    ln_h_nh3 = 178.339 - 15517.91/TR - 25.6767*math.log(TR) + 0.01966*TR + (131.4/TR - 0.1682)*nh3_aq + 0.06*(2*m_CO2 + m_H2S)
    p_nh3 = math.exp(ln_h_nh3) * nh3_aq
    
    # H2S Correlation (Edwards, 1978)
    ln_h_h2s = 8.8 + 0.015 * (TR - 560)
    p_h2s = math.exp(ln_h_h2s + (-0.05*nh3_aq + (0.965 - 486.0/TR)*m_CO2)) * h2s_aq
    
    # CO2 Logic with Threshold and Safety Clamp
    clamped = False
    p_co2 = 0.0
    if m_CO2 > CO2_THRESHOLD:
        # CO2 Correlation (Edwards, 1978)
        ln_h_co2 = 301.68 - 34096.6/TR + 1.2285e8/(TR**2) - 6.4752e10/(TR**3) + 1.1557e13/(TR**4)
        if (m_NH3 + m_CO2 + m_H2S) > 1e-6: ln_h_co2 += -0.09 * (nh3_aq - m_CO2 - m_H2S)
        try:
            val_check = ln_h_co2 + math.log(max(co2_aq, 1e-20))
            if val_check > 7.0: p_co2 = P_CRIT_CO2_PSIA; clamped = True
            else: p_co2 = math.exp(ln_h_co2) * co2_aq
        except: p_co2 = P_CRIT_CO2_PSIA; clamped = True

    # Water Vapor Pressure (Kell Correlation approximation)
    p_w = (10**(5.20389 - 1733.926/(TK - 39.485))) * 14.5038
    p_tot_pa = (p_nh3 + p_h2s + p_co2 + p_w) * PSI_TO_PA

    return {'ph': ph, 'p_total': p_tot_pa, 'clamped': clamped,
            'pp_psia': {'NH3': p_nh3, 'H2S': p_h2s, 'CO2': p_co2, 'H2O': p_w},
            'liq': {'NH3_aq': nh3_aq, 'NH4_ion': nh4_ion, 'H2S_aq': h2s_aq, 'HS_ion': hs_ion, 'CO2_aq': co2_aq, 'HCO3_ion': hco3_ion}}

def Solve_SWEQ():
    try:
        T_in_K = ims1.GetTemperature(); P_in_Pa = ims1.GetPressure()
        f_mol = ims1.GetMolarFlow(); total_mass_flow = ims1.GetMassFlow()
        if f_mol < 1e-9: return
        compositions = ims1.GetOverallComposition(); comp_ids = ims1.ComponentIds
        rho_initial = ims1.GetMassDensity() if hasattr(ims1, 'GetMassDensity') else 0.0
        h_in = ims1.GetMassEnthalpy() if hasattr(ims1, 'GetMassEnthalpy') else 0.0
    except: return

    mol_in = {'NH3': 0.0, 'H2S': 0.0, 'CO2': 0.0, 'H2O': 0.0}
    indices = {'NH3': -1, 'H2S': -1, 'CO2': -1, 'H2O': -1}
    active_idx = []; feed_data = []; mw_list = [0.0] * len(comp_ids)
    
    for i, name in enumerate(comp_ids):
        n_clean = name.lower().strip(); mol_flow_i = f_mol * compositions[i]
        mw_i = MW['DEFAULT']
        for key, synonyms in SYNONYMS.items():
            if any(s in n_clean for s in synonyms):
                mw_i = MW[key]; mol_in[key] = mol_flow_i; indices[key] = i; active_idx.append(i); break
        mw_list[i] = mw_i
        mass_flow_i = mol_flow_i * (mw_i / 1000.0)
        feed_data.append({'name': name, 'mol_frac': compositions[i], 'mass_flow_kgs': mass_flow_i})
    
    if mol_in['H2O'] < 1e-7: return
    calc_total_mass_kgs = sum(d['mass_flow_kgs'] for d in feed_data)

    kg_water = mol_in['H2O'] * 0.018015
    m = {k: v/kg_water for k,v in mol_in.items() if k != 'H2O'}
    res = run_equilibrium(T_in_K, m['NH3'], m['H2S'], m['CO2'])

    P_calc_Pa = min(res['p_total'], 1500 * ATM_TO_PA)
    pp = res['pp_psia']
    y_raw = [0.0] * len(comp_ids)
    for key, idx in indices.items():
        if idx != -1: y_raw[idx] = (pp[key] * PSI_TO_PA) / P_calc_Pa
    for i in range(len(comp_ids)):
        if i not in active_idx: y_raw[i] = compositions[i]
    sum_y = sum(y_raw)
    final_y = [val/sum_y for val in y_raw] if sum_y > 0 else compositions
    mw_vap_avg = sum(final_y[i] * mw_list[i] for i in range(len(comp_ids)))

    wt_tot = (1000 + (m['NH3']*17.03 + m['H2S']*34.08 + m['CO2']*44.01))
    wt_H2O, wt_NH3, wt_H2S, wt_CO2 = 1000/wt_tot*100, (m['NH3']*17.03)/wt_tot*100, (m['H2S']*34.08)/wt_tot*100, (m['CO2']*44.01)/wt_tot*100

    oms1.SetTemperature(T_in_K); oms1.SetPressure(P_calc_Pa); oms1.SetOverallComposition(Array[Double](final_y)); oms1.SetMolarFlow(1e-8)
    oms2.SetTemperature(T_in_K); oms2.SetPressure(P_calc_Pa); oms2.SetOverallComposition(Array[Double](list(compositions))); oms2.SetMolarFlow(f_mol - 1e-8)
    oms1.Calculate(); oms2.Calculate()

    try:
        desktop = os.path.join(os.environ['USERPROFILE'], 'Desktop')
        fname = os.path.join(desktop, "SWEQ_Datasheet.txt")
        liq = res['liq']
        def header(t): return "="*70 + "\n" + t.center(70) + "\n" + "="*70 + "\n"
        def sep(): return "-"*70 + "\n"
        
        rho_final = rho_initial
        density_tag = "(DWSIM)"
        if rho_final < 1.0 or rho_final > 2000.0:
            tc = T_in_K - 273.15
            rho_final = (999.83952 + 16.945176*tc - 7.98704e-3*tc**2 - 46.170461e-6*tc**3) / (1 + 16.879850e-3*tc)
            density_tag = "(Estimated)"

        is_flashing = P_calc_Pa > P_in_Pa
        status_msg = "*** WARNING: UNSTABLE (FLASHING) ***" if is_flashing else "STABLE LIQUID"

        r  = header("SOUR WATER EQUILIBRIUM (SWEQ) DATASHEET")
        r += " Date: %s\n" % System.DateTime.Now
        r += " Model: Edwards, Newman, and Prausnitz (1978)\n"
        
        r += sep() + " 1. FEED STREAM SPECIFICATIONS\n" + sep()
        r += " Temperature:      %8.2f C\n" % (T_in_K-273.15)
        r += " Pressure (Op):    %8.3f atm\n" % (P_in_Pa/101325.0)
        r += " Mass Flow Total:  %8.2f kg/h (DWSIM Ref)\n" % (total_mass_flow*3600)
        r += " Spec. Enthalpy:   %8.2f kJ/kg\n" % h_in
        r += " Status:           %s\n" % status_msg
        r += "\n Composition (Input):\n"
        r += " Component       | Mass % | Mass Flow (kg/h) | Mol % \n"
        r += " ----------------+--------+------------------+--------\n"
        for d in feed_data:
            m_pct = (d['mass_flow_kgs'] / calc_total_mass_kgs) * 100
            r += " %-15s | %6.2f | %16.2f | %6.2f\n" % (d['name'], m_pct, d['mass_flow_kgs']*3600, d['mol_frac']*100)
        
        r += sep() + " 2. EQUILIBRIUM OUTPUT (BUBBLE POINT)\n" + sep()
        r += " pH Calculated:    %8.4f\n" % res['ph']
        r += " Bubble Pressure:  %8.3f atm (%8.2f psia)\n" % (P_calc_Pa/101325.0, P_calc_Pa/6894.76)
        if res['clamped']: r += " NOTE: CO2 partial pressure reached physical saturation limit.\n"
        
        r += "\n ENVIRONMENTAL CONCENTRATIONS (Liquid):\n"
        r += "  > H2S Total: %10.2f mg/L (ppm)\n" % ((wt_H2S/100)*1e6)
        r += "  > NH3 Total: %10.2f mg/L (ppm)\n" % ((wt_NH3/100)*1e6)

        r += "\n VAPOR PHASE (Partial Pressures):\n"
        r += "  > H2S: %8.4f psia\n" % pp['H2S']
        r += "  > NH3: %8.4f psia\n" % pp['NH3']
        r += "  > CO2: %8.4f psia\n" % pp['CO2']
        r += "  > H2O: %8.4f psia\n" % pp['H2O']
        
        r += sep() + " 3. VAPOR PHASE SPECS (Estimated)\n" + sep()
        r += " Avg Mol Weight:   %8.2f g/mol\n" % mw_vap_avg
        r += " Composition (Mol %):\n"
        for i, val in enumerate(final_y):
            if val > 1e-4: r += "  > %-15s: %6.2f %%\n" % (comp_ids[i], val*100)
        
        r += sep() + " 4. LIQUID PHASE SPECS (Product)\n" + sep()
        r += " Density:          %8.2f kg/m3 %s\n" % (rho_final, density_tag)
        r += " Composition (Weight %%): Water %5.2f | NH3 %5.2f | H2S %5.2f | CO2 %5.2f\n" % (wt_H2O, wt_NH3, wt_H2S, wt_CO2)
        
        r += "\n Ionic Detail (Speciation):\n"
        r += "  [NH3 Total] -> NH3: %.3f | NH4+: %.3f\n" % (liq['NH3_aq'], liq['NH4_ion'])
        r += "  [H2S Total] -> H2S: %.3f | HS-:  %.3f\n" % (liq['H2S_aq'], liq['HS_ion'])
        if mol_in['CO2'] > CO2_THRESHOLD:
            r += "  [CO2 Total] -> CO2: %.3f | HCO3-: %.3f\n" % (liq['CO2_aq'], liq['HCO3_ion'])
        r += sep()
        
        with open(fname, "w", encoding='utf-8') as f: f.write(r)
        os.startfile(fname)
        ims1.FlowSheet.WriteMessage("Report saved: " + fname)
    except Exception as e:
        ims1.FlowSheet.WriteMessage("Report Error: " + str(e))

Solve_SWEQ()
