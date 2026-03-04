# SWEQ - Sour Water Equilibrium Solver 💧⚙️

**Version:** 8.1.0 (Enterprise Pitzer Edition)  
**Integration:** DWSIM Python Script Unit Operation  
**License:** GPLv3  

SWEQ is an industrial-grade, open-source thermodynamic solver designed to calculate the vapor-liquid equilibrium (VLE) and ionic speciation of highly concentrated sour water systems ($NH_3 - H_2S - CO_2 - H_2O$). Built as a native Python Script for DWSIM, it bridges the gap between academic electrolyte models and heavy-duty refinery simulation.

## 🚀 Key Features

* **Extreme Concentration Engine:** Bypasses the traditional 0.5 - 0.8 molal physical limits of the extended Debye-Hückel/Davies equations. SWEQ v8.1 stably converges at ionic strengths exceeding **6.0+ molal**, handling the harshest stripper bottoms and atmospheric overhead condensates.
* **State-of-the-Art Thermodynamics:** * **Equilibrium:** Edwards, Newman & Prausnitz (1978) weak electrolyte framework.
  * **Activity Coefficients:** Implementation of the specific **Pitzer Activity Coefficient Model** calibrated for sour water (Rumpf et al., 1999) for $NH_4^+/HS^-$ pairs.
  * **Hybrid Fallback:** Seamlessly utilizes the Davies equation as a thermodynamic fallback for minor carbonate species ($HCO_3^-, CO_3^{2-}$).
* **Advanced Non-Ideal Physics:** * Exact modeling of the "salting-out" effect for dissolved neutral gases ($\gamma_{NH_3}$ and $\gamma_{H_2S}$).
  * Osmotic mole fraction correction for water activity ($a_w$) depression in heavy brines.
  * Poynting correction ($v^\infty$) for high-pressure operations.
* **Robust Numerical Solver:** Features a custom charge-balance bisection loop fortified with a **Successive Substitution Damping algorithm** ($w = 0.5$). This acts as a mathematical shock absorber, completely eliminating non-linear oscillation and `Math Domain` errors during extreme stoichiometric imbalances.
* **Automated Engineering Reports:** Generates comprehensive, scannable plain-text datasheets containing environmental metrics (mg/L), pH, phase statuses, bubble point pressures, and detailed chemical speciation with individual Pitzer coefficients.

## 🛠️ Usage in DWSIM

1. Add a **Python Script** Unit Operation to your DWSIM flowsheet.
2. Connect an inlet stream (`ims1`) and two outlet streams (liquid `oms1`, vapor `oms2`).
3. Paste the contents of `SWEQ_v8.1.py` into the script editor.
4. Run the flowsheet. The solver will automatically detect the components (case-insensitive mapping for Ammonia, Hydrogen Sulfide, Carbon Dioxide, and Water), perform the flash calculation, update the streams, and save the detailed `SWEQ_Datasheet.txt` to your Desktop or `C:\Temp\`.

## 📚 References
* Edwards, T. J., Maurer, G., Newman, J., & Prausnitz, J. M. (1978). *Vapor‐liquid equilibria in multicomponent aqueous solutions of volatile weak electrolytes*. AIChE Journal.
* Rumpf, B., Pérez-Salado Kamps, Á., Sing, R., & Maurer, G. (1999). *Simultaneous solubility of ammonia and hydrogen sulfide in water at temperatures from 313 K to 393 K*. Fluid Phase Equilibria.

---
