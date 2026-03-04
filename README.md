# SWEQ - Sour Water Equilibrium Solver

**Version:** 8.1.0  
**Integration:** DWSIM Python Script Unit Operation  
**License:** GPLv3  

SWEQ is a Python-based thermodynamic solver for calculating the vapor-liquid equilibrium (VLE) and ionic speciation of sour water systems ($NH_3 - H_2S - CO_2 - H_2O$). It is designed to be integrated directly into the DWSIM process simulator.

## Features

* **Thermodynamic Framework:** Based on the weak electrolyte model proposed by Edwards, Newman & Prausnitz (1978).
* **High Concentration Support:** Replaces the standard Davies equation with the Pitzer Activity Coefficient Model (calibrated via Rumpf et al., 1999). This allows the solver to converge in highly concentrated brines (ionic strengths > 6.0 molal) typical of stripper bottoms and atmospheric overheads.
* **Hybrid Activity Model:** * Uses Pitzer parameters for major $NH_4^+/HS^-$ interactions and "salting-out" effects on neutral gases ($\gamma_{NH_3}$, $\gamma_{H_2S}$).
  * Uses the Davies equation as a fallback for minor carbonate species ($HCO_3^-, CO_3^{2-}$).
* **Non-Ideal Phase Corrections:** Includes osmotic mole fraction correction for water activity ($a_w$) depression and the Poynting correction ($v^\infty$) for high-pressure density effects.
* **Numerical Stability:** The charge-balance bisection loop includes a successive substitution damping factor ($w = 0.5$) to prevent non-linear oscillation during the calculation of activity coefficients.
* **Automated Output:** Generates a formatted text report (`SWEQ_Datasheet.txt`) containing phase status, calculated pH, bubble pressure, environmental metrics (mg/L), and specific activity coefficients.

## Usage in DWSIM

1. Add a **Python Script** Unit Operation to your DWSIM flowsheet.
2. Connect one inlet stream (`ims1`) and two outlet streams (liquid `oms1`, vapor `oms2`).
3. Paste the contents of `SWEQ_v8.1.py` into the script editor.
4. Run the flowsheet. The script maps standard components (Ammonia, Hydrogen Sulfide, Carbon Dioxide, and Water), performs the flash calculation, updates the streams, and exports the datasheet to your Desktop or `C:\Temp\`.

## References

* Edwards, T. J., Maurer, G., Newman, J., & Prausnitz, J. M. (1978). *Vapor‐liquid equilibria in multicomponent aqueous solutions of volatile weak electrolytes*. AIChE Journal.
* Rumpf, B., Pérez-Salado Kamps, Á., Sing, R., & Maurer, G. (1999). *Simultaneous solubility of ammonia and hydrogen sulfide in water at temperatures from 313 K to 393 K*. Fluid Phase Equilibria.
