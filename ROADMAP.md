# 🗺️ SWEQ Project Roadmap

SWEQ has evolved from a basic Davies-based equilibrium script into a robust Pitzer-driven simulator. The following roadmap outlines the strategic direction for future releases, focusing on expanding the thermodynamic limits, adding new components, and improving interoperability.

### 🟢 Phase 1: Stabilization & Edge Cases (Current - v8.x)
- [x] **Pitzer Model Integration:** Implement Rumpf (1999) parameters for $NH_3-H_2S$ up to 120°C.
- [x] **High Molality Damping:** Add numerical relaxation to prevent solver crashing at $I > 5.0$ molal.
- [x] **Water Activity Suppression:** Osmotic correction for vapor pressure in heavy brines.
- [ ] **Carbonate Matrix Expansion:** Research and implement specific Pitzer interaction parameters for $CO_2 / HCO_3^-$ to replace the current Davies fallback, fully unifying the ternary gas system.
- [ ] **Density Engine Overhaul:** Replace empirical partial molar volumes ($v^\infty$) with rigorous Redlich-Meyer equations of state for high-pressure density predictions.

### 🟡 Phase 2: The "Danish School" Upgrade (v9.0)
- [ ] **Extended UNIQUAC Framework:** Transition the core thermodynamic engine from Pitzer to the **Extended UNIQUAC** model (Thomsen & Rasmussen, 1999). 
    * *Goal:* Achieve flawless calculations from pure water solvent all the way to molten salt limits without the polynomial explosion risks inherent to the Pitzer model.
- [ ] **Newton-Raphson Multivariable Solver:** Replace the nested bisection loops with an analytical Jacobian matrix solver to handle the highly non-linear surface area/volume fractions ($\theta$) of the UNIQUAC model.
- [ ] **Amine Sweeting Integration:** Add parameters for standard alkanolamines (MEA, DEA, MDEA) to allow SWEQ to simulate full gas sweetening absorption units, not just sour water strippers.

### 🔴 Phase 3: Standalone Architecture & Interoperability (v10.0+)
- [ ] **CAPE-OPEN Compliance:** Package the SWEQ thermodynamic engine as a CAPE-OPEN Property Package, allowing it to be natively plugged into Aspen HYSYS, Aspen Plus, PRO/II, and other commercial simulators.
- [ ] **Standalone GUI:** Develop a lightweight desktop application (using PyQt or CustomTkinter) for quick, on-the-fly sour water flash calculations without needing to open a full DWSIM flowsheet.
- [ ] **Kinetic/Reactive Distillation:** Implement kinetic rate equations for $CO_2$ hydration to support rigorous non-equilibrium tray-by-tray distillation modeling.
