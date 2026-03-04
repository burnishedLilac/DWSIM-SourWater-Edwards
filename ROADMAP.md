# SWEQ Project Roadmap

The evolution of SWEQ follows a rigorous chemical engineering progression: from a basic aqueous speciation script to a full-fledged multiphase process simulator. The roadmap below outlines the strategic direction for future releases.

### Phase 1: Core Engine Stabilization (Current - v8.1)
- [x] **Pitzer Model Integration:** Implemented specific interaction parameters for the $NH_3-H_2S$ system at extreme concentrations.
- [x] **Convergence Damping:** Added Successive Substitution relaxation to prevent solver oscillation at high ionic strengths.
- [x] **Activity-Coupled Henry's Law:** Exact integration of liquid non-ideality into vapor partial pressure calculations.

### Phase 2: Thermodynamic Completeness (v8.5)
- [ ] **Density Engine Overhaul:** Replace infinite dilution empirical volumes with the rigorous Laliberté (2007/2009) density model for aqueous electrolytes.
- [ ] **Carbonate Matrix Expansion:** Implement Pitzer interaction parameters for $CO_2 / HCO_3^-$ to replace the current Davies fallback network.
- [ ] **Rigorous Water Activity:** Upgrade the current osmotic mole fraction approximation to full Pitzer osmotic coefficient ($\phi$) calculations.

### Phase 3: Rigorous VLE & Energy Balances (v9.0)
- [ ] **Isothermal Flash (PT-Flash):** Implement the Rachford-Rice algorithm to calculate exact vapor-liquid fractions ($V/F$) and composition when system pressure drops below the bubble point.
- [ ] **Vapor Phase Non-Ideality:** Integrate a cubic Equation of State (e.g., Peng-Robinson or SRK) to calculate vapor fugacity coefficients ($\hat{\phi}_i$), replacing the current ideal gas assumption at high pressures.
- [ ] **Enthalpy & Heat of Reaction:** Formulate partial molar enthalpies and heat of reaction for ionic species to enable Adiabatic Flash (PH-Flash) and precise temperature drop predictions across expansion valves.

### Phase 4: Extreme Flowsheet Scenarios (v10.0)
- [ ] **Solid Phase Precipitation:** Add solubility product ($K_{sp}$) calculations to predict the exact deposition temperature of solid Ammonium Bisulfide ($NH_4HS$) salts in overhead condensers.
- [ ] **Extended Contaminants:** Introduce thermodynamic parameters for typical deep-conversion sour water impurities, including Hydrogen Cyanide (HCN), Hydrochloric Acid (HCl), and organic acids (Formic/Acetic).

### Phase 5: Interoperability & Next-Gen Framework (v11.0+)
- [ ] **CAPE-OPEN Compliance:** Package the SWEQ thermodynamic engine as a CAPE-OPEN Property Package for native plug-and-play integration with commercial simulators (Aspen HYSYS, PRO/II).
- [ ] **Extended UNIQUAC Migration:** Explore transitioning the core activity model from Pitzer to the Extended UNIQUAC framework (Thomsen) for seamless integration of heavy alkanolamines (MEA, DEA, MDEA) in gas sweetening units.
