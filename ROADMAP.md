# Development Roadmap

This document outlines the development status and future technical goals for the SWEQ-DWSIM project.

## Current Capabilities

### Base Thermodynamics
* Implementation of Edwards (1978) correlations for equilibrium constants.
* Robust Bisection method for pH convergence with bracket validation.
* Critical Pressure clamping for CO2 to prevent numerical divergence.

### Advanced Electrolyte Modeling
* Implementation of Davies Equation for activity coefficients in non-ideal liquid phases.
* Dynamic calculation of dielectric parameters based on temperature.
* Iterative solver loop to converge Ionic Strength and Activity Coefficients self-consistently.
* Activation of Carbonate and Sulfide equilibria for high-pH system simulation.
* Rigorous mg/L calculation correcting for fluid density.

---

## Future Development Goals

### Short-Term Objectives
* **Water Activity Correction:** Implementation of Raoult's Law correction to account for vapor pressure depression in concentrated solutions.
* **UI Integration:** Development of a custom GUI window within DWSIM to display results natively.

### Long-Term Objectives
* **Vapor Phase Non-Ideality:** Implementation of Peng-Robinson EOS to calculate fugacity coefficients for high-pressure systems (> 20 atm).
* **High Salinity Support:** Replacement of the Davies Equation with Pitzer or eNRTL models to handle brines and high ionic strength (I > 6m).
* **Flash Algorithms:** Implementation of a full Isothermal Flash algorithm (Rachford-Rice) to calculate Vapor Fraction (V/F) in addition to Incipient Bubble Point.
* **Energy Balance:** Estimation of Heat of Reaction/Absorption for stripper reboiler optimization.
