# SWEQ-DWSIM: Sour Water Equilibrium Solver

![Platform](https://img.shields.io/badge/platform-DWSIM-green.svg) ![License](https://img.shields.io/badge/license-GPLv3-blue.svg)

**A Semi-Rigorous Thermodynamic Unit Operation for DWSIM based on the Edwards, Newman, and Prausnitz (1978) Electrolyte Model with Activity Coefficient Corrections.**

## Overview
Simulating Sour Water systems (NH3 - H2S - CO2 - H2O) poses a significant challenge in process engineering. Standard Cubic Equations of State often fail to predict the pH-dependent behavior of weak electrolytes, leading to inaccurate Bubble Point, pH, and corrosion predictions.

SWEQ-DWSIM addresses this by implementing a chemical equilibrium solver directly within DWSIM as a Custom Python Unit Operation. The solver includes non-ideal liquid phase corrections suitable for industrial screening and process safety analysis.

## Key Features

### Thermodynamic Framework
* **Electrolyte Model:** Implements the Edwards, Newman, and Prausnitz (1978) correlations.
* **Activity Coefficients:** Utilizes the Davies Equation with temperature-dependent dielectric parameters to correct for non-ideality in the liquid phase.
* **Iterative Solver:** Features a self-consistent loop that converges Ionic Strength and Activity Coefficients simultaneously with the charge balance.

### Extended Chemistry
* **Multi-Step Dissociation:** Handles the full equilibrium chain for accurate high-pH simulation, including Ammonia, Hydrogen Sulfide, and Carbon Dioxide dissociation.
* **Alkaline Systems:** Capable of simulating caustic scrubbers (pH > 10) due to the inclusion of Carbonate and Sulfide ions.

### Operational Safety & Compliance
* **Safety Flash:** Calculates the Incipient Bubble Pressure to detect flashing risks in transfer lines.
* **Environmental Reporting:** Automatically calculates Total H2S and NH3 in mg/L (corrected for fluid density) for environmental compliance verification.
* **Stability Warnings:** Alerts the user if the Ionic Strength exceeds the valid range for the Davies equation (I > 0.5m).

## Installation & Usage

1.  Download the `SWEQ_Solver.py` file from this repository.
2.  Open DWSIM and create a new Flowsheet or open an existing one.
3.  Add a **"Custom Unit Operation (Python Script)"** from the object palette.
4.  Connect inlet and outlet material streams.
5.  Open the script editor in DWSIM, paste the code from `SWEQ_Solver.py`, and save.
6.  Run the simulation. The solver will compute the equilibrium and generate a detailed text datasheet.

## Engineering Limitations

Users must be aware of the following thermodynamic boundaries:

1.  **Ionic Strength Limit:** The Davies equation is rigorously valid for ionic strengths below 0.5 m. For high-salinity brines, accuracy may decrease.
2.  **Vapor Phase Assumption:** The model assumes an ideal gas phase (fugacity coefficient = 1). For pressures significantly above 10-15 atm, this may overestimate the bubble point.
3.  **Water Activity:** The model assumes the solvent follows Raoult's Law without activity correction.

## Technical References
1.  Edwards, T. J., Newman, J., & Prausnitz, J. M. (1978). *Thermodynamics of vapor-liquid equilibria in neutral and aqueous solutions of volatile weak electrolytes*. AIChE Journal.
2.  Davies, C. W. (1962). *Ion Association*. Butterworths.

## License
This project is licensed under the **GNU General Public License v3.0** - see the [LICENSE](LICENSE) file for details.
