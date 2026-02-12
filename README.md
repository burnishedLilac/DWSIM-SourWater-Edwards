# SWEQ-DWSIM: Sour Water Equilibrium Solver

**A Rigorous Thermodynamic Unit Operation for DWSIM based on the Edwards, Newman, and Prausnitz (1978) Electrolyte Model.**

## 📌 Overview
Simulating **Sour Water** systems ($NH_3 - H_2S - CO_2 - H_2O$) poses a significant challenge in process engineering. Standard Cubic Equations of State (like Peng-Robinson or SRK) often fail to predict the pH-dependent behavior of weak electrolytes, leading to inaccurate Bubble Point and VLE calculations.

**SWEQ-DWSIM** solves this by implementing a rigorous chemical equilibrium solver directly within DWSIM as a Custom Python Unit Operation.

## 🚀 Key Features

* **Coupled Phase & Chemical Equilibrium:** Simultaneously solves charge balance (pH) and phase equilibrium.
* **Ionic Speciation:** Accurately calculates concentrations of molecular species ($NH_{3(aq)}$, $H_2S_{(aq)}$) vs. ionic species ($NH_4^+$, $HS^-$, $HCO_3^-$).
* **Environmental Metrics:** Automatically reports liquid concentrations in **mg/L (ppm)** for easy compliance checking.
* **Numerical Robustness:**
    * Includes a physical saturation clamp for $CO_2$ (at ~1070 psia) to prevent mathematical divergence in super-saturated streams.
    * Threshold logic to handle trace components without instability.
* **Automated Reporting:** Generates a professional `.txt` datasheet on the user's desktop after every run.

## 🛠️ Installation & Usage

1.  **Download:** Get the `SWEQ_Solver_v6.7.py` file from this repository.
2.  **Open DWSIM:** Create a new Flowsheet or open an existing one.
3.  **Add Unit Op:** Drag and drop a **"Custom Unit Operation (Python Script)"** from the object palette.
4.  **Connect Streams:** Connect inlet and outlet material streams.
5.  **Paste Code:** Open the script editor in DWSIM, paste the code from `SWEQ_Solver_v6.7.py`, and save.
6.  **Run Simulation:** The solver will compute the equilibrium and open the datasheet automatically.

## 📊 Sample Output
*Typical output for a stripped sour water stream at 25°C:*

```text
======================================================================
             SOUR WATER EQUILIBRIUM (SWEQ) DATASHEET - v6.6            
======================================================================
 Date: 02/11/2026 08:30:00
 Model: Edwards, Newman, and Prausnitz (1978)
----------------------------------------------------------------------
 1. FEED STREAM SPECIFICATIONS
----------------------------------------------------------------------
 Temperature:         25.00 C
 Pressure (Op):       3.000 atm
 Status:           STABLE LIQUID

 Composition (Input):
 Component       | Mass % | Mass Flow (kg/h) | Mol % 
 ----------------+--------+------------------+--------
 Water           |  87.07 |          3600.05 |  89.72
 Ammonia         |   5.93 |           245.25 |   6.47
 Hydrogen sulfide|   7.00 |           289.54 |   3.81
----------------------------------------------------------------------
 2. EQUILIBRIUM OUTPUT (BUBBLE POINT)
----------------------------------------------------------------------
 pH Calculated:      9.5657
 Bubble Pressure:     2.328 atm (   34.22 psia)

 ENVIRONMENTAL CONCENTRATIONS (Liquid):
  > H2S Total:   70025.68 mg/L (ppm)
  > NH3 Total:   59308.97 mg/L (ppm)

 VAPOR PHASE (Partial Pressures):
  > H2S:  33.2340 psia
  > NH3:   0.5242 psia
----------------------------------------------------------------------
 4. LIQUID PHASE SPECS (Product)
----------------------------------------------------------------------
 Composition (Weight %): Water 87.07 | NH3  5.93 | H2S  7.00

 Ionic Detail (Speciation):
  [NH3 Total] -> NH3: 1.648 | NH4+: 2.352
  [H2S Total] -> H2S: 0.008 | HS-:  2.352
----------------------------------------------------------------------
