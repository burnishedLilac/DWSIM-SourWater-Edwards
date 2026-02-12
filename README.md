# SWEQ-DWSIM: Sour Water Equilibrium Solver

**A Rigorous Thermodynamic Unit Operation for DWSIM based on the Edwards, Newman, and Prausnitz (1978) Electrolyte Model.**

 Overview
Simulating **Sour Water** systems ($NH_3 - H_2S - CO_2 - H_2O$) poses a significant challenge in process engineering. Standard Cubic Equations of State (like Peng-Robinson or SRK) often fail to predict the pH-dependent behavior of weak electrolytes, leading to inaccurate Bubble Point and VLE calculations.

**SWEQ-DWSIM** solves this by implementing a rigorous chemical equilibrium solver directly within DWSIM as a Custom Python Unit Operation.

 Key Features

* **Coupled Phase & Chemical Equilibrium:** Simultaneously solves charge balance (pH) and phase equilibrium.
* **Ionic Speciation:** Accurately calculates concentrations of molecular species ($NH_{3(aq)}$, $H_2S_{(aq)}$) vs. ionic species ($NH_4^+$, $HS^-$, $HCO_3^-$).
* **Environmental Metrics:** Automatically reports liquid concentrations in **mg/L (ppm)** for easy compliance checking.
* **Numerical Robustness:**
    * Includes a physical saturation clamp for $CO_2$ (at ~1070 psia) to prevent mathematical divergence in super-saturated streams.
    * Threshold logic to handle trace components without instability.
* **Automated Reporting:** Generates a professional `.txt` datasheet on the user's desktop after every run.

 Installation & Usage

1.  **Download:** Get the `SWEQ_Solver_v6.7.py` file from this repository.
2.  **Open DWSIM:** Create a new Flowsheet or open an existing one.
3.  **Add Unit Op:** Drag and drop a **"Custom Unit Operation (Python Script)"** from the object palette.
4.  **Connect Streams:** Connect inlet and outlet material streams.
5.  **Paste Code:** Open the script editor in DWSIM, paste the code from `SWEQ_Solver_v6.7.py`, and save.
6.  **Run Simulation:** The solver will compute the equilibrium and open the datasheet automatically.
