## 🗺️ Development Roadmap

 v6.x (Current Stable)
- [x] Edwards (1978) Electrolyte Model implementation.
- [x] Robust pH solver with bracket validation.
- [x] Environmental reporting (ppm).
- [x] CO2 Critical Pressure clamping.

  v7.0 (Upcoming High pH Update)
- [ ] **Extended pH Range:** Inclusion of Carbonate ion ($CO_3^{2-}$) equilibrium for caustic scrubber simulation (pH > 10).
- [ ] **Activity Coefficients:** Implementation of Davies Equation for ionic strength correction.
- [ ] **Code Refactoring:** Modular functions for better maintainability.

Future Features (v8.0+)
- [ ] **Rigorous PT-Flash:** Vapor fraction calculation ($V/F$) for 2-phase separators.
- [ ] **Non-Ideal Gas Phase:** Peng-Robinson EOS for vapor fugacity coefficients.
- [ ] **High Salinity:** Corrections for brine systems (salting-out effect).
