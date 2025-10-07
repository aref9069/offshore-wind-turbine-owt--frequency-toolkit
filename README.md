## offshore-wind-turbine-owt--frequency-toolkit
OWT Frequency Toolkit provides parallel analytical and finite-element solutions for the natural frequencies of a two-segment offshore wind turbine monopile–tower system. The analytical path solves an 8×8 characteristic matrix (EB with axial load), while the numerical path assembles a 2-node EB FEM with geometric stiffness, lumped masses, hydrodynamic added mass, and optional coupled base springs. The package outputs CSV comparison tables, publication-quality figures (PNG+PDF), mesh-convergence plots, added-mass/base-spring sensitivity studies, and mode-shape overlays to validate and tune models for research and engineering use.

# OWT-Freq: Analytical vs FEM Natural Frequencies Toolkit

**Author:** Aref Aasi  
**Keywords:** Offshore Wind Turbine, Natural Frequency, Euler–Bernoulli Beam, FEM, Analytical Model, Mode Shapes  

---

## Overview

**OWT-Freq** is an open-source learning and research toolkit for computing and comparing the **natural frequencies** and **mode shapes** of an **Offshore Wind Turbine (OWT)** structure using two complementary methods:

1. **Analytical Solution**  
   - Based on a **two-segment Euler–Bernoulli beam** model (monopile + tower).  
   - Includes **axial load**, **interface mass (waterline)**, **top mass (nacelle + rotor)**, and **optional base springs**.  
   - Frequencies are found by scanning and solving the determinant of an 8×8 characteristic matrix.

2. **Finite Element (FEM) Solution**  
   - Implements **2-node Euler–Bernoulli beam elements** with consistent mass and stiffness matrices.  
   - Includes **geometric stiffness**, **added hydrodynamic mass**, **lumped masses**, and **base flexibility**.  
   - Extracts eigenvalues and mode shapes via `scipy.linalg.eigh`.

Both methods produce **matching frequencies** within a few percent—ideal for verifying modeling accuracy or teaching vibration fundamentals.

---
Each parameter is entered interactively through the **config.py** module so learners can:
- Experiment with geometry and material values.
- Observe how frequencies change with depth, stiffness, or mass.
- Visualize mesh convergence and mode-shape differences.

---
| File               | Description                                                             |
| ------------------ | ----------------------------------------------------------------------- |
| `config.py`        | Interactive configuration (geometry, materials, loads, springs, search) |
| `utils.py`         | Helper utilities (save figures, compute A/I, override depth/mass)       |
| `analytical.py`    | Determinant-based analytical solver for two-segment EB beam             |
| `fem.py`           | FEM assembly and eigen-solver for two-segment EB beam                   |
| `run_cases.py`     | Example runner (clamped & spring base) producing CSV and figures        |
| `README.md`        | Documentation                                                           |
| `LICENSE`          | Open-source license (MIT recommended)                                   |






## Citation

If you use this repository in teaching or research, please cite:

Aref Aasi (2025). OWT-Freq: Analytical vs FEM Natural Frequency Toolkit for Offshore Wind Turbine Structures. GitHub Repository.
https://github.com/
<your-username>/owt-freq-analytical-vs-fem


## License

This project is released under the MIT License — free for research, teaching, and engineering use.
See LICENSE
 for details.
