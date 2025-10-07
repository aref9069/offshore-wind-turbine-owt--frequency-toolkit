"""
config.py — Centralized model configuration (SI units)
------------------------------------------------------
All tunables live here so analyses stay reproducible.

UNITS (SI):
- Length: m
- Force: N
- Moment: N·m
- Mass: kg (line mass kg/m)
- Density: kg/m^3
- Stiffness: N/m (KL), N·m/rad (KR)
- Modulus: Pa
- Frequency output: Hz

Sections:
- Geometry (L, L1, diameters, thicknesses)
- Materials (E, rho), hydrodynamic added mass per length (mA)
- Lumped masses (M1 at waterline, M2 at top)
- Axial loads (P1, P2) derived from masses and gravity
- Base support (springs on/off, KL, KR)
- FEM mesh and geometric stiffness
- Analytical search range and resolution
- Output folders and water properties

Author: Aref Aasi

"""

from __future__ import annotations
import os
import numpy as np


# -------------------------------------------------------------------------
# Helper functions for safe typed input
# -------------------------------------------------------------------------
def _get_float(prompt: str) -> float:
    """Ask the user for a float value (required)."""
    while True:
        txt = input(f"{prompt}: ").strip()
        try:
            return float(txt)
        except ValueError:
            print("Please enter a valid numeric value.")


def _get_int(prompt: str) -> int:
    """Ask the user for an integer value (required)."""
    while True:
        txt = input(f"{prompt}: ").strip()
        try:
            return int(txt)
        except ValueError:
            print("Please enter a valid integer value.")


def _get_bool(prompt: str) -> bool:
    """Ask the user for a yes/no choice."""
    while True:
        txt = input(f"{prompt} [y/n]: ").strip().lower()
        if txt in ("y", "yes"):
            return True
        elif txt in ("n", "no"):
            return False
        print("Please type 'y' or 'n'.")


# -------------------------------------------------------------------------
# Geometry
# -------------------------------------------------------------------------
print("\n--- Geometry Parameters ---")
L   = _get_float("Total height of structure L (m)")
L1  = _get_float("Submerged length L1 (m)")
L2  = L - L1

D1  = _get_float("Outer diameter of submerged part D1 (m)")
D2  = _get_float("Outer diameter of tower part D2 (m)")
t1  = _get_float("Wall thickness of submerged part t1 (m)")
t2  = _get_float("Wall thickness of tower part t2 (m)")

# -------------------------------------------------------------------------
# Material Properties
# -------------------------------------------------------------------------
print("\n--- Material Properties ---")
E   = _get_float("Elastic modulus E (Pa)")
rho = _get_float("Steel density ρ (kg/m^3)")
mA  = _get_float("Hydrodynamic added mass per length mA (kg/m)")

# -------------------------------------------------------------------------
# Lumped Masses
# -------------------------------------------------------------------------
print("\n--- Lumped Masses ---")
M1  = _get_float("Waterline lumped mass M1 (kg)")
M2  = _get_float("Top mass (nacelle + rotor) M2 (kg)")
g   = 9.80665  # constant
P1  = (M1 + M2) * g
P2  = M2 * g

# -------------------------------------------------------------------------
# Base Springs
# -------------------------------------------------------------------------
print("\n--- Base Support Model ---")
USE_BASE_SPRINGS = _get_bool("Use base springs instead of clamped base?")
KL  = _get_float("Base translational spring stiffness KL (N/m)")
KR  = _get_float("Base rotational spring stiffness KR (N·m/rad)")

# -------------------------------------------------------------------------
# FEM Settings
# -------------------------------------------------------------------------
print("\n--- FEM Mesh Settings ---")
n1 = _get_int("Number of elements in submerged segment n1")
n2 = _get_int("Number of elements in tower segment n2")
USE_AXIAL_GEO = _get_bool("Include geometric stiffness (axial load)?")
N_axial_1 = P1
N_axial_2 = P2

# -------------------------------------------------------------------------
# Analytical Frequency Search
# -------------------------------------------------------------------------
print("\n--- Analytical Frequency Search Settings ---")
NMODES_PRINT = _get_int("Number of modes to print")
SCAN_FMIN = _get_float("Frequency scan min (Hz)")
SCAN_FMAX = _get_float("Frequency scan max (Hz)")
NSCAN = _get_int("Number of scan points (resolution)")

# -------------------------------------------------------------------------
# Output Paths
# -------------------------------------------------------------------------
FIG_DIR = "figs"
RES_DIR = "results"
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(RES_DIR, exist_ok=True)

# -------------------------------------------------------------------------
# Water Properties
# -------------------------------------------------------------------------
print("\n--- Water Properties ---")
RHO_WATER = _get_float("Water density ρ_water (kg/m^3)")
CA_BASELINE = mA / (RHO_WATER * (np.pi * D1**2) / 4.0)

# -------------------------------------------------------------------------
# Summary
# -------------------------------------------------------------------------
print("\nConfiguration complete!")
print(f"Geometry: L = {L} m, L1 = {L1} m, D1 = {D1} m, D2 = {D2} m")
print(f"Material: E = {E} Pa, ρ = {rho} kg/m³, mA = {mA} kg/m")
print(f"Masses: M1 = {M1} kg, M2 = {M2} kg")
print(f"Base springs: KL = {KL} N/m, KR = {KR} N·m/rad, enabled = {USE_BASE_SPRINGS}")
print(f"FEM mesh: n1 = {n1}, n2 = {n2}, geometric stiffness = {USE_AXIAL_GEO}")
print(f"Figures -> {FIG_DIR}/, Results -> {RES_DIR}/")
