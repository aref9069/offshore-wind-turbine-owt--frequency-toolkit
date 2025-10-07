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


def _get_float(prompt: str, default: float) -> float:
    """Helper: ask for a float input with default fallback."""
    txt = input(f"{prompt} [{default}]: ")
    return float(txt) if txt.strip() else default


def _get_int(prompt: str, default: int) -> int:
    """Helper: ask for an integer input with default fallback."""
    txt = input(f"{prompt} [{default}]: ")
    return int(txt) if txt.strip() else default


def _get_bool(prompt: str, default: bool) -> bool:
    """Helper: ask for yes/no input with default fallback."""
    txt = input(f"{prompt} [y/n, default={'y' if default else 'n'}]: ").strip().lower()
    if not txt:
        return default
    return txt in ("y", "yes", "true", "1")


# ---------------- Geometry ----------------
print("\n--- Geometry Parameters ---")
L   = _get_float("Total height of structure L (m)", 143.6)
L1  = _get_float("Submerged length L1 (m)", 56.0)
L2  = L - L1

D1  = _get_float("Outer diameter of submerged part D1 (m)", 6.0)
D2  = _get_float("Outer diameter of tower part D2 (m)", 3.87)
t1  = _get_float("Wall thickness of submerged part t1 (m)", 0.060)
t2  = _get_float("Wall thickness of tower part t2 (m)", 0.035)

# ---------------- Material Properties ----------------
print("\n--- Material Properties ---")
E   = _get_float("Elastic modulus E (Pa)", 210e9)
rho = _get_float("Steel density ρ (kg/m^3)", 8500.0)
mA  = _get_float("Hydrodynamic added mass per length mA (kg/m)", 9837.2)

# ---------------- Lumped Masses ----------------
print("\n--- Lumped Masses ---")
M1  = _get_float("Waterline lumped mass M1 (kg)", 1.0e4)
M2  = _get_float("Top mass (nacelle+rotor) M2 (kg)", 3.5e5)
g   = 9.80665  # constant
P1  = (M1 + M2) * g
P2  = M2 * g

# ---------------- Base Springs ----------------
print("\n--- Base Support Model ---")
USE_BASE_SPRINGS = _get_bool("Use base springs instead of clamped base?", True)
KL  = _get_float("Base translational spring stiffness KL (N/m)", 1e12)
KR  = _get_float("Base rotational spring stiffness KR (N·m/rad)", 1e14)

# ---------------- FEM Settings ----------------
print("\n--- FEM Mesh Settings ---")
n1 = _get_int("Number of elements in submerged segment n1", 30)
n2 = _get_int("Number of elements in tower segment n2", 50)
USE_AXIAL_GEO = _get_bool("Include geometric stiffness (axial load)?", True)
N_axial_1 = P1
N_axial_2 = P2

# ---------------- Analytical Frequency Search ----------------
print("\n--- Analytical Frequency Search Settings ---")
NMODES_PRINT = _get_int("Number of modes to print", 6)
SCAN_FMIN = _get_float("Frequency scan min (Hz)", 0.05)
SCAN_FMAX = _get_float("Frequency scan max (Hz)", 5.0)
NSCAN = _get_int("Number of scan points (resolution)", 2400)

# ---------------- Output Paths ----------------
FIG_DIR = "figs"
RES_DIR = "results"
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(RES_DIR, exist_ok=True)

# ---------------- Water Properties ----------------
RHO_WATER = _get_float("Water density ρ_water (kg/m^3)", 1025.0)
CA_BASELINE = mA / (RHO_WATER * (np.pi * D1**2) / 4.0)

print("\nConfiguration complete!")
print(f"L = {L:.2f} m, L1 = {L1:.2f} m, D1 = {D1:.2f} m, D2 = {D2:.2f} m")
print(f"E = {E:.2e} Pa, ρ = {rho:.1f} kg/m³, mA = {mA:.1f} kg/m")
print(f"M1 = {M1:.1f} kg, M2 = {M2:.1f} kg")
print(f"Base springs: KL={KL:.2e}, KR={KR:.2e}, use={USE_BASE_SPRINGS}")
print(f"FEM mesh: n1={n1}, n2={n2}, geometric stiffness={USE_AXIAL_GEO}")
print(f"Figures -> {FIG_DIR}/, Results -> {RES_DIR}/")
