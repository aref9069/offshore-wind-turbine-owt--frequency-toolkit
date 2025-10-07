"""
utils.py — Shared helpers
-------------------------
- save_png_pdf(path): save both PNG and vector PDF (journal-friendly)
- pad_modes(arr, m): pad/trim frequency arrays to exactly m modes
- tube_A_I_exact(D, t): exact area & second moment for thin-walled tube
- temp_depth_and_added_mass(...): context manager to temporarily override
  L1/L2 and submerged line mass (rho*A1 + mA_tmp) during sweeps.
- mA_from_Ca(Ca): convert added-mass coefficient to kg/m using D1, rho_water



Author: Aref Aasi

"""


from __future__ import annotations
import os
import numpy as np
import matplotlib.pyplot as plt
from contextlib import contextmanager
from .config import (L, L1 as L1_cfg, D1, E, rho, mA as mA_cfg, RHO_WATER)

def save_png_pdf(path_png: str):
    plt.savefig(path_png, dpi=200, bbox_inches="tight")
    root, _ = os.path.splitext(path_png)
    plt.savefig(root + ".pdf", bbox_inches="tight")

def pad_modes(freqs: np.ndarray, modes: int) -> np.ndarray:
    out = np.full(modes, np.nan, dtype=float)
    m = min(modes, len(freqs))
    if m > 0:
        out[:m] = freqs[:m]
    return out

def tube_A_I_exact(D: float, t: float):
    Di = D - 2.0*t
    if Di <= 0:
        raise ValueError("Wall thickness too large; inner diameter <= 0.")
    A = 0.25*np.pi*(D**2 - Di**2)
    I = (np.pi/64.0)*(D**4 - Di**4)
    return A, I

# Globals held here so both modules can share consistent overrides
A1, I1 = tube_A_I_exact(D1, 0.060)  # default thickness; may be recomputed by callers
m1_line = rho*A1 + mA_cfg
m2_line = None  # set by caller

@contextmanager
def temp_depth_and_added_mass(L1_new: float, mA_new: float, A1_local: float):
    """
    Temporarily override L1, L2, and m1_line (per-length mass submerged).
    """
    from . import config
    global m1_line
    L1_old, m1_old = config.L1, m1_line
    config.L1 = float(L1_new)
    config.L2 = config.L - config.L1
    m1_line = rho * A1_local + mA_new
    try:
        yield
    finally:
        config.L1 = L1_old
        config.L2 = L - L1_old
        m1_line = m1_old

def mA_from_Ca(Ca_val: float) -> float:
    return Ca_val * RHO_WATER * (np.pi * D1**2) / 4.0
