"""
analytical.py — Two-segment EB analytical solution with axial load
------------------------------------------------------------------
Builds an 8×8 coefficient matrix enforcing:
- Base boundary: clamp OR translational/rotational springs (KL, KR)
- Interface at x=L1: continuity of y, y' and jumps/continuity in shear/moment
  (includes waterline lumped mass M1)
- Tip at x=L: free moment, shear with tip mass M2
The determinant (via signed log-det) is scanned over frequency to bracket roots,
then refined with Brent's method. Mode shapes are recovered from the null
vector (last right singular vector) of the 8×8 matrix at each root.


Author: Aref Aasi

"""


from __future__ import annotations
import numpy as np
from numpy import cos, sin, cosh, sinh
from scipy.optimize import brentq
from .config import (E, L, L1, L2, I1 as I1_cfg, I2 as I2_cfg, P1, P2, M1, M2,
                     SCAN_FMIN, SCAN_FMAX, NSCAN)
from .utils import pad_modes

# If you want I1/I2 from config.t thicknesses, set them there and import.
I1, I2 = I1_cfg, I2_cfg

def alpha_beta(EI: float, P: float, mline: float, omega: float):
    a = 1.0
    b = P / EI
    c = -(mline * omega**2) / EI
    disc = max(b*b - 4*a*c, 0.0)
    r2p  = (-b + np.sqrt(disc)) / 2.0
    r2m  = (-b - np.sqrt(disc)) / 2.0
    alpha = np.sqrt(max(r2p, 0.0))
    beta  = np.sqrt(max(-r2m, 0.0))
    return alpha, beta

def W_vec(x, a, b):   return np.array([cosh(a*x), sinh(a*x), cos(b*x), sin(b*x)])
def dW_vec(x, a, b):  return np.array([a*sinh(a*x), a*cosh(a*x), -b*sin(b*x), b*cos(b*x)])
def d2W_vec(x, a, b): return np.array([a*a*cosh(a*x), a*a*sinh(a*x), -b*b*cos(b*x), -b*b*sin(b*x)])
def d3W_vec(x, a, b): return np.array([a**3*sinh(a*x), a**3*cosh(a*x), b**3*sin(b*x), -b**3*cos(b*x)])

def coeff_matrix(omega: float, use_base_springs: bool, KL_val: float, KR_val: float,
                 m1_line_val: float, m2_line_val: float):
    a1, b1 = alpha_beta(E*I1, P1, m1_line_val, omega)
    a2, b2 = alpha_beta(E*I2, P2, m2_line_val, omega)
    M8 = np.zeros((8, 8))

    if not use_base_springs:
        M8[0, 0:4] = W_vec(0.0, a1, b1)
        M8[1, 0:4] = dW_vec(0.0, a1, b1)
    else:
        M8[0, 0:4] = (E*I1)*d3W_vec(0.0, a1, b1) - KL_val*W_vec(0.0, a1, b1)
        M8[1, 0:4] = (E*I1)*d2W_vec(0.0, a1, b1) + KR_val*dW_vec(0.0, a1, b1)

    # Tip conditions at x=L (with M2)
    M8[2, 4:8] = (E*I2)*d3W_vec(L, a2, b2) + (M2*(omega**2))*W_vec(L, a2, b2)
    M8[3, 4:8] = d2W_vec(L, a2, b2)

    # Interface at L1 (with M1 jump in shear)
    M8[4, 0:4] =  W_vec(L1, a1, b1)
    M8[4, 4:8] = -W_vec(L1, a2, b2)
    M8[5, 0:4] =  dW_vec(L1, a1, b1)
    M8[5, 4:8] = -dW_vec(L1, a2, b2)
    M8[6, 0:4] =  (E*I1)*d3W_vec(L1, a1, b1) + (M1*(omega**2))*W_vec(L1, a1, b1)
    M8[6, 4:8] = -(E*I2)*d3W_vec(L1, a2, b2)
    M8[7, 0:4] =  (E*I1)*d2W_vec(L1, a1, b1)
    M8[7, 4:8] = -(E*I2)*d2W_vec(L1, a2, b2)
    return M8

def signed_logdet(A: np.ndarray) -> float:
    sgn, logabs = np.linalg.slogdet(A)
    return 0.0 if sgn == 0 else sgn * logabs

def _detf(f_hz: float, use_base_springs: bool, KL_val: float, KR_val: float,
          m1_line_val: float, m2_line_val: float):
    omega = 2.0*np.pi*f_hz
    return signed_logdet(coeff_matrix(omega, use_base_springs, KL_val, KR_val,
                                      m1_line_val, m2_line_val))

def analytical_frequencies(use_base_springs: bool, KL_val: float, KR_val: float,
                           m1_line_val: float, m2_line_val: float,
                           fmin: float = SCAN_FMIN, fmax: float = SCAN_FMAX,
                           nscan: int = NSCAN, nroot: int = 6) -> np.ndarray:
    fs = np.linspace(fmin, fmax, nscan)
    vals = np.array([_detf(f, use_base_springs, KL_val, KR_val, m1_line_val, m2_line_val) for f in fs])
    roots = []
    for i in range(len(fs)-1):
        v1, v2 = vals[i], vals[i+1]
        if np.isnan(v1) or np.isnan(v2): continue
        if np.sign(v1) == 0:
            roots.append(fs[i])
        elif np.sign(v1) != np.sign(v2):
            try:
                r = brentq(lambda f: _detf(f, use_base_springs, KL_val, KR_val, m1_line_val, m2_line_val),
                           fs[i], fs[i+1], maxiter=200, xtol=1e-10)
                if all(abs(r - rr) > 1e-4 for rr in roots):
                    roots.append(r)
            except ValueError:
                pass
    roots.sort()
    return np.array(roots[:nroot])

def analytical_mode_shape(f_hz: float, use_base_springs: bool,
                          KL_val: float, KR_val: float, nx: int,
                          m1_line_val: float, m2_line_val: float):
    omega = 2*np.pi*f_hz
    M8 = coeff_matrix(omega, use_base_springs, KL_val, KR_val, m1_line_val, m2_line_val)
    U,S,Vt = np.linalg.svd(M8)
    c = Vt[-1,:]
    a1,b1 = alpha_beta(E*I1, P1, m1_line_val, omega)
    a2,b2 = alpha_beta(E*I2, P2, m2_line_val, omega)
    x1 = np.linspace(0, L1, nx//2)
    x2 = np.linspace(L1, L,  nx - nx//2)
    y1 = np.array([W_vec(x, a1,b1) @ c[0:4] for x in x1])
    y2 = np.array([W_vec(x, a2,b2) @ c[4:8] for x in x2])
    x  = np.concatenate([x1, x2]); y = np.concatenate([y1, y2])
    y /= np.max(np.abs(y))
    return x, y
