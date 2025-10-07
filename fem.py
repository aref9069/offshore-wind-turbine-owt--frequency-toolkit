"""
fem.py — 2-node EB FEM (two segments) with geometric stiffness
--------------------------------------------------------------
- ke_EB(EI, le): 4×4 beam bending stiffness
- me_EB(mline, le): 4×4 consistent mass
- kg_beam_column(N, le): 4×4 geometric stiffness from axial force N
- assemble_matrices(...): assembles K, M, (optionally Kg) across segments,
  adds base springs and lumped masses (M1 at waterline node, M2 at top),
  and returns plotting maps (physical x, w-DOF indices).
- fem_frequencies(...): symmetric GEP solve Kt φ = λ M φ, return first n modes.
- fem_mode_shape(...): extract and normalize lateral DOFs for plotting.
Notes:
- Units must match config (SI). Frequencies returned in Hz.
- If base springs are off, clamped base DOFs are removed (“static condensation”).


Author: Aref Aasi

"""


from __future__ import annotations
import numpy as np
from scipy.linalg import eigh
from .config import (E, L, L1, L2, I1 as I1_cfg, I2 as I2_cfg, M1, M2)
I1, I2 = I1_cfg, I2_cfg

def ke_EB(EI: float, le: float):
    l = le; c = EI/(l**3)
    return c*np.array([
        [ 12,    6*l,  -12,    6*l],
        [ 6*l,  4*l*l, -6*l,  2*l*l],
        [-12,   -6*l,   12,   -6*l],
        [ 6*l,  2*l*l, -6*l,  4*l*l]
    ], dtype=float)

def me_EB(mline: float, le: float):
    l = le; c = mline*l/420.0
    return c*np.array([
        [156,     22*l,   54,   -13*l],
        [22*l,   4*l*l,  13*l,  -3*l*l],
        [54,      13*l,  156,   -22*l],
        [-13*l, -3*l*l, -22*l,   4*l*l]
    ], dtype=float)

def kg_beam_column(N: float, le: float):
    l = le; c = N/(30.0*l)
    return c*np.array([
        [ 36,   3*l,  -36,   3*l],
        [ 3*l,  4*l*l, -3*l, -1*l*l],
        [-36,  -3*l,   36,  -3*l],
        [ 3*l, -1*l*l, -3*l,  4*l*l]
    ], dtype=float)

def assemble_matrices(n1_: int, n2_: int, use_springs: bool, KL_val: float, KR_val: float,
                      use_axial_geo: bool, N1: float, N2: float,
                      m1_line_val: float, m2_line_val: float):
    nnode = (n1_ + n2_) + 1
    ndof  = 2*nnode
    K  = np.zeros((ndof, ndof))
    M  = np.zeros((ndof, ndof))
    Kg = np.zeros((ndof, ndof)) if (use_axial_geo and (N1 != 0 or N2 != 0)) else None

    def add(A, Ae, dofs):
        for i, di in enumerate(dofs):
            A[di, dofs] += Ae[i, :]

    def dofpair(n): return [2*n, 2*n+1]

    # seg1
    le1 = L1/n1_
    for e in range(n1_):
        nL = e; nR = e+1
        dofs = dofpair(nL) + dofpair(nR)
        add(K, ke_EB(E*I1, le1), dofs)
        add(M, me_EB(m1_line_val, le1), dofs)
        if Kg is not None: add(Kg, kg_beam_column(N1, le1), dofs)

    # seg2
    off = n1_; le2 = L2/n2_
    for e in range(n2_):
        nL = off + e; nR = off + e + 1
        dofs = dofpair(nL) + dofpair(nR)
        add(K, ke_EB(E*I2, le2), dofs)
        add(M, me_EB(m2_line_val, le2), dofs)
        if Kg is not None: add(Kg, kg_beam_column(N2, le2), dofs)

    # base springs or clamp-eliminate
    base_w, base_th = 0, 1
    if use_springs:
        K[base_w,  base_w]  += KL_val
        K[base_th, base_th] += KR_val

    # lumped masses
    node_wl  = n1_
    node_top = nnode - 1
    M[2*node_wl,  2*node_wl]  += M1
    M[2*node_top, 2*node_top] += M2

    Kt = K + (Kg if Kg is not None else 0.0)

    if not use_springs:
        keep = np.arange(ndof, dtype=int)
        keep = keep[(keep != base_w) & (keep != base_th)]
        Kt = Kt[np.ix_(keep, keep)]
        M  = M[np.ix_(keep, keep)]
        map_w_dofs = np.arange(0, 2*nnode, 2)[1:]  # all w except base
        map_x = np.linspace(0, L, nnode)[1:]       # no base node
    else:
        map_w_dofs = np.arange(0, 2*nnode, 2)      # all w DOFs
        map_x = np.linspace(0, L, nnode)

    return Kt, M, map_x, map_w_dofs

def fem_frequencies(n1_: int, n2_: int, use_springs: bool, KL_val: float, KR_val: float,
                    use_axial_geo: bool, N1: float, N2: float, nmode: int,
                    m1_line_val: float, m2_line_val: float):
    Kt, M, map_x, map_w = assemble_matrices(n1_, n2_, use_springs, KL_val, KR_val,
                                            use_axial_geo, N1, N2, m1_line_val, m2_line_val)
    Kt = 0.5*(Kt + Kt.T)
    M  = 0.5*(M  + M.T)
    evals, evecs = eigh(Kt, M)
    lam = np.real(evals[evals > 1e-12]); lam.sort()
    omegas = np.sqrt(lam); freqs = omegas/(2.0*np.pi)
    return freqs[:nmode], (Kt, M, evecs, map_x, map_w)

def fem_mode_shape(mode_index: int, evecs: np.ndarray, map_x, map_w_dofs):
    phi = np.real(evecs[map_w_dofs, mode_index])
    phi /= np.max(np.abs(phi))
    return map_x, phi
