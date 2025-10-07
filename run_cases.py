"""
run_cases.py — Reproducible example runs
----------------------------------------
- Executes two scenarios:
  A) Clamped base
  B) Coupled-spring base (KL, KR) from config
- Saves:
  results/freq_compare_<case>.csv
  figs/freq_compare_<case>.png|pdf
  figs/mode_shape_overlay_<case>_mode1.png|pdf
Edit config.py to change geometry/materials, spring values, mesh, and search ranges.



Author: Aref Aasi

"""


from __future__ import annotations
import os
import numpy as np
import matplotlib.pyplot as plt
from .config import *
from .utils import save_png_pdf, pad_modes, tube_A_I_exact
from .analytical import analytical_frequencies, analytical_mode_shape
from .fem import fem_frequencies, fem_mode_shape

# Precompute section properties & line masses
A1, I1 = tube_A_I_exact(D1, t1)
A2, I2 = tube_A_I_exact(D2, t2)
m1_line = rho*A1 + mA
m2_line = rho*A2

def compare_table(f_anal, f_fem, n=6, save_csv_path=None, save_fig_path=None, title=None):
    nshow = min(n, len(f_anal), len(f_fem))
    lines = ["Mode,Analytical(Hz),FEM(Hz),Delta%(Anal-FEM)\n"]
    rows = []
    for i in range(nshow):
        fa, ff = f_anal[i], f_fem[i]
        d = 100.0*(fa-ff)/ff
        lines.append(f"{i+1},{fa:.6f},{ff:.6f},{d:.2f}\n")
        rows.append([i+1, f"{fa:.6f}", f"{ff:.6f}", f"{d:.2f}%"])
    if save_csv_path:
        with open(save_csv_path, "w", encoding="utf-8") as f:
            f.writelines(lines)
    if save_fig_path:
        fig, ax = plt.subplots(figsize=(6.5, 1.0 + 0.4*nshow)); ax.axis("off")
        table = ax.table(cellText=rows,
                         colLabels=["Mode","Analytical (Hz)","FEM (Hz)","Δ% (Anal–FEM)"],
                         cellLoc="center", loc="center")
        table.auto_set_font_size(False); table.set_fontsize(9); table.scale(1.1, 1.2)
        if title: ax.set_title(title, pad=10)
        save_png_pdf(save_fig_path); plt.close(fig)

def run_case(case_name: str, use_springs: bool, KL_val: float, KR_val: float):
    # Analytical
    f_anal = analytical_frequencies(
        use_base_springs=use_springs, KL_val=KL_val, KR_val=KR_val,
        m1_line_val=m1_line, m2_line_val=m2_line,
        fmin=SCAN_FMIN, fmax=SCAN_FMAX, nscan=NSCAN, nroot=NMODES_PRINT
    )
    # FEM
    f_fem, (_Kt,_M,evecs,map_x,map_w) = fem_frequencies(
        n1, n2, use_springs, KL_val, KR_val,
        USE_AXIAL_GEO, N_axial_1, N_axial_2, NMODES_PRINT,
        m1_line, m2_line
    )
    # Save side-by-side table
    compare_table(
        f_anal, f_fem, n=NMODES_PRINT,
        save_csv_path=os.path.join(RES_DIR, f"freq_compare_{case_name}.csv"),
        save_fig_path=os.path.join(FIG_DIR, f"freq_compare_{case_name}.png"),
        title=f"Natural Frequencies: Analytical vs FEM ({case_name})"
    )
    # Example mode shape overlay for mode 1
    xA, yA = analytical_mode_shape(f_anal[0], use_springs, KL_val, KR_val, nx=400,
                                   m1_line_val=m1_line, m2_line_val=m2_line)
    xF, yF = fem_mode_shape(0, evecs, map_x, map_w)
    plt.figure(); plt.plot(xA, yA, label="Analytical", linewidth=2)
    plt.plot(xF, yF, label="FEM", linestyle="--")
    plt.xlabel("Height x (m)"); plt.ylabel("Normalized lateral deflection")
    plt.title(f"Mode shape overlay (mode 1) – {case_name}"); plt.legend(); plt.tight_layout()
    save_png_pdf(os.path.join(FIG_DIR, f"mode_shape_overlay_{case_name}_mode1.png")); plt.close()
    return f_anal, f_fem

def main():
    print(f"[Info] Geometry: L={L} m (L1={L1} m), D1={D1} m, D2={D2} m, t1={t1} m, t2={t2} m")
    print(f"[Info] Masses: M1={M1} kg, M2={M2} kg; added mass mA={mA} kg/m")
    print(f"[Info] FEM mesh: n1={n1}, n2={n2}; axial_geo={USE_AXIAL_GEO}")
    print(f"[Info] E={E:.3e} Pa, rho={rho} kg/m^3")

    print("\n=== CASE A: CLAMPED BASE ===")
    run_case("clamped", use_springs=False, KL_val=0.0, KR_val=0.0)

    print("\n=== CASE B: COUPLED-SPRING BASE ===")
    run_case("cs", use_springs=True, KL_val=KL, KR_val=KR)

if __name__ == "__main__":
    main()
