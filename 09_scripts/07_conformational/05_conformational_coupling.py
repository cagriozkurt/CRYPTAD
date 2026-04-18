"""
CRYPTAD — W16 Apo conformational analysis on metadynamics trajectories.

Three analyses:

  1. DCCM (Dynamic Cross-Correlation Matrix)
     Both S1_BIN1_BAR and S3_PICALM_ANTH.
     Cα cross-correlations after backbone alignment.
     Highlights correlated/anti-correlated domain motions around cryptic pockets.

  2. BAR Curvature Analysis (S1 BIN1)
     Tracks whether site688/site680 pocket opening (CV1 = cross-chain Cα distance)
     correlates with BAR crescent geometry change.
     Metric: inter-tip distance = Cα centroid of N-terminal helix chain A (resids 1–15)
             vs Cα centroid of N-terminal helix chain B (resids 192–206).
     A longer inter-tip distance = straighter (less curved) BAR domain.

  3. PIP2-Loop Analysis (S3 PICALM)
     Tracks whether site473 pocket opening (cv_site4 = ASP203–SER248 distance)
     correlates with displacement of the PIP2-binding helix.
     Metric: RMSD of PIP2-binding region (resids 1–22, incl. KKK basic patch at 20–22)
             after aligning to the stable ANTH core (resids 50–200).

Outputs (03_pocket_analysis/conformational/):
  dccm_S1_BIN1_BAR.png          — DCCM heatmap for BIN1 BAR (apo)
  dccm_S3_PICALM_ANTH.png       — DCCM heatmap for PICALM ANTH (apo)
  dccm_S3_PICALM_ANTH_holo.png  — DCCM heatmap for PICALM ANTH (Ibufenac holo, rep1+rep5)
  dccm_delta_S3_PICALM.png      — 3-panel ΔDCCM: apo | holo | holo−apo
  bar_curvature.png             — CV1 vs inter-tip distance (BIN1)
  pip2_loop.png                 — cv_site4 vs PIP2-helix RMSD (PICALM)
  conformational_summary.txt    — statistics and Pearson r values
  conformational_manifest.json  — run parameters and output paths

Usage:
  python3 09_scripts/07_conformational/05_conformational_coupling.py
  python3 09_scripts/07_conformational/05_conformational_coupling.py --project-root /path/to/CRYPTAD
"""

import argparse
import json
import logging
import time
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
logging.getLogger("MDAnalysis").setLevel(logging.ERROR)

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
import MDAnalysis as mda
from MDAnalysis.analysis import align

_script       = Path(__file__).resolve()
_default_root = _script.parents[2]

matplotlib.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 10,
    "axes.linewidth": 0.8,
    "pdf.fonttype": 42,
    "figure.dpi": 150,
})


# ── Analysis functions ─────────────────────────────────────────────────────────

def compute_dccm(tpr, xtc, sel="protein and name CA", stride=1):
    """
    Compute Cα dynamic cross-correlation matrix after backbone alignment.
    Returns (dccm, resids, resnames, segids) arrays.
    """
    u = mda.Universe(str(tpr), str(xtc))
    ca = u.select_atoms(sel)

    # Align all frames to first frame (backbone Cα)
    align.AlignTraj(u, u, select=sel, in_memory=True).run()

    # Collect positions
    positions = []
    for i, ts in enumerate(u.trajectory):
        if i % stride == 0:
            positions.append(ca.positions.copy())
    positions = np.array(positions)    # (n_frames, n_atoms, 3)
    n_frames  = positions.shape[0]

    # Mean positions and displacements
    mean_pos  = positions.mean(axis=0)           # (n_atoms, 3)
    disp      = positions - mean_pos[np.newaxis]  # (n_frames, n_atoms, 3)

    # Cross-correlation: C(i,j) = <Δri·Δrj> / sqrt(<|Δri|²><|Δrj|²>)
    dot  = np.einsum('fid,fjd->ij', disp, disp) / n_frames   # (n, n)
    norm = np.sqrt(np.diag(dot))                               # (n,)
    dccm = dot / np.outer(norm, norm)
    np.clip(dccm, -1.0, 1.0, out=dccm)

    resids   = np.array([a.resid   for a in ca])
    resnames = np.array([a.resname for a in ca])
    segids   = np.array([a.segid   for a in ca])
    return dccm, resids, resnames, segids


def plot_dccm(dccm, resids, segids, title, outpng, pocket_resids=None):
    """
    Plot DCCM heatmap with residue axis, chain separators, and pocket markers.
    """
    n = len(resids)
    fig, ax = plt.subplots(figsize=(7, 6))
    im = ax.imshow(dccm, cmap="RdBu_r", vmin=-1, vmax=1,
                   origin="lower", aspect="auto", interpolation="nearest")
    plt.colorbar(im, ax=ax, label="Cross-correlation", fraction=0.046, pad=0.04)

    # Chain separator
    unique_segs = list(dict.fromkeys(segids))
    for seg in unique_segs[1:]:
        boundary = np.where(segids == seg)[0][0] - 0.5
        ax.axhline(boundary, color="white", lw=1.0, ls="--", alpha=0.7)
        ax.axvline(boundary, color="white", lw=1.0, ls="--", alpha=0.7)

    # Pocket residue markers
    if pocket_resids:
        for resid in pocket_resids:
            idx_arr = np.where(resids == resid)[0]
            if len(idx_arr):
                idx = idx_arr[0]
                ax.axhline(idx, color="#FFDD00", lw=0.7, alpha=0.6)
                ax.axvline(idx, color="#FFDD00", lw=0.7, alpha=0.6)

    # Axis labels — tick every ~12 positions
    tick_idx = np.arange(0, n, max(1, n // 12))
    ax.set_xticks(tick_idx)
    ax.set_xticklabels([str(resids[i]) for i in tick_idx], fontsize=7, rotation=45)
    ax.set_yticks(tick_idx)
    ax.set_yticklabels([str(resids[i]) for i in tick_idx], fontsize=7)
    ax.set_xlabel("Residue index", fontsize=10)
    ax.set_ylabel("Residue index", fontsize=10)
    ax.set_title(title, fontsize=11, fontweight="bold")

    # Chain labels — inside the top and left edges to avoid collision
    for seg in unique_segs:
        mask   = segids == seg
        center = np.where(mask)[0].mean()
        label  = seg.split("_")[-1]
        ax.text(center, n - 4, label, ha="center", va="top",
                fontsize=7, color="#555555",
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.7, pad=1.5))
        ax.text(4, center, label, ha="left", va="center",
                fontsize=7, color="#555555",
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.7, pad=1.5))

    plt.tight_layout()
    plt.savefig(str(outpng), dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {outpng}")


# ══════════════════════════════════════════════════════════════════════════════
def main():
    parser = argparse.ArgumentParser(
        description="CRYPTAD: W16 apo conformational coupling analysis")
    parser.add_argument("--project-root", type=Path, default=_default_root,
                        help=f"Path to CRYPTAD project root (default: {_default_root})")
    args = parser.parse_args()

    project_root = args.project_root.resolve()
    outdir = project_root / "03_pocket_analysis" / "conformational"
    outdir.mkdir(parents=True, exist_ok=True)

    S1_TPR = project_root / "02_md_simulations/BIN1/metadynamics/S1_BAR/meta.tpr"
    S1_XTC = project_root / "02_md_simulations/BIN1/metadynamics/S1_BAR/meta_nojump.xtc"
    S1_COL = project_root / "02_md_simulations/BIN1/metadynamics/S1_BAR/COLVAR"

    S3_TPR = project_root / "02_md_simulations/PICALM/metadynamics/S3_ANTH/meta.tpr"
    S3_XTC = project_root / "02_md_simulations/PICALM/metadynamics/S3_ANTH/meta.xtc"
    S3_COL = project_root / "02_md_simulations/PICALM/metadynamics/S3_ANTH/COLVAR"

    summary_lines = ["CRYPTAD — W16 Conformational Analysis Summary", "=" * 60, ""]
    manifest_data = {}

    # ══════════════════════════════════════════════════════════════════════════
    # 1. DCCM — S1 BIN1 BAR
    # ══════════════════════════════════════════════════════════════════════════
    print("\n" + "="*60)
    print("1. DCCM — S1 BIN1 BAR")
    print("="*60)
    CACHE_S1 = outdir / "dccm_s1_cache.npz"
    if CACHE_S1.exists():
        print("  Loading cached DCCM (S1) …")
        _c = np.load(str(CACHE_S1), allow_pickle=True)
        dccm_s1 = _c["dccm"]; res_s1 = _c["res"]; rn_s1 = _c["rn"]; seg_s1 = _c["seg"]
    else:
        dccm_s1, res_s1, rn_s1, seg_s1 = compute_dccm(S1_TPR, S1_XTC)
        np.savez(str(CACHE_S1), dccm=dccm_s1, res=res_s1, rn=rn_s1, seg=seg_s1)
        print("  DCCM (S1) cached.")
    print(f"  {len(res_s1)} Cα atoms, {dccm_s1.shape} matrix")
    plot_dccm(dccm_s1, res_s1, seg_s1,
              "DCCM — BIN1 BAR (S1, 200 ns metadynamics)",
              outdir / "dccm_S1_BIN1_BAR.png",
              pocket_resids=None)

    # Summary stats: mean |correlation| within chain vs between chains
    mask_A   = seg_s1 == "seg_0_PROA"
    mask_B   = seg_s1 == "seg_1_PROB"
    intra_A  = dccm_s1[np.ix_(mask_A, mask_A)]
    intra_B  = dccm_s1[np.ix_(mask_B, mask_B)]
    inter_AB = dccm_s1[np.ix_(mask_A, mask_B)]
    np.fill_diagonal(intra_A, np.nan)
    np.fill_diagonal(intra_B, np.nan)
    summary_lines += [
        "BIN1 BAR DCCM:",
        f"  Mean |C| intra-chain A:  {np.nanmean(np.abs(intra_A)):.3f}",
        f"  Mean |C| intra-chain B:  {np.nanmean(np.abs(intra_B)):.3f}",
        f"  Mean |C| inter-chain AB: {np.nanmean(np.abs(inter_AB)):.3f}",
        f"  Anti-corr fraction (inter, C<-0.3): "
        f"{(inter_AB < -0.3).mean():.3f}",
        "",
    ]

    # ══════════════════════════════════════════════════════════════════════════
    # 2. DCCM — S3 PICALM ANTH
    # ══════════════════════════════════════════════════════════════════════════
    print("\n" + "="*60)
    print("2. DCCM — S3 PICALM ANTH")
    print("="*60)
    CACHE_S3 = outdir / "dccm_s3_cache.npz"
    if CACHE_S3.exists():
        print("  Loading cached DCCM (S3) …")
        _c = np.load(str(CACHE_S3), allow_pickle=True)
        dccm_s3 = _c["dccm"]; res_s3 = _c["res"]; rn_s3 = _c["rn"]; seg_s3 = _c["seg"]
    else:
        dccm_s3, res_s3, rn_s3, seg_s3 = compute_dccm(S3_TPR, S3_XTC, stride=10)
        np.savez(str(CACHE_S3), dccm=dccm_s3, res=res_s3, rn=rn_s3, seg=seg_s3)
        print("  DCCM (S3) cached.")
    print(f"  {len(res_s3)} Cα atoms, {dccm_s3.shape} matrix")
    plot_dccm(dccm_s3, res_s3, seg_s3,
              "DCCM — PICALM ANTH (S3, 200 ns metadynamics)",
              outdir / "dccm_S3_PICALM_ANTH.png",
              pocket_resids=[194, 219, 235, 236, 246])

    # Split DCCM into PIP2-loop (1-22) vs pocket (190-250) vs core (50-180)
    pip2_mask   = (res_s3 >= 1)   & (res_s3 <= 22)
    pocket_mask = (res_s3 >= 190) & (res_s3 <= 260)
    core_mask   = (res_s3 >= 50)  & (res_s3 <= 180)
    pip2_pocket = dccm_s3[np.ix_(pip2_mask, pocket_mask)]
    pip2_core   = dccm_s3[np.ix_(pip2_mask, core_mask)]
    summary_lines += [
        "PICALM ANTH DCCM:",
        f"  Mean C (PIP2-loop ↔ site473 pocket): {pip2_pocket.mean():.3f}",
        f"  Mean C (PIP2-loop ↔ core):           {pip2_core.mean():.3f}",
        f"  Anti-corr fraction (PIP2 ↔ pocket, C<-0.2): "
        f"{(pip2_pocket < -0.2).mean():.3f}",
        "",
    ]

    # ══════════════════════════════════════════════════════════════════════════
    # 3. BAR Curvature Analysis (S1 BIN1)
    # ══════════════════════════════════════════════════════════════════════════
    print("\n" + "="*60)
    print("3. BAR Curvature — S1 BIN1")
    print("="*60)

    u_s1  = mda.Universe(str(S1_TPR), str(S1_XTC))
    # N-terminal helix of each chain (first 15 residues per chain = crescent tips)
    tip_A = u_s1.select_atoms("protein and name CA and segid seg_0_PROA and resid 1:15")
    tip_B = u_s1.select_atoms("protein and name CA and segid seg_1_PROB and resid 192:206")
    print(f"  Tip A: {len(tip_A)} atoms (resids {tip_A.resids[0]}–{tip_A.resids[-1]})")
    print(f"  Tip B: {len(tip_B)} atoms (resids {tip_B.resids[0]}–{tip_B.resids[-1]})")

    # Load COLVAR — every 250th row matches each 500 ps trajectory frame
    col_s1 = np.loadtxt(str(S1_COL), comments="#")
    # col_s1[:,0] = time (ps), col_s1[:,1] = cv1 (nm), col_s1[:,2] = cv2 (rad)
    n_traj = len(u_s1.trajectory)
    cv1_traj = col_s1[::250, 1][:n_traj]  # nm

    tip_tip_dist = []
    for ts in u_s1.trajectory:
        centA = tip_A.positions.mean(axis=0)
        centB = tip_B.positions.mean(axis=0)
        tip_tip_dist.append(np.linalg.norm(centA - centB) / 10.0)  # Å → nm

    tip_tip_dist = np.array(tip_tip_dist)
    time_ns = col_s1[::250, 0][:n_traj] / 1000.0  # ps → ns

    # Correlation
    r, p = stats.pearsonr(cv1_traj, tip_tip_dist)
    print(f"  Pearson r(CV1, tip-tip dist) = {r:.3f}  p = {p:.2e}")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    # Panel A: timeseries
    ax = axes[0]
    ax2 = ax.twinx()
    ln1 = ax.plot(time_ns, cv1_traj, color="#1A6FA8", lw=0.8, alpha=0.8, label="CV1 (cross-chain)")
    ln2 = ax2.plot(time_ns, tip_tip_dist, color="#D95F3B", lw=0.8, alpha=0.8, label="Tip-tip distance")
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("CV1 (nm)", color="#1A6FA8")
    ax2.set_ylabel("N-term tip-to-tip distance (nm)", color="#D95F3B")
    ax.tick_params(axis='y', labelcolor="#1A6FA8")
    ax2.tick_params(axis='y', labelcolor="#D95F3B")
    lns = ln1 + ln2
    ax.legend(lns, [l.get_label() for l in lns], fontsize=8, loc="upper left")
    ax.set_title("A   BAR pocket opening vs crescent geometry", fontweight="bold", fontsize=10)
    for spine in ["top"]:
        ax.spines[spine].set_visible(False)
        ax2.spines[spine].set_visible(False)

    # Panel B: scatter
    ax = axes[1]
    sc = ax.scatter(cv1_traj, tip_tip_dist, c=time_ns, cmap="viridis",
                    s=10, alpha=0.6, linewidths=0)
    plt.colorbar(sc, ax=ax, label="Time (ns)", fraction=0.046, pad=0.04)
    m, b = np.polyfit(cv1_traj, tip_tip_dist, 1)
    xfit = np.linspace(cv1_traj.min(), cv1_traj.max(), 100)
    ax.plot(xfit, m * xfit + b, "r--", lw=1.5, alpha=0.8)
    ax.set_xlabel("CV1 — cross-chain distance (nm)")
    ax.set_ylabel("N-term tip-to-tip distance (nm)")
    ax.set_title(f"B   Pearson r = {r:.3f}  (p = {p:.1e})", fontweight="bold", fontsize=10)
    ax.text(0.97, 0.05, f"r = {r:.3f}", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=11, fontweight="bold", color="firebrick")
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

    plt.suptitle("BIN1 BAR — Site688/680 Pocket Opening vs Crescent Curvature\n"
                 "200 ns metadynamics", fontsize=11, y=1.01)
    plt.tight_layout()
    plt.savefig(str(outdir / "bar_curvature.png"), dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {outdir}/bar_curvature.png")

    summary_lines += [
        "BIN1 BAR Curvature:",
        f"  CV1 range: {cv1_traj.min():.2f}–{cv1_traj.max():.2f} nm (mean {cv1_traj.mean():.2f})",
        f"  Tip-tip distance: {tip_tip_dist.min():.2f}–{tip_tip_dist.max():.2f} nm (mean {tip_tip_dist.mean():.2f})",
        f"  Pearson r(CV1, tip-tip) = {r:.3f}  p = {p:.2e}",
        "  Interpretation: " + (
            "pocket opening (↑CV1) correlates with crescent straightening (↑tip-tip)"
            if r > 0.3 else
            "pocket opening and BAR curvature are independent in this simulation"
            if abs(r) < 0.2 else
            "pocket opening (↑CV1) correlates with crescent tightening (↓tip-tip)"
        ),
        "",
    ]
    manifest_data["bar_curvature"] = {
        "pearson_r": float(round(r, 4)),
        "pearson_p": float(f"{p:.2e}"),
        "cv1_mean_nm": float(round(cv1_traj.mean(), 3)),
        "tip_tip_mean_nm": float(round(tip_tip_dist.mean(), 3)),
    }

    # ══════════════════════════════════════════════════════════════════════════
    # 4. PIP2-Loop Analysis (S3 PICALM)
    # ══════════════════════════════════════════════════════════════════════════
    print("\n" + "="*60)
    print("4. PIP2-Loop — S3 PICALM")
    print("="*60)

    u_s3 = mda.Universe(str(S3_TPR), str(S3_XTC))
    col_s3 = np.loadtxt(str(S3_COL), comments="#")
    n_traj_s3 = len(u_s3.trajectory)

    # cv_site4 = col_s3[:,2], every 10th row matches each 20 ps frame
    cv_site4   = col_s3[::10, 2][:n_traj_s3]    # nm
    time_ns_s3 = col_s3[::10, 0][:n_traj_s3] / 1000.0  # ps → ns

    # PIP2-binding helix (resids 1-22, KKK basic patch)
    pip2_helix = u_s3.select_atoms("protein and name CA and resid 1:22")
    # Stable ANTH core for alignment reference
    anth_core  = u_s3.select_atoms("protein and name CA and resid 50:200")
    print(f"  PIP2-helix:  {len(pip2_helix)} Cα (resids 1–22)")
    print(f"  ANTH core:   {len(anth_core)} Cα (resids 50–200) [alignment reference]")

    # Capture reference positions from frame 0 (Kabsch RMSD vs frame 0)
    u_s3.trajectory[0]
    ref_pip2 = pip2_helix.positions.copy()
    ref_core = anth_core.positions.copy()

    # PIP2-helix RMSD after aligning each frame's core to reference core (Å)
    pip2_rmsd_A = []
    for ts in u_s3.trajectory:
        core_now = anth_core.positions.copy()
        pip2_now = pip2_helix.positions.copy()
        c_ref = ref_core.mean(axis=0);  c_now = core_now.mean(axis=0)
        A = (core_now - c_now).T @ (ref_core - c_ref)
        U, S, Vt = np.linalg.svd(A)
        d = np.linalg.det(Vt.T @ U.T)
        R = Vt.T @ np.diag([1, 1, d]) @ U.T
        pip2_rot = (pip2_now - c_now) @ R.T + c_ref
        pip2_rmsd_A.append(np.sqrt(np.mean(np.sum((pip2_rot - ref_pip2)**2, axis=1))))
    pip2_rmsd_A = np.array(pip2_rmsd_A)   # Å (MDAnalysis positions are in Å)

    r3, p3 = stats.pearsonr(cv_site4, pip2_rmsd_A)
    print(f"  Pearson r(cv_site4, PIP2-helix RMSD) = {r3:.3f}  p = {p3:.2e}")
    print(f"  PIP2-helix RMSD: {pip2_rmsd_A.min():.2f}–{pip2_rmsd_A.max():.2f} Å  (mean {pip2_rmsd_A.mean():.2f})")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    # Panel A: timeseries
    ax = axes[0]
    ax2 = ax.twinx()
    step = max(1, len(time_ns_s3) // 2000)
    ln1 = ax.plot(time_ns_s3[::step], cv_site4[::step],
                  color="#2E9E9A", lw=0.6, alpha=0.7, label="cv_site4 (pocket opening)")
    ln2 = ax2.plot(time_ns_s3[::step], pip2_rmsd_A[::step],
                   color="#E07020", lw=0.6, alpha=0.7, label="PIP2-helix RMSD")
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("cv_site4 (nm)", color="#2E9E9A")
    ax2.set_ylabel("PIP2-helix RMSD (Å)", color="#E07020")
    ax.tick_params(axis='y', labelcolor="#2E9E9A")
    ax2.tick_params(axis='y', labelcolor="#E07020")
    lns = ln1 + ln2
    ax.legend(lns, [l.get_label() for l in lns], fontsize=8, loc="upper left")
    ax.set_title("A   PICALM site473 opening vs PIP2-helix displacement", fontweight="bold", fontsize=10)
    for spine in ["top"]:
        ax.spines[spine].set_visible(False); ax2.spines[spine].set_visible(False)

    # Panel B: scatter (subsample to ~2000 points)
    ax = axes[1]
    step_sc = max(1, len(cv_site4) // 2000)
    sc = ax.scatter(cv_site4[::step_sc], pip2_rmsd_A[::step_sc],
                    c=time_ns_s3[::step_sc], cmap="plasma",
                    s=6, alpha=0.4, linewidths=0)
    plt.colorbar(sc, ax=ax, label="Time (ns)", fraction=0.046, pad=0.04)
    m3, b3 = np.polyfit(cv_site4, pip2_rmsd_A, 1)
    xfit3 = np.linspace(cv_site4.min(), cv_site4.max(), 100)
    ax.plot(xfit3, m3 * xfit3 + b3, "r--", lw=1.5, alpha=0.9)
    ax.set_xlabel("cv_site4 — site473 opening (nm)")
    ax.set_ylabel("PIP2-binding helix RMSD (Å)")
    ax.set_title(f"B   Pearson r = {r3:.3f}  (p = {p3:.1e})", fontweight="bold", fontsize=10)
    ax.text(0.97, 0.05, f"r = {r3:.3f}", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=11, fontweight="bold", color="firebrick")
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

    plt.suptitle("PICALM ANTH — Site473 Pocket Opening vs PIP2-Binding Helix Displacement\n"
                 "200 ns metadynamics", fontsize=11, y=1.01)
    plt.tight_layout()
    plt.savefig(str(outdir / "pip2_loop.png"), dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {outdir}/pip2_loop.png")

    summary_lines += [
        "PICALM ANTH PIP2-Loop:",
        f"  cv_site4 range: {cv_site4.min():.2f}–{cv_site4.max():.2f} nm (mean {cv_site4.mean():.2f})",
        f"  PIP2-helix RMSD: {pip2_rmsd_A.min():.2f}–{pip2_rmsd_A.max():.2f} Å (mean {pip2_rmsd_A.mean():.2f})",
        f"  Pearson r(cv_site4, RMSD) = {r3:.3f}  p = {p3:.2e}",
        "  Interpretation: " + (
            "site473 opening correlates with PIP2-helix displacement — mechanistic link to PIP2 binding"
            if r3 > 0.3 else
            "site473 opening does not strongly displace the PIP2-binding helix in this simulation"
            if abs(r3) < 0.2 else
            "site473 opening anti-correlates with PIP2-helix displacement"
        ),
        "",
    ]
    manifest_data["pip2_loop"] = {
        "pearson_r": float(round(r3, 4)),
        "pearson_p": float(f"{p3:.2e}"),
        "cv_site4_mean_nm": float(round(cv_site4.mean(), 3)),
        "pip2_rmsd_mean_A": float(round(pip2_rmsd_A.mean(), 3)),
    }

    # ══════════════════════════════════════════════════════════════════════════
    # 5. Holo DCCM — S3 PICALM ANTH (Ibufenac complex, rep1 + rep5, 100% engaged)
    # ══════════════════════════════════════════════════════════════════════════
    print("\n" + "="*60)
    print("5. Holo DCCM — S3 PICALM ANTH (Ibufenac complex)")
    print("="*60)

    COMPLEX_BASE = project_root / "02_md_simulations/complex_md/S3_site473_CHEMBL341812"
    HOLO_REPS    = ["rep1", "rep5"]   # 100% engaged reps only

    holo_dccms = []
    for rep in HOLO_REPS:
        tpr_h = COMPLEX_BASE / rep / "ref_nw.pdb"
        xtc_h = COMPLEX_BASE / rep / "complexmd_fit_nw.xtc"
        print(f"  Computing holo DCCM: {rep} ...")
        dccm_h, res_h, rn_h, seg_h = compute_dccm(tpr_h, xtc_h, sel="protein and name CA")
        print(f"    {len(res_h)} Cα atoms, {dccm_h.shape} matrix")
        holo_dccms.append((dccm_h, res_h, rn_h, seg_h))

    # Match residues between apo (metadynamics) and holo (complex MD)
    n_apo  = len(res_s3)
    n_holo = len(holo_dccms[0][1])
    if n_apo != n_holo:
        print(f"  NOTE: apo {n_apo} Cα vs holo {n_holo} Cα — restricting to common residues")
        common_res = np.intersect1d(res_s3, holo_dccms[0][1])
        apo_mask   = np.isin(res_s3, common_res)
        holo_mask  = np.isin(holo_dccms[0][1], common_res)
        dccm_apo_use  = dccm_s3[np.ix_(apo_mask, apo_mask)]
        dccm_holo_avg = np.mean([d[np.ix_(holo_mask, holo_mask)] for d, *_ in holo_dccms], axis=0)
        res_cmn = res_s3[apo_mask]
        seg_cmn = seg_s3[apo_mask]
    else:
        dccm_apo_use  = dccm_s3
        dccm_holo_avg = np.mean([d for d, *_ in holo_dccms], axis=0)
        res_cmn = res_s3
        seg_cmn = seg_s3

    np.clip(dccm_holo_avg, -1.0, 1.0, out=dccm_holo_avg)

    # Use HOLO residue numbers for all labeling (matches contact analysis numbering).
    # Apo is offset by -18 (apo resid = holo resid - 18), same 268-residue sequence.
    res_holo_ref = holo_dccms[0][1]   # holo residue numbers (19–286)
    seg_holo_ref = holo_dccms[0][3]   # segids from holo
    print(f"  Holo DCCM ({len(res_holo_ref)} Cα, resids {res_holo_ref[0]}–{res_holo_ref[-1]}): "
          f"mean|C|={np.abs(dccm_holo_avg).mean():.3f}")

    # Plot holo DCCM — label axes in holo residue numbers
    plot_dccm(dccm_holo_avg, res_holo_ref, seg_holo_ref,
              "DCCM — PICALM ANTH holo (Ibufenac, rep1+rep5)",
              outdir / "dccm_S3_PICALM_ANTH_holo.png",
              pocket_resids=[219, 237, 240, 267, 268, 269, 270])

    # Compute ΔDCCM (positional subtraction — sequences are identical, offset=18)
    delta_dccm = dccm_holo_avg - dccm_apo_use
    print(f"  ΔDCCM: min={delta_dccm.min():.3f}  max={delta_dccm.max():.3f}  "
          f"mean|Δ|={np.abs(delta_dccm).mean():.3f}")

    # ══════════════════════════════════════════════════════════════════════════
    # 6. ΔDCCM three-panel figure: apo | holo | holo−apo
    # ══════════════════════════════════════════════════════════════════════════
    print("\n" + "="*60)
    print("6. ΔDCCM — apo vs holo")
    print("="*60)

    # Residue index bands — use HOLO residue numbers for lookup.
    # PIP2-helix: apo resids 1–22 = holo resids 19–40.
    pip2_idx   = np.where((res_holo_ref >= 19)  & (res_holo_ref <= 40))[0]
    pocket_res = [219, 237, 240, 267, 268, 269, 270]   # core Ibufenac contacts (holo numbering)
    pocket_idx = np.array([np.where(res_holo_ref == r)[0][0]
                           for r in pocket_res if r in res_holo_ref])

    fig, axes = plt.subplots(1, 3, figsize=(16, 5.5))
    n_cmn = len(res_holo_ref)
    tick_idx = np.arange(0, n_cmn, max(1, n_cmn // 12))

    panels = [
        (axes[0], dccm_apo_use,  "A   Apo (metadynamics)",           -1.0, 1.0, "RdBu_r", "C"),
        (axes[1], dccm_holo_avg, "B   Holo (Ibufenac, rep1+rep5)",   -1.0, 1.0, "RdBu_r", "C"),
        (axes[2], delta_dccm,    "C   ΔDCCM (holo − apo)",           -0.7, 0.7, "RdBu_r", "ΔC"),
    ]

    for ax, mat, title, vmin, vmax, cmap, cb_label in panels:
        im = ax.imshow(mat, cmap=cmap, vmin=vmin, vmax=vmax,
                       origin="lower", aspect="auto", interpolation="nearest")
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label=cb_label)

        # PIP2-helix band (cyan dotted lines)
        if len(pip2_idx):
            for edge in [pip2_idx[0] - 0.5, pip2_idx[-1] + 0.5]:
                ax.axhline(edge, color="cyan", lw=0.9, ls=":", alpha=0.8)
                ax.axvline(edge, color="cyan", lw=0.9, ls=":", alpha=0.8)

        # Pocket residue lines (yellow)
        for idx in pocket_idx:
            ax.axhline(idx, color="#FFDD00", lw=0.7, alpha=0.65)
            ax.axvline(idx, color="#FFDD00", lw=0.7, alpha=0.65)

        ax.set_xticks(tick_idx)
        ax.set_xticklabels([str(res_holo_ref[i]) for i in tick_idx], fontsize=6, rotation=45)
        ax.set_yticks(tick_idx)
        ax.set_yticklabels([str(res_holo_ref[i]) for i in tick_idx], fontsize=6)
        ax.set_xlabel("Residue", fontsize=9)
        ax.set_title(title, fontsize=10, fontweight="bold")

    axes[0].set_ylabel("Residue", fontsize=9)
    axes[0].text(0.02, 0.98,
                 "Cyan: PIP2-helix (19–40)\nYellow: site473 contacts",
                 transform=axes[0].transAxes, fontsize=7, va="top",
                 bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.75))

    plt.suptitle("PICALM ANTH — ΔDCCM: Apo vs Ibufenac Holo\n"
                 "Apo: 200 ns metadynamics | Holo: complex MD rep1+rep5 (100% engaged)",
                 fontsize=11, y=1.01)
    plt.tight_layout()
    out_delta = outdir / "dccm_delta_S3_PICALM.png"
    plt.savefig(str(out_delta), dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out_delta}")

    # Key statistic: ΔDCCM in PIP2 ↔ pocket block
    delta_pip2_pocket = None
    if len(pip2_idx) and len(pocket_idx):
        delta_pip2_pocket = delta_dccm[np.ix_(pip2_idx, pocket_idx)]
        print(f"  ΔDCCM PIP2-helix ↔ site473 pocket: "
              f"mean={delta_pip2_pocket.mean():.3f}  "
              f"max|Δ|={np.abs(delta_pip2_pocket).max():.3f}")

    summary_lines += [
        "PICALM ANTH ΔDCCM (Holo − Apo):",
        f"  Holo reps: {', '.join(HOLO_REPS)} (100% engaged)",
        f"  Residues compared: {len(res_holo_ref)} Cα",
        f"  ΔDCCM range: {delta_dccm.min():.3f} to {delta_dccm.max():.3f}",
        f"  Mean |ΔDCCM|: {np.abs(delta_dccm).mean():.3f}",
    ]
    if delta_pip2_pocket is not None:
        summary_lines.append(
            f"  ΔDCCM PIP2 ↔ pocket: mean={delta_pip2_pocket.mean():.3f}  "
            f"max|Δ|={np.abs(delta_pip2_pocket).max():.3f}"
        )
    summary_lines.append("")

    manifest_data["delta_dccm"] = {
        "holo_reps": HOLO_REPS,
        "n_residues": int(len(res_holo_ref)),
        "delta_min": float(round(delta_dccm.min(), 4)),
        "delta_max": float(round(delta_dccm.max(), 4)),
        "mean_abs_delta": float(round(np.abs(delta_dccm).mean(), 4)),
        **({"pip2_pocket_mean_delta": float(round(delta_pip2_pocket.mean(), 4)),
            "pip2_pocket_max_abs_delta": float(round(np.abs(delta_pip2_pocket).max(), 4))}
           if delta_pip2_pocket is not None else {}),
    }

    # ── Write summary ──────────────────────────────────────────────────────────
    summ_path = outdir / "conformational_summary.txt"
    with open(str(summ_path), "w") as f:
        f.write("\n".join(summary_lines))
    print(f"\nSummary: {summ_path}")
    print("All outputs in:", outdir)

    # ── Write manifest ─────────────────────────────────────────────────────────
    manifest = {
        "generated_at":  time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "script":        "09_scripts/07_conformational/05_conformational_coupling.py",
        "project_root":  str(project_root),
        "outputs": {
            "dccm_s1":     str(outdir / "dccm_S1_BIN1_BAR.png"),
            "dccm_s3_apo": str(outdir / "dccm_S3_PICALM_ANTH.png"),
            "dccm_s3_holo":str(outdir / "dccm_S3_PICALM_ANTH_holo.png"),
            "delta_dccm":  str(outdir / "dccm_delta_S3_PICALM.png"),
            "bar_curvature": str(outdir / "bar_curvature.png"),
            "pip2_loop":   str(outdir / "pip2_loop.png"),
            "summary_txt": str(summ_path),
        },
        "results": manifest_data,
    }
    manifest_path = outdir / "conformational_manifest.json"
    with open(str(manifest_path), "w") as fh:
        json.dump(manifest, fh, indent=2)
    print(f"Manifest → {manifest_path}")


if __name__ == "__main__":
    main()
