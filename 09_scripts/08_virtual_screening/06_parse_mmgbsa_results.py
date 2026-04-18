"""
CRYPTAD — Step 4.4b + 4.5 pre-processing
Parses all FINAL_RESULTS_MMPBSA.csv files, computes LE, CNS-MPO score,
and outputs ranked tables for ADMET filtering.

Outputs (04_virtual_screening/mmgbsa_results/):
  mmgbsa_all.csv          All results: ΔG, SD, n_frames, LE, CNS-MPO
  mmgbsa_le_filtered.csv  LE ≥ 0.3 compounds only
  top20_per_pocket.csv    Top 20 per pocket by LE (input to ADMET)
  pkcsm_input.csv         Name + SMILES for pkCSM web submission
  parse_mmgbsa_manifest.json

Usage:
  python3 09_scripts/08_virtual_screening/06_parse_mmgbsa_results.py
  python3 09_scripts/08_virtual_screening/06_parse_mmgbsa_results.py --project-root /path/to/CRYPTAD
"""

import argparse
import csv
import json
import statistics
import time
from pathlib import Path

from rdkit import Chem

_script       = Path(__file__).resolve()
_default_root = _script.parents[2]


# ── MM-GBSA delta parser ──────────────────────────────────────────────────────
def parse_delta(csv_path: Path) -> list:
    """Return list of per-frame TOTAL delta values (skips overflow frames)."""
    totals, in_delta = [], False
    with open(csv_path) as f:
        for line in f:
            line = line.strip()
            if line == "Delta Energy Terms":
                in_delta = True; continue
            if in_delta and line.startswith("Frame #"): continue
            if in_delta and line and line[0].isdigit():
                parts = line.split(",")
                try:
                    val = float(parts[-1])
                    if abs(val) < 5000:
                        totals.append(val)
                except (ValueError, IndexError):
                    break
    return totals


# ── CNS-MPO score (Wager et al. 2010, 6 desirability functions) ──────────────
def _d(val, low, high):
    """Linear desirability: 1 if val ≤ low, 0 if val ≥ high."""
    if val <= low:  return 1.0
    if val >= high: return 0.0
    return (high - val) / (high - low)


def cns_mpo(mw, logp, tpsa, hbd, pka=None):
    """
    Compute CNS-MPO score (0–6). pKa is the most basic center.
    If pKa is unknown, that term is omitted and score is scaled to /5.
    """
    d_mw   = _d(mw,   360, 500)    # 1 if ≤360, 0 if ≥500
    d_logp = _d(logp,   3,   5)    # 1 if ≤3,   0 if ≥5
    d_tpsa = _d(tpsa,  60, 120)    # 1 if ≤60,  0 if ≥120
    d_hbd  = _d(hbd,  0.5, 3.5)   # 1 if HBD=0, 0 if HBD≥4
    d_logd = _d(logp,   2,   5)   # LogD proxy (use LogP); 1 if ≤2, 0 if ≥5

    terms = [d_mw, d_logp, d_tpsa, d_hbd, d_logd]

    if pka is not None:
        d_pka = _d(pka, 8, 10)    # 1 if pKa≤8, 0 if pKa≥10
        terms.append(d_pka)
        return sum(terms), 6
    return sum(terms), 5


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="CRYPTAD MM-GBSA results parser (Step 4.4b)")
    parser.add_argument("--project-root", type=Path, default=_default_root,
                        help=f"Path to CRYPTAD project root (default: {_default_root})")
    parser.add_argument("--le-threshold", type=float, default=0.3,
                        help="Ligand efficiency pass threshold (default: 0.3)")
    parser.add_argument("--top-n", type=int, default=20,
                        help="Top N per pocket for ADMET list (default: 20)")
    args = parser.parse_args()

    project_root = args.project_root.resolve()
    mmgbsa_dir   = project_root / "04_virtual_screening" / "mmgbsa"
    manifest_csv = project_root / "04_virtual_screening" / "ligands" / "library_manifest.csv"
    outdir       = project_root / "04_virtual_screening" / "mmgbsa_results"
    outdir.mkdir(parents=True, exist_ok=True)

    # Auto-discover pockets from mmgbsa directory
    pockets = sorted(p.name for p in mmgbsa_dir.iterdir()
                     if p.is_dir() and p.name != "mdp")
    if not pockets:
        print(f"[ERROR] No pocket directories found in {mmgbsa_dir}")
        return
    print(f"Pockets: {pockets}")

    # Load manifest (id → properties)
    manifest = {}
    with open(manifest_csv) as f:
        for row in csv.DictReader(f):
            manifest[row["id"]] = row

    # Parse MM-GBSA results
    rows = []
    for pocket in pockets:
        pdir = mmgbsa_dir / pocket
        for lig_dir in sorted(pdir.iterdir()):
            if not lig_dir.is_dir():
                continue
            csv_path = lig_dir / "mmgbsa" / "FINAL_RESULTS_MMPBSA.csv"
            if not csv_path.exists():
                continue
            totals = parse_delta(csv_path)
            if not totals:
                continue
            dg_mean = statistics.mean(totals)
            dg_sd   = statistics.stdev(totals) if len(totals) > 1 else 0.0

            lig = lig_dir.name
            m = manifest.get(lig, {})
            smiles = m.get("smiles", "")
            name   = m.get("name", lig)
            try:
                mw   = float(m.get("mw",   0))
                logp = float(m.get("logp", 0))
                tpsa = float(m.get("tpsa", 0))
                hbd  = float(m.get("hbd",  0))
                hba  = float(m.get("hba",  0))
            except ValueError:
                mw = logp = tpsa = hbd = hba = 0.0

            mol = Chem.MolFromSmiles(smiles) if smiles else None
            hac = mol.GetNumHeavyAtoms() if mol else 0

            le = (-dg_mean / hac) if hac > 0 else 0.0

            mpo_sum, mpo_max = cns_mpo(mw, logp, tpsa, hbd)
            mpo_score = round(mpo_sum, 2)
            mpo_norm  = round(mpo_sum / mpo_max * 6, 2)

            rows.append({
                "pocket":    pocket,
                "chembl_id": lig,
                "name":      name,
                "smiles":    smiles,
                "dG_mean":   round(dg_mean, 2),
                "dG_sd":     round(dg_sd,   2),
                "n_frames":  len(totals),
                "MW":        mw,
                "LogP":      logp,
                "TPSA":      tpsa,
                "HBD":       hbd,
                "HBA":       hba,
                "HAC":       hac,
                "LE":        round(le, 3),
                "CNS_MPO_5": mpo_score,
                "CNS_MPO_6": mpo_norm,
            })

    print(f"Parsed {len(rows)} MM-GBSA results across {len(pockets)} pockets.")

    fields = ["pocket", "chembl_id", "name", "dG_mean", "dG_sd", "n_frames",
              "LE", "CNS_MPO_5", "CNS_MPO_6", "MW", "LogP", "TPSA", "HBD", "HBA",
              "HAC", "smiles"]
    rows_sorted = sorted(rows, key=lambda r: r["dG_mean"])

    # Write mmgbsa_all.csv
    out_all = outdir / "mmgbsa_all.csv"
    with open(out_all, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader(); w.writerows(rows_sorted)
    print(f"  → mmgbsa_all.csv  ({len(rows_sorted)} rows)")

    # LE filter
    le_pass = [r for r in rows_sorted if r["LE"] >= args.le_threshold]
    le_fail = [r for r in rows_sorted if r["LE"] <  args.le_threshold]
    out_le = outdir / "mmgbsa_le_filtered.csv"
    with open(out_le, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader(); w.writerows(le_pass)
    print(f"  → mmgbsa_le_filtered.csv  "
          f"({len(le_pass)} pass LE≥{args.le_threshold}, {len(le_fail)} fail)")

    # Top N per pocket by LE
    top_rows = []
    print(f"\nTop {args.top_n} per pocket (LE ≥ {args.le_threshold}, ranked by LE):")
    for pocket in pockets:
        pocket_rows = sorted(
            [r for r in le_pass if r["pocket"] == pocket],
            key=lambda r: r["LE"], reverse=True
        )[:args.top_n]
        top_rows.extend(pocket_rows)
        print(f"  {pocket}: {len(pocket_rows)} compounds")

    out_top = outdir / f"top{args.top_n}_per_pocket.csv"
    with open(out_top, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader(); w.writerows(top_rows)
    print(f"  → {out_top.name}  ({len(top_rows)} total)")

    # pkCSM input (unique compounds across pockets)
    seen = {}
    for r in top_rows:
        cid = r["chembl_id"]
        if cid not in seen or r["LE"] > seen[cid]["LE"]:
            seen[cid] = r
    pkcsm_rows = sorted(seen.values(), key=lambda r: r["LE"], reverse=True)
    out_pkcsm = outdir / "pkcsm_input.csv"
    with open(out_pkcsm, "w", newline="") as f:
        w = csv.DictWriter(f,
                           fieldnames=["chembl_id", "name", "smiles",
                                       "MW", "LogP", "TPSA", "HBD"],
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(pkcsm_rows)
    print(f"  → pkcsm_input.csv  ({len(pkcsm_rows)} unique compounds for ADMET submission)")

    # Console summary
    print("\n── Top 10 overall by LE ─────────────────────────────────────────────")
    print(f"{'Rank':>4}  {'Compound':<20}  {'Pocket':<15}  "
          f"{'dG':>7}  {'LE':>6}  {'MPO':>5}  {'MW':>5}")
    print("─" * 75)
    for i, r in enumerate(sorted(le_pass, key=lambda r: r["LE"], reverse=True)[:10], 1):
        print(f"{i:>4}  {r['chembl_id']:<20}  {r['pocket']:<15}  "
              f"{r['dG_mean']:>7.1f}  {r['LE']:>6.3f}  "
              f"{r['CNS_MPO_6']:>5.2f}  {r['MW']:>5.0f}")

    print(f"\n── LE distribution ──────────────────────────────────────────────────")
    les = [r["LE"] for r in rows]
    for (lo, hi), lab in [
        ((0.5, 99),  "≥0.5 (excellent)"),
        ((0.4, 0.5), "0.4–0.5 (good)"),
        ((0.3, 0.4), "0.3–0.4 (acceptable)"),
        ((0.0, 0.3), "<0.3 (fail)"),
    ]:
        n = sum(1 for v in les if lo <= v < hi)
        print(f"  {lab}: {n}")

    # Manifest
    manifest_out = {
        "generated_at":  time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "script":        "09_scripts/08_virtual_screening/06_parse_mmgbsa_results.py",
        "project_root":  str(project_root),
        "pockets":       pockets,
        "le_threshold":  args.le_threshold,
        "top_n":         args.top_n,
        "n_total":       len(rows),
        "n_le_pass":     len(le_pass),
        "n_unique_admet": len(pkcsm_rows),
        "outputs": {
            "all":     str(out_all),
            "le":      str(out_le),
            "top":     str(out_top),
            "pkcsm":   str(out_pkcsm),
        },
    }
    manifest_path = outdir / "parse_mmgbsa_manifest.json"
    with open(manifest_path, "w") as fh:
        json.dump(manifest_out, fh, indent=2)
    print(f"\nManifest → {manifest_path}")


if __name__ == "__main__":
    main()
