"""
CRYPTAD — W13 ADMET filtering and final hit list.

Input:  04_virtual_screening/mmgbsa_results/mmgbsa_all.csv
        (MM-GBSA ΔG, CNS_MPO_5/6, MW, LogP, TPSA, HBD, HBA, SMILES per compound)

Adds via RDKit:
  - RotBonds         : rotatable bond count (Veber rule: ≤ 10)
  - Lipinski_viol    : count of Lipinski RO5 violations (MW>500, cLogP>5, HBD>5, HBA>10)
  - logBBB_pred      : predicted log(brain/blood) via Young et al. (1988)
                       logBBB = 0.152·cLogP − 0.0148·TPSA + 0.139 ; BBB+ if > −1.0
  - PAINS_alert      : 1 if any PAINS A/B/C alert present (Baell & Holloway, 2010)
  - PAINS_detail     : name of PAINS pattern(s) hit, or ""
  - Brenk_alert      : 1 if Brenk structural alert present (reactive/unstable groups)
  - pass_admet       : 1 if compound passes all filters below

Filters for final hit list (pass_admet = 1):
  - CNS_MPO_5  ≥ 3          (CNS drug-likeness; Wager et al., 2010)
  - MW         ≤ 500         (Lipinski)
  - LogP       ≤ 5 and ≥ -1 (CNS range)
  - TPSA       ≤ 90          (BBB penetration)
  - HBD        ≤ 3           (CNS penetration)
  - RotBonds   ≤ 10          (Veber oral bioavailability)
  - Lipinski_viol ≤ 1        (allow one violation)
  - logBBB_pred > -1.0       (BBB permeable)
  - PAINS_alert == 0         (no structural interference flags)
  - dG_mean    ≤ -20 kcal/mol (minimum binding energy threshold)

Outputs (04_virtual_screening/mmgbsa_results/):
  - admet_results.csv   : all compounds + new columns
  - final_hits.csv      : passing compounds, sorted
  - admet_summary.txt   : funnel summary for paper
  - admet_manifest.json

Usage:
  python3 09_scripts/08_virtual_screening/07_admet_filter.py
  python3 09_scripts/08_virtual_screening/07_admet_filter.py --project-root /path/to/CRYPTAD
"""

import argparse
import json
import time
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

_script       = Path(__file__).resolve()
_default_root = _script.parents[2]


def main():
    parser = argparse.ArgumentParser(
        description="CRYPTAD ADMET filtering (Step 4.5)")
    parser.add_argument("--project-root", type=Path, default=_default_root,
                        help=f"Path to CRYPTAD project root (default: {_default_root})")
    parser.add_argument("--dg-max", type=float, default=-20.0,
                        help="MM-GBSA ΔG upper bound (default: -20.0 kcal/mol)")
    parser.add_argument("--cns-mpo-min", type=float, default=3.0,
                        help="Minimum CNS_MPO_5 score (default: 3.0)")
    args = parser.parse_args()

    project_root = args.project_root.resolve()
    mmgbsa_dir   = project_root / "04_virtual_screening" / "mmgbsa_results"
    input_csv    = mmgbsa_dir / "mmgbsa_all.csv"

    # ── Thresholds ────────────────────────────────────────────────────────────
    CNS_MPO5_MIN      = args.cns_mpo_min
    MW_MAX            = 500
    LOGP_MAX          = 5
    LOGP_MIN          = -1
    TPSA_MAX          = 90
    HBD_MAX           = 3
    ROTBONDS_MAX      = 10
    LIPINSKI_VIOL_MAX = 1
    LOGBBB_MIN        = -1.0
    DG_MAX            = args.dg_max

    # ── Build PAINS + Brenk filter catalogs ──────────────────────────────────
    params_pains = FilterCatalogParams()
    params_pains.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
    params_pains.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_B)
    params_pains.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_C)
    pains_catalog = FilterCatalog(params_pains)

    params_brenk = FilterCatalogParams()
    params_brenk.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
    brenk_catalog = FilterCatalog(params_brenk)

    # ── Load data ─────────────────────────────────────────────────────────────
    df = pd.read_csv(input_csv)
    print(f"Loaded {len(df)} compounds from {input_csv.relative_to(project_root)}")

    # ── Compute RDKit properties ──────────────────────────────────────────────
    rot_bonds   = []
    lip_viols   = []
    logbbb      = []
    pains_flag  = []
    pains_names = []
    brenk_flag  = []

    for _, row in df.iterrows():
        smi = row["smiles"]
        mol = Chem.MolFromSmiles(smi) if isinstance(smi, str) else None

        if mol is None:
            rot_bonds.append(None)
            lip_viols.append(None)
            logbbb.append(None)
            pains_flag.append(None)
            pains_names.append("")
            brenk_flag.append(None)
            continue

        rot_bonds.append(rdMolDescriptors.CalcNumRotatableBonds(mol))

        viols = sum([
            row["MW"]   > 500,
            row["LogP"] > 5,
            row["HBD"]  > 5,
            row["HBA"]  > 10,
        ])
        lip_viols.append(viols)

        lbbb = 0.152 * row["LogP"] - 0.0148 * row["TPSA"] + 0.139
        logbbb.append(round(lbbb, 3))

        pains_matches = pains_catalog.GetMatches(mol)
        if pains_matches:
            pains_flag.append(1)
            pains_names.append("; ".join(m.GetDescription() for m in pains_matches))
        else:
            pains_flag.append(0)
            pains_names.append("")

        brenk_matches = brenk_catalog.GetMatches(mol)
        brenk_flag.append(1 if brenk_matches else 0)

    df["RotBonds"]      = rot_bonds
    df["Lipinski_viol"] = lip_viols
    df["logBBB_pred"]   = logbbb
    df["PAINS_alert"]   = pains_flag
    df["PAINS_detail"]  = pains_names
    df["Brenk_alert"]   = brenk_flag

    # ── Apply ADMET pass/fail ─────────────────────────────────────────────────
    def passes(row):
        try:
            return int(
                row["CNS_MPO_5"]     >= CNS_MPO5_MIN      and
                row["MW"]            <= MW_MAX             and
                LOGP_MIN <= row["LogP"] <= LOGP_MAX        and
                row["TPSA"]          <= TPSA_MAX           and
                row["HBD"]           <= HBD_MAX            and
                row["RotBonds"]      <= ROTBONDS_MAX       and
                row["Lipinski_viol"] <= LIPINSKI_VIOL_MAX  and
                row["logBBB_pred"]   >  LOGBBB_MIN         and
                row["PAINS_alert"]   == 0                  and
                row["dG_mean"]       <= DG_MAX
            )
        except Exception:
            return 0

    df["pass_admet"] = df.apply(passes, axis=1)

    # ── Write full results ────────────────────────────────────────────────────
    admet_csv = mmgbsa_dir / "admet_results.csv"
    df.to_csv(admet_csv, index=False)
    print(f"Full ADMET results: {admet_csv.name}  ({len(df)} compounds)")

    # ── Final hit list: passing compounds sorted by dG_mean per pocket ────────
    hits = df[df["pass_admet"] == 1].copy()
    hits = hits.sort_values(["pocket", "dG_mean"])

    final_cols = [
        "pocket", "chembl_id", "name",
        "dG_mean", "dG_sd", "LE",
        "CNS_MPO_5", "CNS_MPO_6",
        "MW", "LogP", "TPSA", "HBD", "HBA", "RotBonds",
        "Lipinski_viol", "logBBB_pred",
        "PAINS_alert", "Brenk_alert",
        "smiles",
    ]
    hits_csv = mmgbsa_dir / "final_hits.csv"
    hits[final_cols].to_csv(hits_csv, index=False)
    print(f"Final hit list:     {hits_csv.name}  ({len(hits)} compounds)")

    # ── Funnel summary ────────────────────────────────────────────────────────
    total          = len(df)
    mask_dg        = df["dG_mean"] <= DG_MAX
    mask_mpo       = mask_dg        & (df["CNS_MPO_5"] >= CNS_MPO5_MIN)
    mask_lipinski  = mask_mpo       & (df["MW"] <= MW_MAX) & (df["LogP"] <= LOGP_MAX) \
                                    & (df["LogP"] >= LOGP_MIN) & (df["HBD"] <= HBD_MAX) \
                                    & (df["TPSA"] <= TPSA_MAX)
    mask_veber     = mask_lipinski  & (df["RotBonds"] <= ROTBONDS_MAX)
    mask_bbb       = mask_veber     & (df["logBBB_pred"] > LOGBBB_MIN)

    n_after_dg       = int(mask_dg.sum())
    n_after_cnsmpo   = int(mask_mpo.sum())
    n_after_lipinski = int(mask_lipinski.sum())
    n_after_veber    = int(mask_veber.sum())
    n_after_bbb      = int(mask_bbb.sum())
    n_final          = len(hits)

    summary_path = mmgbsa_dir / "admet_summary.txt"
    with open(summary_path, "w") as f:
        f.write("CRYPTAD — ADMET Filtering Funnel Summary\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Input compounds (MM-GBSA scored)     : {total}\n")
        f.write(f"After ΔG ≤ {DG_MAX} kcal/mol         : {n_after_dg}\n")
        f.write(f"After CNS MPO ≥ {CNS_MPO5_MIN}              : {n_after_cnsmpo}\n")
        f.write(f"After Lipinski/CNS property filters  : {n_after_lipinski}\n")
        f.write(f"After Veber (RotBonds ≤ 10)          : {n_after_veber}\n")
        f.write(f"After predicted logBBB > {LOGBBB_MIN}     : {n_after_bbb}\n")
        f.write(f"After PAINS removal                  : {n_final}\n")
        f.write(f"\nFinal hits                           : {n_final}\n\n")

        f.write("Filters applied:\n")
        f.write(f"  CNS_MPO_5 >= {CNS_MPO5_MIN}  (Wager et al., J Med Chem, 2010)\n")
        f.write(f"  MW <= {MW_MAX} Da\n")
        f.write(f"  {LOGP_MIN} <= cLogP <= {LOGP_MAX}\n")
        f.write(f"  TPSA <= {TPSA_MAX} Å²\n")
        f.write(f"  HBD <= {HBD_MAX}\n")
        f.write(f"  RotBonds <= {ROTBONDS_MAX}  (Veber et al., J Med Chem, 2002)\n")
        f.write(f"  Lipinski violations <= {LIPINSKI_VIOL_MAX}\n")
        f.write(f"  Predicted logBBB > {LOGBBB_MIN}  (Young et al., J Med Chem, 1988)\n")
        f.write(f"  No PAINS alerts  (Baell & Holloway, J Med Chem, 2010)\n")
        f.write(f"  MM-GBSA ΔG_mean <= {DG_MAX} kcal/mol\n\n")

        f.write("Hits by pocket:\n")
        for pocket, grp in hits.groupby("pocket"):
            f.write(f"  {pocket}: {len(grp)} compounds\n")
            for _, r in grp.iterrows():
                f.write(f"    {r['chembl_id']:15s} {r['name'][:30]:30s} "
                        f"dG={r['dG_mean']:7.2f}  LE={r['LE']:.3f}  "
                        f"CNS_MPO5={r['CNS_MPO_5']:.1f}  logBBB={r['logBBB_pred']:.2f}\n")
            f.write("\n")

        pains_flagged = df[df["PAINS_alert"] == 1]
        if len(pains_flagged):
            f.write(f"PAINS-flagged compounds excluded ({len(pains_flagged)}):\n")
            for _, r in pains_flagged.iterrows():
                f.write(f"  {r['chembl_id']:15s} {r['name'][:30]:30s} — {r['PAINS_detail']}\n")

    print(f"Summary written:    {summary_path.name}")
    print(f"\n{'='*60}")
    print(f"ADMET funnel: {total} → {n_after_dg} (ΔG) → {n_after_cnsmpo} (CNS MPO) → {n_final} (final hits)")
    print(f"{'='*60}")

    manifest = {
        "generated_at":  time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "script":        "09_scripts/08_virtual_screening/07_admet_filter.py",
        "project_root":  str(project_root),
        "input":         str(input_csv),
        "thresholds": {
            "dg_max":            DG_MAX,
            "cns_mpo5_min":      CNS_MPO5_MIN,
            "mw_max":            MW_MAX,
            "logp_range":        [LOGP_MIN, LOGP_MAX],
            "tpsa_max":          TPSA_MAX,
            "hbd_max":           HBD_MAX,
            "rotbonds_max":      ROTBONDS_MAX,
            "lipinski_viol_max": LIPINSKI_VIOL_MAX,
            "logbbb_min":        LOGBBB_MIN,
        },
        "n_input":       total,
        "n_final":       n_final,
        "funnel": {
            "after_dg":       n_after_dg,
            "after_cns_mpo":  n_after_cnsmpo,
            "after_lipinski": n_after_lipinski,
            "after_veber":    n_after_veber,
            "after_bbb":      n_after_bbb,
            "after_pains":    n_final,
        },
        "outputs": {
            "admet_results": str(admet_csv),
            "final_hits":    str(hits_csv),
            "summary":       str(summary_path),
        },
    }
    manifest_path = mmgbsa_dir / "admet_manifest.json"
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2)
    print(f"Manifest → {manifest_path.name}")


if __name__ == "__main__":
    main()
