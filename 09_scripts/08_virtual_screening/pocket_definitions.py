"""
pocket_definitions.py  —  CRYPTAD
Confirmed cryptic pocket definitions for virtual screening.

4 pockets confirmed by: metadynamics sampling + fpocket detection +
persistence gate (≥20% in ≥1 unbiased replica).

Centroids are the mean fpocket alpha-sphere positions from the best-score
metadynamics frame (same frame used for the representative PDB).
Grid box size 25×25×25 Å covers the full pocket cavity with margin.
"""

POCKETS = [
    {
        "id":       "S1_site688",
        "system":   "S1_BIN1_BAR",
        "protein":  "BIN1",
        "domain":   "BAR homodimer",
        "priority": 1,                      # 2/3 replicas PASS — most reproducible
        "pdb":      "03_pocket_analysis/cryptic_inspection/S1_BIN1_BAR/site688_frame302.pdb",
        "frame":    302,
        "time_ns":  151.0,
        "cx": 78.03, "cy": 92.58, "cz": 103.52,
        "box_x": 25, "box_y": 25, "box_z": 25,
        "max_persist": 0.228,
        "note": "Lead. Inter-helix cleft, chain B membrane-binding surface. "
                "PHE54/LEU202/LEU206/CYS57/GLU52 lining. P2Rank rank 2.",
    },
    {
        "id":       "S1_site680",
        "system":   "S1_BIN1_BAR",
        "protein":  "BIN1",
        "domain":   "BAR homodimer",
        "priority": 2,
        "pdb":      "03_pocket_analysis/cryptic_inspection/S1_BIN1_BAR/site680_frame318.pdb",
        "frame":    318,
        "time_ns":  159.0,
        "cx": 63.84, "cy": 95.01, "cz": 94.20,
        "box_x": 25, "box_y": 25, "box_z": 25,
        "max_persist": 0.249,
        "note": "Secondary. Single-replica gate pass. Adjacent to site688 on BAR surface.",
    },
    {
        "id":       "S1_site1813",
        "system":   "S1_BIN1_BAR",
        "protein":  "BIN1",
        "domain":   "BAR homodimer",
        "priority": 3,
        "pdb":      "03_pocket_analysis/cryptic_inspection/S1_BIN1_BAR/site1813_frame288.pdb",
        "frame":    288,
        "time_ns":  144.0,
        "cx": 69.40, "cy": 77.34, "cz": 112.73,
        "box_x": 25, "box_y": 25, "box_z": 25,
        "max_persist": 0.363,
        "note": "Highest single-replica persistence. Inspect in PyMOL before committing to VS.",
    },
    {
        "id":       "S3_site473",
        "system":   "S3_PICALM_ANTH",
        "protein":  "PICALM",
        "domain":   "ANTH monomer",
        "priority": 4,
        "pdb":      "03_pocket_analysis/cryptic_inspection/S3_PICALM_ANTH/site473_frame384.pdb",
        "frame":    384,
        "time_ns":  192.0,
        "cx": 65.41, "cy": 64.13, "cz": 59.47,
        "box_x": 25, "box_y": 25, "box_z": 25,
        "max_persist": 0.204,
        "note": "Borderline gate pass. P2Rank not confirmed. Aβ-transcytosis axis.",
    },
]

# Docking parameters (shared across all pockets)
VINA_PARAMS = {
    "exhaustiveness": 16,
    "num_modes":      20,
    "energy_range":   5,
    "cpu":            8,
}
