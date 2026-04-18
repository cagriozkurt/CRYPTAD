"""
CRYPTAD — STRING Protein Interaction Network
=============================================
Fetches interaction data from STRING (v12, Homo sapiens) for the three
CRYPTAD study proteins plus a curated set of AD-relevant neighbours.

Proteins queried:
  Study targets : BIN1, PICALM, CD2AP
  AD pathway    : APP, APOE, MAPT, TREM2, CLU, SORL1, CR1, ABCA7, BACE1,
                  DNM1L, SH3GL2, AAK1, SYNJ1, ITSN1

Score threshold: ≥ 700 (high confidence)

Node colours:
  CRYPTAD targets  — gold  (#F0A500)
  AD GWAS risk     — orange (#E07800)
  Endocytic / APP  — steel  (#4A90D9)
  Other            — grey   (#AAAAAA)

Output: 06_figures/fig_string_network.png  (300 dpi)

Usage:
  python3 09_scripts/10_figures/04_string_network.py
  python3 09_scripts/10_figures/04_string_network.py --project-root /path/to/CRYPTAD
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import requests
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

_script       = Path(__file__).resolve()
_default_root = _script.parents[2]

# ── Protein list ──────────────────────────────────────────────────────────────
STUDY_TARGETS = {"BIN1", "PICALM", "CD2AP"}

SEED_PROTEINS = [
    "BIN1", "PICALM", "CD2AP",
    "APP", "APOE", "TREM2", "CLU", "CR1", "ABCA7", "SORL1", "MAPT",
    "DNM1", "SH3GL2", "AAK1", "SYNJ1", "ITSN1", "EPS15",
    "BACE1", "PSEN1", "PSEN2",
]

ROLE = {
    "BIN1":   "study",   "PICALM": "study",   "CD2AP":  "study",
    "APP":    "gwas",    "APOE":   "gwas",     "TREM2":  "gwas",
    "CLU":    "gwas",    "CR1":    "gwas",     "ABCA7":  "gwas",
    "SORL1":  "gwas",    "MAPT":   "gwas",
    "DNM1":   "endo",    "SH3GL2": "endo",     "AAK1":   "endo",
    "SYNJ1":  "endo",    "ITSN1":  "endo",     "EPS15":  "endo",
    "BACE1":  "app",     "PSEN1":  "app",      "PSEN2":  "app",
}
COLOR_MAP = {
    "study": "#F0A500",
    "gwas":  "#E07800",
    "endo":  "#4A90D9",
    "app":   "#9B59B6",
    "other": "#AAAAAA",
}
SIZE_MAP = {
    "study": 1800,
    "gwas":  900,
    "endo":  700,
    "app":   700,
    "other": 500,
}

STRING_URL = "https://string-db.org/api/json/network"
SPECIES    = 9606
MIN_SCORE  = 700


def edge_color(u: str, v: str) -> str:
    ru = ROLE.get(u, "other")
    rv = ROLE.get(v, "other")
    if ru == "study" and rv == "study":
        return "#F0A500"
    if ru == "study" or rv == "study":
        return "#4A90D9"
    return "#CCCCCC"


def main():
    parser = argparse.ArgumentParser(
        description="CRYPTAD STRING protein interaction network figure")
    parser.add_argument("--project-root", type=Path, default=_default_root,
                        help=f"Path to CRYPTAD project root (default: {_default_root})")
    args = parser.parse_args()

    project_root = args.project_root.resolve()
    outdir       = project_root / "06_figures"
    outdir.mkdir(parents=True, exist_ok=True)

    # ── STRING API query ──────────────────────────────────────────────────────
    identifiers = "%0d".join(SEED_PROTEINS)
    print(f"Querying STRING v12 for {len(SEED_PROTEINS)} proteins "
          f"(min score {MIN_SCORE}) ...", flush=True)

    params = {
        "identifiers":     identifiers,
        "species":         SPECIES,
        "required_score":  MIN_SCORE,
        "caller_identity": "CRYPTAD_PhD",
    }
    try:
        resp = requests.get(STRING_URL, params=params, timeout=30)
        resp.raise_for_status()
        edges_raw = resp.json()
    except Exception as e:
        print(f"  [ERROR] STRING API: {e}")
        sys.exit(1)

    print(f"  {len(edges_raw)} raw interaction records returned")

    # ── Build graph ───────────────────────────────────────────────────────────
    G = nx.Graph()
    for p in SEED_PROTEINS:
        G.add_node(p)

    for e in edges_raw:
        a     = e.get("preferredName_A", "")
        b     = e.get("preferredName_B", "")
        score = float(e.get("score", 0))
        if a not in G.nodes and b not in G.nodes:
            continue
        if a not in G.nodes:
            G.add_node(a)
        if b not in G.nodes:
            G.add_node(b)
        if G.has_edge(a, b):
            G[a][b]["score"] = max(G[a][b]["score"], score)
        else:
            G.add_edge(a, b, score=score)

    isolated = [n for n in list(G.nodes) if G.degree(n) == 0]
    G.remove_nodes_from(isolated)
    print(f"  Graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges "
          f"({len(isolated)} isolated nodes removed)")

    # ── Layout ────────────────────────────────────────────────────────────────
    pos_seed = {
        "BIN1":   np.array([-0.4,  0.3]),
        "PICALM": np.array([ 0.4,  0.3]),
        "CD2AP":  np.array([ 0.0, -0.4]),
    }
    fixed_pos   = {k: v for k, v in pos_seed.items() if k in G.nodes}
    fixed_nodes = list(fixed_pos.keys())

    rng = np.random.default_rng(42)
    pos = nx.spring_layout(
        G,
        k=2.8 / np.sqrt(G.number_of_nodes()),
        iterations=120,
        fixed=fixed_nodes if set(fixed_nodes).issubset(G.nodes) else None,
        pos={**fixed_pos,
             **{n: (rng.uniform(-1, 1), rng.uniform(-1, 1))
                for n in G.nodes if n not in fixed_pos}},
        weight="score",
        seed=42,
    )

    # ── Draw ──────────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(13, 11))
    ax.set_aspect("equal")
    ax.axis("off")

    scores  = [G[u][v]["score"] for u, v in G.edges()]
    max_sc  = max(scores) if scores else 1.0
    min_sc  = min(scores) if scores else 0.0
    norm_sc = [(s - min_sc) / (max_sc - min_sc + 1e-9) for s in scores]

    for (u, v), nsc in zip(G.edges(), norm_sc):
        ec = edge_color(u, v)
        ax.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]],
                color=ec, lw=0.6 + nsc * 3.0, alpha=0.25 + nsc * 0.55, zorder=1)

    node_list   = list(G.nodes())
    node_colors = [COLOR_MAP.get(ROLE.get(n, "other"), COLOR_MAP["other"]) for n in node_list]
    node_sizes  = [SIZE_MAP.get(ROLE.get(n, "other"),  SIZE_MAP["other"])  for n in node_list]

    nc = nx.draw_networkx_nodes(G, pos,
                                nodelist=node_list,
                                node_color=node_colors,
                                node_size=node_sizes,
                                edgecolors="white",
                                linewidths=1.2,
                                ax=ax)
    nc.set_zorder(2)

    for node in G.nodes():
        x, y = pos[node]
        fw   = "bold" if ROLE.get(node) == "study" else "normal"
        fs   = 9.5    if ROLE.get(node) == "study" else 7.5
        ax.text(x, y, node, ha="center", va="center",
                fontsize=fs, fontweight=fw, color="#222222", zorder=3)

    legend_elements = [
        mpatches.Patch(color=COLOR_MAP["study"],
                       label="CRYPTAD study targets (BIN1, PICALM, CD2AP)"),
        mpatches.Patch(color=COLOR_MAP["gwas"],
                       label="AD GWAS risk proteins (Bellenguez 2022)"),
        mpatches.Patch(color=COLOR_MAP["endo"],
                       label="Endocytic / clathrin-mediated trafficking"),
        mpatches.Patch(color=COLOR_MAP["app"],
                       label="APP processing / γ-secretase"),
        Line2D([0], [0], color="#F0A500", lw=2.5,
               label="Study target – study target interaction"),
        Line2D([0], [0], color="#4A90D9", lw=1.5,
               label="Study target – neighbour interaction"),
        Line2D([0], [0], color="#CCCCCC", lw=1.0,
               label="Background interaction"),
    ]
    fig.legend(handles=legend_elements, loc="lower center",
               fontsize=8.5, framealpha=0.9, frameon=True,
               ncol=2, bbox_to_anchor=(0.5, -0.04), borderpad=0.8)

    ax.set_title(
        "STRING Protein Interaction Network — BIN1, PICALM, CD2AP and AD Pathway Context\n"
        "Homo sapiens · STRING v12 · Score ≥ 0.700 (high confidence)",
        fontsize=11, fontweight="bold", pad=14
    )

    sm = plt.cm.ScalarMappable(
        cmap="Greys",
        norm=plt.Normalize(vmin=min_sc / 1000, vmax=max_sc / 1000))
    sm.set_array([])
    cb = plt.colorbar(sm, ax=ax, fraction=0.018, pad=0.01,
                      label="Interaction confidence score")
    cb.ax.yaxis.label.set_fontsize(9)

    plt.tight_layout()
    out = outdir / "fig_string_network.png"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"\nFigure saved: {out.relative_to(project_root)}")


if __name__ == "__main__":
    main()
