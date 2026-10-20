"""RMSD de CA entre tetrameros: modelos AF3 vs 1C26 e WT vs R337H."""
from itertools import combinations, product
from pathlib import Path

import numpy as np
import pandas as pd

from structure_utils import (N_MODELS, REF, SYSTEMS, best_alignment, ca_array,
                             load_chains, model_path)

OUT = Path("results/tables")
OUT.mkdir(parents=True, exist_ok=True)

ref_ca = ca_array(load_chains(REF))
cas = {(s, i): ca_array(load_chains(model_path(s, i)))
       for s in SYSTEMS for i in range(N_MODELS)}

rows = []


def compare(label, name_a, ca_a, name_b, ca_b):
    rmsd, perm, fit = best_alignment(ca_a, ca_b)
    per_chain = np.sqrt(((fit - ca_b) ** 2).sum(axis=2).mean(axis=1))
    rows.append({
        "comparison": label, "a": name_a, "b": name_b,
        "rmsd_tetramer": round(rmsd, 3),
        "rmsd_chain_mean": round(per_chain.mean(), 3),
        "rmsd_chain_max": round(per_chain.max(), 3),
        "chain_permutation": "".join("ABCD"[p] for p in perm),
    })


for s in SYSTEMS:
    for i in range(N_MODELS):
        compare(f"{s}_vs_1C26", f"{s}_{i}", cas[(s, i)], "1C26", ref_ca)

for i, j in product(range(N_MODELS), repeat=2):
    compare("WT_vs_R337H", f"WT_{i}", cas[("WT", i)], f"R337H_{j}", cas[("R337H", j)])

for s in SYSTEMS:  # linha de base: variacao entre modelos do mesmo sistema
    for i, j in combinations(range(N_MODELS), 2):
        compare(f"{s}_vs_{s}", f"{s}_{i}", cas[(s, i)], f"{s}_{j}", cas[(s, j)])

df = pd.DataFrame(rows)
df.to_csv(OUT / "model_comparison.csv", index=False)

summary = df.groupby("comparison")["rmsd_tetramer"].agg(["count", "median", "min", "max"]).round(3)
summary.to_csv(OUT / "model_comparison_summary.csv")
print(summary)