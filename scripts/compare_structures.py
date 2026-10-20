"""RMSD de CA entre tetrameros: AF3 (50 modelos) vs 1C26 e WT vs R337H."""
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

from structure_utils import (REF, SYSTEMS, all_models, best_alignment,
                             ca_array, load_chains)

OUT = Path("results/tables")
OUT.mkdir(parents=True, exist_ok=True)

ref_ca = ca_array(load_chains(REF))
cas = {(s, seed, i): ca_array(load_chains(p)) for s, seed, i, p in all_models()}
keys = {s: [k for k in cas if k[0] == s] for s in SYSTEMS}

rows = []


def compare(label, ka, ca_a, name_b, ca_b, kb=None):
    rmsd, perm, _ = best_alignment(ca_a, ca_b)
    rows.append({
        "comparison": label,
        "a": f"{ka[0]}_s{ka[1]}_m{ka[2]}",
        "b": name_b,
        "rmsd_tetramer": round(rmsd, 3),
    })


for s in SYSTEMS:
    for k in keys[s]:
        compare(f"{s}_vs_1C26", k, cas[k], "1C26", ref_ca)

for ka in keys["WT"]:
    for kb in keys["R337H"]:
        compare("WT_vs_R337H", ka, cas[ka], f"R337H_s{kb[1]}_m{kb[2]}", cas[kb])

for s in SYSTEMS:
    for ka, kb in combinations(keys[s], 2):
        tipo = "mesma_seed" if ka[1] == kb[1] else "seeds_diferentes"
        compare(f"{s}_vs_{s}_{tipo}", ka, cas[ka], f"{s}_s{kb[1]}_m{kb[2]}", cas[kb])
    print(f"comparacoes de {s} concluidas")

df = pd.DataFrame(rows)
df.to_csv(OUT / "model_comparison.csv", index=False)
summary = df.groupby("comparison")["rmsd_tetramer"].agg(
    ["count", "median", "min", "max"]).round(3)
summary.to_csv(OUT / "model_comparison_summary.csv")
print(summary)