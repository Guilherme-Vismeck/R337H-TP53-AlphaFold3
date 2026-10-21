"""Quais atomos formam o contato mais proximo entre a cadeia lateral do 337 e a do Asp352?"""
from collections import Counter

import numpy as np

from structure_utils import START, all_models, load_chains

BACKBONE = {"N", "CA", "C", "O", "OXT"}
I337, I352 = 337 - START, 352 - START

counts, dists = Counter(), {}
for system, seed, i, path in all_models():
    chains = load_chains(path)
    for ci, (_, res_i) in enumerate(chains):
        a = [x for x in res_i[I337] if x.get_id() not in BACKBONE and x.element != "H"]
        best = None
        for cj, (_, res_j) in enumerate(chains):
            if cj == ci:
                continue
            for b in res_j[I352]:
                if b.get_id() in BACKBONE or b.element == "H":
                    continue
                for x in a:
                    d = float(np.linalg.norm(x.coord - b.coord))
                    if best is None or d < best[0]:
                        best = (d, x.get_id(), b.get_id())
        key = (system, best[1], best[2])
        counts[key] += 1
        dists.setdefault(key, []).append(best[0])

print("sistema  atomo_337  atomo_352  n  dist_mediana")
for key, n in sorted(counts.items(), key=lambda kv: (kv[0][0], -kv[1])):
    print(f"{key[0]:6s} {key[1]:6s} {key[2]:6s} {n:4d}  {np.median(dists[key]):.2f}")