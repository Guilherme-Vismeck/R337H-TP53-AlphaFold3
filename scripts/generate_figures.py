"""Figuras: distancia 337-352 por modelo/seed e tipo de atomo mais proximo."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from structure_utils import START, all_models, load_chains

T = Path("results/tables")
F = Path("results/figures")
F.mkdir(parents=True, exist_ok=True)

BACKBONE = {"N", "CA", "C", "O", "OXT"}
I337, I352 = 337 - START, 352 - START
SEED_ORDER = [1, 2, 3, 4, 333]
COLORS = {"WT": "#5E4AE3", "R337H": "#F26CA7"}
CRYSTAL = 2.72  # N...O no 1C26 (medido em interface_analysis.py)

# ---- Dados do painel (a): ja calculados antes ----
pm = pd.read_csv(T / "key_interaction_per_model.csv")

# ---- Dados do painel (b): atomo mais proximo por subunidade ----
rows = []
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
        rows.append({"system": system, "seed": seed, "model": i, "subunit": "ABCD"[ci],
                     "atom_337": best[1], "atom_352": best[2], "dist": round(best[0], 2)})
at = pd.DataFrame(rows)
at["carbon_closest"] = at.atom_337.str.startswith("C")
at.to_csv(T / "closest_atoms_per_subunit.csv", index=False)
frac = at.groupby(["system", "seed"])["carbon_closest"].mean().unstack("system")
print("Fracao de subunidades cujo atomo mais proximo e um carbono:")
print(frac.reindex(SEED_ORDER).round(2))

# ---- Figura ----
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.8))
rng = np.random.default_rng(0)  # jitter fixo, so para visualizacao
off = {"WT": -0.18, "R337H": 0.18}

for system in ["WT", "R337H"]:
    for k, seed in enumerate(SEED_ORDER):
        y = pm[(pm.system == system) & (pm.seed == seed)].dist_polar_N_O.values
        x = k + off[system] + rng.uniform(-0.06, 0.06, len(y))
        ax1.scatter(x, y, s=38, color=COLORS[system], alpha=0.85,
                    label=system if k == 0 else None, edgecolor="white", linewidth=0.5)
ax1.axhline(CRYSTAL, color="gray", ls="--", lw=1, label="1C26 (cristal)")
ax1.set_xticks(range(len(SEED_ORDER)))
ax1.set_xticklabels([str(s) for s in SEED_ORDER])
ax1.set_xlabel("Seed")
ax1.set_ylabel("Distancia N...O, 337-352 (A)\n(media das 4 subunidades de cada modelo)")
ax1.set_title("(a) Distancia N...O por modelo (5 modelos por seed)")
ax1.legend(frameon=False)

w = 0.38
xs = np.arange(len(SEED_ORDER))
for system, shift in [("WT", -w / 2), ("R337H", w / 2)]:
    vals = [frac.loc[s, system] for s in SEED_ORDER]
    ax2.bar(xs + shift, vals, width=w, color=COLORS[system], label=system)
ax2.set_xticks(xs)
ax2.set_xticklabels([str(s) for s in SEED_ORDER])
ax2.set_ylim(0, 1)
ax2.set_xlabel("Seed")
ax2.set_ylabel("Fracao das subunidades")
ax2.set_title("(b) Atomo mais proximo do Asp352 e um carbono")
ax2.legend(frameon=False)
ax2.text(0.02, 0.96, "WT: zero em todas as seeds\n(o atomo mais proximo e sempre um N)",
         transform=ax2.transAxes, va="top", fontsize=9, color="gray")

fig.suptitle("Previsoes AF3 do tetramero p53 (325-356): 5 seeds por sistema", fontsize=11)
fig.tight_layout()
fig.savefig(F / "distance_337_352_by_seed.png", dpi=200)
print("figura salva em results/figures/distance_337_352_by_seed.png")