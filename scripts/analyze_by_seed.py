"""Distancia 337-352 por seed: a diferenca WT vs R337H se repete?"""
from pathlib import Path
import pandas as pd

T = Path("results/tables")
key = pd.read_csv(T / "key_interaction_337_352.csv")

key = key[key.system != "1C26"].copy()
key["seed"] = key.model // 10
key["seed"] = key["seed"].replace({333: 333})
key["model_in_seed"] = key.model % 10
key["gap"] = key.dist_polar_N_O - key.dist_sidechain_min  # >0: par mais proximo nao e o N...O

# media das 4 subunidades por modelo
per_model = (key.groupby(["system", "seed", "model_in_seed"])
             [["dist_polar_N_O", "dist_sidechain_min", "gap"]].mean().reset_index())
per_model.to_csv(T / "key_interaction_per_model.csv", index=False)

# resumo por seed (5 modelos cada)
per_seed = (per_model.groupby(["system", "seed"])["dist_polar_N_O"]
            .agg(["median", "min", "max"]).round(2))
print("Distancia polar N...O (A) por seed:\n", per_seed)
per_seed.to_csv(T / "key_interaction_per_seed.csv")

# resumo geral
geral = per_model.groupby("system")[["dist_polar_N_O", "dist_sidechain_min", "gap"]].agg(
    ["median", "min", "max"]).round(2)
print("\nResumo geral (25 modelos por sistema):\n", geral.T)

# fracao de modelos com N...O <= 3.5 A
frac = (per_model.assign(ok=per_model.dist_polar_N_O <= 3.5)
        .groupby("system")["ok"].mean().round(2))
print("\nFracao de modelos com N...O <= 3,5 A:\n", frac)