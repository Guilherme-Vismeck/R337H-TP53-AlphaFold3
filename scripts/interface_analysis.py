"""Contatos entre cadeias e interacao 337-352 (1C26, WT e R337H)."""
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

from structure_utils import (N_CHAINS, N_MODELS, N_RES, REF, START, SYSTEMS,
                             best_alignment, ca_array, load_chains, model_path)

OUT = Path("results/tables")
OUT.mkdir(parents=True, exist_ok=True)

CUT_MAIN, CUT_STRICT = 5.0, 4.0  # definidos ANTES de olhar os resultados
LABELS = "ABCD"  # ordem das cadeias da referencia 1C26
BACKBONE = {"N", "CA", "C", "O", "OXT"}
POLAR_337 = {"ARG": {"NE", "NH1", "NH2"}, "HIS": {"ND1", "NE2"}}
POLAR_352 = {"OD1", "OD2"}


def min_dist(a, b):
    return np.sqrt(((a[:, None, :] - b[None, :, :]) ** 2).sum(-1)).min()


def atom_data(chains):
    """[cadeia][residuo] -> (nomes_atomos, coords) so de atomos pesados."""
    out = []
    for _, residues in chains:
        chain_data = []
        for r in residues:
            atoms = [a for a in r if a.element != "H"]
            chain_data.append(([a.get_id() for a in atoms],
                               np.array([a.coord for a in atoms], dtype=float)))
        out.append(chain_data)
    return out


def pick(data, names=None, sidechain=False):
    atom_names, coords = data
    keep = [k for k, n in enumerate(atom_names)
            if (names is None or n in names) and not (sidechain and n in BACKBONE)]
    return coords[keep]


def analyze(system, model, chains):
    atoms = atom_data(chains)
    resname = [[r.get_resname() for r in res] for _, res in chains]
    contacts, key = [], []

    for ci, cj in combinations(range(N_CHAINS), 2):
        for ri in range(N_RES):
            for rj in range(N_RES):
                d = min_dist(atoms[ci][ri][1], atoms[cj][rj][1])
                if d < CUT_MAIN:
                    contacts.append({
                        "system": system, "model": model,
                        "chain_i": LABELS[ci], "chain_j": LABELS[cj],
                        "res_i": START + ri, "resname_i": resname[ci][ri],
                        "res_j": START + rj, "resname_j": resname[cj][rj],
                        "min_dist": round(d, 2), "contact_4A": d < CUT_STRICT,
                    })

    i337, i352 = 337 - START, 352 - START
    for ci in range(N_CHAINS):
        name337 = resname[ci][i337]
        sc337 = pick(atoms[ci][i337], sidechain=True)
        pol337 = pick(atoms[ci][i337], names=POLAR_337[name337])
        cands = []
        for cj in range(N_CHAINS):
            if cj == ci:
                continue
            sc352 = pick(atoms[cj][i352], sidechain=True)
            pol352 = pick(atoms[cj][i352], names=POLAR_352)
            cands.append((min_dist(pol337, pol352), min_dist(sc337, sc352), cj))
        d_pol, d_sc, cj = min(cands)  # parceira = Asp352 com atomo polar mais proximo
        key.append({
            "system": system, "model": model,
            "chain_337": LABELS[ci], "resname_337": name337,
            "partner_chain_352": LABELS[cj],
            "dist_polar_N_O": round(d_pol, 2),
            "dist_sidechain_min": round(d_sc, 2),
        })
    return contacts, key


entries = [("1C26", 0, REF)] + [(s, i, model_path(s, i))
                                for s in SYSTEMS for i in range(N_MODELS)]
ref_ca = ca_array(load_chains(REF))

all_contacts, all_key = [], []
for system, model, path in entries:
    chains = load_chains(path)
    if system != "1C26":  # reordena as cadeias para a ordem da referencia
        _, perm, _ = best_alignment(ca_array(chains), ref_ca)
        chains = [chains[p] for p in perm]
    c, k = analyze(system, model, chains)
    all_contacts += c
    all_key += k
    print(f"ok: {system} {model}")

df = pd.DataFrame(all_contacts)
key = pd.DataFrame(all_key)
df.to_csv(OUT / "interface_contacts.csv", index=False)
key.to_csv(OUT / "key_interaction_337_352.csv", index=False)

# --- Contatos por interface ---
pairs = [f"{a}-{b}" for a, b in combinations(LABELS, 2)]
df["interface"] = df.chain_i + "-" + df.chain_j


def count_table(d):
    t = (d.groupby(["system", "model", "interface"]).size()
         .unstack("interface", fill_value=0)
         .reindex(columns=pairs, fill_value=0))
    t["total"] = t.sum(axis=1)
    return t.groupby(level="system").median()


counts5, counts4 = count_table(df), count_table(df[df.contact_4A])
counts5.to_csv(OUT / "interface_counts_5A.csv")
counts4.to_csv(OUT / "interface_counts_4A.csv")
print(f"\nPares de residuos em contato (<{CUT_MAIN} A), mediana dos modelos:\n{counts5}")
print(f"\nPares de residuos em contato (<{CUT_STRICT} A), mediana dos modelos:\n{counts4}")

# --- Interacao 337-352 ---
print("\nDistancias 337-352 (por subunidade):")
print(key.to_string(index=False))
resumo = key.groupby(["system", "model"])[["dist_polar_N_O", "dist_sidechain_min"]].mean()
print("\nMedia das 4 subunidades, por modelo:")
print(resumo.round(2))

# --- Contatos consistentes entre modelos ---
sets = {k: set(zip(g.chain_i, g.chain_j, g.res_i, g.res_j))
        for k, g in df.groupby(["system", "model"])}
wt_core = set.intersection(*[sets[("WT", i)] for i in range(N_MODELS)])
mut_core = set.intersection(*[sets[("R337H", i)] for i in range(N_MODELS)])
wt_any = set.union(*[sets[("WT", i)] for i in range(N_MODELS)])
mut_any = set.union(*[sets[("R337H", i)] for i in range(N_MODELS)])
ref_set = sets[("1C26", 0)]

print(f"\nContatos presentes nos 5 modelos WT: {len(wt_core)}; nos 5 R337H: {len(mut_core)}")
print("Nos 5 WT e em NENHUM R337H:", sorted(wt_core - mut_any))
print("Nos 5 R337H e em NENHUM WT:", sorted(mut_core - wt_any))
print(f"\n1C26: {len(ref_set)} contatos; em comum com o nucleo WT: {len(ref_set & wt_core)}")
print("No 1C26 mas em nenhum WT:", sorted(ref_set - wt_any))
print("Nos 5 WT mas ausentes no 1C26:", sorted(wt_core - ref_set))