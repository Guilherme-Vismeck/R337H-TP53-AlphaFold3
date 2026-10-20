"""Funcoes compartilhadas pelos scripts de analise estrutural."""
from itertools import permutations
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser

START = 325  # primeiro residuo do fragmento (numeracao da TP53)
N_RES = 32  # residuos 325..356
N_CHAINS = 4
N_MODELS = 5
SYSTEMS = ["WT", "R337H"]
PRED = Path("data/predictions")
REF = Path("data/reference/PDB_1C26_biological_assembly.cif")

_parser = MMCIFParser(QUIET=True)


def load_chains(path):
    """Le um mmCIF e retorna [(id_cadeia, [residuos]), ...] sem aguas/ligantes."""
    model = _parser.get_structure(Path(path).stem, str(path))[0]
    chains = [(ch.id, [r for r in ch if r.id[0] == " "]) for ch in model]
    assert len(chains) == N_CHAINS, f"{path}: {len(chains)} cadeias"
    for cid, res in chains:
        assert len(res) == N_RES, f"{path}: cadeia {cid} tem {len(res)} residuos"
    return chains


def model_path(system, i):
    return next((PRED / system).glob(f"*_model_{i}.cif"))


def ca_array(chains):
    """Array (4, 32, 3) com as coordenadas dos carbonos alfa."""
    return np.array([[r["CA"].coord for r in res] for _, res in chains], dtype=float)


def kabsch(P, Q):
    """Sobrepoe P em Q (minimos quadrados). Retorna (rmsd, P_sobreposto)."""
    pc, qc = P.mean(axis=0), Q.mean(axis=0)
    H = (P - pc).T @ (Q - qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))  # evita reflexao (imagem em espelho)
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    P_fit = (P - pc) @ R.T + qc
    rmsd = np.sqrt(((P_fit - Q) ** 2).sum(axis=1).mean())
    return rmsd, P_fit


def best_alignment(model_ca, target_ca):
    """Testa as 24 permutacoes de cadeias e devolve a de menor RMSD.

    Retorna (rmsd, perm, P_sobreposto) onde a cadeia perm[k] do modelo
    corresponde a cadeia k do alvo.
    """
    target_flat = target_ca.reshape(-1, 3)
    best = None
    for perm in permutations(range(N_CHAINS)):
        P = model_ca[list(perm)].reshape(-1, 3)
        rmsd, P_fit = kabsch(P, target_flat)
        if best is None or rmsd < best[0]:
            best = (rmsd, perm, P_fit.reshape(N_CHAINS, N_RES, 3))
    return best