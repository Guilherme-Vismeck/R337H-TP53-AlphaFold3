"""Controle de qualidade dos modelos AF3 (WT e R337H)."""
import json
from pathlib import Path
import pandas as pd
from Bio.PDB import MMCIFParser
from Bio.SeqUtils import seq1

PRED = Path("data/predictions")
SEQ = Path("data/sequences")
OUT = Path("results/tables")
OUT.mkdir(parents=True, exist_ok=True)

START = 325  # o AF3 numera 1-32; somamos START-1 para voltar a numeracao da TP53
IDX_337, IDX_352 = 337 - START + 1, 352 - START + 1  # 13 e 28


def read_fasta(path):
    return "".join(path.read_text().splitlines()[1:])


expected = {
    "WT": read_fasta(SEQ / "TP53_domain_325_356_WT.fasta"),
    "R337H": read_fasta(SEQ / "TP53_domain_325_356_R337H.fasta"),
}

parser = MMCIFParser(QUIET=True)
rows = []

for system, exp_seq in expected.items():
    for i in range(5):
        cif = next((PRED / system).glob(f"*_model_{i}.cif"))
        summ = next((PRED / system).glob(f"*_summary_confidences_{i}.json"))
        conf = json.loads(summ.read_text())
        model = parser.get_structure(f"{system}_{i}", cif)[0]
        chains = list(model)

        seqs_ok = True
        plddt_337, plddt_352, plddt_mean = [], [], []
        res_337 = set()
        for ch in chains:
            residues = [r for r in ch if r.id[0] == " "]
            seq = "".join(seq1(r.get_resname()) for r in residues)
            seqs_ok &= seq == exp_seq
            res_337.add(residues[IDX_337 - 1].get_resname())
            # No AF3, o B-factor do mmCIF guarda o pLDDT
            for r in residues:
                plddt_mean.append(r["CA"].get_bfactor())
            plddt_337.append(residues[IDX_337 - 1]["CA"].get_bfactor())
            plddt_352.append(residues[IDX_352 - 1]["CA"].get_bfactor())

        rows.append({
            "system": system,
            "model": i,
            "n_chains": len(chains),
            "chain_ids": "".join(c.id for c in chains),
            "sequence_matches": seqs_ok,
            "residue_at_337": "/".join(sorted(res_337)),
            "ptm": conf.get("ptm"),
            "iptm": conf.get("iptm"),
            "ranking_score": conf.get("ranking_score"),
            "has_clash": conf.get("has_clash"),
            "fraction_disordered": conf.get("fraction_disordered"),
            "plddt_mean": round(sum(plddt_mean) / len(plddt_mean), 2),
            "plddt_res337_mean": round(sum(plddt_337) / 4, 2),
            "plddt_res352_mean": round(sum(plddt_352) / 4, 2),
        })

df = pd.DataFrame(rows)
df.to_csv(OUT / "structure_quality.csv", index=False)
print(df.to_string(index=False))