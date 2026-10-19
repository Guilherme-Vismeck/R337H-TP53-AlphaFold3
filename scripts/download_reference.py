"""Baixa a montagem biologica do PDB 1C26 e mostra as cadeias."""
from pathlib import Path
import requests
from Bio.PDB import MMCIFParser

OUT = Path("data/reference")
FILES = {
    "PDB_1C26_biological_assembly.cif": "https://files.rcsb.org/download/1C26-assembly1.cif",
    "PDB_1C26_asymmetric_unit.cif": "https://files.rcsb.org/download/1C26.cif",
}

OUT.mkdir(parents=True, exist_ok=True)
for name, url in FILES.items():
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    (OUT / name).write_text(r.text)
    print(f"baixado: {name}")

parser = MMCIFParser(QUIET=True)
for name in FILES:
    s = parser.get_structure(name, OUT / name)
    print(f"\n{name}")
    for chain in s[0]:
        res = [r for r in chain if r.id[0] == " "]
        print(f"  cadeia {chain.id}: {len(res)} residuos, de {res[0].id[1]} a {res[-1].id[1]}")