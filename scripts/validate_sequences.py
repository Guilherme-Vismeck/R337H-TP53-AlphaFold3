"""Baixa TP53 (UniProt P04637), extrai o dominio 325-356 e gera WT e R337H."""
from pathlib import Path
import requests

START, END = 325, 356
UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/P04637.fasta"
OUT_DIR = Path("data/sequences")


def fetch_canonical_sequence() -> str:
    r = requests.get(UNIPROT_URL, timeout=30)
    r.raise_for_status()
    lines = r.text.strip().splitlines()
    return "".join(lines[1:])  # pula o cabecalho


def main():
    full = fetch_canonical_sequence()
    print(f"Comprimento da TP53 completa: {len(full)} (esperado: 393)")
    assert len(full) == 393

    # Numeracao biologica comeca em 1, Python em 0
    domain_wt = full[START - 1 : END]
    assert len(domain_wt) == END - START + 1

    pos_337 = 337 - START  # indice dentro do fragmento
    pos_352 = 352 - START
    assert domain_wt[pos_337] == "R", "Posicao 337 deveria ser Arg"
    assert domain_wt[pos_352] == "D", "Posicao 352 deveria ser Asp"

    domain_mut = domain_wt[:pos_337] + "H" + domain_wt[pos_337 + 1 :]

    # Somente a posicao 337 pode diferir
    diffs = [i for i, (a, b) in enumerate(zip(domain_wt, domain_mut)) if a != b]
    assert diffs == [pos_337]

    print(f"WT    ({START}-{END}): {domain_wt}")
    print(f"R337H ({START}-{END}): {domain_mut}")
    print(f"Arg337 -> indice {pos_337 + 1} no fragmento (1-based)")
    print(f"Asp352 -> indice {pos_352 + 1} no fragmento (1-based)")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "TP53_domain_325_356_WT.fasta").write_text(
        f">TP53_{START}_{END}_WT\n{domain_wt}\n"
    )
    (OUT_DIR / "TP53_domain_325_356_R337H.fasta").write_text(
        f">TP53_{START}_{END}_R337H\n{domain_mut}\n"
    )
    print("FASTAs salvos em data/sequences/")


if __name__ == "__main__":
    main()