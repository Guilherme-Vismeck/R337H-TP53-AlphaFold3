"""Confere que a seed da pasta bate com modelSeeds do job_request.json."""
import json
from pathlib import Path

bad = 0
for req in Path("data/predictions").rglob("*job_request.json"):
    expected = req.parent.name.replace("seed_", "")
    seeds = json.loads(req.read_text())[0]["modelSeeds"]
    ok = seeds == [expected]
    bad += not ok
    print(("OK   " if ok else "ERRO ") + f"{req.parent.parent.name}/{req.parent.name}: modelSeeds={seeds}")
assert bad == 0, f"{bad} pasta(s) com seed inconsistente"