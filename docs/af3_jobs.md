# Registro dos jobs do AlphaFold Server

Todos: entidade *Protein*, 4 cópias, sequência de `data/sequences/`, `useStructureTemplate: true` (padrão do servidor). Sequências: WT `GEYFTLQIRGRERFEMFRELNEALELKDAQAG`; R337H `GEYFTLQIRGREHFEMFRELNEALELKDAQAG`.

| Sistema | Nome do job | Seed | Pasta |
|---|---|---|---|
| WT | p53_WT_325_356_x4 | 333 | WT/seed_333 |
| WT | p53_WT_325_356_x4_seed1 | 1 | WT/seed_1 |
| WT | p53_WT_325_356_x4_seed2 | 2 | WT/seed_2 |
| WT | p53_WT_325_356_x4_seed3 | 3 | WT/seed_3 |
| WT | p53_WT_325_356_x4_seed4 | 4 | WT/seed_4 |
| R337H | p53_R337H_325_356_x4 | 333 | R337H/seed_333 |
| R337H | p53_R337H_325_356_x4_seed1_v2 | 1 | R337H/seed_1 |
| R337H | p53_R337H_325_356_x4_seed2 | 2 | R337H/seed_2 |
| R337H | p53_R337H_325_356_x4_seed3 | 3 | R337H/seed_3 |
| R337H | p53_R337H_325_356_x4_seed4 | 4 | R337H/seed_4 |

Jobs enviados em 2026-10-19 (seed 333 e seeds 1–4; o `seed1_v2` foi reenviado depois).

**Nota:** o job `p53_R337H_325_356_x4_seed1` foi enviado por engano com a seed 2 (verificado em `modelSeeds` do `job_request.json`) e substituído pelo `..._seed1_v2`, com a seed 1. A pasta `R337H/seed_1` contém apenas o job v2. Os modelos da seed 2 saíram idênticos nos dois jobs, o que indica reprodutibilidade do servidor para a mesma seed e sequência.
