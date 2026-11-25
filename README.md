# Early-life infections and later autism and ADHD

Analysis code and aggregate data for:

> Ssentongo, A. E. & Ssentongo, P. Early-life cytomegalovirus disease but not common infections is associated with later autism.

## What this repository contains

This study used the TriNetX Global Collaborative Network, a federated network of deidentified
electronic health records. **Patient-level data cannot be shared.** They are held by TriNetX and
its participating health care organizations, and access requires a data use agreement with
TriNetX, LLC (https://trinetx.com).

What is shared here is everything needed to reproduce the published estimates from the aggregate
counts returned by the platform: the cohort-level counts and event counts, the derived statistics,
and the code that produced every number and figure in the paper.

```
code/              analysis and figure scripts (Python 3.13, SciPy 1.18)
data/              aggregate results returned by the platform, and derived statistics
figures/           figures as published
requirements.txt   pinned dependencies
run_all.sh         reproduces every statistic, table and figure
```

### `data/`

| File | Contents |
|---|---|
| `all_results_derived.csv` | One row per condition-outcome comparison (173 estimable models): matched pairs, numbers at risk, events, absolute risks, risk differences, hazard ratios with 95% CIs, P values, FDR q values, E-values, and reference-adjusted hazard ratios |
| `stats.json` | The same estimates plus the random-effects pooled results by condition group, in the structure the figure scripts read |
| `table1_cmv.json` | Baseline characteristics for the cytomegalovirus comparison, before and after propensity score matching |
| `sens_covid.json` | Same-era COVID-19 sensitivity analysis |
| `unrestricted.json` | Sensitivity analysis without the well-child requirement in exposed infants |
| `trinetx_results_*.txt` | Raw aggregate output as returned by the platform, before any derivation |

### `code/`

| File | Purpose |
|---|---|
| `stats.py` | Derives P values from each hazard ratio and its confidence interval, applies Benjamini-Hochberg FDR correction across the 144 estimable condition-outcome tests, computes absolute risks and risk differences, reference-adjusted hazard ratios, DerSimonian-Laird random-effects pooling, and E-values |
| `table1_cmv.py` | Builds Table 1 |
| `gen_tables.py` | Builds the supplementary tables |
| `fig1.py` to `fig4.py` | Main text figures |
| `efigs.py` | Supplementary figures |

## Reproducing the results

```bash
pip install -r requirements.txt
bash run_all.sh
```

`run_all.sh` assembles a `build/` directory with the code and data side by side, derives every
statistic, rebuilds Table 1, the supplementary tables and all figures, then copies the figures
back into `figures/`. The scripts read their inputs by bare filename, so run them through
`run_all.sh` rather than individually from `code/`.

Python 3.13 with SciPy 1.18 and Matplotlib. Propensity score matching, Cox proportional hazards
models and log-rank tests were run inside the TriNetX platform and are not reproducible from this
repository; everything downstream of the returned estimates is.

## Methods in brief

Infants with a first diagnosis of each of 22 conditions before age 2 years were compared with
infants who had a well-child visit and none of that comparison's exposure codes. Every exposed and
unexposed infant had a well-child visit before age 2 years, so both groups had similar opportunity
for a developmental concern to be recognised. Cohorts were propensity matched 1:1 on 27 covariates,
giving 1,736,035 matched pairs. Diaper dermatitis served as a benign reference exposure and arm
fracture as a negative-control outcome.

Counts below 11 are returned as 10 by the platform. Estimates that could not be released are
marked as not reported.

## Ethics

Data in the network are deidentified by expert determination under the HIPAA Privacy Rule,
section 164.514(b)(1). The analysis did not constitute human-subjects research, and institutional
review board approval and informed consent were not required.

## Funding

National Institute of Mental Health, National Institutes of Health (grant 11456721 to A.E.S.;
MIND Study).

## Licence

Code is released under the MIT Licence (`LICENSE`). Aggregate data files are released under
CC BY 4.0.

## Citation

Please cite the paper. A citable archive of this repository is deposited on Zenodo; the DOI is
shown in the badge above once minted.
