# LZ WS2024 observed SI limit

This is the numerical release associated with the published [LZ analysis](https://arxiv.org/abs/2410.17036v3), [HEPData 155182 version 2, table 1](https://doi.org/10.17182/hepdata.155182.v2/t1). The original YAML is distributed under CC0.

`observed-si-v2.json` selects the **`limit`** column: the observed, power-constrained 90% CL elastic isoscalar spin-independent limit, per nucleon. It does not select the unconstrained observed limit, expected median or discovery curve. Values are stored in cm²; the comparison converts using **1 pb = 10⁻³⁶ cm²**. The model's neutron SI cross section is used for the isoscalar comparison.

There are 26 masses spanning **9–10,000 GeV**. Representative checks are:

| Mass [GeV] | Observed limit [cm²] |
|---|---:|
| 9 | 9.797484060822151e-47 |
| 40 | 2.1816833824484827e-48 |
| 10,000 | 2.93003225261514e-46 |

The 40 GeV minimum agrees with the published approximately 2.2e-48 cm² result. Between nodes the pipeline interpolates in log(mass)–log(cross section). It never extrapolates; DD outside coverage is unassessed, while the point and other results are retained.

Reproduce the normalized table and verify the original-file checksum with:

```sh
python DM/tools/verify_lz2024.py
```

The YAML SHA256 is `cdcb151339937e9bb1e27fe858359e3c4eb5f14ae4ff82c8538eaba1a8d5a870`. The normalized-table SHA256 is included in every scan manifest. The 2026 high-mass approximation remains available only by explicitly passing `--dm-limit-table DM/data/lz2026/lz2026-figs7-highmass-approx.json`.
