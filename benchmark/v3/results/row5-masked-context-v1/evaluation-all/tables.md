# Row-5 Masked-Context Evaluation (`scope=all`)

Chosen context window (bp): `10000`

## 1. Contrasts (Overall and Per Claim Set)

| Contrast | Claim set | Cells | Improved | Harmed | Unchanged | Discordant | Critical k | p (improve) | Harm elsewhere | Verdict |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `C1_stage0_to_gaps` | **all** | 311 | 0 | 0 | 311 | 0 | - | 1 | 0 | no detectable effect |
| `C1_stage0_to_gaps` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C1_stage0_to_gaps` | `diamond` | 119 | 0 | 0 | 119 | 0 | - | 1 | 0 | no detectable effect |
| `C1_stage0_to_gaps` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C1_stage0_to_gaps` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C1_stage0_to_gaps` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C1_stage0_to_gaps` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C1_stage0_to_gaps` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C1b_gaps_to_gapsfix` | **all** | 311 | 69 | 2 | 240 | 71 | 43 | 1.083e-18 | 0 | systematic improvement |
| `C1b_gaps_to_gapsfix` | `clean` | 32 | 6 | 0 | 26 | 6 | 6 | 0.01562 | 0 | systematic improvement |
| `C1b_gaps_to_gapsfix` | `diamond` | 119 | 44 | 0 | 75 | 44 | 28 | 5.684e-14 | 0 | systematic improvement |
| `C1b_gaps_to_gapsfix` | `dpz` | 32 | 4 | 0 | 28 | 4 | - | 0.0625 | 0 | no detectable effect |
| `C1b_gaps_to_gapsfix` | `enzbert` | 32 | 2 | 0 | 30 | 2 | - | 0.25 | 0 | no detectable effect |
| `C1b_gaps_to_gapsfix` | `graphec` | 32 | 8 | 0 | 24 | 8 | 7 | 0.003906 | 0 | systematic improvement |
| `C1b_gaps_to_gapsfix` | `mapred` | 32 | 5 | 2 | 25 | 7 | 7 | 0.2266 | 0 | no detectable effect |
| `C1b_gaps_to_gapsfix` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C1c_stage0_to_gapsfix` | **all** | 311 | 69 | 2 | 240 | 71 | 43 | 1.083e-18 | 0 | systematic improvement |
| `C1c_stage0_to_gapsfix` | `clean` | 32 | 6 | 0 | 26 | 6 | 6 | 0.01562 | 0 | systematic improvement |
| `C1c_stage0_to_gapsfix` | `diamond` | 119 | 44 | 0 | 75 | 44 | 28 | 5.684e-14 | 0 | systematic improvement |
| `C1c_stage0_to_gapsfix` | `dpz` | 32 | 4 | 0 | 28 | 4 | - | 0.0625 | 0 | no detectable effect |
| `C1c_stage0_to_gapsfix` | `enzbert` | 32 | 2 | 0 | 30 | 2 | - | 0.25 | 0 | no detectable effect |
| `C1c_stage0_to_gapsfix` | `graphec` | 32 | 8 | 0 | 24 | 8 | 7 | 0.003906 | 0 | systematic improvement |
| `C1c_stage0_to_gapsfix` | `mapred` | 32 | 5 | 2 | 25 | 7 | 7 | 0.2266 | 0 | no detectable effect |
| `C1c_stage0_to_gapsfix` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal_w5000` | **all** | 311 | 1 | 13 | 297 | 14 | 11 | 0.9999 | 16 | harm |
| `C5_universal_w5000` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal_w5000` | `diamond` | 119 | 1 | 13 | 105 | 14 | 11 | 0.9999 | 16 | harm |
| `C5_universal_w5000` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal_w5000` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal_w5000` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal_w5000` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal_w5000` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_w5000` | **all** | 311 | 1 | 8 | 302 | 9 | 8 | 0.998 | 3 | harm |
| `C5_pathway_w5000` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_w5000` | `diamond` | 119 | 1 | 8 | 110 | 9 | 8 | 0.998 | 3 | harm |
| `C5_pathway_w5000` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_w5000` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_w5000` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_w5000` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_w5000` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_w5000` | **all** | 311 | 1 | 11 | 299 | 12 | 10 | 0.9998 | 14 | harm |
| `C5_organism_w5000` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_w5000` | `diamond` | 119 | 1 | 11 | 107 | 12 | 10 | 0.9998 | 14 | harm |
| `C5_organism_w5000` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_w5000` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_w5000` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_w5000` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_w5000` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5total_universal_w5000` | **all** | 311 | 57 | 2 | 252 | 59 | 37 | 3.072e-15 | 1 | systematic improvement |
| `C5total_universal_w5000` | `clean` | 32 | 6 | 0 | 26 | 6 | 6 | 0.01562 | 0 | systematic improvement |
| `C5total_universal_w5000` | `diamond` | 119 | 32 | 0 | 87 | 32 | 22 | 2.328e-10 | 1 | systematic improvement |
| `C5total_universal_w5000` | `dpz` | 32 | 4 | 0 | 28 | 4 | - | 0.0625 | 0 | no detectable effect |
| `C5total_universal_w5000` | `enzbert` | 32 | 2 | 0 | 30 | 2 | - | 0.25 | 0 | no detectable effect |
| `C5total_universal_w5000` | `graphec` | 32 | 8 | 0 | 24 | 8 | 7 | 0.003906 | 0 | systematic improvement |
| `C5total_universal_w5000` | `mapred` | 32 | 5 | 2 | 25 | 7 | 7 | 0.2266 | 0 | no detectable effect |
| `C5total_universal_w5000` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w5000` | **all** | 311 | 5 | 0 | 306 | 5 | 5 | 0.03125 | 1 | systematic improvement |
| `C5_pathway_vs_universal_w5000` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w5000` | `diamond` | 119 | 5 | 0 | 114 | 5 | 5 | 0.03125 | 1 | systematic improvement |
| `C5_pathway_vs_universal_w5000` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w5000` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w5000` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w5000` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w5000` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w5000` | **all** | 311 | 2 | 0 | 309 | 2 | - | 0.25 | 1 | no detectable effect |
| `C5_organism_vs_universal_w5000` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w5000` | `diamond` | 119 | 2 | 0 | 117 | 2 | - | 0.25 | 1 | no detectable effect |
| `C5_organism_vs_universal_w5000` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w5000` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w5000` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w5000` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w5000` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway_w5000` | **all** | 311 | 0 | 3 | 308 | 3 | - | 1 | 11 | no detectable effect |
| `C5_organism_vs_pathway_w5000` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway_w5000` | `diamond` | 119 | 0 | 3 | 116 | 3 | - | 1 | 11 | no detectable effect |
| `C5_organism_vs_pathway_w5000` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway_w5000` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway_w5000` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway_w5000` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway_w5000` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal_w10000` | **all** | 311 | 1 | 13 | 297 | 14 | 11 | 0.9999 | 16 | harm |
| `C5_universal_w10000` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal_w10000` | `diamond` | 119 | 1 | 13 | 105 | 14 | 11 | 0.9999 | 16 | harm |
| `C5_universal_w10000` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal_w10000` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal_w10000` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal_w10000` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal_w10000` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_w10000` | **all** | 311 | 1 | 8 | 302 | 9 | 8 | 0.998 | 3 | harm |
| `C5_pathway_w10000` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_w10000` | `diamond` | 119 | 1 | 8 | 110 | 9 | 8 | 0.998 | 3 | harm |
| `C5_pathway_w10000` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_w10000` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_w10000` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_w10000` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_w10000` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_w10000` | **all** | 311 | 1 | 11 | 299 | 12 | 10 | 0.9998 | 14 | harm |
| `C5_organism_w10000` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_w10000` | `diamond` | 119 | 1 | 11 | 107 | 12 | 10 | 0.9998 | 14 | harm |
| `C5_organism_w10000` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_w10000` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_w10000` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_w10000` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_w10000` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5total_universal_w10000` | **all** | 311 | 57 | 2 | 252 | 59 | 37 | 3.072e-15 | 1 | systematic improvement |
| `C5total_universal_w10000` | `clean` | 32 | 6 | 0 | 26 | 6 | 6 | 0.01562 | 0 | systematic improvement |
| `C5total_universal_w10000` | `diamond` | 119 | 32 | 0 | 87 | 32 | 22 | 2.328e-10 | 1 | systematic improvement |
| `C5total_universal_w10000` | `dpz` | 32 | 4 | 0 | 28 | 4 | - | 0.0625 | 0 | no detectable effect |
| `C5total_universal_w10000` | `enzbert` | 32 | 2 | 0 | 30 | 2 | - | 0.25 | 0 | no detectable effect |
| `C5total_universal_w10000` | `graphec` | 32 | 8 | 0 | 24 | 8 | 7 | 0.003906 | 0 | systematic improvement |
| `C5total_universal_w10000` | `mapred` | 32 | 5 | 2 | 25 | 7 | 7 | 0.2266 | 0 | no detectable effect |
| `C5total_universal_w10000` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w10000` | **all** | 311 | 5 | 0 | 306 | 5 | 5 | 0.03125 | 1 | systematic improvement |
| `C5_pathway_vs_universal_w10000` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w10000` | `diamond` | 119 | 5 | 0 | 114 | 5 | 5 | 0.03125 | 1 | systematic improvement |
| `C5_pathway_vs_universal_w10000` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w10000` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w10000` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w10000` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w10000` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w10000` | **all** | 311 | 2 | 0 | 309 | 2 | - | 0.25 | 1 | no detectable effect |
| `C5_organism_vs_universal_w10000` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w10000` | `diamond` | 119 | 2 | 0 | 117 | 2 | - | 0.25 | 1 | no detectable effect |
| `C5_organism_vs_universal_w10000` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w10000` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w10000` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w10000` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w10000` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway_w10000` | **all** | 311 | 0 | 3 | 308 | 3 | - | 1 | 11 | no detectable effect |
| `C5_organism_vs_pathway_w10000` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway_w10000` | `diamond` | 119 | 0 | 3 | 116 | 3 | - | 1 | 11 | no detectable effect |
| `C5_organism_vs_pathway_w10000` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway_w10000` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway_w10000` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway_w10000` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway_w10000` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal_w20000` | **all** | 311 | 1 | 11 | 299 | 12 | 10 | 0.9998 | 16 | harm |
| `C5_universal_w20000` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal_w20000` | `diamond` | 119 | 1 | 11 | 107 | 12 | 10 | 0.9998 | 16 | harm |
| `C5_universal_w20000` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal_w20000` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal_w20000` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal_w20000` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal_w20000` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_w20000` | **all** | 311 | 1 | 7 | 303 | 8 | 7 | 0.9961 | 3 | harm |
| `C5_pathway_w20000` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_w20000` | `diamond` | 119 | 1 | 7 | 111 | 8 | 7 | 0.9961 | 3 | harm |
| `C5_pathway_w20000` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_w20000` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_w20000` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_w20000` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_w20000` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_w20000` | **all** | 311 | 1 | 10 | 300 | 11 | 9 | 0.9995 | 14 | harm |
| `C5_organism_w20000` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_w20000` | `diamond` | 119 | 1 | 10 | 108 | 11 | 9 | 0.9995 | 14 | harm |
| `C5_organism_w20000` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_w20000` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_w20000` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_w20000` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_w20000` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5total_universal_w20000` | **all** | 311 | 59 | 2 | 250 | 61 | 38 | 8.205e-16 | 1 | systematic improvement |
| `C5total_universal_w20000` | `clean` | 32 | 6 | 0 | 26 | 6 | 6 | 0.01562 | 0 | systematic improvement |
| `C5total_universal_w20000` | `diamond` | 119 | 34 | 0 | 85 | 34 | 23 | 5.821e-11 | 1 | systematic improvement |
| `C5total_universal_w20000` | `dpz` | 32 | 4 | 0 | 28 | 4 | - | 0.0625 | 0 | no detectable effect |
| `C5total_universal_w20000` | `enzbert` | 32 | 2 | 0 | 30 | 2 | - | 0.25 | 0 | no detectable effect |
| `C5total_universal_w20000` | `graphec` | 32 | 8 | 0 | 24 | 8 | 7 | 0.003906 | 0 | systematic improvement |
| `C5total_universal_w20000` | `mapred` | 32 | 5 | 2 | 25 | 7 | 7 | 0.2266 | 0 | no detectable effect |
| `C5total_universal_w20000` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w20000` | **all** | 311 | 4 | 0 | 307 | 4 | - | 0.0625 | 1 | no detectable effect |
| `C5_pathway_vs_universal_w20000` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w20000` | `diamond` | 119 | 4 | 0 | 115 | 4 | - | 0.0625 | 1 | no detectable effect |
| `C5_pathway_vs_universal_w20000` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w20000` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w20000` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w20000` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w20000` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w20000` | **all** | 311 | 2 | 1 | 308 | 3 | - | 0.5 | 1 | no detectable effect |
| `C5_organism_vs_universal_w20000` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w20000` | `diamond` | 119 | 2 | 1 | 116 | 3 | - | 0.5 | 1 | no detectable effect |
| `C5_organism_vs_universal_w20000` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w20000` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w20000` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w20000` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w20000` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway_w20000` | **all** | 311 | 0 | 3 | 308 | 3 | - | 1 | 11 | no detectable effect |
| `C5_organism_vs_pathway_w20000` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway_w20000` | `diamond` | 119 | 0 | 3 | 116 | 3 | - | 1 | 11 | no detectable effect |
| `C5_organism_vs_pathway_w20000` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway_w20000` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway_w20000` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway_w20000` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway_w20000` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal` | **all** | 311 | 1 | 13 | 297 | 14 | 11 | 0.9999 | 16 | harm |
| `C5_universal` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal` | `diamond` | 119 | 1 | 13 | 105 | 14 | 11 | 0.9999 | 16 | harm |
| `C5_universal` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway` | **all** | 311 | 1 | 8 | 302 | 9 | 8 | 0.998 | 3 | harm |
| `C5_pathway` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway` | `diamond` | 119 | 1 | 8 | 110 | 9 | 8 | 0.998 | 3 | harm |
| `C5_pathway` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism` | **all** | 311 | 1 | 11 | 299 | 12 | 10 | 0.9998 | 14 | harm |
| `C5_organism` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism` | `diamond` | 119 | 1 | 11 | 107 | 12 | 10 | 0.9998 | 14 | harm |
| `C5_organism` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal` | **all** | 311 | 5 | 0 | 306 | 5 | 5 | 0.03125 | 1 | systematic improvement |
| `C5_pathway_vs_universal` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal` | `diamond` | 119 | 5 | 0 | 114 | 5 | 5 | 0.03125 | 1 | systematic improvement |
| `C5_pathway_vs_universal` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal` | **all** | 311 | 2 | 0 | 309 | 2 | - | 0.25 | 1 | no detectable effect |
| `C5_organism_vs_universal` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal` | `diamond` | 119 | 2 | 0 | 117 | 2 | - | 0.25 | 1 | no detectable effect |
| `C5_organism_vs_universal` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway` | **all** | 311 | 0 | 3 | 308 | 3 | - | 1 | 11 | no detectable effect |
| `C5_organism_vs_pathway` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway` | `diamond` | 119 | 0 | 3 | 116 | 3 | - | 1 | 11 | no detectable effect |
| `C5_organism_vs_pathway` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_gaps` | **all** | 311 | 0 | 0 | 311 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_gaps` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_gaps` | `diamond` | 119 | 0 | 0 | 119 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_gaps` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_gaps` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_gaps` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_gaps` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_gaps` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_gapsfix` | **all** | 311 | 11 | 0 | 300 | 11 | 9 | 0.0004883 | 0 | systematic improvement |
| `S1_sham_gaps_to_sham_gapsfix` | `clean` | 32 | 3 | 0 | 29 | 3 | - | 0.125 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_gapsfix` | `diamond` | 119 | 3 | 0 | 116 | 3 | - | 0.125 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_gapsfix` | `dpz` | 32 | 3 | 0 | 29 | 3 | - | 0.125 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_gapsfix` | `enzbert` | 32 | 1 | 0 | 31 | 1 | - | 0.5 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_gapsfix` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_gapsfix` | `mapred` | 32 | 1 | 0 | 31 | 1 | - | 0.5 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_gapsfix` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w5000` | **all** | 311 | 0 | 2 | 309 | 2 | - | 1 | 13 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w5000` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w5000` | `diamond` | 119 | 0 | 2 | 117 | 2 | - | 1 | 13 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w5000` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w5000` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w5000` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w5000` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w5000` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w5000` | **all** | 311 | 9 | 0 | 302 | 9 | 8 | 0.001953 | 0 | systematic improvement |
| `S1_stage0_to_sham_ctx_universal_w5000` | `clean` | 32 | 3 | 0 | 29 | 3 | - | 0.125 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w5000` | `diamond` | 119 | 1 | 0 | 118 | 1 | - | 0.5 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w5000` | `dpz` | 32 | 3 | 0 | 29 | 3 | - | 0.125 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w5000` | `enzbert` | 32 | 1 | 0 | 31 | 1 | - | 0.5 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w5000` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w5000` | `mapred` | 32 | 1 | 0 | 31 | 1 | - | 0.5 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w5000` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w10000` | **all** | 311 | 0 | 2 | 309 | 2 | - | 1 | 13 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w10000` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w10000` | `diamond` | 119 | 0 | 2 | 117 | 2 | - | 1 | 13 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w10000` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w10000` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w10000` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w10000` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w10000` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w10000` | **all** | 311 | 9 | 0 | 302 | 9 | 8 | 0.001953 | 0 | systematic improvement |
| `S1_stage0_to_sham_ctx_universal_w10000` | `clean` | 32 | 3 | 0 | 29 | 3 | - | 0.125 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w10000` | `diamond` | 119 | 1 | 0 | 118 | 1 | - | 0.5 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w10000` | `dpz` | 32 | 3 | 0 | 29 | 3 | - | 0.125 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w10000` | `enzbert` | 32 | 1 | 0 | 31 | 1 | - | 0.5 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w10000` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w10000` | `mapred` | 32 | 1 | 0 | 31 | 1 | - | 0.5 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w10000` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w20000` | **all** | 311 | 0 | 2 | 309 | 2 | - | 1 | 13 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w20000` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w20000` | `diamond` | 119 | 0 | 2 | 117 | 2 | - | 1 | 13 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w20000` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w20000` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w20000` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w20000` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w20000` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w20000` | **all** | 311 | 9 | 0 | 302 | 9 | 8 | 0.001953 | 0 | systematic improvement |
| `S1_stage0_to_sham_ctx_universal_w20000` | `clean` | 32 | 3 | 0 | 29 | 3 | - | 0.125 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w20000` | `diamond` | 119 | 1 | 0 | 118 | 1 | - | 0.5 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w20000` | `dpz` | 32 | 3 | 0 | 29 | 3 | - | 0.125 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w20000` | `enzbert` | 32 | 1 | 0 | 31 | 1 | - | 0.5 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w20000` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w20000` | `mapred` | 32 | 1 | 0 | 31 | 1 | - | 0.5 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w20000` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal` | **all** | 311 | 0 | 2 | 309 | 2 | - | 1 | 13 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal` | `diamond` | 119 | 0 | 2 | 117 | 2 | - | 1 | 13 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `U1_unmasked_gaps_vs_gaps` | **all** | 311 | 0 | 0 | 311 | 0 | - | 1 | 0 | no detectable effect |
| `U1_unmasked_gaps_vs_gaps` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `U1_unmasked_gaps_vs_gaps` | `diamond` | 119 | 0 | 0 | 119 | 0 | - | 1 | 0 | no detectable effect |
| `U1_unmasked_gaps_vs_gaps` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `U1_unmasked_gaps_vs_gaps` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `U1_unmasked_gaps_vs_gaps` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `U1_unmasked_gaps_vs_gaps` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `U1_unmasked_gaps_vs_gaps` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |

## 2. Harm Elsewhere and Cost

| Contrast | Improved | Harmed | Harm elsewhere | Unconfirmed calls created | Delta calls >= 0.5 | Changed scores | Bootstrap 95% CI | Wilcoxon p (1-sided) |
|---|---:|---:|---:|---:|---:|---:|---|---:|
| `C1_stage0_to_gaps` | 0 | 0 | 0 | 0 | 9268 | 226405 | [0.0, 0.0] | 1 |
| `C1b_gaps_to_gapsfix` | 69 | 2 | 0 | 28 | 2957 | 11298 | [46.0, 90.0] | 0.004816 |
| `C1c_stage0_to_gapsfix` | 69 | 2 | 0 | 28 | 12225 | 229428 | [46.0, 90.0] | 0.004816 |
| `C5_universal_w5000` | 1 | 13 | 16 | 0 | -1020 | 1979 | [-23.0, -3.0] | 0.1587 |
| `C5_pathway_w5000` | 1 | 8 | 3 | 0 | -977 | 1890 | [-17.0, 1.0] | 0.1587 |
| `C5_organism_w5000` | 1 | 11 | 14 | 0 | -674 | 1571 | [-20.0, -2.0] | 0.1587 |
| `C5total_universal_w5000` | 57 | 2 | 1 | 25 | 1937 | 12493 | [38.0, 73.0] | 0.003566 |
| `C5_pathway_vs_universal_w5000` | 5 | 0 | 1 | 0 | 43 | 652 | [1.0, 9.0] | 1 |
| `C5_organism_vs_universal_w5000` | 2 | 0 | 1 | 0 | 346 | 681 | [0.0, 5.0] | 1 |
| `C5_organism_vs_pathway_w5000` | 0 | 3 | 11 | 0 | 303 | 825 | [-6.0, 0.0] | 1 |
| `C5_universal_w10000` | 1 | 13 | 16 | 0 | -980 | 1969 | [-23.0, -3.0] | 0.1587 |
| `C5_pathway_w10000` | 1 | 8 | 3 | 0 | -943 | 1891 | [-17.0, 1.0] | 0.1587 |
| `C5_organism_w10000` | 1 | 11 | 14 | 0 | -650 | 1560 | [-20.0, -2.0] | 0.1587 |
| `C5total_universal_w10000` | 57 | 2 | 1 | 25 | 1977 | 12483 | [38.0, 73.0] | 0.003566 |
| `C5_pathway_vs_universal_w10000` | 5 | 0 | 1 | 0 | 37 | 692 | [1.0, 9.0] | 1 |
| `C5_organism_vs_universal_w10000` | 2 | 0 | 1 | 0 | 330 | 723 | [0.0, 5.0] | 1 |
| `C5_organism_vs_pathway_w10000` | 0 | 3 | 11 | 0 | 293 | 850 | [-6.0, 0.0] | 1 |
| `C5_universal_w20000` | 1 | 11 | 16 | 0 | -965 | 1997 | [-21.0, -1.0] | 0.1587 |
| `C5_pathway_w20000` | 1 | 7 | 3 | 0 | -941 | 1899 | [-16.0, 1.0] | 0.1587 |
| `C5_organism_w20000` | 1 | 10 | 14 | 0 | -643 | 1596 | [-19.0, -1.0] | 0.1587 |
| `C5total_universal_w20000` | 59 | 2 | 1 | 25 | 1992 | 12503 | [39.0, 75.0] | 0.003566 |
| `C5_pathway_vs_universal_w20000` | 4 | 0 | 1 | 0 | 24 | 887 | [1.0, 8.0] | 1 |
| `C5_organism_vs_universal_w20000` | 2 | 1 | 1 | 0 | 322 | 820 | [-2.0, 4.0] | 1 |
| `C5_organism_vs_pathway_w20000` | 0 | 3 | 11 | 0 | 298 | 930 | [-6.0, 0.0] | 1 |
| `C5_universal` | 1 | 13 | 16 | 0 | -980 | 1969 | [-23.0, -3.0] | 0.1587 |
| `C5_pathway` | 1 | 8 | 3 | 0 | -943 | 1891 | [-17.0, 1.0] | 0.1587 |
| `C5_organism` | 1 | 11 | 14 | 0 | -650 | 1560 | [-20.0, -2.0] | 0.1587 |
| `C5_pathway_vs_universal` | 5 | 0 | 1 | 0 | 37 | 692 | [1.0, 9.0] | 1 |
| `C5_organism_vs_universal` | 2 | 0 | 1 | 0 | 330 | 723 | [0.0, 5.0] | 1 |
| `C5_organism_vs_pathway` | 0 | 3 | 11 | 0 | 293 | 850 | [-6.0, 0.0] | 1 |
| `S1_stage0_to_sham_gaps` | 0 | 0 | 0 | 0 | 9250 | 226387 | [0.0, 0.0] | 1 |
| `S1_sham_gaps_to_sham_gapsfix` | 11 | 0 | 0 | 20 | 2900 | 11226 | [2.0, 22.0] | 0.05282 |
| `S1_sham_gaps_to_sham_ctx_universal_w5000` | 0 | 2 | 13 | 0 | -1005 | 1938 | [-5.0, 0.0] | 1 |
| `S1_stage0_to_sham_ctx_universal_w5000` | 9 | 0 | 0 | 17 | 11145 | 229368 | [2.0, 19.0] | 0.05282 |
| `S1_sham_gaps_to_sham_ctx_universal_w10000` | 0 | 2 | 13 | 0 | -966 | 1926 | [-5.0, 0.0] | 1 |
| `S1_stage0_to_sham_ctx_universal_w10000` | 9 | 0 | 0 | 17 | 11184 | 229368 | [2.0, 19.0] | 0.05282 |
| `S1_sham_gaps_to_sham_ctx_universal_w20000` | 0 | 2 | 13 | 0 | -953 | 1955 | [-5.0, 0.0] | 1 |
| `S1_stage0_to_sham_ctx_universal_w20000` | 9 | 0 | 0 | 17 | 11197 | 229368 | [2.0, 19.0] | 0.05282 |
| `S1_sham_gaps_to_sham_ctx_universal` | 0 | 2 | 13 | 0 | -966 | 1926 | [-5.0, 0.0] | 1 |
| `U1_unmasked_gaps_vs_gaps` | 0 | 0 | 0 | 0 | 51 | 89 | [0.0, 0.0] | 1 |

## 3. Strata Breakdown

| Contrast | Stratum axis | Stratum | Cells | Improved | Harmed | Unchanged | Discordant |
|---|---|---|---:|---:|---:|---:|---:|
| `C1_stage0_to_gaps` | `domain` | `archaea` | 60 | 0 | 0 | 60 | 0 |
| `C1_stage0_to_gaps` | `domain` | `bacteria` | 251 | 0 | 0 | 251 | 0 |
| `C1_stage0_to_gaps` | `gap_kind` | `already_gap` | 60 | 0 | 0 | 60 | 0 |
| `C1_stage0_to_gaps` | `gap_kind` | `maskable_gap` | 251 | 0 | 0 | 251 | 0 |
| `C1_stage0_to_gaps` | `pathway_robust` | `marginal` | 142 | 0 | 0 | 142 | 0 |
| `C1_stage0_to_gaps` | `pathway_robust` | `robust` | 169 | 0 | 0 | 169 | 0 |
| `C1b_gaps_to_gapsfix` | `domain` | `archaea` | 60 | 23 | 0 | 37 | 23 |
| `C1b_gaps_to_gapsfix` | `domain` | `bacteria` | 251 | 46 | 2 | 203 | 48 |
| `C1b_gaps_to_gapsfix` | `gap_kind` | `already_gap` | 60 | 11 | 0 | 49 | 11 |
| `C1b_gaps_to_gapsfix` | `gap_kind` | `maskable_gap` | 251 | 58 | 2 | 191 | 60 |
| `C1b_gaps_to_gapsfix` | `pathway_robust` | `marginal` | 142 | 31 | 1 | 110 | 32 |
| `C1b_gaps_to_gapsfix` | `pathway_robust` | `robust` | 169 | 38 | 1 | 130 | 39 |
| `C1c_stage0_to_gapsfix` | `domain` | `archaea` | 60 | 23 | 0 | 37 | 23 |
| `C1c_stage0_to_gapsfix` | `domain` | `bacteria` | 251 | 46 | 2 | 203 | 48 |
| `C1c_stage0_to_gapsfix` | `gap_kind` | `already_gap` | 60 | 11 | 0 | 49 | 11 |
| `C1c_stage0_to_gapsfix` | `gap_kind` | `maskable_gap` | 251 | 58 | 2 | 191 | 60 |
| `C1c_stage0_to_gapsfix` | `pathway_robust` | `marginal` | 142 | 31 | 1 | 110 | 32 |
| `C1c_stage0_to_gapsfix` | `pathway_robust` | `robust` | 169 | 38 | 1 | 130 | 39 |
| `C5_universal_w5000` | `domain` | `archaea` | 60 | 0 | 8 | 52 | 8 |
| `C5_universal_w5000` | `domain` | `bacteria` | 251 | 1 | 5 | 245 | 6 |
| `C5_universal_w5000` | `gap_kind` | `already_gap` | 60 | 0 | 2 | 58 | 2 |
| `C5_universal_w5000` | `gap_kind` | `maskable_gap` | 251 | 1 | 11 | 239 | 12 |
| `C5_universal_w5000` | `pathway_robust` | `marginal` | 142 | 0 | 4 | 138 | 4 |
| `C5_universal_w5000` | `pathway_robust` | `robust` | 169 | 1 | 9 | 159 | 10 |
| `C5_pathway_w5000` | `domain` | `archaea` | 60 | 0 | 6 | 54 | 6 |
| `C5_pathway_w5000` | `domain` | `bacteria` | 251 | 1 | 2 | 248 | 3 |
| `C5_pathway_w5000` | `gap_kind` | `already_gap` | 60 | 0 | 1 | 59 | 1 |
| `C5_pathway_w5000` | `gap_kind` | `maskable_gap` | 251 | 1 | 7 | 243 | 8 |
| `C5_pathway_w5000` | `pathway_robust` | `marginal` | 142 | 0 | 1 | 141 | 1 |
| `C5_pathway_w5000` | `pathway_robust` | `robust` | 169 | 1 | 7 | 161 | 8 |
| `C5_organism_w5000` | `domain` | `archaea` | 60 | 0 | 6 | 54 | 6 |
| `C5_organism_w5000` | `domain` | `bacteria` | 251 | 1 | 5 | 245 | 6 |
| `C5_organism_w5000` | `gap_kind` | `already_gap` | 60 | 0 | 2 | 58 | 2 |
| `C5_organism_w5000` | `gap_kind` | `maskable_gap` | 251 | 1 | 9 | 241 | 10 |
| `C5_organism_w5000` | `pathway_robust` | `marginal` | 142 | 0 | 4 | 138 | 4 |
| `C5_organism_w5000` | `pathway_robust` | `robust` | 169 | 1 | 7 | 161 | 8 |
| `C5total_universal_w5000` | `domain` | `archaea` | 60 | 16 | 0 | 44 | 16 |
| `C5total_universal_w5000` | `domain` | `bacteria` | 251 | 41 | 2 | 208 | 43 |
| `C5total_universal_w5000` | `gap_kind` | `already_gap` | 60 | 9 | 0 | 51 | 9 |
| `C5total_universal_w5000` | `gap_kind` | `maskable_gap` | 251 | 48 | 2 | 201 | 50 |
| `C5total_universal_w5000` | `pathway_robust` | `marginal` | 142 | 27 | 1 | 114 | 28 |
| `C5total_universal_w5000` | `pathway_robust` | `robust` | 169 | 30 | 1 | 138 | 31 |
| `C5_pathway_vs_universal_w5000` | `domain` | `archaea` | 60 | 2 | 0 | 58 | 2 |
| `C5_pathway_vs_universal_w5000` | `domain` | `bacteria` | 251 | 3 | 0 | 248 | 3 |
| `C5_pathway_vs_universal_w5000` | `gap_kind` | `already_gap` | 60 | 1 | 0 | 59 | 1 |
| `C5_pathway_vs_universal_w5000` | `gap_kind` | `maskable_gap` | 251 | 4 | 0 | 247 | 4 |
| `C5_pathway_vs_universal_w5000` | `pathway_robust` | `marginal` | 142 | 3 | 0 | 139 | 3 |
| `C5_pathway_vs_universal_w5000` | `pathway_robust` | `robust` | 169 | 2 | 0 | 167 | 2 |
| `C5_organism_vs_universal_w5000` | `domain` | `archaea` | 60 | 2 | 0 | 58 | 2 |
| `C5_organism_vs_universal_w5000` | `domain` | `bacteria` | 251 | 0 | 0 | 251 | 0 |
| `C5_organism_vs_universal_w5000` | `gap_kind` | `already_gap` | 60 | 0 | 0 | 60 | 0 |
| `C5_organism_vs_universal_w5000` | `gap_kind` | `maskable_gap` | 251 | 2 | 0 | 249 | 2 |
| `C5_organism_vs_universal_w5000` | `pathway_robust` | `marginal` | 142 | 0 | 0 | 142 | 0 |
| `C5_organism_vs_universal_w5000` | `pathway_robust` | `robust` | 169 | 2 | 0 | 167 | 2 |
| `C5_organism_vs_pathway_w5000` | `domain` | `archaea` | 60 | 0 | 0 | 60 | 0 |
| `C5_organism_vs_pathway_w5000` | `domain` | `bacteria` | 251 | 0 | 3 | 248 | 3 |
| `C5_organism_vs_pathway_w5000` | `gap_kind` | `already_gap` | 60 | 0 | 1 | 59 | 1 |
| `C5_organism_vs_pathway_w5000` | `gap_kind` | `maskable_gap` | 251 | 0 | 2 | 249 | 2 |
| `C5_organism_vs_pathway_w5000` | `pathway_robust` | `marginal` | 142 | 0 | 3 | 139 | 3 |
| `C5_organism_vs_pathway_w5000` | `pathway_robust` | `robust` | 169 | 0 | 0 | 169 | 0 |
| `C5_universal_w10000` | `domain` | `archaea` | 60 | 0 | 8 | 52 | 8 |
| `C5_universal_w10000` | `domain` | `bacteria` | 251 | 1 | 5 | 245 | 6 |
| `C5_universal_w10000` | `gap_kind` | `already_gap` | 60 | 0 | 2 | 58 | 2 |
| `C5_universal_w10000` | `gap_kind` | `maskable_gap` | 251 | 1 | 11 | 239 | 12 |
| `C5_universal_w10000` | `pathway_robust` | `marginal` | 142 | 0 | 4 | 138 | 4 |
| `C5_universal_w10000` | `pathway_robust` | `robust` | 169 | 1 | 9 | 159 | 10 |
| `C5_pathway_w10000` | `domain` | `archaea` | 60 | 0 | 6 | 54 | 6 |
| `C5_pathway_w10000` | `domain` | `bacteria` | 251 | 1 | 2 | 248 | 3 |
| `C5_pathway_w10000` | `gap_kind` | `already_gap` | 60 | 0 | 1 | 59 | 1 |
| `C5_pathway_w10000` | `gap_kind` | `maskable_gap` | 251 | 1 | 7 | 243 | 8 |
| `C5_pathway_w10000` | `pathway_robust` | `marginal` | 142 | 0 | 1 | 141 | 1 |
| `C5_pathway_w10000` | `pathway_robust` | `robust` | 169 | 1 | 7 | 161 | 8 |
| `C5_organism_w10000` | `domain` | `archaea` | 60 | 0 | 6 | 54 | 6 |
| `C5_organism_w10000` | `domain` | `bacteria` | 251 | 1 | 5 | 245 | 6 |
| `C5_organism_w10000` | `gap_kind` | `already_gap` | 60 | 0 | 2 | 58 | 2 |
| `C5_organism_w10000` | `gap_kind` | `maskable_gap` | 251 | 1 | 9 | 241 | 10 |
| `C5_organism_w10000` | `pathway_robust` | `marginal` | 142 | 0 | 4 | 138 | 4 |
| `C5_organism_w10000` | `pathway_robust` | `robust` | 169 | 1 | 7 | 161 | 8 |
| `C5total_universal_w10000` | `domain` | `archaea` | 60 | 16 | 0 | 44 | 16 |
| `C5total_universal_w10000` | `domain` | `bacteria` | 251 | 41 | 2 | 208 | 43 |
| `C5total_universal_w10000` | `gap_kind` | `already_gap` | 60 | 9 | 0 | 51 | 9 |
| `C5total_universal_w10000` | `gap_kind` | `maskable_gap` | 251 | 48 | 2 | 201 | 50 |
| `C5total_universal_w10000` | `pathway_robust` | `marginal` | 142 | 27 | 1 | 114 | 28 |
| `C5total_universal_w10000` | `pathway_robust` | `robust` | 169 | 30 | 1 | 138 | 31 |
| `C5_pathway_vs_universal_w10000` | `domain` | `archaea` | 60 | 2 | 0 | 58 | 2 |
| `C5_pathway_vs_universal_w10000` | `domain` | `bacteria` | 251 | 3 | 0 | 248 | 3 |
| `C5_pathway_vs_universal_w10000` | `gap_kind` | `already_gap` | 60 | 1 | 0 | 59 | 1 |
| `C5_pathway_vs_universal_w10000` | `gap_kind` | `maskable_gap` | 251 | 4 | 0 | 247 | 4 |
| `C5_pathway_vs_universal_w10000` | `pathway_robust` | `marginal` | 142 | 3 | 0 | 139 | 3 |
| `C5_pathway_vs_universal_w10000` | `pathway_robust` | `robust` | 169 | 2 | 0 | 167 | 2 |
| `C5_organism_vs_universal_w10000` | `domain` | `archaea` | 60 | 2 | 0 | 58 | 2 |
| `C5_organism_vs_universal_w10000` | `domain` | `bacteria` | 251 | 0 | 0 | 251 | 0 |
| `C5_organism_vs_universal_w10000` | `gap_kind` | `already_gap` | 60 | 0 | 0 | 60 | 0 |
| `C5_organism_vs_universal_w10000` | `gap_kind` | `maskable_gap` | 251 | 2 | 0 | 249 | 2 |
| `C5_organism_vs_universal_w10000` | `pathway_robust` | `marginal` | 142 | 0 | 0 | 142 | 0 |
| `C5_organism_vs_universal_w10000` | `pathway_robust` | `robust` | 169 | 2 | 0 | 167 | 2 |
| `C5_organism_vs_pathway_w10000` | `domain` | `archaea` | 60 | 0 | 0 | 60 | 0 |
| `C5_organism_vs_pathway_w10000` | `domain` | `bacteria` | 251 | 0 | 3 | 248 | 3 |
| `C5_organism_vs_pathway_w10000` | `gap_kind` | `already_gap` | 60 | 0 | 1 | 59 | 1 |
| `C5_organism_vs_pathway_w10000` | `gap_kind` | `maskable_gap` | 251 | 0 | 2 | 249 | 2 |
| `C5_organism_vs_pathway_w10000` | `pathway_robust` | `marginal` | 142 | 0 | 3 | 139 | 3 |
| `C5_organism_vs_pathway_w10000` | `pathway_robust` | `robust` | 169 | 0 | 0 | 169 | 0 |
| `C5_universal_w20000` | `domain` | `archaea` | 60 | 0 | 8 | 52 | 8 |
| `C5_universal_w20000` | `domain` | `bacteria` | 251 | 1 | 3 | 247 | 4 |
| `C5_universal_w20000` | `gap_kind` | `already_gap` | 60 | 0 | 2 | 58 | 2 |
| `C5_universal_w20000` | `gap_kind` | `maskable_gap` | 251 | 1 | 9 | 241 | 10 |
| `C5_universal_w20000` | `pathway_robust` | `marginal` | 142 | 0 | 2 | 140 | 2 |
| `C5_universal_w20000` | `pathway_robust` | `robust` | 169 | 1 | 9 | 159 | 10 |
| `C5_pathway_w20000` | `domain` | `archaea` | 60 | 0 | 6 | 54 | 6 |
| `C5_pathway_w20000` | `domain` | `bacteria` | 251 | 1 | 1 | 249 | 2 |
| `C5_pathway_w20000` | `gap_kind` | `already_gap` | 60 | 0 | 1 | 59 | 1 |
| `C5_pathway_w20000` | `gap_kind` | `maskable_gap` | 251 | 1 | 6 | 244 | 7 |
| `C5_pathway_w20000` | `pathway_robust` | `marginal` | 142 | 0 | 0 | 142 | 0 |
| `C5_pathway_w20000` | `pathway_robust` | `robust` | 169 | 1 | 7 | 161 | 8 |
| `C5_organism_w20000` | `domain` | `archaea` | 60 | 0 | 6 | 54 | 6 |
| `C5_organism_w20000` | `domain` | `bacteria` | 251 | 1 | 4 | 246 | 5 |
| `C5_organism_w20000` | `gap_kind` | `already_gap` | 60 | 0 | 2 | 58 | 2 |
| `C5_organism_w20000` | `gap_kind` | `maskable_gap` | 251 | 1 | 8 | 242 | 9 |
| `C5_organism_w20000` | `pathway_robust` | `marginal` | 142 | 0 | 3 | 139 | 3 |
| `C5_organism_w20000` | `pathway_robust` | `robust` | 169 | 1 | 7 | 161 | 8 |
| `C5total_universal_w20000` | `domain` | `archaea` | 60 | 16 | 0 | 44 | 16 |
| `C5total_universal_w20000` | `domain` | `bacteria` | 251 | 43 | 2 | 206 | 45 |
| `C5total_universal_w20000` | `gap_kind` | `already_gap` | 60 | 9 | 0 | 51 | 9 |
| `C5total_universal_w20000` | `gap_kind` | `maskable_gap` | 251 | 50 | 2 | 199 | 52 |
| `C5total_universal_w20000` | `pathway_robust` | `marginal` | 142 | 29 | 1 | 112 | 30 |
| `C5total_universal_w20000` | `pathway_robust` | `robust` | 169 | 30 | 1 | 138 | 31 |
| `C5_pathway_vs_universal_w20000` | `domain` | `archaea` | 60 | 2 | 0 | 58 | 2 |
| `C5_pathway_vs_universal_w20000` | `domain` | `bacteria` | 251 | 2 | 0 | 249 | 2 |
| `C5_pathway_vs_universal_w20000` | `gap_kind` | `already_gap` | 60 | 1 | 0 | 59 | 1 |
| `C5_pathway_vs_universal_w20000` | `gap_kind` | `maskable_gap` | 251 | 3 | 0 | 248 | 3 |
| `C5_pathway_vs_universal_w20000` | `pathway_robust` | `marginal` | 142 | 2 | 0 | 140 | 2 |
| `C5_pathway_vs_universal_w20000` | `pathway_robust` | `robust` | 169 | 2 | 0 | 167 | 2 |
| `C5_organism_vs_universal_w20000` | `domain` | `archaea` | 60 | 2 | 0 | 58 | 2 |
| `C5_organism_vs_universal_w20000` | `domain` | `bacteria` | 251 | 0 | 1 | 250 | 1 |
| `C5_organism_vs_universal_w20000` | `gap_kind` | `already_gap` | 60 | 0 | 0 | 60 | 0 |
| `C5_organism_vs_universal_w20000` | `gap_kind` | `maskable_gap` | 251 | 2 | 1 | 248 | 3 |
| `C5_organism_vs_universal_w20000` | `pathway_robust` | `marginal` | 142 | 0 | 1 | 141 | 1 |
| `C5_organism_vs_universal_w20000` | `pathway_robust` | `robust` | 169 | 2 | 0 | 167 | 2 |
| `C5_organism_vs_pathway_w20000` | `domain` | `archaea` | 60 | 0 | 0 | 60 | 0 |
| `C5_organism_vs_pathway_w20000` | `domain` | `bacteria` | 251 | 0 | 3 | 248 | 3 |
| `C5_organism_vs_pathway_w20000` | `gap_kind` | `already_gap` | 60 | 0 | 1 | 59 | 1 |
| `C5_organism_vs_pathway_w20000` | `gap_kind` | `maskable_gap` | 251 | 0 | 2 | 249 | 2 |
| `C5_organism_vs_pathway_w20000` | `pathway_robust` | `marginal` | 142 | 0 | 3 | 139 | 3 |
| `C5_organism_vs_pathway_w20000` | `pathway_robust` | `robust` | 169 | 0 | 0 | 169 | 0 |
| `C5_universal` | `domain` | `archaea` | 60 | 0 | 8 | 52 | 8 |
| `C5_universal` | `domain` | `bacteria` | 251 | 1 | 5 | 245 | 6 |
| `C5_universal` | `gap_kind` | `already_gap` | 60 | 0 | 2 | 58 | 2 |
| `C5_universal` | `gap_kind` | `maskable_gap` | 251 | 1 | 11 | 239 | 12 |
| `C5_universal` | `pathway_robust` | `marginal` | 142 | 0 | 4 | 138 | 4 |
| `C5_universal` | `pathway_robust` | `robust` | 169 | 1 | 9 | 159 | 10 |
| `C5_pathway` | `domain` | `archaea` | 60 | 0 | 6 | 54 | 6 |
| `C5_pathway` | `domain` | `bacteria` | 251 | 1 | 2 | 248 | 3 |
| `C5_pathway` | `gap_kind` | `already_gap` | 60 | 0 | 1 | 59 | 1 |
| `C5_pathway` | `gap_kind` | `maskable_gap` | 251 | 1 | 7 | 243 | 8 |
| `C5_pathway` | `pathway_robust` | `marginal` | 142 | 0 | 1 | 141 | 1 |
| `C5_pathway` | `pathway_robust` | `robust` | 169 | 1 | 7 | 161 | 8 |
| `C5_organism` | `domain` | `archaea` | 60 | 0 | 6 | 54 | 6 |
| `C5_organism` | `domain` | `bacteria` | 251 | 1 | 5 | 245 | 6 |
| `C5_organism` | `gap_kind` | `already_gap` | 60 | 0 | 2 | 58 | 2 |
| `C5_organism` | `gap_kind` | `maskable_gap` | 251 | 1 | 9 | 241 | 10 |
| `C5_organism` | `pathway_robust` | `marginal` | 142 | 0 | 4 | 138 | 4 |
| `C5_organism` | `pathway_robust` | `robust` | 169 | 1 | 7 | 161 | 8 |
| `C5_pathway_vs_universal` | `domain` | `archaea` | 60 | 2 | 0 | 58 | 2 |
| `C5_pathway_vs_universal` | `domain` | `bacteria` | 251 | 3 | 0 | 248 | 3 |
| `C5_pathway_vs_universal` | `gap_kind` | `already_gap` | 60 | 1 | 0 | 59 | 1 |
| `C5_pathway_vs_universal` | `gap_kind` | `maskable_gap` | 251 | 4 | 0 | 247 | 4 |
| `C5_pathway_vs_universal` | `pathway_robust` | `marginal` | 142 | 3 | 0 | 139 | 3 |
| `C5_pathway_vs_universal` | `pathway_robust` | `robust` | 169 | 2 | 0 | 167 | 2 |
| `C5_organism_vs_universal` | `domain` | `archaea` | 60 | 2 | 0 | 58 | 2 |
| `C5_organism_vs_universal` | `domain` | `bacteria` | 251 | 0 | 0 | 251 | 0 |
| `C5_organism_vs_universal` | `gap_kind` | `already_gap` | 60 | 0 | 0 | 60 | 0 |
| `C5_organism_vs_universal` | `gap_kind` | `maskable_gap` | 251 | 2 | 0 | 249 | 2 |
| `C5_organism_vs_universal` | `pathway_robust` | `marginal` | 142 | 0 | 0 | 142 | 0 |
| `C5_organism_vs_universal` | `pathway_robust` | `robust` | 169 | 2 | 0 | 167 | 2 |
| `C5_organism_vs_pathway` | `domain` | `archaea` | 60 | 0 | 0 | 60 | 0 |
| `C5_organism_vs_pathway` | `domain` | `bacteria` | 251 | 0 | 3 | 248 | 3 |
| `C5_organism_vs_pathway` | `gap_kind` | `already_gap` | 60 | 0 | 1 | 59 | 1 |
| `C5_organism_vs_pathway` | `gap_kind` | `maskable_gap` | 251 | 0 | 2 | 249 | 2 |
| `C5_organism_vs_pathway` | `pathway_robust` | `marginal` | 142 | 0 | 3 | 139 | 3 |
| `C5_organism_vs_pathway` | `pathway_robust` | `robust` | 169 | 0 | 0 | 169 | 0 |
| `S1_stage0_to_sham_gaps` | `domain` | `archaea` | 60 | 0 | 0 | 60 | 0 |
| `S1_stage0_to_sham_gaps` | `domain` | `bacteria` | 251 | 0 | 0 | 251 | 0 |
| `S1_stage0_to_sham_gaps` | `gap_kind` | `already_gap` | 60 | 0 | 0 | 60 | 0 |
| `S1_stage0_to_sham_gaps` | `gap_kind` | `maskable_gap` | 251 | 0 | 0 | 251 | 0 |
| `S1_stage0_to_sham_gaps` | `pathway_robust` | `marginal` | 142 | 0 | 0 | 142 | 0 |
| `S1_stage0_to_sham_gaps` | `pathway_robust` | `robust` | 169 | 0 | 0 | 169 | 0 |
| `S1_sham_gaps_to_sham_gapsfix` | `domain` | `archaea` | 60 | 0 | 0 | 60 | 0 |
| `S1_sham_gaps_to_sham_gapsfix` | `domain` | `bacteria` | 251 | 11 | 0 | 240 | 11 |
| `S1_sham_gaps_to_sham_gapsfix` | `gap_kind` | `already_gap` | 60 | 11 | 0 | 49 | 11 |
| `S1_sham_gaps_to_sham_gapsfix` | `gap_kind` | `maskable_gap` | 251 | 0 | 0 | 251 | 0 |
| `S1_sham_gaps_to_sham_gapsfix` | `pathway_robust` | `marginal` | 142 | 1 | 0 | 141 | 1 |
| `S1_sham_gaps_to_sham_gapsfix` | `pathway_robust` | `robust` | 169 | 10 | 0 | 159 | 10 |
| `S1_sham_gaps_to_sham_ctx_universal_w5000` | `domain` | `archaea` | 60 | 0 | 0 | 60 | 0 |
| `S1_sham_gaps_to_sham_ctx_universal_w5000` | `domain` | `bacteria` | 251 | 0 | 2 | 249 | 2 |
| `S1_sham_gaps_to_sham_ctx_universal_w5000` | `gap_kind` | `already_gap` | 60 | 0 | 2 | 58 | 2 |
| `S1_sham_gaps_to_sham_ctx_universal_w5000` | `gap_kind` | `maskable_gap` | 251 | 0 | 0 | 251 | 0 |
| `S1_sham_gaps_to_sham_ctx_universal_w5000` | `pathway_robust` | `marginal` | 142 | 0 | 1 | 141 | 1 |
| `S1_sham_gaps_to_sham_ctx_universal_w5000` | `pathway_robust` | `robust` | 169 | 0 | 1 | 168 | 1 |
| `S1_stage0_to_sham_ctx_universal_w5000` | `domain` | `archaea` | 60 | 0 | 0 | 60 | 0 |
| `S1_stage0_to_sham_ctx_universal_w5000` | `domain` | `bacteria` | 251 | 9 | 0 | 242 | 9 |
| `S1_stage0_to_sham_ctx_universal_w5000` | `gap_kind` | `already_gap` | 60 | 9 | 0 | 51 | 9 |
| `S1_stage0_to_sham_ctx_universal_w5000` | `gap_kind` | `maskable_gap` | 251 | 0 | 0 | 251 | 0 |
| `S1_stage0_to_sham_ctx_universal_w5000` | `pathway_robust` | `marginal` | 142 | 0 | 0 | 142 | 0 |
| `S1_stage0_to_sham_ctx_universal_w5000` | `pathway_robust` | `robust` | 169 | 9 | 0 | 160 | 9 |
| `S1_sham_gaps_to_sham_ctx_universal_w10000` | `domain` | `archaea` | 60 | 0 | 0 | 60 | 0 |
| `S1_sham_gaps_to_sham_ctx_universal_w10000` | `domain` | `bacteria` | 251 | 0 | 2 | 249 | 2 |
| `S1_sham_gaps_to_sham_ctx_universal_w10000` | `gap_kind` | `already_gap` | 60 | 0 | 2 | 58 | 2 |
| `S1_sham_gaps_to_sham_ctx_universal_w10000` | `gap_kind` | `maskable_gap` | 251 | 0 | 0 | 251 | 0 |
| `S1_sham_gaps_to_sham_ctx_universal_w10000` | `pathway_robust` | `marginal` | 142 | 0 | 1 | 141 | 1 |
| `S1_sham_gaps_to_sham_ctx_universal_w10000` | `pathway_robust` | `robust` | 169 | 0 | 1 | 168 | 1 |
| `S1_stage0_to_sham_ctx_universal_w10000` | `domain` | `archaea` | 60 | 0 | 0 | 60 | 0 |
| `S1_stage0_to_sham_ctx_universal_w10000` | `domain` | `bacteria` | 251 | 9 | 0 | 242 | 9 |
| `S1_stage0_to_sham_ctx_universal_w10000` | `gap_kind` | `already_gap` | 60 | 9 | 0 | 51 | 9 |
| `S1_stage0_to_sham_ctx_universal_w10000` | `gap_kind` | `maskable_gap` | 251 | 0 | 0 | 251 | 0 |
| `S1_stage0_to_sham_ctx_universal_w10000` | `pathway_robust` | `marginal` | 142 | 0 | 0 | 142 | 0 |
| `S1_stage0_to_sham_ctx_universal_w10000` | `pathway_robust` | `robust` | 169 | 9 | 0 | 160 | 9 |
| `S1_sham_gaps_to_sham_ctx_universal_w20000` | `domain` | `archaea` | 60 | 0 | 0 | 60 | 0 |
| `S1_sham_gaps_to_sham_ctx_universal_w20000` | `domain` | `bacteria` | 251 | 0 | 2 | 249 | 2 |
| `S1_sham_gaps_to_sham_ctx_universal_w20000` | `gap_kind` | `already_gap` | 60 | 0 | 2 | 58 | 2 |
| `S1_sham_gaps_to_sham_ctx_universal_w20000` | `gap_kind` | `maskable_gap` | 251 | 0 | 0 | 251 | 0 |
| `S1_sham_gaps_to_sham_ctx_universal_w20000` | `pathway_robust` | `marginal` | 142 | 0 | 1 | 141 | 1 |
| `S1_sham_gaps_to_sham_ctx_universal_w20000` | `pathway_robust` | `robust` | 169 | 0 | 1 | 168 | 1 |
| `S1_stage0_to_sham_ctx_universal_w20000` | `domain` | `archaea` | 60 | 0 | 0 | 60 | 0 |
| `S1_stage0_to_sham_ctx_universal_w20000` | `domain` | `bacteria` | 251 | 9 | 0 | 242 | 9 |
| `S1_stage0_to_sham_ctx_universal_w20000` | `gap_kind` | `already_gap` | 60 | 9 | 0 | 51 | 9 |
| `S1_stage0_to_sham_ctx_universal_w20000` | `gap_kind` | `maskable_gap` | 251 | 0 | 0 | 251 | 0 |
| `S1_stage0_to_sham_ctx_universal_w20000` | `pathway_robust` | `marginal` | 142 | 0 | 0 | 142 | 0 |
| `S1_stage0_to_sham_ctx_universal_w20000` | `pathway_robust` | `robust` | 169 | 9 | 0 | 160 | 9 |
| `S1_sham_gaps_to_sham_ctx_universal` | `domain` | `archaea` | 60 | 0 | 0 | 60 | 0 |
| `S1_sham_gaps_to_sham_ctx_universal` | `domain` | `bacteria` | 251 | 0 | 2 | 249 | 2 |
| `S1_sham_gaps_to_sham_ctx_universal` | `gap_kind` | `already_gap` | 60 | 0 | 2 | 58 | 2 |
| `S1_sham_gaps_to_sham_ctx_universal` | `gap_kind` | `maskable_gap` | 251 | 0 | 0 | 251 | 0 |
| `S1_sham_gaps_to_sham_ctx_universal` | `pathway_robust` | `marginal` | 142 | 0 | 1 | 141 | 1 |
| `S1_sham_gaps_to_sham_ctx_universal` | `pathway_robust` | `robust` | 169 | 0 | 1 | 168 | 1 |
| `U1_unmasked_gaps_vs_gaps` | `domain` | `archaea` | 60 | 0 | 0 | 60 | 0 |
| `U1_unmasked_gaps_vs_gaps` | `domain` | `bacteria` | 251 | 0 | 0 | 251 | 0 |
| `U1_unmasked_gaps_vs_gaps` | `gap_kind` | `already_gap` | 60 | 0 | 0 | 60 | 0 |
| `U1_unmasked_gaps_vs_gaps` | `gap_kind` | `maskable_gap` | 251 | 0 | 0 | 251 | 0 |
| `U1_unmasked_gaps_vs_gaps` | `pathway_robust` | `marginal` | 142 | 0 | 0 | 142 | 0 |
| `U1_unmasked_gaps_vs_gaps` | `pathway_robust` | `robust` | 169 | 0 | 0 | 169 | 0 |

## 4. Learned-Stage Candidate Counts

| Arm | Candidates | Rescored | Unsupported | Unchanged unsupported | Retracted | Scored only |
|---|---:|---:|---:|---:|---:|---:|
| - | 0 | 0 | 0 | 0 | 0 | 0 |

## 5. Masking & Build Audit Summary

| Check | Count |
|---|---:|
| `total_cells` | 311 |
| `complete_cells` | 311 |
| `candidate_universe_preserved_cells` | 311 |
| `masked_gap_present_cells` | 311 |
| `gaps_ctxjar_identical_cells` | 311 |
| `already_gap_unmasked_total` | 60 |
| `already_gap_unmasked_matches` | 60 |
