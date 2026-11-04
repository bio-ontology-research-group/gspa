# Row-5 Masked-Context Evaluation (`scope=heldout`)

Chosen context window (bp): `10000`

## 1. Contrasts (Overall and Per Claim Set)

| Contrast | Claim set | Cells | Improved | Harmed | Unchanged | Discordant | Critical k | p (improve) | Harm elsewhere | Verdict |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `C1_stage0_to_gaps` | **all** | 253 | 0 | 0 | 253 | 0 | - | 1 | 0 | no detectable effect |
| `C1_stage0_to_gaps` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C1_stage0_to_gaps` | `diamond` | 61 | 0 | 0 | 61 | 0 | - | 1 | 0 | no detectable effect |
| `C1_stage0_to_gaps` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C1_stage0_to_gaps` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C1_stage0_to_gaps` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C1_stage0_to_gaps` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C1_stage0_to_gaps` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C1b_gaps_to_gapsfix` | **all** | 253 | 47 | 2 | 204 | 49 | 31 | 2.178e-12 | 0 | systematic improvement |
| `C1b_gaps_to_gapsfix` | `clean` | 32 | 6 | 0 | 26 | 6 | 6 | 0.01562 | 0 | systematic improvement |
| `C1b_gaps_to_gapsfix` | `diamond` | 61 | 22 | 0 | 39 | 22 | 16 | 2.384e-07 | 0 | systematic improvement |
| `C1b_gaps_to_gapsfix` | `dpz` | 32 | 4 | 0 | 28 | 4 | - | 0.0625 | 0 | no detectable effect |
| `C1b_gaps_to_gapsfix` | `enzbert` | 32 | 2 | 0 | 30 | 2 | - | 0.25 | 0 | no detectable effect |
| `C1b_gaps_to_gapsfix` | `graphec` | 32 | 8 | 0 | 24 | 8 | 7 | 0.003906 | 0 | systematic improvement |
| `C1b_gaps_to_gapsfix` | `mapred` | 32 | 5 | 2 | 25 | 7 | 7 | 0.2266 | 0 | no detectable effect |
| `C1b_gaps_to_gapsfix` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C1c_stage0_to_gapsfix` | **all** | 253 | 47 | 2 | 204 | 49 | 31 | 2.178e-12 | 0 | systematic improvement |
| `C1c_stage0_to_gapsfix` | `clean` | 32 | 6 | 0 | 26 | 6 | 6 | 0.01562 | 0 | systematic improvement |
| `C1c_stage0_to_gapsfix` | `diamond` | 61 | 22 | 0 | 39 | 22 | 16 | 2.384e-07 | 0 | systematic improvement |
| `C1c_stage0_to_gapsfix` | `dpz` | 32 | 4 | 0 | 28 | 4 | - | 0.0625 | 0 | no detectable effect |
| `C1c_stage0_to_gapsfix` | `enzbert` | 32 | 2 | 0 | 30 | 2 | - | 0.25 | 0 | no detectable effect |
| `C1c_stage0_to_gapsfix` | `graphec` | 32 | 8 | 0 | 24 | 8 | 7 | 0.003906 | 0 | systematic improvement |
| `C1c_stage0_to_gapsfix` | `mapred` | 32 | 5 | 2 | 25 | 7 | 7 | 0.2266 | 0 | no detectable effect |
| `C1c_stage0_to_gapsfix` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal_w5000` | **all** | 253 | 1 | 6 | 246 | 7 | 7 | 0.9922 | 14 | no detectable effect |
| `C5_universal_w5000` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal_w5000` | `diamond` | 61 | 1 | 6 | 54 | 7 | 7 | 0.9922 | 14 | no detectable effect |
| `C5_universal_w5000` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal_w5000` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal_w5000` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal_w5000` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal_w5000` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_w5000` | **all** | 253 | 1 | 2 | 250 | 3 | - | 0.875 | 2 | no detectable effect |
| `C5_pathway_w5000` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_w5000` | `diamond` | 61 | 1 | 2 | 58 | 3 | - | 0.875 | 2 | no detectable effect |
| `C5_pathway_w5000` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_w5000` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_w5000` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_w5000` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_w5000` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_w5000` | **all** | 253 | 1 | 5 | 247 | 6 | 6 | 0.9844 | 13 | no detectable effect |
| `C5_organism_w5000` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_w5000` | `diamond` | 61 | 1 | 5 | 55 | 6 | 6 | 0.9844 | 13 | no detectable effect |
| `C5_organism_w5000` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_w5000` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_w5000` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_w5000` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_w5000` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5total_universal_w5000` | **all** | 253 | 41 | 2 | 210 | 43 | 28 | 1.077e-10 | 0 | systematic improvement |
| `C5total_universal_w5000` | `clean` | 32 | 6 | 0 | 26 | 6 | 6 | 0.01562 | 0 | systematic improvement |
| `C5total_universal_w5000` | `diamond` | 61 | 16 | 0 | 45 | 16 | 12 | 1.526e-05 | 0 | systematic improvement |
| `C5total_universal_w5000` | `dpz` | 32 | 4 | 0 | 28 | 4 | - | 0.0625 | 0 | no detectable effect |
| `C5total_universal_w5000` | `enzbert` | 32 | 2 | 0 | 30 | 2 | - | 0.25 | 0 | no detectable effect |
| `C5total_universal_w5000` | `graphec` | 32 | 8 | 0 | 24 | 8 | 7 | 0.003906 | 0 | systematic improvement |
| `C5total_universal_w5000` | `mapred` | 32 | 5 | 2 | 25 | 7 | 7 | 0.2266 | 0 | no detectable effect |
| `C5total_universal_w5000` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w5000` | **all** | 253 | 4 | 0 | 249 | 4 | - | 0.0625 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w5000` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w5000` | `diamond` | 61 | 4 | 0 | 57 | 4 | - | 0.0625 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w5000` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w5000` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w5000` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w5000` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w5000` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w5000` | **all** | 253 | 1 | 0 | 252 | 1 | - | 0.5 | 0 | no detectable effect |
| `C5_organism_vs_universal_w5000` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w5000` | `diamond` | 61 | 1 | 0 | 60 | 1 | - | 0.5 | 0 | no detectable effect |
| `C5_organism_vs_universal_w5000` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w5000` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w5000` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w5000` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w5000` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway_w5000` | **all** | 253 | 0 | 3 | 250 | 3 | - | 1 | 11 | no detectable effect |
| `C5_organism_vs_pathway_w5000` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway_w5000` | `diamond` | 61 | 0 | 3 | 58 | 3 | - | 1 | 11 | no detectable effect |
| `C5_organism_vs_pathway_w5000` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway_w5000` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway_w5000` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway_w5000` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway_w5000` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal_w10000` | **all** | 253 | 1 | 6 | 246 | 7 | 7 | 0.9922 | 14 | no detectable effect |
| `C5_universal_w10000` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal_w10000` | `diamond` | 61 | 1 | 6 | 54 | 7 | 7 | 0.9922 | 14 | no detectable effect |
| `C5_universal_w10000` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal_w10000` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal_w10000` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal_w10000` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal_w10000` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_w10000` | **all** | 253 | 1 | 2 | 250 | 3 | - | 0.875 | 2 | no detectable effect |
| `C5_pathway_w10000` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_w10000` | `diamond` | 61 | 1 | 2 | 58 | 3 | - | 0.875 | 2 | no detectable effect |
| `C5_pathway_w10000` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_w10000` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_w10000` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_w10000` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_w10000` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_w10000` | **all** | 253 | 1 | 5 | 247 | 6 | 6 | 0.9844 | 13 | no detectable effect |
| `C5_organism_w10000` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_w10000` | `diamond` | 61 | 1 | 5 | 55 | 6 | 6 | 0.9844 | 13 | no detectable effect |
| `C5_organism_w10000` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_w10000` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_w10000` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_w10000` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_w10000` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5total_universal_w10000` | **all** | 253 | 41 | 2 | 210 | 43 | 28 | 1.077e-10 | 0 | systematic improvement |
| `C5total_universal_w10000` | `clean` | 32 | 6 | 0 | 26 | 6 | 6 | 0.01562 | 0 | systematic improvement |
| `C5total_universal_w10000` | `diamond` | 61 | 16 | 0 | 45 | 16 | 12 | 1.526e-05 | 0 | systematic improvement |
| `C5total_universal_w10000` | `dpz` | 32 | 4 | 0 | 28 | 4 | - | 0.0625 | 0 | no detectable effect |
| `C5total_universal_w10000` | `enzbert` | 32 | 2 | 0 | 30 | 2 | - | 0.25 | 0 | no detectable effect |
| `C5total_universal_w10000` | `graphec` | 32 | 8 | 0 | 24 | 8 | 7 | 0.003906 | 0 | systematic improvement |
| `C5total_universal_w10000` | `mapred` | 32 | 5 | 2 | 25 | 7 | 7 | 0.2266 | 0 | no detectable effect |
| `C5total_universal_w10000` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w10000` | **all** | 253 | 4 | 0 | 249 | 4 | - | 0.0625 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w10000` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w10000` | `diamond` | 61 | 4 | 0 | 57 | 4 | - | 0.0625 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w10000` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w10000` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w10000` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w10000` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w10000` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w10000` | **all** | 253 | 1 | 0 | 252 | 1 | - | 0.5 | 0 | no detectable effect |
| `C5_organism_vs_universal_w10000` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w10000` | `diamond` | 61 | 1 | 0 | 60 | 1 | - | 0.5 | 0 | no detectable effect |
| `C5_organism_vs_universal_w10000` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w10000` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w10000` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w10000` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w10000` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway_w10000` | **all** | 253 | 0 | 3 | 250 | 3 | - | 1 | 11 | no detectable effect |
| `C5_organism_vs_pathway_w10000` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway_w10000` | `diamond` | 61 | 0 | 3 | 58 | 3 | - | 1 | 11 | no detectable effect |
| `C5_organism_vs_pathway_w10000` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway_w10000` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway_w10000` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway_w10000` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway_w10000` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal_w20000` | **all** | 253 | 1 | 4 | 248 | 5 | 5 | 0.9688 | 14 | no detectable effect |
| `C5_universal_w20000` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal_w20000` | `diamond` | 61 | 1 | 4 | 56 | 5 | 5 | 0.9688 | 14 | no detectable effect |
| `C5_universal_w20000` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal_w20000` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal_w20000` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal_w20000` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal_w20000` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_w20000` | **all** | 253 | 1 | 1 | 251 | 2 | - | 0.75 | 2 | no detectable effect |
| `C5_pathway_w20000` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_w20000` | `diamond` | 61 | 1 | 1 | 59 | 2 | - | 0.75 | 2 | no detectable effect |
| `C5_pathway_w20000` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_w20000` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_w20000` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_w20000` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_w20000` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_w20000` | **all** | 253 | 1 | 4 | 248 | 5 | 5 | 0.9688 | 13 | no detectable effect |
| `C5_organism_w20000` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_w20000` | `diamond` | 61 | 1 | 4 | 56 | 5 | 5 | 0.9688 | 13 | no detectable effect |
| `C5_organism_w20000` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_w20000` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_w20000` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_w20000` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_w20000` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5total_universal_w20000` | **all** | 253 | 43 | 2 | 208 | 45 | 29 | 2.944e-11 | 0 | systematic improvement |
| `C5total_universal_w20000` | `clean` | 32 | 6 | 0 | 26 | 6 | 6 | 0.01562 | 0 | systematic improvement |
| `C5total_universal_w20000` | `diamond` | 61 | 18 | 0 | 43 | 18 | 13 | 3.815e-06 | 0 | systematic improvement |
| `C5total_universal_w20000` | `dpz` | 32 | 4 | 0 | 28 | 4 | - | 0.0625 | 0 | no detectable effect |
| `C5total_universal_w20000` | `enzbert` | 32 | 2 | 0 | 30 | 2 | - | 0.25 | 0 | no detectable effect |
| `C5total_universal_w20000` | `graphec` | 32 | 8 | 0 | 24 | 8 | 7 | 0.003906 | 0 | systematic improvement |
| `C5total_universal_w20000` | `mapred` | 32 | 5 | 2 | 25 | 7 | 7 | 0.2266 | 0 | no detectable effect |
| `C5total_universal_w20000` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w20000` | **all** | 253 | 3 | 0 | 250 | 3 | - | 0.125 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w20000` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w20000` | `diamond` | 61 | 3 | 0 | 58 | 3 | - | 0.125 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w20000` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w20000` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w20000` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w20000` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w20000` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w20000` | **all** | 253 | 1 | 1 | 251 | 2 | - | 0.75 | 0 | no detectable effect |
| `C5_organism_vs_universal_w20000` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w20000` | `diamond` | 61 | 1 | 1 | 59 | 2 | - | 0.75 | 0 | no detectable effect |
| `C5_organism_vs_universal_w20000` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w20000` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w20000` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w20000` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w20000` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway_w20000` | **all** | 253 | 0 | 3 | 250 | 3 | - | 1 | 11 | no detectable effect |
| `C5_organism_vs_pathway_w20000` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway_w20000` | `diamond` | 61 | 0 | 3 | 58 | 3 | - | 1 | 11 | no detectable effect |
| `C5_organism_vs_pathway_w20000` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway_w20000` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway_w20000` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway_w20000` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway_w20000` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal` | **all** | 253 | 1 | 6 | 246 | 7 | 7 | 0.9922 | 14 | no detectable effect |
| `C5_universal` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal` | `diamond` | 61 | 1 | 6 | 54 | 7 | 7 | 0.9922 | 14 | no detectable effect |
| `C5_universal` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway` | **all** | 253 | 1 | 2 | 250 | 3 | - | 0.875 | 2 | no detectable effect |
| `C5_pathway` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway` | `diamond` | 61 | 1 | 2 | 58 | 3 | - | 0.875 | 2 | no detectable effect |
| `C5_pathway` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism` | **all** | 253 | 1 | 5 | 247 | 6 | 6 | 0.9844 | 13 | no detectable effect |
| `C5_organism` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism` | `diamond` | 61 | 1 | 5 | 55 | 6 | 6 | 0.9844 | 13 | no detectable effect |
| `C5_organism` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal` | **all** | 253 | 4 | 0 | 249 | 4 | - | 0.0625 | 0 | no detectable effect |
| `C5_pathway_vs_universal` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal` | `diamond` | 61 | 4 | 0 | 57 | 4 | - | 0.0625 | 0 | no detectable effect |
| `C5_pathway_vs_universal` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal` | **all** | 253 | 1 | 0 | 252 | 1 | - | 0.5 | 0 | no detectable effect |
| `C5_organism_vs_universal` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal` | `diamond` | 61 | 1 | 0 | 60 | 1 | - | 0.5 | 0 | no detectable effect |
| `C5_organism_vs_universal` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway` | **all** | 253 | 0 | 3 | 250 | 3 | - | 1 | 11 | no detectable effect |
| `C5_organism_vs_pathway` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway` | `diamond` | 61 | 0 | 3 | 58 | 3 | - | 1 | 11 | no detectable effect |
| `C5_organism_vs_pathway` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_gaps` | **all** | 253 | 0 | 0 | 253 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_gaps` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_gaps` | `diamond` | 61 | 0 | 0 | 61 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_gaps` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_gaps` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_gaps` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_gaps` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_gaps` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_gapsfix` | **all** | 253 | 11 | 0 | 242 | 11 | 9 | 0.0004883 | 0 | systematic improvement |
| `S1_sham_gaps_to_sham_gapsfix` | `clean` | 32 | 3 | 0 | 29 | 3 | - | 0.125 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_gapsfix` | `diamond` | 61 | 3 | 0 | 58 | 3 | - | 0.125 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_gapsfix` | `dpz` | 32 | 3 | 0 | 29 | 3 | - | 0.125 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_gapsfix` | `enzbert` | 32 | 1 | 0 | 31 | 1 | - | 0.5 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_gapsfix` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_gapsfix` | `mapred` | 32 | 1 | 0 | 31 | 1 | - | 0.5 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_gapsfix` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w5000` | **all** | 253 | 0 | 2 | 251 | 2 | - | 1 | 13 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w5000` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w5000` | `diamond` | 61 | 0 | 2 | 59 | 2 | - | 1 | 13 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w5000` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w5000` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w5000` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w5000` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w5000` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w5000` | **all** | 253 | 9 | 0 | 244 | 9 | 8 | 0.001953 | 0 | systematic improvement |
| `S1_stage0_to_sham_ctx_universal_w5000` | `clean` | 32 | 3 | 0 | 29 | 3 | - | 0.125 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w5000` | `diamond` | 61 | 1 | 0 | 60 | 1 | - | 0.5 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w5000` | `dpz` | 32 | 3 | 0 | 29 | 3 | - | 0.125 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w5000` | `enzbert` | 32 | 1 | 0 | 31 | 1 | - | 0.5 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w5000` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w5000` | `mapred` | 32 | 1 | 0 | 31 | 1 | - | 0.5 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w5000` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w10000` | **all** | 253 | 0 | 2 | 251 | 2 | - | 1 | 13 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w10000` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w10000` | `diamond` | 61 | 0 | 2 | 59 | 2 | - | 1 | 13 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w10000` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w10000` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w10000` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w10000` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w10000` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w10000` | **all** | 253 | 9 | 0 | 244 | 9 | 8 | 0.001953 | 0 | systematic improvement |
| `S1_stage0_to_sham_ctx_universal_w10000` | `clean` | 32 | 3 | 0 | 29 | 3 | - | 0.125 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w10000` | `diamond` | 61 | 1 | 0 | 60 | 1 | - | 0.5 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w10000` | `dpz` | 32 | 3 | 0 | 29 | 3 | - | 0.125 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w10000` | `enzbert` | 32 | 1 | 0 | 31 | 1 | - | 0.5 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w10000` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w10000` | `mapred` | 32 | 1 | 0 | 31 | 1 | - | 0.5 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w10000` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w20000` | **all** | 253 | 0 | 2 | 251 | 2 | - | 1 | 13 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w20000` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w20000` | `diamond` | 61 | 0 | 2 | 59 | 2 | - | 1 | 13 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w20000` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w20000` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w20000` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w20000` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w20000` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w20000` | **all** | 253 | 9 | 0 | 244 | 9 | 8 | 0.001953 | 0 | systematic improvement |
| `S1_stage0_to_sham_ctx_universal_w20000` | `clean` | 32 | 3 | 0 | 29 | 3 | - | 0.125 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w20000` | `diamond` | 61 | 1 | 0 | 60 | 1 | - | 0.5 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w20000` | `dpz` | 32 | 3 | 0 | 29 | 3 | - | 0.125 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w20000` | `enzbert` | 32 | 1 | 0 | 31 | 1 | - | 0.5 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w20000` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w20000` | `mapred` | 32 | 1 | 0 | 31 | 1 | - | 0.5 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w20000` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal` | **all** | 253 | 0 | 2 | 251 | 2 | - | 1 | 13 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal` | `diamond` | 61 | 0 | 2 | 59 | 2 | - | 1 | 13 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `U1_unmasked_gaps_vs_gaps` | **all** | 253 | 0 | 0 | 253 | 0 | - | 1 | 0 | no detectable effect |
| `U1_unmasked_gaps_vs_gaps` | `clean` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `U1_unmasked_gaps_vs_gaps` | `diamond` | 61 | 0 | 0 | 61 | 0 | - | 1 | 0 | no detectable effect |
| `U1_unmasked_gaps_vs_gaps` | `dpz` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `U1_unmasked_gaps_vs_gaps` | `enzbert` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `U1_unmasked_gaps_vs_gaps` | `graphec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `U1_unmasked_gaps_vs_gaps` | `mapred` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |
| `U1_unmasked_gaps_vs_gaps` | `topec` | 32 | 0 | 0 | 32 | 0 | - | 1 | 0 | no detectable effect |

## 2. Harm Elsewhere and Cost

| Contrast | Improved | Harmed | Harm elsewhere | Unconfirmed calls created | Delta calls >= 0.5 | Changed scores | Bootstrap 95% CI | Wilcoxon p (1-sided) |
|---|---:|---:|---:|---:|---:|---:|---|---:|
| `C1_stage0_to_gaps` | 0 | 0 | 0 | 0 | 7412 | 111306 | [0.0, 0.0] | 1 |
| `C1b_gaps_to_gapsfix` | 47 | 2 | 0 | 20 | 2889 | 10792 | [30.0, 61.0] | 0.02116 |
| `C1c_stage0_to_gapsfix` | 47 | 2 | 0 | 20 | 10301 | 114261 | [30.0, 61.0] | 0.02116 |
| `C5_universal_w5000` | 1 | 6 | 14 | 0 | -674 | 1551 | [-10.0, 0.0] | 0.1587 |
| `C5_pathway_w5000` | 1 | 2 | 2 | 0 | -653 | 1483 | [-4.0, 2.0] | 0.1587 |
| `C5_organism_w5000` | 1 | 5 | 13 | 0 | -502 | 1342 | [-8.0, 0.0] | 0.1587 |
| `C5total_universal_w5000` | 41 | 2 | 0 | 17 | 2215 | 11741 | [25.0, 54.0] | 0.01648 |
| `C5_pathway_vs_universal_w5000` | 4 | 0 | 0 | 0 | 21 | 489 | [1.0, 8.0] | 1 |
| `C5_organism_vs_universal_w5000` | 1 | 0 | 0 | 0 | 172 | 423 | [0.0, 3.0] | 1 |
| `C5_organism_vs_pathway_w5000` | 0 | 3 | 11 | 0 | 151 | 542 | [-6.0, 0.0] | 1 |
| `C5_universal_w10000` | 1 | 6 | 14 | 0 | -659 | 1549 | [-10.0, 0.0] | 0.1587 |
| `C5_pathway_w10000` | 1 | 2 | 2 | 0 | -641 | 1483 | [-4.0, 2.0] | 0.1587 |
| `C5_organism_w10000` | 1 | 5 | 13 | 0 | -484 | 1337 | [-8.0, 0.0] | 0.1587 |
| `C5total_universal_w10000` | 41 | 2 | 0 | 17 | 2230 | 11737 | [25.0, 54.0] | 0.01648 |
| `C5_pathway_vs_universal_w10000` | 4 | 0 | 0 | 0 | 18 | 502 | [1.0, 8.0] | 1 |
| `C5_organism_vs_universal_w10000` | 1 | 0 | 0 | 0 | 175 | 449 | [0.0, 3.0] | 1 |
| `C5_organism_vs_pathway_w10000` | 0 | 3 | 11 | 0 | 157 | 544 | [-6.0, 0.0] | 1 |
| `C5_universal_w20000` | 1 | 4 | 14 | 0 | -636 | 1561 | [-7.0, 1.0] | 0.1587 |
| `C5_pathway_w20000` | 1 | 1 | 2 | 0 | -637 | 1485 | [-3.0, 3.0] | 0.1587 |
| `C5_organism_w20000` | 1 | 4 | 13 | 0 | -477 | 1346 | [-7.0, 1.0] | 0.1587 |
| `C5total_universal_w20000` | 43 | 2 | 0 | 17 | 2253 | 11736 | [27.0, 56.0] | 0.01648 |
| `C5_pathway_vs_universal_w20000` | 3 | 0 | 0 | 0 | -1 | 585 | [0.0, 6.0] | 1 |
| `C5_organism_vs_universal_w20000` | 1 | 1 | 0 | 0 | 159 | 507 | [-3.0, 3.0] | 1 |
| `C5_organism_vs_pathway_w20000` | 0 | 3 | 11 | 0 | 160 | 563 | [-6.0, 0.0] | 1 |
| `C5_universal` | 1 | 6 | 14 | 0 | -659 | 1549 | [-10.0, 0.0] | 0.1587 |
| `C5_pathway` | 1 | 2 | 2 | 0 | -641 | 1483 | [-4.0, 2.0] | 0.1587 |
| `C5_organism` | 1 | 5 | 13 | 0 | -484 | 1337 | [-8.0, 0.0] | 0.1587 |
| `C5_pathway_vs_universal` | 4 | 0 | 0 | 0 | 18 | 502 | [1.0, 8.0] | 1 |
| `C5_organism_vs_universal` | 1 | 0 | 0 | 0 | 175 | 449 | [0.0, 3.0] | 1 |
| `C5_organism_vs_pathway` | 0 | 3 | 11 | 0 | 157 | 544 | [-6.0, 0.0] | 1 |
| `S1_stage0_to_sham_gaps` | 0 | 0 | 0 | 0 | 7413 | 111307 | [0.0, 0.0] | 1 |
| `S1_sham_gaps_to_sham_gapsfix` | 11 | 0 | 0 | 20 | 2845 | 10749 | [3.0, 21.0] | 0.05282 |
| `S1_sham_gaps_to_sham_ctx_universal_w5000` | 0 | 2 | 13 | 0 | -668 | 1521 | [-5.0, 0.0] | 1 |
| `S1_stage0_to_sham_ctx_universal_w5000` | 9 | 0 | 0 | 17 | 9590 | 114233 | [2.0, 18.0] | 0.05282 |
| `S1_sham_gaps_to_sham_ctx_universal_w10000` | 0 | 2 | 13 | 0 | -654 | 1517 | [-5.0, 0.0] | 1 |
| `S1_stage0_to_sham_ctx_universal_w10000` | 9 | 0 | 0 | 17 | 9604 | 114233 | [2.0, 18.0] | 0.05282 |
| `S1_sham_gaps_to_sham_ctx_universal_w20000` | 0 | 2 | 13 | 0 | -633 | 1530 | [-5.0, 0.0] | 1 |
| `S1_stage0_to_sham_ctx_universal_w20000` | 9 | 0 | 0 | 17 | 9625 | 114233 | [2.0, 18.0] | 0.05282 |
| `S1_sham_gaps_to_sham_ctx_universal` | 0 | 2 | 13 | 0 | -654 | 1517 | [-5.0, 0.0] | 1 |
| `U1_unmasked_gaps_vs_gaps` | 0 | 0 | 0 | 0 | 19 | 57 | [0.0, 0.0] | 1 |

## 3. Strata Breakdown

| Contrast | Stratum axis | Stratum | Cells | Improved | Harmed | Unchanged | Discordant |
|---|---|---|---:|---:|---:|---:|---:|
| `C1_stage0_to_gaps` | `domain` | `archaea` | 6 | 0 | 0 | 6 | 0 |
| `C1_stage0_to_gaps` | `domain` | `bacteria` | 247 | 0 | 0 | 247 | 0 |
| `C1_stage0_to_gaps` | `gap_kind` | `already_gap` | 58 | 0 | 0 | 58 | 0 |
| `C1_stage0_to_gaps` | `gap_kind` | `maskable_gap` | 195 | 0 | 0 | 195 | 0 |
| `C1_stage0_to_gaps` | `pathway_robust` | `marginal` | 122 | 0 | 0 | 122 | 0 |
| `C1_stage0_to_gaps` | `pathway_robust` | `robust` | 131 | 0 | 0 | 131 | 0 |
| `C1b_gaps_to_gapsfix` | `domain` | `archaea` | 6 | 3 | 0 | 3 | 3 |
| `C1b_gaps_to_gapsfix` | `domain` | `bacteria` | 247 | 44 | 2 | 201 | 46 |
| `C1b_gaps_to_gapsfix` | `gap_kind` | `already_gap` | 58 | 11 | 0 | 47 | 11 |
| `C1b_gaps_to_gapsfix` | `gap_kind` | `maskable_gap` | 195 | 36 | 2 | 157 | 38 |
| `C1b_gaps_to_gapsfix` | `pathway_robust` | `marginal` | 122 | 26 | 1 | 95 | 27 |
| `C1b_gaps_to_gapsfix` | `pathway_robust` | `robust` | 131 | 21 | 1 | 109 | 22 |
| `C1c_stage0_to_gapsfix` | `domain` | `archaea` | 6 | 3 | 0 | 3 | 3 |
| `C1c_stage0_to_gapsfix` | `domain` | `bacteria` | 247 | 44 | 2 | 201 | 46 |
| `C1c_stage0_to_gapsfix` | `gap_kind` | `already_gap` | 58 | 11 | 0 | 47 | 11 |
| `C1c_stage0_to_gapsfix` | `gap_kind` | `maskable_gap` | 195 | 36 | 2 | 157 | 38 |
| `C1c_stage0_to_gapsfix` | `pathway_robust` | `marginal` | 122 | 26 | 1 | 95 | 27 |
| `C1c_stage0_to_gapsfix` | `pathway_robust` | `robust` | 131 | 21 | 1 | 109 | 22 |
| `C5_universal_w5000` | `domain` | `archaea` | 6 | 0 | 1 | 5 | 1 |
| `C5_universal_w5000` | `domain` | `bacteria` | 247 | 1 | 5 | 241 | 6 |
| `C5_universal_w5000` | `gap_kind` | `already_gap` | 58 | 0 | 2 | 56 | 2 |
| `C5_universal_w5000` | `gap_kind` | `maskable_gap` | 195 | 1 | 4 | 190 | 5 |
| `C5_universal_w5000` | `pathway_robust` | `marginal` | 122 | 0 | 4 | 118 | 4 |
| `C5_universal_w5000` | `pathway_robust` | `robust` | 131 | 1 | 2 | 128 | 3 |
| `C5_pathway_w5000` | `domain` | `archaea` | 6 | 0 | 0 | 6 | 0 |
| `C5_pathway_w5000` | `domain` | `bacteria` | 247 | 1 | 2 | 244 | 3 |
| `C5_pathway_w5000` | `gap_kind` | `already_gap` | 58 | 0 | 1 | 57 | 1 |
| `C5_pathway_w5000` | `gap_kind` | `maskable_gap` | 195 | 1 | 1 | 193 | 2 |
| `C5_pathway_w5000` | `pathway_robust` | `marginal` | 122 | 0 | 1 | 121 | 1 |
| `C5_pathway_w5000` | `pathway_robust` | `robust` | 131 | 1 | 1 | 129 | 2 |
| `C5_organism_w5000` | `domain` | `archaea` | 6 | 0 | 0 | 6 | 0 |
| `C5_organism_w5000` | `domain` | `bacteria` | 247 | 1 | 5 | 241 | 6 |
| `C5_organism_w5000` | `gap_kind` | `already_gap` | 58 | 0 | 2 | 56 | 2 |
| `C5_organism_w5000` | `gap_kind` | `maskable_gap` | 195 | 1 | 3 | 191 | 4 |
| `C5_organism_w5000` | `pathway_robust` | `marginal` | 122 | 0 | 4 | 118 | 4 |
| `C5_organism_w5000` | `pathway_robust` | `robust` | 131 | 1 | 1 | 129 | 2 |
| `C5total_universal_w5000` | `domain` | `archaea` | 6 | 2 | 0 | 4 | 2 |
| `C5total_universal_w5000` | `domain` | `bacteria` | 247 | 39 | 2 | 206 | 41 |
| `C5total_universal_w5000` | `gap_kind` | `already_gap` | 58 | 9 | 0 | 49 | 9 |
| `C5total_universal_w5000` | `gap_kind` | `maskable_gap` | 195 | 32 | 2 | 161 | 34 |
| `C5total_universal_w5000` | `pathway_robust` | `marginal` | 122 | 22 | 1 | 99 | 23 |
| `C5total_universal_w5000` | `pathway_robust` | `robust` | 131 | 19 | 1 | 111 | 20 |
| `C5_pathway_vs_universal_w5000` | `domain` | `archaea` | 6 | 1 | 0 | 5 | 1 |
| `C5_pathway_vs_universal_w5000` | `domain` | `bacteria` | 247 | 3 | 0 | 244 | 3 |
| `C5_pathway_vs_universal_w5000` | `gap_kind` | `already_gap` | 58 | 1 | 0 | 57 | 1 |
| `C5_pathway_vs_universal_w5000` | `gap_kind` | `maskable_gap` | 195 | 3 | 0 | 192 | 3 |
| `C5_pathway_vs_universal_w5000` | `pathway_robust` | `marginal` | 122 | 3 | 0 | 119 | 3 |
| `C5_pathway_vs_universal_w5000` | `pathway_robust` | `robust` | 131 | 1 | 0 | 130 | 1 |
| `C5_organism_vs_universal_w5000` | `domain` | `archaea` | 6 | 1 | 0 | 5 | 1 |
| `C5_organism_vs_universal_w5000` | `domain` | `bacteria` | 247 | 0 | 0 | 247 | 0 |
| `C5_organism_vs_universal_w5000` | `gap_kind` | `already_gap` | 58 | 0 | 0 | 58 | 0 |
| `C5_organism_vs_universal_w5000` | `gap_kind` | `maskable_gap` | 195 | 1 | 0 | 194 | 1 |
| `C5_organism_vs_universal_w5000` | `pathway_robust` | `marginal` | 122 | 0 | 0 | 122 | 0 |
| `C5_organism_vs_universal_w5000` | `pathway_robust` | `robust` | 131 | 1 | 0 | 130 | 1 |
| `C5_organism_vs_pathway_w5000` | `domain` | `archaea` | 6 | 0 | 0 | 6 | 0 |
| `C5_organism_vs_pathway_w5000` | `domain` | `bacteria` | 247 | 0 | 3 | 244 | 3 |
| `C5_organism_vs_pathway_w5000` | `gap_kind` | `already_gap` | 58 | 0 | 1 | 57 | 1 |
| `C5_organism_vs_pathway_w5000` | `gap_kind` | `maskable_gap` | 195 | 0 | 2 | 193 | 2 |
| `C5_organism_vs_pathway_w5000` | `pathway_robust` | `marginal` | 122 | 0 | 3 | 119 | 3 |
| `C5_organism_vs_pathway_w5000` | `pathway_robust` | `robust` | 131 | 0 | 0 | 131 | 0 |
| `C5_universal_w10000` | `domain` | `archaea` | 6 | 0 | 1 | 5 | 1 |
| `C5_universal_w10000` | `domain` | `bacteria` | 247 | 1 | 5 | 241 | 6 |
| `C5_universal_w10000` | `gap_kind` | `already_gap` | 58 | 0 | 2 | 56 | 2 |
| `C5_universal_w10000` | `gap_kind` | `maskable_gap` | 195 | 1 | 4 | 190 | 5 |
| `C5_universal_w10000` | `pathway_robust` | `marginal` | 122 | 0 | 4 | 118 | 4 |
| `C5_universal_w10000` | `pathway_robust` | `robust` | 131 | 1 | 2 | 128 | 3 |
| `C5_pathway_w10000` | `domain` | `archaea` | 6 | 0 | 0 | 6 | 0 |
| `C5_pathway_w10000` | `domain` | `bacteria` | 247 | 1 | 2 | 244 | 3 |
| `C5_pathway_w10000` | `gap_kind` | `already_gap` | 58 | 0 | 1 | 57 | 1 |
| `C5_pathway_w10000` | `gap_kind` | `maskable_gap` | 195 | 1 | 1 | 193 | 2 |
| `C5_pathway_w10000` | `pathway_robust` | `marginal` | 122 | 0 | 1 | 121 | 1 |
| `C5_pathway_w10000` | `pathway_robust` | `robust` | 131 | 1 | 1 | 129 | 2 |
| `C5_organism_w10000` | `domain` | `archaea` | 6 | 0 | 0 | 6 | 0 |
| `C5_organism_w10000` | `domain` | `bacteria` | 247 | 1 | 5 | 241 | 6 |
| `C5_organism_w10000` | `gap_kind` | `already_gap` | 58 | 0 | 2 | 56 | 2 |
| `C5_organism_w10000` | `gap_kind` | `maskable_gap` | 195 | 1 | 3 | 191 | 4 |
| `C5_organism_w10000` | `pathway_robust` | `marginal` | 122 | 0 | 4 | 118 | 4 |
| `C5_organism_w10000` | `pathway_robust` | `robust` | 131 | 1 | 1 | 129 | 2 |
| `C5total_universal_w10000` | `domain` | `archaea` | 6 | 2 | 0 | 4 | 2 |
| `C5total_universal_w10000` | `domain` | `bacteria` | 247 | 39 | 2 | 206 | 41 |
| `C5total_universal_w10000` | `gap_kind` | `already_gap` | 58 | 9 | 0 | 49 | 9 |
| `C5total_universal_w10000` | `gap_kind` | `maskable_gap` | 195 | 32 | 2 | 161 | 34 |
| `C5total_universal_w10000` | `pathway_robust` | `marginal` | 122 | 22 | 1 | 99 | 23 |
| `C5total_universal_w10000` | `pathway_robust` | `robust` | 131 | 19 | 1 | 111 | 20 |
| `C5_pathway_vs_universal_w10000` | `domain` | `archaea` | 6 | 1 | 0 | 5 | 1 |
| `C5_pathway_vs_universal_w10000` | `domain` | `bacteria` | 247 | 3 | 0 | 244 | 3 |
| `C5_pathway_vs_universal_w10000` | `gap_kind` | `already_gap` | 58 | 1 | 0 | 57 | 1 |
| `C5_pathway_vs_universal_w10000` | `gap_kind` | `maskable_gap` | 195 | 3 | 0 | 192 | 3 |
| `C5_pathway_vs_universal_w10000` | `pathway_robust` | `marginal` | 122 | 3 | 0 | 119 | 3 |
| `C5_pathway_vs_universal_w10000` | `pathway_robust` | `robust` | 131 | 1 | 0 | 130 | 1 |
| `C5_organism_vs_universal_w10000` | `domain` | `archaea` | 6 | 1 | 0 | 5 | 1 |
| `C5_organism_vs_universal_w10000` | `domain` | `bacteria` | 247 | 0 | 0 | 247 | 0 |
| `C5_organism_vs_universal_w10000` | `gap_kind` | `already_gap` | 58 | 0 | 0 | 58 | 0 |
| `C5_organism_vs_universal_w10000` | `gap_kind` | `maskable_gap` | 195 | 1 | 0 | 194 | 1 |
| `C5_organism_vs_universal_w10000` | `pathway_robust` | `marginal` | 122 | 0 | 0 | 122 | 0 |
| `C5_organism_vs_universal_w10000` | `pathway_robust` | `robust` | 131 | 1 | 0 | 130 | 1 |
| `C5_organism_vs_pathway_w10000` | `domain` | `archaea` | 6 | 0 | 0 | 6 | 0 |
| `C5_organism_vs_pathway_w10000` | `domain` | `bacteria` | 247 | 0 | 3 | 244 | 3 |
| `C5_organism_vs_pathway_w10000` | `gap_kind` | `already_gap` | 58 | 0 | 1 | 57 | 1 |
| `C5_organism_vs_pathway_w10000` | `gap_kind` | `maskable_gap` | 195 | 0 | 2 | 193 | 2 |
| `C5_organism_vs_pathway_w10000` | `pathway_robust` | `marginal` | 122 | 0 | 3 | 119 | 3 |
| `C5_organism_vs_pathway_w10000` | `pathway_robust` | `robust` | 131 | 0 | 0 | 131 | 0 |
| `C5_universal_w20000` | `domain` | `archaea` | 6 | 0 | 1 | 5 | 1 |
| `C5_universal_w20000` | `domain` | `bacteria` | 247 | 1 | 3 | 243 | 4 |
| `C5_universal_w20000` | `gap_kind` | `already_gap` | 58 | 0 | 2 | 56 | 2 |
| `C5_universal_w20000` | `gap_kind` | `maskable_gap` | 195 | 1 | 2 | 192 | 3 |
| `C5_universal_w20000` | `pathway_robust` | `marginal` | 122 | 0 | 2 | 120 | 2 |
| `C5_universal_w20000` | `pathway_robust` | `robust` | 131 | 1 | 2 | 128 | 3 |
| `C5_pathway_w20000` | `domain` | `archaea` | 6 | 0 | 0 | 6 | 0 |
| `C5_pathway_w20000` | `domain` | `bacteria` | 247 | 1 | 1 | 245 | 2 |
| `C5_pathway_w20000` | `gap_kind` | `already_gap` | 58 | 0 | 1 | 57 | 1 |
| `C5_pathway_w20000` | `gap_kind` | `maskable_gap` | 195 | 1 | 0 | 194 | 1 |
| `C5_pathway_w20000` | `pathway_robust` | `marginal` | 122 | 0 | 0 | 122 | 0 |
| `C5_pathway_w20000` | `pathway_robust` | `robust` | 131 | 1 | 1 | 129 | 2 |
| `C5_organism_w20000` | `domain` | `archaea` | 6 | 0 | 0 | 6 | 0 |
| `C5_organism_w20000` | `domain` | `bacteria` | 247 | 1 | 4 | 242 | 5 |
| `C5_organism_w20000` | `gap_kind` | `already_gap` | 58 | 0 | 2 | 56 | 2 |
| `C5_organism_w20000` | `gap_kind` | `maskable_gap` | 195 | 1 | 2 | 192 | 3 |
| `C5_organism_w20000` | `pathway_robust` | `marginal` | 122 | 0 | 3 | 119 | 3 |
| `C5_organism_w20000` | `pathway_robust` | `robust` | 131 | 1 | 1 | 129 | 2 |
| `C5total_universal_w20000` | `domain` | `archaea` | 6 | 2 | 0 | 4 | 2 |
| `C5total_universal_w20000` | `domain` | `bacteria` | 247 | 41 | 2 | 204 | 43 |
| `C5total_universal_w20000` | `gap_kind` | `already_gap` | 58 | 9 | 0 | 49 | 9 |
| `C5total_universal_w20000` | `gap_kind` | `maskable_gap` | 195 | 34 | 2 | 159 | 36 |
| `C5total_universal_w20000` | `pathway_robust` | `marginal` | 122 | 24 | 1 | 97 | 25 |
| `C5total_universal_w20000` | `pathway_robust` | `robust` | 131 | 19 | 1 | 111 | 20 |
| `C5_pathway_vs_universal_w20000` | `domain` | `archaea` | 6 | 1 | 0 | 5 | 1 |
| `C5_pathway_vs_universal_w20000` | `domain` | `bacteria` | 247 | 2 | 0 | 245 | 2 |
| `C5_pathway_vs_universal_w20000` | `gap_kind` | `already_gap` | 58 | 1 | 0 | 57 | 1 |
| `C5_pathway_vs_universal_w20000` | `gap_kind` | `maskable_gap` | 195 | 2 | 0 | 193 | 2 |
| `C5_pathway_vs_universal_w20000` | `pathway_robust` | `marginal` | 122 | 2 | 0 | 120 | 2 |
| `C5_pathway_vs_universal_w20000` | `pathway_robust` | `robust` | 131 | 1 | 0 | 130 | 1 |
| `C5_organism_vs_universal_w20000` | `domain` | `archaea` | 6 | 1 | 0 | 5 | 1 |
| `C5_organism_vs_universal_w20000` | `domain` | `bacteria` | 247 | 0 | 1 | 246 | 1 |
| `C5_organism_vs_universal_w20000` | `gap_kind` | `already_gap` | 58 | 0 | 0 | 58 | 0 |
| `C5_organism_vs_universal_w20000` | `gap_kind` | `maskable_gap` | 195 | 1 | 1 | 193 | 2 |
| `C5_organism_vs_universal_w20000` | `pathway_robust` | `marginal` | 122 | 0 | 1 | 121 | 1 |
| `C5_organism_vs_universal_w20000` | `pathway_robust` | `robust` | 131 | 1 | 0 | 130 | 1 |
| `C5_organism_vs_pathway_w20000` | `domain` | `archaea` | 6 | 0 | 0 | 6 | 0 |
| `C5_organism_vs_pathway_w20000` | `domain` | `bacteria` | 247 | 0 | 3 | 244 | 3 |
| `C5_organism_vs_pathway_w20000` | `gap_kind` | `already_gap` | 58 | 0 | 1 | 57 | 1 |
| `C5_organism_vs_pathway_w20000` | `gap_kind` | `maskable_gap` | 195 | 0 | 2 | 193 | 2 |
| `C5_organism_vs_pathway_w20000` | `pathway_robust` | `marginal` | 122 | 0 | 3 | 119 | 3 |
| `C5_organism_vs_pathway_w20000` | `pathway_robust` | `robust` | 131 | 0 | 0 | 131 | 0 |
| `C5_universal` | `domain` | `archaea` | 6 | 0 | 1 | 5 | 1 |
| `C5_universal` | `domain` | `bacteria` | 247 | 1 | 5 | 241 | 6 |
| `C5_universal` | `gap_kind` | `already_gap` | 58 | 0 | 2 | 56 | 2 |
| `C5_universal` | `gap_kind` | `maskable_gap` | 195 | 1 | 4 | 190 | 5 |
| `C5_universal` | `pathway_robust` | `marginal` | 122 | 0 | 4 | 118 | 4 |
| `C5_universal` | `pathway_robust` | `robust` | 131 | 1 | 2 | 128 | 3 |
| `C5_pathway` | `domain` | `archaea` | 6 | 0 | 0 | 6 | 0 |
| `C5_pathway` | `domain` | `bacteria` | 247 | 1 | 2 | 244 | 3 |
| `C5_pathway` | `gap_kind` | `already_gap` | 58 | 0 | 1 | 57 | 1 |
| `C5_pathway` | `gap_kind` | `maskable_gap` | 195 | 1 | 1 | 193 | 2 |
| `C5_pathway` | `pathway_robust` | `marginal` | 122 | 0 | 1 | 121 | 1 |
| `C5_pathway` | `pathway_robust` | `robust` | 131 | 1 | 1 | 129 | 2 |
| `C5_organism` | `domain` | `archaea` | 6 | 0 | 0 | 6 | 0 |
| `C5_organism` | `domain` | `bacteria` | 247 | 1 | 5 | 241 | 6 |
| `C5_organism` | `gap_kind` | `already_gap` | 58 | 0 | 2 | 56 | 2 |
| `C5_organism` | `gap_kind` | `maskable_gap` | 195 | 1 | 3 | 191 | 4 |
| `C5_organism` | `pathway_robust` | `marginal` | 122 | 0 | 4 | 118 | 4 |
| `C5_organism` | `pathway_robust` | `robust` | 131 | 1 | 1 | 129 | 2 |
| `C5_pathway_vs_universal` | `domain` | `archaea` | 6 | 1 | 0 | 5 | 1 |
| `C5_pathway_vs_universal` | `domain` | `bacteria` | 247 | 3 | 0 | 244 | 3 |
| `C5_pathway_vs_universal` | `gap_kind` | `already_gap` | 58 | 1 | 0 | 57 | 1 |
| `C5_pathway_vs_universal` | `gap_kind` | `maskable_gap` | 195 | 3 | 0 | 192 | 3 |
| `C5_pathway_vs_universal` | `pathway_robust` | `marginal` | 122 | 3 | 0 | 119 | 3 |
| `C5_pathway_vs_universal` | `pathway_robust` | `robust` | 131 | 1 | 0 | 130 | 1 |
| `C5_organism_vs_universal` | `domain` | `archaea` | 6 | 1 | 0 | 5 | 1 |
| `C5_organism_vs_universal` | `domain` | `bacteria` | 247 | 0 | 0 | 247 | 0 |
| `C5_organism_vs_universal` | `gap_kind` | `already_gap` | 58 | 0 | 0 | 58 | 0 |
| `C5_organism_vs_universal` | `gap_kind` | `maskable_gap` | 195 | 1 | 0 | 194 | 1 |
| `C5_organism_vs_universal` | `pathway_robust` | `marginal` | 122 | 0 | 0 | 122 | 0 |
| `C5_organism_vs_universal` | `pathway_robust` | `robust` | 131 | 1 | 0 | 130 | 1 |
| `C5_organism_vs_pathway` | `domain` | `archaea` | 6 | 0 | 0 | 6 | 0 |
| `C5_organism_vs_pathway` | `domain` | `bacteria` | 247 | 0 | 3 | 244 | 3 |
| `C5_organism_vs_pathway` | `gap_kind` | `already_gap` | 58 | 0 | 1 | 57 | 1 |
| `C5_organism_vs_pathway` | `gap_kind` | `maskable_gap` | 195 | 0 | 2 | 193 | 2 |
| `C5_organism_vs_pathway` | `pathway_robust` | `marginal` | 122 | 0 | 3 | 119 | 3 |
| `C5_organism_vs_pathway` | `pathway_robust` | `robust` | 131 | 0 | 0 | 131 | 0 |
| `S1_stage0_to_sham_gaps` | `domain` | `archaea` | 6 | 0 | 0 | 6 | 0 |
| `S1_stage0_to_sham_gaps` | `domain` | `bacteria` | 247 | 0 | 0 | 247 | 0 |
| `S1_stage0_to_sham_gaps` | `gap_kind` | `already_gap` | 58 | 0 | 0 | 58 | 0 |
| `S1_stage0_to_sham_gaps` | `gap_kind` | `maskable_gap` | 195 | 0 | 0 | 195 | 0 |
| `S1_stage0_to_sham_gaps` | `pathway_robust` | `marginal` | 122 | 0 | 0 | 122 | 0 |
| `S1_stage0_to_sham_gaps` | `pathway_robust` | `robust` | 131 | 0 | 0 | 131 | 0 |
| `S1_sham_gaps_to_sham_gapsfix` | `domain` | `archaea` | 6 | 0 | 0 | 6 | 0 |
| `S1_sham_gaps_to_sham_gapsfix` | `domain` | `bacteria` | 247 | 11 | 0 | 236 | 11 |
| `S1_sham_gaps_to_sham_gapsfix` | `gap_kind` | `already_gap` | 58 | 11 | 0 | 47 | 11 |
| `S1_sham_gaps_to_sham_gapsfix` | `gap_kind` | `maskable_gap` | 195 | 0 | 0 | 195 | 0 |
| `S1_sham_gaps_to_sham_gapsfix` | `pathway_robust` | `marginal` | 122 | 1 | 0 | 121 | 1 |
| `S1_sham_gaps_to_sham_gapsfix` | `pathway_robust` | `robust` | 131 | 10 | 0 | 121 | 10 |
| `S1_sham_gaps_to_sham_ctx_universal_w5000` | `domain` | `archaea` | 6 | 0 | 0 | 6 | 0 |
| `S1_sham_gaps_to_sham_ctx_universal_w5000` | `domain` | `bacteria` | 247 | 0 | 2 | 245 | 2 |
| `S1_sham_gaps_to_sham_ctx_universal_w5000` | `gap_kind` | `already_gap` | 58 | 0 | 2 | 56 | 2 |
| `S1_sham_gaps_to_sham_ctx_universal_w5000` | `gap_kind` | `maskable_gap` | 195 | 0 | 0 | 195 | 0 |
| `S1_sham_gaps_to_sham_ctx_universal_w5000` | `pathway_robust` | `marginal` | 122 | 0 | 1 | 121 | 1 |
| `S1_sham_gaps_to_sham_ctx_universal_w5000` | `pathway_robust` | `robust` | 131 | 0 | 1 | 130 | 1 |
| `S1_stage0_to_sham_ctx_universal_w5000` | `domain` | `archaea` | 6 | 0 | 0 | 6 | 0 |
| `S1_stage0_to_sham_ctx_universal_w5000` | `domain` | `bacteria` | 247 | 9 | 0 | 238 | 9 |
| `S1_stage0_to_sham_ctx_universal_w5000` | `gap_kind` | `already_gap` | 58 | 9 | 0 | 49 | 9 |
| `S1_stage0_to_sham_ctx_universal_w5000` | `gap_kind` | `maskable_gap` | 195 | 0 | 0 | 195 | 0 |
| `S1_stage0_to_sham_ctx_universal_w5000` | `pathway_robust` | `marginal` | 122 | 0 | 0 | 122 | 0 |
| `S1_stage0_to_sham_ctx_universal_w5000` | `pathway_robust` | `robust` | 131 | 9 | 0 | 122 | 9 |
| `S1_sham_gaps_to_sham_ctx_universal_w10000` | `domain` | `archaea` | 6 | 0 | 0 | 6 | 0 |
| `S1_sham_gaps_to_sham_ctx_universal_w10000` | `domain` | `bacteria` | 247 | 0 | 2 | 245 | 2 |
| `S1_sham_gaps_to_sham_ctx_universal_w10000` | `gap_kind` | `already_gap` | 58 | 0 | 2 | 56 | 2 |
| `S1_sham_gaps_to_sham_ctx_universal_w10000` | `gap_kind` | `maskable_gap` | 195 | 0 | 0 | 195 | 0 |
| `S1_sham_gaps_to_sham_ctx_universal_w10000` | `pathway_robust` | `marginal` | 122 | 0 | 1 | 121 | 1 |
| `S1_sham_gaps_to_sham_ctx_universal_w10000` | `pathway_robust` | `robust` | 131 | 0 | 1 | 130 | 1 |
| `S1_stage0_to_sham_ctx_universal_w10000` | `domain` | `archaea` | 6 | 0 | 0 | 6 | 0 |
| `S1_stage0_to_sham_ctx_universal_w10000` | `domain` | `bacteria` | 247 | 9 | 0 | 238 | 9 |
| `S1_stage0_to_sham_ctx_universal_w10000` | `gap_kind` | `already_gap` | 58 | 9 | 0 | 49 | 9 |
| `S1_stage0_to_sham_ctx_universal_w10000` | `gap_kind` | `maskable_gap` | 195 | 0 | 0 | 195 | 0 |
| `S1_stage0_to_sham_ctx_universal_w10000` | `pathway_robust` | `marginal` | 122 | 0 | 0 | 122 | 0 |
| `S1_stage0_to_sham_ctx_universal_w10000` | `pathway_robust` | `robust` | 131 | 9 | 0 | 122 | 9 |
| `S1_sham_gaps_to_sham_ctx_universal_w20000` | `domain` | `archaea` | 6 | 0 | 0 | 6 | 0 |
| `S1_sham_gaps_to_sham_ctx_universal_w20000` | `domain` | `bacteria` | 247 | 0 | 2 | 245 | 2 |
| `S1_sham_gaps_to_sham_ctx_universal_w20000` | `gap_kind` | `already_gap` | 58 | 0 | 2 | 56 | 2 |
| `S1_sham_gaps_to_sham_ctx_universal_w20000` | `gap_kind` | `maskable_gap` | 195 | 0 | 0 | 195 | 0 |
| `S1_sham_gaps_to_sham_ctx_universal_w20000` | `pathway_robust` | `marginal` | 122 | 0 | 1 | 121 | 1 |
| `S1_sham_gaps_to_sham_ctx_universal_w20000` | `pathway_robust` | `robust` | 131 | 0 | 1 | 130 | 1 |
| `S1_stage0_to_sham_ctx_universal_w20000` | `domain` | `archaea` | 6 | 0 | 0 | 6 | 0 |
| `S1_stage0_to_sham_ctx_universal_w20000` | `domain` | `bacteria` | 247 | 9 | 0 | 238 | 9 |
| `S1_stage0_to_sham_ctx_universal_w20000` | `gap_kind` | `already_gap` | 58 | 9 | 0 | 49 | 9 |
| `S1_stage0_to_sham_ctx_universal_w20000` | `gap_kind` | `maskable_gap` | 195 | 0 | 0 | 195 | 0 |
| `S1_stage0_to_sham_ctx_universal_w20000` | `pathway_robust` | `marginal` | 122 | 0 | 0 | 122 | 0 |
| `S1_stage0_to_sham_ctx_universal_w20000` | `pathway_robust` | `robust` | 131 | 9 | 0 | 122 | 9 |
| `S1_sham_gaps_to_sham_ctx_universal` | `domain` | `archaea` | 6 | 0 | 0 | 6 | 0 |
| `S1_sham_gaps_to_sham_ctx_universal` | `domain` | `bacteria` | 247 | 0 | 2 | 245 | 2 |
| `S1_sham_gaps_to_sham_ctx_universal` | `gap_kind` | `already_gap` | 58 | 0 | 2 | 56 | 2 |
| `S1_sham_gaps_to_sham_ctx_universal` | `gap_kind` | `maskable_gap` | 195 | 0 | 0 | 195 | 0 |
| `S1_sham_gaps_to_sham_ctx_universal` | `pathway_robust` | `marginal` | 122 | 0 | 1 | 121 | 1 |
| `S1_sham_gaps_to_sham_ctx_universal` | `pathway_robust` | `robust` | 131 | 0 | 1 | 130 | 1 |
| `U1_unmasked_gaps_vs_gaps` | `domain` | `archaea` | 6 | 0 | 0 | 6 | 0 |
| `U1_unmasked_gaps_vs_gaps` | `domain` | `bacteria` | 247 | 0 | 0 | 247 | 0 |
| `U1_unmasked_gaps_vs_gaps` | `gap_kind` | `already_gap` | 58 | 0 | 0 | 58 | 0 |
| `U1_unmasked_gaps_vs_gaps` | `gap_kind` | `maskable_gap` | 195 | 0 | 0 | 195 | 0 |
| `U1_unmasked_gaps_vs_gaps` | `pathway_robust` | `marginal` | 122 | 0 | 0 | 122 | 0 |
| `U1_unmasked_gaps_vs_gaps` | `pathway_robust` | `robust` | 131 | 0 | 0 | 131 | 0 |

## 4. Learned-Stage Candidate Counts

| Arm | Candidates | Rescored | Unsupported | Unchanged unsupported | Retracted | Scored only |
|---|---:|---:|---:|---:|---:|---:|
| - | 0 | 0 | 0 | 0 | 0 | 0 |

## 5. Masking & Build Audit Summary

| Check | Count |
|---|---:|
| `total_cells` | 253 |
| `complete_cells` | 253 |
| `candidate_universe_preserved_cells` | 253 |
| `masked_gap_present_cells` | 253 |
| `gaps_ctxjar_identical_cells` | 253 |
| `already_gap_unmasked_total` | 58 |
| `already_gap_unmasked_matches` | 58 |
