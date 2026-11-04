# Row-5 Masked-Context Evaluation (`scope=all`)

Chosen context window (bp): `10000`

## 1. Contrasts (Overall and Per Claim Set)

| Contrast | Claim set | Cells | Improved | Harmed | Unchanged | Discordant | Critical k | p (improve) | Harm elsewhere | Verdict |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `C1_stage0_to_gaps` | **all** | 23 | 0 | 0 | 23 | 0 | - | 1 | 0 | no detectable effect |
| `C1_stage0_to_gaps` | `clean` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C1_stage0_to_gaps` | `diamond` | 5 | 0 | 0 | 5 | 0 | - | 1 | 0 | no detectable effect |
| `C1_stage0_to_gaps` | `dpz` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C1_stage0_to_gaps` | `enzbert` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C1_stage0_to_gaps` | `graphec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C1_stage0_to_gaps` | `mapred` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C1_stage0_to_gaps` | `topec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C2_gaps_to_gaps_learned` | **all** | 23 | 0 | 0 | 23 | 0 | - | 1 | 0 | no detectable effect |
| `C2_gaps_to_gaps_learned` | `clean` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C2_gaps_to_gaps_learned` | `diamond` | 5 | 0 | 0 | 5 | 0 | - | 1 | 0 | no detectable effect |
| `C2_gaps_to_gaps_learned` | `dpz` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C2_gaps_to_gaps_learned` | `enzbert` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C2_gaps_to_gaps_learned` | `graphec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C2_gaps_to_gaps_learned` | `mapred` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C2_gaps_to_gaps_learned` | `topec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C3_universal_vs_gaps_learned` | **all** | 23 | 0 | 0 | 23 | 0 | - | 1 | 0 | no detectable effect |
| `C3_universal_vs_gaps_learned` | `clean` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C3_universal_vs_gaps_learned` | `diamond` | 5 | 0 | 0 | 5 | 0 | - | 1 | 0 | no detectable effect |
| `C3_universal_vs_gaps_learned` | `dpz` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C3_universal_vs_gaps_learned` | `enzbert` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C3_universal_vs_gaps_learned` | `graphec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C3_universal_vs_gaps_learned` | `mapred` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C3_universal_vs_gaps_learned` | `topec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C3_pathway_vs_gaps_learned` | **all** | 23 | 0 | 0 | 23 | 0 | - | 1 | 0 | no detectable effect |
| `C3_pathway_vs_gaps_learned` | `clean` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C3_pathway_vs_gaps_learned` | `diamond` | 5 | 0 | 0 | 5 | 0 | - | 1 | 0 | no detectable effect |
| `C3_pathway_vs_gaps_learned` | `dpz` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C3_pathway_vs_gaps_learned` | `enzbert` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C3_pathway_vs_gaps_learned` | `graphec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C3_pathway_vs_gaps_learned` | `mapred` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C3_pathway_vs_gaps_learned` | `topec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C3_organism_vs_gaps_learned` | **all** | 23 | 0 | 0 | 23 | 0 | - | 1 | 0 | no detectable effect |
| `C3_organism_vs_gaps_learned` | `clean` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C3_organism_vs_gaps_learned` | `diamond` | 5 | 0 | 0 | 5 | 0 | - | 1 | 0 | no detectable effect |
| `C3_organism_vs_gaps_learned` | `dpz` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C3_organism_vs_gaps_learned` | `enzbert` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C3_organism_vs_gaps_learned` | `graphec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C3_organism_vs_gaps_learned` | `mapred` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C3_organism_vs_gaps_learned` | `topec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C4_pathway_vs_universal` | **all** | 23 | 0 | 0 | 23 | 0 | - | 1 | 0 | no detectable effect |
| `C4_pathway_vs_universal` | `clean` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C4_pathway_vs_universal` | `diamond` | 5 | 0 | 0 | 5 | 0 | - | 1 | 0 | no detectable effect |
| `C4_pathway_vs_universal` | `dpz` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C4_pathway_vs_universal` | `enzbert` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C4_pathway_vs_universal` | `graphec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C4_pathway_vs_universal` | `mapred` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C4_pathway_vs_universal` | `topec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C4_organism_vs_universal` | **all** | 23 | 0 | 0 | 23 | 0 | - | 1 | 0 | no detectable effect |
| `C4_organism_vs_universal` | `clean` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C4_organism_vs_universal` | `diamond` | 5 | 0 | 0 | 5 | 0 | - | 1 | 0 | no detectable effect |
| `C4_organism_vs_universal` | `dpz` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C4_organism_vs_universal` | `enzbert` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C4_organism_vs_universal` | `graphec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C4_organism_vs_universal` | `mapred` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C4_organism_vs_universal` | `topec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C4_organism_vs_pathway` | **all** | 23 | 0 | 0 | 23 | 0 | - | 1 | 0 | no detectable effect |
| `C4_organism_vs_pathway` | `clean` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C4_organism_vs_pathway` | `diamond` | 5 | 0 | 0 | 5 | 0 | - | 1 | 0 | no detectable effect |
| `C4_organism_vs_pathway` | `dpz` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C4_organism_vs_pathway` | `enzbert` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C4_organism_vs_pathway` | `graphec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C4_organism_vs_pathway` | `mapred` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C4_organism_vs_pathway` | `topec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal_w5000` | **all** | 23 | 3 | 0 | 20 | 3 | - | 0.125 | 0 | no detectable effect |
| `C5_universal_w5000` | `clean` | 3 | 1 | 0 | 2 | 1 | - | 0.5 | 0 | no detectable effect |
| `C5_universal_w5000` | `diamond` | 5 | 1 | 0 | 4 | 1 | - | 0.5 | 0 | no detectable effect |
| `C5_universal_w5000` | `dpz` | 3 | 1 | 0 | 2 | 1 | - | 0.5 | 0 | no detectable effect |
| `C5_universal_w5000` | `enzbert` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal_w5000` | `graphec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal_w5000` | `mapred` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal_w5000` | `topec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_w5000` | **all** | 23 | 3 | 0 | 20 | 3 | - | 0.125 | 0 | no detectable effect |
| `C5_pathway_w5000` | `clean` | 3 | 1 | 0 | 2 | 1 | - | 0.5 | 0 | no detectable effect |
| `C5_pathway_w5000` | `diamond` | 5 | 1 | 0 | 4 | 1 | - | 0.5 | 0 | no detectable effect |
| `C5_pathway_w5000` | `dpz` | 3 | 1 | 0 | 2 | 1 | - | 0.5 | 0 | no detectable effect |
| `C5_pathway_w5000` | `enzbert` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_w5000` | `graphec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_w5000` | `mapred` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_w5000` | `topec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_w5000` | **all** | 23 | 3 | 0 | 20 | 3 | - | 0.125 | 0 | no detectable effect |
| `C5_organism_w5000` | `clean` | 3 | 1 | 0 | 2 | 1 | - | 0.5 | 0 | no detectable effect |
| `C5_organism_w5000` | `diamond` | 5 | 1 | 0 | 4 | 1 | - | 0.5 | 0 | no detectable effect |
| `C5_organism_w5000` | `dpz` | 3 | 1 | 0 | 2 | 1 | - | 0.5 | 0 | no detectable effect |
| `C5_organism_w5000` | `enzbert` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_w5000` | `graphec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_w5000` | `mapred` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_w5000` | `topec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5total_universal_w5000` | **all** | 23 | 3 | 0 | 20 | 3 | - | 0.125 | 0 | no detectable effect |
| `C5total_universal_w5000` | `clean` | 3 | 1 | 0 | 2 | 1 | - | 0.5 | 0 | no detectable effect |
| `C5total_universal_w5000` | `diamond` | 5 | 1 | 0 | 4 | 1 | - | 0.5 | 0 | no detectable effect |
| `C5total_universal_w5000` | `dpz` | 3 | 1 | 0 | 2 | 1 | - | 0.5 | 0 | no detectable effect |
| `C5total_universal_w5000` | `enzbert` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5total_universal_w5000` | `graphec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5total_universal_w5000` | `mapred` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5total_universal_w5000` | `topec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w5000` | **all** | 23 | 0 | 0 | 23 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w5000` | `clean` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w5000` | `diamond` | 5 | 0 | 0 | 5 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w5000` | `dpz` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w5000` | `enzbert` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w5000` | `graphec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w5000` | `mapred` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w5000` | `topec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w5000` | **all** | 23 | 0 | 0 | 23 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w5000` | `clean` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w5000` | `diamond` | 5 | 0 | 0 | 5 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w5000` | `dpz` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w5000` | `enzbert` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w5000` | `graphec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w5000` | `mapred` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w5000` | `topec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway_w5000` | **all** | 23 | 0 | 0 | 23 | 0 | - | 1 | 1 | no detectable effect |
| `C5_organism_vs_pathway_w5000` | `clean` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway_w5000` | `diamond` | 5 | 0 | 0 | 5 | 0 | - | 1 | 1 | no detectable effect |
| `C5_organism_vs_pathway_w5000` | `dpz` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway_w5000` | `enzbert` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway_w5000` | `graphec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway_w5000` | `mapred` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway_w5000` | `topec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal_w10000` | **all** | 23 | 3 | 0 | 20 | 3 | - | 0.125 | 0 | no detectable effect |
| `C5_universal_w10000` | `clean` | 3 | 1 | 0 | 2 | 1 | - | 0.5 | 0 | no detectable effect |
| `C5_universal_w10000` | `diamond` | 5 | 1 | 0 | 4 | 1 | - | 0.5 | 0 | no detectable effect |
| `C5_universal_w10000` | `dpz` | 3 | 1 | 0 | 2 | 1 | - | 0.5 | 0 | no detectable effect |
| `C5_universal_w10000` | `enzbert` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal_w10000` | `graphec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal_w10000` | `mapred` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal_w10000` | `topec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_w10000` | **all** | 23 | 3 | 0 | 20 | 3 | - | 0.125 | 0 | no detectable effect |
| `C5_pathway_w10000` | `clean` | 3 | 1 | 0 | 2 | 1 | - | 0.5 | 0 | no detectable effect |
| `C5_pathway_w10000` | `diamond` | 5 | 1 | 0 | 4 | 1 | - | 0.5 | 0 | no detectable effect |
| `C5_pathway_w10000` | `dpz` | 3 | 1 | 0 | 2 | 1 | - | 0.5 | 0 | no detectable effect |
| `C5_pathway_w10000` | `enzbert` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_w10000` | `graphec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_w10000` | `mapred` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_w10000` | `topec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_w10000` | **all** | 23 | 3 | 0 | 20 | 3 | - | 0.125 | 0 | no detectable effect |
| `C5_organism_w10000` | `clean` | 3 | 1 | 0 | 2 | 1 | - | 0.5 | 0 | no detectable effect |
| `C5_organism_w10000` | `diamond` | 5 | 1 | 0 | 4 | 1 | - | 0.5 | 0 | no detectable effect |
| `C5_organism_w10000` | `dpz` | 3 | 1 | 0 | 2 | 1 | - | 0.5 | 0 | no detectable effect |
| `C5_organism_w10000` | `enzbert` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_w10000` | `graphec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_w10000` | `mapred` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_w10000` | `topec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5total_universal_w10000` | **all** | 23 | 3 | 0 | 20 | 3 | - | 0.125 | 0 | no detectable effect |
| `C5total_universal_w10000` | `clean` | 3 | 1 | 0 | 2 | 1 | - | 0.5 | 0 | no detectable effect |
| `C5total_universal_w10000` | `diamond` | 5 | 1 | 0 | 4 | 1 | - | 0.5 | 0 | no detectable effect |
| `C5total_universal_w10000` | `dpz` | 3 | 1 | 0 | 2 | 1 | - | 0.5 | 0 | no detectable effect |
| `C5total_universal_w10000` | `enzbert` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5total_universal_w10000` | `graphec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5total_universal_w10000` | `mapred` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5total_universal_w10000` | `topec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w10000` | **all** | 23 | 0 | 0 | 23 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w10000` | `clean` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w10000` | `diamond` | 5 | 0 | 0 | 5 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w10000` | `dpz` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w10000` | `enzbert` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w10000` | `graphec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w10000` | `mapred` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w10000` | `topec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w10000` | **all** | 23 | 0 | 0 | 23 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w10000` | `clean` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w10000` | `diamond` | 5 | 0 | 0 | 5 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w10000` | `dpz` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w10000` | `enzbert` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w10000` | `graphec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w10000` | `mapred` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w10000` | `topec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway_w10000` | **all** | 23 | 0 | 0 | 23 | 0 | - | 1 | 1 | no detectable effect |
| `C5_organism_vs_pathway_w10000` | `clean` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway_w10000` | `diamond` | 5 | 0 | 0 | 5 | 0 | - | 1 | 1 | no detectable effect |
| `C5_organism_vs_pathway_w10000` | `dpz` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway_w10000` | `enzbert` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway_w10000` | `graphec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway_w10000` | `mapred` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway_w10000` | `topec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal_w20000` | **all** | 23 | 3 | 0 | 20 | 3 | - | 0.125 | 0 | no detectable effect |
| `C5_universal_w20000` | `clean` | 3 | 1 | 0 | 2 | 1 | - | 0.5 | 0 | no detectable effect |
| `C5_universal_w20000` | `diamond` | 5 | 1 | 0 | 4 | 1 | - | 0.5 | 0 | no detectable effect |
| `C5_universal_w20000` | `dpz` | 3 | 1 | 0 | 2 | 1 | - | 0.5 | 0 | no detectable effect |
| `C5_universal_w20000` | `enzbert` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal_w20000` | `graphec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal_w20000` | `mapred` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal_w20000` | `topec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_w20000` | **all** | 23 | 3 | 0 | 20 | 3 | - | 0.125 | 0 | no detectable effect |
| `C5_pathway_w20000` | `clean` | 3 | 1 | 0 | 2 | 1 | - | 0.5 | 0 | no detectable effect |
| `C5_pathway_w20000` | `diamond` | 5 | 1 | 0 | 4 | 1 | - | 0.5 | 0 | no detectable effect |
| `C5_pathway_w20000` | `dpz` | 3 | 1 | 0 | 2 | 1 | - | 0.5 | 0 | no detectable effect |
| `C5_pathway_w20000` | `enzbert` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_w20000` | `graphec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_w20000` | `mapred` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_w20000` | `topec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_w20000` | **all** | 23 | 3 | 0 | 20 | 3 | - | 0.125 | 0 | no detectable effect |
| `C5_organism_w20000` | `clean` | 3 | 1 | 0 | 2 | 1 | - | 0.5 | 0 | no detectable effect |
| `C5_organism_w20000` | `diamond` | 5 | 1 | 0 | 4 | 1 | - | 0.5 | 0 | no detectable effect |
| `C5_organism_w20000` | `dpz` | 3 | 1 | 0 | 2 | 1 | - | 0.5 | 0 | no detectable effect |
| `C5_organism_w20000` | `enzbert` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_w20000` | `graphec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_w20000` | `mapred` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_w20000` | `topec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5total_universal_w20000` | **all** | 23 | 3 | 0 | 20 | 3 | - | 0.125 | 0 | no detectable effect |
| `C5total_universal_w20000` | `clean` | 3 | 1 | 0 | 2 | 1 | - | 0.5 | 0 | no detectable effect |
| `C5total_universal_w20000` | `diamond` | 5 | 1 | 0 | 4 | 1 | - | 0.5 | 0 | no detectable effect |
| `C5total_universal_w20000` | `dpz` | 3 | 1 | 0 | 2 | 1 | - | 0.5 | 0 | no detectable effect |
| `C5total_universal_w20000` | `enzbert` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5total_universal_w20000` | `graphec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5total_universal_w20000` | `mapred` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5total_universal_w20000` | `topec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w20000` | **all** | 23 | 0 | 0 | 23 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w20000` | `clean` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w20000` | `diamond` | 5 | 0 | 0 | 5 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w20000` | `dpz` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w20000` | `enzbert` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w20000` | `graphec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w20000` | `mapred` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal_w20000` | `topec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w20000` | **all** | 23 | 0 | 0 | 23 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w20000` | `clean` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w20000` | `diamond` | 5 | 0 | 0 | 5 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w20000` | `dpz` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w20000` | `enzbert` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w20000` | `graphec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w20000` | `mapred` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal_w20000` | `topec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway_w20000` | **all** | 23 | 0 | 0 | 23 | 0 | - | 1 | 1 | no detectable effect |
| `C5_organism_vs_pathway_w20000` | `clean` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway_w20000` | `diamond` | 5 | 0 | 0 | 5 | 0 | - | 1 | 1 | no detectable effect |
| `C5_organism_vs_pathway_w20000` | `dpz` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway_w20000` | `enzbert` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway_w20000` | `graphec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway_w20000` | `mapred` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway_w20000` | `topec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal` | **all** | 23 | 3 | 0 | 20 | 3 | - | 0.125 | 0 | no detectable effect |
| `C5_universal` | `clean` | 3 | 1 | 0 | 2 | 1 | - | 0.5 | 0 | no detectable effect |
| `C5_universal` | `diamond` | 5 | 1 | 0 | 4 | 1 | - | 0.5 | 0 | no detectable effect |
| `C5_universal` | `dpz` | 3 | 1 | 0 | 2 | 1 | - | 0.5 | 0 | no detectable effect |
| `C5_universal` | `enzbert` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal` | `graphec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal` | `mapred` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_universal` | `topec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway` | **all** | 23 | 3 | 0 | 20 | 3 | - | 0.125 | 0 | no detectable effect |
| `C5_pathway` | `clean` | 3 | 1 | 0 | 2 | 1 | - | 0.5 | 0 | no detectable effect |
| `C5_pathway` | `diamond` | 5 | 1 | 0 | 4 | 1 | - | 0.5 | 0 | no detectable effect |
| `C5_pathway` | `dpz` | 3 | 1 | 0 | 2 | 1 | - | 0.5 | 0 | no detectable effect |
| `C5_pathway` | `enzbert` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway` | `graphec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway` | `mapred` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway` | `topec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism` | **all** | 23 | 3 | 0 | 20 | 3 | - | 0.125 | 0 | no detectable effect |
| `C5_organism` | `clean` | 3 | 1 | 0 | 2 | 1 | - | 0.5 | 0 | no detectable effect |
| `C5_organism` | `diamond` | 5 | 1 | 0 | 4 | 1 | - | 0.5 | 0 | no detectable effect |
| `C5_organism` | `dpz` | 3 | 1 | 0 | 2 | 1 | - | 0.5 | 0 | no detectable effect |
| `C5_organism` | `enzbert` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism` | `graphec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism` | `mapred` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism` | `topec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal` | **all** | 23 | 0 | 0 | 23 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal` | `clean` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal` | `diamond` | 5 | 0 | 0 | 5 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal` | `dpz` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal` | `enzbert` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal` | `graphec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal` | `mapred` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_pathway_vs_universal` | `topec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal` | **all** | 23 | 0 | 0 | 23 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal` | `clean` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal` | `diamond` | 5 | 0 | 0 | 5 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal` | `dpz` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal` | `enzbert` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal` | `graphec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal` | `mapred` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_universal` | `topec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway` | **all** | 23 | 0 | 0 | 23 | 0 | - | 1 | 1 | no detectable effect |
| `C5_organism_vs_pathway` | `clean` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway` | `diamond` | 5 | 0 | 0 | 5 | 0 | - | 1 | 1 | no detectable effect |
| `C5_organism_vs_pathway` | `dpz` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway` | `enzbert` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway` | `graphec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway` | `mapred` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `C5_organism_vs_pathway` | `topec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_gaps` | **all** | 23 | 0 | 0 | 23 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_gaps` | `clean` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_gaps` | `diamond` | 5 | 0 | 0 | 5 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_gaps` | `dpz` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_gaps` | `enzbert` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_gaps` | `graphec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_gaps` | `mapred` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_gaps` | `topec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w5000` | **all** | 23 | 1 | 0 | 22 | 1 | - | 0.5 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w5000` | `clean` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w5000` | `diamond` | 5 | 0 | 0 | 5 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w5000` | `dpz` | 3 | 1 | 0 | 2 | 1 | - | 0.5 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w5000` | `enzbert` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w5000` | `graphec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w5000` | `mapred` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w5000` | `topec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w5000` | **all** | 23 | 1 | 0 | 22 | 1 | - | 0.5 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w5000` | `clean` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w5000` | `diamond` | 5 | 0 | 0 | 5 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w5000` | `dpz` | 3 | 1 | 0 | 2 | 1 | - | 0.5 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w5000` | `enzbert` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w5000` | `graphec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w5000` | `mapred` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w5000` | `topec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w10000` | **all** | 23 | 1 | 0 | 22 | 1 | - | 0.5 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w10000` | `clean` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w10000` | `diamond` | 5 | 0 | 0 | 5 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w10000` | `dpz` | 3 | 1 | 0 | 2 | 1 | - | 0.5 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w10000` | `enzbert` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w10000` | `graphec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w10000` | `mapred` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w10000` | `topec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w10000` | **all** | 23 | 1 | 0 | 22 | 1 | - | 0.5 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w10000` | `clean` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w10000` | `diamond` | 5 | 0 | 0 | 5 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w10000` | `dpz` | 3 | 1 | 0 | 2 | 1 | - | 0.5 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w10000` | `enzbert` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w10000` | `graphec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w10000` | `mapred` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w10000` | `topec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w20000` | **all** | 23 | 1 | 0 | 22 | 1 | - | 0.5 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w20000` | `clean` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w20000` | `diamond` | 5 | 0 | 0 | 5 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w20000` | `dpz` | 3 | 1 | 0 | 2 | 1 | - | 0.5 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w20000` | `enzbert` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w20000` | `graphec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w20000` | `mapred` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal_w20000` | `topec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w20000` | **all** | 23 | 1 | 0 | 22 | 1 | - | 0.5 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w20000` | `clean` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w20000` | `diamond` | 5 | 0 | 0 | 5 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w20000` | `dpz` | 3 | 1 | 0 | 2 | 1 | - | 0.5 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w20000` | `enzbert` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w20000` | `graphec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w20000` | `mapred` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `S1_stage0_to_sham_ctx_universal_w20000` | `topec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal` | **all** | 23 | 1 | 0 | 22 | 1 | - | 0.5 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal` | `clean` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal` | `diamond` | 5 | 0 | 0 | 5 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal` | `dpz` | 3 | 1 | 0 | 2 | 1 | - | 0.5 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal` | `enzbert` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal` | `graphec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal` | `mapred` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `S1_sham_gaps_to_sham_ctx_universal` | `topec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `U1_unmasked_gaps_vs_gaps` | **all** | 23 | 0 | 0 | 23 | 0 | - | 1 | 0 | no detectable effect |
| `U1_unmasked_gaps_vs_gaps` | `clean` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `U1_unmasked_gaps_vs_gaps` | `diamond` | 5 | 0 | 0 | 5 | 0 | - | 1 | 0 | no detectable effect |
| `U1_unmasked_gaps_vs_gaps` | `dpz` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `U1_unmasked_gaps_vs_gaps` | `enzbert` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `U1_unmasked_gaps_vs_gaps` | `graphec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `U1_unmasked_gaps_vs_gaps` | `mapred` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |
| `U1_unmasked_gaps_vs_gaps` | `topec` | 3 | 0 | 0 | 3 | 0 | - | 1 | 0 | no detectable effect |

## 2. Harm Elsewhere and Cost

| Contrast | Improved | Harmed | Harm elsewhere | Unconfirmed calls created | Delta calls >= 0.5 | Changed scores | Bootstrap 95% CI | Wilcoxon p (1-sided) |
|---|---:|---:|---:|---:|---:|---:|---|---:|
| `C1_stage0_to_gaps` | 0 | 0 | 0 | 0 | 621 | 4177 | [0.0, 0.0] | 1 |
| `C2_gaps_to_gaps_learned` | 0 | 0 | 0 | 0 | 0 | 0 | [0.0, 0.0] | 1 |
| `C3_universal_vs_gaps_learned` | 0 | 0 | 0 | 0 | 0 | 0 | [0.0, 0.0] | 1 |
| `C3_pathway_vs_gaps_learned` | 0 | 0 | 0 | 0 | 0 | 0 | [0.0, 0.0] | 1 |
| `C3_organism_vs_gaps_learned` | 0 | 0 | 0 | 0 | 0 | 0 | [0.0, 0.0] | 1 |
| `C4_pathway_vs_universal` | 0 | 0 | 0 | 0 | 0 | 0 | [0.0, 0.0] | 1 |
| `C4_organism_vs_universal` | 0 | 0 | 0 | 0 | 0 | 0 | [0.0, 0.0] | 1 |
| `C4_organism_vs_pathway` | 0 | 0 | 0 | 0 | 0 | 0 | [0.0, 0.0] | 1 |
| `C5_universal_w5000` | 3 | 0 | 0 | 2 | 157 | 532 | [1.0, 5.0] | 1 |
| `C5_pathway_w5000` | 3 | 0 | 0 | 2 | 161 | 531 | [1.0, 5.0] | 1 |
| `C5_organism_w5000` | 3 | 0 | 0 | 2 | 174 | 527 | [1.0, 5.0] | 1 |
| `C5total_universal_w5000` | 3 | 0 | 0 | 2 | 157 | 532 | [1.0, 5.0] | 1 |
| `C5_pathway_vs_universal_w5000` | 0 | 0 | 0 | 0 | 4 | 40 | [0.0, 0.0] | 1 |
| `C5_organism_vs_universal_w5000` | 0 | 0 | 0 | 0 | 17 | 36 | [0.0, 0.0] | 1 |
| `C5_organism_vs_pathway_w5000` | 0 | 0 | 1 | 0 | 13 | 47 | [0.0, 0.0] | 1 |
| `C5_universal_w10000` | 3 | 0 | 0 | 2 | 158 | 532 | [1.0, 5.0] | 1 |
| `C5_pathway_w10000` | 3 | 0 | 0 | 2 | 162 | 531 | [1.0, 5.0] | 1 |
| `C5_organism_w10000` | 3 | 0 | 0 | 2 | 175 | 527 | [1.0, 5.0] | 1 |
| `C5total_universal_w10000` | 3 | 0 | 0 | 2 | 158 | 532 | [1.0, 5.0] | 1 |
| `C5_pathway_vs_universal_w10000` | 0 | 0 | 0 | 0 | 4 | 42 | [0.0, 0.0] | 1 |
| `C5_organism_vs_universal_w10000` | 0 | 0 | 0 | 0 | 17 | 39 | [0.0, 0.0] | 1 |
| `C5_organism_vs_pathway_w10000` | 0 | 0 | 1 | 0 | 13 | 46 | [0.0, 0.0] | 1 |
| `C5_universal_w20000` | 3 | 0 | 0 | 2 | 160 | 533 | [1.0, 5.0] | 1 |
| `C5_pathway_w20000` | 3 | 0 | 0 | 2 | 163 | 531 | [1.0, 5.0] | 1 |
| `C5_organism_w20000` | 3 | 0 | 0 | 2 | 175 | 528 | [1.0, 5.0] | 1 |
| `C5total_universal_w20000` | 3 | 0 | 0 | 2 | 160 | 533 | [1.0, 5.0] | 1 |
| `C5_pathway_vs_universal_w20000` | 0 | 0 | 0 | 0 | 3 | 49 | [0.0, 0.0] | 1 |
| `C5_organism_vs_universal_w20000` | 0 | 0 | 0 | 0 | 15 | 43 | [0.0, 0.0] | 1 |
| `C5_organism_vs_pathway_w20000` | 0 | 0 | 1 | 0 | 12 | 51 | [0.0, 0.0] | 1 |
| `C5_universal` | 3 | 0 | 0 | 2 | 158 | 532 | [1.0, 5.0] | 1 |
| `C5_pathway` | 3 | 0 | 0 | 2 | 162 | 531 | [1.0, 5.0] | 1 |
| `C5_organism` | 3 | 0 | 0 | 2 | 175 | 527 | [1.0, 5.0] | 1 |
| `C5_pathway_vs_universal` | 0 | 0 | 0 | 0 | 4 | 42 | [0.0, 0.0] | 1 |
| `C5_organism_vs_universal` | 0 | 0 | 0 | 0 | 17 | 39 | [0.0, 0.0] | 1 |
| `C5_organism_vs_pathway` | 0 | 0 | 1 | 0 | 13 | 46 | [0.0, 0.0] | 1 |
| `S1_stage0_to_sham_gaps` | 0 | 0 | 0 | 0 | 620 | 4176 | [0.0, 0.0] | 1 |
| `S1_sham_gaps_to_sham_ctx_universal_w5000` | 1 | 0 | 0 | 2 | 157 | 531 | [0.0, 3.0] | 1 |
| `S1_stage0_to_sham_ctx_universal_w5000` | 1 | 0 | 0 | 2 | 777 | 4398 | [0.0, 3.0] | 1 |
| `S1_sham_gaps_to_sham_ctx_universal_w10000` | 1 | 0 | 0 | 2 | 158 | 531 | [0.0, 3.0] | 1 |
| `S1_stage0_to_sham_ctx_universal_w10000` | 1 | 0 | 0 | 2 | 778 | 4398 | [0.0, 3.0] | 1 |
| `S1_sham_gaps_to_sham_ctx_universal_w20000` | 1 | 0 | 0 | 2 | 160 | 531 | [0.0, 3.0] | 1 |
| `S1_stage0_to_sham_ctx_universal_w20000` | 1 | 0 | 0 | 2 | 780 | 4398 | [0.0, 3.0] | 1 |
| `S1_sham_gaps_to_sham_ctx_universal` | 1 | 0 | 0 | 2 | 158 | 531 | [0.0, 3.0] | 1 |
| `U1_unmasked_gaps_vs_gaps` | 0 | 0 | 0 | 0 | 2 | 2 | [0.0, 0.0] | 1 |

## 3. Strata Breakdown

| Contrast | Stratum axis | Stratum | Cells | Improved | Harmed | Unchanged | Discordant |
|---|---|---|---:|---:|---:|---:|---:|
| `C1_stage0_to_gaps` | `domain` | `bacteria` | 23 | 0 | 0 | 23 | 0 |
| `C1_stage0_to_gaps` | `gap_kind` | `already_gap` | 7 | 0 | 0 | 7 | 0 |
| `C1_stage0_to_gaps` | `gap_kind` | `maskable_gap` | 16 | 0 | 0 | 16 | 0 |
| `C1_stage0_to_gaps` | `pathway_robust` | `marginal` | 7 | 0 | 0 | 7 | 0 |
| `C1_stage0_to_gaps` | `pathway_robust` | `robust` | 16 | 0 | 0 | 16 | 0 |
| `C2_gaps_to_gaps_learned` | `domain` | `bacteria` | 23 | 0 | 0 | 23 | 0 |
| `C2_gaps_to_gaps_learned` | `gap_kind` | `already_gap` | 7 | 0 | 0 | 7 | 0 |
| `C2_gaps_to_gaps_learned` | `gap_kind` | `maskable_gap` | 16 | 0 | 0 | 16 | 0 |
| `C2_gaps_to_gaps_learned` | `pathway_robust` | `marginal` | 7 | 0 | 0 | 7 | 0 |
| `C2_gaps_to_gaps_learned` | `pathway_robust` | `robust` | 16 | 0 | 0 | 16 | 0 |
| `C3_universal_vs_gaps_learned` | `domain` | `bacteria` | 23 | 0 | 0 | 23 | 0 |
| `C3_universal_vs_gaps_learned` | `gap_kind` | `already_gap` | 7 | 0 | 0 | 7 | 0 |
| `C3_universal_vs_gaps_learned` | `gap_kind` | `maskable_gap` | 16 | 0 | 0 | 16 | 0 |
| `C3_universal_vs_gaps_learned` | `pathway_robust` | `marginal` | 7 | 0 | 0 | 7 | 0 |
| `C3_universal_vs_gaps_learned` | `pathway_robust` | `robust` | 16 | 0 | 0 | 16 | 0 |
| `C3_pathway_vs_gaps_learned` | `domain` | `bacteria` | 23 | 0 | 0 | 23 | 0 |
| `C3_pathway_vs_gaps_learned` | `gap_kind` | `already_gap` | 7 | 0 | 0 | 7 | 0 |
| `C3_pathway_vs_gaps_learned` | `gap_kind` | `maskable_gap` | 16 | 0 | 0 | 16 | 0 |
| `C3_pathway_vs_gaps_learned` | `pathway_robust` | `marginal` | 7 | 0 | 0 | 7 | 0 |
| `C3_pathway_vs_gaps_learned` | `pathway_robust` | `robust` | 16 | 0 | 0 | 16 | 0 |
| `C3_organism_vs_gaps_learned` | `domain` | `bacteria` | 23 | 0 | 0 | 23 | 0 |
| `C3_organism_vs_gaps_learned` | `gap_kind` | `already_gap` | 7 | 0 | 0 | 7 | 0 |
| `C3_organism_vs_gaps_learned` | `gap_kind` | `maskable_gap` | 16 | 0 | 0 | 16 | 0 |
| `C3_organism_vs_gaps_learned` | `pathway_robust` | `marginal` | 7 | 0 | 0 | 7 | 0 |
| `C3_organism_vs_gaps_learned` | `pathway_robust` | `robust` | 16 | 0 | 0 | 16 | 0 |
| `C4_pathway_vs_universal` | `domain` | `bacteria` | 23 | 0 | 0 | 23 | 0 |
| `C4_pathway_vs_universal` | `gap_kind` | `already_gap` | 7 | 0 | 0 | 7 | 0 |
| `C4_pathway_vs_universal` | `gap_kind` | `maskable_gap` | 16 | 0 | 0 | 16 | 0 |
| `C4_pathway_vs_universal` | `pathway_robust` | `marginal` | 7 | 0 | 0 | 7 | 0 |
| `C4_pathway_vs_universal` | `pathway_robust` | `robust` | 16 | 0 | 0 | 16 | 0 |
| `C4_organism_vs_universal` | `domain` | `bacteria` | 23 | 0 | 0 | 23 | 0 |
| `C4_organism_vs_universal` | `gap_kind` | `already_gap` | 7 | 0 | 0 | 7 | 0 |
| `C4_organism_vs_universal` | `gap_kind` | `maskable_gap` | 16 | 0 | 0 | 16 | 0 |
| `C4_organism_vs_universal` | `pathway_robust` | `marginal` | 7 | 0 | 0 | 7 | 0 |
| `C4_organism_vs_universal` | `pathway_robust` | `robust` | 16 | 0 | 0 | 16 | 0 |
| `C4_organism_vs_pathway` | `domain` | `bacteria` | 23 | 0 | 0 | 23 | 0 |
| `C4_organism_vs_pathway` | `gap_kind` | `already_gap` | 7 | 0 | 0 | 7 | 0 |
| `C4_organism_vs_pathway` | `gap_kind` | `maskable_gap` | 16 | 0 | 0 | 16 | 0 |
| `C4_organism_vs_pathway` | `pathway_robust` | `marginal` | 7 | 0 | 0 | 7 | 0 |
| `C4_organism_vs_pathway` | `pathway_robust` | `robust` | 16 | 0 | 0 | 16 | 0 |
| `C5_universal_w5000` | `domain` | `bacteria` | 23 | 3 | 0 | 20 | 3 |
| `C5_universal_w5000` | `gap_kind` | `already_gap` | 7 | 1 | 0 | 6 | 1 |
| `C5_universal_w5000` | `gap_kind` | `maskable_gap` | 16 | 2 | 0 | 14 | 2 |
| `C5_universal_w5000` | `pathway_robust` | `marginal` | 7 | 1 | 0 | 6 | 1 |
| `C5_universal_w5000` | `pathway_robust` | `robust` | 16 | 2 | 0 | 14 | 2 |
| `C5_pathway_w5000` | `domain` | `bacteria` | 23 | 3 | 0 | 20 | 3 |
| `C5_pathway_w5000` | `gap_kind` | `already_gap` | 7 | 1 | 0 | 6 | 1 |
| `C5_pathway_w5000` | `gap_kind` | `maskable_gap` | 16 | 2 | 0 | 14 | 2 |
| `C5_pathway_w5000` | `pathway_robust` | `marginal` | 7 | 1 | 0 | 6 | 1 |
| `C5_pathway_w5000` | `pathway_robust` | `robust` | 16 | 2 | 0 | 14 | 2 |
| `C5_organism_w5000` | `domain` | `bacteria` | 23 | 3 | 0 | 20 | 3 |
| `C5_organism_w5000` | `gap_kind` | `already_gap` | 7 | 1 | 0 | 6 | 1 |
| `C5_organism_w5000` | `gap_kind` | `maskable_gap` | 16 | 2 | 0 | 14 | 2 |
| `C5_organism_w5000` | `pathway_robust` | `marginal` | 7 | 1 | 0 | 6 | 1 |
| `C5_organism_w5000` | `pathway_robust` | `robust` | 16 | 2 | 0 | 14 | 2 |
| `C5total_universal_w5000` | `domain` | `bacteria` | 23 | 3 | 0 | 20 | 3 |
| `C5total_universal_w5000` | `gap_kind` | `already_gap` | 7 | 1 | 0 | 6 | 1 |
| `C5total_universal_w5000` | `gap_kind` | `maskable_gap` | 16 | 2 | 0 | 14 | 2 |
| `C5total_universal_w5000` | `pathway_robust` | `marginal` | 7 | 1 | 0 | 6 | 1 |
| `C5total_universal_w5000` | `pathway_robust` | `robust` | 16 | 2 | 0 | 14 | 2 |
| `C5_pathway_vs_universal_w5000` | `domain` | `bacteria` | 23 | 0 | 0 | 23 | 0 |
| `C5_pathway_vs_universal_w5000` | `gap_kind` | `already_gap` | 7 | 0 | 0 | 7 | 0 |
| `C5_pathway_vs_universal_w5000` | `gap_kind` | `maskable_gap` | 16 | 0 | 0 | 16 | 0 |
| `C5_pathway_vs_universal_w5000` | `pathway_robust` | `marginal` | 7 | 0 | 0 | 7 | 0 |
| `C5_pathway_vs_universal_w5000` | `pathway_robust` | `robust` | 16 | 0 | 0 | 16 | 0 |
| `C5_organism_vs_universal_w5000` | `domain` | `bacteria` | 23 | 0 | 0 | 23 | 0 |
| `C5_organism_vs_universal_w5000` | `gap_kind` | `already_gap` | 7 | 0 | 0 | 7 | 0 |
| `C5_organism_vs_universal_w5000` | `gap_kind` | `maskable_gap` | 16 | 0 | 0 | 16 | 0 |
| `C5_organism_vs_universal_w5000` | `pathway_robust` | `marginal` | 7 | 0 | 0 | 7 | 0 |
| `C5_organism_vs_universal_w5000` | `pathway_robust` | `robust` | 16 | 0 | 0 | 16 | 0 |
| `C5_organism_vs_pathway_w5000` | `domain` | `bacteria` | 23 | 0 | 0 | 23 | 0 |
| `C5_organism_vs_pathway_w5000` | `gap_kind` | `already_gap` | 7 | 0 | 0 | 7 | 0 |
| `C5_organism_vs_pathway_w5000` | `gap_kind` | `maskable_gap` | 16 | 0 | 0 | 16 | 0 |
| `C5_organism_vs_pathway_w5000` | `pathway_robust` | `marginal` | 7 | 0 | 0 | 7 | 0 |
| `C5_organism_vs_pathway_w5000` | `pathway_robust` | `robust` | 16 | 0 | 0 | 16 | 0 |
| `C5_universal_w10000` | `domain` | `bacteria` | 23 | 3 | 0 | 20 | 3 |
| `C5_universal_w10000` | `gap_kind` | `already_gap` | 7 | 1 | 0 | 6 | 1 |
| `C5_universal_w10000` | `gap_kind` | `maskable_gap` | 16 | 2 | 0 | 14 | 2 |
| `C5_universal_w10000` | `pathway_robust` | `marginal` | 7 | 1 | 0 | 6 | 1 |
| `C5_universal_w10000` | `pathway_robust` | `robust` | 16 | 2 | 0 | 14 | 2 |
| `C5_pathway_w10000` | `domain` | `bacteria` | 23 | 3 | 0 | 20 | 3 |
| `C5_pathway_w10000` | `gap_kind` | `already_gap` | 7 | 1 | 0 | 6 | 1 |
| `C5_pathway_w10000` | `gap_kind` | `maskable_gap` | 16 | 2 | 0 | 14 | 2 |
| `C5_pathway_w10000` | `pathway_robust` | `marginal` | 7 | 1 | 0 | 6 | 1 |
| `C5_pathway_w10000` | `pathway_robust` | `robust` | 16 | 2 | 0 | 14 | 2 |
| `C5_organism_w10000` | `domain` | `bacteria` | 23 | 3 | 0 | 20 | 3 |
| `C5_organism_w10000` | `gap_kind` | `already_gap` | 7 | 1 | 0 | 6 | 1 |
| `C5_organism_w10000` | `gap_kind` | `maskable_gap` | 16 | 2 | 0 | 14 | 2 |
| `C5_organism_w10000` | `pathway_robust` | `marginal` | 7 | 1 | 0 | 6 | 1 |
| `C5_organism_w10000` | `pathway_robust` | `robust` | 16 | 2 | 0 | 14 | 2 |
| `C5total_universal_w10000` | `domain` | `bacteria` | 23 | 3 | 0 | 20 | 3 |
| `C5total_universal_w10000` | `gap_kind` | `already_gap` | 7 | 1 | 0 | 6 | 1 |
| `C5total_universal_w10000` | `gap_kind` | `maskable_gap` | 16 | 2 | 0 | 14 | 2 |
| `C5total_universal_w10000` | `pathway_robust` | `marginal` | 7 | 1 | 0 | 6 | 1 |
| `C5total_universal_w10000` | `pathway_robust` | `robust` | 16 | 2 | 0 | 14 | 2 |
| `C5_pathway_vs_universal_w10000` | `domain` | `bacteria` | 23 | 0 | 0 | 23 | 0 |
| `C5_pathway_vs_universal_w10000` | `gap_kind` | `already_gap` | 7 | 0 | 0 | 7 | 0 |
| `C5_pathway_vs_universal_w10000` | `gap_kind` | `maskable_gap` | 16 | 0 | 0 | 16 | 0 |
| `C5_pathway_vs_universal_w10000` | `pathway_robust` | `marginal` | 7 | 0 | 0 | 7 | 0 |
| `C5_pathway_vs_universal_w10000` | `pathway_robust` | `robust` | 16 | 0 | 0 | 16 | 0 |
| `C5_organism_vs_universal_w10000` | `domain` | `bacteria` | 23 | 0 | 0 | 23 | 0 |
| `C5_organism_vs_universal_w10000` | `gap_kind` | `already_gap` | 7 | 0 | 0 | 7 | 0 |
| `C5_organism_vs_universal_w10000` | `gap_kind` | `maskable_gap` | 16 | 0 | 0 | 16 | 0 |
| `C5_organism_vs_universal_w10000` | `pathway_robust` | `marginal` | 7 | 0 | 0 | 7 | 0 |
| `C5_organism_vs_universal_w10000` | `pathway_robust` | `robust` | 16 | 0 | 0 | 16 | 0 |
| `C5_organism_vs_pathway_w10000` | `domain` | `bacteria` | 23 | 0 | 0 | 23 | 0 |
| `C5_organism_vs_pathway_w10000` | `gap_kind` | `already_gap` | 7 | 0 | 0 | 7 | 0 |
| `C5_organism_vs_pathway_w10000` | `gap_kind` | `maskable_gap` | 16 | 0 | 0 | 16 | 0 |
| `C5_organism_vs_pathway_w10000` | `pathway_robust` | `marginal` | 7 | 0 | 0 | 7 | 0 |
| `C5_organism_vs_pathway_w10000` | `pathway_robust` | `robust` | 16 | 0 | 0 | 16 | 0 |
| `C5_universal_w20000` | `domain` | `bacteria` | 23 | 3 | 0 | 20 | 3 |
| `C5_universal_w20000` | `gap_kind` | `already_gap` | 7 | 1 | 0 | 6 | 1 |
| `C5_universal_w20000` | `gap_kind` | `maskable_gap` | 16 | 2 | 0 | 14 | 2 |
| `C5_universal_w20000` | `pathway_robust` | `marginal` | 7 | 1 | 0 | 6 | 1 |
| `C5_universal_w20000` | `pathway_robust` | `robust` | 16 | 2 | 0 | 14 | 2 |
| `C5_pathway_w20000` | `domain` | `bacteria` | 23 | 3 | 0 | 20 | 3 |
| `C5_pathway_w20000` | `gap_kind` | `already_gap` | 7 | 1 | 0 | 6 | 1 |
| `C5_pathway_w20000` | `gap_kind` | `maskable_gap` | 16 | 2 | 0 | 14 | 2 |
| `C5_pathway_w20000` | `pathway_robust` | `marginal` | 7 | 1 | 0 | 6 | 1 |
| `C5_pathway_w20000` | `pathway_robust` | `robust` | 16 | 2 | 0 | 14 | 2 |
| `C5_organism_w20000` | `domain` | `bacteria` | 23 | 3 | 0 | 20 | 3 |
| `C5_organism_w20000` | `gap_kind` | `already_gap` | 7 | 1 | 0 | 6 | 1 |
| `C5_organism_w20000` | `gap_kind` | `maskable_gap` | 16 | 2 | 0 | 14 | 2 |
| `C5_organism_w20000` | `pathway_robust` | `marginal` | 7 | 1 | 0 | 6 | 1 |
| `C5_organism_w20000` | `pathway_robust` | `robust` | 16 | 2 | 0 | 14 | 2 |
| `C5total_universal_w20000` | `domain` | `bacteria` | 23 | 3 | 0 | 20 | 3 |
| `C5total_universal_w20000` | `gap_kind` | `already_gap` | 7 | 1 | 0 | 6 | 1 |
| `C5total_universal_w20000` | `gap_kind` | `maskable_gap` | 16 | 2 | 0 | 14 | 2 |
| `C5total_universal_w20000` | `pathway_robust` | `marginal` | 7 | 1 | 0 | 6 | 1 |
| `C5total_universal_w20000` | `pathway_robust` | `robust` | 16 | 2 | 0 | 14 | 2 |
| `C5_pathway_vs_universal_w20000` | `domain` | `bacteria` | 23 | 0 | 0 | 23 | 0 |
| `C5_pathway_vs_universal_w20000` | `gap_kind` | `already_gap` | 7 | 0 | 0 | 7 | 0 |
| `C5_pathway_vs_universal_w20000` | `gap_kind` | `maskable_gap` | 16 | 0 | 0 | 16 | 0 |
| `C5_pathway_vs_universal_w20000` | `pathway_robust` | `marginal` | 7 | 0 | 0 | 7 | 0 |
| `C5_pathway_vs_universal_w20000` | `pathway_robust` | `robust` | 16 | 0 | 0 | 16 | 0 |
| `C5_organism_vs_universal_w20000` | `domain` | `bacteria` | 23 | 0 | 0 | 23 | 0 |
| `C5_organism_vs_universal_w20000` | `gap_kind` | `already_gap` | 7 | 0 | 0 | 7 | 0 |
| `C5_organism_vs_universal_w20000` | `gap_kind` | `maskable_gap` | 16 | 0 | 0 | 16 | 0 |
| `C5_organism_vs_universal_w20000` | `pathway_robust` | `marginal` | 7 | 0 | 0 | 7 | 0 |
| `C5_organism_vs_universal_w20000` | `pathway_robust` | `robust` | 16 | 0 | 0 | 16 | 0 |
| `C5_organism_vs_pathway_w20000` | `domain` | `bacteria` | 23 | 0 | 0 | 23 | 0 |
| `C5_organism_vs_pathway_w20000` | `gap_kind` | `already_gap` | 7 | 0 | 0 | 7 | 0 |
| `C5_organism_vs_pathway_w20000` | `gap_kind` | `maskable_gap` | 16 | 0 | 0 | 16 | 0 |
| `C5_organism_vs_pathway_w20000` | `pathway_robust` | `marginal` | 7 | 0 | 0 | 7 | 0 |
| `C5_organism_vs_pathway_w20000` | `pathway_robust` | `robust` | 16 | 0 | 0 | 16 | 0 |
| `C5_universal` | `domain` | `bacteria` | 23 | 3 | 0 | 20 | 3 |
| `C5_universal` | `gap_kind` | `already_gap` | 7 | 1 | 0 | 6 | 1 |
| `C5_universal` | `gap_kind` | `maskable_gap` | 16 | 2 | 0 | 14 | 2 |
| `C5_universal` | `pathway_robust` | `marginal` | 7 | 1 | 0 | 6 | 1 |
| `C5_universal` | `pathway_robust` | `robust` | 16 | 2 | 0 | 14 | 2 |
| `C5_pathway` | `domain` | `bacteria` | 23 | 3 | 0 | 20 | 3 |
| `C5_pathway` | `gap_kind` | `already_gap` | 7 | 1 | 0 | 6 | 1 |
| `C5_pathway` | `gap_kind` | `maskable_gap` | 16 | 2 | 0 | 14 | 2 |
| `C5_pathway` | `pathway_robust` | `marginal` | 7 | 1 | 0 | 6 | 1 |
| `C5_pathway` | `pathway_robust` | `robust` | 16 | 2 | 0 | 14 | 2 |
| `C5_organism` | `domain` | `bacteria` | 23 | 3 | 0 | 20 | 3 |
| `C5_organism` | `gap_kind` | `already_gap` | 7 | 1 | 0 | 6 | 1 |
| `C5_organism` | `gap_kind` | `maskable_gap` | 16 | 2 | 0 | 14 | 2 |
| `C5_organism` | `pathway_robust` | `marginal` | 7 | 1 | 0 | 6 | 1 |
| `C5_organism` | `pathway_robust` | `robust` | 16 | 2 | 0 | 14 | 2 |
| `C5_pathway_vs_universal` | `domain` | `bacteria` | 23 | 0 | 0 | 23 | 0 |
| `C5_pathway_vs_universal` | `gap_kind` | `already_gap` | 7 | 0 | 0 | 7 | 0 |
| `C5_pathway_vs_universal` | `gap_kind` | `maskable_gap` | 16 | 0 | 0 | 16 | 0 |
| `C5_pathway_vs_universal` | `pathway_robust` | `marginal` | 7 | 0 | 0 | 7 | 0 |
| `C5_pathway_vs_universal` | `pathway_robust` | `robust` | 16 | 0 | 0 | 16 | 0 |
| `C5_organism_vs_universal` | `domain` | `bacteria` | 23 | 0 | 0 | 23 | 0 |
| `C5_organism_vs_universal` | `gap_kind` | `already_gap` | 7 | 0 | 0 | 7 | 0 |
| `C5_organism_vs_universal` | `gap_kind` | `maskable_gap` | 16 | 0 | 0 | 16 | 0 |
| `C5_organism_vs_universal` | `pathway_robust` | `marginal` | 7 | 0 | 0 | 7 | 0 |
| `C5_organism_vs_universal` | `pathway_robust` | `robust` | 16 | 0 | 0 | 16 | 0 |
| `C5_organism_vs_pathway` | `domain` | `bacteria` | 23 | 0 | 0 | 23 | 0 |
| `C5_organism_vs_pathway` | `gap_kind` | `already_gap` | 7 | 0 | 0 | 7 | 0 |
| `C5_organism_vs_pathway` | `gap_kind` | `maskable_gap` | 16 | 0 | 0 | 16 | 0 |
| `C5_organism_vs_pathway` | `pathway_robust` | `marginal` | 7 | 0 | 0 | 7 | 0 |
| `C5_organism_vs_pathway` | `pathway_robust` | `robust` | 16 | 0 | 0 | 16 | 0 |
| `S1_stage0_to_sham_gaps` | `domain` | `bacteria` | 23 | 0 | 0 | 23 | 0 |
| `S1_stage0_to_sham_gaps` | `gap_kind` | `already_gap` | 7 | 0 | 0 | 7 | 0 |
| `S1_stage0_to_sham_gaps` | `gap_kind` | `maskable_gap` | 16 | 0 | 0 | 16 | 0 |
| `S1_stage0_to_sham_gaps` | `pathway_robust` | `marginal` | 7 | 0 | 0 | 7 | 0 |
| `S1_stage0_to_sham_gaps` | `pathway_robust` | `robust` | 16 | 0 | 0 | 16 | 0 |
| `S1_sham_gaps_to_sham_ctx_universal_w5000` | `domain` | `bacteria` | 23 | 1 | 0 | 22 | 1 |
| `S1_sham_gaps_to_sham_ctx_universal_w5000` | `gap_kind` | `already_gap` | 7 | 1 | 0 | 6 | 1 |
| `S1_sham_gaps_to_sham_ctx_universal_w5000` | `gap_kind` | `maskable_gap` | 16 | 0 | 0 | 16 | 0 |
| `S1_sham_gaps_to_sham_ctx_universal_w5000` | `pathway_robust` | `marginal` | 7 | 0 | 0 | 7 | 0 |
| `S1_sham_gaps_to_sham_ctx_universal_w5000` | `pathway_robust` | `robust` | 16 | 1 | 0 | 15 | 1 |
| `S1_stage0_to_sham_ctx_universal_w5000` | `domain` | `bacteria` | 23 | 1 | 0 | 22 | 1 |
| `S1_stage0_to_sham_ctx_universal_w5000` | `gap_kind` | `already_gap` | 7 | 1 | 0 | 6 | 1 |
| `S1_stage0_to_sham_ctx_universal_w5000` | `gap_kind` | `maskable_gap` | 16 | 0 | 0 | 16 | 0 |
| `S1_stage0_to_sham_ctx_universal_w5000` | `pathway_robust` | `marginal` | 7 | 0 | 0 | 7 | 0 |
| `S1_stage0_to_sham_ctx_universal_w5000` | `pathway_robust` | `robust` | 16 | 1 | 0 | 15 | 1 |
| `S1_sham_gaps_to_sham_ctx_universal_w10000` | `domain` | `bacteria` | 23 | 1 | 0 | 22 | 1 |
| `S1_sham_gaps_to_sham_ctx_universal_w10000` | `gap_kind` | `already_gap` | 7 | 1 | 0 | 6 | 1 |
| `S1_sham_gaps_to_sham_ctx_universal_w10000` | `gap_kind` | `maskable_gap` | 16 | 0 | 0 | 16 | 0 |
| `S1_sham_gaps_to_sham_ctx_universal_w10000` | `pathway_robust` | `marginal` | 7 | 0 | 0 | 7 | 0 |
| `S1_sham_gaps_to_sham_ctx_universal_w10000` | `pathway_robust` | `robust` | 16 | 1 | 0 | 15 | 1 |
| `S1_stage0_to_sham_ctx_universal_w10000` | `domain` | `bacteria` | 23 | 1 | 0 | 22 | 1 |
| `S1_stage0_to_sham_ctx_universal_w10000` | `gap_kind` | `already_gap` | 7 | 1 | 0 | 6 | 1 |
| `S1_stage0_to_sham_ctx_universal_w10000` | `gap_kind` | `maskable_gap` | 16 | 0 | 0 | 16 | 0 |
| `S1_stage0_to_sham_ctx_universal_w10000` | `pathway_robust` | `marginal` | 7 | 0 | 0 | 7 | 0 |
| `S1_stage0_to_sham_ctx_universal_w10000` | `pathway_robust` | `robust` | 16 | 1 | 0 | 15 | 1 |
| `S1_sham_gaps_to_sham_ctx_universal_w20000` | `domain` | `bacteria` | 23 | 1 | 0 | 22 | 1 |
| `S1_sham_gaps_to_sham_ctx_universal_w20000` | `gap_kind` | `already_gap` | 7 | 1 | 0 | 6 | 1 |
| `S1_sham_gaps_to_sham_ctx_universal_w20000` | `gap_kind` | `maskable_gap` | 16 | 0 | 0 | 16 | 0 |
| `S1_sham_gaps_to_sham_ctx_universal_w20000` | `pathway_robust` | `marginal` | 7 | 0 | 0 | 7 | 0 |
| `S1_sham_gaps_to_sham_ctx_universal_w20000` | `pathway_robust` | `robust` | 16 | 1 | 0 | 15 | 1 |
| `S1_stage0_to_sham_ctx_universal_w20000` | `domain` | `bacteria` | 23 | 1 | 0 | 22 | 1 |
| `S1_stage0_to_sham_ctx_universal_w20000` | `gap_kind` | `already_gap` | 7 | 1 | 0 | 6 | 1 |
| `S1_stage0_to_sham_ctx_universal_w20000` | `gap_kind` | `maskable_gap` | 16 | 0 | 0 | 16 | 0 |
| `S1_stage0_to_sham_ctx_universal_w20000` | `pathway_robust` | `marginal` | 7 | 0 | 0 | 7 | 0 |
| `S1_stage0_to_sham_ctx_universal_w20000` | `pathway_robust` | `robust` | 16 | 1 | 0 | 15 | 1 |
| `S1_sham_gaps_to_sham_ctx_universal` | `domain` | `bacteria` | 23 | 1 | 0 | 22 | 1 |
| `S1_sham_gaps_to_sham_ctx_universal` | `gap_kind` | `already_gap` | 7 | 1 | 0 | 6 | 1 |
| `S1_sham_gaps_to_sham_ctx_universal` | `gap_kind` | `maskable_gap` | 16 | 0 | 0 | 16 | 0 |
| `S1_sham_gaps_to_sham_ctx_universal` | `pathway_robust` | `marginal` | 7 | 0 | 0 | 7 | 0 |
| `S1_sham_gaps_to_sham_ctx_universal` | `pathway_robust` | `robust` | 16 | 1 | 0 | 15 | 1 |
| `U1_unmasked_gaps_vs_gaps` | `domain` | `bacteria` | 23 | 0 | 0 | 23 | 0 |
| `U1_unmasked_gaps_vs_gaps` | `gap_kind` | `already_gap` | 7 | 0 | 0 | 7 | 0 |
| `U1_unmasked_gaps_vs_gaps` | `gap_kind` | `maskable_gap` | 16 | 0 | 0 | 16 | 0 |
| `U1_unmasked_gaps_vs_gaps` | `pathway_robust` | `marginal` | 7 | 0 | 0 | 7 | 0 |
| `U1_unmasked_gaps_vs_gaps` | `pathway_robust` | `robust` | 16 | 0 | 0 | 16 | 0 |

## 4. Learned-Stage Candidate Counts

| Arm | Candidates | Rescored | Unsupported | Unchanged unsupported | Retracted | Scored only |
|---|---:|---:|---:|---:|---:|---:|
| `gaps_learned` | 36657 | 0 | 36657 | 36657 | 0 | 0 |
| `organism` | 36657 | 0 | 36657 | 36657 | 0 | 0 |
| `pathway` | 36657 | 0 | 36657 | 36657 | 0 | 0 |
| `universal` | 36657 | 0 | 36657 | 36657 | 0 | 0 |

## 5. Masking & Build Audit Summary

| Check | Count |
|---|---:|
| `total_cells` | 23 |
| `complete_cells` | 23 |
| `candidate_universe_preserved_cells` | 23 |
| `masked_gap_present_cells` | 23 |
| `gaps_ctxjar_identical_cells` | 23 |
| `already_gap_unmasked_total` | 7 |
| `already_gap_unmasked_matches` | 7 |
