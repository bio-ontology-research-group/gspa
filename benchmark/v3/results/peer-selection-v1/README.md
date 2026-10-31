# Close-relative peer genome selection v1 (2026-10-28)

Produced by `benchmark/v3/select_peer_genomes.py` (NCBI Datasets v2 family reports, complete RefSeq
assemblies, atypical excluded; species resolved through NCBI Taxonomy EFETCH batches; raw responses
cached under `~/tmp/gspa-peer-selection-v1/out/ncbi-cache`, 280 family pages + 17 taxonomy batches).
For each of the 18 label-rich cohort targets: up to 6 close peers (other species of the same family,
same genus first, one assembly per species) and one near-duplicate (another assembly of the target's
species). Both members of every cohort GenBank/RefSeq pair are excluded, so a target's own RefSeq twin
never counts as a peer (Aeropyrum pernix and Thermoplasma acidophilum therefore have no
near-duplicate; Paramagnetospirillum has no other complete genome in its family). Taxonomy drift is
recorded (Aeropyrum: Desulfurococcaceae → Acidilobaceae in current NCBI). 111 peer rows. No genome
package is downloaded by this step; `peers.tsv` is the input to the recurrence experiment (WP4).
