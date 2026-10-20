# Empty Quarter genome annotation with the full GSPA pipeline

Annotates every high-quality (completeness >= 90 %, contamination <= 5 %) and medium-quality
(>= 50 %, <= 10 %) genome of the Empty Quarter project, from two sources:

- the soil bins of the single-sample SemiBin run: 16,333 bins (4,030 HQ, 12,303 MQ) from 542 samples, as
  listed by `rob/metagenomics/results/checkm2_single_sample_full.tsv` on ibex (dataset `soil_bin`);
- the culture genomes of the RAK-PR catalogue (Alejandra's Rub' al-Khali pangenome reference,
  `alelopezv/eqmgc/catalogue.tsv`): 722 HQ/MQ bacterial genomes, all of them from cultures: 520 bins of
  enrichment cultures (Flye `--meta` + MetaBAT2), 131 isolate assemblies (hifiasm) and 71 genomes of the
  thermophilic cultures, which the catalogue files under batch `soil` (datasets `culture_enrichment`,
  `culture_isolate`, `culture_thermophilic`). Assembly paths come from the catalogue's input registry
  (`inputs.tsv`), Bakta annotations from `results/annotation/bakta` (646 of the 722; the other 76 get
  pyrodigal). Genome ids and sample names of the two sources are disjoint.

Every genome gets gene calling, the three Stage 0 channels (DIAMOND, Pfam, InterProScan), a gapsmith
metabolic model and GSPA Stages 0, 1, 2A and 2B. Stage 3 (the LLM agent) is not part of this run. The
culture genomes join the soil bins in the same protein clustering, so a culture protein shares its
channel results with its near-identical soil relatives and the reverse.

Host: KAUST ibex. Work area `/encrypted/e3008/gspa/eq_mags` (`code/`, `ref/`, `lists/`, `work/<mag>/`,
`cluster/`, `chunks/`, `claims_by_sample/`, `logs/`, `results/`). Every step is a Slurm job; `run.sh`
only submits and copies.

## Design: channels once per protein cluster

The 542 samples are site replicates, so the bins carry near-identical genomes many times over. The
channels therefore run once per cluster of proteins rather than once per bin: all proteomes are pooled
(ids `<mag>|<protein_id>`), clustered with mmseqs linclust at 90 % identity and 90 % bidirectional
coverage (`--min-seq-id 0.9 -c 0.9 --cov-mode 0`), and DIAMOND, Pfam and InterProScan run on the
cluster representatives in chunks of 25,000. The resulting claims (InterPro entries and their
InterPro2GO terms, Pfam domains with pfam2go, DIAMOND hits with experimental GO) are then copied to every
member of the cluster. Domain architecture is stable at 90 % identity with 90 % coverage in both
directions, and GSPA consumes the entries and GO terms, not match coordinates, so a member receives the
same claims a direct run would give it; the coverage condition keeps fragments from inheriting the
domains of a full-length relative or the reverse. The linclust of the project's protein catalog
(`04_structure`) suggests about 2.3 times fewer representatives than proteins.

| Step | Script | What it does | Resources |
| --- | --- | --- | --- |
| select | `00_select.py` | HQ/MQ filter of the soil bins, every RAK-PR catalogue genome, GTDB lineage -> NCBI taxon of GSPA's taxon-constraint table (most specific rank whose name the table knows; GTDB synonyms translated, placeholder names skipped), biomass template (gram: monoderm phyla and archaea -> pos), essential profile, input paths | debug, 1 CPU |
| 01 proteins | `01_proteins.sbatch` (array per genome) | `annot/<mag>.{faa,gff3,fna}` from the Bakta annotation when present (15,571 soil bins from the metagenomics tree, 646 culture genomes from RAK-PR) or pyrodigal 3.7.1 otherwise | 2 CPU, 8 GB |
| 02 cluster | `02_cluster.sbatch` (one job, after 01) | pool all proteomes, linclust, per-bin member maps, representative chunks; then submits 03, 04 and 06 | 32 CPU, 180 GB |
| 03 chunk channels | `03_chunk_channels.sbatch` (array per chunk) | DIAMOND 2.1.16 vs experimental-GO Swiss-Prot (`--evalue 1e-5 --max-target-seqs 5`); HMMER 3.3.2 Pfam-A `--cut_ga`; InterProScan 5.61 (Pfam, TIGRFAM, Hamap, PIRSF, SMART, CDD, ProSite, PRINTS, SFLD; `-goterms -iprlookup`); claims via `bgc_benchmark.py build-claims`, keyed by representative | 8 CPU, 24 GB, <= 24 h |
| 04 scatter | `04_scatter.py` (one job, afterok 03) | claims of every representative copied to its members; one JSONL per sample; `coverage.tsv` per bin | 4 CPU, 64 GB |
| 05 model | `05_gapsmith.sbatch` (array per bin, aftercorr 01) | gapsmith 0.1.0 `doall` (find, find-transport, draft, medium, fill; DIAMOND aligner; `-t Archaea` for archaeal bins) on the bin's proteome, uniformly for all bins; a gapseq model of the HQ run is linked under `gapseq_ref/` for comparison only | 4 CPU, 12 GB |
| 06 integrate | `06_integrate.sbatch` (array per bin, afterok 04 and aftercorr 05) | the bin's claims from its sample file; `eq_prep.py` (layout and operon chains from the GFF3, gap reactions from the model via P14 `prep`); Stage 0 (`--enable-priors none`); Stage 1 + 2A + 2B (`--enable-priors` all six, `--taxon`, `--essential-profile`, operons, gaps, RLGC suggester, `--enable-fba-refiner` with the ModelSEED universal model, `--enforce-stage2 --enforce-consistency-mode relax`); P14 `score` | 4 CPU, 20 GB |
| aggregate | `04_aggregate.py` | one row per genome and a summary by class, by dataset and overall | debug |

Outputs per bin under `work/<mag>/`: `annot/`, `gapseq/` (gapsmith model as `<mag>.xml`, `<mag>-draft.xml`,
tables), `gapseq_ref/` (gapseq model when one exists), `<mag>_claims.jsonl.gz`, `<mag>_{layout,operons}.tsv`,
`<mag>_gaps.jsonl`, `<mag>_stage0.tsv.gz`, `<mag>_stage12.tsv.gz`, `<mag>_stage2_actions.tsv`,
`<mag>_rlgc_suggestions.tsv`, `<mag>_score.tsv`, `timing*.tsv`, `DONE_proteins`, `DONE_gapseq`, `DONE`.
Channel results per chunk under `chunks/chunk_NNNNN/` (`diamond.tsv.gz`, `pfam.domtbl.gz`, `ips.tsv.gz`,
`claims.jsonl`); cluster under `cluster/` (`clu_rep_seq.fasta`, `clu_cluster.tsv`, `members_by_bin/`).

## Why gapsmith for every bin

The HQ run's gapseq models (3,794 bins) were built from contigs with tblastn; gapsmith works from
proteins. On three bacterial HQ bins the two agree only partly (draft reaction Jaccard 0.63-0.66,
gap-filled 0.64-0.68, about half as many predicted pathways from gapsmith; `test/compare_hq_bacteria.tsv`).
Mixing the two sources would make the gap prior depend on which tool a bin happened to get, so every bin
gets a gapsmith model under one setting, at about a minute per bin; the gapseq models stay available
through `run.sh compare`. The gapseq fallback script is in `legacy/`.

## Scheduling

The per-user limits on ibex are 1,300 CPUs and 2,000 submitted jobs, and every array task counts as a
job, so the per-genome steps (01, 05, 06) run on blocks of 50 consecutive rows of `bins.tsv` per array
task (`GENOMES_PER_TASK`; 342 tasks for 17,055 genomes) and the channel chunks hold 50,000
representatives (about 6 h on 8 CPUs; a few hundred tasks). `run.sh launch <JOBID>` submits the
proteins array with `--dependency=afterany:<JOBID>` (the AllTheBacteria Stage D array), the gapsmith
array with `aftercorr` on it (same block index), and the cluster job with `afterok` on it; the cluster
job submits the chunk array (throttle 140 x 8 CPUs), the scatter job and the integration array
(`afterok` scatter, `aftercorr` gapsmith). `run.sh rebalance THR ARRAY` changes a throttle. Expected
cost: about 22,000 CPU-h for the channels on representatives (InterProScan dominates; the 722 culture
genomes add about 2.5 million proteins, most of them clustering with soil proteins), about 1,500 for
gapsmith and 3,000 for integration; about a day of wall time at the quota once the cluster job has run.

## Validation

`run.sh test` runs the whole chain on five genomes (an HQ archaeal soil bin with gapseq model and Bakta
annotation, an MQ soil bin with Bakta annotation, an HQ soil bin without Bakta annotation, an isolate with
Bakta annotation, a thermophilic culture genome without) under `test/` on ibex. The run of 2026-10-20
(jobs 54067424-6, 54071845-7, 54078956) completed every step for all five; the numbers are in
`test/validation_2026-10-20.tsv` (fetched copy: `results/test_validation_2026-10-20.tsv`).

Clustered against direct channels: for the three soil bins that also have a direct per-bin run
(`legacy/`), the GO claims of the clustered run are the same set to within a few DIAMOND claims
(Jaccard 0.9985, 1.0 and 0.9987 over all claims; Pfam and InterProScan claims identical for two bins and
3 of 2,808 and 4,143 apart for the third). The member of a cluster receives the DIAMOND hits of its
representative, which differ from its own hits when the representative's best five subjects differ; the
Pfam and InterProScan differences are single domains at the 90 % identity boundary.

EC claims and Stage 2B: `add_ec_claims.py` adds about 10 % EC claims on top of the GO claims (ec2go on GO
terms with one complete EC, Swiss-Prot EC transfer from the DIAMOND hits); with them the
StoichiometricLocusRefiner finds a core of 494 to 2,933 reactions per genome and repairs 1 to 25.

Gap prior and RLGC with gapsmith models: gapsmith writes no EC annotations into its SBML, which left the
P14 prep without a single fillable gap (`n_gap_fillable = 0`; RLGC skipped). `mag_pathway_gaps.py prep`
now takes the EC numbers of a gap reaction also from the gapsmith Reactions.tbl (column `ec` via `dbhit`)
and from the ModelSEED Enzyme Class alias table (`--ec-from-reactions-tbl`, `--reaction-ec-aliases`; both
off by default so P14 is unchanged), and maps gapsmith gene products and hits, which are protein ids, to
loci directly. On the test archaeon this makes 172 of 237 gap reactions fillable (gapseq's own model of
the same bin: 94 of 214 through its SBML ECs); over the five genomes 60 to 172 gaps are fillable, RLGC
emits 60 to 166 suggestions, and the P14 score reports 46 to 109 new fills at q = 0.30.

Timing per genome on 4 CPUs: gapsmith 1.0 to 1.7 min, Stage 0 0.2 to 2.7 min, Stage 1 + 2A + 2B 1.1 to
2.0 min, prep and score seconds. Channel chunks of 3,000 representatives took 18 to 25 min on 8 CPUs
(InterProScan dominates), i.e. about 4 CPU-seconds per representative.

Observations on the inputs: isolate `MR59-1` (RAK-PR, Staphylococcus epidermidis, HQ by CheckM2 and
BUSCO) is an 8.1 Mb assembly in 528 contigs with 11,133 proteins, about three times the size of the
species; it is annotated as delivered, but its catalogue entry deserves a look.

## Reused code and data

`bgc_benchmark.py build-claims` (claims), `mag_pathway_gaps.py prep` and `score` (gap reactions and
gap-fill scoring, as P14), `stoichiometric_locus_refiner.py` (Stage 2B sidecar), `gspa-unified.jar`
v2.1.0, the P14 reference files (GO 2026-03-25, ec2go, pfam2go, KEGG, SEED tables, essential profiles,
experimental-GO Swiss-Prot DIAMOND database) and the universal-model data of P07. Checksums are in
`ref/sha256.txt` on ibex.
