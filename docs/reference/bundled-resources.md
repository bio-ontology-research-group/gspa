# Bundled resources

GSPA packages reference ontology structures, taxonomic hierarchies, anatomical
cross-products, complex classifications, and domain essentiality profiles
directly inside `gspa-core/src/main/resources/`. When running `gspa integrate`,
`gspa annotate`, or `gspa evaluate`, these bundled files allow full symbolic
consistency, coherence, and completeness checking without requiring external
database downloads or network access.

## Summary of bundled files

| Directory | File | Total lines | Data records | Primary source ontologies |
|---|---|---|---|---|
| `taxon-constraints/` | `go-taxon-constraints.tsv` | 26,912 | 26,911 | Gene Ontology (`go.owl`, `go-computed-taxon-constraints.obo`) |
| `taxon-constraints/` | `ncbi-taxon-hierarchy.tsv` | 1,819 | 1,818 | NCBI Taxonomy (`ncbitaxon.owl`) |
| `taxon-constraints/` | `ncbi-taxon-labels.tsv` | 1,290 | 1,289 | NCBI Taxonomy (`ncbitaxon.owl`) |
| `taxon-constraints/` | `PROVENANCE.md` | 21 | N/A | Documentation on upstream repository provenance |
| `anatomical-constraints/` | `go-anatomical-constraints.tsv` | 2,716 | 2,715 | GO, Uberon, Cell Ontology (CL), Plant Ontology (PO) |
| `anatomical-constraints/` | `uberon-po-taxon-constraints.tsv` | 8,332 | 8,331 | Uberon, Plant Ontology, Cell Ontology, GO-Plus |
| `anatomical-constraints/` | `uberon-po-hierarchy.tsv` | 7,723 | 7,722 | Uberon, Plant Ontology, Cell Ontology |
| `complex-terms/` | `complex_terms.tsv` | 159 | 158 | Gene Ontology (`GO:0032991`), EBI Complex Portal |
| `essential-functions/` | `bacteria.tsv` | 35 | 32 | JCVI Syn1.0 minimal genome, GAEF benchmark |
| `essential-functions/` | `archaea.tsv` | 29 | 26 | Archaeal core functional profile |
| `essential-functions/` | `eukaryote.tsv` | 34 | 32 | *S. cerevisiae* essential functional profile |
| `essential-functions/` | `virus.tsv` | 27 | 24 | Curated viral replication and virion hallmark functions |
| `config/` | `bacteria.yaml` | 31 | N/A | Bacterial isolate pipeline configuration |
| `config/` | `archaea.yaml` | 31 | N/A | Archaeal isolate pipeline configuration |
| `config/` | `eukaryote.yaml` | 29 | N/A | Eukaryotic pipeline configuration |
| `config/` | `mag.yaml` | 24 | N/A | Metagenome-assembled genome configuration |

## Taxon constraints (`taxon-constraints/`)

Loaded by `gspa.ontology.TaxonConstraints` and `gspa.ontology.SatConsistencyChecker`.
These tables encode which Gene Ontology terms are biologically permitted or
prohibited in specific taxonomic lineages.

### `go-taxon-constraints.tsv`

- **Record count**: 26,911 data rows (plus 1 header row).
- **Columns**: `GO_ID`, `Constraint_Type`, `Taxon_ID`.
- **Constraint types**:
  - `only_in_taxon`: 13,819 rows. Indicates that the GO term may only appear in
    the specified taxon lineage (e.g. `GO:0000001` only in `NCBITaxon_2759` *Eukaryota*).
  - `never_in_taxon`: 13,092 rows. Indicates that the GO term is forbidden in
    the specified taxon lineage (e.g. `GO:0000011` never in `NCBITaxon_4896` *Schizosaccharomyces pombe*).
- **Format example**:

```text
GO_ID	Constraint_Type	Taxon_ID
GO:0000001	only_in_taxon	NCBITaxon_2759
GO:0000011	never_in_taxon	NCBITaxon_4896
GO:0000015	only_in_taxon	NCBITaxon_131567
GO:0000022	only_in_taxon	NCBITaxon_2759
```

### `ncbi-taxon-hierarchy.tsv`

- **Record count**: 1,818 data rows (plus 1 header row).
- **Columns**: `Term`, `Relationship`, `Parent/Disjoint_From_Term`.
- **Relationships**:
  - `is_a`: 1,329 rows. Encodes taxonomic lineage subsumption.
  - `disjoint_from`: 454 rows. Encodes mutual exclusivity between taxa (for
    example, *Bacteria* and *Eukaryota* are mutually disjoint).
  - `union_of`: 35 rows. Encodes taxonomic grouping definitions.
- **Role**: Provides the backbone for SAT4J consistency verification. When two
  predicted GO terms mandate conflicting taxon lineages that are disjoint in this
  backbone, a taxonomic violation core is identified.

### `ncbi-taxon-labels.tsv`

- **Record count**: 1,289 data rows (plus 1 header row).
- **Columns**: `Taxon_ID`, `Label`.
- **Role**: Provides human-readable scientific names for NCBI Taxonomy identifiers
  (e.g. `NCBITaxon_2` $\to$ *Bacteria*, `NCBITaxon_2759` $\to$ *Eukaryota*),
  used in quality reports and provenance logs.

## Anatomical constraints (`anatomical-constraints/`)

Loaded by `gspa.ontology.AnatomicalConstraints` to enforce developmental and
structural consistency on multicellular eukaryotic annotations.

### `go-anatomical-constraints.tsv`

- **Record count**: 2,715 data rows (plus 1 header row).
- **Columns**: `GO_ID`, `GO_Label`, `Relation_ID`, `Relation_Label`, `Anatomy_ID`, `Anatomy_Label`.
- **Anatomy namespaces**:
  - `UBERON`: 1,333 axioms (Uberon multi-species anatomy ontology).
  - `CL`: 1,196 axioms (Cell Ontology).
  - `PO`: 186 axioms (Plant Ontology).
- **Primary relations**:
  - `results in development of` (`RO:0002296`): 756 axioms.
  - `results in morphogenesis of` (`RO:0002298`): 382 axioms.
  - `results in acquisition of features of` (`RO:0002315`): 340 axioms.
  - `occurs in` (`BFO:0000066`): 220 axioms.
  - `results in formation of` (`RO:0002297`): 205 axioms.
  - `part of` (`BFO:0000050`): 122 axioms.
- **Format example**:

```text
GO_ID	GO_Label	Relation_ID	Relation_Label	Anatomy_ID	Anatomy_Label
GO:0000902	cell morphogenesis	RO:0002298	results in morphogenesis of	CL:0000000	cell
GO:0001410	chlamydospore formation	RO:0002297	results in formation of	CL:0000726	chlamydospore
GO:0001501	skeletal system development	RO:0002296	results in development of	UBERON:0001434	skeletal system
GO:0001503	ossification	RO:0002297	results in formation of	UBERON:0001474	bone element
```

### `uberon-po-taxon-constraints.tsv`

- **Record count**: 8,331 data rows (plus 1 header row).
- **Columns**: `Anatomy_ID`, `Anatomy_Label`, `Constraint_Type`, `Taxon_ID`, `Taxon_Label`, `Source`.
- **Constraint types**: `only_in_taxon` (4,767 rows), `never_in_taxon` (3,564 rows).
- **Role**: Links anatomical structures to their phylogenetic distribution. For
  instance, *neural crest cell* (`CL:0000011`) is constrained to `NCBITaxon_33213`
  (*Bilateria*), preventing its annotation in plants, fungi, or non-bilaterian metazoans.

### `uberon-po-hierarchy.tsv`

- **Record count**: 7,722 data rows (plus 1 header row).
- **Columns**: `Child_ID`, `Relationship`, `Parent_ID`.
- **Relationships**: `is_a` (4,445 rows), `part_of` (2,130 rows), `develops_from` (1,147 rows).
- **Role**: Structural closure graph over Uberon, Cell Ontology, and Plant Ontology
  terms, allowing constraints asserted on anatomical superstructures or tissues
  to propagate downward to specific sub-structures and cell types.

## Protein complex terms (`complex-terms/`)

Loaded by `gspa.metrics.Coherence` and `gspa.integration.prior.CoherencePrior`.

### `complex_terms.tsv`

- **Record count**: 158 data rows (plus 1 header row).
- **Columns**: `GO_term`, `classification`, `definition`.
- **Classification values**:
  - `n` (obligate heteromeric complex): 70 terms. Requires two or more distinct
    gene products to form a functional quaternary complex (e.g. ribonucleases,
    DNA polymerase holoenzymes, SAGA complexes).
  - `h` (homomeric complex): 70 terms. Formed by multiple copies of a single gene
    product.
  - `a` (ambiguous / variable stoichiometry): 18 terms. Complex composition varies
    across taxa or conditions.
- **Role**: Used during coherence assessment to detect complex singleton false
  positives: if an obligate heteromeric complex (`n`) is triggered by only a
  single protein in the genome without partner subunits, it is flagged as an
  incoherence anomaly and subjected to demotion or partner recruitment.

## Essential function profiles (`essential-functions/`)

Loaded by `gspa.metrics.Completeness` and `gspa.integration.prior.EssentialityPrior`.
These TSVs define domain-wide core biological processes and molecular functions
that must be present in a self-sustaining organism.

- **Columns**: `GO_term`, `category`, `description`.
- **Profiles**:
  - `bacteria.tsv` (32 terms): Derived from the JCVI Syn1.0 synthetic minimal
    bacterial genome. Categories include `Core` (DNA replication, transcription,
    translation, protein folding), `Nutrient Uptake`, `Defense`, `Environmental Adaptation`,
    and `Glucose Metabolism`.
  - `archaea.tsv` (26 terms): Based on universal archaeal core processes. Categories
    include `Core`, `Archaeal-specific`, and `Energy Metabolism`.
  - `eukaryote.tsv` (32 terms): Derived from *Saccharomyces cerevisiae* essential
    gene knockout screens. Categories include `Core`, `Energy Metabolism`, and
    `Eukaryote-specific` (nuclear transport, vesicle trafficking, chromatin organization).
  - `virus.tsv` (24 terms): Hallmark viral functions covering `Core`, `Replication`,
    `Virion_Structure`, `Virion_Assembly`, `Maturation_Lysis`, and `Host_Interaction`.

## Regeneration of bundled resources

All bundled resources are extracted directly from official upstream OBO Foundry
ontologies using reproducible Python extraction scripts located in `scripts/`.

### Regenerating anatomical and taxon constraints

The script `scripts/extract_anatomical_and_taxon_constraints.py` (also referenced
as `scripts/extract_anatomical_constraints.py`) parses `go-plus.owl`, `uberon.obo`,
and `po.obo` to rebuild the anatomical constraints and NCBI hierarchy.

```bash
# Download upstream ontology files
wget http://purl.obolibrary.org/obo/go/extensions/go-plus.owl
wget http://purl.obolibrary.org/obo/uberon.obo
wget http://purl.obolibrary.org/obo/po.obo

# Run extraction script
python3 scripts/extract_anatomical_and_taxon_constraints.py \
  --go-plus go-plus.owl \
  --uberon-obo uberon.obo \
  --po-obo po.obo \
  --existing-taxon-hier gspa-core/src/main/resources/taxon-constraints/ncbi-taxon-hierarchy.tsv \
  --existing-taxon-labels gspa-core/src/main/resources/taxon-constraints/ncbi-taxon-labels.tsv \
  --out-dir gspa-core/src/main/resources
```

The script extracts cross-product definitions across 30 relation types, parses
UBERON and PO `is_a`, `disjoint_from`, and `relationship` stanzas, and updates
the output TSVs under `anatomical-constraints/` and `taxon-constraints/`.

### Regenerating complex terms

The complex terms classification table is derived from the Gene Ontology cellular
component hierarchy under `GO:0032991` (*protein-containing complex*) combined with
stoichiometric curation from the EBI Complex Portal.

```bash
# Extract complex descendants and align with Complex Portal curation
python3 scripts/extract_complex_terms.py \
  --go-obo go.obo \
  --complex-portal complextab.tsv \
  --out gspa-core/src/main/resources/complex-terms/complex_terms.tsv
```

The extraction script verifies all GO identifiers against `go.obo`, extracts
complex descriptions, checks stoichiometry evidence for multi-subunit requirements,
and stamps each complex with its classification (`n`, `h`, or `a`).
