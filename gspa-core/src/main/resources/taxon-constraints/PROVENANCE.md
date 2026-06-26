# Taxon constraint data (vendored)

Source: https://github.com/bio-ontology-research-group/genome-scale-pfp-adjust
(commit 9861e8d), companion code for *Genome-scale protein function
adjustment using constraint optimization* (A. Toonsi et al.).

- `go-taxon-constraints.tsv` — `data/go_taxon_constraints_updated.tsv`:
  GO_ID, Constraint_Type (only_in_taxon|never_in_taxon), Taxon_ID (NCBITaxon_N).
  ~26.9k constraints extracted from the Gene Ontology.
- `ncbi-taxon-hierarchy.tsv` — `data/taxon_hierarchy.tsv`:
  Term, Relationship (is_a|disjoint_from|union_of), Parent/Disjoint/Member.
  The NCBI-taxonomy disjointness backbone (kingdom-level + major eukaryote
  groups) used to decide whether two taxon requirements can co-occur.

Loaded by gspa.ontology.TaxonConstraints (constraints) and
gspa.ontology.SatConsistencyChecker (hierarchy + explicit disjointness).
The backbone covers kingdom-level taxa (Bacteria 2, Archaea 2157,
Eukaryota 2759, Viruses 10239, cellular organisms 131567, Metazoa, Fungi,
Viridiplantae, ...). Constraints on finer taxa not in the backbone are
loaded but only enforced when their lineage reaches a backbone node.
