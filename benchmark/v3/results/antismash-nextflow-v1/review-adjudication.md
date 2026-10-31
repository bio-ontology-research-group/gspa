Sonnet 5.5 read-only review found no correctness defects. Accepted the new importer
performance finding: replaced per-CDS and per-protein linear contig searches with one
ID index. Existing shared region projection still uses linear contig lookup; this is
a separate scalability limitation and this change does not claim a metagenome-scale
performance result. Full runtime validation is required; review is not execution evidence.
