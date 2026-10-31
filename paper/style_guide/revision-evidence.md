# Evidence from Robert's manuscript revisions

Initial local analysis, 2026-09-30. This reference supports the scientific-writing
skill with observed revision choices. Current explicit preferences take precedence
over historical wording. The evidence concerns recurring editorial choices; it does
not establish that these features uniquely identify Robert against other writers.

## Sampling and attribution

- Inspected 603 immediate Git repositories under `~/Documents/papers`; 290 had commits before 2022. These counts include duplicates and non-paper projects, not 290 independent papers.
- Extracted eight historical manuscript snapshots across conceptual, resource, and machine-learning work. Selected source files rather than latexdiff files, templates, or cover letters.
- Following Robert's suggestion, shifted the main evidence to his edits of coauthor text. Screened 334 moderate-sized revisions in five repositories: `machine-learning-with-ontologies`, `flopo`, `el-embedding`, `function-prediction-paper`, and `self-normalizing-ontologies`.
- A coarse screen retained 297 changes with at least 25 before/after whitespace-delimited tokens and non-identical whitespace-normalized text. 190 touched lines last attributed to another identity. These are candidate counts: they include tables, comments, technical fixes, and near-reflow, not 190 demonstrated stylistic rewrites.
- Closely inspected selected diffs and retained six examples from four coauthored manuscripts. The selected source passages are in [revision-evidence.json](revision-evidence.json), including exact before/after lines, full commit and parent IDs, dates, and last-touch author counts.
- Robert identities were recognized by his name or known leechuck addresses; predecessor identities were checked in parent blame. Student/coauthor status is not inferred solely from a Git identity. Anonymous/import commits and large deletions were not treated as clean authorship proof.
- Use `git blame -w -M` to reduce whitespace and within-file movement artifacts. Git records who last touched lines, not who first formulated the ideas. Robert may commit collaboratively supplied wording. A rewrite can also change scientific content; separate those changes from style.
- Initial-import snapshots, unfinished drafts, and formatting changes have lower evidential weight. In particular, self-normalizing commit `07b52abfce924672400a9d7875adfa9fa3f436c0` is mostly line wrapping, pronoun changes, and removal of comments; it is not counted as a substantive reasoning rewrite.
- This sample is intentionally rich in ontology/representation-learning writing. It is a first profile, not an exhaustive survey or a blind predictive test.

## 1. Replace a catalogue with an explanatory sequence

Source: `machine-learning-with-ontologies/main.tex`, 2020-04-26, commit `ae7a76a85ff866da781ea1b4ef734a8781a89aa8`.
Parent `dc8eba71a253297d8d54f26bd73135a6288bcb50`. Selected before lines 730–755; after lines 766–771.

The earlier passage enumerates Walking RDF and OWL, OWL2Vec, and graph walking under bold run-in labels. Robert removes those labels from active prose, moves graph-construction material to the earlier representation discussion, and opens with the property graph embeddings preserve. He then describes graph construction, neighbourhood walks, and embedding in dependency order before introducing DeepWalk as an instance.

Learn the organization by representation and operation, with methods introduced as examples of that organization. Removing bold formatting alone would miss the substantive edit.

## 2. Explain why a representation requires the next operation

Source: `machine-learning-with-ontologies/main.tex`, 2020-04-27, commit `8096c0bfa5ca8e6b0ab724e0a3316d76af8e0bc9`.
Parent `344a7047f9526bada89abc106ae01f801cb4336c`. Selected before lines 1000–1020; after lines 986–1004.

The earlier OPA2Vec account lists formal/informal content, annotation axioms, sentences, and a pretrained language model. Robert explains the annotation-to-sentence mapping with a class/label example, then distinguishes where class identifiers obtain meaning from where natural-language terms obtain meaning. That distinction motivates biomedical-language pretraining.

Describe an operation through its input, transformation, and consequence. When moving to the next operation, state the specific requirement it satisfies. Preserve the meaning of the scientific objects across representations.

## 3. Define the operation before explaining its benefit

Source: `self-normalizing-ontologies/main.tex`, 2020-04-19, commit `c41cb6a32436540d6fa69cc450c849e2637c13c1`.
Parent `4f962cabef05a22eca62dad0e99ba1641db426b1`. Selected before lines 503–533; after lines 515–535.

The earlier Results text mixes normalization, task specificity, prior work, and expected prediction benefits. Robert defines normalization as making the symbols referring to an ontology class identical, then describes replacing mentions with the class IRI and how the recognition dictionary is constructed. His adjacent comment requests an explanation of the three embedding corpora and a worked biological example.

Build the account from a precise operation through its inputs and outputs to the expected effect. Use an example to expose the operation. This draft still contains incomplete citations and TODOs; those are evidence of his requested reasoning, not prose to imitate.

## 4. Give a comparison a question and explain its consequences

Source: `self-normalizing-ontologies/main.tex`, 2020-04-20, commit `f14b02d95e980d01fa589bb1aeae7a2d7d3f43d5`.
Parent `46840f2d0c0d84b649412d595cac347a0d818f2f`. Selected before lines 853–864; after lines 867–879.

The earlier passage says that a comparison is also performed. Robert makes its purpose explicit: testing improvement over a relevant baseline in toxicological-effect prediction. He connects the absence of comparator confidence scores to the reported metrics and adds an author comment asking for the finding after the table reference. Elsewhere in the same revision he removes an unsupported mechanism claim about converting noise into useful information.

A comparison should identify the question, relevant comparator, evidence, and finding. Explain why a design constraint affects the comparison. State what was learned; do not stop at a pointer to results or supply an untested explanation.

## 5. Integrate evaluation into the argument

Source: `function-prediction-paper/function.tex`, 2017-08-14, commit `46fe61f570c4be3cff9a8268e1e227352a0ffe7e`.
Parent `43826cd0eb9e58b8f3bcb39b540d2229a13730a2`. Selected before lines 644–648; after lines 550–576.

The earlier text isolates method comparison and sequence clustering under separate subsections. Robert moves acquisition details to Methods and integrates evaluation paragraphs into the Results flow. For the sequence split, he explains the similarity issue, the grouping rule, the resulting partition, and the scenario being evaluated. He also distinguishes metric-specific wins and losses in the baseline comparison.

Organize around what each evaluation establishes. Explain the design property needed to interpret its result. This is also a factual revision to the split description; it cannot be classified as style alone.

## 6. Preserve the premises of a compact derivation

Source: `el-embedding/main.tex`, 2019-05-19, commit `2320308d7481b265f923721942d2d78d647d4478`.
Parent `3962a46e7b245cbe422bb35fa558db358cb7211b`. Selected before lines 590–602; after lines 590–608.

The student passage already derives runtime from normalization, per-axiom cost, and training epochs. Robert keeps that chain, clarifies symbols and asymptotic notation, and connects the conclusion with Consequently and Therefore.

Do not claim Robert introduced the reasoning here. The evidence supports preserving and sharpening an existing derivation, and permitting explicit consequence markers when premises justify them. This is corroborating evidence, weaker than the structural rewrites above.

## Reproduce and extend the analysis

For any case, read the parent and revised manuscript at the recorded commits,
including neighbouring paragraphs; current working-tree line numbers may differ.

```bash
git -C ~/Documents/papers/machine-learning-with-ontologies diff --word-diff=plain ae7a76a85ff866da781ea1b4ef734a8781a89aa8^ ae7a76a85ff866da781ea1b4ef734a8781a89aa8 -- main.tex
git -C ~/Documents/papers/machine-learning-with-ontologies blame -w -M -L 730,755 ae7a76a85ff866da781ea1b4ef734a8781a89aa8^ -- main.tex
```

For further calibration, sample unseen Robert revisions from additional papers.
Keep each paper wholly in either the development or evaluation sample. Give the
writer only the before passage, its context, and the scientific evidence it needs;
withhold Robert's after version. Compare the resulting revision with his edit by
argument order, explicit premises, referent continuity, example placement, and
scientific fidelity. Compare against the previous skill on the same passages.
Use Robert's assessment to resolve tradeoffs. Exact phrase matching would reward
copying vocabulary and miss the feature this profile is intended to capture.
This paired forward evaluation has not yet been run.

## What should not be learned from this sample

- Older drafts contain `novel`, `unsurprisingly`, negated contrasts, draft labels, missing citations, and grammatical errors. Preserve the explanatory structure while applying Robert's current preferences.
- Preserve supported qualifications and actual scientific history when relevant. Do not convert the history of editing a draft into manuscript content.
- Neither maximal expansion nor maximal compression describes the edits: some explain a mechanism at greater length; others remove duplication and move details. Retain the premises needed for the inference.
- The revised runtime paragraph confirms that a student's already coherent structure can be retained. Do not force a rewrite merely to make a passage look different.
- Before/after examples can reveal preferred interventions; they do not establish that every scientific change made in an old draft is correct. Do not reuse the examples' factual claims in other papers.
