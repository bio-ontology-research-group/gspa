# Robert Hoehndorf scientific writing style: critical adaptation

Source consulted:
`/home/leechuck/Documents/papers/hott/WRITINGSTYLE.md`.

This reference preserves the useful house style while resolving its internal
conflicts and avoiding rules that can reduce accuracy or readability.

## Primary criterion: logical flow of the argument

Robert's explicit priority is the logical flow of the argument; shape comes later
and is secondary. Determine what the passage needs to establish, which premises
and observations support it, and how the conclusion follows. Select and order that
content before deciding on paragraph boundaries or polishing sentences.

For each statement, ask what necessary work it does in this argument. An accurate
fact about the same topic may still be a digression. If removing it leaves the
inference intact and loses no necessary qualification or reproducibility detail,
consider deleting or relocating it. Do not invent a motivation or build an extra
transition simply to retain peripheral material. Preserve necessary intermediate
steps; make missing ones explicit only when the evidence supports them.

Use explanation proportionately to the audience and inference. Defining objects
and operations should make the reasoning follow, without turning every local edit
into a tutorial. A concise deletion can improve logical flow more than several
new explanatory sentences. Evaluate the resulting argument before evaluating its
paragraph shape, rhythm, or resemblance to preferred wording.

## Apply directly

- Prefer active voice when it identifies the responsible actor or makes the
  action clearer. Use passive voice when the actor is irrelevant or the object
  of study should remain the grammatical focus.
- Give readers the context needed to interpret a claim before presenting the
  claim. Do not add a dependent clause mechanically.
- Prefer direct, technically precise sentences. Use a longer multi-clause
  sentence only when the clauses form one inference.
- Calibrate certainty: assert established facts and direct observations;
  use `may`, `can`, or `suggests` for limited or uncertain inferences.
- State what data and results show. Never state what they do not show, mean,
  measure, or imply, and never use "X, not Y" or "X rather than Y" to fence
  off an interpretation the reader has not been given. Write the positive
  statement that carries the information, or delete the sentence.
- Use sentence case for titles and headings unless the venue requires another
  convention.
- Never use em dashes. Replace them with commas, parentheses, colons,
  semicolons, or separate sentences according to the logical relation. In
  LaTeX, do not use `---`.
- Do not use `thus`; use `therefore`.
- Prefix semantic version tags with `v` in prose, for example `PhysioMap
  v1.1.1`.
- Use `record` only for an actual row or serialized data object. For scientific
  content, name the object precisely, for example `axiom`, `assertion`,
  `relation`, `variable`, `value`, or `derivative-sign value`.
- Define every symbol before its first use. When a compact list of formal
  constraints would require readers to decode several new objects at once,
  use a table with one concrete example for each constraint.
- Refer to supplementary material by the shortest specific item, such as a
  subsection, figure, equation, or table. Never use an unqualified `the
  Supplement`. Prefix supplementary item numbers with `S`, and order their
  first citations numerically, for example `Suppl. Fig. S1` before `Suppl.
  Fig. S2`.
- Prefer finite verbs over unnecessary gerunds and nominalizations.
- Preserve domain terminology such as ontology, axiom, embedding, semantic
  similarity, annotation, phenotype, and genotype.
- Use plain verbs where they are equally precise: `use` rather than `utilize`,
  `show` rather than `demonstrate`, `apply` rather than `employ`, `allow` rather
  than `enable`, `reduce` rather than `mitigate`, and `rank` rather than
  `prioritize`.
- Remove editorializers (`fortunately`, `interestingly`, `surprisingly`),
  assertions of obviousness (`clearly`, `obviously`), vague intensifiers
  (`very`, `quite`, `really`), and throat-clearing phrases such as `it is
  important to note`, `in order to`, `a number of`, and `due to the fact that`.
- Use transitions only when they state a real relationship: addition,
  contrast, specification, cause, or consequence.
- Support factual claims with appropriate citations or empirical evidence.
  Report data sources, versions, evaluation splits, sampling, metrics, and
  repository locations when relevant to reproduction.
- Do not print cryptographic digests in a human-facing manuscript or PDF.
  Link to a versioned, machine-readable checksum manifest or metadata file
  that identifies the relevant artifact instead.
- Treat figures as designed scientific arguments rather than collections of
  default plots. In a multi-panel figure, align the panels and use a consistent
  visual grammar for typography, line weights, spacing, panel labels, and
  colour. Give each panel a clear evidential role, remove dead space and
  draft-like annotations, and redesign or separate panels that do not form one
  coherent visual narrative.
- Treat a repository cited by a paper as part of the publication. Its default
  branch and landing page should present the current released resource, its
  semantics, supported operations, limitations, installation, reproducible
  examples, release artifacts, and citation. Keep development history in Git
  history or a changelog rather than making it the repository's public
  narrative. Before claiming availability, verify every link without
  authentication and ensure that the tagged release matches the manuscript.

## Paragraph form and connected reasoning

Write manuscript body text as ordinary prose paragraphs. Do not introduce each
paragraph with a bold label followed by text, such as `**Evaluation.** We ...`,
`\textbf{Evaluation.} We ...`, or `\paragraph{Evaluation}`. Removing the bold
alone is insufficient: turn the label into a sentence that advances the argument.
Use genuine section and subsection headings for sustained divisions of the paper.
Do not promote every paragraph label into a subsection. Preserve bold mathematical
notation, and follow a venue's required run-in format or Robert's explicit request.

Build each paragraph around one question, claim, or explanatory step. Its opening
sentence should locate that step in the argument. Develop it with the necessary
premises, mechanism, evidence, or example in an order the reader can follow. Introduce
an object before using it, and make the subject of a pronoun or comparison clear.
Interpret the evidence when an interpretation adds meaning; do not append a
restatement merely to give every paragraph a concluding sentence. Paragraphs can
vary in length and structure with their scientific purpose.

A shared topic is insufficient to connect paragraphs. Identify what the preceding
paragraph establishes and what the next needs to explain or test. Carry the relevant
object, result, or unresolved question into the next paragraph's opening sentence.
If the logical connection is missing, repair the order or add the missing premise.
Words such as "Furthermore" and "Moreover" cannot supply that reasoning. Do not
invent a causal link, motivation, or experiment merely to make a transition smooth.

Use a private one-sentence-per-paragraph outline to check the progression. Planning
labels such as question, method, finding, and meaning stay in that outline; they
must not appear as run-in labels in the manuscript. For a small revision, inspect
the neighbouring paragraphs so that the replacement still fits the passage.

Illustrative example, assuming these methods and results are supported:

Weak: "**Evaluation.** We divided the data by patient. Accuracy was 0.81.
**Robustness.** We tested a second hospital. Accuracy was 0.73."

Connected prose: "We divided the data by patient so that no patient contributed
samples to both training and testing. The model achieved an accuracy of 0.81 on
the held-out patients."

"To assess whether performance extended beyond the first hospital, we evaluated
the model on patients from a second hospital. Accuracy was 0.73 in this cohort."

The purpose and evaluation design must come from the study evidence; this example
is a pattern for relating supported facts, not permission to supply missing facts.

## Corrections belong in the corrected account

Treat Robert's correction as an instruction to revise the account. Write the
scientific statement that remains true after the correction. Update dependent
sentences, definitions, captions, and cross-references within the requested scope;
flag affected material outside that scope separately. Keep editorial history in
an author-facing change note, a response to reviewers, or version control when
appropriate. The paper's reader has not participated in the editing conversation.

Do not insert an old error, abandoned formulation, failed draft, or explanation of
why wording changed into the manuscript. This includes sentences such as "Previous
versions conflated A with B, whereas the current version separates them", "We now
correctly account for X", and "Following feedback, we clarified Y" when their only
basis is an editing exchange. Do not keep the same comparison after deleting just
"previously" or "now"; replace the whole claim with the supported scientific account.

Illustrative correction: Robert says, "The split was by patient, not by image."

Manuscript: "We partitioned the dataset by patient, with all images from each
patient assigned to the same partition."

Optional author-facing note: "Corrected the split description to patient-level."

The added detail about image assignment is appropriate only if the correction or
source material establishes it. Otherwise state only "We partitioned the dataset
by patient" and resolve any material uncertainty separately.

Retain real scientific history when it is part of the argument: cited prior work,
a documented comparison of released versions, the sequence of a reported experiment,
a correction notice, or a requested reviewer response. Distinguish those evidence-
based comparisons from the history of editing this manuscript. A factual correction
may expose a material methodological flaw; report it to Robert and describe its
scientific consequences accurately rather than concealing it as a prose change.

## Argument patterns

- **Introduction:** problem from general to specific → state-of-the-art
  approaches and successive limitations → solution mapped to those limitations
  → evaluation and impact → availability.
- **Methods:** purpose and output → fully specified inputs → reproducible
  protocol, tools, versions, parameters, and decision rules → provenance,
  quality control, evaluation design, and analyses.
- **Results:** hypothesis or question → brief test → finding with a figure or
  table pointer → interpretation → next hypothesis.
- **Discussion or conclusion:** contribution → relation to prior work →
  limitations → specific future tests.

The Results pattern above supersedes the source guide's weaker instruction to
restate a method and then list metrics.

## Preserve Robert's explanatory style

The revision-based profile below is supported by six inspected pre-2022 edits
in four coauthored manuscripts. It describes preferred explanatory moves, not a
fixed paragraph template or a reason to imitate dated vocabulary. Detailed before/
after evidence is recorded in this skill's `references/revision-evidence.md` and
`references/revision-evidence.json`.

- Preserve the premises of an inference even when shortening. Compress repeated
  setup or peripheral detail before removing the definition, reason, or intermediate
  step that makes the conclusion follow. Retain an already coherent explanation.
- Organize methods and related work around the objects represented, operations
  performed, and properties preserved. Introduce named methods as instances of
  that explanation; avoid a catalogue of names and features.
- Make the reason for each operation recoverable: identify its input, what it
  changes or preserves, and why that matters for the next step. Define a technical
  operation before claiming its benefit. Explain the same scientific object as
  it passes through different representations.
- Introduce technical operations in plain terms at first use, for example
  prediction as computing the consequences of an intervention and abduction as
  identifying an intervention that can explain observations.
- Narrate a contribution through concrete author actions in their scientific
  order: `we collected`, `we classified`, `we encoded`, `we projected`, and
  `we evaluated`. Do not replace this sequence with a dense inventory of
  artifact properties.
- Build conceptual distinctions one at a time and preserve their dependency
  order. State what an object is, why the distinction matters, and only then
  introduce the notation or operation that depends on it.
- In related work, first state what an approach represents or computes, then
  state the limitation that matters for the paper. Avoid compressed contrasts
  that require the reader to reconstruct the missing premise.
- Repeat the precise technical noun when a pronoun or stylistic synonym could
  blur its referent. In particular, keep ontology semantics, causal semantics,
  representation, trait, and intervention distinct.
- Separate the content and semantics of a scientific resource from the
  evidence and provenance used to construct it. Describe admission evidence
  and provenance in Methods; do not present evidence classes as content types
  or as inputs to an application that does not read them. Organize application
  Results around the represented content that actually changes inference.
- When a formalism has several semantic or computational layers, introduce
  them in dependency order. First state what the primary variables denote and
  what the relations constrain. Then define any derivative, sign, symbolic, or
  other abstraction of that semantics. Finally state exactly which abstraction
  an implemented solver reads. Do not let an operational graph or calculus
  stand in for the underlying formalism.
- Keep explanatory prose at the level of scientific objects, formal
  parameters, and method concepts. Do not expose internal source-code
  identifiers, function arguments, schema field paths, enum labels, or fixture
  names as if they were scientific terminology. Use human-readable names
  instead. Retain public artifact names, stable content identifiers, external
  accessions, and exact commands only where they are needed for audit or
  computational reproduction.
- Use examples to perform explanatory work. Show how the defined operation acts
  on a concrete instance, then use that instance to explain its consequence.
  Carry a useful example across adjacent paragraphs rather than replacing it at
  every step. Supply only examples and mechanisms supported by the sources.
- Give each evaluation a question and state what the result establishes. Explain
  why a constraint changes the design or the metric. A table reference alone is
  not a finding. When one result raises a further question, make that question
  explicit before introducing the next comparison; do not invent experimental
  motivation or causation to complete this pattern.
- Use simple signposting such as `First`, `We then`, `Additionally`, `Among
  these`, and `However` when it exposes the actual argumentative sequence.
  Avoid varying transitions merely for style.
- Treat an application such as rare disease as an application of the general
  representation unless the evidence justifies making the application itself
  the paper's identity.
- Preserve this explanatory pacing while correcting draft-level grammar,
  repetition that adds no conceptual clarity, placeholders, and uncalibrated
  words such as `trivially`.

## Apply with judgment

- Do not target an average sentence length of 25–35 words. Vary length according
  to the inference and split sentences when this improves comprehension.
- Do not meet a numerical citation quota. Cite claims at the density required
  by the field and avoid citation padding.
- Do not add decimal places unless the measurement, sample size, and source
  justify that precision.
- Vocabulary preferences are defaults, not substitutions. For example,
  `workflow`, `pipeline`, `limit`, `hinder`, `allow`, and `enable` have
  different meanings. Choose the accurate word.
- Avoid promotional uses of `comprehensive`, `novel`, `unique`, `rigorous`, and
  `robust`, but retain them when they denote a defined property or metric, such
  as a unique top-ranked candidate.
- Do not use stock novelty formulas such as “to the best of our knowledge” or
  “the first” unless a current literature search supports the claim.
- Never leave invented numerical or citation placeholders in a final
  manuscript. Placeholders are acceptable only in an explicitly marked draft.
- Do not explain every correction unless the user requests a teaching edit.
  Report material structural and scientific changes instead.
- Technical precision and readability are joint requirements. When they appear
  to conflict, revise the structure or add necessary context rather than
  sacrificing either one by rule.
