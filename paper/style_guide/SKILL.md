---
name: write-edit-scientific-paper
description: Write, rewrite, or structurally edit scientific papers and manuscript sections, including abstracts, introductions, methods, results, discussions, and supplements. Use when asked to draft, revise, shorten, reorganize, polish, or review scientific prose; improve argument flow; remove markers of AI-written text; or make methods reproducible and results evidence-led. Applies to LaTeX, Markdown, plain text, and word-processing source.
---

# Write and edit scientific papers

Use this skill as both a scientific writer and a structural editor. Preserve the
science while making the paper's reasoning explicit, proportionate to the
evidence, and easy to follow.

Robert's primary editing criterion is the logical flow of the argument. First
decide which claims, premises, and intermediate steps belong and how they support
the conclusion. Paragraph shape, transitions, and sentence polish come afterward.
A smoother passage is not an improvement if it retains an irrelevant detour.

Before drafting or editing, read
`references/robert-scientific-style.md` completely. It is a critical adaptation
of the house style, not a list of mechanical substitutions. In particular, apply
its rules for ordinary prose paragraphs, connected reasoning, and corrections
without revision history. The headings and checklists in this skill organize
agent instructions; they are not a template for manuscript formatting.

## Scope and related skills

This skill owns manuscript prose, argument structure, and scientific writing style.
Use `grounded-research` for reference/identifier audits or requested project setup;
its scaffolder exports this skill's style guide. A prose edit does not require a
new project layout or a full audit of unchanged identifiers. Use `write-paper-review`
when the requested deliverable is a referee report, and `new-paper` for recording
and announcing an accepted or published paper.

## Establish the source of truth

1. Identify the audience, venue, length limit, manuscript role, and requested
   scope. Inspect venue instructions when they are available.
2. Locate the authoritative data, tables, figures, result files, citations, and
   terminology. Treat generated values and manuscript macros as fixed unless the
   user asks to recompute them.
3. Separate established facts, results shown by the current study, plausible
   interpretations, and speculation.
4. Never invent a result, citation, method detail, data version, or novelty
   claim. Use an explicit placeholder only in a working draft. If a final claim
   requires unavailable evidence or a scientific choice, ask the user.

## Learning and applying Robert's style

Apply the revision-derived explanatory guidance in the shared style reference.
For substantial structural rewrites or further style calibration, consult
[the historical revision evidence](references/revision-evidence.md). Prefer
before/after edits to coauthor prose over surface imitation of finished papers.
Preserve necessary premises, define operations, and make the reason for the next
step explicit. Historical examples guide structure; current instructions govern
formatting and wording. Do not claim a stylistic improvement from metadata checks
alone.

## Build the argument before polishing sentences

For the passage being drafted or revised, identify its role in the paper's
contribution. Scale this planning to the request: a local edit needs the surrounding
paragraphs, not a new outline of the entire paper. For each relevant section or
subsection, identify:

- the question or claim it advances;
- the evidence or method that supports it;
- why that evidence changes what the reader should believe; and
- why the next paragraph or section follows.

Remove material that does not support this chain, or move necessary technical
detail to the Methods or Supplement. Do not move details whose absence would
prevent reproduction.

Before polishing, outline the argument independently of the existing paragraph
boundaries: the question, necessary premises or observations, intermediate steps,
and supported conclusion. Test whether each piece of the draft contributes to
that argument. Remove a detour when its removal leaves the reasoning intact;
do not manufacture a bridge to retain it merely because it is accurate and related
to the topic. Preserve scientific qualifications and information required for
reproduction. Add a missing premise only when the sources support it; otherwise
flag the gap. Then arrange paragraphs around the resulting progression and check
their internal dependencies. These planning notes are not manuscript headings.
An explicit argument need not be longer: explain only the steps the intended
reader needs to follow this particular inference.

## Use section-specific structures

### Introduction

Build the Introduction as a problem-driven argument:

1. State the problem from general to specific. Define it precisely for the
   paper type. For a computational problem, identify the inputs, outputs, and
   constraints governing how the outputs relate to the inputs. Do not force
   this input--output formulation onto resource, observational, or other papers
   for which it does not fit.
2. Explain how the problem has been addressed. Organize the state of the art
   around methods or approaches, not names. Explain how each approach works at
   the level needed for the argument, which part of the problem it solves,
   where it fails, and how later approaches build on it. Continue until the
   limitations of the current state of the art are explicit.
3. State how the new work addresses those limitations. Introduce the work
   directly, then map each established limitation to the corresponding
   component or advance. Do not claim capabilities unrelated to the gap just
   established.
4. State how the work is evaluated and what scientific or practical impact
   that evaluation can support. Keep the impact claims within the scope of the
   evidence.
5. State where the method, resource, software, and data are available, as
   applicable.

Remove prior work that does not support this argument or move it to the point
where its role becomes clear. Do not turn the Introduction into an exhaustive
review, and state novelty only when the literature search supports it.

### Methods

Write the Methods so an informed reader can reproduce every reported result
without guessing. For each method, specify, as applicable, its purpose and
output; data or other inputs, sources, versions, access dates, sampling, and
inclusion or exclusion criteria; preprocessing and curation; protocol,
algorithm, or formal definition; software, instruments, databases, and models
with exact versions; parameters, prompts, seeds, configurations, and stopping
criteria; manual and automated steps with their decision rules; provenance
from inputs to outputs; and quality control, baselines, evaluation design,
metrics, and statistical analyses. State missing details of the reported
experiment when needed for reproducibility. Describe the method actually used;
keep the history of correcting its description out of the manuscript. Keep the essential design in
the main Methods; put commands, checksums, long parameter listings, and
diagnostics in the Supplement only when the main text and Supplement together
remain sufficient for reproduction.

Define every formal object unambiguously at first use. If several objects could
share a generic label such as `model`, give each a distinct name and state the
relationships or mappings among them. Do not use an unqualified generic term
when it could refer to more than one object.

### Results

Use the following reasoning sequence across a Results passage. Express it in
connected prose; the labels below are planning prompts, not visible headings or
a requirement for five sentences in every paragraph:

1. **Hypothesis or question:** state what the analysis tests or aims to show.
2. **How:** describe the comparison, experiment, or analysis briefly.
3. **Finding:** point to the relevant figure or table and state the pattern.
4. **Meaning:** state what the pattern supports, in positive terms (see
   "State what results show" below).
5. **Bridge:** when appropriate, identify the next uncertainty and use it to
   motivate the next hypothesis.

Do not narrate every value already visible in a table or figure. Report numbers
in prose only when the reader needs them to follow the inference. Distinguish
an accuracy test from an illustration, sensitivity analysis, negative control,
or mechanistic example.

### Discussion

Return to the paper's main claim, relate it to the closest prior work, delimit
what the evaluation establishes, state material limitations, and identify
specific next tests. Do not repeat the Results section as a list of values.

## Edit in five passes

1. **Structure:** repair the claim–evidence sequence and paragraph order.
2. **Evidence:** verify that claims match figures, tables, citations, and the
   evaluation design; weaken claims that exceed their controls.
3. **Reproducibility:** ensure that Methods and Supplement jointly specify the
   data, versions, protocol, metrics, exclusions, and code availability.
4. **Paragraphs:** check that each sentence develops the paragraph's main point
   and follows from information already introduced. Read adjacent paragraphs
   together: the next must take up a result, unresolved question, or dependency
   established by the previous one. Repair the reasoning before adding a transition
   word. Use a closing interpretation only when it advances the argument.
5. **Prose:** remove redundancy, nominalizations, vague intensifiers, inflated
   verbs, unnecessary hedges, and formulaic transitions.

Preserve labels, citations, macros, terminology, and numerical values during a
prose edit. Flag any edit that changes a scientific endpoint, comparison,
interpretation, or claim strength.

When Robert corrects a mistake, replace the mistaken account with the supported
account and update affected explanations within scope. The correction is an editing
instruction, not a new historical fact to narrate. Report any material change to
Robert separately from the manuscript. See the style guide's correction examples.

## State what results show

State positively what data, results, and analyses show. Never write a sentence
whose job is to tell the reader what the data or a result does *not* show,
mean, measure, or imply. Sentences of this kind pre-empt a misreading that the
writer has imagined for the reader; the house style does not make assumptions
about what the reader is thinking, so such sentences are removed rather than
softened.

Forbidden forms, in any section, in the main text, the supplement, and a data
paper alike:

- contrastive negation: "These shares describe the marker-gene profile, not
  cell counts"; "This is a selection outcome rather than a rate"
- negated inference: "does not show", "does not mean", "does not imply",
  "cannot be read as", "should not be interpreted as", "is not evidence of"
- disclaimers of scope attached to a result: "This test does not establish
  causation"; "These values are not comparable across trips"

Replace each with the positive statement that carries the same information:

- "These shares describe the relative abundance of the 16S marker gene."
- "Communities in this treatment were selected on nitrate."
- "The comparison is within trip; each trip was sequenced on one run."

If the positive statement is already in the text, delete the negated sentence.
If the positive statement conveys nothing the reader needs, delete the sentence
without replacement. A limitation belongs in the Discussion, phrased as what the
design does establish and what a named further experiment would add.

Apply this rule during the Evidence and Prose passes, and check for it
explicitly before delivering any draft or revision.

## Remove markers of AI-written prose

Look for generic scene-setting, repeated transition words, uniform paragraph
shapes, symmetrical lists without scientific need, restated captions,
unsupported claims of novelty or robustness, and sentences that sound polished
but do not advance the argument. Replace them with domain-specific subjects,
concrete actions, calibrated claims, and transitions that express the actual
logical relation. Vary sentence length naturally; do not introduce quirks or
informality merely to appear human.

## Check the actual output before delivery

Read the revised passage, not only the diff. Check these three recurring failures:

- Remove paragraph-opening bold labels and run-in headings, including Markdown
  `**Label.**`, LaTeX `\paragraph{...}` / `\paragraph*{...}`, and label-like
  `\textbf{...}`. Use ordinary prose paragraphs. Keep genuine section/subsection
  headings; use run-in headings only when Robert explicitly requests them or the
  venue requires them. Preserve mathematical bold notation.
- State each paragraph's main point and its relationship to its neighbour. If the
  relationship is only "another fact about the topic", rebuild the progression.
  Check sentence dependencies within paragraphs as well as paragraph order.
- Inspect "previously", "earlier version", "now", "corrected", and similar wording
  in context. Remove references whose only source is an editing exchange or a
  superseded draft. Preserve supported prior-work comparisons, reported version
  comparisons, and relevant scientific chronology.

For a local edit, also read the preceding and following paragraph. Do not claim
these checks passed for material outside the inspected scope.

## Deliver the revision

Lead with the revised text or completed artifact. Briefly report structural or
scientific changes that affect how the argument is read in a separate author-facing
note, never inserted into the manuscript. List unresolved
evidence gaps or decisions only when they require the user's input. Do not
explain every line edit unless requested.
