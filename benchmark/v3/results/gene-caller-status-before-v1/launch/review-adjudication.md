# Sonnet 5.5 review adjudication

Read-only review completed. Native observation 10507 supersedes unobserved assumptions.

Accepted: output-before-wait timeout flaw; ignored wait return; no caller exception boundary; permissive partial mutation/order join; genuine zero/status ambiguity; reset typed status per invocation; missing-sequence coverage; possible exact-ID/header-metadata alias collision warrants native reproduction before shared-reader extraction.

Native facts: Pyrodigal3.7.1 meta short input yields header-only GFF+empty FAA exit0; single mode same input exit1; empty meta input exits0 with empty GFF (so reviewer assertion that all output always has a header is false). All-N meta input produces one partial all-X protein, not zero. No native Prodigal executable on worker PATH.

Not accepted as required scope: automatic runtime fallback after failure (changes semantics and can hide primary failures); different ID sorting (preserve valid established IDs/order); broad eukaryotic multi-row reader support (this helper is for single-row caller CDS models); arbitrary strict protein length relation without partial/stop/translation-table observations. NumberFormatException already extends IllegalArgumentException, so the alleged exception-type distinction is not a defect. Caller is used sequentially by the existing pipeline; no observed concurrent-sharing defect.

No production edits have yet been made.
