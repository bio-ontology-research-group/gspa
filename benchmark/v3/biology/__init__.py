"""Biology evaluation suite for real GSPA whole-proteome predictions.

Streams Stage 0 (raw predictor evidence without priors) and Stage 1 (full GSPA
with default priors) outputs, writes compact per-cell summaries with input
SHA-256 receipts, evaluates fixed and cited biological rules, and renders
side-by-side Markdown tables with plain-word prose.
"""
