# Format fixtures

These files support parser and injected-process tests. They lack verifiable
input, executable and run receipts and must not be described as observed tool
executions or biological validation. Some contain placeholder protein identifiers
and literal escaped tabs expanded by their tests.

Versioned HMMER output captured from an actual execution is separately stored in
`../observed-output/hmmer/`, with input sequences, model, commands and hashes.
Those sequences/model are explicitly synthetic; the capture proves a tool-format
and adapter contract, not accuracy on a biological protein family.
