# Source-only publication scope

This repository is a sanitized snapshot, prepared separately from the working research directory. It includes source code, synthetic tests, configurations, plans, public upstream attribution and aggregate Markdown reports.

Excluded material includes raw and processed datasets, cell-level expression, target/gene data files, statistical banks, trained model binaries, predictions, submission packages, raw server responses, authentication stores, temporary HTTP configuration, signed URLs, logs, images, caches and nested repository histories.

Personal names, email addresses, local usernames, machine-specific filesystem paths, team identifiers and submission identifiers have been removed from the included text. Commit metadata uses a GitHub no-reply email. Public upstream author attribution and scientific dataset names remain.

Executable paths now use the project-local `.venv` or standard commands. H1 input paths use `data/h1_2025/`; optional legacy comparator paths use `data/legacy_model3/`. Numerical algorithms and model parameters were not intentionally changed. The full-axis writer tests use 18,533 synthetic gene identifiers instead of loading an external gene-list CSV. The included tests are run against this exported copy.

Historical plans and reports are records of their respective experiments, not promises that all artifacts needed for rerunning them are distributed here. Links to omitted artifacts are converted into explicit local-artifact references. Original freeze hashes remain historical evidence and are not hashes of privacy-edited documents.

The `.gitignore` and staged-content scanner reduce accidental disclosure in future commits. They cannot recognize every possible sensitive value; review new material before publishing it.

Publication validation: all **181 synthetic tests passed** in the exported copy. All Python sources and JSON configurations parsed, all retained local Markdown links resolved, the staged-content scanner passed, and comparisons against the locally stored VCC and GitHub credentials found no matches. Credential values and raw validation logs are not distributed.
