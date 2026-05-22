# Data Privacy

This repository is intended for public GitHub portfolio use. Public modules must not expose private research data, sensitive paths, real sample-level measurements, or manually prepared original result tables.

## Public Data Policy

All public data should be one of the following:

- simulated data;
- toy data;
- desensitized data;
- small synthetic examples designed to preserve analysis logic without exposing original measurements.

## Private Data Policy

The `_private_original/` directory is a git-ignored private reference area. It is excluded from the public project and should not be used as a public data source.

Do not publish:

- original sample metadata;
- original environmental measurements;
- original abundance tables;
- original differential analysis tables;
- original enrichment results;
- manually edited final plotting tables;
- private absolute paths;
- reference PDFs or third-party course materials copied into the private folder.

## Demo Input Policy

Each demo module should start from raw-like toy tables, such as sample metadata, environmental variables, abundance tables, taxonomy tables, or functional annotation tables. Any differential testing, enrichment, correlation, normalization, filtering, or model fitting required by a figure should be performed inside the module.
