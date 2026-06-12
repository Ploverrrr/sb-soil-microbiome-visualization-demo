# Codex Working Rules

These rules define how future work should be done in this repository.

## Private Reference Boundary

- Do not copy files from `_private_original/` into public modules.
- Do not add `_private_original/` to git.
- Do not use private original files as public demo inputs.
- Treat `_private_original/` as read-only private reference material.

## Public Module Rules

- Each numbered folder must be independently runnable.
- Each module must use relative paths.
- Each module must start from raw-like toy input data in `data/toy/`.
- Each module must compute intermediate statistics inside the module.
- Precomputed plotting tables should not be used as primary inputs.
- Scripts should write reproducible outputs to `results/` and `figures/`.

## Implementation Rules

- Do not create a global all-in-one pipeline.
- Do not require modules to run in project-wide order.
- Preserve publication-style colors, grouping logic, and scientific narrative where appropriate.
- Prefer clear R scripts with a user-editable parameter section and explanatory comments.
