# Workflow Notes

This module intentionally starts from raw-like toy input tables instead of a final ternary plotting table.

## Workflow

1. Read shared toy input tables:
   - `sample_metadata.csv`
   - `abundance_table.csv`
   - `taxonomy_table.csv`

2. Run sanity checks:
   - input files exist;
   - required columns exist;
   - sample IDs in metadata match abundance sample columns;
   - feature IDs in abundance match taxonomy;
   - the three selected axis groups exist in metadata.

3. Calculate relative abundance:
   - abundance values are divided by each sample total;
   - the resulting table is written to `results/relative_abundance_by_feature.csv`.

4. Aggregate to the selected taxonomy level:
   - default level is `Genus`;
   - features sharing the same taxon are summed.

5. Calculate group mean abundance:
   - sample-level relative abundance is averaged within each group;
   - the table is written to `results/group_mean_abundance_by_<level>.csv`.

6. Select top taxa:
   - taxa are ranked by overall mean relative abundance;
   - the top `top_n` taxa are retained.

7. Calculate ternary proportions:
   - mean abundance in the three selected groups is extracted;
   - the three values are normalized so each taxon sums to 1 across the ternary axes.

8. Draw the ternary plot:
   - the plot uses `ggtern`;
   - point size represents overall mean abundance;
   - point color identifies taxa;
   - PDF and PNG outputs are written to `figures/`.

## Why Not Use a Final `Ter.txt` Table?

A final ternary plotting table is convenient for manual figure production, but it hides important workflow steps. This demo recomputes the ternary plotting table from raw-like abundance and taxonomy inputs so users can understand and modify the full workflow.
