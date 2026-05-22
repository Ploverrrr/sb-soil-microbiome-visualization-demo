# Workflow Notes

This module builds a FAPROTAX-style functional profile from raw-like shared toy functional annotation data. It does not use manually edited group summary tables.

## Workflow

1. Read shared toy inputs:
   - `sample_metadata.csv`;
   - `functional_annotation_table.csv`.

2. Run sanity checks:
   - input files exist;
   - required columns are present;
   - functional sample IDs match metadata sample IDs;
   - selected annotation sources exist;
   - `top_n_functions` is valid.

3. Filter functional annotations:
   - keep selected annotation sources;
   - default sources are `FAPROTAX`, `SulfurCycle`, `NitrogenCycle`, and `CarbonCycle`.

4. Aggregate counts:
   - counts are summed by sample and pathway;
   - sample-wise relative abundance is calculated from the aggregated counts.

5. Select top functions:
   - functions are ranked by overall mean relative abundance;
   - the top functions are retained for plotting.

6. Summarize by group:
   - mean, standard deviation, standard error, and sample count are calculated by group and function.

7. Draw figures:
   - bubble profile: sample-by-function view of relative abundance;
   - grouped barplot: group mean relative abundance with error bars.

## Interpretation

The plots summarize simulated functional patterns across groups. They are meant to demonstrate workflow structure and publication-style visualization, not to report real biological results.
