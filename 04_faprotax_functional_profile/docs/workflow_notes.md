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

4. Aggregate function profiles:
   - counts are summed by sample and pathway;
   - sample-wise count-weighted relative abundance is calculated from the aggregated counts;
   - unweighted feature-assignment percentage is calculated as the percentage of active features assigned to each function in each sample.

5. Select top functions:
   - functions are ranked by overall mean unweighted feature-assignment percentage;
   - the top functions are retained for plotting.

6. Summarize by group:
   - mean, standard deviation, standard error, and sample count are calculated by group and function.

7. Draw figures:
   - bubble profile: sample-by-function view of unweighted feature-assignment percentage;
   - scaled bubble profile: count-weighted function profiles are z-score scaled across samples, following the original script style;
   - grouped barplot: group mean unweighted feature-assignment percentage with error bars.

## Interpretation

The plots summarize simulated functional patterns across groups. They are meant to demonstrate workflow structure and publication-style visualization, not to report real biological results.
