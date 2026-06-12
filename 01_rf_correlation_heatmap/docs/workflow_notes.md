# Workflow Notes

This module recomputes all analysis results from raw-like shared toy inputs. It does not use hand-edited plotting tables or precomputed RF importance tables.

## Workflow

1. Read shared toy input tables:
   - sample metadata;
   - environmental variables;
   - feature abundance;
   - taxonomy;
   - functional annotation.

2. Run sanity checks:
   - input files exist;
   - required columns exist;
   - sample IDs match across metadata, environmental variables, and abundance columns;
   - feature IDs match between abundance and taxonomy;
   - user-selected environmental variables exist;
   - `top_n_features` is valid.

3. Prepare feature abundance:
   - if `feature_source = "taxonomy"`, calculate relative abundance from the feature-by-sample abundance table and aggregate to `target_taxonomic_level`;
   - if `feature_source = "function"`, aggregate functional counts by pathway and sample, then calculate relative abundance.

4. Select top features:
   - rank features by overall mean abundance;
   - retain the top `top_n_features`.

5. Calculate feature-environment correlations:
   - run `cor.test` for each selected feature and each selected environmental variable;
   - adjust p values using the selected method, default `BH`;
   - save the full correlation table.

6. Train random forest:
   - use selected feature abundances as predictors;
   - use `target_environmental_variable_for_rf`, default `Sb_total`, as the response;
   - extract feature importance from the fitted model;
   - report out-of-bag pseudo R2 and MSE from the model object.

7. Draw figures:
   - correlation heatmap with cyan-white-pink diverging colors;
   - RF importance bubbles overlaid on the target environmental variable column;
   - a separate RF importance bar plot;
   - a combined portfolio figure with the correlation heatmap on the left and RF importance on the right.

The combined figure is the primary display output. The standalone heatmap and standalone RF plot are auxiliary outputs that make it easier to inspect each analysis component separately.

Both panels in the combined figure use the same selected features and the same feature order. Feature labels are shown on the heatmap side only to avoid duplicate labels.

## Interpretation

The correlation heatmap summarizes pairwise associations. The RF model evaluates which selected features are useful predictors of the target environmental variable when considered together. These are complementary views and should not be interpreted as causal evidence.
