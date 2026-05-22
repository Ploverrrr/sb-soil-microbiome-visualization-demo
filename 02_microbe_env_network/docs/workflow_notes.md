# Workflow Notes

This module builds a microbe-environment association network from raw-like shared toy inputs. It does not start from a precomputed network edge list.

## Workflow

1. Read shared toy inputs:
   - sample metadata;
   - environmental variables;
   - abundance table;
   - taxonomy table.

2. Run sanity checks:
   - input files exist;
   - required columns are present;
   - sample IDs match across tables;
   - abundance feature IDs are present in taxonomy;
   - selected environmental variables exist and are numeric.

3. Calculate relative abundance:
   - feature pseudo-counts are divided by each sample total.

4. Aggregate microbes:
   - features are grouped by the selected taxonomic level, default `Genus`;
   - taxon-level relative abundance is calculated per sample.

5. Select taxa:
   - taxa are ranked by overall mean abundance;
   - the top `top_n_taxa` taxa are retained.

6. Calculate associations:
   - each selected taxon is correlated with each selected environmental variable;
   - p values are adjusted with the selected method, default `BH`;
   - edges are retained if they pass both the correlation and adjusted p-value thresholds.

7. Build network tables:
   - retained associations become edges;
   - taxa and environmental variables become nodes;
   - taxonomic metadata and node type are added to the node table.

8. Draw network:
   - edge color shows positive or negative association;
   - edge width scales with absolute correlation;
   - node color separates taxa and environmental variables;
   - node size reflects degree.

## Interpretation

The network is an association map, not a causal model. It is useful for identifying candidate relationships among microbial taxa, contaminants, and soil chemistry variables.
