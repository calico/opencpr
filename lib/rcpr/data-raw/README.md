Order of operations for generating example data supporting network contextualization (`netcontextr`)

1. Create Human Protein Atlas data (don't trim since this is a small dataset)

```r
source("data-raw/example_protein_subcellular_localizations.R")
```

2. Create example string graph (trim heavily since this is a large dataset). Trim to interactions involving GWAS hits from NAFLD so we have a coherent example dataset

```r
source("data-raw/example_string_graph.R")
```

3. Create example GTEX gene expression data (subset since this is a large dataset). Trim based to genes present in example_string_graph and subset to 8 tissues

```r
source("data-raw/example_gtex_data.R")
```