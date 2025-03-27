# code to download and prepare `example_string_graph` dataset
# This is a subset of "full" protein links between human proteins, from string.db
# NB full network means data includes distinction of direct vs. interologs
# NB interologs = interolog is a conserved interaction between a pair of
# proteins which have interacting homologs in another organism.
# NB For most types of evidence, there are two types of scores:
# the 'normal' score, and the 'transferred' score.
# The latter is computed from data that is not originally observed in the
# organism of interest, but instead in some other organism and then transferred
# via homology/orthology. All potential source organisms are searched for
# evidence, but the actual transfers to the receiving organism are made
# non-redundant (according to 'clades' of closely related organisms in the tree
# of life).
# see also documentation for "example_string_graph" in "/R/data.R"

suppressPackageStartupMessages(library(tidyverse))
library(igraph)

STRING_NET_URL <- "https://stringdb-static.org/download/protein.links.full.v11.5/9606.protein.links.full.v11.5.txt.gz"
STRING_NET_REL_PATH <- "data-raw/9606.protein.links.full.v11.5.txt.gz"
STRING_ALIASES_URL <- "https://stringdb-static.org/download/protein.aliases.v11.5/9606.protein.aliases.v11.5.txt.gz"
STRING_ALIASES_REL_PATH <- "data-raw/9606.protein.aliases.v11.5.txt.gz"
DEGREE_CUTOFF <- 4000

# download data from string
# this is protein network data (full network, incl. distinction: direct vs. interologs)
if(!file.exists(STRING_NET_REL_PATH)){
  download.file(STRING_NET_URL, STRING_NET_REL_PATH)
}

# next download and open key for converting ensembl protein ID to ensembl gene ID
if(!file.exists(STRING_ALIASES_REL_PATH)){
  download.file(STRING_ALIASES_URL, destfile = STRING_ALIASES_REL_PATH)
}

# now convert string flat text file to igraph object
example_string_graph <- createStringGraph(
  stringNetwork = STRING_NET_REL_PATH,
  primaryID = "Ensembl_gene",
  aliases = STRING_ALIASES_REL_PATH
  )

# now trim string data
example_string_graph <- trim_network_by_nodes(
  example_string_graph,
  rcpr::example_protein_subcellular_localizations$ensembl_gene_id,
  requireBothNodes = TRUE
  )

# filter by degree and then take the major component (to create a small, example graph)

string_ig <- igraph::graph_from_edgelist(
  as.matrix(example_string_graph$interactions[, c("protein1", "protein2")]),
  directed = FALSE
)

valid_nodes_by_degree <- igraph::V(string_ig)[igraph::degree(string_ig) > DEGREE_CUTOFF]$name
high_degree_subgraph <- igraph::subgraph(string_ig, valid_nodes_by_degree)

example_string_graph <- trim_network_by_nodes(
  example_string_graph,
  valid_nodes_by_degree,
  requireBothNodes = TRUE
  )

# update file in /data
usethis::use_data(example_string_graph, overwrite = TRUE)
