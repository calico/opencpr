#' Create igraph object from String network data
#'
#' `createStringGraph` uses a flat text file of gene/protein interactions,
#' downloaded from string.db to create netcontextr graph object
#'
#' @param stringNetwork a flat text file of interactions downloaded from string
#' @param primaryID desired primary ID to use in created igraph object,
#'  by default this is Ensembl gene ID
#' @param aliases key for converting between string protein ID and primary ID
#'
#' @return A netcontextr graph object, specifically, a list comprised of
#' "genes" (a tibble of the unique set of genes in the graph,
#' eventually this will also contain metadata for the genes), and
#' "interactions" (a multi-column tibble of the full set of interactions along
#' with the various weight metrics derived from string).
#'
#' @examples
#' if(interactive()){
#'   stringNetwork = "data-raw/9606.protein.links.full.v11.5.txt.gz"
#'   primaryID = "Ensembl_gene"
#'   aliases = "data-raw/9606.protein.aliases.v11.5.txt.gz"
#'
#'   createStringGraph(
#'     stringNetwork = stringNetwork,
#'     primaryID = primaryID,
#'     aliases = aliases
#'    )
#' }
#' @export
createStringGraph <- function(
  stringNetwork = "data-raw/9606.protein.links.full.v11.5.txt.gz",
  primaryID = "Ensembl_gene",
  aliases = "data-raw/9606.protein.aliases.v11.5.txt.gz"
  ){
  
  # Check for existence of string network file, throw error if doesn't exist
  checkmate::assertFileExists(stringNetwork)
  # Check for existence of string aliases file, throw error if doesn't exist
  checkmate::assertFileExists(aliases)
  
  # Load string data
  cli::cli_alert_info("Loading string network, this may take a little time...")
  string_protein_links <- readr::read_delim(stringNetwork, show_col_types = FALSE)
  
  # Load protein aliases from string metadata file
  # Note: would eventually like to find a way to silence the message from name_repair
  cli::cli_alert_info("Loading protein aliases")
  protein_aliases <- readr::read_tsv(aliases, name_repair = "universal", show_col_types = FALSE)
  
  # Filter to just ensembl gene conversions (data is in long format and there are multiple entries per string protein).
  # Default option is Ensembl_gene but in future versions could list the myriad of other options.
  cli::cli_alert_info(paste("Converting primary ID to", primaryID, sep = " "))
  
  # Confirm that primaryID exists
  checkmate::assertChoice(primaryID, protein_aliases$source)
  # Now filter to just these IDs
  protein_aliases <- protein_aliases %>% dplyr::filter(source == primaryID)
  # Check that there are not duplicated entries
  checkmate::assertCharacter(protein_aliases$.string_protein_id, unique = TRUE)
  # Use named vector as lookup table
  protein_alias_lookup <- protein_aliases$alias
  names(protein_alias_lookup) <- protein_aliases$.string_protein_id
  
  # Now extract Ensembl gene id and exchange with protein id
  string_protein_links <- string_protein_links %>%
    dplyr::mutate(
      protein1 = unname(protein_alias_lookup[protein1]),
      protein2 = unname(protein_alias_lookup[protein2])
      )
  
  # Assert that no NA were introduced
  checkmate::assertCharacter(string_protein_links$protein1,
                             any.missing = FALSE)
  checkmate::assertCharacter(string_protein_links$protein2,
                             any.missing = FALSE)
  
  # Create simplified tibble of interactions and tibble with corresponding weights
  cli::cli_alert_info("Creating edge tibble and adding weights")
  
  # Next create tibble containing edge weights and convert to a fraction by dividing by 1000
  interactions <- dplyr::mutate(string_protein_links, dplyr::across(.data$neighborhood:.data$combined_score, ~.x/1000))
  
  # Create tibble containing genes
  # Eventually will add metadata from protein_aliases to this (but see below)
  genes <- tibble::tibble(gene = unique(c(interactions$protein1, interactions$protein2)))
  
  # create graph object
  string_graph_full <- list(
    genes = genes,
    interactions = interactions
  )
  
  return(string_graph_full)
}

#' Load and format Human Protein Atlas subcellular localization data
#'
#' `load_and_clean_hpa_data` loads subcellular localization data from the
#' Human Protein Atlas and then creates a "compartment" column which is a
#' list of tibbles (all of the localisations of each protein are denoted in the
#' single column tibble for that protein). The "compartment" column can be used
#' to annotate the genes in a netcontextr graph object using the
#' `annotate_genes` function. After annotation, the netcontextr graph can then
#' be trimmed using the `trim_network_by_gene_attribute` function
#'
#' @param hpa_localization_data path to HPA subcellular localization data
#'
#' @returns a tibble containing all of the information from the HPA flat text
#'  file plus an addition "compartment" column, which is a list of tibbles, with
#'  each tibble containing a single column denoting all of the annotated
#'  subcellular locations (or a one row one column tibble with NA if none)
#'
#' @examples
#' if(interactive()){
#'   hpa_file <- "data-raw/subcellular_location.tsv.zip"
#'   hpa_localization_data <- load_and_clean_hpa_data(hpa_file)
#' }
#' @export
load_and_clean_hpa_data <- function(hpa_localization_data){
  # check for existence of human protein atlas file, throw error if doesn't exist
  checkmate::assertFileExists(hpa_localization_data)
  
  # load HPA data
  cli::cli_alert_info("loading Human Protein Atlas subcellular localization data")
  protein_subcellular_localizations <- readr::read_tsv(hpa_localization_data, name_repair = "universal", show_col_types = FALSE)
  # change name of Gene column to be more informative
  protein_subcellular_localizations <- dplyr::rename(protein_subcellular_localizations, ensembl_gene_id = .data$Gene)
  
  # convert GO.id column to a list of tibbles
  compartments <- stringr::str_split(protein_subcellular_localizations$GO.id, ";", n = Inf)
  compartments_list_tbl <- lapply(compartments, function(x) tibble::tibble(compartment = unlist(x)))
  protein_subcellular_localizations$compartment <- compartments_list_tbl
  
  return(protein_subcellular_localizations)
}

#' Load GTEX tissue specific expression data
#'
#' `load_and_clean_gtex_data` loads tissue specific expression data from GTEX
#' (median value per gene per tissue).
#'
#' @param gtex_data path to GTEX tissue specific expression data (medians)
#'
#' @returns a tibble containing all the information from the GTEX flat text file,
#' without modification (unlike related function `load_and_clean_hpa_data`)
#'
#' @examples
#' if(interactive()){
#'   # load gtex expression data
#'   gtex_file <- "data-raw/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"
#'   gtex_expression_data <- load_and_clean_gtex_data(gtex_file)
#' }
#' @export
load_and_clean_gtex_data <- function(gtex_data){
  # check for existence of gtex file, throw error if doesn't exist
  checkmate::assertFileExists(gtex_data)
  
  # load gtex data
  cli::cli_alert_info("loading GTEX tissue specific expression data")
  gtex_expression_data <- suppressMessages(readr::read_tsv(
    gtex_data,
    skip = 2,
    name_repair = "universal",
    show_col_types = FALSE
  ))
  
  # format gtex data
  # create ensembl_gene_id column by turning "hybrid" Ensembl gene IDs
  # in Name column into true Ensembl gene IDs
  gtex_expression_data <- dplyr::mutate(gtex_expression_data, ensembl_gene_id = stringr::str_remove_all(.data$Name, "\\.[0-9]+"))
  # also rename Name column to be more informative
  gtex_expression_data <- dplyr::rename(gtex_expression_data, ensembl_geneTranscript_id = .data$Name)
  # reorder columns so that ensembl_gene_id comes first
  gtex_expression_data <- dplyr::select(
    gtex_expression_data,
    .data$ensembl_gene_id,
    .data$ensembl_geneTranscript_id,
    .data$Description,
    tidyselect::everything()
    )
  
  return(gtex_expression_data)
}

#' Determine gene by tissue expression
#'
#' `gene_expression_by_tissue` uses GTEX tissue-specific gene expression data
#' loaded by `load_and_clean_gtex_data` along with the zFPKM normalization
#' algorithm/package to determine within which tissue (if any) each gene in
#' "expression_data" is expressed. More specifically, a node is considered to be
#' expressed if the corresponding gene's zFPKM value is greater than or equal to
#' -3. This cutoff is based on empirical data from ENCODE and corresponds to the
#' approximate threshold where the the ratio of active to repressed promoters
#' begins to be greater than 1. In other words, this algorithm attempts to
#' separate biologically relevant genes (associated with active promoters) from
#' ultralow-expression genes (associated with repressed promoters). -3 is a
#' liberal threshold, erring on the side of capturing some noisy genes rather
#' than missing active ones. After determining expression TRUE/FALSE this
#' function creates a "tissue" column which is a list of tibbles (all of the
#' tissues where a given gene is expressed are denoted in the single column
#' tibble for that gene). The "tissue" column can be used to annotate the genes
#' in a netcontextr graph object using the `annotate_genes` function. After
#' annotation, the netcontextr graph can then be trimmed using the
#' `trim_network_by_gene_attribute` function, preferably by specifying a
#' specific tissue of interest, since the filtering to all edges where the two
#' nodes overlap in one or more tissues is very memory intensive and will likely
#' crash (and is also not particularly meaningful)
#'
#' @param expression_data GTEX data loaded by `load_and_clean_gtex_data`
#' @inheritParams detect_unexpressed
#'
#' @returns a tibble containing all of the information from the GTEX flat text
#'  file plus an addition "tissue" column, which is a list of tibbles, with
#'  each tibble containing a single column denoting all of the tissues within
#'  which a gene is expressed (or a one row one column tibble with NA if none)
#'
#' @examples
#' if(interactive()){
#'   # load gtex expression data
#'   gtex_file <- "data-raw/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"
#'   gtex_expression_data <- load_and_clean_gtex_data(gtex_file)
#'   # determine gene expression by tissue using zFPKM
#'   gtex_expression_data_zFPKM <- gene_expression_by_tissue(gtex_expression_data)
#' }
#' @export
gene_expression_by_tissue <- function(
    expression_data,
    feature_defining_var = "ensembl_gene_id",
    invalid_vars = c("ensembl_geneTranscript_id", "Description")
){
  
  detected_expression <- detect_unexpressed(
    expression_data,
    feature_defining_var,
    invalid_vars
  )
  
  gene_by_tissue_expression <- detected_expression$is_expressed
  
  # now create a list of tibbles, one tibble per gene
  expressed_tissues <- apply(gene_by_tissue_expression, 1, function(x) colnames(gene_by_tissue_expression)[which(x)])
  expressed_tissues_list_tbl <- lapply(expressed_tissues, function(x) tibble::tibble(tissue = unlist(x)))
  expression_data$tissue <- expressed_tissues_list_tbl
  
  return(expression_data)
}

#' Annotate netcontextr graph with data from a list-col of tibbles
#'
#' `annotate_genes` modifies a netcontextr graph by adding annotation data to
#' the "genes" object. Currently this function assumes that the "field_name" in
#' "annotation_data" is a list of tibbles. This function calls the helper
#' function `extract_gene_annotations` which assumes that "annotation_data" has
#' a column called "ensembl_gene_id".
#'
#' @param netcontextr_graph a netcontextr graph object with "genes" and
#'  "interactions" elements, Gene ID needs to be in Ensembl format and specified
#'  in a column called "ensembl_gene_id"
#' @param annotation_data annotation data containing "ensembl_gene_id" and
#'  "field_name"
#' @param field_name the specific column in the annotation data to use for
#' annotating the genes in the netcontextr graph object. This column should be a
#' list of tibbles
#'
#' @returns modified netcontextr graph object, with a list of tibbles named
#' "field_name" added on to the "genes" element (each tibble containing a single
#' column denoting all of the annotated values in "annotation_data$field_name"
#' (or a one row one column tibble with NA if none)). The modified netcontextr
#' graph can now be trimmed using the `trim_network_by_gene_attribute` function
#'
#' @examples
#' # annotate the genes in a netcontextr graph with subcellular localization(s)
#' # add this information to a new field called "compartment"
#' annotated_string_graph <- annotate_genes(
#'   example_string_graph,
#'   example_protein_subcellular_localizations,
#'   "compartment"
#'   )
#'
#' # annotate the genes with expression data from GTEX
#' # first determine gene expression by tissue using zFPKM
#' gtex_expression_data_zFPKM <- gene_expression_by_tissue(example_gtex_data)
#' # second annotate genes with tissue expression(s)
#' # add this information to a new field called "tissue"
#' annotated_string_graph <- annotate_genes(example_string_graph,
#'   gtex_expression_data_zFPKM, "tissue")
#' @export
annotate_genes <- function(netcontextr_graph, annotation_data, field_name){
  
  # annotate each gene in the network with annotation data
  cli::cli_alert_info(
    "annotating {.field {field_name}} for each gene with annotation_data of shape {.field {nrow(annotation_data)}} x {.field {ncol(annotation_data)}}"
  )
  
  gene_annotations <- apply(
    netcontextr_graph$genes,
    1,
    function(x) extract_gene_annotations(x["gene"], annotation_data, field_name)
    )
  
  netcontextr_graph$genes <- netcontextr_graph$genes %>%
    dplyr::mutate("{field_name}" := gene_annotations)
  
  return(netcontextr_graph)
}


#' Extract gene annotations from annotation data
#'
#' `extract_gene_annotations` matches a gene to the corresponding tibble in
#' "annotation_data$field_name". This is a helper function for `annotate_genes`,
#' and assumes that "annotation_data" has a column called "ensembl_gene_id".
#'
#' @param gene gene to be annotated, in ensembl ID format
#' @param annotation_data annotation data containing "ensembl_gene_id" and
#'  "field_name"
#' @param field_name the specific column in the annotation data to use for
#' annotating the genes in the netcontextr graph object. This column should be a
#' list of tibbles
#'
#' @returns a single column tibble containing all of the annotations for gene
#' in "annotation_data$field_name"
#'
#' @examples
#' # add "compartment" annotation to a gene
#' example_gene_annotation <- extract_gene_annotations(example_string_graph$genes$gene[1],
#'   example_protein_subcellular_localizations, "compartment")
#' @export
extract_gene_annotations <- function(gene, annotation_data, field_name){
  # find gene in annotation data
  myIndex <- match(gene, annotation_data$ensembl_gene_id)
  
  if(is.na(myIndex)){
    field_annotation <- tibble::tibble("{field_name}" := NA)
    return(field_annotation)
  } else {
    # extract annotations for gene/protein
    field_annotation <- annotation_data[[myIndex, field_name]][[1]]
    
    return(field_annotation)
  }
}
