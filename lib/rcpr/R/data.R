#' EFO Diseases
#'
#' The Experimental Factor Ontology used for diseases and phenotypes by
#'   Open Targets.
#'
#' @format A tibble with 18,497 rows and 2 columns:
#' \describe{
#'   \item{id}{The EFO ID}
#'   \item{name}{The name of the EFO phenotype}
#' }
"efo_diseases"


#' Example protein network data from string.db
#'
#'
#' A subset of "full" protein links between human proteins, from string.db
#' NB full network means data includes distinction of direct vs. interologs
#' NB interologs = interolog is a conserved interaction between a pair of
#' proteins which have interacting homologs in another organism.
#' NB For most types of evidence, there are two types of scores:
#' the 'normal' score, and the 'transferred' score. The latter is computed from
#' data that is not originally observed in the organism of interest, but instead
#' in some other organism and then transferred via homology/orthology. All
#' potential source organisms are searched for evidence, but the actual transfers
#' to the receiving organism are made non-redundant (according to 'clades' of
#' closely related organisms in the tree of life).
#'
#' From String paper at https://academic.oup.com/nar/article/49/D1/D605/6006194
#'
#' All interaction evidence that contributes to a given network is benchmarked
#' and scored, and the scores are integrated into a final ‘combined score’. This
#' score is scaled between zero and one and provides an estimate of STRING’s
#' confidence on whether a proposed association is biologically meaningful given
#' all the contributing evidence. More specifically, The combined score is
#' computed by combining the probabilities from the different evidence channels
#' and corrected for the probability of randomly observing an interaction.
#' NB in the downloadable protein links files, the scores are multiplied by 1000
#' to make them integers.
#'
#' Confidence limits are as follows:
#' low confidence - 0.15 (or better) | or 150 or better in protein links files
#' medium confidence - 0.4 | or 400 or better in protein links files
#' high confidence - 0.7 | or 700 or better in protein links files
#' highest confidence - 0.9 | or 900 or better in protein links files
#'
#' Each association is provided as a connection
#' between two non-identical proteins, each from a different protein-coding gene
#' locus. STRING does not differentiate between splicing variants or
#' post-translationally modified protein isoforms encoded from the same locus—
#' instead, all such isoforms are collapsed and represented by a single,
#' canonical protein (i.e. a single protein per gene locus).
#'
#' The various evidence types that contribute to STRING associations are first
#' benchmarked and scored separately, in seven distinct evidence channels.
#'
#' The first three of the channels (neighborhood, fusion and co-occurrence)
#' contain computational association predictions based on whole-genome
#' comparisons. These so-called ‘genomic context’ channels can be computed for
#' all organisms for which a completely sequenced genome is available, and do
#' not depend on any further lab assays or measurements.
#' Neighborhood: two proteins are given an association score when their encoding
#' genes are in close proximity to each other on the chromosome. This channel is
#' applicable mostly for Bacteria and Archaea.
#' Fusion: STRING scans all genomes for open reading frames that appear to be
#' the result of gene-fusion events. For all inferred fusion events, the
#' constituent, non-fused genes in other genomes/organisms are given an
#' association score; the score is higher the better the fusion event can be
#' delineated in terms of the orthology of the participating genes. The
#' detection of gene fusion events across genomes can be used for predicting
#' functional associations of proteins, including physical interaction or
#' complex formation. See https://pubmed.ncbi.nlm.nih.gov/11820254/
#' Co-occurrence: STRING searches for pairs of genes whose occurrence patterns
#' throughout evolution show similarities. Such similarities can arise when
#' genes have been transferred, lost, or duplicated together during evolution,
#' which in turn can signify a shared function.
#'
#' The next two channels are dealing with functional genomics experiments or
#' direct lab assays.
#' Co-expression: STRING is collecting gene expression evidence from a number
#' of sources; this is then normalized, pruned, and the expression profiles over
#' a large variety of conditions are compared. Pairs of genes that show
#' consistent similarities between their expression profiles are assigned
#' association scores; the majority of the expression data is RNA-based, but we
#' also import proteome expression data, from the ProteomeHD database.
#' Experiments: STRING collects protein–protein interaction evidence from
#' experiments and assays in the lab. This includes biochemical, biophysical and
#' genetic experiments; all such interaction evidence is imported from the
#' curated interaction database organized in the iMEX consortium, plus BioGRID.
#'
#' The final two evidence channels deal with prior, consolidated knowledge on
#' protein–protein associations.
#' Database (AKA knowledge): STRING parses association evidence from curated
#' pathway databases, where it has been collected and consolidated manually by
#' expert curators. These include pathways annotated in KEGG, Reactome, and
#' MetaCyc, as well as protein complexes defined at the EBI Complex Portal or by
#' the Gene Ontology Consortium.
#' Textmining: STRING conducts statistical co-occurrence analysis across the
#' scientific literature. As of version 11.5 of STRING, the text-mining channel
#' is based on PubMed abstracts, articles from the PMC open access subset and
#' text from OMIM and SGD entry descriptions. Pairs of proteins mentioned
#' together in the same sentence, the same paragraph, or merely the same
#' publication are assigned a benchmarked association score, the calculation of
#' which is described in detail in (33).
#'
#' NB the string interactions were trimmed to edges where both genes have
#' subcellular localization data (see example_protein_subcellular_localizations)
#' and edges involving a gene with genetic association data (see
#' example_genetic_association_data).
#'
#' @format ## `example_string_graph`
#' A netcontextr graph object with 128 genes and 13,050 interactions.
#' Specifically, a list comprised of "genes" (a tibble of the unique set of
#' genes in the graph, eventually this will also contain metadata for the genes),
#' and "interactions" (a multi-column tibble of the full set of interactions
#' along with the various weight metrics derived from string (ranging from 0-1)).
#' \describe{
#'   \item{genes}{one column tibble of all genes in graph in ensembl gene id format}
#'   \item{interactions}{16 column tibble of all interactions in graph, first two
#'   columns are "gene 1" and "gene 2" in Ensembl gene id format, followed by
#'   14 different measures of interaction "score"}
#' }
#' @source <https://string-db.org/cgi/download.pl>
"example_string_graph"


#' Example RNA-Seq data (median TPM per tissue) from GTEX V8
#'
#' These are Analysis V8 data (dbGaP Accession phs000424.v8.p2) and correspond to
#' gene level median TPMs per tissue (54 tissues) summarized from ~17,000 samples.
#'
#' NB these data were trimmed to genes present in example_string_graph as well
#' as to 8 tissue to reduce example dataset size.
#'
#' @format ## `example_gtex_data`
#' A tibble with 11,471 rows and 11 columns:
#' \describe{
#'   \item{ensembl_gene_id}{Ensembl gene id, e.g. ENSG00000223972}
#'   \item{ensembl_geneTranscript_id}{GTEX's hybrid gene/transcript ID}
#'   \item{Decription}{human gene symbol, e.g. DDX11L1}
#'   \item{Liver}{median TPM for Liver}
#'   \item{Spleen}{median TPM for Spleen}
#'   \item{Lung}{median TPM for Lung}
#'   \item{Testis}{median TPM for Testis}
#'   \item{Ovary}{median TPM for Ovary}
#'   \item{Whole.Blood}{median TPM for Whole.Blood}
#'   \item{Kidney...Cortex}{median TPM for Kidney...Cortex}
#'   \item{Kidney...Medulla}{median TPM for Kidney...Medulla}
#' }
#' @source <https://www.gtexportal.org/home/datasets>
"example_gtex_data"


#' Subcellular localization data from The Protein Atlas
#'
#'
#' Subcellular location of proteins based on immunofluorescently stained cells.
#' The data is based on The Human Protein Atlas version 21.1 and Ensembl version 103.38.
#' See https://www.proteinatlas.org/about/assays+annotation#ihk for more
#' information on how these data were annotated / derived. Some specifics are
#' included below for reference.
#'
#' The Subcellular section of the Human Protein Atlas provides high-resolution
#' insights into the expression and spatiotemporal distribution of proteins
#' encoded by 13041 genes (65% of the human protein-coding genes). For each gene,
#' the subcellular distribution of the protein has been investigated by
#' immunofluorescence (ICC-IF) and confocal microscopy in up to three different
#' cell lines, selected from a subset of 36 of the cell lines found in the Cell
#' Line Section. Upon image analysis, the subcellular localization of the protein
#' has been classified into one or more of 35 different organelles and fine
#' subcellular structures. In addition, the section includes an annotation of
#' genes that display single-cell variation in protein expression levels and/or
#' subcellular distribution, as well as an extended analysis of cell cycle
#' dependency of such variations. See also key publication:
#' Thul PJ et al. (2017) “A subcellular map of the human proteome” Science 356
#' https://pubmed.ncbi.nlm.nih.gov/28495876/
#'
#' Annotation
#' In order to provide an interpretation of the staining patterns, all images
#' generated by ICC-IF are manually annotated. For each cell line and antibody,
#' the staining is described in terms of intensity, subcellular location and
#' single-cell variability (SCV). The staining intensity is classified as
#' negative, weak, moderate or strong based on the detector gain settings used
#' for image acquisition in combination with the visual appearance of the image.
#' SCVs within an immunofluorescence image are classified as intensity variation
#' (variation in their expression level) or as spatial variation (variation in
#' the spatial distribution).
#'
#' Knowledge-based Annotation
#' The knowledge-based annotation aims to provide an interpretation of the
#' detected subcellular localization of a protein. In the first step, stainings
#' in different cell lines with the same antibody are reviewed and the results
#' are compared with external experimental protein/gene characterization data
#' for subcellular localization, available in the UniProtKB/Swiss-Prot database.
#' In the second step, all antibodies targeting the same protein are taken in
#' consideration for a final annotation of the subcellular distribution of the
#' protein.
#'
#' Reliability Score
#' Each location is separately given one of the four reliability scores
#' (Enhanced, Supported, Approved, or Uncertain) based on available
#' protein/RNA/gene characterization data from both HPA and the
#' UniProtKB/Swiss-Prot database. The reliability score also encompasses several
#' additional factors, including reproducibility of the antibody staining in
#' different cell lines, correlation between staining intensity and RNA
#' expression levels, and assays for enhanced antibody validation. Enhanced
#' validation in achieved by using antibodies binding to different epitopes on
#' the same target protein (independent antibody validation), by assessing
#' staining intensity upon knockdown/knockout of the target protein
#' (genetic validation) and/or by matching of the signal with a GFP-tagged
#' protein (recombinant expression validation), and experimental evidence for
#' subcellular location described in literature. The individual location
#' relibility scores are summarized in an overall gene reliability score.
#'
#' There are four different reliability scores:
#' Enhanced - The antibody has enhanced validation and there is no contradicting
#' data, such as literature describing experimental evidence for a different
#' location.
#' Supported - There is no enhanced validation of the antibody, but the
#' annotated localization is reported in literature.
#' Approved - The localization of the protein has not been previously described
#' and was detected by only one antibody without additional antibody validation.
#' Uncertain - The antibody-staining pattern contradicts experimental data or
#' expression is not detected at RNA level.
#'
#' From Thul PJ et al. (2017) “A subcellular map of the human proteome” Science
#' main and additional locations are generated for each protein on the basis of
#' a clear difference either in the signal strength or in the occurrence across
#' the tested cell lines.
#'
#' NB these data were not trimmed at all
#'
#' @format ## `example_protein_subcellular_localizations`
#' A tibble with 13,105 rows and 15 columns:
#' \describe{
#'   \item{Gene}{Ensembl gene identifier}
#'   \item{Gene name}{name of gene}
#'   \item{Reliability}{gene reliability score}
#'   \item{Main Location}{Primary location(s) of the protein}
#'   \item{Additional Location}{Secondary location(s) of the protein}
#'   \item{Extracellular Location}{i.e. is the protein secreted}
#'   \item{Enhanced}{enhanced locations (strongest evidence)}
#'   \item{Supported}{supported locations (second strongest evidence)}
#'   \item{Approved}{approved locations (third strongest evidence)}
#'   \item{Uncertain}{uncertain locations (weakest evidence)}
#'   \item{Single-cell variation intensity}{locations with single-cell variation in intensity}
#'   \item{Single-cell variation spatial}{locations with spatial single-cell variation}
#'   \item{Cell cycle dependency}{locations with observed cell cycle dependency (type can be one or more of biological definition, custom data or correlation)}
#'   \item{GO id}{Gene Ontology Cellular Component term identifier}
#'   \item{compartment}{list-col containing a one column tibble with all of the subcompartments a given gene is observed in}
#' }
#' @source <https://www.proteinatlas.org/about/downloads>
"example_protein_subcellular_localizations"
