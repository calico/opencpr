test_that("Test Open Targets GraphQL API", {
  target_ensembl_gene_ids <- c("ENSG00000119718", "ENSG00000175354")
  open_targets_data <- summarize_open_targets_targets(target_ensembl_gene_ids)

  expect_equal(nrow(open_targets_data$targets), 2)
})
