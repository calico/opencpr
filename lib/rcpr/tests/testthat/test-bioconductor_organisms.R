test_that("Species name to R package", {
  expect_equal(bioconductor_org_package_prefix("Homo sapiens"), "org.Hs.eg")
  expect_equal(bioconductor_org_package_prefix("Saccharomyces cerevisae"), "org.Sc.sgd")
})

test_that("Calling functions from org package works", {
  results <- bioconductor_org_function("ENZYME", "Saccharomyces cerevisae") # org.Sc.sgd is a suggested dependency
  expect_equal(colnames(results), c("systematic_name", "ec_number"))
})
