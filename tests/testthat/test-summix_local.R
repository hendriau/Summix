test_that("summix_local() works on ancestryData", {
  data("ancestryData")
  testval <- summix_local(data = ancestryData,
                          reference = c("reference_AF_afr",
                                        "reference_AF_eas",
                                        "reference_AF_eur",
                                        "reference_AF_iam",
                                        "reference_AF_sas"),
                          NSimRef = c(704,787,741,47,545),
                          observed="gnomad_AF_afr",
                          goodness.of.fit = TRUE,
                          type = "variants",
                          algorithm = "fastcatch",
                          minVariants = 150,
                          maxVariants = 250,
                          maxStepSize = 1000,
                          diffThreshold = .02,
                          override_fit = FALSE,
                          override_removeSmallAnc = TRUE,
                          selection_scan = FALSE,
                          position_col = "POS") 
  
  expect_equal(length(colnames(testval$results)), expected = 12)
})