test_that("adjAF() works on ancestryData", {
  expect <- c(0.045000811, 0.219423999, 0.179859133, 0.006041453, 0.064062087)
  data("ancestryData")
  testval <- adjusted_data<-adjAF(data   = ancestryData,
                                  reference  = c("reference_AF_afr", "reference_AF_eur"),
                                  observed    = "gnomad_AF_afr",
                                  pi.target   = c(1, 0),
                                  pi.observed = c(.85, .15),
                                  adj_method = 'average',
                                  N_reference = c(704,741),
                                  N_observed = 20744,
                                  filter = TRUE)
  expect_equal(testval$adjusted.AF[1:5,"adjustedAF"], expected = expect)
})


