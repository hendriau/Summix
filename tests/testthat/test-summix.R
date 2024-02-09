test_that("summix() works on ancestryData", {
  data("ancestryData")
  testval <- summix(data = ancestryData,
                    reference=c("reference_AF_afr",
                                "reference_AF_eas",
                                "reference_AF_eur",
                                "reference_AF_iam",
                                "reference_AF_sas"),
                    observed="gnomad_AF_afr",
                    pi.start = c(.2, .2, .2, .2, .2),
                    goodness.of.fit=TRUE) 
  
  expect_equal(length(testval), expected = 7)
})

