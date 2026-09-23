test_that("run_app is available", {
  expect_true(is.function(run_app))
})

test_that("stacked PEAKS PTMs both convert", {
  expect_equal(
    normalize_peptidoform("M(+42.01)(+15.99)NSLSEANTKF"),
    "M(Acetyl)(Oxidation)NSLSEANTKF"   # adjust to your real expected output
  )
})
