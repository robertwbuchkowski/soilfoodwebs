test_that("Comana returns reasonable values", {
  expect_equal(attributes(comana(intro_comm))$names, c("consumption", "Cmin", "Nmin", "Nminmat","fmat", "Nfmat", "usin"))
  expect_type(comana(intro_comm)$consumption, "double")
  expect_type(comana(intro_comm)$Cmin, "double")
  expect_type(comana(intro_comm)$Nmin, "double")
  expect_type(comana(intro_comm)$Nminmat, "double")
  expect_type(comana(intro_comm)$fmat, "double")
  expect_type(comana(intro_comm)$usin, "list")
})
