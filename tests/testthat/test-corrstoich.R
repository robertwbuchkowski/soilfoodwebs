test_that("Correct stoichiometry gives correct properties.", {
  expect_equal(corrstoich(intro_comm)$prop$ID, intro_comm$prop$ID)
  expect_equal(dim(corrstoich(intro_comm)$imat), dim(intro_comm$imat))
  expect_type(corrstoich(intro_comm), "list")
  expect_type(corrstoich(intro_comm)$imat, "double")
  expect_s3_class(corrstoich(intro_comm)$prop, "data.frame")
})
