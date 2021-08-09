test_that("Correct dimensions", {
  expect_equal(dim(whomineralizes(intro_comm)), c(7,5))
  expect_equal(whomineralizes(intro_comm)$ID, intro_comm$prop$ID)
})
