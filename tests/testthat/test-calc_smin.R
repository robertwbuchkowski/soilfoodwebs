test_that("multiplication works", {
  expect_equal(calc_smin(intro_comm), 4e-04)
  expect_type(calc_smin(intro_comm), "double")
})
