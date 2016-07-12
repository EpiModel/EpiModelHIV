context("Numeric attributes for C++ code")

test_that("Numeric for treatment trajectory", {
    expect_equal(sum(sim$attr$tt.traj.W.prob), 1)
    expect_error(sum(sim$attr$tt.traj.B.prob) != 1, "Proportion of men in treatment categories must sum to 1")
})

test_that("Numeric for disease stage", {
    expect_true(all(!(is.na(sim$attr$inf.stage)) %in% c(1, 2, 3, 4)))
})