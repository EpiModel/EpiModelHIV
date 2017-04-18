context("Aging Module")

data(est)
data(st)

test_that("Aging module", {

  param <- param_msm(nwstats = st)
  init <- init_msm(nwstats = st)
  control <- control_msm()

  at <- 1
  dat <- initialize_msm(est, param, init, control, s = 1)

  pre.age <- dat$attr$age
  expect_true(all(!is.na(pre.age)))
  expect_true(sum(pre.age) > 0)

  dat <- aging_msm(dat, at = 2)
  post.age <- dat$attr$age
  expect_true(all(post.age == pre.age + st$time.unit/365))

})
