
## Paper 1 estimation file

suppressMessages(library("EpiModelHIV"))
rm(list = ls())

load("est/nwstats.rda")


# 1. Main Model -----------------------------------------------------------

# Initialize network
nw.main <- base_nw.msm(st)

# Assign degree
nw.main <- assign_degree(nw.main, deg.type = "pers", nwstats = st)

# Formulas
formation.m <- ~edges +
                nodefactor("deg.pers") +
                absdiff("sqrt.age") +
                offset(nodematch("role.class", diff = TRUE, keep = 1:2))

# Fit model
fit.m <- netest(nw.main,
                formation = formation.m,
                coef.form = c(-Inf, -Inf),
                target.stats = st$stats.m,
                coef.diss = st$coef.diss.m,
                constraints = ~bd(maxout = 1),
                set.control.ergm = control.ergm(MPLE.max.dyad.types = 1e10,
                                                init.method = "zeros",
                                                MCMLE.maxit = 250))


# 2. Casual Model ---------------------------------------------------------

# Initialize network
nw.pers <- nw.main

# Assign degree
nw.pers <- assign_degree(nw.pers, deg.type = "main", nwstats = st)

# Formulas
formation.p <- ~edges +
                nodefactor("deg.main") +
                concurrent +
                absdiff("sqrt.age") +
                offset(nodematch("role.class", diff = TRUE, keep = 1:2))

# Fit model
fit.p <- netest(nw.pers,
                formation = formation.p,
                coef.form = c(-Inf, -Inf),
                target.stats = st$stats.p,
                coef.diss = st$coef.diss.p,
                constraints = ~bd(maxout = 2),
                set.control.ergm = control.ergm(MPLE.max.dyad.types = 1e9,
                                                init.method = "zeros",
                                                MCMLE.maxit = 250))


# Fit inst model ----------------------------------------------------------

# Initialize network
nw.inst <- nw.main

# Assign degree
nw.inst <- set.vertex.attribute(nw.inst, "deg.main", nw.pers %v% "deg.main")
nw.inst <- set.vertex.attribute(nw.inst, "deg.pers", nw.main %v% "deg.pers")
table(nw.inst %v% "deg.main", nw.inst %v% "deg.pers")

# Formulas
formation.i <- ~edges +
                nodefactor(c("deg.main", "deg.pers")) +
                nodefactor("riskg", base = 3) +
                absdiff("sqrt.age") +
                offset(nodematch("role.class", diff = TRUE, keep = 1:2))

# Fit model
fit.i <- netest(nw.inst,
                formation = formation.i,
                target.stats = st$stats.i,
                coef.form = c(-Inf, -Inf),
                coef.diss = dissolution_coefs(~offset(edges), 1),
                set.control.ergm = control.ergm(MPLE.max.dyad.types = 1e9,
                                                MCMLE.maxit = 250))

# Save data
est <- list(fit.m, fit.p, fit.i)
save(est, file = "est/fit.rda")


# Diagnostics -------------------------------------------------------------

load("est/fit.rda")
dx <- netdx(est[[3]], nsims = 10, dynamic = FALSE)
dx

dx <- simulate(est[[3]]$fit, statsonly = TRUE, nsim = 1000)
cbind(colMeans(dx)[1:11], est[[3]]$target.stats)
