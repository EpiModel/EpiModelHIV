## About this branch

This branch is a frozen snapshot of the EpiModelHIV engine used in:

> Onwubiko UN, Le Guillou A, Chamberlain AT, Benkeser D, Holland DP, Baral SD, Jenness SM. Evaluating mental health integration in PrEP programs to reduce HIV incidence among US gay and bisexual men: a network modeling study. _BMC Health Serv Res._ 2026; 26(1): 1381. https://doi.org/10.1186/s12913-026-15074-0

- Analysis code: https://github.com/EpiModel/StigmaDepressionHIV
- Release: https://github.com/EpiModel/EpiModelHIV/releases/tag/MentalHealth-BMCHSR
- Install: `remotes::install_github("EpiModel/EpiModelHIV@MentalHealth-BMCHSR")`

It is not maintained. It is a copy of EpiModelHIV-p `StigmaDepressionHIV_real` at `85b8284` (May 15, 2024). The private ARTnetData package is listed only under Suggests, so the package installs without it.

EpiModelHIV
===============

Modules for simulating HIV/STI transmission dynamics among men who have sex with men and heterosexual populations, developed as an extension to our general network-based epidemic modeling platform, [EpiModel](http://www.epimodel.org).

`EpiModel` and `EpiModelHIV` use the statistical framework of temporal exponential-family random graph models to fit and simulate models of dynamic networks. These [statistical methods](https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/rssb.12014) have been developed and implemented as open-source software, building on the extensive efforts of the [Statnet](https://statnet.org/) research group to build software tools for the representation, analysis, and visualization of complex network data.

These packages combine these Statnet methods with an agent-based epidemic modeling engine to simulate HIV transmission over networks, allowing for complex dependencies between the network, epidemiological, and demographic changes in the simulated populations. Readers new to these methods are recommended to consult our [EpiModel](http://www.epimodel.org) resources, including our main methods paper [Vignette](https://doi.org/10.18637/jss.v084.i08) describing the theory and implementation.

## Installation

You can install `EpiModelHIV` and package dependencies in R using `pak` in the project root directory:
```r
install.packages("pak")
pak::pak()
```

Note that `EpiModelHIV` depends on the `ARTnetData` package, which is private. See details about access [here](https://github.com/EpiModel/ARTnet).
