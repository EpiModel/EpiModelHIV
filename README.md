## About this branch

This branch is a frozen snapshot of the EpiModelHIV engine used in:

> Chandra C, Marcus JL, Jones J, Crawford ND, Galvin HK, Mayer KH, Sinha P, Le Guillou A, Jenness SM. PrEP reinitiation, uptake, and persistence as targets for reducing HIV incidence and racial and ethnic inequities among men who have sex with men in Boston, USA: a network modeling study. _Lancet Reg Health Am._ 2026; 63: 101620. https://doi.org/10.1016/j.lana.2026.101620

- Analysis code: https://github.com/EpiModel/PrEP-Reinitiation
- Release: https://github.com/EpiModel/EpiModelHIV/releases/tag/PrEPReinit-LRHA
- Install: `remotes::install_github("EpiModel/EpiModelHIV@PrEPReinit-LRHA")`

It is not maintained. It is a copy of EpiModelHIV-p `PrEP_reinitiation` at `5c63f1f` (May 12, 2025). The private ARTnetData package is listed only under Suggests, so the package installs without it.

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
