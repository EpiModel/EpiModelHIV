## About this branch

This branch is a frozen snapshot of the EpiModelHIV engine used in:

> Chandra C, Maloney KM, Le Guillou A, Hoover KW, Jenness SM. Modeling the Potential Impact of Scaling Up Event-Driven PrEP Among Gay, Bisexual, and Other Men Who Have Sex With Men. _J Acquir Immune Defic Syndr._ 2026; 101(11): 1219-1230. https://doi.org/10.1097/QAI.0000000000003934

- Analysis code: https://github.com/EpiModel/EventDrivenPrEP
- Release: https://github.com/EpiModel/EpiModelHIV/releases/tag/EDP-JAIDS
- Install: `remotes::install_github("EpiModel/EpiModelHIV@EDP-JAIDS")`

It is not maintained. Its engine code is a copy of EpiModelHIV-p `EDP` at `096e28d` (December 1, 2023), the version pinned in the project's `renv.lock`. The private ARTnetData package is listed only under Suggests, so the package installs without it.

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
