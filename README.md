## About this branch

This branch is a frozen snapshot of the EpiModelHIV engine used in these two papers:

> Jones J, Le Guillou A, Gift TL, Chesson H, Bernstein KT, Delaney KP, Lyles C, Berruti A, Sullivan PS, Jenness SM. Effect of Screening and Treatment for Gonorrhea and Chlamydia on HIV Incidence Among Men Who Have Sex With Men in the United States: A Modeling Analysis. _Sex Transm Dis._ 2022; 49(10): 669-676. https://doi.org/10.1097/OLQ.0000000000001685
>
> Jones J, Jenness SM, Le Guillou A, Sullivan PS, Gift TL, Delaney KP, Chesson H. Estimated Number of Incident HIV Infections in Men Who Have Sex With Men Attributable to Gonorrhea and Chlamydia, Per Gonococcal or Chlamydial Infection, in the United States. _Sex Transm Dis._ 2023; 50(2): 83-85. https://doi.org/10.1097/OLQ.0000000000001724

- Analysis code: https://github.com/EpiModel/PAFGuidelines
- Release: https://github.com/EpiModel/EpiModelHIV/releases/tag/GCCT-STD
- Install: `remotes::install_github("EpiModel/EpiModelHIV@GCCT-STD")`

It is not maintained. It is a copy of EpiModelHIV-p tag `PAFGuidelines` at `e96a668` (May 17, 2022). The private ARTnetData package is listed only under Suggests, so the package installs without it.

EpiModelHIV
===============

[![Build Status](https://travis-ci.org/statnet/EpiModelHIV.svg?branch=master)](https://travis-ci.org/statnet/EpiModelHIV)

Modules for simulating HIV/STI transmission dynamics among men who have sex with men and heterosexual populations, developed as an extension to our general network-based epidemic modeling platform, [EpiModel](http://epimodel.org).

`EpiModel` and `EpiModelHIV` use the statistical framework of temporal exponential-family random graph models to fit and simulate models of dynamic networks. These [statistical methods](http://onlinelibrary.wiley.com/doi/10.1111/rssb.12014/abstract) have been developed and implemented as open-source software, building on the extensive efforts of the [Statnet](https://statnet.org/) research group to build software tools for the representation, analysis, and visualization of complex network data.

These packages combine these Statnet methods with an agent-based epidemic modeling engine to simulate HIV transmission over networks, allowing for complex dependencies between the network, epidemiological, and demographic changes in the simulated populations. Readers new to these methods are recommended to consult our [EpiModel](http://epimodel.org) resources, including our main methods paper [Vignette](http://doi.org/10.18637/jss.v084.i08) describing the theory and implementation.

## Installation

You can install `EpiModelHIV` in R using `remotes`:
```
install.packages("EpiModel", dependencies = TRUE)
remotes::install_github("statnet/tergmLite")
remotes::install_github("statnet/EpiModelHPC")
remotes::install_github("statnet/EpiModelHIV")
```

Documentation on using this software package is forthcoming, although limited function documentation is provided within the package and available with the `help(package = "EpiModelHIV")` command.
