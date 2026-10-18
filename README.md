## About this branch

This branch is a frozen snapshot of the EpiModelHIV engine used in:

> Jenness SM, Le Guillou A, Lyles C, Bernstein KT, Krupinsky K, Enns EA, Sullivan PS, Delaney KP. The Role of HIV Partner Services in the Modern Biomedical HIV Prevention Era: A Network Modeling Study. _Sex Transm Dis._ 2022; 49(12): 801-807. https://doi.org/10.1097/OLQ.0000000000001711

- Analysis code: https://github.com/EpiModel/CombPrevNet
- Release: https://github.com/EpiModel/EpiModelHIV/releases/tag/PartServ-STD
- Install: `remotes::install_github("EpiModel/EpiModelHIV@PartServ-STD")`

It is not maintained. It is a copy of EpiModelHIV-p tag `CombPrevNet` at `d3dda21` (June 16, 2021). The private ARTnetData package is listed only under Suggests, so the package installs without it.

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
