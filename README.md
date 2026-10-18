## About this branch

This branch is a frozen snapshot of the EpiModelHIV engine used in:

> Wheatley MM, Knowlton G, Kao SY, Jenness SM, Enns EA. Cost-Effectiveness of Interventions to Improve HIV Pre-exposure Prophylaxis Initiation, Adherence, and Persistence Among Men Who Have Sex With Men. _J Acquir Immune Defic Syndr._ 2022; 90(1): 41-49. https://doi.org/10.1097/QAI.0000000000002921

- Analysis code: https://github.com/EpiModel/PrEP-CEA
- Release: https://github.com/EpiModel/EpiModelHIV/releases/tag/PrEPCEA-JAIDS
- Install: `remotes::install_github("EpiModel/EpiModelHIV@PrEPCEA-JAIDS")`

It is not maintained. It is a copy of the EpiModelHIV-p `PrEP-CEA` branch at `92c88cb` (July 20, 2021), which is kept in EpiModelHIV-p as tag `PrEP-CEA-PSA`. It has no private package dependencies.

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
