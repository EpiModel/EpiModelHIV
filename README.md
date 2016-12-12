EpiModelHIV
===============

[![Build Status](https://travis-ci.org/statnet/EpiModelHIV.svg?branch=master)](https://travis-ci.org/statnet/EpiModelHIV)

An R package for simulating HIV transmission dynamics among men who have sex with men and heterosexual populations, developed as an extension to our general network-based epidemic modeling platform, [EpiModel](http://epimodel.org).

`EpiModel` and `EpiModelHIV` use the statistical framework of temporal exponential-family random graph models to fit and simulate models of dynamic networks. These [statistical methods](http://onlinelibrary.wiley.com/doi/10.1111/rssb.12014/abstract) have been developed and implemented as open-source software, building on the extensive efforts of the [Statnet](https://statnet.org/) research group to build software tools for the representation, analysis, and visualization of complex network data.

These packages combine these Statnet methods with an agent-based epidemic modeling engine to simulate HIV transmission over networks, allowing for complex dependencies between the network, epidemiological, and demographic changes in the simulated populations. Readers new to these methods are recommended to consult our [EpiModel](http://epimodel.org) resources, including our main [Vignette](http://statnet.github.io/tut/EpiModelVignette.pdf) describe the theory and implementation.

## Installation

You can install `EpiModelHIV` in R using `devtools`:
```R
install.packages("EpiModel", dependencies = TRUE)
devtools::install_github("statnet/tergmLite", subdir = "tergmLite")
devtools::install_github("statnet/EpiModelHPC", build_vignettes = TRUE)
devtools::install_github("statnet/EpiModelHIV")
```

Documentation on using this software package is forthcoming, although limited function documentation is provided within the package and available with the `help(package = "EpiModelHIV")` command.

### Note on Repository and Package Name
This Github repository `EpiModelHIV` and the R package within it were previously named `EpiModelHIVmsm`. On 6/24/2016, we merged that MSM package with our `EpiModelHIVhet` (HIV models for heterosexuals) into this combined repository and package. All scripts from those separate packages should still function with this current version after changing the input to `library`. 

## Software Development Team

Author                                                        | Department                 | Institution
------------------------------------------------------------- | -------------------------- | ------------------------
[Samuel M. Jenness](http://samueljenness.org/)                | Department of Epidemiology | Emory University
[Steven M. Goodreau](http://faculty.washington.edu/goodreau/) | Department of Anthropology | University of Washington


## Literature

`EpiModelHIV` has been used in the following scientific articles:

1. Jenness SM, Goodreau SM, Morris M, Cassels S. Effectiveness of Combination Packages for HIV-1 Prevention in Sub-Saharan Africa Depends on Partnership Network Structure. _Sexually Transmitted Infections._ [DOI: 10.1136/sextrans-2015-052476.](http://sti.bmj.com/content/early/2016/06/09/sextrans-2015-052476.abstract)

2. Jenness SM, Goodreau SM, Rosenberg E, Beylerian EN, Hoover KW, Smith DK, Sullivan P. Impact of CDC’s HIV Preexposure Prophylaxis Guidelines among MSM in the United States. _Journal of Infectious Diseases._ [DOI: 10.1093/infdis/jiw223.](http://jid.oxfordjournals.org/content/early/2016/07/12/infdis.jiw223.full)

