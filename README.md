EpiModelHIV
===============

[![Build Status](https://travis-ci.org/statnet/EpiModelHIV.svg?branch=master)](https://travis-ci.org/statnet/EpiModelHIV)

An R package for simulating HIV transmission dynamics among men who have sex with men and heterosexual populations, developed as an extension to our general network-based epidemic modeling platform, [EpiModel](http://epimodel.org).

`EpiModel` and `EpiModelHIV` use the statistical framework of temporal exponential-family random graph models to fit and simulate models of dynamic networks. These [statistical methods](http://onlinelibrary.wiley.com/doi/10.1111/rssb.12014/abstract) have been developed and implemented as open-source software, building on the extensive efforts of the [Statnet](https://statnet.org/) research group to build software tools for the representation, analysis, and visualization of complex network data.

These packages combine these Statnet methods with an agent-based epidemic modeling engine to simulate HIV transmission over networks, allowing for complex dependencies between the network, epidemiological, and demographic changes in the simulated populations. Readers new to these methods are recommended to consult our [EpiModel](http://epimodel.org) resources, including our main [Vignette](http://statnet.github.io/tut/EpiModelVignette.pdf) describe the theory and implementation.

## Installation

You can install `EpiModelHIV` in R using `devtools`:
```
install.packages("EpiModel", dependencies = TRUE)
devtools::install_github("statnet/tergmLite")
devtools::install_github("statnet/EpiModelHPC")
devtools::install_github("statnet/EpiModelHIV")
```

Documentation on using this software package is forthcoming, although limited function documentation is provided within the package and available with the `help(package = "EpiModelHIV")` command.

## EpiModelHIV Software Development Team

Author                                                        | Department                 | Institution
------------------------------------------------------------- | -------------------------- | ------------------------
[Samuel M. Jenness](http://samueljenness.org/)                | Department of Epidemiology | Emory University
[Steven M. Goodreau](http://faculty.washington.edu/goodreau/) | Department of Anthropology | University of Washington
Deven T. Hamilton                                             | CSDE                       | University of Washington
Kevin M. Weiss                                                | Department of Epidemiology | Emory University


## Literature

`EpiModelHIV` has been used in the following scientific articles:

1. Jenness SM, Goodreau SM, Morris M, Cassels S. Effectiveness of Combination Packages for HIV-1 Prevention in Sub-Saharan Africa Depends on Partnership Network Structure. _Sexually Transmitted Infections._ 2016; 92(8): 619-624. [LINK](http://sti.bmj.com/content/early/2016/06/09/sextrans-2015-052476.abstract)

2. Jenness SM, Goodreau SM, Rosenberg E, Beylerian EN, Hoover KW, Smith DK, Sullivan P. Impact of CDC’s HIV Preexposure Prophylaxis Guidelines among MSM in the United States. _Journal of Infectious Diseases._ 2016; 214(12): 1800-1807. [LINK](http://jid.oxfordjournals.org/content/early/2016/07/12/infdis.jiw223.full)

3. Jenness SM, Sharma A, Goodreau SM, Rosenberg ES, Weiss KM, Hoover KW, Smith DK, Sullivan P. Individual HIV Risk versus Population Impact of Risk Compensation after HIV Preexposure Prophylaxis Initiation among Men Who Have Sex with Men. _PLoS One._ 2017; 12(1): e0169484. [LINK](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0169484)

4. Goodreau SM, Rosenberg ES, Jenness SM, Luisi N, Stansfield SE, Millett G, Sullivan P. Sources of Racial Disparities in HIV Prevalence among Men Who Have Sex with Men in Atlanta, GA: A Modeling Study. _Lancet HIV._ Epub ahead of print. DOI: 10.1016/ S2352-3018(17)30067. [LINK](https://www.ncbi.nlm.nih.gov/pubmed/28431923)

5. Jenness SM, Weiss KM, Goodreau SM, Rosenberg E, Gift T, Chesson H, Hoover KW, Smith DK, Liu AY, Sullivan P. Incidence of Gonorrhea and Chlamydia Following HIV Preexposure Prophylaxis among Men Who Have Sex with Men: A Modeling Study. _Clinical Infectious Diseases._ In Press.
