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
