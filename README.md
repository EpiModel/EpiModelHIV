# EpiModelHIV
An R package for simulating HIV transmission dynamics among men who have sex with men and heterosexual populations, developed as an extension to our general network-based epidemic modeling platform, [EpiModel](http://epimodel.org).

`EpiModel` and `EpiModelHIV` use the statistical framework of temporal exponential-family random graph models to fit and simulate models of dynamic networks. These statistical methods have been developed by [Krivitsky and Handcock](http://onlinelibrary.wiley.com/doi/10.1111/rssb.12014/abstract) and implemented within the [Statnet](https://cran.r-project.org/web/packages/statnet/index.html) suite of statistical software. 

These packages combine these Statnet methods with an agent-based epidemic modeling engine to simulate HIV transmission over networks, allowing for complex dependencies between the network, epidemiological, and demographic changes in the simulated populations. Readers new to these methods are recommended to consult our [EpiModel](http://epimodel.org) resources, including our main [Vignette](http://statnet.github.io/tut/EpiModelVignette.pdf) describe the theory and implementation.

## Installation

You can install `EpiModelHIV` in R using `devtools`:
```
install.packages("EpiModel", dependencies = TRUE)
devtools::install_github("statnet/tergmLite")
devtools::install_github("statnet/EpiModelHIV")
```

Documentation on using this software package is forthcoming, although limited function documentation is provided within the package and available with the `help(package = "EpiModelHIV")` command.

### Note on Repository and Package Name
This Github repository `EpiModelHIV` and the R package within it were previously named `EpiModelHIVmsm`. On 6/24/2016, we merged that MSM package with our `EpiModelHIVhet` (HIV models for heterosexuals) into this combined repository and package. All scripts from those separate packages should still function with this current version after changing the input to `library`. 

## Authors
<table>
  <tr>
    <td><a href="http://samueljenness.org/" target="_blank">Samuel M. Jenness</a></th>
    <td>Department of Epidemiology</th>
    <td>Emory University</th>
  </tr>
  <tr>
    <td><a href="http://faculty.washington.edu/goodreau/" target="_blank">Steven M. Goodreau</a></td>
    <td>Department of Anthropology</td>
    <td>University of Washington</td>
  </tr>
</table>

## Literature

`EpiModelHIV` has been used in the following scientific articles published or in press:

1. Jenness SM, Goodreau SM, Morris M, Cassels S. Effectiveness of Combination Packages for HIV-1 Prevention in Sub-Saharan Africa Depends on Partnership Network Structure. Sexually Transmitted Infections. Epub ahead of print. DOI: 10.1136/sextrans-2015-052476. [LINK](http://sti.bmj.com/content/early/2016/06/09/sextrans-2015-052476.abstract)

2. Jenness SM, Goodreau SM, Rosenberg E, Beylerian EN, Hoover KW, Smith DK, Sullivan P. Impact of CDC’s HIV Preexposure Prophylaxis Guidelines among MSM in the United States. In Press, Journal of Infectious Diseases.

