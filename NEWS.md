## EpiModelHIV 3.2.0

### NEW FEATURES

- The Syphilis processes now are closer to the Tuite et al 2020 paper

### OTHER

- Refactor of the PrEP module. Update the indication to the 2021 guidelines and enforce consistent naming conventions
- The "Partner Services" modules is now disable by default. See "R/mod.partident.R" for information on how to enable it.
- The modules have been re-ordered to simplify the handling of competitiveness between processes.
- The STI modules have been deprecated in favor of a single `sti.FUN` module that calls the Gonorrhea, Chlamydia and Syphilis process functions.
- The attributes and parameters names for all STIs have been standardized

## EpiModelHIV 3.1.0

### NEW FEATURES

- default parameters of the model as well as there definition and sources are
now defined in "inst/model_parameters.xlsx"
- ability to overwrite some node attributes with a `data.frame`. See
`?EpiModelHIV:::overwrite_attributes`

### DEPRECATION

- `param_msm` has been deprecated in favor of the more generic `EpiModel::param.net`.

### OTHER

- `init_msm` default values for STIs prevalence were removed. In order to make
no STI models a conscious choice.
