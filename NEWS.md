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

