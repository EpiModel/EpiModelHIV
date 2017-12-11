
# Update stiPrEP stack

# devtools::install_github("statnet/EpiModel")
# devtools::install_github("statnet/EpiModelHPC")
# devtools::install_github("statnet/tergmLite", subdir = "tergmLite")
devtools::install_github("statnet/EpiModelHIV", ref = "recal-sti")


## interface with hyak

# upload scripts
system("scp scripts/burnin/*.burn.[Rs]* hyak:/gscratch/csde/sjenness/stia")

# upload inputs
# system("scp est/*.rda hyak:/gscratch/csde/sjenness/stia/est")
