#Cohab

load("F:/GitHub/EpiModelHIV/data/ego.obj_c.rda")

## Sample weights
numegos <- length(ego.obj_c$egoWt)
sumwt <- sum(ego.obj_c$egoWt)
wtscale <- numegos/sumwt
ego.obj_c$egoWt <- ego.obj_c$egoWt*wtscale
ego.obj_c$egos$weight <- NULL
save(ego.obj_c, file='data/sampwt.ego.obj_c')

## Unwtd
ego.obj_c$egoWt <- rep(1, length(ego.obj_c$egoWt))
ego.obj_c$egos$weight <- NULL
save(ego.obj_c, file='data/unwtd.ego.obj_c')

# Persistent
load("F:/GitHub/EpiModelHIV/data/ego.obj_p.rda")

## Sample weights
numegos <- length(ego.obj_p$egoWt)
sumwt <- sum(ego.obj_p$egoWt)
wtscale <- numegos/sumwt
ego.obj_p$egoWt <- ego.obj_p$egoWt*wtscale
ego.obj_p$egos$weight <- NULL
save(ego.obj_p, file='data/sampwt.ego.obj_p')

## Unwtd
ego.obj_p$egoWt <- rep(1, length(ego.obj_p$egoWt))
ego.obj_p$egos$weight <- NULL
save(ego.obj_p, file='data/unwtd.ego.obj_p')

# One time

load("F:/GitHub/EpiModelHIV/data/ego.obj_i.rda")

## Sample weights
numegos <- length(ego.obj_i$egoWt)
sumwt <- sum(ego.obj_i$egoWt)
wtscale <- numegos/sumwt
ego.obj_i$egoWt <- ego.obj_i$egoWt*wtscale
ego.obj_i$egos$weight <- NULL
save(ego.obj_i, file='data/sampwt.ego.obj_i')

## Unwtd
ego.obj_i$egoWt <- rep(1, length(ego.obj_i$egoWt))
ego.obj_i$egos$weight <- NULL
save(ego.obj_i, file='data/unwtd.ego.obj_i')




