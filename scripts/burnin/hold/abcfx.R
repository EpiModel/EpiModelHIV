
## ABC fx

simParmSamp <- function(n = NULL, parms = NULL) {

  if (is.null(parms)) {
    rgc.tprob <- runif(n, 0.35, 0.60)
    ugc.tprob <- runif(n, 0.20, 0.40)
    rct.tprob <- runif(n, 0.35, 0.60)
    uct.tprob <- runif(n, 0.20, 0.40)

    rgc.dur.asympt <- runif(n, 26, 52)
    ugc.dur.asympt <- runif(n, 26, 52)
    rct.dur.asympt <- runif(n, 39, 65)
    uct.dur.asympt <- runif(n, 39, 65)

    rgc.sympt.prob <- runif(n, 0.05, 0.20)
    ugc.sympt.prob <- runif(n, 0.60, 0.95)
    rct.sympt.prob <- runif(n, 0.05, 0.20)
    uct.sympt.prob <- runif(n, 0.60, 0.95)

    hiv.rect.rr <- runif(n, 2, 3)
    hiv.ureth.rr <- runif(n, 1, 2)

    df <- data.frame(rgc.tprob, ugc.tprob, rct.tprob, uct.tprob,
                     rgc.dur.asympt, ugc.dur.asympt, rct.dur.asympt, uct.dur.asympt,
                     rgc.sympt.prob, ugc.sympt.prob, rct.sympt.prob, uct.sympt.prob,
                     hiv.rect.rr, hiv.ureth.rr)
  } else {
    df <- parms
    df$rgc.tprob <- dunif(parms$rgc.tprob, 0.35, 0.60)
    df$ugc.tprob <- dunif(parms$ugc.tprob, 0.20, 0.40)
    df$rct.tprob <- dunif(parms$rct.tprob, 0.35, 0.60)
    df$uct.tprob <- dunif(parms$uct.tprob, 0.20, 0.40)

    df$rgc.dur.asympt <- dunif(parms$rgc.dur.asympt, 26, 52)
    df$ugc.dur.asympt <- dunif(parms$ugc.dur.asympt, 26, 52)
    df$rct.dur.asympt <- dunif(parms$rct.dur.asympt, 39, 65)
    df$uct.dur.asympt <- dunif(parms$uct.dur.asympt, 39, 65)

    df$rgc.sympt.prob <- dunif(parms$rgc.sympt.prob, 0.05, 0.20)
    df$ugc.sympt.prob <- dunif(parms$ugc.sympt.prob, 0.60, 0.95)
    df$rct.sympt.prob <- dunif(parms$rct.sympt.prob, 0.05, 0.20)
    df$uct.sympt.prob <- dunif(parms$uct.sympt.prob, 0.60, 0.95)

    df$hiv.rect.rr <- dunif(parms$hiv.rect.rr, 2, 3)
    df$hiv.ureth.rr <- dunif(parms$hiv.ureth.rr, 1, 2)
  }
  return(df)
}

perturbParticle <- function(parms,
                            pnames,
                            from = NULL,
                            sds) {

  old.parms <- new.parms <- parms[, pnames, drop = FALSE]

  if (is.null(from)) {

    for(pp in 1:ncol(old.parms)) {
      new.parms[, pp] <- old.parms[, pp, drop = FALSE] + runif(nrow(old.parms), -sds[pp], sds[pp])
    }

    PriorDensities <- simParmSamp(parms = new.parms)
    PriorDensity <- apply(PriorDensities, 1, prod)
    new.parms <- new.parms[PriorDensity > 0, , drop = FALSE]

    return(new.parms)

  } else {

    Lparms <- as.matrix(parms)
    dprior <- Lfrom <- as.matrix(from)

    for(pp in 1:length(Lparms)) {
      dprior[,pp] <- dunif(Lparms[,pp], Lfrom[,pp] - sds[pp], Lfrom[,pp] + sds[pp])
    }

    dpriorProd <- as.numeric(apply(dprior, 1, prod))

    return(dpriorProd)
  }
}

weightParticles <- function(pnames, currentBatch, lastBatch, sdsUse) {

  dpriors <- simParmSamp(parms = currentBatch)[, pnames, drop = FALSE]
  numerator <- apply(dpriors, 1, prod)
  denominator <- rep(NA, length(numerator))

  ## for each particle, calculate the probability it could have been gotten to
  ## from all previous particles, weighted by their weights
  for(jj in 1:length(numerator)) {
    Ks <- perturbParticle(parms = currentBatch[jj, pnames, drop = FALSE],
                          from = lastBatch[, pnames, drop = FALSE],
                          sds = sdsUse)
    denominator[jj] <- sum(lastBatch$weight * Ks)
    if (denominator[jj] == 0) browser()
  }
  return(numerator/denominator)
}

pal <- RColorBrewer::brewer.pal(3, "Set1")
nhist <- function(...) hist(..., col = pal[3], border = "white")

plotStats <- function(sim, targets) {
  par(mfrow = c(2,3), mar = c(3,3,3,1), mgp = c(2,1,0))
  stats <- 15:19
  for (i in seq_along(stats)) {
    nhist(sim[[stats[i]]], main = names(sim[stats[i]]))
    abline(v = mean(sim[[stats[i]]]), col = "blue", lwd = 2)
    abline(v = median(sim[[stats[i]]]), col = "blue", lty = 2, lwd = 2)
    abline(v = targets[i], col = "red", lwd = 2)
  }
}

plotParam <- function(sim) {
  par(mfrow = c(4,4))
  for (i in 1:14) {
    nhist(sim[[i]], main = names(sim[i]))
  }
}

rejection <- function(sim, targets, threshold) {
  edist <- sapply(1:nrow(sim), function(x) sqrt(sum((targets - sim[x, 15:19])^2)))
  edist.quant <- quantile(edist, threshold)
  accepted <- which(edist <= edist.quant)
  post <- sim[accepted, ]
  cat("\n Accepted n:", nrow(post))
  return(post)
}

abc_pre <- function(batch, batch.size) {
  fn <- paste0("scripts/burnin/simChosen.b", batch-1, ".rda")
  load(fn)

  init.samp.size = 5e5
  lastParticleSamp <- simChosen[sample(x = 1:nrow(simChosen),
                                       size = init.samp.size,
                                       prob = simChosen$weight,
                                       replace = TRUE), ]

  # restrict new sample to initial priors
  simParms <- perturbParticle(parms = lastParticleSamp,
                              pnames = names(simChosen[1:14]),
                              from = NULL,
                              sds = sds[1:14])
  simParms <- simParms[sample(1:nrow(simParms), size = batch.size), ]
  row.names(simParms) <- NULL

  cat("\n simParms dimensions: \n"); print(dim(simParms))
  cat("\n\n simParms lengths: \n"); print(sapply(simParms, function(x) length(unique(x))))
  cat("\n\n simParms head: \n"); print(head(simParms))
  return(simParms)
}

abc_post <- function(simChosen, batch) {

  sdsNew <- apply(simChosen, 2, sd)

  simChosenNew <- simChosen

  old.fn <- paste0("scripts/burnin/simChosen.b", batch-1, ".rda")
  load(old.fn)
  simChosenOld <- simChosen
  sdsOld <- sds

  simChosenNew$weight <- weightParticles(pnames = names(simChosen[1:14]),
                                         currentBatch = simChosenNew,
                                         lastBatch = simChosenOld,
                                         sdsUse = sdsOld[1:14])
  simChosenNew$weight[simChosenNew$weight == Inf] <- 0
  simChosenNew$weight <- simChosenNew$weight/sum(simChosenNew$weight)

  simChosen <- simChosenNew
  sds <- sdsNew

  out <- list(simChosen, sds)
  return(out)
}
