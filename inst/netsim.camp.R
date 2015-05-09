
#' @export
netsim.camp <- function(x, param, init, control){

  dat <- initialize.camp(x, param, init, control, s = 1)

  for (at in max(2, control$start):control$nsteps) {

    # Age
    dat <- aging.camp(dat, at)

    # Deaths
    dat <- deaths.camp(dat, at)

    # Births
    dat <- births.camp(dat, at)

    # Testing
    dat <- test.camp(dat, at)

    # Treatment
    dat <- tx.camp(dat, at)

    # Progress
    dat <- progress.camp(dat, at)

    # Viral load
    dat <- update_vl.camp(dat, at)

    # Change inst.ai.class
    dat <- update_aiclass.camp(dat, at)

    # Change role.class
    dat <- update_roleclass.camp(dat, at)

    # Update edges terms
    dat <- edges_correct.camp(dat, at)

    # Relational dynamics
    dat <- resim_nets.camp(dat, at)

    # Disclosure
    dat <- disclose.camp(dat, at)

    # Discordant AI: Main / Pers / Inst
    dat <- acts.camp(dat, at)

    # Condom use: Main / Pers / Inst
    dat <- condoms.camp(dat, at)

    # Position: Main / Pers / Inst
    dat <- position.camp(dat, at)

    # Transmission
    dat <- trans.camp(dat, at)

    # Summary statistics
    dat <- getprev.camp(dat, at)

    # Console output
    verbose.camp(dat, type = "progress", s = 1, at)

  }

  return(dat)
}
