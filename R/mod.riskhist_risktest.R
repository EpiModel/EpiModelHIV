
#' @title Risk History for Risk Based Testers Module
#'
#' @description Module function to track the risk history persons
#'              who test selectively based on recent risk.
#'
#' @inheritParams aging_msm
#'
#' @keywords module msm
#'
#' @export
#'
riskhist_risktest_msm <- function(dat, at) {


  ## Attributes
  uid <- dat$attr$uid
  dx <- dat$attr$diag.status

  ind.uai.nonmain <- dat$attr$ind.uai.nonmain 
  ind.uai.known.sd <- dat$attr$ind.uai.known.sd
  ind.newmain <- dat$attr$ind.newmain
  main.num <- dat$attr$main.num <- get_degree(dat$el[[1]])
  main.num.previous <- dat$attr$main.num.previous
  
  dat$attr$ind.uai.nonmain.clock <- ifelse(is.na(dat$attr$ind.uai.nonmain.clock) == FALSE, dat$attr$ind.uai.nonmain.clock + 1, dat$attr$ind.uai.nonmain.clock) 
  dat$attr$ind.uai.known.sd.clock <- ifelse(is.na(dat$attr$ind.uai.known.sd.clock) == FALSE, dat$attr$ind.uai.known.sd.clock + 1, dat$attr$ind.uai.known.sd.clock)
  dat$attr$ind.newmain.clock <- ifelse(is.na(dat$attr$ind.newmain.clock) == FALSE, dat$attr$ind.newmain.clock + 1, dat$attr$ind.newmain.clock)



  ## Parameters
  time.unit <- dat$param$time.unit

  ## Edgelist, adds uai summation per partnership from act list
  pid <- NULL # For R CMD Check
  al <- as.data.frame(dat$temp$al)
  by_pid <- group_by(al, pid)
  uai <- summarise(by_pid, uai = sum(uai))[, 2]
  el <- as.data.frame(cbind(dat$temp$el, uai))

  # Remove concordant positive edges
  #  el2 <- el[el$st2 == 0, ]
  # Keep discordant pairs for testing counts
  el2<-el
  
  ## Degree ##
  main.deg <- get_degree(dat$el[[1]])
  casl.deg <- get_degree(dat$el[[2]])
  inst.deg <- get_degree(dat$el[[3]])


  ## Condition 1: UAI in non-main partnerships
  uai.nmain <- unique(c(el2$p1[el2$st1 == 0 & el2$uai > 0 & el2$ptype %in% 2:3],
                        el2$p2[el2$uai > 0 & el2$ptype %in% 2:3]))
 
  ##DON"T UPDATE UNLESS NA
  uai.main.elig <- which(is.na(dat$attr$ind.uai.nonmain) == TRUE)
  uai.nmain <- intersect(uai.main.elig,uai.nmain)
  
  dat$attr$ind.uai.nonmain[uai.nmain] <- 1
  dat$attr$ind.uai.nonmain.clock[uai.nmain] <- 1
  
  ## Condition 2: AI within known serodiscordant partnerships
  el2.cond3 <- el2[el2$st1 == 1 & el2$ptype %in% 1:3 & el2$uai > 0, ]

  # Disclosure
  discl.list <- dat$temp$discl.list
  disclose.cdl <- discl.list[, 1] * 1e7 + discl.list[, 2]
  delt.cdl <- uid[el2.cond3[, 1]] * 1e7 + uid[el2.cond3[, 2]]
  discl <- (delt.cdl %in% disclose.cdl)
  ai.sd <- el2.cond3$p2[discl == TRUE]
 
  ##DON"T UPDATE UNLESS NA 
  uai.known.sd.elig <- which(is.na(dat$attr$ind.uai.known.sd) ==TRUE)
  ai.sd <- intersect(uai.known.sd.elig,ai.sd)
 
  dat$attr$ind.uai.known.sd[ai.sd] <- 1
  dat$attr$ind.uai.known.sd.clock[ai.sd] <- 1

  ## Condition 3: a new main partner
  idsNM <- which(main.deg > dat$attr$main.num.previous)

  ##DON"T UPDATE UNLESS NA
  NM.elig <- which(is.na(dat$attr$ind.newmain) ==TRUE)
  idsNM <- intersect(NM.elig,idsNM)
  
  dat$attr$ind.newmain[idsNM] <- 1
  dat$attr$ind.newmain.clock[idsNM] <- 1
  
  dat$attr$main.num.previous <- get_degree(dat$el[[1]])

  
  return(dat)
}
