
#' @title Disclosure Module for up to 5 race groups and heterosexuals and MSM.
#'
#' @description Module function for disclosure of HIV status to partners given
#'              non-disclosure in the past.
#'
#' @inheritParams aging_msm
#'
#' @details
#' Persons who are infected may disclose their status to partners at three
#' distinct time points: at relationship onset for newly formed discordant
#' pairs; at diagnosis for pairs starting as both negative but with one newly
#' infected; or post diagnosis for one recently infected. The rates of disclosure
#' vary at these three points, and also by the partnership type race and sexual identity.
#'
#' @return
#' This function returns the \code{dat} object with the updated disclosure list,
#' on \code{temp$discl.list}.
#'
#' @keywords module SHAMP msm het
#' @export
#'
disclose_shamp <- function(dat, at){

  for (type in c("main", "pers", "inst")) {

    # Variables --------------------------------------------------------------

    # Attributes
    status <- dat$attr$status
    uid <- dat$attr$uid
    diag.status <- dat$attr$diag.status
    diag.time <- dat$attr$diag.time
    race <- dat$attr$race
    sex <- dat$attr$sex
    sex.ident <- dat$attr$sex.ident
    immig.loc <- dat$attr$immig.loc

    # Parameters and network
    if (type == "main") {
      #Female
      disc.outset.B.f.prob <- dat$param$disc.outset.main.B.f.prob
      disc.at.diag.B.f.prob <- dat$param$disc.at.diag.main.B.f.prob
      disc.post.diag.B.f.prob <- dat$param$disc.post.diag.main.B.f.prob
      
      disc.outset.BI.f.prob <- dat$param$disc.outset.main.BI.f.prob
      disc.at.diag.BI.f.prob <- dat$param$disc.at.diag.main.BI.f.prob
      disc.post.diag.BI.f.prob <- dat$param$disc.post.diag.main.BI.f.prob
      
      disc.outset.H.f.prob <- dat$param$disc.outset.main.H.f.prob
      disc.at.diag.H.f.prob <- dat$param$disc.at.diag.main.H.f.prob
      disc.post.diag.H.f.prob <- dat$param$disc.post.diag.main.H.f.prob
      
      disc.outset.HI.f.prob <- dat$param$disc.outset.main.HI.f.prob
      disc.at.diag.HI.f.prob <- dat$param$disc.at.diag.main.HI.f.prob
      disc.post.diag.HI.f.prob <- dat$param$disc.post.diag.main.HI.f.prob
      
      disc.outset.W.f.prob <- dat$param$disc.outset.main.W.f.prob
      disc.at.diag.W.f.prob <- dat$param$disc.at.diag.main.W.f.prob
      disc.post.diag.W.f.prob <- dat$param$disc.post.diag.main.W.f.prob
      
      #Male het msf
      disc.outset.B.msf.prob <- dat$param$disc.outset.main.B.msf.prob
      disc.at.diag.B.msf.prob <- dat$param$disc.at.diag.main.B.msf.prob
      disc.post.diag.B.msf.prob <- dat$param$disc.post.diag.main.B.msf.prob
      
      disc.outset.BI.msf.prob <- dat$param$disc.outset.main.BI.msf.prob
      disc.at.diag.BI.msf.prob <- dat$param$disc.at.diag.main.BI.msf.prob
      disc.post.diag.BI.msf.prob <- dat$param$disc.post.diag.main.BI.msf.prob

      disc.outset.H.msf.prob <- dat$param$disc.outset.main.H.msf.prob
      disc.at.diag.H.msf.prob <- dat$param$disc.at.diag.main.H.msf.prob
      disc.post.diag.H.msf.prob <- dat$param$disc.post.diag.main.H.msf.prob
      
      disc.outset.HI.msf.prob <- dat$param$disc.outset.main.HI.msf.prob
      disc.at.diag.HI.msf.prob <- dat$param$disc.at.diag.main.HI.msf.prob
      disc.post.diag.HI.msf.prob <- dat$param$disc.post.diag.main.HI.msf.prob
      
      disc.outset.W.msf.prob <- dat$param$disc.outset.main.W.msf.prob
      disc.at.diag.W.msf.prob <- dat$param$disc.at.diag.main.W.msf.prob
      disc.post.diag.W.msf.prob <- dat$param$disc.post.diag.main.W.msf.prob
      
      ## Males MSM
      disc.outset.B.msm.prob <- dat$param$disc.outset.main.B.msm.prob
      disc.at.diag.B.msm.prob <- dat$param$disc.at.diag.main.B.msm.prob
      disc.post.diag.B.msm.prob <- dat$param$disc.post.diag.main.B.msm.prob
      
      disc.outset.BI.msm.prob <- dat$param$disc.outset.main.BI.msm.prob
      disc.at.diag.BI.msm.prob <- dat$param$disc.at.diag.main.BI.msm.prob
      disc.post.diag.BI.msm.prob <- dat$param$disc.post.diag.main.BI.msm.prob
      
      disc.outset.H.msm.prob <- dat$param$disc.outset.main.H.msm.prob
      disc.at.diag.H.msm.prob <- dat$param$disc.at.diag.main.H.msm.prob
      disc.post.diag.H.msm.prob <- dat$param$disc.post.diag.main.H.msm.prob
      
      disc.outset.HI.msm.prob <- dat$param$disc.outset.main.HI.msm.prob
      disc.at.diag.HI.msm.prob <- dat$param$disc.at.diag.main.HI.msm.prob
      disc.post.diag.HI.msm.prob <- dat$param$disc.post.diag.main.HI.msm.prob
      
      disc.outset.W.msm.prob <- dat$param$disc.outset.main.W.msm.prob
      disc.at.diag.W.msm.prob <- dat$param$disc.at.diag.main.W.msm.prob
      disc.post.diag.W.msm.prob <- dat$param$disc.post.diag.main.W.msm.prob
      
      disc.outset.B.msmf.prob <- dat$param$disc.outset.main.B.msmf.prob
      disc.at.diag.B.msmf.prob <- dat$param$disc.at.diag.main.B.msmf.prob
      disc.post.diag.B.msmf.prob <- dat$param$disc.post.diag.main.B.msmf.prob
      
      disc.outset.BI.msmf.prob <- dat$param$disc.outset.main.BI.msmf.prob
      disc.at.diag.BI.msmf.prob <- dat$param$disc.at.diag.main.BI.msmf.prob
      disc.post.diag.BI.msmf.prob <- dat$param$disc.post.diag.main.BI.msmf.prob
      
      disc.outset.H.msmf.prob <- dat$param$disc.outset.main.H.msmf.prob
      disc.at.diag.H.msmf.prob <- dat$param$disc.at.diag.main.H.msmf.prob
      disc.post.diag.H.msmf.prob <- dat$param$disc.post.diag.main.H.msmf.prob
      
      disc.outset.HI.msmf.prob <- dat$param$disc.outset.main.HI.msmf.prob
      disc.at.diag.HI.msmf.prob <- dat$param$disc.at.diag.main.HI.msmf.prob
      disc.post.diag.HI.msmf.prob <- dat$param$disc.post.diag.main.HI.msmf.prob
      
      disc.outset.W.msmf.prob <- dat$param$disc.outset.main.W.msmf.prob
      disc.at.diag.W.msmf.prob <- dat$param$disc.at.diag.main.W.msmf.prob
      disc.post.diag.W.msmf.prob <- dat$param$disc.post.diag.main.W.msmf.prob
      
      
      el <- dat$el[[1]]
    }

    if (type == "pers") {
      #Females
      disc.outset.B.f.prob <- dat$param$disc.outset.pers.B.f.prob
      disc.at.diag.B.f.prob <- dat$param$disc.at.diag.pers.B.f.prob
      disc.post.diag.B.f.prob <- dat$param$disc.post.diag.pers.B.f.prob
      
      disc.outset.BI.f.prob <- dat$param$disc.outset.pers.BI.f.prob
      disc.at.diag.BI.f.prob <- dat$param$disc.at.diag.pers.BI.f.prob
      disc.post.diag.BI.f.prob <- dat$param$disc.post.diag.pers.BI.f.prob
 
      disc.outset.H.f.prob <- dat$param$disc.outset.pers.H.f.prob
      disc.at.diag.H.f.prob <- dat$param$disc.at.diag.pers.H.f.prob
      disc.post.diag.H.f.prob <- dat$param$disc.post.diag.pers.H.f.prob
      
      disc.outset.HI.f.prob <- dat$param$disc.outset.pers.HI.f.prob
      disc.at.diag.HI.f.prob <- dat$param$disc.at.diag.pers.HI.f.prob
      disc.post.diag.HI.f.prob <- dat$param$disc.post.diag.pers.HI.f.prob
      
      disc.outset.W.f.prob <- dat$param$disc.outset.pers.W.f.prob
      disc.at.diag.W.f.prob <- dat$param$disc.at.diag.pers.W.f.prob
      disc.post.diag.W.f.prob <- dat$param$disc.post.diag.pers.W.f.prob
      
      #Males het msf
      disc.outset.B.msf.prob <- dat$param$disc.outset.pers.B.msf.prob
      disc.at.diag.B.msf.prob <- dat$param$disc.at.diag.pers.B.msf.prob
      disc.post.diag.B.msf.prob <- dat$param$disc.post.diag.pers.B.msf.prob
      
      disc.outset.BI.msf.prob <- dat$param$disc.outset.pers.BI.msf.prob
      disc.at.diag.BI.msf.prob <- dat$param$disc.at.diag.pers.BI.msf.prob
      disc.post.diag.BI.msf.prob <- dat$param$disc.post.diag.pers.BI.msf.prob
      
      disc.outset.H.msf.prob <- dat$param$disc.outset.pers.H.msf.prob
      disc.at.diag.H.msf.prob <- dat$param$disc.at.diag.pers.H.msf.prob
      disc.post.diag.H.msf.prob <- dat$param$disc.post.diag.pers.H.msf.prob
      
      disc.outset.HI.msf.prob <- dat$param$disc.outset.pers.HI.msf.prob
      disc.at.diag.HI.msf.prob <- dat$param$disc.at.diag.pers.HI.msf.prob
      disc.post.diag.HI.msf.prob <- dat$param$disc.post.diag.pers.HI.msf.prob
      
      disc.outset.W.msf.prob <- dat$param$disc.outset.pers.W.msf.prob
      disc.at.diag.W.msf.prob <- dat$param$disc.at.diag.pers.W.msf.prob
      disc.post.diag.W.msf.prob <- dat$param$disc.post.diag.pers.W.msf.prob
      
      # Males MSM
      disc.outset.B.msm.prob <- dat$param$disc.outset.pers.B.msm.prob
      disc.at.diag.B.msm.prob <- dat$param$disc.at.diag.pers.B.msm.prob
      disc.post.diag.B.msm.prob <- dat$param$disc.post.diag.pers.B.msm.prob
      
      disc.outset.BI.msm.prob <- dat$param$disc.outset.pers.BI.msm.prob
      disc.at.diag.BI.msm.prob <- dat$param$disc.at.diag.pers.BI.msm.prob
      disc.post.diag.BI.msm.prob <- dat$param$disc.post.diag.pers.BI.msm.prob
      
      disc.outset.H.msm.prob <- dat$param$disc.outset.pers.H.msm.prob
      disc.at.diag.H.msm.prob <- dat$param$disc.at.diag.pers.H.msm.prob
      disc.post.diag.H.msm.prob <- dat$param$disc.post.diag.pers.H.msm.prob
      
      disc.outset.HI.msm.prob <- dat$param$disc.outset.pers.HI.msm.prob
      disc.at.diag.HI.msm.prob <- dat$param$disc.at.diag.pers.HI.msm.prob
      disc.post.diag.HI.msm.prob <- dat$param$disc.post.diag.pers.HI.msm.prob
      
      disc.outset.W.msm.prob <- dat$param$disc.outset.pers.W.msm.prob
      disc.at.diag.W.msm.prob <- dat$param$disc.at.diag.pers.W.msm.prob
      disc.post.diag.W.msm.prob <- dat$param$disc.post.diag.pers.W.msm.prob
      
      # Males MSMF
      disc.outset.B.msmf.prob <- dat$param$disc.outset.pers.B.msmf.prob
      disc.at.diag.B.msmf.prob <- dat$param$disc.at.diag.pers.B.msmf.prob
      disc.post.diag.B.msmf.prob <- dat$param$disc.post.diag.pers.B.msmf.prob
      
      disc.outset.BI.msmf.prob <- dat$param$disc.outset.pers.BI.msmf.prob
      disc.at.diag.BI.msmf.prob <- dat$param$disc.at.diag.pers.BI.msmf.prob
      disc.post.diag.BI.msmf.prob <- dat$param$disc.post.diag.pers.BI.msmf.prob
      
      disc.outset.H.msmf.prob <- dat$param$disc.outset.pers.H.msmf.prob
      disc.at.diag.H.msmf.prob <- dat$param$disc.at.diag.pers.H.msmf.prob
      disc.post.diag.H.msmf.prob <- dat$param$disc.post.diag.pers.H.msmf.prob
      
      disc.outset.HI.msmf.prob <- dat$param$disc.outset.pers.HI.msmf.prob
      disc.at.diag.HI.msmf.prob <- dat$param$disc.at.diag.pers.HI.msmf.prob
      disc.post.diag.HI.msmf.prob <- dat$param$disc.post.diag.pers.HI.msmf.prob
      
      disc.outset.W.msmf.prob <- dat$param$disc.outset.pers.W.msmf.prob
      disc.at.diag.W.msmf.prob <- dat$param$disc.at.diag.pers.W.msmf.prob
      disc.post.diag.W.msmf.prob <- dat$param$disc.post.diag.pers.W.msmf.prob
      
      el <- dat$el[[2]]
    }

    if (type == "inst") {
      #Females
      disc.inst.B.f.prob <- dat$param$disc.inst.B.f.prob
      disc.inst.BI.f.prob <- dat$param$disc.inst.BI.f.prob
      disc.inst.H.f.prob <- dat$param$disc.inst.H.f.prob
      disc.inst.HI.f.prob <- dat$param$disc.inst.HI.f.prob
      disc.inst.W.f.prob <- dat$param$disc.inst.W.f.prob
      
      #Males Het msf
      disc.inst.B.msf.prob <- dat$param$disc.inst.B.msf.prob
      disc.inst.BI.msf.prob <- dat$param$disc.inst.BI.msf.prob
      disc.inst.H.msf.prob <- dat$param$disc.inst.H.msf.prob
      disc.inst.HI.msf.prob <- dat$param$disc.inst.HI.msf.prob
      disc.inst.W.msf.prob <- dat$param$disc.inst.W.msf.prob

      #MSM
      disc.inst.B.msm.prob <- dat$param$disc.inst.B.msm.prob
      disc.inst.BI.msm.prob <- dat$param$disc.inst.BI.msm.prob
      disc.inst.H.msm.prob <- dat$param$disc.inst.H.msm.prob
      disc.inst.HI.msm.prob <- dat$param$disc.inst.HI.msm.prob
      disc.inst.W.msm.prob <- dat$param$disc.inst.W.msm.prob
      
      #MSMF
      disc.inst.B.msmf.prob <- dat$param$disc.inst.B.msmf.prob
      disc.inst.BI.msmf.prob <- dat$param$disc.inst.BI.msmf.prob
      disc.inst.H.msmf.prob <- dat$param$disc.inst.H.msmf.prob
      disc.inst.HI.msmf.prob <- dat$param$disc.inst.HI.msmf.prob
      disc.inst.W.msmf.prob <- dat$param$disc.inst.W.msmf.prob
      
      
      el <- dat$el[[3]]
    }


    # Processes --------------------------------------------------------------

    # Check for discordant rels in location
    posneg <- el[which((status[el[, 1]] - status[el[, 2]] == 1) & (immig.loc[el[, 1]]==0 & immig.loc[el[, 2]]==0)) , , drop = FALSE]
    negpos <- el[which((status[el[, 2]] - status[el[, 1]] == 1) & (immig.loc[el[, 2]]==0 & immig.loc[el[, 1]]==0)), , drop = FALSE]
    disc.el <- rbind(posneg, negpos[, 2:1])

    # Check for not already disclosed
    discl.list <- dat$temp$discl.list
    disclose.cdl <- discl.list[, 1] * 1e7 + discl.list[, 2]
    discord.cdl <- uid[disc.el[, 1]] * 1e7 + uid[disc.el[, 2]]
    notdiscl <- !(discord.cdl %in% disclose.cdl)


    # data frame of non-disclosed pairs
    nd <- disc.el[notdiscl, , drop = FALSE]

    # Check for positive diagnosis
    notdiscl.dx <- which(diag.status[nd[, 1]] == 1)

    # data frame of non-disclosed pairs where infected is dx'ed
    nd.dx <- nd[notdiscl.dx, , drop = FALSE]

    # If there are any eligible pairs
    if (nrow(nd.dx) > 0) {

      # Split by race, sex and sex identity of pos node
      pos.race <- race[nd.dx[, 1]]
      pos.sex <- sex[nd.dx[,1]]
      pos.sex.ident <- sex.ident[nd.dx[,1]]

      if (type %in% c("main", "pers")) {

        # Check that rel is new
        # new.edges matrix is expressed in uid, so need to transform nd.dx
        new.edges <- dat$temp$new.edges
        new.rel <- ((uid[nd.dx[, 1]] * 1e7 + uid[nd.dx[, 2]]) %in%
                      (new.edges[, 1] * 1e7 + new.edges[, 2])) |
                   ((uid[nd.dx[, 2]] * 1e7 + uid[nd.dx[, 1]]) %in%
                      (new.edges[, 1] * 1e7 + new.edges[, 2]))

        # Check if diag is new
        new.dx <- diag.time[nd.dx[, 1]] == at

        # Assign disclosure probs
        dl.prob <- vector("numeric", length = nrow(nd.dx))
        
        #Females
        dl.prob[pos.race == "B" & pos.sex == "F" & new.rel == TRUE] <- disc.outset.B.f.prob
        dl.prob[pos.race == "B" & pos.sex == "F"  & new.rel == FALSE & new.dx == TRUE] <- disc.at.diag.B.f.prob
        dl.prob[pos.race == "B" & pos.sex == "F"  & new.rel == FALSE & new.dx == FALSE] <- disc.post.diag.B.f.prob

        dl.prob[pos.race == "BI" & pos.sex == "F" & new.rel == TRUE] <- disc.outset.BI.f.prob
        dl.prob[pos.race == "BI" & pos.sex == "F"  & new.rel == FALSE & new.dx == TRUE] <- disc.at.diag.BI.f.prob
        dl.prob[pos.race == "BI" & pos.sex == "F"  & new.rel == FALSE & new.dx == FALSE] <- disc.post.diag.BI.f.prob
        
        dl.prob[pos.race == "H" & pos.sex == "F" & new.rel == TRUE] <- disc.outset.H.f.prob
        dl.prob[pos.race == "H" & pos.sex == "F"  & new.rel == FALSE & new.dx == TRUE] <- disc.at.diag.H.f.prob
        dl.prob[pos.race == "H" & pos.sex == "F"  & new.rel == FALSE & new.dx == FALSE] <- disc.post.diag.H.f.prob
        
        dl.prob[pos.race == "HI" & pos.sex == "F" & new.rel == TRUE] <- disc.outset.HI.f.prob
        dl.prob[pos.race == "HI" & pos.sex == "F"  & new.rel == FALSE & new.dx == TRUE] <- disc.at.diag.HI.f.prob
        dl.prob[pos.race == "HI" & pos.sex == "F"  & new.rel == FALSE & new.dx == FALSE] <- disc.post.diag.HI.f.prob
        
        dl.prob[pos.race == "W" & pos.sex == "F" & new.rel == TRUE] <- disc.outset.W.f.prob
        dl.prob[pos.race == "W" & pos.sex == "F"  & new.rel == FALSE & new.dx == TRUE] <- disc.at.diag.W.f.prob
        dl.prob[pos.race == "W" & pos.sex == "F"  & new.rel == FALSE & new.dx == FALSE] <- disc.post.diag.W.f.prob
        
        #Males Het msf
        dl.prob[pos.race == "B" & pos.sex == "M" & pos.sex.ident == "msf" & new.rel == TRUE] <- disc.outset.B.msf.prob
        dl.prob[pos.race == "B" & pos.sex == "M" & pos.sex.ident == "msf"  & new.rel == FALSE & new.dx == TRUE] <- disc.at.diag.B.msf.prob
        dl.prob[pos.race == "B" & pos.sex == "M" & pos.sex.ident == "msf"  & new.rel == FALSE & new.dx == FALSE] <- disc.post.diag.B.msf.prob
        
        dl.prob[pos.race == "BI" & pos.sex == "M" & pos.sex.ident == "msf" & new.rel == TRUE] <- disc.outset.BI.msf.prob
        dl.prob[pos.race == "BI" & pos.sex == "M" & pos.sex.ident == "msf"  & new.rel == FALSE & new.dx == TRUE] <- disc.at.diag.BI.msf.prob
        dl.prob[pos.race == "BI" & pos.sex == "M" & pos.sex.ident == "msf"  & new.rel == FALSE & new.dx == FALSE] <- disc.post.diag.BI.msf.prob
        
        dl.prob[pos.race == "H" & pos.sex == "M" & pos.sex.ident == "msf" & new.rel == TRUE] <- disc.outset.H.msf.prob
        dl.prob[pos.race == "H" & pos.sex == "M" & pos.sex.ident == "msf"  & new.rel == FALSE & new.dx == TRUE] <- disc.at.diag.H.msf.prob
        dl.prob[pos.race == "H" & pos.sex == "M" & pos.sex.ident == "msf"  & new.rel == FALSE & new.dx == FALSE] <- disc.post.diag.H.msf.prob
        
        dl.prob[pos.race == "HI" & pos.sex == "M" & pos.sex.ident == "msf" & new.rel == TRUE] <- disc.outset.HI.msf.prob
        dl.prob[pos.race == "HI" & pos.sex == "M" & pos.sex.ident == "msf"  & new.rel == FALSE & new.dx == TRUE] <- disc.at.diag.HI.msf.prob
        dl.prob[pos.race == "HI" & pos.sex == "M" & pos.sex.ident == "msf"  & new.rel == FALSE & new.dx == FALSE] <- disc.post.diag.HI.msf.prob
        
        dl.prob[pos.race == "W" & pos.sex == "M" & pos.sex.ident == "msf" & new.rel == TRUE] <- disc.outset.W.msf.prob
        dl.prob[pos.race == "W" & pos.sex == "M" & pos.sex.ident == "msf"  & new.rel == FALSE & new.dx == TRUE] <- disc.at.diag.W.msf.prob
        dl.prob[pos.race == "W" & pos.sex == "M" & pos.sex.ident == "msf"  & new.rel == FALSE & new.dx == FALSE] <- disc.post.diag.W.msf.prob
  
        #MSM
        dl.prob[pos.race == "B" & pos.sex == "M" & pos.sex.ident == "msm" & new.rel == TRUE] <- disc.outset.B.msm.prob
        dl.prob[pos.race == "B" & pos.sex == "M" & pos.sex.ident == "msm"  & new.rel == FALSE & new.dx == TRUE] <- disc.at.diag.B.msm.prob
        dl.prob[pos.race == "B" & pos.sex == "M" & pos.sex.ident == "msm"  & new.rel == FALSE & new.dx == FALSE] <- disc.post.diag.B.msm.prob
        
        dl.prob[pos.race == "BI" & pos.sex == "M" & pos.sex.ident == "msm" & new.rel == TRUE] <- disc.outset.BI.msm.prob
        dl.prob[pos.race == "BI" & pos.sex == "M" & pos.sex.ident == "msm"  & new.rel == FALSE & new.dx == TRUE] <- disc.at.diag.BI.msm.prob
        dl.prob[pos.race == "BI" & pos.sex == "M" & pos.sex.ident == "msm"  & new.rel == FALSE & new.dx == FALSE] <- disc.post.diag.BI.msm.prob
        
        dl.prob[pos.race == "H" & pos.sex == "M" & pos.sex.ident == "msm" & new.rel == TRUE] <- disc.outset.H.msm.prob
        dl.prob[pos.race == "H" & pos.sex == "M" & pos.sex.ident == "msm"  & new.rel == FALSE & new.dx == TRUE] <- disc.at.diag.H.msm.prob
        dl.prob[pos.race == "H" & pos.sex == "M" & pos.sex.ident == "msm"  & new.rel == FALSE & new.dx == FALSE] <- disc.post.diag.H.msm.prob
        
        dl.prob[pos.race == "HI" & pos.sex == "M" & pos.sex.ident == "msm" & new.rel == TRUE] <- disc.outset.HI.msm.prob
        dl.prob[pos.race == "HI" & pos.sex == "M" & pos.sex.ident == "msm"  & new.rel == FALSE & new.dx == TRUE] <- disc.at.diag.HI.msm.prob
        dl.prob[pos.race == "HI" & pos.sex == "M" & pos.sex.ident == "msm"  & new.rel == FALSE & new.dx == FALSE] <- disc.post.diag.HI.msm.prob
        
        dl.prob[pos.race == "W" & pos.sex == "M" & pos.sex.ident == "msm" & new.rel == TRUE] <- disc.outset.W.msm.prob
        dl.prob[pos.race == "W" & pos.sex == "M" & pos.sex.ident == "msm"  & new.rel == FALSE & new.dx == TRUE] <- disc.at.diag.W.msm.prob
        dl.prob[pos.race == "W" & pos.sex == "M" & pos.sex.ident == "msm"  & new.rel == FALSE & new.dx == FALSE] <- disc.post.diag.W.msm.prob
        
        #MSMF
        dl.prob[pos.race == "B" & pos.sex == "M" & pos.sex.ident == "msmf" & new.rel == TRUE] <- disc.outset.B.msmf.prob
        dl.prob[pos.race == "B" & pos.sex == "M" & pos.sex.ident == "msmf"  & new.rel == FALSE & new.dx == TRUE] <- disc.at.diag.B.msmf.prob
        dl.prob[pos.race == "B" & pos.sex == "M" & pos.sex.ident == "msmf"  & new.rel == FALSE & new.dx == FALSE] <- disc.post.diag.B.msmf.prob
        
        dl.prob[pos.race == "BI" & pos.sex == "M" & pos.sex.ident == "msmf" & new.rel == TRUE] <- disc.outset.BI.msmf.prob
        dl.prob[pos.race == "BI" & pos.sex == "M" & pos.sex.ident == "msmf"  & new.rel == FALSE & new.dx == TRUE] <- disc.at.diag.BI.msmf.prob
        dl.prob[pos.race == "BI" & pos.sex == "M" & pos.sex.ident == "msmf"  & new.rel == FALSE & new.dx == FALSE] <- disc.post.diag.BI.msmf.prob
        
        dl.prob[pos.race == "H" & pos.sex == "M" & pos.sex.ident == "msmf" & new.rel == TRUE] <- disc.outset.H.msmf.prob
        dl.prob[pos.race == "H" & pos.sex == "M" & pos.sex.ident == "msmf"  & new.rel == FALSE & new.dx == TRUE] <- disc.at.diag.H.msmf.prob
        dl.prob[pos.race == "H" & pos.sex == "M" & pos.sex.ident == "msmf"  & new.rel == FALSE & new.dx == FALSE] <- disc.post.diag.H.msmf.prob
        
        dl.prob[pos.race == "HI" & pos.sex == "M" & pos.sex.ident == "msmf" & new.rel == TRUE] <- disc.outset.HI.msmf.prob
        dl.prob[pos.race == "HI" & pos.sex == "M" & pos.sex.ident == "msmf"  & new.rel == FALSE & new.dx == TRUE] <- disc.at.diag.HI.msmf.prob
        dl.prob[pos.race == "HI" & pos.sex == "M" & pos.sex.ident == "msmf"  & new.rel == FALSE & new.dx == FALSE] <- disc.post.diag.HI.msmf.prob
        
        dl.prob[pos.race == "W" & pos.sex == "M" & pos.sex.ident == "msmf" & new.rel == TRUE] <- disc.outset.W.msmf.prob
        dl.prob[pos.race == "W" & pos.sex == "M" & pos.sex.ident == "msmf"  & new.rel == FALSE & new.dx == TRUE] <- disc.at.diag.W.msmf.prob
        dl.prob[pos.race == "W" & pos.sex == "M" & pos.sex.ident == "msmf"  & new.rel == FALSE & new.dx == FALSE] <- disc.post.diag.W.msmf.prob
        
      }
      
      

      if (type == "inst") {
        dl.prob <- vector("numeric", length = nrow(nd.dx))
        #Females
        dl.prob[pos.race == "B" & pos.sex == "F"] <- disc.inst.B.f.prob
        dl.prob[pos.race == "BI" & pos.sex == "F"] <- disc.inst.BI.f.prob
        dl.prob[pos.race == "H" & pos.sex == "F"] <- disc.inst.H.f.prob
        dl.prob[pos.race == "HI" & pos.sex == "F"] <- disc.inst.HI.f.prob
        dl.prob[pos.race == "W" & pos.sex == "F"] <- disc.inst.W.f.prob
        
        #Males het
        dl.prob[pos.race == "B" & pos.sex == "M" & pos.sex.ident == "msf"] <- disc.inst.B.msf.prob
        dl.prob[pos.race == "BI" & pos.sex == "M" & pos.sex.ident == "msf"] <- disc.inst.BI.msf.prob
        dl.prob[pos.race == "H" & pos.sex == "M" & pos.sex.ident == "msf"] <- disc.inst.H.msf.prob
        dl.prob[pos.race == "HI" & pos.sex == "M" & pos.sex.ident == "msf"] <- disc.inst.HI.msf.prob
        dl.prob[pos.race == "W" & pos.sex == "M" & pos.sex.ident == "msf"] <- disc.inst.W.msf.prob
        
        #MSM
        #Males 
        dl.prob[pos.race == "B" & pos.sex == "M" & pos.sex.ident == "msm"] <- disc.inst.B.msm.prob
        dl.prob[pos.race == "BI" & pos.sex == "M" & pos.sex.ident == "msm"] <- disc.inst.BI.msm.prob
        dl.prob[pos.race == "H" & pos.sex == "M" & pos.sex.ident == "msm"] <- disc.inst.H.msm.prob
        dl.prob[pos.race == "HI" & pos.sex == "M" & pos.sex.ident == "msm"] <- disc.inst.HI.msm.prob
        dl.prob[pos.race == "W" & pos.sex == "M" & pos.sex.ident == "msm"] <- disc.inst.W.msm.prob
        
        #MSMF
        #Males 
        dl.prob[pos.race == "B" & pos.sex == "M" & pos.sex.ident == "msmf"] <- disc.inst.B.msmf.prob
        dl.prob[pos.race == "BI" & pos.sex == "M" & pos.sex.ident == "msmf"] <- disc.inst.BI.msmf.prob
        dl.prob[pos.race == "H" & pos.sex == "M" & pos.sex.ident == "msmf"] <- disc.inst.H.msmf.prob
        dl.prob[pos.race == "HI" & pos.sex == "M" & pos.sex.ident == "msmf"] <- disc.inst.HI.msmf.prob
        dl.prob[pos.race == "W" & pos.sex == "M" & pos.sex.ident == "msmf"] <- disc.inst.W.msmf.prob
      }

      # Determine disclosers
      discl <- which(rbinom(length(dl.prob), 1, dl.prob) == 1)

      # Write output
      if (length(discl) > 0) {
        discl.mat <- cbind(pos = uid[nd.dx[discl, 1]],
                           neg = uid[nd.dx[discl, 2]],
                           discl.time = at)
        dat$temp$discl.list <- rbind(dat$temp$discl.list, discl.mat)
      }
    }
  }

  if (at > 2) {
    discl.list <- dat$temp$discl.list
    master.el <- rbind(dat$el[[1]], dat$el[[2]], dat$el[[3]])
    m <- which(match(discl.list[, 1] * 1e7 + discl.list[, 2],
                     uid[master.el[, 1]] * 1e7 + uid[master.el[, 2]]) |
               match(discl.list[, 2] * 1e7 + discl.list[, 1],
                     uid[master.el[, 1]] * 1e7 + uid[master.el[, 2]]))
    dat$temp$discl.list <- discl.list[m, ]
  }

  return(dat)
}
