
#' @title Disclosure Module for up to 5 race groups and heterosexuals and MSM.
#'
#' @description Module function for disclosure of HIV status to partners given
#'              non-disclosure in the past.
#'
#' @inheritParams aging_shamp
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


    # Parameters and network
    if (type == "main") {
     
      
      disc.outset.prob <- dat$param$disc.outset.main.prob
      disc.at.diag.prob <- dat$param$disc.at.diag.main.prob
      disc.post.diag.prob <- dat$param$disc.post.diag.main.prob

      
      el <- dat$el[[1]]
    }

    if (type == "pers") {
      #Females
      disc.outset.prob <- dat$param$disc.outset.pers.prob
      disc.at.diag.prob <- dat$param$disc.at.diag.pers.prob
      disc.post.diag.prob <- dat$param$disc.post.diag.pers.prob
      
      el <- dat$el[[2]]
    }

    if (type == "inst") {
      #Females
      disc.inst.prob <- dat$param$disc.inst.prob

      
      el <- dat$el[[3]]
    }


    # Processes --------------------------------------------------------------

    # Check for discordant rels in location
    posneg <- el[which((status[el[, 1]] - status[el[, 2]] == 1)) , , drop = FALSE]
    negpos <- el[which((status[el[, 2]] - status[el[, 1]] == 1)), , drop = FALSE]
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

      # Split by sex of pos node
      pos.sex <- sex[nd.dx[,1]]

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
        #dl.prob[pos.sex == "F" & new.rel == TRUE] <- disc.outset.prob
        #dl.prob[pos.sex == "F"  & new.rel == FALSE & new.dx == TRUE] <- disc.at.diag.prob
        #dl.prob[pos.sex == "F"  & new.rel == FALSE & new.dx == FALSE] <- disc.post.diag.prob

        #Males
        #dl.prob[pos.sex == "M" & new.rel == TRUE] <- disc.outset.prob
        #dl.prob[pos.sex == "M" & new.rel == FALSE & new.dx == TRUE] <- disc.at.diag.prob
        #dl.prob[pos.sex == "M" & new.rel == FALSE & new.dx == FALSE] <- disc.post.diag.prob
        
        dl.prob[new.rel == TRUE] <- disc.outset.prob
        dl.prob[new.rel == FALSE & new.dx == TRUE] <- disc.at.diag.prob
        dl.prob[new.rel == FALSE & new.dx == FALSE] <- disc.post.diag.prob
        
       
      }
      
      

      if (type == "inst") {
        dl.prob <- vector("numeric", length = nrow(nd.dx))
        #dl.prob[nd.dx] <- disc.inst.prob
        dl.prob <- rep(disc.inst.prob,length = nrow(nd.dx))
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
    
    if(length(dat$temp$discl.list) < 4){
    dat$temp$discl.list <- rbind(dat$temp$discl.list,c(1,2,1))
    dat$temp$discl.list <- rbind(dat$temp$discl.list,c(3,4,1))}
  }

  return(dat)
}
