
#' @title Disclosure Module
#'
#' @description Module function for disclosure of HIV status to partners given
#'              non-disclosure in the past.
#'
#' @inheritParams aging.mard
#'
#' @details
#' Persons who are infected may disclose their status to partners at three
#' distinct time points: at relationship onset for newly formed discordant
#' pairs; at diagnosis for pairs starting as both negative but with one newly
#' infected; or post diagnosis for one recently infected. The rates of disclosure
#' vary at these three points, and also by the partnership type.
#'
#' @return
#' This function returns the \code{dat} object with the updated disclosure list,
#' on \code{temp$discl.list}.
#'
#' @keywords module
#' @export
#'
disclose.mard <- function(dat, at){

  for (type in c("main", "pers", "inst")) {

    # Variables --------------------------------------------------------------

    # Attributes
    status <- dat$attr$status
    uid <- dat$attr$uid
    diag.status <- dat$attr$diag.status
    diag.time <- dat$attr$diag.time
    race <- dat$attr$race

    # Parameters and network
    if (type == "main") {
      disc.outset.prob.B <- dat$param$disc.main.outset.prob.B
      disc.at.diag.prob.B <- dat$param$disc.main.at.diag.prob.B
      disc.post.diag.prob.B <- dat$param$disc.main.post.diag.prob.B
      disc.outset.prob.W <- dat$param$disc.main.outset.prob.W
      disc.at.diag.prob.W <- dat$param$disc.main.at.diag.prob.W
      disc.post.diag.prob.W <- dat$param$disc.main.post.diag.prob.W
      net <- dat$nw$m
      discl.type <- "M"
    }

    if (type == "pers") {
      disc.outset.prob.B <- dat$param$disc.pers.outset.prob.B
      disc.at.diag.prob.B <- dat$param$disc.pers.at.diag.prob.B
      disc.post.diag.prob.B <- dat$param$disc.pers.post.diag.prob.B
      disc.outset.prob.W <- dat$param$disc.pers.outset.prob.W
      disc.at.diag.prob.W <- dat$param$disc.pers.at.diag.prob.W
      disc.post.diag.prob.W <- dat$param$disc.pers.post.diag.prob.W
      net <- dat$nw$p
      discl.type <- "P"
    }

    if (type == "inst") {
      disc.inst.prob.B <- dat$param$disc.inst.prob.B
      disc.inst.prob.W <- dat$param$disc.inst.prob.W
      net <- dat$nw$i
      discl.type <- "I"
    }


    # Processes --------------------------------------------------------------

    # Check for discordant rels
    if (dat$control$delete.nodes == TRUE) {
      el <- matrix(as.edgelist(net), ncol = 2)
    } else {
      el <- get.dyads.active(net, at = at)
    }

    posneg <- el[which(status[el[, 1]] - status[el[, 2]] == 1), , drop = FALSE]
    negpos <- el[which(status[el[, 2]] - status[el[, 1]] == 1), , drop = FALSE]
    disc.el <- rbind(posneg, negpos[, 2:1])
    disc.el <- as.data.frame(disc.el)
    names(disc.el) <- c("pos", "neg")

    # Check for not already disclosed
    discl.list <- dat$temp$discl.list
    notdiscl <- sapply(1:nrow(disc.el), function(x) {
                       length(intersect(
                         which(uid[disc.el[x, 1]] == discl.list$pos),
                         which(uid[disc.el[x, 2]] == discl.list$neg))) == 0
                       })

    # data frame of non-disclosed pairs
    nd <- disc.el[notdiscl, , drop = FALSE]

    # Check for positive diagnosis
    notdiscl.dx <- which(diag.status[nd[, 1]] == 1)

    # data frame of non-disclosed pairs where infected is dx'ed
    nd.dx <- nd[notdiscl.dx, , drop = FALSE]

    # If there are any eligible pairs
    if (nrow(nd.dx) > 0) {

      # Split by race of pos node
      nd.dx$pos.race <- race[nd.dx[, 1]]

      if (type %in% c("main", "pers")) {

        # Check that rel is new
        # new.edges matrix is expressed in uid, so need to transform nd.dx
        new.edges <- dat$temp$new.edges
        nd.dx$new.rel <- ((uid[nd.dx[, 1]] * 1e12 + uid[nd.dx[, 2]]) %in%
                            (new.edges[, 1] * 1e12 + new.edges[, 2])) |
                         ((uid[nd.dx[, 2]] * 1e12 + uid[nd.dx[, 1]]) %in%
                            (new.edges[, 1] * 1e12 + new.edges[, 2]))

        # Check if diag is new
        nd.dx$new.dx <- diag.time[nd.dx[, 1]] == at

        # Assign disclosure probs
        dl.prob <- vector("numeric", length = nrow(nd.dx))
        dl.prob[nd.dx$pos.race == "B" &
                nd.dx$new.rel == TRUE] <- disc.outset.prob.B
        dl.prob[nd.dx$pos.race == "B" &
                nd.dx$new.rel == FALSE &
                nd.dx$new.dx == TRUE] <- disc.at.diag.prob.B
        dl.prob[nd.dx$pos.race == "B" &
                nd.dx$new.rel == FALSE &
                nd.dx$new.dx == FALSE] <- disc.post.diag.prob.B

        dl.prob[nd.dx$pos.race == "W" &
                nd.dx$new.rel == TRUE] <- disc.outset.prob.W
        dl.prob[nd.dx$pos.race == "W" &
                nd.dx$new.rel == FALSE &
                nd.dx$new.dx == TRUE] <- disc.at.diag.prob.W
        dl.prob[nd.dx$pos.race == "W" &
                nd.dx$new.rel == FALSE &
                nd.dx$new.dx == FALSE] <- disc.post.diag.prob.W
      }

      if (type == "inst") {
        dl.prob <- vector("numeric", length = nrow(nd.dx))
        dl.prob[nd.dx$pos.race == "B"] <- disc.inst.prob.B
        dl.prob[nd.dx$pos.race == "W"] <- disc.inst.prob.W
      }

      # Determine disclosers
      discl <- which(rbinom(length(dl.prob), 1, dl.prob) == 1)

      # Write output
      if (length(discl) > 0) {
        discl.df <- data.frame(pos = uid[nd.dx[discl, 1]],
                               neg = uid[nd.dx[discl, 2]],
                               discl.time = at,
                               discl.type = discl.type)
        dat$temp$discl.list <- rbind(dat$temp$discl.list, discl.df)
      }
    }
  }

  return(dat)

}
