
#' @title Partnership tracking Module
#'
#' @description Module function for tracking partnerships for STD testing
#'              and EPT.
#'
#' @inheritParams aging_msm
#'
#' @details
#' Partnerships are tracked in a persistent edge list that allows for easy
#' reference to determine if a participant has been in a particular type of
#' relationship within a defined time frame infected; or post diagnosis for
#' one recently infected. The rates of disclosure vary at these three points,
#' and also by the partnership type.
#'
#' @return
#' This function returns the \code{dat} object with the updated master
#' partnership list, on \code{temp$part.list}.
#'
#' @keywords module msm
#' @export
#'
part_msm <- function(dat, at){

  if (at < dat$param$partlist.start) {
    return(dat)
  }

  # Cycle through three partnership types
  for (type in 1:3) {

    # Variables -----------------------------------------------------------

    # Attributes
    uid <- dat$attr$uid

    # Parameters and network
    part.int <- dat$param$riskhist.int

    # pull edgelist, expressed as uid
    el <- dat$el[[type]]
    el <- matrix(uid[el], ncol = 2)


    # Processes -----------------------------------------------------------

    # STI tracking - start with existing edge list
    highlow <- el[which(el[, 1] > el[, 2]), , drop = FALSE]
    lowhigh <- el[which(el[, 1] < el[, 2]), , drop = FALSE]
    part.el <- rbind(highlow[, 2:1], lowhigh)

    # Check for not already in partnership list
    part.list <- dat$temp$part.list
    part.list <- part.list[which(part.list[, "ptype"] == type), , drop = FALSE]

    exist.partel.ids <- part.list[, 1] * 1e7 + part.list[, 2]
    check.partel.ids <- part.el[, 1] * 1e7 + part.el[, 2]
    new.part.ids <- which(!(check.partel.ids %in% exist.partel.ids))

    # matrix of dyads not yet in cumulative edgelist
    new.part.el <- part.el[new.part.ids, , drop = FALSE]

    # Write output
    if (nrow(new.part.el) > 0) {
      new.part <- cbind(uid1 = new.part.el[, 1],
                        uid2 = new.part.el[, 2],
                        ptype = type,
                        start.time = at,
                        last.active.time = at,
                        end.time = NA)

      if (type %in% 1:2) {
        # Dissolved dyads: in part.list but not in part.el *that have not already ended*
        diss.part.ids <- which(!(exist.partel.ids %in% check.partel.ids))
        toUpdate <- intersect(diss.part.ids, which(is.na(part.list[, "end.time"])))
        part.list[toUpdate, "end.time"] <- at

        # Active dyads: end.time is now or have no end.time yet
        # For those, set last.active.time to now
        last.active.now <- which(part.list[, "end.time"] == at |
                                 is.na(part.list[, "end.time"]))
        part.list[last.active.now, "last.active.time"] <- at
      }

      if (type == 3) {
        # Set end.time for all one-offs to now
        new.part[, "end.time"] <- at

        # Newly re-active one-offs: of those in current EL, also in existing PL
        # For those, reset last.active.time and end.time to now
        update.oneoff.ids <- (exist.partel.ids %in% check.partel.ids)
        if (sum(update.oneoff.ids) > 0) {
          part.list[update.oneoff.ids, c("last.active.time", "end.time")] <- at
        }
      }

      # Bind old PL and new PL
      part.list <- rbind(part.list, new.part)
    }

    # Update PL on dat$temp
    toRemove <- dat$temp$part.list[, "ptype"] == type
    dat$temp$part.list <- dat$temp$part.list[!toRemove, ]
    dat$temp$part.list <- rbind(dat$temp$part.list, part.list)
  }

  # Subset PL to current observation window
  if (at > (dat$param$partlist.start)) {
    toKeep <- which((at - (dat$temp$part.list[, "last.active.time"]) <= part.int))
    # toDrop <- which((at - (dat$temp$part.list[, "last.active.time"]) > part.int))
    dat$temp$part.list <- dat$temp$part.list[toKeep, , drop = FALSE]
  }

  return(dat)
}
