
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
    part.int <- dat$param$sti.highrisktest.int

    # pull edgelist
    el <- dat$el[[type]]


    # Processes -----------------------------------------------------------
    # STI tracking - start with existing edge list
    # Order with lowest uid value first - just to create matrix
    highlow <- el[which(uid[el[, 1]] > uid[el[, 2]]), , drop = FALSE]
    lowhigh <- el[which(uid[el[, 1]] < uid[el[, 2]]), , drop = FALSE]
    part.el <- rbind(highlow[, 2:1], lowhigh)

    # Check for not already in partnership list
    part.list <- dat$temp$part.list
    exist.partel <- part.list[, 1] * 1e7 + part.list[, 2]
    check.partel <- uid[part.el[, 1]] * 1e7 + uid[part.el[, 2]]
    notpartlist <- !(check.partel %in% exist.partel)

    # data frame of pairs not yet in edgelist
    notyet <- part.el[notpartlist, , drop = FALSE]

    # If there are any eligible pairs to add
    # TODO: why checking against new.edges matrix? Shouldn't notyet be it?
    #       this also doesn't work if at == 1
    if (nrow(notyet) > 0) {

      # TODO: is that conditional check necessary?
      if (type %in% 1:3) {
          # Assign partnership type
          parttype <- type

          # new.edges matrix is expressed in uid, notyet vs new.edges
          new.edges <- dat$temp$new.edges # only includes types 1 and 2
          new.rel <- ((uid[notyet[, 1]] * 1e7 + uid[notyet[, 2]]) %in%
                        (new.edges[, 1] * 1e7 + new.edges[, 2])) |
                     ((uid[notyet[, 2]] * 1e7 + uid[notyet[, 1]]) %in%
                        (new.edges[, 1] * 1e7 + new.edges[, 2]))
      }
    }

    # Write output
    if (length(new.rel) > 0) {
      new.part <- cbind(uid1 = uid[notyet[, 1]],
                        uid2 = uid[notyet[, 2]],
                        ptype = parttype, # TODO: can you just input type here, instead of creating new parttype?
                        start.time = at,
                        last.active.time = at,
                        end.time = NA)
      dat$temp$part.list <- rbind(dat$temp$part.list, new.part)

      if (type == 3) {

        # One-off: last.active.time and end.time columns get value of start.time
        dat$temp$part.list[which(dat$temp$part.list[, 3] == 3), 5:6] <-
                      dat$temp$part.list[which(dat$temp$part.list[, 3] == 3), 4]
      }
    }
  }

  # if partlist already exists, update it
  if (at > (dat$param$partlist.start)) {

    # Existing edges to reference against partnership list
    master.el <- rbind(dat$el[[1]], dat$el[[2]], dat$el[[3]])

    # Partnership tracking - last x months
    part.list <- dat$temp$part.list

    # Add partnership end dates for non-instantaneous
    dead.edges.m <- attributes(dat$el[[1]])$changes
    dead.edges.m <- dead.edges.m[dead.edges.m[, "to"] == 0, 1:2, drop = FALSE]
    dead.edges.p <- attributes(dat$el[[2]])$changes
    dead.edges.p <- dead.edges.p[dead.edges.p[, "to"] == 0, 1:2, drop = FALSE]
    dead.edges <- rbind(dead.edges.m, dead.edges.p)

    dead.rel <- (
      (uid[dead.edges[, 1]] * 1e7 + uid[dead.edges[, 2]]) %in% (part.list[, 1] * 1e7 + part.list[, 2])) |
        ((uid[dead.edges[, 2]] * 1e7 + uid[dead.edges[, 1]]) %in% (part.list[, 1] * 1e7 + part.list[, 2]))

    # Set dead edges to have ended at this timepoint
    if (length(dead.rel) > 0) {
        part.list[which(
          (match(part.list[, 1] * 1e7 + part.list[, 2], uid[dead.edges[, 1]] * 1e7 + uid[dead.edges[, 2]]) |
           match(part.list[, 2] * 1e7 + part.list[, 1], uid[dead.edges[, 1]] * 1e7 + uid[dead.edges[, 2]]))), 6] <- at
    }

    # Select matching ((currently in both edgelist and existing
    # dat$temp$part.list) and no end date yet) partnerships to update
    # last active date of partnership
    part.list[which(
      (match(part.list[, 1] * 1e7 + part.list[, 2], uid[master.el[, 1]] * 1e7 + uid[master.el[, 2]]) |
       match(part.list[, 2] * 1e7 + part.list[, 1], uid[master.el[, 1]] * 1e7 + uid[master.el[, 2]])) &
       is.na(part.list[, 6])), 5] <- at

    # Subset part.list to include only partnerships active in last x months
    part.list <- part.list[which((at - (part.list[, 5]) <= part.int)), , drop = FALSE]
    dat$temp$part.list <- part.list
  }

  return(dat)
}
