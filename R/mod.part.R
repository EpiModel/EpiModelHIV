
#' @title Partnership tracking Module
#'
#' @description Module function for tracking partnerships for STD testing and EPT.
#'
#' @inheritParams aging_msm
#'
#' @details
#' Partnerships are tracked in a persistent edge list that allows for easy
#' reference to determine if a participant has been in a particular type of
#' relationship within a defined time frame
#' infected; or post diagnosis for one recently infected. The rates of disclosure
#' vary at these three points, and also by the partnership type.
#'
#' @return
#' This function returns the \code{dat} object with the updated master partnership list,
#' on \code{temp$part.list}.
#'
#' @keywords module msm
#' @export
#'
part_msm <- function(dat, at){
    for (type in c("main", "pers", "inst")) {
        
        # Variables --------------------------------------------------------------
        
        # Attributes
        uid <- dat$attr$uid
        
        # Parameters and network
        part.int <- dat$param$sti.highrisktest.int
        if (type == "main") {
            el <- dat$el[[1]]
        }
        
        if (type == "pers") {
            el <- dat$el[[2]]
        }
        
        if (type == "inst") {
            el <- dat$el[[3]]
        }
        
        
        # Processes --------------------------------------------------------------
        # STI tracking
        # Order with lowest uid value first - just to create matrix
        highlow <- el[which(uid[el[, 1]] > uid[el[, 2]]), , drop = FALSE]
        lowhigh <- el[which(uid[el[, 1]] < uid[el[, 2]]), , drop = FALSE]
        part.el <- rbind(highlow[, 2:1], lowhigh)
        
        # Check for not already in partnership list
        part.list <- dat$temp$part.list
        exist.cdl <- part.list[, 1] * 1e7 + part.list[, 2]
        check.cdl <- uid[part.el[, 1]] * 1e7 + uid[part.el[, 2]]
        notpartlist <- !(check.cdl %in% exist.cdl)
        
        # data frame of pairs not yet in edgelist
        notyet <- part.el[notpartlist, , drop = FALSE]
        
        # If there are any eligible pairs to add
        if (nrow(notyet) > 0) {
            
            if (type %in% c("main", "pers", "inst")) {
                
                # Assign partnership type
                if (type == "main") {
                    parttype = 1
                }
                if (type == "pers") {
                    parttype = 2
                }
                if (type == "inst") {
                    parttype = 3
                }
                
                # Check that rel is new
                # new.edges matrix is expressed in uid, so need to transform notyet
                new.edges <- dat$temp$new.edges #only includes types 1 and 2 so far
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
                              ptype = parttype,
                              start.time = at,
                              last.active.time = at,
                              end.time = NA)
            dat$temp$part.list <- rbind(dat$temp$part.list, new.part)
            
            if (type %in% "inst") {
                
            # Instantaneous - last.active.time and end.time columnns get value of start.time
                dat$temp$part.list[which(dat$temp$part.list[, 3] == 3), 5:6] <- dat$temp$part.list[which(dat$temp$part.list[, 3] == 3), 4]
        }
        }
    }
    
    if (at > 2) {
        
        # Existing edges to reference against partnership list
        master.el <- rbind(dat$el[[1]], dat$el[[2]], dat$el[[3]])

        # Partnership tracking - last x months
        part.list <- dat$temp$part.list
        
        # Select matching (currently in edgelist and existing temp$part.list) partnerships to update last active date of partnership
        m2 <- which(match(part.list[, 1] * 1e7 + part.list[, 2],
                          uid[master.el[, 1]] * 1e7 + uid[master.el[, 2]]) |
                        match(part.list[, 2] * 1e7 + part.list[, 1],
                              uid[master.el[, 1]] * 1e7 + uid[master.el[, 2]]))
        dat$temp$part.list[existing, 5] <- at
        
        # Update for edges that do not have an end date (is.na(partlist[, "end.time"]))
        # part.list2 <- part.list[is.na(part.list[ , 6]), ]
        
        # Add partnership end dates for non-instantaneous
        
        # Subset part.list to include only partnerships active in last x months
        dat$temp$part.list <- dat$temp$part.list[(at-(dat$temp$part.list[, 5]) <= part.int), ]
    }
return(dat)    
}

    