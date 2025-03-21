library(dplyr)
library(tidyr)

# 1. Add the following lines to mod.departure.R:
dat <- EpiModel::record_attr_history(dat, at, "tx.status", which(status == 1), tx.status[which(status == 1)])
dat <- EpiModel::record_attr_history(dat, at, "stage", which(status == 1), stage[which(status == 1)])

# 2. Run BigNets 04-netsim_interactive.R and save output:
save(sim, file = "out/sim.RData")

# 3. Load output
attr_history <- EpiModel::get_attr_history(sim)

# 4. Process output - overall mortality rate
raw <- cbind(attr_history[["stage"]][, c(2, 4, 5)], attr_history[["tx.status"]][, 5])
names(raw) <- c("ts", "uid", "stage", "tx")

inf.counts <- raw %>% group_by(ts) %>% summarise(n(), .groups = "drop_last")
aids.dep <- cbind(inf.counts, (sim[["epi"]][["dep.AIDS.on.tx"]][-1, 1] +
                                 sim[["epi"]][["dep.AIDS.off.tx"]][-1, 1]))
names(aids.dep) <- c("ts", "num", "dep")
aids.dep$mr <- aids.dep$dep / aids.dep$num

disease.mr.wk <- mean(aids.dep$mr, na.rm = TRUE)
disease.mr.yr <- disease.mr.wk * 52

# 5. Compare results
GDPH_val <- (6 + 23 + 51 + 10) / (3452 + 9889 + 16106 + 3551)
print(paste0("Actual annual disease-related mr: ", round(disease.mr.yr, 3)))
print(paste0("Expected (GDPH data): ", round(GDPH_val, 3)))

# 6. Average and max time spent in each stage/tx status
# 6a. Remove censored nodes
uncensored_ids <- raw %>%
  group_by(uid) %>%
  summarise(maxTS = max(ts), minTS = min(ts)) %>%
  filter(maxTS < max(raw$ts) & minTS > min(raw$ts))
hist_df <- raw %>% filter(uid %in% uncensored_ids$uid)

# 6b. For each node, calculate number of time steps spent in each state
hist_df <- hist_df %>%
  group_by(tx, stage, uid) %>%
  summarise(time = n(), .groups = "drop_last")

# 6c. Summarize across nodes
stageTimeSummary <- hist_df %>%
  group_by(tx, stage) %>%
  summarise(quant25Weeks = quantile(time, 0.25), avgWeeks = mean(time),
            quant75Weeks = quantile(time, 0.75), maxWeeks = max(time),
            .groups = "drop_last")

# 6d. Format and convert weeks to years
stageTimeSummary$quant25Yrs <- round(stageTimeSummary$quant25Weeks / 52, 2)
stageTimeSummary$avgYrs <- round(stageTimeSummary$avgWeeks / 52, 2)
stageTimeSummary$quant75Yrs <- round(stageTimeSummary$quant75Weeks / 52, 2)
stageTimeSummary$maxYrs <- round(stageTimeSummary$maxWeeks / 52, 2)
stageTimeSummary[, 3:6] <- data.frame(lapply(stageTimeSummary[, 3:6],
                                             function(x) round(x, 2)))
print(stageTimeSummary[, c(1, 2, 8)])
