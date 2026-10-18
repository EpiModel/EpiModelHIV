# 2X faster than plogis and qlogis
prob_to_log_odds <- function(p) log(p / (1 - p))
log_odds_to_prob <- function(x) 1 / (1 + exp(-x))
apply_odds_ratio <- function(p, or) log_odds_to_prob(prob_to_log_odds(p) + log(or))
