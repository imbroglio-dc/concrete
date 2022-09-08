source("./scripts/sim-data/sim_functions.R")

# true risks for 3 competing risks
# X1_full
risks1 <- getTrueRisks(assign_A = function(W, n) return(rep_len(1, n)), 
                       n = 1e3)
risks0 <- getTrueRisks(assign_A = function(W, n) return(rep_len(0, n)), 
                       n = 1e3)