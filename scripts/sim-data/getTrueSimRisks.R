library(tidyverse); library(data.table)
devtools::load_all(".")
source("./scripts/sim-data/sim_functions.R")

# true risks for 3 competing risks
# X1_full
risks1 <- getTrueRisks(n = 1e6, assign_A = function(W, n) return(rep_len(1, n)))
gc()
risks0 <- getTrueRisks(n = 1e6, assign_A = function(W, n) return(rep_len(0, n)))
gc()

risks <- rbind(cbind("trt" = "A=1", 
                     rbind(data.table("Event" = 1, "Time" = 1:nrow(risks1), True = risks1[[1]]), 
                           data.table("Event" = 2, "Time" = 1:nrow(risks1), True = risks1[[2]]),
                           data.table("Event" = 3, "Time" = 1:nrow(risks1), True = risks1[[3]]))), 
               cbind("trt" = "A=0", 
                     rbind(data.table("Event" = 1, "Time" = 1:nrow(risks1), True = risks0[[1]]), 
                           data.table("Event" = 2, "Time" = 1:nrow(risks1), True = risks0[[2]]),
                           data.table("Event" = 3, "Time" = 1:nrow(risks1), True = risks0[[3]]))))

# ggplot(risks, aes(x = time, y = Risk)) + geom_line(aes(colour = A, shape = ))

set.seed(0)
seeds <- sample(0:12345678, size = 240)
library(foreach)
library(doParallel)
registerDoParallel(cl = makeCluster(8))

output <- foreach(i = seeds) %dopar% {
    devtools::load_all(".")
    
    data <- simConCR(n = 5e2, random_seed = i)
    
    concreteArgs <- formatArguments(DataTable = data, 
                                    EventTime = "TIME", 
                                    EventType = "EVENT", 
                                    Treatment = "ARM",
                                    ID = "id",
                                    Intervention = 0:1, 
                                    TargetTime = 1e3)
    concreteArgs$Model$`1`$model2 <- 
        concreteArgs$Model$`2`$model2 <- 
        concreteArgs$Model$`3`$model2 <- NULL
    concreteArgs <- formatArguments(concreteArgs)
    concreteEst <- doConcrete(concreteArgs)
    concreteOut <- getOutput(concreteEst, "Risk")$Risk
    out <- bind_rows(lapply(seq_along(concreteOut), function(j) {
        cbind("trt" = paste0("A=", as.numeric(j == 1)), concreteOut[[j]])
    }))
    return(out)
}
# stopImplicitCluster()

out <- as.data.table(bind_rows(output))
out[Event == 1 & trt == "A=1", TrueRisk := risks[Time == 1000 & trt == "A=1" & Event == 1, True]]
out[Event == 2 & trt == "A=1", TrueRisk := risks[Time == 1000 & trt == "A=1" & Event == 2, True]]
out[Event == 3 & trt == "A=1", TrueRisk := risks[Time == 1000 & trt == "A=1" & Event == 3, True]]
out[Event == 1 & trt == "A=0", TrueRisk := risks[Time == 1000 & trt == "A=0" & Event == 1, True]]
out[Event == 2 & trt == "A=0", TrueRisk := risks[Time == 1000 & trt == "A=0" & Event == 2, True]]
out[Event == 3 & trt == "A=0", TrueRisk := risks[Time == 1000 & trt == "A=0" & Event == 3, True]]
result.obsgstar <- out[, list("Truth" = unique(TrueRisk), 
           "MSE" = mean((Risk - TrueRisk)^2), 
           "Bias" = mean(Risk - TrueRisk),
           "Oracle.se" = sd(Risk), 
           "Cover" = mean(Risk + 1.96*se > TrueRisk & Risk - 1.96*se < TrueRisk)), 
    by = c("trt", "Estimator", "Event")]


# half leader ---------------------------------------------------------------------------------
set.seed(0)
data <- simConCR(n = 5e2, random_seed = sample(1:1e8, 1))

concreteArgs <- formatArguments(DataTable = data, 
                                EventTime = "TIME", 
                                EventType = "EVENT", 
                                Treatment = "ARM",
                                ID = "id",
                                Intervention = 0:1, 
                                TargetTime = 1e3)

concreteEst <- doConcrete(concreteArgs)

concreteOut <- getOutput(concreteEst, "RD")$RD

