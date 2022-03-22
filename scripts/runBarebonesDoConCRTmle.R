library(survival); library(zoo); library(prodlim)
library(foreach); library(doRNG); library(doParallel)
library(tidyverse); library(data.table)
library(tlverse); library(sl3); library(origami)
setwd(dir = "/Shared/Projects/ConCR-TMLE/")
lapply(paste0("R/", list.files("R/")), source)
source("./R/doConCRTmle.R")
lapply(paste0("../contTMLE/R/", list.files("../contTMLE/R/")), source)
source("./R/contmle.R")


data <- as.data.table(survival::pbc)
set.seed(12345)
data[is.na(trt), trt := sample(1:2, 1)][, trt := trt - 1]
data[, status := as.numeric(status > 1)]
EventTime = data$time
EventType = data$status
Treatment = data$trt
CovDataTable = data[, c("age", "sex")]
CovTrtTime = NULL
ID = data$id
TargetTimes = 1200 #500*1:4
TargetEvents = sort(unique(data[status > 0, status]))

CVArgs = NULL
NumUpdateSteps = 100
OneStepEps = 0.5
MinNuisanceDenom = 0.05
PropScoreBackend = "sl3"
Verbose = T
GComp = TRUE

logreg <- make_learner(Lrnr_glm)
lasso <- make_learner(Lrnr_glmnet) # alpha default is 1
ridge <- Lrnr_glmnet$new(alpha = 0)
e_net <- make_learner(Lrnr_glmnet, alpha = 0.5)
a_lrnrs <- logreg # make_learner(Stack, logreg, lasso, ridge, e_net)
Intervention <- list(
    "A == 1" = function(a, L) {
        regime <- rep_len(1, length(a))
        attr(regime, "g.star") <- as.numeric(a == regime)
        return(regime)
    },
    "A == 0" = function(a, L) {
        regime <- rep_len(0, length(a))
        attr(regime, "g.star") <- as.numeric(a == regime)
        return(regime)
    })

Models <- list("Trt" = a_lrnrs,
               "0" = list(mod1 = Surv(Time, Event == 0) ~ Trt,
                          mod3 = Surv(Time, Event == 0) ~ Trt*age + Trt*sex),
               "1" = list(mod1 = Surv(Time, Event == 1) ~ Trt,
                          mod2 = Surv(Time, Event == 1) ~ Trt + age + sex,
                          mod3 = Surv(Time, Event == 1) ~ Trt*age + Trt*sex))



