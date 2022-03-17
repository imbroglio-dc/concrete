## for testing / developing the doConCRTmle function


# setup -------------------------------------------------------------------

library(MOSS); library(survtmle); library(survival); library(zoo)
library(tidyverse); library(foreach); library(doRNG); library(doParallel)
library(data.table); library(tlverse); library(sl3); library(origami)
setwd(dir = "/Shared/Projects/ConCR-TMLE/")
lapply(paste0("R/", list.files("R/")), source)
source("./R/doConCRTmle.R")

# data cleaning -------------------------------------------------------------------------------
base_data <- readxl::read_excel("./data/test_leader.xlsx") %>%
  mutate_if(is_character, as_factor) %>%
  mutate_if(~length(levels(.)) == 2, ~as.logical(as.numeric(.)-1)) %>%
  mutate(ARM = as.numeric(ARM), TIME = time_days, EVENT = event,
         SMOKER = case_when(SMOKER == "NEVER SMOKED" ~ 0,
                            SMOKER == "PREVIOUS SMOKER" ~ 1,
                            T ~ 2),
         BMIBL = case_when(is.na(BMIBL) ~ mean(BMIBL, na.rm = T),
                           T ~ BMIBL)) %>%
  dplyr::select(ARM, TIME, everything(), -subjid, -time_days, -event) %>%
  as.data.table()

data <- simulate_data(n = 1e3, base_data = base_data)


# sl parameters ----------------------------------------------------------------------

EventTime = data$TIME
EventType = data$EVENT
Treatment = data$ARM
CovDataTable = data[, -c("EVENT", "TIME", "ARM", "ID")]
CovTrtTime = NULL
ID = data$ID
TargetTimes = 500*1:3
TargetEvents = sort(unique(data[EVENT > 0, get("EVENT")]))

CVArgs = NULL
NumUpdateSteps = 100
OneStepEps = 0.5
PropScoreCutoff = 0.05
Verbose = T

logreg <- make_learner(Lrnr_glm)
lasso <- make_learner(Lrnr_glmnet) # alpha default is 1
ridge <- Lrnr_glmnet$new(alpha = 0)
e_net <- make_learner(Lrnr_glmnet, alpha = 0.5)
a_lrnrs <- make_learner(Stack, logreg, lasso, ridge, e_net)

Models <- list("A" = a_lrnrs,
               "0" = list(mod1 = Surv(Time, Event == 0) ~ Trt,
                          mod2 = Surv(Time, Event == 0) ~ Trt + AGE,
                          mod3 = Surv(Time, Event == 0) ~ Trt + AGE + SMOKER + STROKSFL,
                          mod4 = Surv(Time, Event == 0) ~ .),
               "1" = list(mod1 = Surv(Time, Event == 1) ~ Trt,
                          mod2 = Surv(Time, Event == 1) ~ Trt + SMOKER + BMIBL,
                          # mod3 = Surv(Time, Event == 1) ~ Trt*SMOKER + I(BMIBL>30)*Trt,
                          mod4 = Surv(Time, Event == 1) ~ .),
               "2" = list(mod1 = Surv(Time, Event == 2) ~ Trt,
                          mod2 = Surv(Time, Event == 2) ~ Trt + STROKSFL + MIFL,
                          # mod3 = Surv(Time, Event == 2) ~ Trt*STROKSFL + Trt*MIFL + MIFL:STROKSFL,
                          mod4 = Surv(Time, Event == 2) ~ .),
               "3" = list(mod1 = Surv(Time, Event == 3) ~ Trt,
                          mod2 = Surv(Time, Event == 3) ~ Trt + SMOKER + BMIBL,
                          mod3 = Surv(Time, Event == 3) ~ Trt + SMOKER + BMIBL + STROKSFL + MIFL,
                          mod4 = Surv(Time, Event == 3) ~ .))


# run doConCRTmle -----------------------------------------------------------------------------

debugonce(doConCRTmle)
doConCRTmle(EventTime, EventType, Treatment, CovDataTable, CovTrtTime, ID, TargetTimes,
            TargetEvents, Models, CVArgs, NumUpdateSteps, OnestepEps, PropScoreCutoff, Verbose)


# Helene's data sim ---------------------------------------------------------
library(tidyverse); library(data.table); library(zoo); library(survival); library(prodlim)
lapply(paste0("../contTMLE/R/", list.files("../contTMLE/R/")), source)
source("./R/contmle.R")

formatContmle <- function(contmleOutput) {
  tmleOutput <-
    data.table(
      "J" = rep(names(contmleOutput$tmle), each = 2),
      "val" = c("ATE", "se"),
      do.call(rbind, contmleOutput$tmle)
    )
  tmleOutput <- melt(tmleOutput, id.vars = c("J", "val"))
  setnames(tmleOutput, "variable", "time")
  tmleOutput[, J := gsub("F", "", J)]
  tmleOutput[, time := as.numeric(gsub("tau=", "", time))]
  tmleOutput <-
    dcast(tmleOutput, J + time ~ val, value.var = "value")
}

B <- 100
seeds <- seq(123456, length.out = B)
results <- vector("list", length = B)
for (i in 1:B) {
  set.seed(seeds[i])
  dt3 <- sim.data2(1e3, setting = 2, no.cr = 1, competing.risk = TRUE)
  dt3$delta[dt3$delta == 2] <- 1

  logreg <- make_learner(Lrnr_glm)
  lasso <- make_learner(Lrnr_glmnet) # alpha default is 1
  ridge <- Lrnr_glmnet$new(alpha = 0)
  e_net <- make_learner(Lrnr_glmnet, alpha = 0.5)
  a_lrnrs <- make_learner(Stack, logreg, lasso, ridge, e_net)

  Models <- list("A" = a_lrnrs,
                 "0" = list(mod1 = Surv(Time, Event == 0) ~ Trt*L1 + L1 + L2 + L3),
                 "1" = list(mod1 = Surv(Time, Event == 1) ~ Trt + L1 + L3),
                 "2" = list(mod1 = Surv(Time, Event == 2) ~ Trt + L1 + L2 + L3),
                 "3" = list(mod1 = Surv(Time, Event == 3) ~ Trt + L1 + L2))

  concreteOutput <- doConCRTmle(EventTime = dt3$time, EventType = dt3$delta, Treatment = dt3$A,
                                CovDataTable = dt3[, c("L1", "L2", "L3")], CovTrtTime = NULL,
                                ID = dt3$id, TargetTimes = 0.3*1:3, Models = Models,
                                TargetEvents = sort(unique(dt3[delta > 0, get("delta")])),
                                CVArgs = NULL, NumUpdateSteps = NumUpdateSteps,
                                OneStepEps = OneStepEps, PropScoreCutoff = 0.05, Verbose = T)
  results[[i]] <- cbind(fn = "doConCRTmle", getATE(concreteOutput))

  run <- contmle(
    dt3, #-- dataset
    target = 1:3, #-- go after cause 1 and cause 2 specific risks
    iterative = FALSE, #-- use one-step tmle to target F1 and F2 simultaneously
    treat.effect = "ate", #-- target the ate directly
    tau = 0.3 * 1:3, #-- time-point of interest
    estimation = list(
      "cens" = list(fit = "cox",
                    model = Surv(time, delta == 0) ~ L1 + L2 + L3 + A*L1),
      "cause1" = list(fit = "cox",
                      model = Surv(time, delta == 1) ~ A + L1 + L3),
      "cause2" = list(fit = "cox",
                      model = Surv(time, delta == 2) ~ A + L1 + L2 + L3),
      "cause3" = list(fit = "cox",
                      model = Surv(time, delta == 3) ~ A + L1 + L2))
  )

  results[[i]] <- rbind(results[[i]],
                        cbind(fn = "contmle",
                              formatContmle(run)))
}

bind_rows(results) %>%
  group_by(J, time, fn) %>%
  summarise(mean_ATE = mean(ATE), mean_se = mean(se))

