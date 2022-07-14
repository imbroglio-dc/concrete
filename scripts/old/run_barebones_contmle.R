## running barebones contmle (Helene's repo)

# packages -------------------------------------------------------------------
# helene's repo
contmle.dir <- c("/Shared/Projects/continuousTMLE/R", "~/research/SoftWare/continuousTMLE/R")
x <- lapply(contmle.dir, function(dir) lapply(list.files(dir, full.names = TRUE), source))

library(survival)
library(zoo)
library(tidyverse)
library(data.table)
library(nleqslv)

dt3 <- sim.data2(1e3, setting = 2, no.cr = 1, competing.risk = TRUE)
# dt3$delta[dt3$delta == 2] <- 1

run <- contmle(dt = dt3, #-- dataset
               target = 1:3, #-- go after cause 1 and cause 2 specific risks
               iterative = TRUE, #-- use one-step tmle to target F1 and F2 simultaneously
               treat.effect = "ate", #-- target the ate directly
               tau = 0.3 * 1:3, #-- time-point of interest
               estimation = list(
                 "cens" = list(fit = "cox",
                               model = Surv(time, delta == 0) ~ L1 + L2 + L3 + A*L1),
                 "cause1" = list(fit = "cox",
                                 model = Surv(time, delta == 1) ~ A + L1 + L3),
                 "cause2" = list(fit = "cox",
                                 model = Surv(time, delta == 2) ~ A + L1 + L2 + L3))
)


# concrete - contmle comparison sim  ---------------------------------------------------------
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

