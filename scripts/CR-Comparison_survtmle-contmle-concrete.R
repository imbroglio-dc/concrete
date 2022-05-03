## compare contmle, concrete, and survtmle
# setup -------------------------------------------------------------------
library(tidyverse); library(data.table); library(foreach); library(doParallel)
library(doRNG)

# helene's repo
contmle.dir <- c("/Shared/Projects/continuousTMLE/R", "~/research/SoftWare/continuousTMLE/R")
x <- lapply(contmle.dir, function(dir) lapply(list.files(dir, full.names = TRUE), source))

# concrete
try(setwd(dir = "/Shared/Projects/ConCR-TMLE/"), silent = TRUE)
try(setwd("~/research/SoftWare/devel-tmle-survival/ConCR-TMLE/"), silent = TRUE)
x <- lapply(list.files(path = "R", full.names = TRUE), source)
source("scripts/packages.R")
source("scripts/prepare-pbc.R")
source("scripts/sim_functions.R")
intervention <- list("A == 1" = list("intervention" = function(a, L) {rep_len(1, length(a))},
                                     "g.star" = function(a, L) {as.numeric(a == 1)}),
                     "A == 0" = list("intervention" = function(a, L) {rep_len(0, length(a))},
                                     "g.star" = function(a, L) {as.numeric(a == 0)}))

# data cleaning -------------------------------------------------------------------------------
load("scripts/PseudoLEADER.RData")
PseudoLEADER <- PseudoLEADER %>%
  mutate_if(is_character, as_factor) %>%
  mutate_if(~length(levels(.)) == 2, ~as.logical(as.numeric(.)-1)) %>%
  mutate(ARM = as.numeric(ARM), TIME = time_days, EVENT = event,
         SMOKER = case_when(SMOKER == "NEVER SMOKED" ~ 0,
                            SMOKER == "PREVIOUS SMOKER" ~ 1,
                            T ~ 2),
         BMIBL = case_when(is.na(BMIBL) ~ mean(BMIBL, na.rm = T),
                           T ~ BMIBL)) %>%
  dplyr::select(ARM, TIME, everything(), -time_days, -event) %>%
  as.data.table()

data <- simulate_data(n = 1e3, base_data = PseudoLEADER)


# simulation ----------------------------------------------------------------------

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
  # dt <- sim.data2(1e3, setting = 2, no.cr = 3, competing.risk = TRUE)
  dt <- simulate_data(n = 1e3, base_data = PseudoLEADER)
  setnames(dt, c("TIME", 'EVENT', 'ARM', 'AGE', 'HBA1CBL', 'EGFMDRBC'),
           c("time", "delta", 'A', "L1", 'L2', 'L3'))

  logreg <- make_learner(Lrnr_glm)
  lasso <- make_learner(Lrnr_glmnet) # alpha default is 1
  # ridge <- Lrnr_glmnet$new(alpha = 0)
  e_net <- make_learner(Lrnr_glmnet, alpha = 0.5)
  a_lrnrs <- make_learner(Stack, logreg, lasso, e_net)

  models <- list("Trt" = a_lrnrs,
                 "0" = list(mod1 = Surv(time, delta == 0) ~ A + L1 + L2 + L3),
                 "1" = list(mod1 = Surv(time, delta == 1) ~ A + L1 + L2 + L3),
                 "2" = list(mod1 = Surv(time, delta == 2) ~ A + L1 + L2 + L3),
                 "3" = list(mod1 = Surv(time, delta == 3) ~ A + L1 + L2 + L3))
  intervention <- list("A == 1" = list("intervention" = function(a, L) {rep_len(1, length(a))},
                                       "g.star" = function(a, L) {as.numeric(a == 1)}),
                       "A == 0" = list("intervention" = function(a, L) {rep_len(0, length(a))},
                                       "g.star" = function(a, L) {as.numeric(a == 0)}))
  target.time <- quantile(dt[["time"]][dt[["delta"]] > 0], 1:3 / 4)
  target.event <- sort(unique(dt[delta > 0, delta]))

  concrete.args <- formatArguments(DataTable = dt[, c("time", "delta", "A", "id", "L1", "L2", 'L3')],
                                   EventTime = "time", EventType = "delta",
                                   Treatment = "A", ID = "id", Intervention = intervention,
                                   TargetTime = target.time, TargetEvent = target.event,
                                   Model = models, Verbose = TRUE)
  concrete

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

