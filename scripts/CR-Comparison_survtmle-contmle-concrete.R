## compare contmle, concrete, and survtmle
# setup -------------------------------------------------------------------
library(tidyverse); library(data.table); library(foreach); library(doParallel)
library(doRNG); library(SuperLearner); library(survtmle); library(MOSS)


# surv_tmle ---------------------------------------------------------------
surv_tmle.dir <- c("/Shared/Projects/Roadmap_CVOT/R/functions/")
x <- lapply(surv_tmle.dir, function(dir) lapply(list.files(dir, full.names = TRUE),
                                                function(x) try(source(x), silent = TRUE)))

# helene's contmle repo
contmle.dir <- c("/Shared/Projects/continuousTMLE/R", "~/research/SoftWare/continuousTMLE/R")
x <- lapply(contmle.dir, function(dir) lapply(list.files(dir, full.names = TRUE), source))

# concrete
try(setwd(dir = "/Shared/Projects/ConCR-TMLE/"), silent = TRUE)
try(setwd("~/research/SoftWare/devel-tmle-survival/ConCR-TMLE/"), silent = TRUE)
x <- lapply(list.files(path = "R", full.names = TRUE), source)
source("scripts/packages.R")
source("scripts/prepare-pbc.R")
source("scripts/sim_functions.R")

n_cores <- 10
registerDoParallel(n_cores)

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


# get true risks ----------------------------------------------------------
if (file.exists("./data/true_risks.csv")) {
  true_risks <- as.data.table(read.csv("./data/true_risks.csv"))
} else {
  true_risks <- list("A=1" = NULL, "A=0" = NULL)
  for (a in 1:0) { # for binary treatment only
    obs <- as.data.table(bind_rows(lapply(1:8000, function(b) PseudoLEADER)))
    n <- nrow(obs)
    A <- rep(a, nrow(obs))
    outcomes <- data.table("T1" = T1_fn(A, obs[["SMOKER"]], obs[["BMIBL"]], t1_coefs,
                                        output = "F_inv.u", u = runif(nrow(obs), 0, 1))$F_inv.u,
                           "T2" = T2_fn(A, obs[["STROKSFL"]], obs[["MIFL"]], t2_coefs,
                                        output = "F_inv.u", u = runif(nrow(obs), 0, 1))$F_inv.u,
                           "T3" = T3_fn(t3_coefs, output = "F_inv.u",
                                        u = runif(nrow(obs), 0, 1))$F_inv.u)
    outcomes <- cbind("A" = A,
                      t(apply(outcomes, 1,function(r) c("T" = ceiling(min(r)),
                                                        "J" = which.min(r)))))
    rm(obs); rm(A); gc()

    true_risks[[paste0("A=", a)]] <- foreach(t = interval,
                                             .combine = rbind,
                                             .inorder = T) %dopar% {
                                               tabulate(outcomes[["J"]][outcomes[["T"]] <= t])
                                             }
    true_risks[[paste0("A=", a)]] <- as.data.table(true_risks[[paste0("A=", a)]] / n) %>%
      rename_all(~paste0("F.j", 1:3, ".a", a))
  }
  rm(outcomes);
  true_risks <- rbind(
    data.table(A = 1, "time" = 1:nrow(true_risks[["A=1"]]), true_risks[["A=1"]]),
    data.table(A = 0, "time" = 1:nrow(true_risks[["A=0"]]), true_risks[["A=0"]]),
    use.names=F)
  setnames(true_risks, 3:5, paste0("F.j", 1:3))
  true_risks[, "S.t" := 1 - F.j1 - F.j2 - F.j3]
  write_csv(true_risks, "data/true_risks.csv")
}


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
  tmleOutput <- dcast(tmleOutput, J + time ~ val, value.var = "value")
  setnames(tmleOutput, c("J", "time", 'ATE'), c("Event", "Time", "RD"))
}

B <- 200
seeds <- sample(0:1e9, size = B)
results <- vector("list", length = B)
target.time <- 1:4 * 400
target.event <- 1:3

results <- foreach(i = 1:B, .combine = rbind) %dopar% {
  set.seed(seeds[i])
  # dt <- sim.data2(1e3, setting = 2, no.cr = 3, competing.risk = TRUE)
  dt <- simulate_data(n = 500, base_data = PseudoLEADER)
  setnames(dt, c("TIME", 'EVENT', 'ARM', 'AGE', 'HBA1CBL', 'EGFMDRBC'),
           c("time", "delta", 'A', "L1", 'L2', 'L3'))

  # concrete ----------------------------------------------------------------

  logreg <- make_learner(Lrnr_glm)
  # lasso <- make_learner(Lrnr_glmnet) # alpha default is 1
  # ridge <- Lrnr_glmnet$new(alpha = 0)
  # e_net <- make_learner(Lrnr_glmnet, alpha = 0.5)
  a_lrnrs <- make_learner(Stack, logreg)

  models <- list("Trt" = a_lrnrs,
                 "0" = list(mod1 = Surv(time, delta == 0) ~ A + L1 + L2 + L3),
                 "1" = list(mod1 = Surv(time, delta == 1) ~ A + L1 + L2 + L3),
                 "2" = list(mod1 = Surv(time, delta == 2) ~ A + L1 + L2 + L3),
                 "3" = list(mod1 = Surv(time, delta == 3) ~ A + L1 + L2 + L3))
  intervention <- list("A == 1" = list("intervention" = function(a, L) {rep_len(1, length(a))},
                                       "g.star" = function(a, L) {as.numeric(a == 1)}),
                       "A == 0" = list("intervention" = function(a, L) {rep_len(0, length(a))},
                                       "g.star" = function(a, L) {as.numeric(a == 0)}))


  concrete.args <- formatArguments(DataTable = dt[, c("time", "delta", "A", "id", "L1", "L2", 'L3')],
                                   EventTime = "time", EventType = "delta",
                                   Treatment = "A", ID = "id", Intervention = intervention,
                                   TargetTime = target.time, TargetEvent = target.event,
                                   Model = models, Verbose = TRUE)
  concrete.est <- doConcrete(ConcreteArgs = concrete.args)

  concrete.ate <- getOutput(Estimate = concrete.est, Estimand = c("rd"), TargetTime = target.time,
                            TargetEvent = target.event, GComp = TRUE)$RD

  result.i <- cbind(fn = "concrete", concrete.ate)

  # contmle -----------------------------------------------------------------

  run <- contmle(
    dt, #-- dataset
    target = target.event, #-- go after cause 1 and cause 2 specific risks
    iterative = FALSE, #-- use one-step tmle to target F1 and F2 simultaneously
    treat.effect = "ate", #-- target the ate directly
    tau = target.time, #-- time-point of interest
    estimation = list(
      "cens" = list(fit = "cox",
                    model = Surv(time, delta == 0) ~ A + L1 + L2 + L3),
      "cause1" = list(fit = "cox",
                      model = Surv(time, delta == 1) ~ A + L1 + L2 + L3),
      "cause2" = list(fit = "cox",
                      model = Surv(time, delta == 2) ~ A + L1 + L2 + L3),
      "cause3" = list(fit = "cox",
                      model = Surv(time, delta == 3) ~ A + L1 + L2 + L3))
  )

  result.i <- rbind(result.i,
                    cbind(fn = "contmle", 'Estimator' = 'tmle',
                          formatContmle(run)))


  # survtmle ----------------------------------------------------------------

  screeners <- "All"

  sl_lib_g <- expand.grid(c("SL.glm"), screeners)
  sl_lib_g <- lapply(1:nrow(sl_lib_g),
                     function(i) as.character(unlist(sl_lib_g[i,])))

  sl_lib_censor <-
    expand.grid(c("SL.glm", "SL.glmnet"), screeners)
  sl_lib_censor <- lapply(1:nrow(sl_lib_censor),
                          function(i) as.character(unlist(sl_lib_censor[i,])))

  sl_lib_failure <-
    expand.grid(c("SL.glm", "SL.glmnet"), screeners)
  sl_lib_failure <- lapply(1:nrow(sl_lib_failure),
                           function(i) as.character(unlist(sl_lib_failure[i,])))

  sl_fit <- my_init_sl_fit(
    T_tilde = ceiling(dt$time/80),
    Delta = as.numeric(dt$delta),
    A = as.numeric(dt$A),
    W = as.data.frame(dt[, list(L1, L2, L3)]),
    t_max = max(ceiling(target.time/80)),
    sl_failure = sl_lib_failure,
    sl_censoring = sl_lib_censor,
    sl_treatment = "SL.glm",
    cv.Control = list(V = 10)
  )

  sl_fit$models$A$env <- sl_fit$models$C$env <- sl_fit$models$Y$J1$env <- NULL

  haz_sl <- list(sl_fit$density_failure_1$clone(),
                 sl_fit$density_failure_0$clone())
  haz_sl[[1]]$haz2surv()
  haz_sl[[2]]$haz2surv()
  names(haz_sl) <- c("A = 1", "A = 0")

  SL_ftime <- sl_fit$models$Y
  sl_G_dC <- sl_fit$G_dC
  # glm_trt <- paste0(colnames(adjust_vars), collapse = " + ")
  rm(sl_fit)

  tmle_sl <- surv_tmle(
    ftime = ceiling(dt$time/80),
    ftype = dt$delta,
    targets = ceiling(target.time/80),
    trt = dt$A,
    t0 = max(ceiling(target.time/80)),
    adjustVars = as.data.frame(dt[, list(L1, L2, L3)]),
    SL.ftime = SL_ftime,
    SL.ctime = sl_G_dC,
    SL.trt = sl_lib_g,
    # glm.trt = glm_trt,
    returnIC = TRUE,
    returnModels = TRUE,
    ftypeOfInterest = target.event,
    trtOfInterest = c(1, 0),
    maxIter = 20,
    method = "hazard"
  )

  survtmle.out <- cbind(A = rep(0:1, times = length(target.event)),
                        Event = rep(target.event, each = length(target.time)),
                        tmle_sl$est) %>% as.data.table()
  survtmle.out <- melt(data = survtmle.out, id.vars = c("A", "Event"),
                       variable.name = "Time", value.name = "Risk")
  survtmle.out[["Time"]] <- as.numeric(str_extract(survtmle.out[["Time"]], '\\d+')) * 80
  survtmle.out <- full_join(survtmle.out,
                            data.frame(A = rep(0:1, each = length(target.time) * length(target.event)),
                                       Time = rep(target.time, times = length(target.event) * 2),
                                       Event = rep(1:3, each = length(target.time)),
                                       'se' = sqrt(diag(tmle_sl$var))))
  survtmle.out <- dcast(survtmle.out, ... ~ A, value.var = c("Risk", "se"))
  survtmle.out <- survtmle.out[, list(Event = Event, Time = Time,
                                      RD = Risk_1 - Risk_0, se = sqrt(se_1^2 + se_0^2))]
  result.i <- rbind(result.i,
                    cbind(fn = "survtmle", Estimator = "tmle",
                          survtmle.out))

  return(result.i)
}
stopCluster()

# obsolete ----------------------------------------------------------------


# for (i in 1:B) {
#   set.seed(seeds[i])
#   # dt <- sim.data2(1e3, setting = 2, no.cr = 3, competing.risk = TRUE)
#   dt <- simulate_data(n = 400, base_data = PseudoLEADER)
#   setnames(dt, c("TIME", 'EVENT', 'ARM', 'AGE', 'HBA1CBL', 'EGFMDRBC'),
#            c("time", "delta", 'A', "L1", 'L2', 'L3'))
#
#   logreg <- make_learner(Lrnr_glm)
#   # lasso <- make_learner(Lrnr_glmnet) # alpha default is 1
#   # ridge <- Lrnr_glmnet$new(alpha = 0)
#   # e_net <- make_learner(Lrnr_glmnet, alpha = 0.5)
#   a_lrnrs <- make_learner(Stack, logreg)
#
#   models <- list("Trt" = a_lrnrs,
#                  "0" = list(mod1 = Surv(time, delta == 0) ~ A + L1 + L2 + L3),
#                  "1" = list(mod1 = Surv(time, delta == 1) ~ A + L1 + L2 + L3),
#                  "2" = list(mod1 = Surv(time, delta == 2) ~ A + L1 + L2 + L3),
#                  "3" = list(mod1 = Surv(time, delta == 3) ~ A + L1 + L2 + L3))
#   intervention <- list("A == 1" = list("intervention" = function(a, L) {rep_len(1, length(a))},
#                                        "g.star" = function(a, L) {as.numeric(a == 1)}),
#                        "A == 0" = list("intervention" = function(a, L) {rep_len(0, length(a))},
#                                        "g.star" = function(a, L) {as.numeric(a == 0)}))
#
#
#   concrete.args <- formatArguments(DataTable = dt[, c("time", "delta", "A", "id", "L1", "L2", 'L3')],
#                                    EventTime = "time", EventType = "delta",
#                                    Treatment = "A", ID = "id", Intervention = intervention,
#                                    TargetTime = target.time, TargetEvent = target.event,
#                                    Model = models, Verbose = TRUE)
#   concrete.est <- doConcrete(ConcreteArgs = concrete.args)
#
#   concrete.ate <- getOutput(Estimate = concrete.est, Estimand = c("rd"), TargetTime = target.time,
#                             TargetEvent = target.event, GComp = TRUE)$RD
#
#   results[[i]] <- cbind(fn = "concrete", concrete.ate)
#
#   run <- contmle(
#     dt, #-- dataset
#     target = target.event, #-- go after cause 1 and cause 2 specific risks
#     iterative = FALSE, #-- use one-step tmle to target F1 and F2 simultaneously
#     treat.effect = "ate", #-- target the ate directly
#     tau = target.time, #-- time-point of interest
#     estimation = list(
#       "cens" = list(fit = "cox",
#                     model = Surv(time, delta == 0) ~ A + L1 + L2 + L3),
#       "cause1" = list(fit = "cox",
#                       model = Surv(time, delta == 1) ~ A + L1 + L2 + L3),
#       "cause2" = list(fit = "cox",
#                       model = Surv(time, delta == 2) ~ A + L1 + L2 + L3),
#       "cause3" = list(fit = "cox",
#                       model = Surv(time, delta == 3) ~ A + L1 + L2 + L3))
#   )
#
#   results[[i]] <- rbind(results[[i]],
#                         cbind(fn = "contmle", 'Estimator' = 'tmle',
#                               formatContmle(run)))
# }


# results -----------------------------------------------------------------
truth <- melt(true_risks, id.vars = c('time', 'A'),
              measure.vars = c("F.j1", "F.j2", "F.j3"),
              variable.name = 'Event', value.name = "Risk")
truth[["Event"]] <- as.character(str_extract(truth[["Event"]], "\\d+"))
truth <- dcast(truth, time + Event ~ A, value.var = c("Risk"), sep = ".a")
truth <- truth[, list(Time = time, Event = Event, trueRD = `1` - `0`)]

truth %>% ggplot(aes(x = Time, y = trueRD, colour = Event)) + geom_line()

truth <- truth[Time %in% target.time, ]
bind_rows(results) %>% filter(Event != "S") %>%
  mutate(Event = as.character(Event),
         fn = paste(fn, Estimator, sep = "-")) %>%
  group_by(Event, Time, fn, Estimator) %>%
  full_join(., truth) %>%
  summarise(mean_ATE = mean(RD), mean_se = mean(se),
            cov = mean(RD + 1.96*se >= trueRD &
                         RD - 1.96*se <= trueRD)) %>%
  mutate(Time = as.factor(Time),
         lower = mean_ATE - 1.96 * mean_se,
         upper = mean_ATE + 1.96 * mean_se) %>%
  ggplot(aes(x = Time, y = mean_ATE, colour = fn)) +
  facet_wrap(~Event) + theme_minimal() +
  geom_errorbar(
    aes(
      ymin = lower,
      ymax = upper,
      colour = fn
    ),
    width = .5,
    position = position_dodge(.5)
  ) +
  geom_point(aes(y = mean_ATE, colour = fn), position = position_dodge(0.5)) +
  geom_label(
    aes(
      y = upper,
      label = paste0(round(cov, 2) * 100, "%"),
      group = fn
    ),
    position = position_dodge(0.5),
    vjust = -.2
  )

