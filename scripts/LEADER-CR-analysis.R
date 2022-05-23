# load surv_tmle
surv_tmle.dir <- c("/Shared/Projects/Roadmap_CVOT/R/functions/")
x <- lapply(surv_tmle.dir, function(dir) lapply(list.files(dir, full.names = TRUE),
                                                function(x) {
                                                    if (!grepl("LEADER", x))
                                                        try(source(x), silent = TRUE)
                                                    }))

# helene's contmle repo
contmle.dir <- c("/Shared/Projects/continuousTMLE/R", "~/research/SoftWare/continuousTMLE/R")
x <- lapply(contmle.dir, function(dir) lapply(list.files(dir, full.names = TRUE), source))

# concrete
try(setwd(dir = "/Shared/Projects/ConCR-TMLE/"), silent = TRUE)
try(setwd("~/research/SoftWare/devel-tmle-survival/ConCR-TMLE/"), silent = TRUE)
x <- lapply(list.files(path = "R", full.names = TRUE), source)
source("scripts/packages.R")

try(setwd("/Shared/Projects/novo_nordisk/"))
library(tidyverse); library(here); library(survival); library(mice);
library(SuperLearner); library(survtmle); library(MOSS)
i_am("R/competing_risks/LEADER_CR_discrete.R")

target.event <- 1:3
target.time <- 3 * 1:8 / 12 * 365.25
timescale <- 3 # 3 month interval

set.seed(0)

source(file = here("R/functions/LEADER_W_clean_new.R"))

out <- as.data.table(outcomes)
setnames(out, colnames(out), c("id", 'j', 't', 'delta'))
out[, t.min := min(t), by = id]
out[delta == 0, j := 'CNSR']
out <- distinct(out)[t == t.min, ]
out[, n := .N, by = id]
out[j %in% c("NFMI", "NFSTROKE"), j := "NFMACE"]
out <- dplyr::filter(distinct(out), !(delta == 0 & n > 1))
out[, n := .N, by = id]
out[j == "CNSR", j := 0][j == "CVDEATH", j := 1][j == "NFMACE", j := 2][j == "NONCVDEATH", j := 3]
out[, j := as.numeric(j)]
out[n > 1, j := j[which.min(j)], by = id]
out <- dplyr::select(distinct(out), -c("t.min", 'n', 'delta'))
out[, t := t * 365.25 / 12]

set.seed(123456)
dt <- full_join(out, W_imputed, by = c("id" = "USUBJID")) %>% sample_n(size = 1000) %>%
    dplyr::select(t, j, ARM, everything(), -id) %>% setDT()
setnames(dt, c("t", 'j', 'ARM'), c("time", 'delta', 'A'))
dt <- as.data.frame(dt)

# estimates ----------------------------------------------------------------

# aalen-johansen ------------------------------------------------------------------------------

AJ <- survfit(Surv(time = time, event = as.factor(delta))
              ~ A, data = dt, ctype = 1) %>%
    summary(., times = target.time)

result.i <- tibble('Risk' = unlist(as.data.frame(AJ$pstate[, target.event + 1])),
                   se = unlist(as.data.frame(AJ$std.err[, target.event + 1])),
                   A = rep(rep(0:1, each = length(target.time)), times = length(target.event)),
                   Event = rep(target.event, each = length(target.time)*2),
                   Time = rep(target.time, times = length(target.event)*2)) %>%
    pivot_wider(names_from = "A", values_from = c(Risk, se), names_sep = "_") %>%
    mutate(Estimator = "Aalen-Johansen") %>%
    dplyr::select(Estimator, Event, Time, everything())

# concrete ----------------------------------------------------------------

# W <- dplyr::select(dt, -c(time, delta, A))
W <- dplyr::select(dt, c("AGE", 'AGE.imp', "SEX", "MIFL", 'HBA1CBL', 'RACE', "LDL1BL", 'LDL1BL.imp',
                         'EGFR.CKD.BL', "BMIBL", "BMIBL.imp", 'HDL1BL.imp', 'HYPFL', "HBA1CBL",
                         "WGTTBL", 'SMOKER', 'STROKSFL', "HDL1BL", 'TRIG1BL.imp', "TRIG1BL"))
mod_cens <- eval(paste0("Surv(time, delta == 0) ~ A + ", paste0(colnames(W), collapse = " + ")))
mod_1 <- eval(paste0("Surv(time, delta == 1) ~ A + ", paste0(colnames(W), collapse = " + ")))
mod_2 <- eval(paste0("Surv(time, delta == 2) ~ A + ", paste0(colnames(W), collapse = " + ")))
mod_3 <- eval(paste0("Surv(time, delta == 3) ~ A + ", paste0(colnames(W), collapse = " + ")))

# logreg <- make_learner(Lrnr_glm)
lasso <- make_learner(Lrnr_glmnet) # alpha default is 1
ridge <- Lrnr_glmnet$new(alpha = 0)
# e_net <- make_learner(Lrnr_glmnet, alpha = 0.5)
a_lrnrs <- make_learner(Stack, lasso, ridge)

models <- list("Trt" = a_lrnrs,
               "0" = list(mod1 = formula(mod_cens)),
               "1" = list(mod1 = formula(mod_1)),
               "2" = list(mod1 = formula(mod_2)),
               "3" = list(mod1 = formula(mod_3)))
intervention <- list("A == 1" = list("intervention" = function(a, L) {rep_len(1, length(a))},
                                     "g.star" = function(a, L) {as.numeric(a == 1)}),
                     "A == 0" = list("intervention" = function(a, L) {rep_len(0, length(a))},
                                     "g.star" = function(a, L) {as.numeric(a == 0)}))

concrete.args <- formatArguments(DataTable = dt,
                                 EventTime = "time", EventType = "delta",
                                 Treatment = "A", ID = NULL, Intervention = intervention,
                                 TargetTime = target.time, TargetEvent = target.event,
                                 Model = models, Verbose = TRUE)

concrete.est <- doConcrete(ConcreteArgs = concrete.args)

concrete.out <- getOutput(Estimate = concrete.est, Estimand = c("risk"), TargetTime = target.time,
                          TargetEvent = target.event, GComp = TRUE)$Risk
setnames(concrete.out$`A == 1`, c("Risk", "se"), c("Risk_1", "se_1"))
setnames(concrete.out$`A == 0`, c("Risk", "se"), c("Risk_0", "se_0"))
concrete.out <- setDT(full_join(concrete.out$`A == 0`, concrete.out$`A == 1`))
concrete.out[, Estimator := as.character(Estimator)]
concrete.out[Estimator == "tmle", Estimator := 'concrete']
concrete.out[Estimator == "gcomp", Estimator := 'c.gcomp']

result.i <- rbind(result.i, concrete.out)
saveRDS(result.i, "../ConCR-TMLE/output/leader-est.RDS")
result.i <- readRDS("../ConCR-TMLE/output/leader-est.RDS")

# # contmle -----------------------------------------------------------------
#
# run <- contmle(
#   dt, #-- dataset
#   target = target.event, #-- go after cause 1 and cause 2 specific risks
#   iterative = FALSE, #-- use one-step tmle to target F1 and F2 simultaneously
#   treat.effect = "ate", #-- target the ate directly
#   tau = target.time, #-- time-point of interest
#   estimation = list(
#     "cens" = list(fit = "cox",
#                   model = Surv(time, delta == 0) ~ A + L1 + L2 + L3 + L4 + L5),
#     "cause1" = list(fit = "cox",
#                     model = Surv(time, delta == 1) ~ A + L1 + L2 + L3 + L4 + L5),
#     "cause2" = list(fit = "cox",
#                     model = Surv(time, delta == 2) ~ A + L1 + L2 + L3 + L4 + L5),
#     "cause3" = list(fit = "cox",
#                     model = Surv(time, delta == 3) ~ A + L1 + L2 + L3 + L4 + L5)),
#   treat.model = A ~ L1 + L2 + L3 + L4 + L5,
#   verbose = FALSE,
#   sl.models = list(model = Surv(time, delta == 0) ~ A + L1 + L2 + L3 + L4 + L5))
#
# result.i <- rbind(result.i,
#                   cbind(fn = "contmle", 'Estimator' = 'tmle',
#                         formatContmle(run)))
#

# survtmle ----------------------------------------------------------------


sl_trt <- c("SL.glm", glmnets$names, "SL.bayesglm")
sl_censor <- c("SL.glm", glmnets$names, "SL.bayesglm")
sl_events <- c("SL.glm", glmnets$names, "SL.bayesglm")

timescale <- 1
tmle_sl <- surv_tmle(
    ftime = ceiling(dt$time * 12 / 365.25 / timescale),
    ftype = dt$delta,
    targets = ceiling(target.time * 12 / 365.25 / timescale),
    trt = dt$A,
    t0 = max(ceiling(target.time * 12 / 365.25 / timescale)),
    adjustVars = W,
    SL.ftime = sl_events,
    SL.ctime = sl_censor,
    SL.trt = sl_trt,
    # glm.trt = glm_trt,
    returnIC = TRUE,
    returnModels = FALSE,
    ftypeOfInterest = target.event,
    trtOfInterest = c(1, 0),
    maxIter = 50,
    method = "hazard",
)

timescale <- 3
tmle_sl3 <- surv_tmle(
    ftime = ceiling(dt$time * 12 / 365.25 / timescale),
    ftype = dt$delta,
    targets = ceiling(target.time * 12 / 365.25 / timescale),
    trt = dt$A,
    t0 = max(ceiling(target.time * 12 / 365.25 / timescale)),
    adjustVars = W,
    SL.ftime = sl_events,
    SL.ctime = sl_censor,
    SL.trt = sl_trt,
    # glm.trt = glm_trt,
    returnIC = TRUE,
    returnModels = FALSE,
    ftypeOfInterest = target.event,
    trtOfInterest = c(1, 0),
    maxIter = 50,
    method = "hazard",
)

survtmle.est <- rbind(cbind("Estimator" = "d.gcomp",
                            "A" = rep(0:1, times = length(target.event)),
                            "Event" = rep(target.event, each = 2),
                            as.data.frame(tmle_sl$init_est)),
                      cbind(Estimator = "survtmle",
                            A = rep(0:1, times = length(target.event)),
                            Event = rep(target.event, each = 2),
                            as.data.frame(tmle_sl$est)))%>% as.data.table() %>%
    melt(., id.vars = c("Estimator", "A", "Event"), value.name = "Risk", variable.name = "Time")
survtmle.est[["Time"]] <- rep(target.time, each = length(target.event)*2)
survtmle.est <- full_join(survtmle.est,
                          cbind(Estimator = "survtmle",
                                A = rep(0:1, each = length(target.event) * length(target.time)),
                                Event = rep(target.event, each = length(target.time), times = 2),
                                Time = rep(target.time, times = 2 * length(target.event)),
                                data.table(se = sqrt(diag(tmle_sl$var))))
)

survtmle.est <- dcast(survtmle.est, ... ~ A, value.var = c("Risk", "se"))
survtmle.est <- survtmle.est[, list(Event = Event, Time = Time,
                                    RD = Risk_1 - Risk_0, se = sqrt(se_1^2 + se_0^2))]
result.i <- rbind(result.i, survtmle.est)


# plot --------------------------------------------------------------------



estimates %>%
    transmute(Estimator = Estimator, t = `t`, delta = delta,
              "Risk.treated" = F1,
              "Risk.control" = F0,
              "Risk.ratio" = F1 / F0,
              "Risk.diff" = F1 - F0,
              "se.treated" = se1, "se.control"= se0,
              "se.ratio" = sqrt(1/F0^2 * se1^2 + (F1 / F0^2)^2 * se0^2),
              "se.diff" = sqrt(se1^2 + se0^2)) %>%
    pivot_longer(names_to = c("dummy", "Estimand"), names_sep = "\\.",
                 values_to = "Estimate", cols = contains("Risk")) %>%
    pivot_longer(names_to = c("dummy2", "estimand"), names_sep = "\\.",
                 values_to = "se", cols = contains("se")) %>%
    filter(Estimand == estimand) %>% dplyr::select(-c(dummy, dummy2, estimand)) %>%
    mutate(Estimand = case_when(Estimand == "control" ~ 1,
                                Estimand == "treated" ~ 2,
                                Estimand == "diff" ~ 3,
                                T ~ 4),
           delta = case_when(delta == "CVDeath" ~ 1,
                             delta == "nonfatalMACE" ~ 2,
                             T ~ 3),
           Estimand = factor(Estimand, labels = c("Control Risk", "Treated Risk", "Risk Difference",
                                                  "Relative Risk")),
           delta = factor(delta, labels = c("CV Death", "nonfatal MACE", "nonCV Death")),
           Estimator = case_when(Estimator == "Aalen-Johansen" ~ 1,
                                 Estimator == "TMLE" ~ 2,
                                 T ~ 3),
           Estimator = factor(Estimator, labels = c("Aalen-Johansen", "TMLE",
                                                    "G-Comp")),
           t = t/4) -> est

est %>% ggplot(aes(x = t, y = `Estimate`, colour = Estimator)) +
    facet_wrap(delta~Estimand, scales = "free", ncol = 4) +
    geom_errorbar(aes(ymin = Estimate - 1.96*se,
                      ymax = Estimate + 1.96*se), width = .5,
                  position = position_dodge(.5)) +
    geom_point(size = 2, position = position_dodge(.5)) + theme_bw() +
    labs(x = "Years")

estimates <- readRDS("R/competing_risks/")
estimates %>% mutate(t = t/4) %>%
    ggplot(aes(x = t, y = `Estimate`, colour = Estimator)) +
    facet_wrap(delta~Estimand, scales = "free", ncol = 4) +
    geom_errorbar(aes(ymin = Estimate - 1.96*se,
                      ymax = Estimate + 1.96*se), width = .5,
                  position = position_dodge(.5)) +
    geom_point(size = 2, position = position_dodge(.5)) + theme_bw() +
    labs(x = "Years")

est %>% arrange(`t`, `delta`, `Estimand`) %>% group_by(`t`, `delta`, `Estimand`) %>%
    mutate(Eff = se / head(se, 1)) -> est

est %>% filter(t == 3, Estimator != "G-Comp", Estimand == "Risk Difference") %>%
    mutate(delta = factor(as.numeric(delta), levels = 1:3,
                          labels = c("CV Death", "nonfatal MI or Stroke", "nonCV Death"))) %>%
    ggplot(aes(x = as.character(`t`), y = Estimate, colour = Estimator)) +
    facet_wrap(~delta, ncol = 3) +
    geom_hline(aes(yintercept = 0), alpha = 0.3, size = 1) +
    geom_errorbar(aes(ymin = Estimate - 1.96*se,
                      ymax = Estimate + 1.96*se), width = .5,
                  position = position_dodge(.5)) +
    geom_point(size = 2, position = position_dodge(.5)) + theme_bw() +
    labs(x = "Year", y = "Risk Difference") # +
# geom_label(aes(x = as.character(`t`), y = Estimate - 1.96 * `se` - 0.001,
#                label = paste0("Eff = ", round(Eff, 3)*100, "%")),
#            data = filter(est, Estimand == "Risk Difference", Estimator == "TMLE"),
#            colour = 'black')
ggsave(filename = "leader_riskdiff_yr3.png", path = "R/competing_risks/LEADER/",
       device = "png", width = 10.5, height = 6, units = "in")

