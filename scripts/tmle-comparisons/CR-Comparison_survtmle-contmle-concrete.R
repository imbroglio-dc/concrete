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

# n_cores <- 10
# registerDoParallel(n_cores)

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
if (file.exists("./output/true_risks.csv")) {
    true_risks <- as.data.table(read.csv("./output/true_risks.csv"))
} else {
    true_risks <- list("A=1" = NULL, "A=0" = NULL)
    for (a in 1:0) { # for binary treatment only
        # obs <- as.data.table(bind_rows(lapply(1:4000, function(b) PseudoLEADER)))
        obs <- PseudoLEADER
        n <- nrow(obs)
        N <- 5e3
        A <- rep(a, nrow(obs))
        set.seed(12345678)
        seeds <- sample(1:1e9, N)
        outcomes <- foreach(j = seq_along(seeds),
                            .combine = rbind) %do% {
                                seed <- seeds[j]
                                outcomes <- data.table("T1" = T1_fn(A, obs[["SMOKER"]], obs[["BMIBL"]], t1_coefs,
                                                                    output = "F_inv.u", u = runif(nrow(obs), 0, 1))$F_inv.u,
                                                       "T2" = T2_fn(A, obs[["STROKSFL"]], obs[["MIFL"]], t2_coefs,
                                                                    output = "F_inv.u", u = runif(nrow(obs), 0, 1))$F_inv.u,
                                                       "T3" = T3_fn(t3_coefs, output = "F_inv.u",
                                                                    u = runif(nrow(obs), 0, 1))$F_inv.u)
                                outcomes <- cbind("A" = A,
                                                  t(apply(outcomes, 1,function(r) c("T" = min(r),
                                                                                    "J" = which.min(r))))) %>%
                                    as.data.frame()
                                colnames(outcomes) <- c("A", "T", "J")
                                return(outcomes)
                            }

        true_risks[[paste0("A=", a)]] <- foreach(t = interval,
                                                 .combine = rbind,
                                                 .inorder = T) %do% {
                                                     tabulate(outcomes[["J"]][outcomes[["T"]] <= t]) / (n * N)
                                                 }
        colnames(true_risks[[paste0("A=", a)]]) <- paste0("F.j", 1:3, ".a", a)
    }
    rm(outcomes)
    true_risks <- rbind(
        data.table(A = 1, "time" = 1:nrow(true_risks[["A=1"]]), true_risks[["A=1"]]),
        data.table(A = 0, "time" = 1:nrow(true_risks[["A=0"]]), true_risks[["A=0"]]),
        use.names=F)
    setnames(true_risks, 3:5, paste0("F.j", 1:3))
    true_risks[, "S.t" := 1 - F.j1 - F.j2 - F.j3]
    write_csv(true_risks, "./output/true_risks.csv")
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
    tmleOutput <- cbind("Estimator" = "contmle", tmleOutput)
    return(tmleOutput)
}

B <- 1000
n <- 400
set.seed(123456)
seeds <- sample(0:1e9, size = B)
results <- vector("list", length = B)
target.time <- 2:8 * 200
target.event <- 1:3
# results <- read_rds("output/concrete-survtmle-aj.rds")

# results <- foreach(i = 1:B) %do% {
for (i in 1:B) {
    set.seed(seeds[i])
    # dt <- sim.data2(1e3, setting = 2, no.cr = 3, competing.risk = TRUE)
    dt <- simulate_data(n = n, base_data = PseudoLEADER)
    setnames(dt, c("TIME", 'EVENT', 'ARM', 'AGE', 'STROKSFL', 'SMOKER', 'BMIBL', 'MIFL'),
             c("time", "delta", 'A', "L1", 'L2', 'L3', 'L4', 'L5'))


    # aalen-johansen ----------------------------------------------------------
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

    logreg <- make_learner(Lrnr_glm)
    # lasso <- make_learner(Lrnr_glmnet) # alpha default is 1
    # ridge <- Lrnr_glmnet$new(alpha = 0)
    # e_net <- make_learner(Lrnr_glmnet, alpha = 0.5)
    a_lrnrs <- make_learner(Stack, logreg)

    models <- list("Trt" = a_lrnrs,
                   "0" = list(mod1 = Surv(time, delta == 0) ~ A + L1 + L2 + L3 + L4 + L5),
                   "1" = list(mod1 = Surv(time, delta == 1) ~ A + L1 + L2 + L3 + L4 + L5),
                   "2" = list(mod1 = Surv(time, delta == 2) ~ A + L1 + L2 + L3 + L4 + L5),
                   "3" = list(mod1 = Surv(time, delta == 3) ~ A + L1 + L2 + L3 + L4 + L5))
    intervention <- list("A == 1" = list("intervention" = function(a, L) {rep_len(1, length(a))},
                                         "g.star" = function(a, L) {as.numeric(a == 1)}),
                         "A == 0" = list("intervention" = function(a, L) {rep_len(0, length(a))},
                                         "g.star" = function(a, L) {as.numeric(a == 0)}))

    concrete.args <- formatArguments(DataTable = dt[, c("time", "delta", "A", "id",
                                                        "L1", "L2", 'L3', 'L4', "L5")],
                                     EventTime = "time", EventType = "delta",
                                     Treatment = "A", ID = "id", Intervention = intervention,
                                     TargetTime = target.time, TargetEvent = target.event,
                                     Model = models, Verbose = FALSE)

    concrete.est <- doConcrete(ConcreteArgs = concrete.args)

    concrete.out <- getOutput(Estimate = concrete.est, Estimand = c("risk"), TargetTime = target.time,
                              TargetEvent = target.event, GComp = TRUE)$Risk
    setnames(concrete.out$`A == 1`, c("Risk", "se"), c("Risk_1", "se_1"))
    setnames(concrete.out$`A == 0`, c("Risk", "se"), c("Risk_0", "se_0"))
    concrete.out <- full_join(concrete.out$`A == 0`, concrete.out$`A == 1`)
    concrete.out[, Estimator := as.character(Estimator)]
    concrete.out[Estimator == "tmle", Estimator := 'concrete']
    concrete.out[Estimator == "gcomp", Estimator := 'c.gcomp']

    result.i <- rbind(result.i, concrete.out)

    # contmle -----------------------------------------------------------------

    run <- contmle(
      dt, #-- dataset
      target = target.event, #-- go after cause 1 and cause 2 specific risks
      iterative = FALSE, #-- use one-step tmle to target F1 and F2 simultaneously
      treat.effect = "ate", #-- target the ate directly
      tau = target.time, #-- time-point of interest
      estimation = list(
        "cens" = list(fit = "cox",
                      model = Surv(time, delta == 0) ~ A + L1 + L2 + L3 + L4 + L5),
        "cause1" = list(fit = "cox",
                        model = Surv(time, delta == 1) ~ A + L1 + L2 + L3 + L4 + L5),
        "cause2" = list(fit = "cox",
                        model = Surv(time, delta == 2) ~ A + L1 + L2 + L3 + L4 + L5),
        "cause3" = list(fit = "cox",
                        model = Surv(time, delta == 3) ~ A + L1 + L2 + L3 + L4 + L5)),
      treat.model = A ~ L1 + L2 + L3 + L4 + L5,
      verbose = FALSE,
      sl.models = list(model = Surv(time, delta == 0) ~ A + L1 + L2 + L3 + L4 + L5))

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

    tmle_sl <- list()
    class(tmle_sl) <- "try-error"
    while(inherits(tmle_sl, "try-error")) {
        tmle_sl <- try(surv_tmle(
            ftime = ceiling(dt$time/200),
            ftype = dt$delta,
            targets = ceiling(target.time/200),
            trt = dt$A,
            t0 = max(ceiling(target.time/200)),
            adjustVars = as.data.frame(dt[, list(L1, L2, L3, L4, L5)]),
            SL.ftime = sl_lib_failure,
            SL.ctime = sl_lib_censor,
            SL.trt = sl_lib_g,
            # glm.trt = glm_trt,
            returnIC = TRUE,
            returnModels = TRUE,
            ftypeOfInterest = target.event,
            trtOfInterest = c(1, 0),
            maxIter = 50,
            method = "hazard",
        ))
    }


    survtmle.est <- rbind(cbind("Estimator" = "d.gcomp",
                                "A" = rep(0:1, times = length(target.event)),
                                "Event" = rep(target.event, each = 2),
                                as.data.frame(tmle_sl$init_est)),
                          cbind(Estimator = "survtmle",
                                A = rep(0:1, times = length(target.event)),
                                Event = rep(target.event, each = 2),
                                as.data.frame(tmle_sl$est))) %>% as.data.table() %>%
        melt(., id.vars = c("Estimator", "A", "Event"), value.name = "Risk", variable.name = "Time")
    survtmle.est[["Time"]] <- rep(target.time, each = length(target.event)*4)
    survtmle.est <- full_join(survtmle.est,
                              cbind(Estimator = "survtmle",
                                    A = rep(0:1, each = length(target.event) * length(target.time)),
                                    Event = rep(target.event, each = length(target.time), times = 2),
                                    Time = rep(target.time, times = 2 * length(target.event)),
                                    data.table(se = sqrt(diag(tmle_sl$var))))
    )

    survtmle.est <- dcast(survtmle.est, ... ~ A, value.var = c("Risk", "se"))
    # survtmle.est <- survtmle.est[, list(Event = Event, Time = Time,
    #                                     RD = Risk_1 - Risk_0, se = sqrt(se_1^2 + se_0^2))]
    result.i <- rbind(result.i,
                      survtmle.est)
    results[[i]] <- result.i
    saveRDS(results, "output/concrete-survtmle-aj-correctedtime.rds", compress = FALSE)
    #   return(result.i)
}

# stopImplicitCluster()
results <- read_rds("output/concrete-survtmle-aj-correctedtime.rds")


# results -----------------------------------------------------------------
truth <- melt(true_risks, id.vars = c('time', 'A'),
              measure.vars = c("F.j1", "F.j2", "F.j3"),
              variable.name = 'Event', value.name = "Risk")
truth[["Event"]] <- as.character(str_extract(truth[["Event"]], "\\d+"))
truth <- dcast(truth, time + Event ~ A, value.var = c("Risk"), sep = ".a")
truth <- truth[, list(Time = time, Event = Event, RR = `1`/`0`,
                      RD = `1` - `0`, F1 = `1`, F0 = `0`)] %>%
    melt(id.vars = c("Time", "Event"), variable.name = "Estimand", value.name = "Truth")

truth %>% ggplot(aes(x = Time, y = Truth, colour = Event)) + geom_line() +
    facet_wrap(~Estimand)

truth <- truth[Time %in% target.time, ]

res.tbl <- results[sapply(results, function(r) !is.null(results))] %>%
    bind_rows() %>% filter(Event != "S") %>%
    transmute(Estimator = Estimator, Event = Event, Time = Time,
              Event = as.character(Event),
              F1 = Risk_1, F0 = Risk_0, RR = F1/F0, RD = F1 - F0,
              "se.RR" = sqrt(1/Risk_0^2 * se_1^2 + (Risk_1 / Risk_0^2)^2 * se_0^2),
              "se.RD" = sqrt(se_1^2 + se_0^2),
              se1 = se_1, se0 = se_0) %>%
    melt(id.vars = c("Estimator", "Event", "Time"),
         measure.vars = list("Estimate" = c("F1", "F0", "RD", "RR"),
                             "se" = c("se1", "se0", "se.RD", "se.RR")),
         variable.name = "Estimand") %>%
    mutate(Estimand = case_when(Estimand == 1 ~ "F1",
                                Estimand == 2 ~ "F0",
                                Estimand == 3 ~ "RD",
                                Estimand == 4 ~ "RR")) %>%
    left_join(., truth) %>%
    group_by(Estimator, Event, Time, Estimand) %>%
    summarise(mean = mean(Estimate),
              lower = quantile(Estimate, 0.025),
              upper = quantile(Estimate, 0.975),
              cov95 = mean(Estimate + 1.96*se >= Truth &
                               Estimate - 1.96*se <= Truth),
              Truth = mean(Truth),
              MSE = mean((Truth - Estimate)^2),
              var = var(Estimate),
              bias = mean(Truth - Estimate)) %>%
    mutate(Time = as.factor(Time))

res.tbl %>%
    filter(Estimator %in% c("Aalen-Johansen", "concrete", "survtmle"),
           Estimand %in% c("RD")) %>%
    mutate(Estimand = "Absolute Risk Difference",
           Event = case_when(Event == 1 ~ "Event 1",
                             Event == 2 ~ "Event 2",
                             Event == 3 ~ "Event 3")) %>%
    ggplot(aes(x = Time, y = mean, colour = Estimator)) +
    facet_wrap(Estimand ~ Event, nrow = 3, scales = "free") + theme_minimal() +
    geom_errorbar(
        aes(
            ymin = lower,
            ymax = upper,
            colour = Estimator
        ),
        width = .75,
        position = position_dodge(.85)
    ) +
    geom_point(aes(y = mean, colour = Estimator), position = position_dodge(.85)) +
    geom_point(aes(y = upper + 0.03, colour = Estimator), alpha = 0, position = position_dodge(.85)) +
    geom_label(
      aes(
        y = upper,
        label = paste0(round(cov95, 2) * 100, "%"),
        group = Estimator
      ),
      position = position_dodge(.85),
      vjust = -.2
    ) +
    geom_segment(
        aes(
            y = Truth,
            yend = Truth,
            x = as.numeric(Time) - .4,
            xend = as.numeric(Time) + .4
        )
    ) +
    labs(x = "\nTime", y = "Simulation Mean\n")
ggsave(filename = "../ConCR-TMLE/output/concrete-survtmle-sim-plot.png", device = "png",
       height = 8, width = 14, units = "in")
