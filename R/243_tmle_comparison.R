
# setup ---------------------------------------------------------------------------------------

library(tidyverse); library(sl3); library(tmle3); library(survtmle); library(MOSS)
library(data.table); library(foreach); library(doParallel); library(doRNG); library(parallel)
library(survival); library(zoo)

## Edit to suit local working directory structure ---------------------------------------------
setwd("/Shared/Projects/ConCR-TMLE/")

# source("./R/contmle-competing-risks-simulation.R")
## just need base_data and true_risks from above
## get base_data ----
set.seed(0)
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

## get true_risks ----
## True risks approximated by computing risks in a very large sample
if (file.exists("./data/true_risks.csv")) {
    true_risks <- read.csv("./data/true_risks.csv")
} else {
    true_risks <- list("A=1" = NULL, "A=0" = NULL)
    for (a in 1:0) { # for binary treatment only
        obs <- as.data.table(bind_rows(lapply(1:5000, function(b) base_data)))
        A <- rep(a, nrow(obs))
        outcomes <- data.table("T1" = T1_fn(A, obs[["SMOKER"]], obs[["BMIBL"]], t1_coefs,
                                            output = "F_inv.u", u = runif(nrow(obs), 0, 1))$F_inv.u,
                               "T2" = T2_fn(A, obs[["STROKSFL"]], obs[["MIFL"]], t2_coefs,
                                            output = "F_inv.u", u = runif(nrow(obs), 0, 1))$F_inv.u,
                               "T3" = T3_fn(t3_coefs, output = "F_inv.u",
                                            u = runif(nrow(obs), 0, 1))$F_inv.u)
        outcomes <- outcomes %>%
            mutate("1>2" = (T1 - T2) > 0,
                   "2>3" = (T2 - T3) > 0,
                   "3>1" = (T3 - T1) > 0,
                   "T" = T1*`3>1`*(1 - `1>2`) + T2*`1>2`*(1 - `2>3`) + T3*`2>3`*(1 - `3>1`),
                   "J" = `3>1`*(1 - `1>2`) + 2*`1>2`*(1 - `2>3`) + 3*`2>3`*(1 - `3>1`)) %>%
            dplyr::select(`T`, J)
        true_risks[[paste0("A=", a)]] <- foreach(t = interval,
                                                 .combine = rbind,
                                                 .inorder = T) %dopar% {
                                                     tabulate(outcomes[["J"]][outcomes[["T"]] <= t])
                                                 }
        true_risks[[paste0("A=", a)]] <- as.data.table(true_risks[[paste0("A=", a)]] / nrow(obs)) %>%
            rename_all(~paste0("F.j", 1:3, ".a", a))
    }
    rm(outcomes); rm(obs); rm(A); gc()
    true_risks <- rbind(
        data.table(A = 1, "time" = 1:nrow(true_risks[["A=1"]]), true_risks[["A=1"]]), 
        data.table(A = 0, "time" = 1:nrow(true_risks[["A=0"]]), true_risks[["A=0"]]), 
        use.names=F)
    setnames(true_risks, 3:5, paste0("F.j", 1:3))
    true_risks[, "S.t" := 1 - F.j1 - F.j2 - F.j3]
    write_csv(true_risks, "data/true_risks.csv")
}

source(file = "R/functions/my_sl_functions.R")
source(file = "R/functions/my_MOSS_hazard_methods.R")
source(file = "R/functions/my_survtmle_functions.R")
source(file = "R/functions/sim_functions.R")
source(file = "R/contmle.R")
source(file = "R/barebones_contmle.R")
set.seed(0)

# simulation parameters -----------------------------------------------------------------------

B <- 200
n_cores <- 4
registerDoParallel(n_cores)
registerDoRNG(0)
target_times_cont <- 1:4 * 360
target_events <- 1:3
sim_results <- vector("list", B)
b <- 1

# true risks ----------------------------------------------------------------------------------

psi0 <- NULL

# generate data -------------------------------------------------------------------------------
sim_data <- foreach(i = 1:B) %dopar% {
    return(simulate_data(n = 1e3, base_data = base_data) %>%
               mutate(ARM = as.numeric(ARM)))
}


# estimation ----------------------------------------------------------------------------------


while (b <= B) {
    indices <- b:min(B, b + 5*(n_cores))
    cat("Simulations", b, "-", tail(indices, 1), "\n")
    b <- b + 5*(n_cores)
    
    sim_results[indices] <- foreach(i = indices) %dopar% {
        results <- list()
        
        obs <- sim_data[[1]] %>% 
            mutate(TIME_1mo = ceiling(TIME / 30), 
                   TIME_3mo = ceiling(TIME / (30 * 3)), 
                   TIME_6mo = ceiling(TIME / (30 * 6)))
        
        ftimes_cont <- obs$TIME
        
        adjust_vars <- dplyr::select(obs, -c(grep("TIME", colnames(obs), value = T), 
                                             "EVENT", "ARM", "id"))
        
        
        # survtmle ----------------------------------------------------------------------------
        
        screeners <- "All"
        sl_glmnets <- create.Learner("SL.glmnet", tune = list("alpha" = c(0.5, 1)), 
                                     name_prefix = "glmnet", detailed_names = T)
        
        sl_lib_g <- expand.grid(c("SL.glm"), screeners)
        sl_lib_g <- lapply(1:nrow(sl_lib_g), function(i) as.character(unlist(sl_lib_g[i, ])))
        
        sl_lib_censor <- expand.grid(c("SL.glm", sl_glmnets$names), screeners)
        sl_lib_censor <- lapply(1:nrow(sl_lib_censor),
                                function(i) as.character(unlist(sl_lib_censor[i, ])))
        
        sl_lib_failure <- expand.grid(c("SL.glm", sl_glmnets$names), screeners)
        sl_lib_failure <- lapply(1:nrow(sl_lib_failure),
                                 function(i) as.character(unlist(sl_lib_failure[i, ])))
        
        for (timescale in c(1, 3, 6)) {
            target_times <- target_times_cont / (30 * timescale)
            ftimes <- ceiling(obs$TIME / (30 * timescale))
            
            sl_fit <- my_init_sl_fit(
                T_tilde = ftimes,
                Delta = as.numeric(obs$EVENT),
                A = as.numeric(obs$ARM),
                W = adjust_vars,
                t_max = max(target_times),
                sl_failure = sl_lib_failure,
                sl_censoring = sl_lib_censor,
                sl_treatment = "SL.glm",
                cv.Control = list(V = 10))
            
            haz_sl <- list(sl_fit$density_failure_1$clone(),
                           sl_fit$density_failure_0$clone())
            haz_sl[[1]]$haz2surv()
            haz_sl[[2]]$haz2surv()
            names(haz_sl) <- c("A = 1", "A = 0")
            
            SL_ftime <- sl_fit$models$Y
            sl_G_dC <- sl_fit$G_dC
            # glm_trt <- paste0(colnames(adjust_vars), collapse = " + ")
            rm(sl_fit); gc()
            
            tmle_sl <- surv_tmle(ftime = ftimes,
                                 ftype = obs$EVENT,
                                 targets = target_times,
                                 trt = obs$ARM,
                                 t0 = max(target_times), adjustVars = adjust_vars,
                                 SL.ftime = SL_ftime, SL.ctime = sl_G_dC,
                                 SL.trt = sl_lib_g, # glm.trt = glm_trt,
                                 returnIC = T, returnModels = T,
                                 ftypeOfInterest = target_events, 
                                 trtOfInterest = c(1, 0),
                                 maxIter = 10, method = "hazard")
            
            tmle_sl_out <- suppressWarnings(
                cbind(A = rep(0:1, times = length(target_events)), 
                      J = rep(target_events, each = 2), 
                      tmle_sl$est) %>% as.data.table() %>% 
                    melt(., id.var = c("A", "J"), variable = "time") %>% 
                    .[, time := rep(target_times_cont, 
                                    each = length(target_events) * 2)] %>% 
                    as_tibble() %>%
                    cbind(Estimator = "survtmle", 
                          timescale = paste0(timescale, "mo"), 
                          .)
            )
            setnames(tmle_sl_out, "value", "RiskEst")
            tmle_sl_out <- suppressWarnings(
                cbind("A" = rep(0:1, each = length(target_times_cont) * length(target_events)), 
                      "J" = rep(target_events, each = length(target_times_cont), times = 2), 
                      "time" = rep(target_times_cont, times = length(target_events) * 2),
                      "se" = sqrt(diag(tmle_sl$var))) %>%
                    as.data.frame() %>% full_join(tmle_sl_out, .)
            )
            if (is.null(results$estimates)) {
                results$estimates <- tmle_sl_out
            } else {
                results$estimates <- rbind(results$estimates, tmle_sl_out)
            }
        }
        
        
        # # tmle3 -------------------------------------------------------------------------------
        # k_grid <- 1:max(obs$TIME)
        # 
        # all_times <- lapply(k_grid, function(t_current) {
        #     df_time <- copy(obs)
        #     
        #     df_time$N <- as.numeric(t_current == obs$TIME & obs$EVENT == 1)
        #     df_time$A_c <- as.numeric(t_current == obs$TIME & obs$EVENT == 0)
        #     df_time$pre_failure <- as.numeric(t_current <= obs$TIME)
        #     df_time$t <- t_current
        #     df_time$X <- 1:nrow(df_time)
        #     
        #     return(df_time)
        # })
        # 
        # df_long <- rbindlist(all_times)
        # 
        # node_list <- list(
        #     W = c("GEOGR1", "SEX", "AGE", "CREATBL", "HBA1CBL", "MIFL",
        #           "SMOKER", "STROKSFL", "BMIBL", "ETHNIC", "EGFMDRBC"),
        #     A = "ARM",
        #     T_tilde = "TIME",
        #     Delta = "EVENT",
        #     time = "t",
        #     N = "N",
        #     A_c = "A_c",
        #     id = "X",
        #     pre_failure = "pre_failure"
        # )
        # 
        # learners <- list(
        #     glmnet = make_learner(Lrnr_glmnet),
        #     glm = make_learner(Lrnr_glm),
        #     gam = make_learner(Lrnr_gam),
        #     rf = make_learner(Lrnr_ranger)
        # )
        # lrnr_glm <- make_learner(Lrnr_glm)
        # 
        # sl_A <- Lrnr_sl$new(lrnr_glm)
        # sl_Y <- Lrnr_sl$new(learners)
        # 
        # learner_list <- list(A = sl_A, N = sl_Y, A_c = sl_Y)
        # var_types <- list(T_tilde = Variable_Type$new("continuous"),
        #                   t = Variable_Type$new("continuous"),
        #                   Delta = Variable_Type$new("binomial"))
        # 
        # survival_spec1 <- tmle_survival(
        #     treatment_level = 1, control_level = 0,
        #     variable_types = var_types, target_times = target_times
        # )
        # survival_spec0 <- tmle_survival(
        #     treatment_level = 0, control_level = 1,
        #     variable_types = var_types, target_times = target_times
        # )
        # 
        # tmle_task1 <- survival_spec1$make_tmle_task(df_long, node_list)
        # tmle_task0 <- survival_spec0$make_tmle_task(df_long, node_list)
        # 
        # initial_likelihood1 <- survival_spec1$make_initial_likelihood(tmle_task1, learner_list)
        # initial_likelihood0 <- survival_spec0$make_initial_likelihood(tmle_task0, learner_list)
        # 
        # up1 <- tmle3_Update_survival$new(
        #     maxit = 25,
        #     cvtmle = TRUE,
        #     convergence_type = "scaled_var",
        #     delta_epsilon = 1e-2,
        #     fit_method = "l2",
        #     use_best = TRUE,
        #     verbose = TRUE
        # )
        # up0 <- tmle3_Update_survival$new(
        #     maxit = 25,
        #     cvtmle = TRUE,
        #     convergence_type = "scaled_var",
        #     delta_epsilon = 1e-2,
        #     fit_method = "l2",
        #     use_best = TRUE,
        #     verbose = TRUE
        # )
        # 
        # targeted_likelihood1 <- Targeted_Likelihood$new(initial_likelihood1, updater = up1)
        # targeted_likelihood0 <- Targeted_Likelihood$new(initial_likelihood0, updater = up0)
        # 
        # tmle_params1 <- survival_spec1$make_params(tmle_task1, targeted_likelihood1)
        # tmle_params0 <- survival_spec0$make_params(tmle_task0, targeted_likelihood0)
        # 
        # # max(abs(colMeans(tmle_params1[[1]]$estimates(tmle_task, "validation")$IC[, 1:10])))
        # # debugonce(tmle_params[[1]]$estimates)
        # 
        # tmle_fit_manual1 <- fit_tmle3(
        #     tmle_task1, targeted_likelihood1, tmle_params1,
        #     targeted_likelihood1$updater
        # )
        # 
        # tmle_fit_manual0 <- fit_tmle3(
        #     tmle_task0, targeted_likelihood0, tmle_params0,
        #     targeted_likelihood0$updater
        # )
        # 
        # # conv <- apply(abs(do.call(rbind,up$EDs)),1,max)
        # results$estimates <-
        #     data.frame("Estimator" = "tmle3",
        #                "t" = target_times,
        #                "s0" = tmle_fit_manual0$estimates[[1]]$psi[target_times],
        #                "s1" = tmle_fit_manual1$estimates[[1]]$psi[target_times],
        #                "se0" = sqrt(diag(var(tmle_fit_manual0$estimates[[1]]$IC)) /
        #                                 nrow(obs))[target_times],
        #                "se1" = sqrt(diag(var(tmle_fit_manual1$estimates[[1]]$IC)) /
        #                                 nrow(obs))[target_times]) %>%
        #     bind_rows(results$estimates, .)
        # 
        # 
        # contmle discretized scale ------------------------------------------------------------
        
        for (timescale in c("cont", 1, 3, 6)) {
            if (timescale == "cont") {
                obs$TIME <- ftimes_cont
                target_times <- target_times_cont
            } else {
                obs$TIME <- ftimes_cont / (30 * timescale)
                target_times <- target_times_cont / (30 * timescale)
            }
            est <- lapply(1:0, function(a) {
                contmle(obs, 
                        target = target_events, 
                        iterative = F, one.step = T, 
                        treat.effect = as.character(a), 
                        tau = target_times, 
                        estimation = list("cause1" = list(fit = "cox",
                                                          model = Surv(TIME, EVENT == 1) ~
                                                              ARM + SEX + AGE + BMIBL +
                                                              SMOKER + STROKSFL + MIFL),
                                          "cause2" = list(fit = "cox",
                                                          model = Surv(TIME, EVENT == 2) ~
                                                              ARM + SEX + AGE + BMIBL +
                                                              SMOKER + STROKSFL + MIFL), 
                                          "cause3" = list(fit = "cox",
                                                          model = Surv(TIME, EVENT == 3) ~
                                                              ARM + SEX + AGE + BMIBL +
                                                              SMOKER + STROKSFL + MIFL),
                                          "cens" = list(fit = "cox",
                                                        model = Surv(TIME, EVENT == 0) ~
                                                            ARM + SEX + AGE + BMIBL +
                                                            SMOKER + STROKSFL + MIFL)
                        ),
                        treat.model = ARM ~ SEX + AGE + BMIBL + SMOKER + STROKSFL + MIFL, 
                        sl.models = list(mod1 = list(Surv(TIME, EVENT == 1) ~ ARM),
                                         mod2 = list(Surv(TIME, EVENT == 1) ~
                                                         ARM + SEX + AGE + BMIBL),
                                         mod3 = list(Surv(TIME, EVENT == 1) ~
                                                         ARM + SEX + AGE + BMIBL +
                                                         SMOKER + STROKSFL + MIFL))
                )
            })
            timescale_contmle_out <- ifelse(timescale == "cont", "cont", paste0(timescale, "mo"))
            contmle_out <- as.data.frame(t(est[[1]]$tmle)) %>%
                cbind("Estimator" = "contmle", 
                      timescale = timescale_contmle_out, 
                      "t" = target_times_cont, .) %>%
                rename_all(~c("Estimator", "t", "s1", "se1")) %>%
                cbind(t(est[[2]]$tmle)) %>%
                rename_all(~c("Estimator", "t", "s1", "se1", "s0", "se0")) %>%
                mutate_at(c("s1", "s0"), ~ 1 - .) %>%
                dplyr::select(colnames(results$estimates))
            
            if(is.null(results$estimates)) {
                results$estimates <- contmle_out
            } else {
                results$estimates <- bind_rows(results$estimates, contmle_out)
            }
            
        }
    }
}

sim_out <- sim_results[!sapply(sim_results, is.null)]
saveRDS(sim_out, file = "R/sim_est.RDS")
sim_out <- lapply(1:length(sim_out), function(i) cbind("iter" = i, sim_out[[i]]$estimates)) %>%
    bind_rows() %>%
    mutate(RD = (1 - s1) - (1 - s0),
           RR = (1 - s1) / (1 - s0),
           SR = s1 / s0,
           s0_se = se0, s1_se = se1,
           RR_se = sqrt(se1^2 / (1 - s0)^2 + se0^2 * ((1 - s1) / (1 - s0)^2)^2),
           RD_se = sqrt(se1^2 + se0^2),
           SR_se = sqrt(se1^2 / s0^2 + se0^2 * s1^2 / s0^4),
           Estimator = as.character(Estimator)) %>% select(-c(se0, se1)) %>%
    pivot_longer(cols = c(`s0`, `s1`, `RD`, `RR`, `SR`), names_to = "Estimand",
                 values_to = "Estimate") %>%
    pivot_longer(cols = contains("se"), names_to = c("se_est", "tmp"),
                 names_sep = "_", values_to = "se") %>% dplyr::select(-`tmp`) %>%
    filter(se_est == Estimand) %>% dplyr::select(-se_est)

sim_tbl <- left_join(sim_out, psi0) %>%
    mutate(`t` = as_factor(`t`), Estimand = as_factor(Estimand),
           bias = Truth - Estimate,
           MSE = bias^2,
           coverage = (Estimate + 1.96*se > Truth) & (Estimate - 1.96*se < Truth)) %>%
    dplyr::select(-iter) %>% group_by(Estimator, t, Estimand) %>%
    mutate(upper = quantile(Estimate, 0.975), lower = quantile(Estimate, 0.025)) %>%
    summarise_all(mean) %>% ungroup()

sim_tbl %>% ggplot(aes(x = `t`)) + facet_wrap(~Estimand, scales = "free", nrow = 5) +
    geom_errorbar(aes(ymin = lower, ymax = upper, colour = Estimator),
                  width = .5, position = position_dodge(.5)) +
    geom_point(aes(y = Estimate, colour = Estimator), position = position_dodge(0.5)) +
    geom_label(aes(y = upper, label = paste0(round(coverage, 2)*100, "%"),
                   group = Estimator), position = position_dodge(0.5), vjust = -.2) +
    geom_segment(aes(y = Truth, yend = Truth,
                     x=as.numeric(`t`)-.3, xend = as.numeric(`t`)+.3),
                 data = distinct(dplyr::select(sim_tbl, `t`, Estimand, Truth))) +
    theme_minimal()

ggsave(filename = "output/contmle-tmle3-survtmle", device = "png",
       width = 12, height = 8, units = "in")
