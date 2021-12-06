### compare-tmle-methods-leader.R --- 
#----------------------------------------------------------------------
## Author: David Chen & Thomas Alexander Gerds
## Created: Dec  6 2021 (16:59) 
## Version: 
## Last-Updated: Dec  6 2021 (17:42) 
##           By: Thomas Alexander Gerds
##     Update #: 3
#----------------------------------------------------------------------
## 
### Commentary:
##
## Here we discuss the input methods and the results of
## tmle3, survtmle, contmle
## in a head-to-head comparison
##
##
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

# setup ---------------------------------------------------------------------------------------

library(readxl)
library(tidyverse); library(sl3); library(tmle3); library(survtmle); library(MOSS)
library(data.table); library(foreach); library(doParallel); library(doRNG); library(parallel)
library(glmnet)
library(SuperLearner)

## Add your pwd to suit local working directory structure -------------------------------------

try(setwd("/Shared/Projects/ConCR-TMLE/"),silent=TRUE)
try(setwd("~/research/SoftWare/test-tmle-surival/ConCR-TMLE/"),silent=TRUE)
source("R/functions/sim_functions.R")
source(file = "R/functions/my_sl_functions.R")
source(file = "R/functions/my_MOSS_hazard_methods.R")
source(file = "R/functions/my_survtmle_functions.R")
source(file = "R/contmle.R")

# ---------------------------------------------------------------------------------------------
# get base_data and discretize time 
# ---------------------------------------------------------------------------------------------
set.seed(0)
base_data <- read_excel("./data/test_leader.xlsx") %>%
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

obs <- base_data
obs_cont <- obs
obs$TIME <- ceiling(obs$TIME / 30 / 3) # convert to 3 month intervals
target_times <- (1:5) * 4 # years 1:5 for time discretized to 3 month intervals
target_times_cont <- target_times * 30.45 * 3
adjust_vars <- dplyr::select(obs, -c("TIME", "EVENT", "ARM"))

# ---------------------------------------------------------------------------------------------
# 1. survtmle
# ---------------------------------------------------------------------------------------------

screeners <- "All"
sl_lib_g <- expand.grid(c("SL.glm"), screeners)
sl_lib_g <- lapply(1:nrow(sl_lib_g), function(i) as.character(unlist(sl_lib_g[i, ])))

## sl_lib_censor <- expand.grid(c("SL.glm", "SL.glmnet", "SL.gam", "SL.ranger"), screeners)
sl_lib_censor <- expand.grid(c("SL.glm"), screeners)
sl_lib_censor <- lapply(1:nrow(sl_lib_censor), function(i) as.character(unlist(sl_lib_censor[i, ])))

## sl_lib_failure <- expand.grid(c("SL.glm", "SL.glmnet", "SL.gam", "SL.ranger"), screeners)
sl_lib_failure <- expand.grid(c("SL.glm"), screeners)
sl_lib_failure <- lapply(1:nrow(sl_lib_failure),function(i) as.character(unlist(sl_lib_failure[i, ])))

sl_fit <- my_init_sl_fit(
    T_tilde = obs$TIME,
    Delta = as.numeric(obs$EVENT),
    A = as.numeric(obs$ARM),
    W = adjust_vars,
    t_max = max(target_times),
    sl_failure = sl_lib_failure,
    sl_censoring = sl_lib_censor,
    sl_treatment = "SL.glm",
    cv.Control = list(V = 10))

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

        tmle_sl <- surv_tmle(ftime = obs$TIME,
                             ftype = obs$EVENT,
                             targets = target_times,
                             trt = obs$ARM,
                             t0 = max(target_times), adjustVars = adjust_vars,
                             SL.ftime = SL_ftime, SL.ctime = sl_G_dC,
                             SL.trt = sl_lib_g, # glm.trt = glm_trt,
                             returnIC = T, returnModels = T,
                             ftypeOfInterest = 1, trtOfInterest = c(1, 0),
                             maxIter = 25, method = "hazard")

        tmle_sl_out <- suppressWarnings(
            t(1 - tmle_sl$est) %>% unname() %>% as_tibble() %>%
                cbind(Estimator = "survtmle", t = target_times, .) %>%
                rename("s0" = V1, "s1" = V2))
        results$estimates <- suppressWarnings(
            matrix(sqrt(diag(tmle_sl$var)),
                   nrow = length(target_times), byrow = F,
                   dimnames = list(NULL, c("se0", "se1"))) %>%
                as.data.frame() %>% cbind(tmle_sl_out, .))
#------ -------------------------------------------------------------------------------
# tmle3
#------ -------------------------------------------------------------------------------
k_grid <- 1:max(obs$TIME)

        all_times <- lapply(k_grid, function(t_current) {
            df_time <- copy(obs)

            df_time$N <- as.numeric(t_current == obs$TIME & obs$EVENT == 1)
            df_time$A_c <- as.numeric(t_current == obs$TIME & obs$EVENT == 0)
            df_time$pre_failure <- as.numeric(t_current <= obs$TIME)
            df_time$t <- t_current
            df_time$X <- 1:nrow(df_time)

            return(df_time)
        })

        df_long <- rbindlist(all_times)

        node_list <- list(
            W = c("GEOGR1", "SEX", "AGE", "CREATBL", "HBA1CBL", "MIFL",
                  "SMOKER", "STROKSFL", "BMIBL", "ETHNIC", "EGFMDRBC"),
            A = "ARM",
            T_tilde = "TIME",
            Delta = "EVENT",
            time = "t",
            N = "N",
            A_c = "A_c",
            id = "X",
            pre_failure = "pre_failure"
        )

        learners <- list(
            glmnet = make_learner(Lrnr_glmnet),
            glm = make_learner(Lrnr_glm),
            gam = make_learner(Lrnr_gam),
            rf = make_learner(Lrnr_ranger)
        )
        lrnr_glm <- make_learner(Lrnr_glm)

        sl_A <- Lrnr_sl$new(lrnr_glm)
        sl_Y <- Lrnr_sl$new(learners)

        learner_list <- list(A = sl_A, N = sl_Y, A_c = sl_Y)
        var_types <- list(T_tilde = Variable_Type$new("continuous"),
                          t = Variable_Type$new("continuous"),
                          Delta = Variable_Type$new("binomial"))

        survival_spec1 <- tmle_survival(
            treatment_level = 1, control_level = 0,
            variable_types = var_types, target_times = target_times
        )
        survival_spec0 <- tmle_survival(
            treatment_level = 0, control_level = 1,
            variable_types = var_types, target_times = target_times
        )

        tmle_task1 <- survival_spec1$make_tmle_task(df_long, node_list)
        tmle_task0 <- survival_spec0$make_tmle_task(df_long, node_list)

        initial_likelihood1 <- survival_spec1$make_initial_likelihood(tmle_task1, learner_list)
        initial_likelihood0 <- survival_spec0$make_initial_likelihood(tmle_task0, learner_list)

        up1 <- tmle3_Update_survival$new(
            maxit = 25,
            cvtmle = TRUE,
            convergence_type = "scaled_var",
            delta_epsilon = 1e-2,
            fit_method = "l2",
            use_best = TRUE,
            verbose = TRUE
        )
        up0 <- tmle3_Update_survival$new(
            maxit = 25,
            cvtmle = TRUE,
            convergence_type = "scaled_var",
            delta_epsilon = 1e-2,
            fit_method = "l2",
            use_best = TRUE,
            verbose = TRUE
        )

        targeted_likelihood1 <- Targeted_Likelihood$new(initial_likelihood1, updater = up1)
        targeted_likelihood0 <- Targeted_Likelihood$new(initial_likelihood0, updater = up0)

        tmle_params1 <- survival_spec1$make_params(tmle_task1, targeted_likelihood1)
        tmle_params0 <- survival_spec0$make_params(tmle_task0, targeted_likelihood0)

        # max(abs(colMeans(tmle_params1[[1]]$estimates(tmle_task, "validation")$IC[, 1:10])))
        # debugonce(tmle_params[[1]]$estimates)

        tmle_fit_manual1 <- fit_tmle3(
            tmle_task1, targeted_likelihood1, tmle_params1,
            targeted_likelihood1$updater
        )

        tmle_fit_manual0 <- fit_tmle3(
            tmle_task0, targeted_likelihood0, tmle_params0,
            targeted_likelihood0$updater
        )

        # conv <- apply(abs(do.call(rbind,up$EDs)),1,max)
        results$estimates <-
            data.frame("Estimator" = "tmle3",
                       "t" = target_times,
                       "s0" = tmle_fit_manual0$estimates[[1]]$psi[target_times],
                       "s1" = tmle_fit_manual1$estimates[[1]]$psi[target_times],
                       "se0" = sqrt(diag(var(tmle_fit_manual0$estimates[[1]]$IC)) /
                                        nrow(obs))[target_times],
                       "se1" = sqrt(diag(var(tmle_fit_manual1$estimates[[1]]$IC)) /
                                        nrow(obs))[target_times]) %>%
            bind_rows(results$estimates, .)


        # contmle discretized scale ------------------------------------------------------------

        est_disc <- lapply(c("1", "0"), function(a) {
            contmle(obs, #-- dataset
                    target = 1, #-- target competing events
                    iterative = FALSE, #-- one-step tmle to target \simultaneously
                    treat.effect = a,
                    tau = target_times, #-- time-point of interest
                    estimation = list("cause1" = list(fit = "cox",
                                                      model = Surv(TIME, EVENT == 1) ~
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
        results$estimates <- est_disc[[1]]$tmle %>% t() %>% as.data.frame() %>%
            data.frame("Estimator" = "contmle_d", "t" = target_times, .) %>%
            rename_all(~c("Estimator", "t", "s1", "se1")) %>%
            cbind(t(est_disc[[2]]$tmle)) %>%
            rename_all(~c("Estimator", "t", "s1", "se1", "s0", "se0")) %>%
            mutate_at(c("s1", "s0"), ~ 1 - .) %>%
            dplyr::select(colnames(results$estimates)) %>%
            bind_rows(results$estimates, .)


        # contmle - continuous time -----------------------------------------------------------

        est_cont <- lapply(c("1", "0"), function(a) {
            contmle(obs_cont, #-- dataset
                    target = 1, #-- target competing events
                    iterative = FALSE, #-- one-step tmle to target \simultaneously
                    treat.effect = a,
                    tau = target_times_cont, #-- time-point of interest
                    estimation = list("cause1" = list(fit = "cox",
                                                      model = Surv(TIME, EVENT == 1) ~
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

        results$estimates <- est_cont[[1]]$tmle %>% t() %>% as.data.frame() %>%
            data.frame("Estimator" = "contmle_c", "t" = target_times, .) %>%
            rename_all(~c("Estimator", "t", "s1", "se1")) %>%
            cbind(t(est_cont[[2]]$tmle)) %>%
            rename_all(~c("Estimator", "t", "s1", "se1", "s0", "se0")) %>%
            mutate_at(c("s1", "s0"), ~ 1 - .) %>%
            dplyr::select(colnames(results$estimates)) %>%
            bind_rows(results$estimates, .)
        return(results)
    }
}

#----------------------------------------------------------------------
### compare-tmle-methods-leader.R ends here
