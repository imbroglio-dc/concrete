
# setup ---------------------------------------------------------------------------------------

library(tidyverse); library(sl3); library(tmle3); library(survtmle); library(MOSS)
library(data.table); library(foreach); library(doParallel); library(doRNG); library(parallel)
library(survival); library(zoo)

source(file = "R/functions/my_sl_functions.R")
source(file = "R/functions/my_MOSS_hazard_methods.R")
source(file = "R/functions/my_survtmle_functions.R")
source(file = "R/functions/sim_functions.R")
source(file = "R/functions/contmle.R")
source(file = "R/barebones_contmle.R")


# simulation parameters -----------------------------------------------------------------------

B <- 40
n_cores <- 8
cl <- makeForkCluster(n_cores)
registerDoParallel(cl)
target_times_cont <- 1:4 * 360
target_events <- 1:3
sim_results <- vector("list", B)
b <- 1


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
    true_risks <- as.data.table(read.csv("./data/true_risks.csv"))
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
        rm(obs); gc()
        outcomes[, "1>2" := T1 > T2]
        outcomes[, "2>3" := T2 > T3]
        outcomes[, "3>1" := T3 > T1]
        outcomes[, "J" := `3>1`*(1 - `1>2`) + 2*`1>2`*(1 - `2>3`) + 3*`2>3`*(1 - `3>1`)]
        outcomes[, "T" := T1*(J == 1) + T2*(J == 2) + T3*(J == 3)]
        outcomes <- outcomes[, c("J", "T")]
        risk <- as.data.table(foreach(t = interval,
                                      .combine = rbind,
                                      .inorder = T) %dopar% {
                                          return(tabulate(outcomes[["J"]][outcomes[["T"]] <= t]))
                                      })
        if (!exists("true_risks")) {
            true_risks <- risk[, lapply(.SD, function(x) x / length(A))]
            setnames(true_risks, 1:3, do.call(paste0, expand.grid("F.j", 1:3, ".a", a)))
        } else {
            risk <- risk[, lapply(.SD, function(x) x / length(A))]
            setnames(risk, 1:3, do.call(paste0, expand.grid("F.j", 1:3, ".a", a)))
            true_risks <- cbind(true_risks, risk)
        }
        rm(risk); rm(outcomes); rm(A); gc()
    }
    true_risks <- 
        rbind(cbind(J = 1, "time" = interval,
                    true_risks[, mget(grep("j1", colnames(true_risks), value = T))]), 
              cbind(J = 2,  "time" = interval,
                    true_risks[, mget(grep("j2", colnames(true_risks), value = T))]), 
              cbind(J = 3,  "time" = interval,
                    true_risks[, mget(grep("j3", colnames(true_risks), value = T))]), 
              use.names = FALSE)
    setnames(true_risks, 2:3, c("F.a1", "F.a0"))
    
    write_csv(true_risks, "data/true_risks.csv")
}


melt(true_risks, id.vars = c("J", "time"), variable.name = "Parameter") %>% 
    filter(Parameter %in% c("F.a1", "F.a0")) %>% 
    mutate(J = as.character(J)) %>% ggplot() + 
    geom_line(aes(x = time, y = value, colour = Parameter, linetype = J)) + 
    theme_minimal()


# true risks ----------------------------------------------------------------------------------

psi0 <- true_risks[time %in% target_times_cont,]

# generate data -------------------------------------------------------------------------------

registerDoRNG(123456789)
sim_data <- foreach(i = 1:B) %dopar% {
    return(simulate_data(n = 1e3, base_data = base_data) %>%
               mutate(ARM = as.numeric(ARM)))
}


# estimation ----------------------------------------------------------------------------------

screeners <- "All"
sl_glmnets <- create.Learner("SL.glmnet", tune = list("alpha" = c(1)), 
                             name_prefix = "glmnet", detailed_names = T)

sl_lib_g <- expand.grid(c("SL.glm"), screeners)
sl_lib_g <- lapply(1:nrow(sl_lib_g), function(i) as.character(unlist(sl_lib_g[i, ])))

sl_lib_censor <- expand.grid(c("SL.glm", "SL.glmnet"), screeners)
sl_lib_censor <- lapply(1:nrow(sl_lib_censor),
                        function(i) as.character(unlist(sl_lib_censor[i, ])))

sl_lib_failure <- expand.grid(c("SL.glm", "SL.glmnet"), screeners)
sl_lib_failure <- lapply(1:nrow(sl_lib_failure),
                         function(i) as.character(unlist(sl_lib_failure[i, ])))

while (b <= B) {
    indices <- b:min(B, b + 1*(n_cores))
    cat("Simulations", b, "-", tail(indices, 1), "\n")
    b <- b + 1*(n_cores)
    
    sim_results[indices] <- foreach(i = indices) %dopar% {
        obs <- sim_data[[i]]
        ftimes_cont <- obs$TIME
        adjust_vars <- dplyr::select(obs, -c(grep("TIME", colnames(obs), value = T), 
                                             "EVENT", "ARM", "id"))
        
        
        # survtmle ----------------------------------------------------------------------------
        
        for (timescale in c(1, 3, 6)) {
            target_times <- target_times_cont / (30 * timescale)
            ftimes <- ceiling(obs$TIME / (30 * timescale))
            
            sl_fit <- suppressMessages(suppressWarnings(
                my_init_sl_fit(
                    T_tilde = ftimes,
                    Delta = as.numeric(obs$EVENT),
                    A = as.numeric(obs$ARM),
                    W = adjust_vars,
                    t_max = max(target_times),
                    sl_failure = sl_lib_failure,
                    sl_censoring = sl_lib_censor,
                    sl_treatment = "SL.glm",
                    cv.Control = list(V = 10))
            ))
            
            haz_sl <- list(sl_fit$density_failure_1$clone(),
                           sl_fit$density_failure_0$clone())
            haz_sl[[1]]$haz2surv()
            haz_sl[[2]]$haz2surv()
            names(haz_sl) <- c("A = 1", "A = 0")
            
            SL_ftime <- sl_fit$models$Y
            sl_G_dC <- sl_fit$G_dC
            # glm_trt <- paste0(colnames(adjust_vars), collapse = " + ")
            rm(sl_fit); gc()
            
            survtmle_fit <- suppressMessages(suppressWarnings(
                surv_tmle(ftime = ftimes,
                          ftype = obs$EVENT,
                          targets = target_times,
                          trt = obs$ARM, 
                          t0 = max(target_times), adjustVars = adjust_vars,
                          SL.ftime = SL_ftime, SL.ctime = sl_G_dC,
                          SL.trt = sl_lib_g, # glm.trt = glm_trt,
                          returnIC = T, returnModels = T,
                          ftypeOfInterest = target_events, 
                          trtOfInterest = c(1, 0), verbose = FALSE,
                          maxIter = 10, method = "hazard")
            ))
            
            survtmle_out <- suppressWarnings(
                cbind(estimand = rep(paste0("F", target_events), each = 2), 
                      A = rep(0:1, times = length(target_events)), 
                      as.data.table(survtmle_fit$est)) %>% 
                    melt(., id.var = c("A", "estimand"), variable = "time") %>% 
                    .[, time := rep(target_times_cont, 
                                    each = length(target_events) * 2)] %>% 
                    cbind(estimator = "survtmle", 
                          timescale = paste0(timescale, "mo"), 
                          .)
            )
            
            survtmle_out <- 
                rbind(cbind(estimator = "glm_gcomp", as.data.table(survtmle_fit$init_est)), 
                      cbind(estimator = "survtmle", as.data.table(survtmle_fit$est))) %>% 
                cbind(estimand = rep(paste0("F", target_events), each = 2), 
                      A = rep(0:1, times = length(target_events)), 
                      .) %>% 
                melt(., id.var = c("estimator", "estimand", "A"), variable = "time") %>% 
                .[, time := rep(target_times_cont, 
                                each = length(target_events) * 4)] %>% 
                cbind(timescale = paste0(timescale, "mo"), 
                      .)
            
            setnames(survtmle_out, "value", "estimate")
            setcolorder(survtmle_out, c('estimator', 'timescale', 'estimand', 'A', 'time'))
            
            if (!exists("results")) {
                results <- copy(survtmle_out)
            } else {
                results <- rbind(results, copy(survtmle_out))
            }
        }
        
        # contmle discretized scale ------------------------------------------------------------
        
        
        for (timescale in c("cont", 3, 6)) {
            if (timescale == "cont") {
                obs$TIME <- ftimes_cont
                target_times <- target_times_cont
            } else {
                obs$TIME <- ftimes_cont / (30 * as.numeric(timescale))
                target_times <- target_times_cont / (30 * as.numeric(timescale))
            }
            
            logreg <- make_learner(Lrnr_glm)
            lasso <- make_learner(Lrnr_glmnet) # alpha default is 1
            ridge <- Lrnr_glmnet$new(alpha = 0)
            e_net <- make_learner(Lrnr_glmnet, alpha = 0.5)
            a_lrnrs <- make_learner(Stack, logreg, lasso, ridge, e_net)
            
            models <- list("A" = a_lrnrs, 
                           "0" = list(mod1 = Surv(TIME, EVENT == 0) ~ ARM,
                                      mod2 = Surv(TIME, EVENT == 0) ~ ARM + AGE,
                                      mod3 = Surv(TIME, EVENT == 0) ~ ARM + AGE + SMOKER + STROKSFL,
                                      mod4 = Surv(TIME, EVENT == 0) ~ .), 
                           "1" = list(mod1 = Surv(TIME, EVENT == 1) ~ ARM,
                                      mod2 = Surv(TIME, EVENT == 1) ~ ARM + SMOKER + BMIBL,
                                      # mod3 = Surv(TIME, EVENT == 1) ~ ARM*SMOKER + I(BMIBL>30)*ARM, 
                                      mod4 = Surv(TIME, EVENT == 1) ~ .), 
                           "2" = list(mod1 = Surv(TIME, EVENT == 2) ~ ARM,
                                      mod2 = Surv(TIME, EVENT == 2) ~ ARM + STROKSFL + MIFL,
                                      # mod3 = Surv(TIME, EVENT == 2) ~ ARM*STROKSFL + ARM*MIFL + MIFL:STROKSFL, 
                                      mod4 = Surv(TIME, EVENT == 2) ~ .), 
                           "3" = list(mod1 = Surv(TIME, EVENT == 3) ~ ARM,
                                      mod2 = Surv(TIME, EVENT == 3) ~ ARM + SMOKER + BMIBL,
                                      mod3 = Surv(TIME, EVENT == 3) ~ ARM + SMOKER + BMIBL + STROKSFL + MIFL,
                                      mod4 = Surv(TIME, EVENT == 3) ~ .))
            timescale_contmle_out <- ifelse(timescale == "cont", "cont", 
                                            paste0(timescale, "mo"))
            
            contmle_est <- concr_tmle(obs, target_times, target_events, models)
            tmle_out <- rbind(cbind('estimator' = 'contmle', 
                                    as.data.table(contmle_est$estimates$tmle[, -"S"])), 
                              cbind('estimator' = 'cox_gcomp', 
                                    as.data.table(contmle_est$estimates$`g-comp`[, -"S"]))) %>% 
                melt(., id.vars = c("estimator", "A", "time")) %>%
                .[, timescale := timescale_contmle_out]
            setnames(tmle_out, c("variable", "value"), c("estimand", "estimate"))
            setcolorder(tmle_out, c('estimator', 'timescale', 'estimand', 'A', 'time'))
            
            if (!exists("results")) {
                results <- copy(tmle_out)
            } else {
                results <- rbind(results, copy(tmle_out))
            }
        }
        
        # return ------------------------------------------------------------------
        wide_results <- dcast(results, ... ~ A, value.var = "estimate") %>%
            .[, RR := `1` / `0`] %>%
            .[, RD := `1` - `0`] %>%
            .[, SR := (1 - `1`) / (1 - `0`)]
        setnames(wide_results, c("estimand", "0", "1"), c("J", "F.a0", "F.a1"))
        
        # sim_est <- rbind(sim_est, copy(results))
        
        return(wide_results)
    }
}

sim_results

sim_out <- bind_rows(sim_results[!sapply(sim_results, is.null)])
write_csv(sim_out, file = "output/survtmle_contmle_40.csv")



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
