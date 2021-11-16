
# 0. Preparation ------------------------------------------------------------------------------

## 0.1 Packages and Base Data -----------------------------------------------------------------

library(tidyverse); library(readxl); library(skimr); library(data.table); library(here)
library(doParallel); library(foreach); library(survival); library(zoo)
setwd("/Shared/Projects/ConCR-TMLE/")
i_am("./R/contmle-competing-risks-simulation.R")
source("./R/contmle.R")

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

# skimr::skim(base_data)
# base_data covariates: 2 factor, 4 logical, 5 numeric
#   ARM, TIME, and EVENT will be replaced


# 0.2 Variables -------------------------------------------------------------------------------

interval <- 1:2e3
n <- 1e3

# 0.3 simulation functions ------------------------------------------------

source("R/functions/sim_functions.R")


# 3. Estimation -------------------------------------------------------------------------------

# Estimation targets
tau <- 720
target <- 1:3
B <- 800
n_cores <- 8
registerDoParallel(n_cores)

# 3.1 True Psi --------------------------------------------------------------------------------

if (file.exists("./data/true_risks.RDS")) {
    true_risks <- readRDS("./data/true_risks.RDS")
} else {
    for (a in 1:0) { # for binary treatment only
        obs <- as.data.table(bind_rows(lapply(1:5000, function(b) base_data)))
        true_risks <- list("A=1" = NULL, "A=0" = NULL)
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
    saveRDS(true_risks, "./data/true_risks.RDS")
}

# plot survival curves
lapply(true_risks, function(r) rename_all(r, ~c("J=1", "J=2", "J=3"))) %>% 
    bind_rows() %>% mutate(`A` = rep(1:0, each = nrow(true_risks[[1]])), 
                           `t` = rep(1:nrow(true_risks[[1]]), times = 2)) %>% 
    pivot_longer(cols = c(`J=1`, `J=2`, `J=3`), names_to = "event", 
                 values_to = "risk") %>% 
    ggplot(aes(x = `t`, y = risk, 
               linetype = as.character(A), colour = event)) + 
    geom_line()+ theme_minimal()


# 3.2 contmle Estimation ----------------------------------------------------------------------

if (file.exists("./output/contmle_estimates.RDS")) {
    estimates <- readRDS("./output/contmle_estimates.RDS")
} else {
    estimates <- foreach(i=1:B,
                         .combine = rbind,
                         .packages = c("data.table", "tidyverse", "survival", "zoo")) %dopar% {
                             
                             obs <- simulate_data(n = n, base_data = base_data)
                             est <- contmle(obs, #-- dataset
                                            target = target, #-- target competing events
                                            iterative = FALSE, #-- one-step tmle to target \simultaneously
                                            treat.effect = "ate", #-- target the ate directly
                                            tau = tau, #-- time-point of interest
                                            estimation = list("cause1" = list(fit = "cox",
                                                                              model = Surv(TIME, EVENT == 1) ~ 
                                                                                  ARM + SEX + AGE + BMIBL + 
                                                                                  SMOKER + STROKSFL + MIFL),
                                                              "cens" = list(fit = "cox",
                                                                            model = Surv(TIME, EVENT == 0) ~ 
                                                                                ARM + SEX + AGE + BMIBL + 
                                                                                SMOKER + STROKSFL + MIFL), 
                                                              "cause2" = list(fit = "cox",
                                                                              model = Surv(TIME, EVENT == 2) ~ 
                                                                                  ARM + SEX + AGE + BMIBL + 
                                                                                  SMOKER + STROKSFL + MIFL), 
                                                              "cause3" = list(fit = "cox",
                                                                              model = Surv(TIME, EVENT == 3) ~ 
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
                             estimates <- bind_rows(unlist(bind_cols(est$init)[1, ]),
                                                    unlist(bind_cols(est$tmle)[1, ]))
                             estimates <- as.data.table(cbind("run" = i, 
                                                              "Estimator" = c("init", "TMLE"), 
                                                              estimates))
                             estimates
                         }
    saveRDS(estimates, "./output/contmle_estimates.RDS")
}


# 4. Evaluation / Visualization ---------------------------------------------------------------


Psi0 <- cbind("time" = tau, (true_risks[["A=1"]] - true_risks[["A=0"]])[tau, ])
estimates %>% group_by(Estimator) %>% summarise_all(list(mean = mean, se = var)) %>% 
    mutate_at(vars(ends_with("se")), sqrt)



