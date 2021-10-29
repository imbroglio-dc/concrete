
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

# first two terms in coefs are B and k respectively for weibull hazard = B * k * x^(k - 1)
ltfu_coefs <- c(8e-5, 1, 1.5, 1.1)
eos_coefs <- c("eos_start_time" = 1460, "eos_end_time" = 2000)
t1_coefs <- c(1e-4, 1, 1.2, 1.3, 1.2, 1.2)
t2_coefs <- c(9e-6, 1.3, 1.5, 1.5, 1.2)
t3_coefs <- c(1.8e-5, 1.2)

# Estimation targets
tau <- 720
target <- 1:3
B <- 400
n_cores <- 8
registerDoParallel(n_cores)

# 1. Event Process Models ---------------------------------------------------------------------

return_weibull_outputs <- function(phi, B, k, output = c("h.t", "S.t", "F_inv.u"), t, u) {
    out <- list()
    if ("h.t" %in% output)
        out[["h.t"]] <- phi * B * k * t^(k - 1)
    if ("S.t" %in% output)
        out[["S.t"]] <- exp(-phi * B * t^k)
    if ("F_inv.u" %in% output)
        out[["F_inv.u"]] <- ( -log(1 - u) / (phi * B) )^(1/k)
    return(out)
}

## 1.1 Censoring Model ------------------------------------------------------------------------

# lost-to-followup
ltfu_fn <- function(ARM, AGE, params, output = c("h.t", "S.t", "F_inv.u"), t = NULL, u = NULL) {
    B <- params[1]
    k <- params[2]
    b_A <- log(params[3]) * as.numeric(ARM)
    b_AGE <- log(params[4]) * as.numeric(scale(AGE, scale = max(abs(AGE - mean(AGE))) / 3))
    phi <- exp(b_A + b_AGE)
    
    return(return_weibull_outputs(phi, B, k, output, t, u))
}

# end of study
eos_fn <- function(params, output = c("S.t", "F_inv.u"), t = NULL, u = NULL) {
    out <- list()
    if ("S.t" %in% output)
        out[["S.t"]] <- (t < params[1]) - (t > params[1]) * (t <= params[2]) / diff(params)
    if ("F_inv.u" %in% output)
        out[["F_inv.u"]] <- params[1] + (1 - u) * diff(params)
    return(out)
}


## 1.2 Hazard Models --------------------------------------------------------------------------

### 1.2.1 Event 1 -----------------------------------------------------------------------------

T1_fn <- function(ARM, SMOKER, BMIBL, params,
                  output = c("h.t", "S.t", "F_inv.u"), t = NULL, u = NULL) {
    B <- params[1]
    k <- params[2]
    b1 <- log(params[3]) * as.numeric(SMOKER)
    b2 <- log(params[4]) * as.numeric(ARM == 0) * as.numeric(SMOKER)
    b3 <- log(params[5]) * as.numeric(BMIBL > 30)
    b4 <- log(params[6]) * as.numeric(ARM == 0) * as.numeric(BMIBL > 30)
    phi <- exp(b1 + b2 + b3 + b4)
    
    return(return_weibull_outputs(phi, B, k, output, t, u))
}

### 1.2.2 Event 2 -----------------------------------------------------------------------------
T2_fn <- function(ARM, STROKSFL, MIFL, params,
                  output = c("h.t", "S.t", "F_inv.u"), t = NULL, u = NULL) {
    B <- params[1]
    k <- params[2]
    b1 <- log(params[3]) * as.numeric(ARM == 0) * as.numeric(STROKSFL)
    b2 <- log(params[4]) * as.numeric(ARM == 0) * as.numeric(MIFL)
    b3 <- log(params[5]) * as.numeric(STROKSFL) * as.numeric(MIFL)
    phi <- exp(b1 + b2 + b3)
    
    return(return_weibull_outputs(phi, B, k, output, t, u))
}

### 1.2.3 Event 3 -----------------------------------------------------------------------------

## just a Weibull, shape lambda, scale p
T3_fn <- function(params, output = c("h.t", "S.t", "F_inv.u"), t = NULL, u = NULL) {
    B <- params[1]
    k <- params[2]
    phi <- 1
    
    return(return_weibull_outputs(phi, B, k, output, t, u))
}


# 2. Simulate Data ----------------------------------------------------------------------------

simulate_data <- function(n = 1e3, assign_A = function(W, n) rbinom(n, 1, 0.5), base_data) {
    obs <- dplyr::select(sample_n(base_data, size = n), -ARM, -TIME, -EVENT)
    A <- assign_A(obs, n)
    outcomes <- data.table("C_ltfu" = ltfu_fn(A, obs[["AGE"]], ltfu_coefs,
                                              output = "F_inv.u", u = runif(n, 0, 1))$F_inv.u,
                           "C_eos" = eos_fn(eos_coefs, "F_inv.u", u = runif(n, 0, 1))$F_inv.u,
                           "T1" = T1_fn(A, obs[["SMOKER"]], obs[["BMIBL"]],
                                        t1_coefs, output = "F_inv.u", u = runif(n, 0, 1))$F_inv.u,
                           "T2" = T2_fn(A, obs[["STROKSFL"]], obs[["MIFL"]],
                                        t2_coefs, output = "F_inv.u", u = runif(n, 0, 1))$F_inv.u,
                           "T3" = T3_fn(t3_coefs, output = "F_inv.u", u = runif(n, 0, 1))$F_inv.u)
    
    obs <- cbind(cbind(t(apply(outcomes, 1, function(r)
        c("TIME" = ceiling(min(r)), "EVENT" = max(0, which.min(r) - 2)))),
        "ARM" = A), obs)
    obs <- as.data.table(cbind(obs, "id" = 1:nrow(obs)))
    return(obs)
}

# 3. Estimation -------------------------------------------------------------------------------

# 3.1 True Psi --------------------------------------------------------------------------------

obs <- as.data.table(bind_rows(lapply(1:1200, function(b) base_data)))
true_risks <- list("A=1" = NULL, "A=0" = NULL)
for (a in 1:0) {
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

# plot survival curves
lapply(true_risks, function(r) rename_all(r, ~c("J=1", "J=2", "J=3"))) %>% 
    bind_rows() %>% mutate(`A` = rep(1:0, each = nrow(true_risks[[1]])), 
                           `t` = rep(1:nrow(true_risks[[1]]), times = 2)) %>% 
    pivot_longer(cols = c(`J=1`, `J=2`, `J=3`), names_to = "event", 
                 values_to = "risk") %>% 
    ggplot(aes(x = `t`, y = risk, 
               linetype = as.character(A), colour = event)) + 
    geom_line()+ theme_minimal()

Psi0 <- cbind("time" = tau, (true_risks[["A=1"]] - true_risks[["A=0"]])[tau, ])


# 3.2 contmle Estimation ----------------------------------------------------------------------

estimates <- foreach(i=1:B,
                     .combine = rbind,
                     .packages = c("data.table", "tidyverse", "survival", "zoo")) %dopar% {
                         
                         obs <- simulate_data(n = n, base_data = base_data)
                         est <- contmle(obs, #-- dataset
                                        target = target, #-- go after cause 1 and cause 2 specific risks
                                        iterative = FALSE, #-- use one-step tmle to target F1 and F2 simultaneously
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
                                                         mod2 = list(Surv(TIME, EVENT == 1) ~ ARM + SEX + AGE + BMIBL),
                                                         mod3 = list(Surv(TIME, EVENT == 1) ~ ARM + SEX + AGE + BMIBL + 
                                                                         SMOKER + STROKSFL + MIFL))
                         )
                         estimates <- bind_rows(unlist(bind_cols(est$init)[1, ]),
                                                unlist(bind_cols(est$tmle)[1, ]))
                         estimates <- as.data.table(cbind("run" = i, "Estimator" = c("init", "TMLE"), estimates))
                         estimates
                     }

estimates %>% group_by(Estimator) %>% summarise_all(mean)