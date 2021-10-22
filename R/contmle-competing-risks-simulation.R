
# 0. Preparation ------------------------------------------------------------------------------

## 0.1 Packages and Base Data -----------------------------------------------------------------

library(tidyverse); library(readxl); library(skimr); library(data.table)
library(doParallel); library(foreach); library(survival); library(zoo)

set.seed(0)
base_data <- read_excel("data/test_leader.xlsx") %>%
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


# 1. Event Process Models ---------------------------------------------------------------------

t <- 1:2e3

# first two terms in coefs are B and k respectively for weibull hazard = B * k * x^(k - 1)
ltfu_coefs <- c(1e-4, 1, 1.2, 1.1)
eos_coefs <- c("eos_start_time" = 1460, "eos_end_time" = 2000)
t1_coefs <- c(1e-4, 1, 1.2, 1.3, 1.2, 1.2)
t2_coefs <- c(9e-6, 1.3, 1.5, 1.5, 1.2)
t3_coefs <- c(1.8e-5, 1.2)

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

f1_fn <- function(ARM, SMOKER, BMIBL, params,
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
f2_fn <- function(ARM, STROKSFL, MIFL, params,
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
f3_fn <- function(params, output = c("h.t", "S.t", "F_inv.u"), t = NULL, u = NULL) {
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
                           "T1" = f1_fn(A, obs[["SMOKER"]], obs[["BMIBL"]],
                                        t1_coefs, output = "F_inv.u", u = runif(n, 0, 1))$F_inv.u,
                           "T2" = f2_fn(A, obs[["STROKSFL"]], obs[["MIFL"]],
                                        t2_coefs, output = "F_inv.u", u = runif(n, 0, 1))$F_inv.u,
                           "T3" = f3_fn(t3_coefs, output = "F_inv.u", u = runif(n, 0, 1))$F_inv.u)

    obs <- cbind(cbind(t(apply(outcomes, 1, function(r)
        c("TIME" = ceiling(min(r)), "EVENT" = max(0, which.min(r) - 2)))),
        "ARM" = A), obs)
    obs <- as.data.table(cbind(obs, "id" = 1:nrow(obs)))
    return(obs)
}

# 3. Estimation -------------------------------------------------------------------------------

tau <- 720

# 3.1 True Psi --------------------------------------------------------------------------------

obs <- as.data.table(bind_rows(lapply(1:1000, function(b) base_data)))
Psi0 <- data.frame("F1" = 1:2, "F2" = 1:2, "F3" = 1:2)
n <- nrow(obs)
for (a in 1:0) {
    A <- rep(a, nrow(obs))
    outcomes <- data.table("T1" = f1_fn(A, obs[["SMOKER"]], obs[["BMIBL"]],
                                        t1_coefs, output = "F_inv.u", u = runif(n, 0, 1))$F_inv.u,
                           "T2" = f2_fn(A, obs[["STROKSFL"]], obs[["MIFL"]],
                                        t2_coefs, output = "F_inv.u", u = runif(n, 0, 1))$F_inv.u,
                           "T3" = f3_fn(t3_coefs, output = "F_inv.u", u = runif(n, 0, 1))$F_inv.u)
    outcomes <- as.data.table(t(apply(outcomes, 1, function(r) {
        J = which.min(r)
        c("EVENT" = r[J] <= tau, "J" = J)
    })))
    Psi0[a+1, ] <- (count(outcomes[outcomes[["V1"]] == 1, ], V2) / n)[["n"]]
}
rm(outcomes); rm(obs); gc()
Psi0[3, ] <- Psi0[2, ] - Psi0[1, ]
Psi0 <- cbind("A" = c("1", "0", "ATE"), Psi0[c(2, 1, 3), ])



# 3.2 contmle Estimation ----------------------------------------------------------------------

source("../continuousTMLE/R/contmle.R")
B <- 400
registerDoParallel(cores = 8)
est <- foreach(i=1:B,
               .combine = rbind,
               .packages = c("data.table", "tidyverse", "survival", "zoo")) %dopar% {

                   obs <- simulate_data(base_data = base_data)
                   est <- contmle(obs, #-- dataset
                                  target=1:3, #-- go after cause 1 and cause 2 specific risks
                                  iterative=FALSE, #-- use one-step tmle to target F1 and F2 simultaneously
                                  treat.effect="ate", #-- target the ate directly
                                  tau=720, #-- time-point of interest
                                  estimation=list("cause1"=list(fit="cox",
                                                                model=Surv(TIME, EVENT==1)~ARM+SEX+AGE+BMIBL+SMOKER+STROKSFL+MIFL),
                                                  "cens"=list(fit="cox",
                                                              model=Surv(TIME, EVENT==0)~ARM+SEX+AGE+BMIBL+SMOKER+STROKSFL+MIFL),
                                                  "cause2"=list(fit="cox",
                                                                model=Surv(TIME, EVENT==2)~ARM+SEX+AGE+BMIBL+SMOKER+STROKSFL+MIFL),
                                                  "cause3"=list(fit="cox",
                                                                model=Surv(TIME, EVENT==3)~ARM+SEX+AGE+BMIBL+SMOKER+STROKSFL+MIFL)
                                  ),
                                  treat.model = ARM ~ SEX+AGE+BMIBL+SMOKER+STROKSFL+MIFL,
                                  sl.models=list(mod1=list(Surv(TIME, EVENT==1)~ARM),
                                                 mod2=list(Surv(TIME, EVENT==1)~ARM+SEX+AGE+BMIBL),
                                                 mod3=list(Surv(TIME, EVENT==1)~ARM+SEX+AGE+BMIBL+SMOKER+STROKSFL+MIFL))
                   )
                   estimates <- bind_rows(unlist(bind_cols(est$init)[1, ]),
                                          unlist(bind_cols(est$tmle)[1, ]))
                   estimates <- as.data.table(cbind("run" = i, "Estimator" = c("init", "TMLE"), estimates))
                   estimates
               }
