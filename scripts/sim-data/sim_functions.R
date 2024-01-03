# all models are simulated from underlying Weibull distribution

simConCR <- function(interval = 1:2e3,
                     ltfu_coefs = c(7.5e-5, 1, 4, 4),
                     eos_coefs = c("eos_start_time" = 1460, "eos_end_time" = 2000),
                     t1_coefs = c(7.5e-5, 1, 2, 1.2, 2, 1.2),
                     t2_coefs = c(1.5e-5, 1.3, 2, 2, 1.4),
                     t3_coefs = c(0, 1),
                     n = 1e3,
                     assign_A = function(W, n) rbinom(n, 1, 0.5),
                     test_leader.xlsx_path = "/Shared/Projects/concrete/scripts/sim-data/test_leader.xlsx",
                     random_seed = 12345678) {
    library(data.table)
    library(tidyverse)
    base_data <- try(readxl::read_excel(test_leader.xlsx_path))
    
    if (inherits(base_data, "try-error"))
        stop("could not find the 'test_leader.xlsx' document at the specified filepath")
    
    base_data <- base_data %>%
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
    
    
    # 0. parameters -----------------------------------------------------------
    # first two terms in coefs are B and k respectively for weibull hazard = B * k * x^(k - 1)
    
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
    ltfu_fn <- function(ARM, BMIBL, MIFL, params, output = c("h.t", "S.t", "F_inv.u"), t = NULL, u = NULL) {
        B <- params[1]
        k <- params[2]
        b_1 <- log(params[3]) * as.numeric(BMIBL > 30) # * as.numeric(ARM == 0) 
        b_2 <- log(params[4]) * as.numeric(MIFL) # * as.numeric(ARM == 1) 
        phi <- exp(b_1 + b_2)
        
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
        b2 <- log(params[4]) * as.numeric(ARM == 0) # * as.numeric(SMOKER)
        b3 <- log(params[5]) * as.numeric(BMIBL > 30)
        b4 <- log(params[6]) * as.numeric(ARM == 0) # * as.numeric(BMIBL > 30)
        phi <- exp(b1 + b2 + b3 + b4)
        
        return(return_weibull_outputs(phi, B, k, output, t, u))
    }
    
    ### 1.2.2 Event 2 -----------------------------------------------------------------------------
    T2_fn <- function(ARM, STROKSFL, MIFL, params,
                      output = c("h.t", "S.t", "F_inv.u"), t = NULL, u = NULL) {
        B <- params[1]
        k <- params[2]
        b1 <- log(params[3]) * as.numeric(STROKSFL) # * as.numeric(ARM == 0)
        b2 <- log(params[4]) * as.numeric(MIFL) # * as.numeric(ARM == 1)
        b3 <- log(params[5]) * as.numeric(ARM == 1) # as.numeric(STROKSFL) * as.numeric(MIFL)
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
    
    simulate_data <- function(n = 1e3, assign_A = function(W, n) rbinom(n, 1, 0.5), base_data = NULL,
                              random_seed = random_seed) {
        set.seed(random_seed)
        obs <- dplyr::select(sample_n(base_data, size = n, replace = T), -ARM, -TIME, -EVENT)
        A <- assign_A(obs, n)
        outcomes <- data.table("C_ltfu" = ltfu_fn(A, obs[["BMIBL"]], obs[["MIFL"]], ltfu_coefs,
                                                  output = "F_inv.u", u = runif(n, 0, 1))$F_inv.u,
                               "C_eos" = eos_fn(eos_coefs, "F_inv.u", u = runif(n, 0, 1))$F_inv.u,
                               "T1" = T1_fn(A, obs[["SMOKER"]], obs[["BMIBL"]],
                                            t1_coefs, output = "F_inv.u", u = runif(n, 0, 1))$F_inv.u,
                               "T2" = T2_fn(A, obs[["STROKSFL"]], obs[["MIFL"]],
                                            t2_coefs, output = "F_inv.u", u = runif(n, 0, 1))$F_inv.u,
                               "T3" = T3_fn(t3_coefs, output = "F_inv.u", u = runif(n, 0, 1))$F_inv.u)
        
        obs <- 
            cbind(
                cbind(
                    t(apply(outcomes, 1, 
                            function(r) c("TIME" = min(r), 
                                          "EVENT" = max(0, which(r == min(r, na.rm = T)) - 2)))),
                    "ARM" = A), 
                obs)
        obs <- as.data.table(cbind("id" = 1:nrow(obs), obs))
        return(obs)
    }
    
    return(simulate_data(n = n, assign_A = assign_A, base_data = base_data, random_seed = random_seed))
}

getTrueRisks <- function(target_time = 1:2000,
                         n = 1e5,
                         assign_A = function(W, n) return(rep_len(1, n)),
                         ltfu_coefs = c(0, 1, 1, 1),
                         eos_coefs = c(max(target_time) + 1, max(target_time) + 1), 
                         rep = 5, 
                         parallel = FALSE) {
    seeds <- sample(0:999999999, size = rep)
    if (all(parallel, require(doParallel), require(foreach))) {
        n_cores <- min(parallel::detectCores() - 1, 10)
        cl <- makeCluster(n_cores, type = "FORK")
        registerDoParallel(cl)
        risks <- data.table(Reduce("+", foreach(i = seq_along(seeds)) %dopar% {
            outcomes <- simConCR(assign_A = assign_A,
                                 ltfu_coefs = ltfu_coefs,
                                 eos_coefs = eos_coefs,
                                 n = n, 
                                 random_seed = seeds[i])
            Js <- sort(unique(outcomes$EVENT))
            risks <- data.table()
            risks[, (as.character(Js)) := lapply(Js, function(j) {
                    unlist(colSums(outer(
                        outcomes[EVENT == j, TIME], target_time, `<=`
                    )))   
            })]
            return(as.matrix(risks))
        }))
        registerDoSEQ()
        stopCluster(cl)
        rm(cl); gc()
    } else {
        outcomes <- simConCR(assign_A = assign_A,
                             ltfu_coefs = ltfu_coefs,
                             eos_coefs = eos_coefs,
                             n = n, 
                             random_seed = seeds[1])
        Js <- sort(unique(outcomes$EVENT))
        risks <- data.table(matrix(0, nrow = length(target_time), ncol = length(Js)))
        setnames(risks, as.character(Js))
        i <- 1
        while (i <= rep) {
            risks[, (as.character(Js)) := lapply(Js, function(j) { 
                risks[[as.character(j)]] +
                    unlist(colSums(outer(
                        outcomes[EVENT == j, TIME], target_time, `<=`)))   
            })]
            i <- i + 1
            if (i < rep) {
                outcomes <- simConCR(assign_A = assign_A,
                                     ltfu_coefs = ltfu_coefs,
                                     eos_coefs = eos_coefs,
                                     n = n, 
                                     random_seed = seeds[i])
            }
        }
    }
    risks <- subset(risks, select = apply(risks, 2, function(j) max(j) > 0))
    return(cbind("Time" = target_time, risks / (n * rep)))
}
