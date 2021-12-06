
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

n <- 1e3

# 0.3 simulation functions ------------------------------------------------

source("R/functions/sim_functions.R")


# 3. Estimation -------------------------------------------------------------------------------

# Estimation targets
tau <- 720
target <- 1:3
B <- 1000
n_cores <- 10
registerDoParallel(n_cores)

# 3.1 True Psi --------------------------------------------------------------------------------
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
    rm(obs); rm(A); gc()
    
    true_risks[[paste0("A=", a)]] <- foreach(t = interval,
                                             .combine = rbind,
                                             .inorder = T) %dopar% {
                                               tabulate(outcomes[["J"]][outcomes[["T"]] <= t])
                                             }
    true_risks[[paste0("A=", a)]] <- as.data.table(true_risks[[paste0("A=", a)]] / nrow(obs)) %>%
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

# plot survival curves
# lapply(true_risks, function(r) rename_all(r, ~c("J=1", "J=2", "J=3"))) %>%
#   bind_rows() %>% mutate(`A` = rep(1:0, each = nrow(true_risks[[1]])),
#                          `t` = rep(1:nrow(true_risks[[1]]), times = 2)) %>%
#   pivot_longer(cols = c(`J=1`, `J=2`, `J=3`), names_to = "event",
#                values_to = "risk") %>%
#   ggplot(aes(x = `t`, y = risk,
#              linetype = as.character(A), colour = event)) +
#   geom_line()+ theme_minimal()

melt(true_risks, id.vars = c("A", "time"), variable.name = "Parameter") %>% 
  mutate(A = as.character(A)) %>% ggplot() + 
  geom_line(aes(x = time, y = value, colour = Parameter, linetype = A)) + 
  theme_minimal()


# 3.2 contmle Estimation ----------------------------------------------------------------------

if (file.exists("./output/contmle_estimates.RDS")) {
  estimates <- readRDS("./output/contmle_estimates.RDS")
} else { 
  estimates <- foreach(i=1:B,
                       .combine = rbind) %dopar% 
    {
      obs <- simulate_data(n = n, base_data = base_data)
      est <- lapply(1:0, function(a) {
        contmle(obs, #-- dataset
                target = target, #-- target competing events
                iterative = FALSE, #-- one-step tmle to target \simultaneously
                treat.effect = a, #-- target the ate directly
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
      })
      
      estimates <- lapply(c("tmle", "init"), function(estimator) {
        lapply(1:0, function(a) {
          est[[2-a]][[estimator]] %>% do.call(cbind, .) %>% t() %>% as.data.frame() %>% 
            rename_all(~paste0(c("Estimate.", "se."), a)) 
        }) %>% bind_cols() %>% mutate("Estimand" = c(paste0("F", 1:3), "S"), 
                                      Estimator = estimator) %>% 
          dplyr::select(Estimator, Estimand, everything())
      }) %>% bind_rows() %>% cbind("run" = i, .)
      return(estimates)
    }
  saveRDS(estimates, "./output/contmle_estimates.RDS")
}


# 4. Evaluation / Visualization ---------------------------------------------------------------


Psi0 <- cbind("time" = tau, 
              rbind(cbind("A" = 1, unname(true_risks[["A=1"]][tau, ])), 
                    cbind("A" = 0, unname(true_risks[["A=0"]][tau, ])))) %>% 
  rename_all(~c("time", "A", "J1", "J2", "J3")) %>% 
  mutate(S = 1 - J1 - J2 - J3) %>% 
  pivot_longer(c("J1", "J2", "J3", "S"), names_to = "Event", values_to = "Truth") %>% 
  pivot_wider(names_from = A, values_from = Truth, names_prefix = "F") %>% 
  mutate(RD = F1 - F0, 
         RR = F1 / F0, 
         SR = (1 - F1) / (1 - F0)) %>% 
  pivot_longer(cols = c(F1, F0, RD, RR, SR), names_to = "Estimand", values_to = "Truth")

result_tbl <- estimates %>% dplyr::select(-run) %>% as_tibble() %>% 
  rename(F1 = Estimate.1, F0 = Estimate.0, se.F1 = se.1, se.F0 = se.0) %>% 
  mutate(RD = F1 - F0, 
         RR = F1 / F0, 
         SR = (1 - F1) / (1 - F0), 
         se.RD = sqrt(se.F1^2 + se.F0^2), 
         se.RR = sqrt(se.F1^2 / (F1)^2 + se.F0^2 * ((F1) / (F0)^2)^2), 
         se.SR = sqrt(se.F1^2 / (1 - F1)^2 + se.F0^2 * (1 - F1)^2 / (1 - F0)^4)) %>% 
  mutate(Estimand = case_when(Estimand == "F1" ~ "J1", 
                              Estimand == "F2" ~ "J2", 
                              Estimand == "F3" ~ "J3", 
                              Estimand == "S" ~ "S")) %>% rename(Event = Estimand) %>% 
  pivot_longer(cols = c(`F0`, `F1`, `RD`, `RR`, `SR`), names_to = "Estimand",
               values_to = "Estimate") %>%
  pivot_longer(cols = contains("se"), names_to = c("tmp", "se_est"),
               names_sep = "\\.", values_to = "se") %>% dplyr::select(-`tmp`) %>%
  filter(se_est == Estimand) %>% dplyr::select(-se_est) %>% left_join(., Psi0) %>% 
  mutate(Bias = Truth - Estimate, MSE = (Truth - Estimate)^2) %>%
  group_by(Estimator, Estimand, Event) %>%
  mutate(upper95 = quantile(Estimate, 0.975),
         lower95 = quantile(Estimate, 0.025), 
         cover = (Estimate + 1.96 * se > Truth) & (Estimate - 1.96 * se < Truth)) %>% 
  group_by(Estimator, Event, Estimand, time) %>% summarise_all(mean)

result_plot <- result_tbl %>%
  ggplot(aes(x = factor(`time`), y = Estimate, colour = Estimator)) +
  facet_wrap(Event~Estimand, scales = "free", ncol = 5) +
  geom_errorbar(aes(ymin = lower95, ymax = upper95), width = .5,
                position = position_dodge(.5)) +
  geom_point(position = position_dodge(.5)) + theme_bw() +
  labs(x = "Days", y = "Estimate", title = "contmle Simulation") +
  geom_label(aes(y = upper95, label = paste0(round(cover, 2)*100, "%")),
             colour = 'black', vjust = 1.2, position = position_dodge2(.5))
