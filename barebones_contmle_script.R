
# setup -------------------------------------------------------------------

library(MOSS); library(survtmle); library(survival); library(zoo)
library(tidyverse); library(foreach); library(doRNG); library(doParallel)
library(data.table); library(tlverse); library(sl3); library(origami)
setwd(dir = "/Shared/Projects/ConCR-TMLE/")
lapply(paste0("R/functions/", list.files("R/functions/")), source)
source("./R/barebones_contmle.R")

# data cleaning -------------------------------------------------------------------------------
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

data <- simulate_data(n = 1e3, base_data = base_data)



# sl ----------------------------------------------------------------------





target_events <- 1:3
target_times <- seq(from = 100, to = 1600, by = 100)

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

psi0 <- do.call(cbind, readRDS("./data/true_risks.RDS"))[target_times, ]


# generate data -------------------------------------------------------------------------------

cl8 <- makeForkCluster(8)
registerDoParallel(cl8)
registerDoRNG(123456789)
B <- 50
sim_data <- foreach(i = 1:B) %dopar% {
  return(simulate_data(n = 1e3, base_data = base_data) %>%
           mutate(ARM = as.numeric(ARM)))
}

target_times_cont <- 1:4 * 360
target_events <- 1:3

cont_fit <- foreach(data = sim_data[1:8]) %dopar% {
  out <- concr_tmle(data, target_times_cont, target_events, models)
  tmle_out <- out$estimates$tmle[, -"S"] %>% melt(., id.vars = c("A", "time")) %>%
    dcast(., time + variable ~ A, value.var = "value") %>%
    .[, RR := `1` / `0`] %>%
    .[, RD := `1` - `0`] %>%
    .[, SR := (1 - `1`) / (1 - `0`)]
  setnames(tmle_out, c("variable", "0", "1"), c("J", "F.a0", "F.a1"))
  setcolorder(tmle_out, c("J", "time", "F.a1", "F.a0"))
  return(tmle_out)
}

fit_1mo <- c(fit_1mo, lapply(sim_data[1:10], function(data) {
  dat <- copy(data)[, TIME := ceiling(TIME / 30)]
  out <- concr_tmle(dat, target_times_cont / 30, target_events, models)
  tmle_out <- out$estimates$tmle[, -"S"] %>% melt(., id.vars = c("A", "time")) %>%
    dcast(., time + variable ~ A, value.var = "value") %>%
    .[, RR := `1` / `0`] %>%
    .[, RD := `1` - `0`] %>%
    .[, SR := (1 - `1`) / (1 - `0`)]
  setnames(tmle_out, c("variable", "0", "1"), c("J", "F.a0", "F.a1"))
  setcolorder(tmle_out, c("J", "time", "F.a1", "F.a0"))
  return(tmle_out)
}))

fit_3mo <- c(fit_3mo, lapply(sim_data[1:30], function(data) {
  dat <- copy(data)[, TIME := ceiling(TIME / 30 / 3)]
  out <- concr_tmle(dat, target_times_cont / 30 / 3, target_events, models)
  tmle_out <- out$estimates$tmle[, -"S"] %>% melt(., id.vars = c("A", "time")) %>%
    dcast(., time + variable ~ A, value.var = "value") %>%
    .[, RR := `1` / `0`] %>%
    .[, RD := `1` - `0`] %>%
    .[, SR := (1 - `1`) / (1 - `0`)]
  setnames(tmle_out, c("variable", "0", "1"), c("J", "F.a0", "F.a1"))
  setcolorder(tmle_out, c("J", "time", "F.a1", "F.a0"))
  return(tmle_out)
}))

fit_6mo <- c(fit_6mo, lapply(sim_data[1:30], function(data) {
  dat <- copy(data)[, TIME := ceiling(TIME / 30 / 6)]
  out <- concr_tmle(dat, target_times_cont / 30 / 6, target_events, models)
  tmle_out <- out$estimates$tmle[, -"S"] %>% melt(., id.vars = c("A", "time")) %>%
    dcast(., time + variable ~ A, value.var = "value") %>%
    .[, RR := `1` / `0`] %>%
    .[, RD := `1` - `0`] %>%
    .[, SR := (1 - `1`) / (1 - `0`)]
  setnames(tmle_out, c("variable", "0", "1"), c("J", "F.a0", "F.a1"))
  setcolorder(tmle_out, c("J", "time", "F.a1", "F.a0"))
  return(tmle_out)
}))

