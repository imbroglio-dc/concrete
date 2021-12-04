
# setup -------------------------------------------------------------------

library(MOSS); library(survtmle); library(survival); library(zoo)
library(tidyverse); library(foreach); library(doRNG); library(doParallel)
library(data.table); library(tlverse); library(sl3); library(origami)
setwd(dir = "/Shared/Projects/ConCR-TMLE/")
lapply(paste0("R/functions/", list.files("R/functions/")), source)


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

data <- simulate_data(base_data = base_data)



# sl ----------------------------------------------------------------------



logreg <- make_learner(Lrnr_glm)
lasso <- make_learner(Lrnr_glmnet) # alpha default is 1
ridge <- Lrnr_glmnet$new(alpha = 0)
e_net <- make_learner(Lrnr_glmnet, alpha = 0.5)
a_lrnrs <- make_learner(Stack, logreg, lasso, ridge, e_net)

target_events <- 1:3
target_times <- seq(from = 100, to = 1600, by = 100)

models <- list("A" = a_lrnrs, 
               "0" = list(mod1 = Surv(TIME, EVENT == 0) ~ ARM,
                           mod2 = Surv(TIME, EVENT == 0) ~ ARM + AGE,
                           mod3 = Surv(TIME, EVENT == 0) ~ ARM + AGE + SMOKER + STROKSFL, 
                           mod4 = Surv(TIME, EVENT == 0) ~ .), 
               "1" = list(mod1 = Surv(TIME, EVENT == 1) ~ ARM,
                           mod2 = Surv(TIME, EVENT == 1) ~ ARM + SMOKER + BMIBL,
                           mod3 = Surv(TIME, EVENT == 1) ~ ARM*SMOKER + I(BMIBL>30)*ARM, 
                           mod4 = Surv(TIME, EVENT == 1) ~ .), 
               "2" = list(mod1 = Surv(TIME, EVENT == 2) ~ ARM,
                           mod2 = Surv(TIME, EVENT == 2) ~ ARM + STROKSFL + MIFL,
                           mod3 = Surv(TIME, EVENT == 2) ~ ARM*STROKSFL + ARM*MIFL + MIFL:STROKSFL, 
                           mod4 = Surv(TIME, EVENT == 2) ~ .), 
               "3" = list(mod1 = Surv(TIME, EVENT == 3) ~ ARM,
                           mod2 = Surv(TIME, EVENT == 3) ~ ARM + SMOKER + BMIBL,
                           mod3 = Surv(TIME, EVENT == 3) ~ ARM+ SMOKER + BMIBL + STROKSFL + MIFL, 
                           mod4 = Surv(TIME, EVENT == 3) ~ .))

true_risks <- readRDS("./data/true_risks.RDS")
psi0 <- do.call(cbind, true_risks)[target_times, ]
