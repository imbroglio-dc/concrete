
# setup -------------------------------------------------------------------

library(MOSS); library(survtmle); library(survival); library(zoo); library(tlverse)
library(tidyverse); library(foreach); library(doRNG); library(doParallel)
library(data.table)
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


models <- list("0"  = list(mod1 = Surv(TIME, EVENT == 0) ~ ARM,
                           mod2 = Surv(TIME, EVENT == 0) ~ 
                                         ARM + SMOKER + HBA1CBL + BMIBL + AGE,
                           mod3 = Surv(TIME, EVENT == 0) ~ .), 
               "1"  = list(mod1 = Surv(TIME, EVENT == 1) ~ ARM,
                           mod2 = Surv(TIME, EVENT == 1) ~ 
                                         ARM + SMOKER + HBA1CBL + BMIBL + AGE,
                           mod3 = Surv(TIME, EVENT == 1) ~ .), 
               "2"  = list(mod1 = Surv(TIME, EVENT == 2) ~ ARM,
                           mod2 = Surv(TIME, EVENT == 2) ~ 
                                         ARM + SMOKER + HBA1CBL + BMIBL + AGE,
                           mod3 = Surv(TIME, EVENT == 2) ~ .), 
               "3"  = list(mod1 = Surv(TIME, EVENT == 3) ~ ARM,
                           mod2 = Surv(TIME, EVENT == 3) ~ 
                                         ARM + SMOKER + HBA1CBL + BMIBL + AGE,
                           mod3 = Surv(TIME, EVENT == 3) ~ .))


