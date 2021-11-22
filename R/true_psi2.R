library(tidyverse)
library(readxl)
library(skimr)
library(data.table)
library(here)
library(doParallel)
library(foreach)
library(survival)
library(zoo)
setwd("../ConCR-TMLE/")
i_am("./R/contmle-competing-risks-simulation.R")
source("./R/contmle.R")

base_data <- read_excel("./data/test_leader.xlsx") %>%
  mutate_if(is_character, as_factor) %>%
  mutate_if( ~ length(levels(.)) == 2, ~ as.logical(as.numeric(.) - 1)) %>%
  mutate(
    ARM = as.numeric(ARM),
    TIME = time_days,
    EVENT = event,
    SMOKER = case_when(SMOKER == "NEVER SMOKED" ~ 0,
                       SMOKER == "PREVIOUS SMOKER" ~ 1,
                       T ~ 2),
    BMIBL = case_when(is.na(BMIBL) ~ mean(BMIBL, na.rm = T),
                      T ~ BMIBL)
  ) %>%
  dplyr::select(ARM, TIME, everything(),-subjid,-time_days,-event) %>%
  as.data.table()

source("R/functions/sim_functions.R")

set.seed(1245837)
dat =  simulate_data(n = 5e3, base_data = base_data)

colnames(base_data)

ARM <- base_data$ARM
SMOKER <- base_data$SMOKER
BMIBL <- base_data$BMIBL
STROKSFL <- base_data$STROKSFL
MIFL <- base_data$MIFL

t1_coefs <- c(1e-4, 1, 1.2, 1.3, 1.2, 1.2)
t2_coefs <- c(9e-6, 1.3, 1.5, 1.5, 1.2)
t3_coefs <- c(1.8e-5, 1.2)

T1 = T1_fn(
  ARM,
  SMOKER,
  BMIBL,
  t1_coefs,
  output = c("h.t", "S.t", "F_inv.u"),
  t = NULL,
  u = NULL
)

T2 = T2_fn(
  ARM,
  STROKSFL,
  MIFL,
  t2_coefs,
  output = c("h.t", "S.t", "F_inv.u"),
  t = NULL,
  u = NULL
)

T3 = T3_fn(
  t3_coefs,
  output = c("h.t", "S.t", "F_inv.u"),
  t = NULL,
  u = NULL
)


# g = lambda_1( t, a, w) + lambda(j, t, a, w) + lambda(j, t, a, w)

# lambda_1(s, a, w) = T1_fn(ARM = a, SMOKER = w[1], BMIBL = w[2],
#                           params = t1_coefs, output = "h.t", t = s)

w <-
  dplyr::select(base_data, ARM, SMOKER, BMIBL, STROKSFL, MIFL) %>%
  distinct()


# 3d array
result <- data.frame(matrix(0, nrow = nrow(w), ncol = 3))
names(result) =  c("F1", "F2", "F3")

w_all <-
  dplyr::select(base_data, ARM, SMOKER, BMIBL, STROKSFL, MIFL) %>%
  distinct()

for (i in 1:nrow(w_all)) {
  w = w_all[i]
  time_var = 5
  g <- function(s) {
    T1_fn(
      ARM = w[['ARM']],
      SMOKER = w[['SMOKER']],
      BMIBL = w[['BMIBL']],
      params = t1_coefs,
      output = "h.t",
      t = s
    )$h.t +
      T2_fn(
        ARM = w[['ARM']],
        STROKSFL = w[['STROKSFL']],
        MIFL = w[['MIFL']],
        params = t2_coefs,
        output = "h.t",
        t = s
      )$h.t +
      T3_fn(params = t3_coefs,
            output = "h.t",
            t = s)$h.t
  }
  
  InnerIntegral = Vectorize(function(t) {
    integrate(g, lower = 0, upper = t)$value
  })
  # integrate(InnerIntegral, 0, 5)
  
  
  f_1 = Vectorize(function(t) {
    T1_fn(
      ARM = w[['ARM']],
      SMOKER = w[['SMOKER']],
      BMIBL = w[['BMIBL']],
      params = t1_coefs,
      output = "h.t",
      t = t
    )$h.t *
      exp(-integrate(InnerIntegral, 0, t)$value)
  })
  
  f_2 = Vectorize(function(t) {
    T2_fn(
      ARM = w[['ARM']],
      STROKSFL = w[['STROKSFL']],
      MIFL = w[['MIFL']],
      params = t2_coefs,
      output = "h.t",
      t = t
    )$h.t *
      exp(-integrate(InnerIntegral, 0, t)$value)
  })
  
  f_3 = Vectorize(function(t) {
    T3_fn(params = t3_coefs,
          output = "h.t",
          t = t)$h.t *
      exp(-integrate(InnerIntegral, 0, t)$value)
  })
  
  F_1 = integrate(f = f_1, lower = 0, upper = 1)$value
  
  F_2 = integrate(f = f_2, lower = 0, upper = time_var)$value
  
  F_3 = integrate(f = f_3, lower = 0, upper = time_var)$value
  
  result[i, ] = c(F_1, F_2, F_3)
}

write.csv(result,'data/true_psi_calc.csv')
