# one step contmle for ATE for competing events

# read_csv("data/true_psi_calc.csv")

# 0.0.1 Packages and Base Data ------------------------------------------------------------------

library(tidyverse)
library(readxl)
library(skimr)
library(data.table)
library(here)
library(doParallel)
library(foreach)
library(survival)
library(zoo)
setwd("Research/Projects/ConCR-TMLE/")
i_am("./R/contmle-competing-risks-simulation.R")
source("./R/contmle.R")

set.seed(0)
base_data <- read_excel("./data/test_leader.xlsx") %>%
  mutate_if(is_character, as_factor) %>%
  mutate_if(~ length(levels(.)) == 2, ~ as.logical(as.numeric(.) - 1)) %>%
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
  dplyr::select(ARM, TIME, everything(), -subjid, -time_days, -event) %>%
  as.data.table()

# 0. set up ---------------------------------------------------------------

## specifications

# for data
n <- 1e3
# run sim_functions.R
obs <- simulate_data(n = n, base_data = base_data)
# Estimation targets
tau <- 720 # time point of interest
target <- 1:3 # for each target
B <- 800 # repeat estimations

# let us just use target = 1 for starters
# target = 1
dt = obs = simulate_data(n = n, base_data = base_data)

estimation = list(
  "cens" = list(
    fit = "cox",
    model = Surv(TIME, EVENT == 0) ~
      ARM + SEX + AGE + BMIBL +
      SMOKER + STROKSFL + MIFL
  ),
  "cause1" = list(
    fit = "cox",
    model = Surv(TIME, EVENT == 1) ~
      ARM + SEX + AGE + BMIBL +
      SMOKER + STROKSFL + MIFL
  ),
  "cause2" = list(
    fit = "cox",
    model = Surv(TIME, EVENT == 2) ~
      ARM + SEX + AGE + BMIBL +
      SMOKER + STROKSFL + MIFL
  ),
  "cause3" = list(
    fit = "cox",
    model = Surv(TIME, EVENT == 3) ~
      ARM + SEX + AGE + BMIBL +
      SMOKER + STROKSFL + MIFL
  )
)

iterative = FALSE #-- one-step tmle to target \simultaneously
treat.model = ARM ~ SEX + AGE + BMIBL + SMOKER + STROKSFL + MIFL
sl.models = list(
  mod1 = list(Surv(TIME, EVENT == 1) ~ ARM),
  mod2 = list(Surv(TIME, EVENT == 1) ~
                ARM + SEX + AGE + BMIBL),
  mod3 = list(Surv(TIME, EVENT == 1) ~
                ARM + SEX + AGE + BMIBL +
                SMOKER + STROKSFL + MIFL)
)
treat.effect = "ate"


# 0.1 other specifications ------------------------------------------------

hal.screening = FALSE
one.step = FALSE
deps.size = 0.1
no.small.steps = 200
push.criterion = FALSE
no.update.if.solved = FALSE
separate.cr = FALSE
simultaneous.ci = FALSE
cut.off.cens = 1e-3
weighted.norm = c(FALSE, "sigma", "Sigma")
pi.star.fun = function(L)
  L$L1 * 0.2
use.observed.times = FALSE
length.times = length(tau)
check.times.size = NULL # 100,
delta.min = 0.01
check.min = min(tau) + 0.1
check.max = max(tau) - 0.1
sl.method = 3
V = 5
lambda.cv = NULL
maxit = 1e3
init.lambda.cvs = c(sapply(1:5, function(jjj)
  (9:1) / (10 ^ jjj)))
lambda.cvs = seq(0.0000001, 0.01, length = 50) #seq(0, 0.008, length=51)[-1],
lambda.grid.size = 50
penalize.time = FALSE
cut.covars = 8
cut.time = 10
cut.time.A = 10
cut.L.A = 8
cut.L.interaction = 3
maxIter = 10
verbose = FALSE
verbose.sl = FALSE
check.sup = FALSE
output.km = FALSE
only.km = FALSE
only.cox.sl = FALSE
output.RR = FALSE
sl.change.points = (0:12) / 10


# 0.2. barebones set up ---------------------------------------------------

## 0.1 names of time variable and event (delta) variable ----------------------------------
time.var <- "TIME"
delta.var <- "EVENT"

## 0.2 add event value to list of estimation, and only include those observed -------------
estimation <- lapply(estimation, function(x) {
  event <-
    as.numeric(gsub("\\D", "", unlist(strsplit(
      as.character(x[["model"]])[2], ","
    ))[2]))
  x[[length(x) + 1]] <- event
  names(x)[length(x)] <- "event"
  return(x)
})

# grab out index of which ones are competing events (excluding censoring)
outcome.index <-
  unlist(lapply(1:length(estimation), function(each) {
    event <- estimation[[each]][["event"]]
    each[event > 0]
  }))

target.S <- FALSE

# are there competing risks? assumed true so
cr <- TRUE
# also ignoring check for this
# if (length(dt[get(delta.var) > 0, unique(get(delta.var))]) > 1)
#   cr <- TRUE
# else
#   cr <- FALSE
# if (!cr)
#   target <-
#   1
# else if (length(dt[get(delta.var) > 0, unique(get(delta.var))]) < 3)
#   target <- target[target < 3]

## 0.6 get number of subjects -------------------------------------------------------------
n <- length(dt[, unique(id)])

## 0.7 get treatment colname -------------------------------------------------------------
A.name <- as.character(treat.model)[2]
A.name <- "ARM"

## 0.8 list of covariates ------------------------------------------------------------------

# instead of this we can just like... list the covariates as input
covars <- c("SEX",
            "AGE",
            "BMIBL",
            "SMOKER",
            "STROKSFL",
            "MIFL")

unique.times <- sort(unique(dt[, get(time.var)]))
unique.times2 <- sort(unique(dt[get(delta.var) > 0, get(time.var)]))


a <- c(1, 0)

dt2 <- NULL
bhaz.cox <-
  do.call("rbind", lapply(a, function(aa)
    data.table(
      time = c(0, unique.times), A = aa
    )))
setnames(bhaz.cox, "A", A.name)


# we can probably add to estimate propensity
# 1. Estimate Propensity ------------------------------------------------------------------
prob.A <-
  predict(glm(as.formula(deparse(treat.model)), data = dt), type = "response")

for (each in 1:length(estimation)) {
  fit <- estimation[[each]][["fit"]][1]
  fit.model <- estimation[[each]][["model"]]
  fit.delta <- estimation[[each]][["event"]]
  fit.name <- names(estimation)[each]
  fit.cox <-
    coxph(as.formula(paste0(deparse(fit.model), collapse = "")), data = dt)
  
  estimation[[each]]$fit.cox <- fit.cox
  
  ### Baseline hazard if event is modeled using Cox YES
  bhaz.cox <- merge(bhaz.cox,
                    rbind(data.table(time = 0, hazard = 0),
                          suppressWarnings(setDT(
                            basehaz(# Compute the predicted survival curve for a Cox model.
                              fit.cox, centered = TRUE)
                          ))),
                    by = "time",
                    all.x = TRUE)
  
  
  bhaz.cox[, hazard := na.locf(hazard), by = A.name]
  bhaz.cox[, (paste0("dhaz", fit.delta)) := c(0, diff(hazard)), by = A.name]
  setnames(bhaz.cox, "hazard", paste0("chaz", fit.delta))
}

setnames(bhaz.cox, c("time", A.name), c(time.var, A.name))

dt2 <- copy(dt)
dt2.a <- do.call("rbind", lapply(a, function(aa) {
  dt.tmp <- copy(dt2)
  dt.tmp[, (A.name) := aa]
}))

# TMLE update step

bhaz.cox[, chaz0.1 := c(0, chaz0[-.N])]
mat <- do.call("rbind", lapply(1:n, function(i) {
  tmp <-
    cbind(id = i, bhaz.cox) # attach subject id=i to matrix of hazards
  tmp[, time.obs := dt[id == i, get(time.var)]] # get failure time for subject i
  tmp[, delta.obs := dt[id == i, get(delta.var)]] # get event type for subject i
  tmp[, A.obs := dt[id == i, get(A.name)]] # get treatment assignment for person i
  tmp[, prob.A := prob.A[i]] # get propensity for person i
  
  ## Clever covariate (inverse intervention weights)
  # for ATE the length(a) == 2
  tmp[, Ht := -((get(A.name) == 1) - (get(A.name) == 0)) / # treatment and censoring weights
        ((prob.A ^ get(A.name) * (1 - prob.A) ^ (1 - get(A.name))))]
  
  for (each in 1:length(estimation)) {
    dt2.a[, (paste0("fit.cox", estimation[[each]][["event"]])) :=
            predict(estimation[[each]][["fit.cox"]], newdata = dt2.a, type = "risk")]
    dt2.a[get(paste0("fit.cox", estimation[[each]][["event"]])) > 500,
          (paste0("fit.cox", estimation[[each]][["event"]])) := 500]
    tmp <- merge(tmp, unique(dt2.a[id == i,
                                   c(A.name, paste0("fit.cox", estimation[[each]][["event"]])),
                                   with = FALSE]), by = c(A.name))
  }
  
  # censoring survival denominator
  if (any(unlist(lapply(estimation, function(x)
    x[["event"]])) == 0)) {
    tmp[, surv.C1 := exp(-fit.cox0 * chaz0.1)]
    if (tmp[, min(surv.C1) < cut.off.cens]) {
      tmp[surv.C1 < cut.off.cens, surv.C1 := cut.off.cens]
    }
    tmp[, Ht := Ht / surv.C1]
  }
}))

mat[, surv.t := 1]
for (each in outcome.index) {
  fit.delta <- estimation[[each]][["event"]]
  mat[, surv.t := surv.t * exp(-cumsum(get(paste0("dhaz", fit.delta)) *
                                         get(paste0(
                                           "fit.cox", fit.delta
                                         )))),
      by = c("id", A.name)]
}

# append surv = 0 to the start of every subject's survival prob vector
mat[, surv.t1 := c(0, surv.t[-.N]), by = c("id", A.name)]

for (kk in 1:length(tau)) {
  mat[, (paste0("surv.tau", kk)) :=
        surv.t[get(time.var) == max(get(time.var)[get(time.var) <= tau[kk]])],
      by = c("id", A.name)]
}

# cause-specific, conditional CIFs
if (cr) {
  for (each in outcome.index) {
    # looping over each target event
    fit.delta <- estimation[[each]][["event"]]
    mat[, (paste0("F", fit.delta, ".t")) := cumsum(surv.t * get(paste0("dhaz", fit.delta)) *
                                                     get(paste0("fit.cox", fit.delta))),
        by = c("id", A.name)]
    for (kk in 1:length(tau)) {
      # cause-specific, conditional CIF at target times
      mat[, (paste0("F", fit.delta, ".tau", kk)) :=
            get(paste0("F", fit.delta, ".t"))[get(time.var) ==
                                                max(get(time.var)[get(time.var) <= tau[kk]])],
          by = c("id", A.name)]
      for (each2 in outcome.index) {
        ### 4.3.3 the competing events part of clever covariates ----
        fit.delta2 <- estimation[[each2]][["event"]]
        if (fit.delta == fit.delta2) {
          mat[surv.t > 0, (paste0("Ht", fit.delta, ".lambda", fit.delta2, ".", kk)) :=
                -(1 - (get(
                  paste0("F", fit.delta, ".tau", kk)
                ) -
                  get(paste0(
                    "F", fit.delta, ".t"
                  ))) / surv.t)]
          mat[surv.t == 0, (paste0("Ht", fit.delta, ".lambda", fit.delta2, ".", kk)) := -1]
        } else {
          mat[surv.t > 0, (paste0("Ht", fit.delta, ".lambda", fit.delta2, ".", kk)) :=
                (get(paste0("F", fit.delta, ".tau", kk)) - get(paste0("F", fit.delta, ".t"))) / surv.t]
          mat[surv.t == 0, (paste0("Ht", fit.delta, ".lambda", fit.delta2, ".", kk)) := 0]
        }
      }
    }
  }
}

# Initial estimate of marginal cause-specific CIFs (subdistributions)
if (cr) {
  # Mean cause-specific risk, at target time(s)
  init.fit <- lapply(target, function(each) {
    sapply(1:length(tau), function(kk) {
      mean(rowSums(sapply(a, function(aa)
        (2 * (
          aa == a[1]
        ) - 1) * (mat[get(A.name) == aa, get(paste0("F", estimation[[outcome.index[each]]][["event"]],
                                                    ".tau", kk))[1],
                      by = "id"][, 2][[1]]))))
    })
  })
  names(init.fit) <-
    paste0("F", sapply(outcome.index[target], function(each)
      estimation[[each]][["event"]]))
}

# EIC
# function to evaluate influence function with competing risks
if (cr) {
  eval.ic <-
    function(mat,
             fit,
             target.index = outcome.index,
             Sigma = FALSE,
             tau.values = tau,
             survival = FALSE,
             tau.all = tau) {
      outer <-
        lapply(target.index, function(each) {
          # loop over target times
          fit.delta <- estimation[[each]][["event"]]
          each.index <-
            (1:length(target.index))[target.index == each]
          sapply(1:length(tau.values), function(kk) {
            k2 <-
              (1:length(tau.all))[tau.all == max(tau.all[tau.all <= tau.values[kk]])]
            out <- 0
            for (each2 in outcome.index) {
              # loop over events
              fit.delta2 <- estimation[[each2]][["event"]]
              mat[get(paste0("Ht", fit.delta, ".lambda", fit.delta2, ".", k2)) <= -500,
                  (paste0("Ht", fit.delta, ".lambda", fit.delta2, ".", k2)) := -500]
              out2 <- mat[, sum(
                (get(A.name) == A.obs) * (get(time.var) <= tau.values[kk]) *
                  (get(time.var) <= time.obs) * Ht *
                  get(
                    paste0("Ht", fit.delta, ".lambda", fit.delta2, ".", k2)
                  ) *
                  ((delta.obs == fit.delta2 &
                      get(time.var) == time.obs) -
                     get(paste0("dhaz", fit.delta2)) * get(paste0("fit.cox", fit.delta2))
                  ),
                na.rm = TRUE
              ), by = "id"]
              out <- out + out2[, 2][[1]]
            }
            
            ic.squared <- (out + rowSums(sapply(a, function(aa)
              (
                2 * (aa == a[1]) - 1
              ) * (mat[get(A.name) == aa, get(paste0("F", fit.delta,
                                                     ".tau", k2))[1],
                       by = "id"][, 2][[1]]))) -
                (length(a) == 2) * fit[[paste0("F", fit.delta)]][kk] -
                (length(a) == 1) * fit[[paste0("F", fit.delta)]][kk])
            
            
            return(sqrt(mean(ic.squared ^ 2) / n))
          })
        })
      
      names(outer) <-
        paste0("F", sapply(target.index, function(each)
          estimation[[each]]["event"]))
      return(outer)
      
    }
}

# # compute initial eic
# init.ic <-
#   eval.ic(mat, fit = init.fit, target.index = outcome.index[target])
# 
# init.list <- rbind(init.est = init.fit, init.se = init.ic)
# # colnames(init.list) <- paste0("tau=", tau)
# tmle.list <- list(init = init.list)

# create initial eic for tmle storage
if (cr) {
  init.list <-  lapply(1:length(init.fit), function(each.index) {
    out <- rbind(init.est = init.fit[[each.index]],
                 init.se = init.ic[[each.index]])
    colnames(out) <- paste0("tau=", tau)
    return(out)
  })
  names(init.list) <-
    paste0("F", sapply(outcome.index[target], function(each)
      estimation[[each]]["event"]))
  tmle.list <- list(init = init.list)
}

# evaluate equation
if (cr) {
  eval.equation <- function(mat,
                            eps = 0,
                            target.index = outcome.index,
                            cr.index = outcome.index,
                            tau.values = tau,
                            tau.all = tau) {
    outer <- lapply(target.index, function(each) { # loop over events
      fit.delta <- estimation[[each]][["event"]]
      sapply(1:length(tau.values), function(kk) { # loop over target times
        k2 <- (1:length(tau.all))[tau.all == max(tau.all[tau.all <= tau.values[kk]])]
        # k2 <- which(tau.all  ==  max(tau.all  <=  tau.values[kk]))
        # or just sort the tau's and we don't have to be so complex
        # sum(tau.all <= tau.values[kk])
        
        out <- 0
        for (each2 in cr.index) {
          fit.delta2 <- estimation[[each2]][["event"]]
          out2 <- mat[(get(time.var) <= tau.values[kk]) & (get(time.var) <= time.obs),
                      sum((get(A.name) == A.obs) * (get(time.var) <= time.obs) * Ht *
                            get(paste0("Ht", fit.delta, ".lambda", fit.delta2, ".", k2)) *
                            ((delta.obs == fit.delta2 & get(time.var) == time.obs) -
                               exp(eps * Ht * get(paste0("Ht", fit.delta, 
                                                         ".lambda", fit.delta2, ".", k2))) *
                               get(paste0("dhaz", fit.delta2)) * get(paste0("fit.cox", fit.delta2))
                            )), by = "id"]
          out <- out + out2[, 2][[1]]
        }
        return(mean(out))
      })
    })
    names(outer) <-
      paste0("F", sapply(target.index, function(each)
        estimation[[each]]["event"]))
    return(outer)
  }
}

# one step tmle
second.round <- FALSE
if (cr) {
  Pn.eic.fun <- function(mat) {
    eval <- eval.equation(mat, eps = 0, target.index = outcome.index[target])
    return(eval)
  }
}

Pn.eic <- Pn.eic.fun(mat)
Pn.eic.norm.fun <- function(x2, x) {
  return(sqrt(sum(unlist(x2) * unlist(x))))
}

Pn.eic2.fun <- function(Pn.eic)
  Pn.eic

Pn.eic2 <- Pn.eic2.fun(Pn.eic)

criterion <- 1 / (sqrt(n) * log(n))
Pn.eic.norm.prev <- Pn.eic.norm <- Pn.eic.norm.fun(Pn.eic2, Pn.eic)

for (step in 1:no.small.steps) { 
  if (cr) { ## making duplicate columns for c-s hazards and clever covs
    for (each in outcome.index) { # loop over events
      fit.delta <- estimation[[each]][["event"]]
      mat[, (paste0("fit.cox", fit.delta, ".tmp")) := get(paste0("fit.cox", fit.delta))]
      for (kk in 1:length(tau)) { # loop over target times
        if (target.S)
          mat[, (paste0("Ht", ".lambda", fit.delta, ".", kk, ".tmp")) :=
                get(paste0("Ht", ".lambda", fit.delta, ".", kk))]
        for (each2 in outcome.index[target]) { # nested loop over events
          fit.delta2 <- estimation[[each2]][["event"]]
          mat[, (paste0("Ht", fit.delta2, ".lambda", fit.delta, ".", kk, ".tmp")) :=
                get(paste0("Ht", fit.delta2, ".lambda", fit.delta, ".", kk))]
        }
      }
    }
  } 
  
  if (cr) {
    if (target.S) {
      for (each in outcome.index) {
        fit.delta <- estimation[[each]][["event"]]
        mat[, (paste0("delta", fit.delta, ".dx")) := 0]
        for (kk in 1:length(tau)) {
          mat[, (paste0("delta", fit.delta, ".dx")) :=
                get(paste0("delta", fit.delta, ".dx")) +
                (get(time.var) <= tau[kk]) * Ht * ((get(
                  paste0("Ht", ".lambda", fit.delta, ".", kk)
                )) *
                  Pn.eic2[kk]) / Pn.eic.norm]
        }
      }
    } else { ## 5.1 calculate update step direction? -------------------------------------
      if (!separate.cr) {
        for (each in outcome.index) { # loop over event hazards
          fit.delta <- estimation[[each]][["event"]]
          mat[, (paste0("delta", fit.delta, ".dx")) := 0] # fluctuation of cs-haz
          for (each2 in outcome.index[target]) { # loop over target events
            fit.delta2 <- estimation[[each2]][["event"]]
            for (kk in 1:length(tau)) { # loop over target times
              mat[, (paste0("delta", fit.delta, ".dx")) :=
                    get(paste0("delta", fit.delta, ".dx")) + 
                    (get(time.var) <= tau[kk]) * Ht * 
                    ((get( paste0( "Ht", fit.delta2, ".lambda", fit.delta, ".", kk ))) *
                       Pn.eic2[each2 == outcome.index[target]][[1]][kk]) / Pn.eic.norm]
            }
          }
        }
      } else {
        Pn.eic.separate <- Pn.eic.fun.separate(mat)
        for (each in outcome.index) {
          fit.delta <- estimation[[each]][["event"]]
          mat[, (paste0("delta", fit.delta, ".dx")) := 0]
          for (each2 in outcome.index[target]) {
            fit.delta2 <- estimation[[each2]][["event"]]
            for (kk in 1:length(tau)) {
              mat[, (paste0("delta", fit.delta, ".dx")) :=
                    get(paste0("delta", fit.delta, ".dx")) +
                    (get(time.var) <= tau[kk]) * Ht *
                    ((get(paste0("Ht", fit.delta2, ".lambda", fit.delta, ".", kk))) * 
                       Pn.eic.separate[each2 == outcome.index[target]][[1]][
                         each == outcome.index[target]][[1]][kk]) / Pn.eic.norm]
            }
          }
        }
      }
    }
  }
  
  ## 5.2 set update epsilon down-scaling factor ----------------------------------------
  # why not just use deps.size?
  deps <- deps.size
  
  if (cr) { ## 5.3 update cause specific hazards ----------------------------------------
    mat[, surv.t := 1]
    for (each in outcome.index) { # loop over target events
      fit.delta <- estimation[[each]][["event"]]
      mat[, (paste0("fit.cox", fit.delta)) :=
            get(paste0("fit.cox", fit.delta)) *
            exp(deps * get(paste0("delta", fit.delta, ".dx")))]
      mat[get(paste0("fit.cox", fit.delta)) > 500, (paste0("fit.cox", fit.delta)) := 500]
      # Does changing the order of event updates matter? -------------------------------
      mat[, surv.t := surv.t * exp(-cumsum(get(paste0("dhaz", fit.delta)) *
                                             get(paste0("fit.cox", fit.delta)))),
          by = c("id", A.name)]
      mat[get(paste0("fit.cox", fit.delta)) == Inf, surv.t := 0]
    }
    mat[, surv.t1 := c(0, surv.t[-.N]), by = c("id", A.name)]
    for (kk in 1:length(tau)) {
      mat[, (paste0("surv.tau", kk)) :=
            surv.t[get(time.var) == max(get(time.var)[get(time.var) <= tau[kk]])],
          by = c("id", A.name)]
    }
    for (each in outcome.index[target]) {
      fit.delta <- estimation[[each]][["event"]]
      mat[, (paste0("F", fit.delta, ".t")) := cumsum(surv.t * get(paste0("dhaz", fit.delta)) *
                                                       get(paste0("fit.cox", fit.delta))),
          by = c("id", A.name)]
      for (kk in 1:length(tau)) {
        mat[, (paste0("F", fit.delta, ".tau", kk)) :=
              get(paste0("F", fit.delta, ".t"))[get(time.var) == max(get(time.var)[
                get(time.var) <= tau[kk]])], by = c("id", A.name)]
        for (each2 in outcome.index) {
          fit.delta2 <- estimation[[each2]][["event"]]
          if (fit.delta == fit.delta2) {
            mat[surv.t > 0, (paste0("Ht", fit.delta, ".lambda", fit.delta2, ".", kk)) :=
                  -(1 - (get(paste0("F", fit.delta, ".tau", kk)) - 
                           get(paste0("F", fit.delta, ".t"))) / surv.t)]
            mat[round(surv.t, 8) == 0, 
                (paste0("Ht", fit.delta, ".lambda", fit.delta2, ".", kk)) := -1]
          } else {
            mat[surv.t > 0, (paste0("Ht", fit.delta, ".lambda", fit.delta2, ".", kk)) :=
                  (get(paste0("F", fit.delta, ".tau", kk)) - 
                     get(paste0("F", fit.delta, ".t"))) / surv.t]
            mat[round(surv.t, 8) == 0, 
                (paste0("Ht", fit.delta, ".lambda", fit.delta2, ".", kk)) := 0]
          }
        }
      }
    }
  } 
  
  Pn.eic <- Pn.eic.fun(mat)
  Pn.eic2 <- Pn.eic2.fun(Pn.eic)
  Pn.eic.norm <- Pn.eic.norm.fun(Pn.eic2, Pn.eic)
  
  print(paste0("step = ", step))
  
  if (Pn.eic.norm.prev <= Pn.eic.norm) {
    # reset all columns to '.tmp'
    
    if (verbose) {
      print("----")
      print(step)
      print("----")
    }
    
    if (cr) {
      for (each in outcome.index) {
        fit.delta <- estimation[[each]][["event"]]
        mat[, (paste0("fit.cox", fit.delta)) := get(paste0("fit.cox", fit.delta, ".tmp"))]
        for (kk in 1:length(tau)) {
          if (target.S)
            mat[, (paste0("Ht", ".lambda", fit.delta, ".", kk)) :=
                  get(paste0("Ht", ".lambda", fit.delta, ".", kk, ".tmp"))]
          for (each2 in outcome.index[target]) {
            fit.delta2 <- estimation[[each2]][["event"]]
            mat[, (paste0("Ht", fit.delta2, ".lambda", fit.delta, ".", kk)) :=
                  get(paste0("Ht", fit.delta2, ".lambda", fit.delta, ".", kk, ".tmp"))]
          }
        }
      }
    } else {
      mat[, fit.cox1.tmp := fit.cox1]
      for (kk in 1:length(tau)) {
        mat[, (paste0("Ht.lambda.", kk)) := get((paste0("Ht.lambda.", kk, ".tmp")))]
      }
    }
    
    Pn.eic <- Pn.eic.fun(mat)
    Pn.eic2 <- Pn.eic2.fun(Pn.eic)
    
    Pn.eic.norm <- Pn.eic.norm.fun(Pn.eic2, Pn.eic)
    deps.size <- 0.5 * deps.size # 0.1 * deps.size
    
  } else {
    Pn.eic.norm.prev <- Pn.eic.norm
   
    Pn.eic3 <-
      lapply(1:length(Pn.eic), function(kk)
        Pn.eic[[kk]] / (ifelse(any(unlist(init.ic) == 0), 
                               init.ic[[kk]] + 0.001, 
                               init.ic[[kk]]
        ) * sqrt(n)))
  }
  
  
  Pn.eic3 <-
    lapply(1:length(Pn.eic), function(kk)
      Pn.eic[[kk]] / (ifelse(
        any(unlist(init.ic) == 0), init.ic[[kk]] + 0.001, init.ic[[kk]]
      ) * sqrt(n)))
  
  check.sup.norm <- max(abs(unlist(Pn.eic3))) <= (criterion)
  if (check.sup.norm |
      step == no.small.steps) {
    # (left.criterion <= criterion) {
    if (step == no.small.steps) {
      message("Warning: Algorithm did not converge")
    }
    
    if (verbose)
      print(paste0("converged", " at ", step, "th step"))
    if (verbose)
      print(paste0("eic = ", Pn.eic.fun(mat)))
    
    #-- 12c -- compute sd:
    if (cr) {
      if (treat.effect[1] == "stochastic") {
        final.fit <- lapply(target, function(each) {
          sapply(1:length(tau), function(kk) {
            mean(rowSums(sapply(a, function(aa)
              (
                mat[get(A.name) == aa, pi.star[1] *
                      get(paste0("F", estimation[[outcome.index[each]]][["event"]],
                                 ".tau", kk))[1],
                    by = "id"][, 2][[1]]
              ))))
          })
        })
      } else {
        final.fit <- lapply(outcome.index[target], function(each) {
          sapply(1:length(tau), function(kk) {
            mean(rowSums(sapply(a, function(aa)
              (2 * (aa == a[1]) - 1) * (mat[get(A.name) == aa,
                                            get(paste0("F", estimation[[each]][["event"]],
                                                       ".tau", kk))[1],
                                            by = "id"][, 2][[1]]))))
          })
        })
      }
      names(final.fit) <-
        paste0("F", sapply(outcome.index[target], function(each)
          estimation[[each]][["event"]]))
      final.ic <-
        eval.ic(mat, final.fit, target.index = outcome.index[target])
    } else {
      if (treat.effect[1] == "stochastic") {
        final.fit <- list(sapply(1:length(tau), function(kk) {
          mean(rowSums(sapply(a, function(aa)
            (
              mat[get(A.name) == aa, pi.star[1] * (1 - get(paste0("surv.tau", kk))[1]), 
                  by = "id"][, 2][[1]]
            ))))
        }))
      } else {
        final.fit <- list(sapply(1:length(tau), function(kk) {
          mean(rowSums(sapply(a, function(aa)
            (2 * (aa == a[1]) - 1) * (mat[get(A.name) == aa, 1 - get(paste0("surv.tau", kk))[1], 
                                          by = "id"][, 2][[1]]))))
        }))
      }
      final.ic <- list(eval.ic(mat, final.fit[[1]]))
    }
    
    final.list <-
      lapply(1:length(final.fit), function(each.index) {
        out <- rbind(tmle.est = final.fit[[each.index]],
                     tmle.se = final.ic[[each.index]])
        colnames(out) <- paste0("tau=", tau)
        return(out)
      })
    
    if (cr)
      names(final.list) <-
      paste0("F", sapply(outcome.index[target], function(each)
        estimation[[each]]["event"]))
    else
      final.list <- final.list[[1]]
    
    tmle.list$tmle <- final.list

    if (!second.round) {
      if (check.sup)
        tmle.list$check.sup.norm <- list(
          check.sup.norm = check.sup.norm,
          lhs = max(abs(unlist(Pn.eic3))),
          rhs = criterion / sqrt(length(target) * length(tau))
        )
      tmle.list$convergenced.at.step <- step
      
      if (!push.criterion)
        break
      else
        criterion <- criterion / sqrt(n)
      second.round <- TRUE
    } else {
      if (target.S) {
        tmle.list$tmle.second.round <- final.list[names(final.list) == "S"]
      } else {
        tmle.list$tmle.second.round <- final.list
      }
      if (check.sup)
        tmle.list$check.sup.norm.second.round <- list(
          check.sup.norm =
            max(abs(unlist(Pn.eic3))) <= (criterion),
          lhs = max(abs(unlist(Pn.eic3))),
          rhs = criterion
        )
      tmle.list$convergenced.at.step.second.round <- step
      break
    }
  }
}

return(tmle.list)
