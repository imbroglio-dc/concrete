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
# time.var <- gsub("Surv\\(",
#                  "",
#                  unlist(strsplit(as.character(estimation[[1]][["model"]])[2], ","))[1]
# )
# delta.var <- gsub(" ", "",
#                   gsub(" == ", "",
#                        gsub("[0-9]+\\)", "",
#                             unlist(
#                               strsplit(as.character(estimation[[1]][["model"]])[2], ",")
#                             )[2])
#                   )
# )

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
# include fit i.e. "cox",
# model i.e. Surv(TIME, EVENT == 3) ~
# ARM + SEX + AGE + BMIBL + SMOKER + STROKSFL + MIFL,
# event i.e. 0, 1, 2, 3 where 0 is censoring event all extracted from input

# grab out index of which ones are competing events (excluding censoring)
outcome.index <-
  unlist(lapply(1:length(estimation), function(each) {
    event <- estimation[[each]][["event"]]
    each[event > 0]
  }))

target.S <- FALSE
# # to target the censoring event then target.S is true
# # will be ignored for now
# if (all(target == 0)) {
#   target <- 1:length(outcome.index)
#   target.S <- TRUE
# } else {
#   target.S <- FALSE
# }

# # also don't think this is needed
# # Check if multiple models are for a single event type
# events <- unlist(lapply(estimation, function(each) each[["event"]]))
# if (any(table(events) > 1))
#   stop(paste0("multiple models specified for ", delta.var, "=",
#               paste0(names(table(events))[table(events) > 1], collapse = ",")))
# also ignoring the check for missing event types

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

# covars <- NULL
# for (mod in c(lapply(sl.models, function(x) x[[1]]),
#               unlist(lapply(estimation, function(x) x[["model"]])))) {
#   mod3 <- as.character(mod)[3]
#   covars <- unique(c(covars, unlist(strsplit(gsub("\\+", " ", mod3), " "))))
#   covars <- covars[!covars %in% c(A.name, "", "*")]
#   if (length(grep(".squared", mod3)) > 0) {
#     names.squared <- unique(gsub(".squared", "",
#                                  grep(".squared",
#                                       unlist(strsplit(gsub("\\+", " ", mod3), " ")),
#                                       value = TRUE
#                                  )))
#     for (col in names.squared)
#       dt[, (paste0(col, ".squared")) := get(col)^2]
#   }
#   if (length(grep(".log", mod3)) > 0) {
#     names.log <- unique(gsub(".log", "",
#                              grep(".log", unlist(strsplit(gsub("\\+", " ", mod3), " ")),
#                                   value = TRUE)))
#     for (col in names.log)
#       dt[, (paste0(col, ".log")) := log(get(col))]
#   }
# }

# also ignoring these :)
# if (length(grep(".squared", covars)) > 0)
#   covars <- covars[-grep(".squared", covars)]
# if (length(grep(".log", covars)) > 0)
#   covars <- covars[-grep(".log", covars)]
# if (length(grep("cut.", covars)) > 0)
#   covars <- covars[-grep("cut.", covars)]

# # also assuming the data is "nice"
# # such that every variable has more than one value
# # aka ignoring this
# covars <- covars[covars != "1"]
# for (covar in covars) {
#   if (dt[, length(unique(get(covar)))] == 1)
#     covars <- covars[covars != covar]
# }

## 0.9 get unique times in dataset -----------------------------------------------------
# every time recorded in the dataset
unique.times <- sort(unique(dt[, get(time.var)]))
# every time that doesn't have a censoring event (aka delta.var (EVENT = 0))
unique.times2 <- sort(unique(dt[get(delta.var) > 0, get(time.var)]))

# also assuming this is nice :) i.e. intervals in tau has observations
# ## 0.10 truncate tau if there intervals in tau without obs --------------------------------
# if (use.observed.times) {
#   n <- nrow(dt)
#   unique.times3 <- unique.times[unique.times <= max(tau)]
#   tau <-
#     unique.times3[(1:length(unique.times3)) %% floor(length(unique.times3) / length.times) == 1]
#   #set.seed(1)
#   sort(unique.times3[sample(length(unique.times3), length.times)])
# } else {
#   test.tau <- findInterval(tau, unique.times2)
#   # this gets out the tau index in unique times2
#   if (FALSE & !length(test.tau) == length(unique(test.tau))) {
#     tau <- na.omit(tau[(1:length(tau))[unique(findInterval(unique.times2, tau)) + 1]])
#     warning("no observations between tau as specified, truncating tau")
#   }
# }

# ## 0.11 specify treatment levels of interest for binary treatment ------------------------
# # stochastic treatment check/specified later?
# if (treat.effect[1] == "1")
#   a <- 1
# else if (treat.effect[1] == "0")
#   a <- 0
# else
#   a <- c(1, 0)

# since we have ate
a <- c(1, 0)

## 0.12 initialize duplicate dataset to be used later for? -------------------------------
dt2 <- NULL
bhaz.cox <-
  do.call("rbind", lapply(a, function(aa)
    data.table(
      time = c(0, unique.times), A = aa
    )))
setnames(bhaz.cox, "A", A.name)

# here we have bhaz.cox as time and
# treatment ARM (all set to 1 then all set to 0)
# not sure what to use it for
# this next part is not useful for our code because it returns the exact same

# ## 0.13 prepare variables for if outcome models use coxnet --------------------------------
# sl.models.tmp <- sl.models
# sl.models <- list()
#
# ## 0.14 add separate sl models when specified with, e.g., multiple changepoints -----------
# for (k1 in 1:length(sl.models.tmp)) {
#   if (length(sl.models.tmp[[k1]]) > 1) {
#     for (k2 in 2:length(sl.models.tmp[[k1]])) {
#       sl.models[[length(sl.models) + 1]] <- c(sl.models.tmp[[k1]][1], sl.models.tmp[[k1]][k2])
#       if (length(sl.models.tmp[[k1]]) >= 3) {
#         names(sl.models)[length(sl.models)] <- paste0(names(sl.models.tmp)[k1], k2)
#       } else {
#         names(sl.models)[length(sl.models)] <- names(sl.models.tmp)[k1]
#       }
#     }
#   } else {
#     sl.models[[length(sl.models) + 1]] <- c(sl.models.tmp[[k1]][1])
#     names(sl.models)[length(sl.models)] <- names(sl.models.tmp)[k1]
#   }
# }

# we can probably add to estimate propensity
# 1. Estimate Propensity ------------------------------------------------------------------
prob.A <-
  predict(glm(as.formula(deparse(treat.model)), data = dt), type = "response")

# if (verbose)
#   print(summary(glm(as.formula(deparse(treat.model)), data = dt)))
# if (verbose)
#   print(summary(prob.A))

# 2. Estimate Hazards ---------------------------------------------------------------------
for (each in 1:length(estimation))
{
  fit <- estimation[[each]][["fit"]][1]
  fit.model <- estimation[[each]][["model"]]
  ## ignoring changepoint for now
  # if (any(names(estimation[[each]]) == "changepoint"))
  #   fit.changepoint <- estimation[[each]][["changepoint"]]
  # else
  #   fit.changepoint <- NULL
  fit.delta <- estimation[[each]][["event"]]
  fit.name <- names(estimation)[each]
  
  
  ## 2.1 Cox-SL ---------------------------------------------------------------------------
  if (fit[1] == "sl" | fit[1] == "cox.hal.sl") {
    set.seed(123849348)
    sl.pick <- suppressWarnings(
      cox.sl(
        #???????? where is this from
        loss.fun = cox.loss.fun,
        dt = dt,
        delta.var = delta.var,
        treatment = A.name,
        V = V,
        delta.value = fit.delta,
        cox.models = sl.models,
        change.points = sl.change.points
      )
    )
    cve.sl.pick <- sl.pick$picked.cox.model$cve
    fit.model <- sl.pick$picked.cox.model$form
    
    estimation[[each]]$model <- fit.model
    
  } else {
    ## 2.2 Kaplan-Meier ------------------------------------------------------------
    sl.pick <- ""
    cve.sl.pick <- ""
    
    if (fit[1] %in% c("km")) {
      #-- later for hal or if uses km
      tmp.model <- as.character(fit.model)
      if (fit[1] == "km")
        tmp.model[3] <- "strata(A)"
      else
        tmp.model[3] <- "1"
      estimation[[each]]$model <-
        fit.model <-
        formula(paste0(tmp.model[2], tmp.model[1], tmp.model[3]))
      estimation[[each]]$changepoint <- fit.changepoint <- NULL
    }
  }
  
  estimation[[each]]$sl.pick <- fit.model
  estimation[[each]]$cve.sl.pick <- cve.sl.pick
  
  
  ## 2.3 Cox -------------------------------------------------------------------------------
  # if (length(fit.changepoint) > 0) {
  #   if (length(fit.changepoint) > 0 & length(dt2) == 0) {
  #     dt2 <- rbind(dt, dt)[order(id)]
  #   }
  #   ### 2.3.1 Cox with change-point -----------------------
  #   delta1 <- abs(fit.delta - 1)
  #   dt2[, time.indicator := (get(time.var) <= fit.changepoint)]
  #   dt2[, (paste0("period", fit.delta)) := 1:.N, by = "id"]
  #   dt2[get(paste0("period", fit.delta)) == 1, (paste0("tstart", fit.delta)) := 0]
  #   dt2[get(paste0("period", fit.delta)) == 1, (paste0("tstop", fit.delta)) :=
  #         (get(time.var) <= fit.changepoint) * get(time.var) +
  #         (get(time.var) > fit.changepoint) * fit.changepoint]
  #   dt2[get(paste0("period", fit.delta)) == 1, (paste0("tstart", fit.delta)) := 0]
  #   dt2[get(paste0("period", fit.delta)) == 1, (paste0("tstop", fit.delta)) :=
  #         (get(time.var) <= fit.changepoint) * get(time.var) +
  #         (get(time.var) > fit.changepoint) * fit.changepoint]
  #   dt2[get(paste0("period", fit.delta)) == 2, (paste0("tstart", fit.delta)) := fit.changepoint]
  #   dt2[get(paste0("period", fit.delta)) == 2, (paste0("tstop", fit.delta)) := get(time.var)]
  #   dt2[get(paste0("period", fit.delta)) == 1 &
  #         !time.indicator, (delta.var) := delta1]
  #   mod1 <- as.character(fit.model)
  #   mod2 <-
  #     paste0(gsub(
  #       substr(mod1[2],
  #              which(strsplit(mod1[2], "")[[1]] == "(") + 1,
  #              which(strsplit(mod1[2], "")[[1]] == ",") - 1
  #       ),
  #       paste0("tstart", fit.delta, ", tstop", fit.delta),
  #       mod1[2]),
  #       "~",
  #       gsub(paste0("\\+", A.name, "\\+"), "",
  #            gsub(
  #              paste0("\\+", A.name, " "), "",
  #              gsub(" ", "",
  #                   paste0(
  #                     "I((period", fit.delta, " == 1)&(",
  #                     A.name, " == 1))",
  #                     " + I((period", fit.delta, " == 2)&(",
  #                     A.name, " == 1))",
  #                     " + ", paste0("+", mod1[3])
  #                   )
  #              )
  #            )))
  #   fit.cox <- coxph(formula(mod2), data = dt2[!time.indicator |
  #                                                get(paste0("period", fit.delta)) == 1])
  # } else {
  ### 2.3.2 Cox without change-point ---------------------------------------------
  if (fit[1] == "sl" & length(grep("coxnet", sl.pick)) > 0) {
    X <- model.matrix(as.formula(deparse(fit.model)), data = dt)
    y <- dt[, Surv(get(time.var), get(delta.var) == fit.delta)]
    fit.cox <- glmnet(
      x = X,
      y = y,
      family = "cox",
      maxit = 1000,
      lambda = fit.penalty
    )
  } else {
    fit.cox <-
      coxph(as.formula(paste0(deparse(fit.model), collapse = "")), data = dt)
  }
  # }
  
  estimation[[each]]$fit.cox <- fit.cox
  
  ## 2.4 Baseline Hazards ------------------------------------------------------------------
  
  if (fit[1] == "km") {
    ### 2.4.1 Baseline hazard if event is modeled using Kaplan-Meier -----------------------
    tmp <-
      suppressWarnings(setDT(basehaz(fit.cox, centered = TRUE)))
    setnames(tmp, "strata", "A")
    tmp[, A := as.numeric(gsub("A=", "", A))]
    bhaz.cox <- merge(bhaz.cox,
                      rbind(do.call(
                        "rbind", lapply(a, function(aa)
                          data.table(
                            time = 0,
                            hazard = 0,
                            A = aa
                          ))
                      ),
                      tmp),
                      by = c("time", "A"),
                      all.x = TRUE)
    bhaz.cox[, hazard := na.locff(hazard), by = "A"]
    bhaz.cox[, (paste0("dhaz.", fit.delta)) := c(0, diff(hazard)), by = "A"]
    setnames(bhaz.cox, "hazard", paste0("chaz", fit.delta))
  } else {
    ### 2.4.2 Baseline hazard if event is modeled using SL or coxnet -------------------
    if (fit[1] == "sl" & length(grep("coxnet", sl.pick)) > 0) {
      basehaz <- glmnet_basesurv(dt[, get(time.var)],
                                 dt[, get(delta.var) == fit.delta], X, centered = TRUE)
      bhaz.cox <-
        merge(bhaz.cox,
              rbind(
                data.table(time = 0, hazard = 0),
                data.table(
                  time = basehaz$time,
                  hazard = basehaz$cumulative_base_hazard
                )
              ),
              by = "time",
              all.x = TRUE)
    } else {
      ### 2.4.3 Baseline hazard if event is modeled using Cox YES -----------------------
      bhaz.cox <- merge(
        bhaz.cox,
        rbind(data.table(
          time = 0, hazard = 0
        ),
        suppressWarnings(setDT(
          basehaz(# Compute the predicted survival curve for a Cox model.
            fit.cox, centered = TRUE)
        ))),
        by = "time",
        all.x = TRUE
      )
    }
    bhaz.cox[, hazard := na.locf(hazard), by = A.name]
    bhaz.cox[, (paste0("dhaz", fit.delta)) := c(0, diff(hazard)), by = A.name]
    setnames(bhaz.cox, "hazard", paste0("chaz", fit.delta))
  }
}
## set names of bhaz.cox to match observed data
# just the time variable
setnames(bhaz.cox, c("time", A.name), c(time.var, A.name))
# define chaz0.1 - the censoring survival one time-point back: S^c(t- | A, W)
# this line does NOT work:
# bhaz.cox[, chaz0.1 := c(0, chaz0[-.N])] ??????
# bhaz.cox[, chaz0.1 := c(0, chaz0[-.N])] unsure what this is at all.
# Error in eval(jsub, SDenv, parent.frame()) : object 'chaz0' not found

# ignoring output.km because its False
# ignoring variable prep for HAL
# ignoring ??.?? because dt2 only used with change-point information
dt2 <- copy(dt)
dt2.a <- do.call("rbind", lapply(a, function(aa) {
  dt.tmp <- copy(dt2)
  dt.tmp[, (A.name) := aa]
}))
# 4. TMLE update step --------------------------------------------------------

bhaz.cox[, chaz0.1 := c(0, chaz0[-.N])]
mat <- do.call("rbind", lapply(1:n, function(i) {
  tmp <-
    cbind(id = i, bhaz.cox) # attach subject id=i to matrix of hazards
  tmp[, time.obs := dt[id == i, get(time.var)]] # get failure time for subject i
  tmp[, delta.obs := dt[id == i, get(delta.var)]] # get event type for subject i
  tmp[, A.obs := dt[id == i, get(A.name)]] # get treatment assignment for person i
  tmp[, prob.A := prob.A[i]] # get propensity for person i
  
  ## 4.1 Clever covariate (inverse intervention weights) --------------------------------
  
  # for ATE the length(a) == 2
  if (length(a) == 2) {
    tmp[, Ht := -((get(A.name) == 1) - (get(A.name) == 0)) / # treatment and censoring weights
          ((prob.A ^ get(A.name) * (1 - prob.A) ^ (1 - get(A.name))))]
  } else {
    tmp[, Ht := -((get(A.name) == a)) / # treatment and censoring weights
          ((prob.A ^ get(A.name) * (1 - prob.A) ^ (1 - get(A.name))))]
  }
  
  # for (each in 1:length(estimation)) {
  # again ignoring changepoint stuff
  # if (length(estimation[[each]][["changepoint"]]) > 0) {
  #   tmp[, (paste0("period", estimation[[each]][["event"]])) :=
  #         (get(time.var)  <=  estimation[[each]][["changepoint"]]) *
  #         1 + (get(time.var) > estimation[[each]][["changepoint"]]) * 2]
  #   tmp <- merge(tmp, dt2.a[id == i, c(paste0("period", estimation[[each]][["event"]]),
  #                                      A.name, paste0("fit.cox", estimation[[each]][["event"]])),
  #                           with = FALSE], by = c(paste0("period", estimation[[each]][["event"]]), A.name))
  # } else {
  # tmp <- merge(tmp, unique(dt2.a[id == i,
  # c(A.name, paste0("fit.cox", estimation[[each]][["event"]])),
  # with = FALSE]), by = c(A.name))
  # }
  # }
  
  for (each in 1:length(estimation)) {
    if (estimation[[each]][["fit"]][1] == "sl" &
        length(grep("coxnet", estimation[[each]][["sl.pick"]])) > 0) {
      # if using sl coxnet
      X2.a <-
        model.matrix(as.formula(deparse(estimation[[each]][["model"]])), data = dt2.a)
      dt2.a[, (paste0("fit.cox", estimation[[each]][["event"]])) :=
              exp(predict(estimation[[each]][["fit.cox"]], newx = X2.a, type = "link"))]
    } else {
      dt2.a[, (paste0("fit.cox", estimation[[each]][["event"]])) :=
              predict(estimation[[each]][["fit.cox"]], newdata = dt2.a, type = "risk")]
    }
    dt2.a[get(paste0("fit.cox", estimation[[each]][["event"]])) > 500,
          (paste0("fit.cox", estimation[[each]][["event"]])) := 500]
    tmp <- merge(tmp, unique(dt2.a[id == i,
                                   c(A.name, paste0("fit.cox", estimation[[each]][["event"]])),
                                   with = FALSE]), by = c(A.name))
  }
  # ???????
  # censoring survival denominator
  if (any(unlist(lapply(estimation, function(x)
    x[["event"]])) == 0)) {
    # where is fit.cox0 from? i cannot find.
    # chaz0.1 doesn't work.
    tmp[, surv.C1 := exp(-fit.cox0 * chaz0.1)] # what if censoring is estimated not with cox?
    # if (tmp[, any(is.na(surv.C1))]) browser()
    if (tmp[, min(surv.C1) < cut.off.cens]) {
      tmp[surv.C1 < cut.off.cens, surv.C1 := cut.off.cens]
      # warning(paste0("there seems to be positivity issues; truncated at level ",
      #               cut.off.cens))
    }
    tmp[, Ht := Ht / surv.C1]
  }
}))

# event free survivial at time t
mat[, surv.t := 1]
for (each in outcome.index) {
  fit.delta <- estimation[[each]][["event"]]
  mat[, surv.t := surv.t * exp(-cumsum(get(paste0("dhaz", fit.delta)) *
                                         get(paste0(
                                           "fit.cox", fit.delta
                                         )))),
      by = c("id", A.name)]
}
# ignore poisson HAL initial estimate
## within HAL ignore truncating weights because idk how to get surv.C1
## also ignore 4.3 b/c HAL

##
## 4.4 EIC ---------------------------------------------------------------------------
if (cr) {
  ### 4.4.1 function to evaluate influence function with competing risks ---------
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
              out <- out + out2[, 2][[1]] #EIF 3.1
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
# mat[, surv.t := 1]
# for (each in outcome.index) {
#   fit.delta <- estimation[[each]][["event"]]
#   mat[, surv.t := surv.t * exp(-cumsum(get(paste0("dhaz", fit.delta)) *
#                                          get(paste0(
#                                            "fit.cox", fit.delta
#                                          )))),
#       by = c("id", A.name)]
# }


if (cr) {
  eval.equation <- function(mat,
                            eps = 0,
                            target.index = outcome.index,
                            cr.index = outcome.index,
                            tau.values = tau,
                            tau.all = tau) {
    outer <- lapply(target.index, function(each) {
      # loop over events
      fit.delta <- estimation[[each]][["event"]]
      sapply(1:length(tau.values), function(kk) {
        # loop over target times
        k2 <-
          (1:length(tau.all))[tau.all == max(tau.all[tau.all <= tau.values[kk]])]
        # k2 <- which(tau.all  ==  max(tau.all  <=  tau.values[kk]))
        # or just sort the tau's and we don't have to be so complex
        # sum(tau.all <= tau.values[kk])
        
        out <- 0
        for (each2 in cr.index) {
          fit.delta2 <- estimation[[each2]][["event"]]
          out2 <-
            mat[(get(time.var) <= tau.values[kk]) &
                  (get(time.var) <= time.obs),
                sum((get(A.name) == A.obs) * (get(time.var) <= time.obs) * Ht *
                      get(
                        paste0("Ht", fit.delta, ".lambda", fit.delta2, ".", k2)
                      ) *
                      ((delta.obs == fit.delta2 &
                          get(time.var) == time.obs) -
                         exp(eps * Ht * get(
                           paste0("Ht", fit.delta,
                                  ".lambda", fit.delta2, ".", k2)
                         )) *
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

#### automatically true: run this one step
second.round <- FALSE

if (cr) {
  Pn.eic.fun <- function(mat) {
    eval <-
      eval.equation(mat, eps = 0, target.index = outcome.index[target])
    # if (target.S) {
    #   #-- if only want to target survival in competing risks setting
    #   eval <-
    #     sapply(tau, function(tt)
    #       - sum(sapply(eval, function(xx)
    #         xx[tau == tt])))
    # }
    return(eval)
  }
  
}

Pn.eic <- Pn.eic.fun(mat)
Pn.eic.norm.fun <- function(x2, x) {
  return(sqrt(sum(unlist(x2) * unlist(x))))
}

# useless ...?
Pn.eic2.fun <- function(Pn.eic)
  Pn.eic
Pn.eic2 <- Pn.eic2.fun(Pn.eic)

criterion <- 1 / (sqrt(n) * log(n))
Pn.eic.norm.prev <-
  Pn.eic.norm <- Pn.eic.norm.fun(Pn.eic2, Pn.eic)

init.fit <- lapply(target, function(each) {
  sapply(1:length(tau), function(kk) {
    mean(rowSums(sapply(a, function(aa)
      (2 * (
        aa == a[1]
      ) - 1) * (mat[get(A.name) == aa, 
                    get(paste0("F", estimation[[outcome.index[each]]][["event"]],
                               ".tau", kk))[1],
                    by = "id"][, 2][[1]]))))
  })
})

names(init.fit) <-
  paste0("F", sapply(outcome.index[target], function(each)
    estimation[[each]][["event"]]))

# insert init.list
init.ic <- eval.ic(mat, fit = init.fit, target.index = outcome.index[target])
init.list <-  lapply(1:length(init.fit), function(each.index) {
  out <- rbind(init.est = init.fit[[each.index]],
               init.se = init.ic[[each.index]])
  colnames(out) <- paste0("tau=", tau)
  return(out)
})
names(init.list) <- paste0("F", sapply(outcome.index[target], function(each)
  estimation[[each]]["event"]))

if (length(target) == length(outcome.index)) { ## what if not all events are targeted? ------------------------------------------------
  S.se.init <- sapply(tau, function(tt) {
    sqrt(
      mean(
        rowSums(
          sapply(1:length(outcome.index), function(target11) {
            unlist(
              eval.ic(mat,
                      unlist(lapply(init.list, function(xx) xx["init.est", paste0("tau=", tt)])),
                      target.index = outcome.index[target11],
                      tau.values = tt,
                      survival = TRUE
              )
            )
          })
        )^2) / n)
  })
  S.fit.init <-
    sapply(tau, function(tt)
      (treat.effect != "ate") - sum(sapply(init.list, function(fl)
        fl["init.est", paste0("tau=", tt)])))
  init.list$S <- rbind(init.est = S.fit.init, init.se = S.se.init)
  colnames(init.list$S) <- paste0("tau=", tau)
} 

tmle.list <- list(init = init.list)

init.ic <-
  eval.ic(mat, fit = init.fit, target.index = outcome.index[target])


for (step in 1:no.small.steps) {
  if (cr) {
    ## making duplicate columns for c-s hazards and clever covs
    for (each in outcome.index) {
      # loop over events
      fit.delta <- estimation[[each]][["event"]]
      mat[, (paste0("fit.cox", fit.delta, ".tmp")) := get(paste0("fit.cox", fit.delta))]
      for (kk in 1:length(tau)) {
        # loop over target times
        if (target.S)
          mat[, (paste0("Ht", ".lambda", fit.delta, ".", kk, ".tmp")) :=
                get(paste0("Ht", ".lambda", fit.delta, ".", kk))]
        for (each2 in outcome.index[target]) {
          # nested loop over events
          fit.delta2 <- estimation[[each2]][["event"]]
          mat[, (paste0(
            "Ht",
            fit.delta2,
            ".lambda",
            fit.delta,
            ".",
            kk,
            ".tmp"
          )) :=
            get(paste0("Ht", fit.delta2, ".lambda", fit.delta, ".", kk))]
        }
      }
    }
  } else {
    mat[, fit.cox1.tmp := fit.cox1]
    for (kk in 1:length(tau)) {
      mat[, (paste0("Ht.lambda.", kk, ".tmp")) := get((paste0("Ht.lambda.", kk)))]
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
    } else {
      ## 5.1 calculate update step direction? -------------------------------------
      if (!separate.cr) {
        for (each in outcome.index) {
          # loop over event hazards
          fit.delta <- estimation[[each]][["event"]]
          mat[, (paste0("delta", fit.delta, ".dx")) := 0] # fluctuation of cs-haz
          for (each2 in outcome.index[target]) {
            # loop over target events
            fit.delta2 <- estimation[[each2]][["event"]]
            for (kk in 1:length(tau)) {
              # loop over target times
              mat[, (paste0("delta", fit.delta, ".dx")) :=
                    get(paste0("delta", fit.delta, ".dx")) +
                    (get(time.var) <= tau[kk]) * Ht *
                    ((get(
                      paste0("Ht", fit.delta2, ".lambda", fit.delta, ".", kk)
                    )) *
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
                    ((get(
                      paste0("Ht", fit.delta2, ".lambda", fit.delta, ".", kk)
                    )) *
                      Pn.eic.separate[each2 == outcome.index[target]][[1]][each == outcome.index[target]][[1]][kk]) / Pn.eic.norm]
            }
          }
        }
      }
    }
  } else {
    mat[, delta.dx := 0]
    for (kk in 1:length(tau)) {
      mat[, delta.dx := delta.dx +
            (get(time.var) <= tau[kk]) *
            Ht * get(paste0("Ht.lambda.", kk)) * Pn.eic2[kk] / Pn.eic.norm]
    }
  }
  
  ## 5.2 set update epsilon down-scaling factor ----------------------------------------
  # why not just use deps.size?
  deps <- deps.size
  
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
  if (cr) {
    ## 5.3 update cause specific hazards ----------------------------------------
    mat[, surv.t := 1]
    for (each in outcome.index) {
      # loop over target events
      fit.delta <- estimation[[each]][["event"]]
      mat[, (paste0("fit.cox", fit.delta)) :=
            get(paste0("fit.cox", fit.delta)) *
            exp(deps * get(paste0("delta", fit.delta, ".dx")))]
      mat[get(paste0("fit.cox", fit.delta)) > 500, (paste0("fit.cox", fit.delta)) := 500]
      # Does changing the order of event updates matter? -------------------------------
      mat[, surv.t := surv.t * exp(-cumsum(get(paste0(
        "dhaz", fit.delta
      )) *
        get(paste0(
          "fit.cox", fit.delta
        )))),
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
              get(paste0("F", fit.delta, ".t"))[get(time.var) == max(get(time.var)[get(time.var) <= tau[kk]])], by = c("id", A.name)]
        for (each2 in outcome.index) {
          fit.delta2 <- estimation[[each2]][["event"]]
          if (fit.delta == fit.delta2) {
            mat[surv.t > 0, (paste0("Ht", fit.delta, ".lambda", fit.delta2, ".", kk)) :=
                  -(1 - (get(
                    paste0("F", fit.delta, ".tau", kk)
                  ) -
                    get(paste0(
                      "F", fit.delta, ".t"
                    ))) / surv.t)]
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
    if (target.S) {
      for (each2 in outcome.index) {
        fit.delta2 <- estimation[[each2]][["event"]]
        for (kk in 1:length(tau)) {
          mat[, (paste0("Ht", ".lambda", fit.delta2, ".", kk)) := 0]
          for (each in outcome.index) {
            fit.delta <- estimation[[each]][["event"]]
            mat[, (paste0("Ht", ".lambda", fit.delta2, ".", kk)) :=
                  get(paste0("Ht", ".lambda", fit.delta2, ".", kk)) -
                  get(paste0("Ht", fit.delta, ".lambda", fit.delta2, ".", kk))]
          }
        }
      }
    }
  } else {
    mat[, fit.cox1 := fit.cox1 * exp(deps * delta.dx)]
    mat[fit.cox1 > 500, fit.cox1 := 500]
    mat[, surv.t := exp(-cumsum(dhaz1 * fit.cox1)), by = c("id", A.name)]
    for (kk in 1:length(tau)) {
      mat[, (paste0("surv.tau", kk)) :=
            surv.t[get(time.var) == max(get(time.var)[get(time.var) <= tau[kk]])],
          by = c("id", A.name)]
      mat[surv.t > 0, (paste0("Ht.lambda.", kk)) := get(paste0("surv.tau", kk)) / surv.t]
      mat[surv.t == 0, (paste0("Ht.lambda.", kk)) := 1]
    }
  }
  
  Pn.eic <- Pn.eic.fun(mat)
  Pn.eic2 <- Pn.eic2.fun(Pn.eic)
  Pn.eic.norm <- Pn.eic.norm.fun(Pn.eic2, Pn.eic)
  
  
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
  
  # reset all columns to '.tmp'
  
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
                get(paste0(
                  "Ht",
                  fit.delta2,
                  ".lambda",
                  fit.delta,
                  ".",
                  kk,
                  ".tmp"
                ))]
        }
      }
    }
  }
  
  Pn.eic <- Pn.eic.fun(mat)
  Pn.eic2 <- Pn.eic2.fun(Pn.eic)
  
  Pn.eic.norm <- Pn.eic.norm.fun(Pn.eic2, Pn.eic)
  deps.size <- 0.5 * deps.size # 0.1 * deps.size
  
  
  
  
  Pn.eic3 <-
    lapply(1:length(Pn.eic), function(kk)
      Pn.eic[[kk]] / (ifelse(
        any(unlist(init.ic) == 0), init.ic[[kk]] + 0.001, init.ic[[kk]]
      ) * sqrt(n)))
  
  
  
  check.sup.norm <- max(abs(unlist(Pn.eic3))) <= (criterion)
  
  
  
  
  
  #-- 12c -- compute sd:
  final.fit <- lapply(outcome.index[target], function(each) {
    sapply(1:length(tau), function(kk) {
      mean(rowSums(sapply(a, function(aa)
        (2 * (
          aa == a[1]
        ) - 1) * (mat[get(A.name) == aa,
                      get(paste0("F", estimation[[each]][["event"]],
                                 ".tau", kk))[1],
                      by = "id"][, 2][[1]]))))
    })
  })
  
  
  names(final.fit) <-
    paste0("F", sapply(outcome.index[target], function(each)
      estimation[[each]][["event"]]))
  final.ic <-
    eval.ic(mat, final.fit, target.index = outcome.index[target])
  
  
  final.list <-
    lapply(1:length(final.fit), function(each.index) {
      out <- rbind(tmle.est = final.fit[[each.index]],
                   tmle.se = final.ic[[each.index]])
      colnames(out) <- paste0("tau=", tau)
      return(out)
    })
  
  names(final.list) <-
    paste0("F", sapply(outcome.index[target], function(each)
      estimation[[each]]["event"]))
  
  
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



# return(tmle.list)

