library(survival)
library(survtmle)
library(data.table)
library(tmle3)
library(sl3)
data(pbc,package="survival")
setDT(pbc)
d <- pbc[!is.na(trt)]
d[,trt:=trt-1]
d[,ctime:=time]
d[,time:=round(time/365.25)+1]
d[,event:=1*(status!=0)]

x <- with(d,
          survtmle(ftime = time,
                   ftype = event,
                   trt = trt,
                   t0 = 5, adjustVars = data.frame(age,sex),
                   glm.ftime="trt+age+sex",
                   glm.ctime="trt+age+sex",
                   SL.ftime = NULL, SL.ctime = NULL, SL.trt = NULL,
                   glm.trt = "age+sex",
                   returnIC = TRUE, returnModels = TRUE,
                   ftypeOfInterest = 1, trtOfInterest = c(1, 0),
                   maxIter = 25, method = "hazard"))

y <- myTMLE3(data=d,A="trt",W=c("age","sex"),TIME="time",EVENT="event",tau=5)
d[,treatment:=factor(trt)]
z <- riskRegression::ate(event=Surv(ctime,event)~treatment+age+sex,
                         treatment=treatment~age+sex,
                         censor=Surv(ctime,event==0)~treatment+age+sex,
                         times=5*365.25,
                         data=d)
source(file = "~/research/SoftWare/test-tmle-surival/ConCR-TMLE/R/contmle.R")
library(zoo)
library(nleqslv)
u <- contmle(d,target = 1,iterative = FALSE,treat.effect = 0,tau = 5*365.25,estimation = list("cause1" = list(fit = "cox",model = Surv(ctime, event == 1) ~ trt + sex + age),"cens" = list(fit = "cox",model = Surv(ctime, event == 0) ~ trt + sex + age)),treat.model = trt ~ sex + age,sl.models=NULL)


source(file = "~/research/SoftWare/test-tmle-surival/ConCR-TMLE/R/functions/doConCRTmle.R")
logreg <- make_learner(Lrnr_glm)
Models <- list("A" = logreg,
               "0" = list(mod1 = Surv(Time, Event == 0) ~ Trt+age+sex,
                          mod2 = Surv(Time, Event == 0) ~ Trt + age),
               "1" = list(mod1 = Surv(Time, Event == 1) ~ Trt,
                          mod2 = Surv(Time, Event == 1) ~ Trt+age+sex))

v <- with(d,doConCRTmle(EventTime=ctime, EventType=event, Treatment=trt, CovDataTable=d[,.(age,sex)], CovTrtTime = NULL,
                        ID = NULL, TargetTimes = 5*365.25,
                        TargetEvents = 1, Models=Models, CVArgs = NULL, NumUpdateSteps = 25,
                        OneStepEps = 0.1, PropScoreCutoff = 0.05))

myTMLE3 <- function(data,A,W,TIME,EVENT,tau){
    k_grid <- 1:max(data[[TIME]])
    all_times <- lapply(k_grid, function(t_current) {
        df_time <- copy(data)
        df_time$N <- as.numeric(t_current == data[[TIME]] & data[[EVENT]] == 1)
        df_time$A_c <- as.numeric(t_current == data[[TIME]] & data[[EVENT]] == 0)
        df_time$pre_failure <- as.numeric(t_current <= data[[TIME]])
        df_time$t <- t_current
        df_time$X <- 1:nrow(df_time)
        return(df_time)
    })
    df_long <- rbindlist(all_times)
    node_list <- list(
        W=W,
        A =A,
        T_tilde = TIME,
        Delta = EVENT,
        time = "t",
        N = "N",
        A_c = "A_c",
        id = "X",
        pre_failure = "pre_failure"
    )
    learners <- list(
        ## glmnet = make_learner(Lrnr_glmnet),
        glm = make_learner(Lrnr_glm)
        ## gam = make_learner(Lrnr_gam),
        ## rf = make_learner(Lrnr_ranger)
    )
    lrnr_glm <- make_learner(Lrnr_glm)
    sl_A <- Lrnr_sl$new(lrnr_glm)
    sl_Y <- Lrnr_sl$new(learners)
    learner_list <- list(A = sl_A, N = sl_Y, A_c = sl_Y)
    # FIXME
    var_types <- list(T_tilde = Variable_Type$new("continuous"),
                      t = Variable_Type$new("continuous"),
                      Delta = Variable_Type$new("binomial"))

    survival_spec1 <- tmle_survival(
        treatment_level = 1,
        control_level = 0,
        variable_types = var_types,
        target_times = tau
    )
    survival_spec0 <- tmle_survival(
        treatment_level = 0,
        control_level = 1,
        variable_types = var_types,
        target_times = tau
    )

    tmle_task1 <- survival_spec1$make_tmle_task(df_long, node_list)
    tmle_task0 <- survival_spec0$make_tmle_task(df_long, node_list)

    initial_likelihood1 <- survival_spec1$make_initial_likelihood(tmle_task1, learner_list)
    initial_likelihood0 <- survival_spec0$make_initial_likelihood(tmle_task0, learner_list)

    up1 <- tmle3_Update_survival$new(
                                     maxit = 25,
                                     cvtmle = TRUE,
                                     convergence_type = "scaled_var",
                                     delta_epsilon = 1e-2,
                                     fit_method = "l2",
                                     use_best = TRUE,
                                     verbose = TRUE
                                 )
    up0 <- tmle3_Update_survival$new(
                                     maxit = 25,
                                     cvtmle = TRUE,
                                     convergence_type = "scaled_var",
                                     delta_epsilon = 1e-2,
                                     fit_method = "l2",
                                     use_best = TRUE,
                                     verbose = TRUE
                                 )

    targeted_likelihood1 <- Targeted_Likelihood$new(initial_likelihood1, updater = up1)
    targeted_likelihood0 <- Targeted_Likelihood$new(initial_likelihood0, updater = up0)

    tmle_params1 <- survival_spec1$make_params(tmle_task1, targeted_likelihood1)
    tmle_params0 <- survival_spec0$make_params(tmle_task0, targeted_likelihood0)

    # max(abs(colMeans(tmle_params1[[1]]$estimates(tmle_task, "validation")$IC[, 1:10])))
    # debugonce(tmle_params[[1]]$estimates)

    tmle_fit_manual1 <- fit_tmle3(
        tmle_task1, targeted_likelihood1, tmle_params1,
        targeted_likelihood1$updater
    )

    tmle_fit_manual0 <- fit_tmle3(
        tmle_task0, targeted_likelihood0, tmle_params0,
        targeted_likelihood0$updater
    )

    # conv <- apply(abs(do.call(rbind,up$EDs)),1,max)
    pos.tau <- prodlim::sindex(eval.times=tau,jump.times=k_grid)
    results <-         data.frame("Estimator" = "tmle3",
                   "t" = tau,
                   "s0" = tmle_fit_manual0$estimates[[1]]$psi[pos.tau],
                   "s1" = tmle_fit_manual1$estimates[[1]]$psi[pos.tau],
                   "se0" = sqrt(diag(var(tmle_fit_manual0$estimates[[1]]$IC)) /
                                nrow(data))[pos.tau],
                   "se1" = sqrt(diag(var(tmle_fit_manual1$estimates[[1]]$IC)) /
                                nrow(data))[pos.tau]) 
    return(results)
}
