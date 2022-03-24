
my_init_sl_fit <- function (T_tilde, Delta, A, W, t_max, sl_failure = c("SL.glm"),
                            sl_censoring = c("SL.glm"), sl_treatment = c("SL.glm"),
                            cv.Control = NULL, gtol = 0.001, cl = NULL, tt_fun = NULL)
{
    ftime <- T_tilde
    ftype <- Delta
    trt <- A
    adjustVars <- W
    t_0 <- t_max
    trtOfInterest <- 0:1
    SL.ftime <- sl_failure
    SL.ctime <- sl_censoring
    SL.trt <- sl_treatment
    adjustVars <- data.frame(adjustVars)
    ftypeOfInterest <- unique(ftype)
    n <- length(ftime)
    id <- seq_len(n)
    dat <- data.frame(id = id, ftime = ftime, ftype = ftype, trt = trt)
    if (!is.null(adjustVars))
        dat <- cbind(dat, adjustVars)
    nJ <- length(ftypeOfInterest)
    allJ <- sort(unique(ftype[ftype != 0]))
    ofInterestJ <- sort(ftypeOfInterest)
    ntrt <- length(trtOfInterest)
    uniqtrt <- sort(trtOfInterest)
    message("estimating treatment")
    trtOut <- my_estimate_treatment(dat = dat, ntrt = ntrt, cv.Control = cv.Control,
                                    uniqtrt = uniqtrt, adjustVars = adjustVars, SL.trt = SL.trt,
                                    verbose = TRUE, returnModels = TRUE, gtol = gtol)
    dat <- trtOut$dat
    trtMod <- trtOut$trtMod
    dataList <- survtmle::makeDataList(dat = dat, J = allJ, ntrt = ntrt,
                                       uniqtrt = uniqtrt, t0 = t_0, bounds = NULL)

    # insert time-varying covariates ------------------------------------------
    if (!is.null(tt_fun)) {
        dataList <- tryCatch({
            lapply(c("obs", "0", "1"), function(i) {
                tt_fun[["update_data"]](dataList[[i]])
            })
        }, error = function(cond) {
            message("tt_fun update_data error")
            dataList
        })
        names(dataList) <- c("obs", 0, 1)
        adjustVars <- tryCatch({
            tt_fun[["update_adjust"]](adjustVars)
        }, error = function(cond) {
            message("tt_fun update_adjust error")
            adjustVars
        })
    }

    # st = c(798, 1099)
    # dataList$obs <- dplyr::mutate(dataList$obs, st1 = t*7 > st[1], st2 = t*7 > st[2])
    # dataList$'0' <- dplyr::mutate(dataList$'0', st1 = t*7 > st[1], st2 = t*7 > st[2])
    # dataList$'1' <- dplyr::mutate(dataList$'1', st1 = t*7 > st[1], st2 = t*7 > st[2])
    # adjustVars <- cbind(adjustVars, st1 = 1, st2 = 1)

    # back to init_sl_fit -----------------------------------------------------


    message("estimating censoring")
    censOut <- tryCatch({
        my_estimate_censoring(dataList = dataList, ntrt = ntrt, cv.Control = cv.Control,
                              uniqtrt = uniqtrt, t0 = t_0, verbose = TRUE, adjustVars = adjustVars,
                              SL.ctime = SL.ctime, glm.family = "binomial", returnModels = TRUE,
                              gtol = gtol)
    }, error = function(cond) {
        message("censoring sl error")
        NULL
    })
    if (is.null(censOut)) {
        censOut <- list()
        censOut$dataList <- dataList
        censOut$dataList$obs[, "G_dC"] <- 1
        censOut$dataList$"0"[, "G_dC"] <- 1
        censOut$dataList$"1"[, "G_dC"] <- 1
        is_sl_censoring_converge <- FALSE
        dataList <- censOut$dataList
        ctimeMod <- NULL
    }
    else {
        dataList <- censOut$dataList
        ctimeMod <- censOut$ctimeMod
        is_sl_censoring_converge <- TRUE
    }
    message("estimating hazards")
    estOut <- my_estimate_hazards(dataList = dataList, cv.Control = cv.Control,
                                  J = allJ, verbose = TRUE, bounds = NULL,
                                  adjustVars = adjustVars, returnModels = TRUE,
                                  SL.ftime = SL.ftime, glm.family = "binomial",
                                  cl = cl)
    dataList <- estOut$dataList
    ftimeMod <- estOut$ftimeMod
    suppressWarnings(if (all(dataList[[1]] == "convergence failure")) {
        return("estimation convergence failure")
    })

    # del_fitLibrary <- function(x) {
    #     if ("env" %in% names(x)) {
    #         x$fitLibrary <- NULL
    #         x$env <- NULL
    #     } else {
    #         for (i in seq_along(x)) {
    #             value <- x[[i]]
    #             if (class(x[[i]]) %in% c("survtmle" ,"SuperLearner", "list")) {
    #                 if ("env" %in% names(x[[i]])) {
    #                     x[[i]]$fitLibrary <- NULL
    #                     x[[i]]$env <- NULL
    #                 } else {
    #                     x[[i]] <- del_fitLibrary(value)
    #                 }
    #             }
    #         }
    #     }
    #
    #     x
    # }
    # ftimeMod <- del_fitLibrary(ftimeMod)
    # ctimeMod <- del_fitLibrary(ctimeMod)
    # trtMod <- del_fitLibrary(trtMod)

    g_1 <- dat$g_1
    g_0 <- dat$g_0
    d1 <- dataList$`1`
    d0 <- dataList$`0`
    haz1 <- d1[, c("id", "t", "Q1Haz")]
    haz1 <- tidyr::spread(haz1, t, Q1Haz)
    haz1$id <- NULL
    haz0 <- d0[, c("id", "t", "Q1Haz")]
    haz0 <- tidyr::spread(haz0, t, Q1Haz)
    haz0$id <- NULL
    S_Ac_1 <- d1[, c("id", "t", "G_dC")]
    S_Ac_1 <- tidyr::spread(S_Ac_1, t, G_dC)
    S_Ac_1 <- S_Ac_1[, -1]
    S_Ac_0 <- d0[, c("id", "t", "G_dC")]
    S_Ac_0 <- tidyr::spread(S_Ac_0, t, G_dC)
    S_Ac_0 <- S_Ac_0[, -1]
    G_dC <- lapply(dataList, function(d) d[, c("id", "t", "trt", "G_dC")])
    density_failure_1 <- survival_curve$new(t = seq(range(ftime)[1],
                                                    range(ftime)[2]), hazard = haz1)
    density_failure_0 <- survival_curve$new(t = seq(range(ftime)[1],
                                                    range(ftime)[2]), hazard = haz0)
    density_censor_1 <- survival_curve$new(t = seq(range(ftime)[1],
                                                   range(ftime)[2]), survival = S_Ac_1)
    density_censor_0 <- survival_curve$new(t = seq(range(ftime)[1],
                                                   range(ftime)[2]), survival = S_Ac_0)
    return(list(density_failure_1 = density_failure_1, density_failure_0 = density_failure_0,
                density_censor_1 = density_censor_1, density_censor_0 = density_censor_0,
                g1W = g_1[, 1], models = list("A" = trtMod, "C" = ctimeMod,
                                              "Y" = ftimeMod),
                G_dC = G_dC))
}
environment(my_init_sl_fit) <- asNamespace('MOSS')


sl_glm <- function (Y, X, newX, family, obsWeights, model = TRUE, ...)
{
    if (is.matrix(X)) {
        X = as.data.frame(X)
    }
    if (is.matrix(newX)) {
        newX = as.data.frame(newX)
    }
    fit.glm <- glm(Y ~ factor(t) + ., data = X, family = family, weights = obsWeights,
                   model = model)
    pred <- predict(fit.glm, newdata = newX, type = "response")
    fit <- list(object = fit.glm)
    class(fit) <- "SL.glm"
    out <- list(pred = pred, fit = fit)
    return(out)
}
environment(sl_glm) <- asNamespace("SuperLearner")


sl_xgboost <- function (Y, X, newX, family, obsWeights, id, ntrees = 1500,
                        max_depth = 4, shrinkage = 0.1, minobspernode = 10, params = list(),
                        nthread = 8, verbose = FALSE, save_period = NULL, ...)
{
    SuperLearner:::.SL.require("xgboost")
    if (packageVersion("xgboost") < 0.6)
        stop("SL.xgboost requires xgboost version >= 0.6, try help('SL.xgboost') for details")
    if (!is.matrix(X)) {
        X = model.matrix(~. - 1, X)
    }
    if (!is.matrix(newX)) {
        newX = model.matrix(~. - 1, newX)
    }
    xgmat = xgboost::xgb.DMatrix(data = X, label = Y, weight = obsWeights)
    if (family$family == "gaussian") {
        model = xgboost::xgboost(data = xgmat, objective = "reg:linear",
                                 nrounds = ntrees, max_depth = max_depth,
                                 min_child_weight = minobspernode,
                                 eta = shrinkage, verbose = verbose, nthread = nthread,
                                 params = params, save_period = save_period)
    }
    if (family$family == "binomial") {
        model = xgboost::xgboost(data = xgmat, objective = "binary:logistic",
                                 nrounds = ntrees, max_depth = max_depth,
                                 min_child_weight = minobspernode,
                                 eta = shrinkage, verbose = verbose, nthread = nthread,
                                 params = params, save_period = save_period)
    }
    if (family$family == "multinomial") {
        model = xgboost::xgboost(data = xgmat, objective = "multi:softmax",
                                 nrounds = ntrees, max_depth = max_depth,
                                 min_child_weight = minobspernode,
                                 eta = shrinkage, verbose = verbose, num_class = length(unique(Y)),
                                 nthread = nthread, params = params, save_period = save_period)
    }
    pred = predict(model, newdata = newX)
    fit = list(object = model)
    class(fit) = c("SL.xgboost")
    out = list(pred = pred, fit = fit)
    return(out)
}
environment(sl_xgboost) <- asNamespace("SuperLearner")

sl_bayesglm <- function (Y, X, newX, family, obsWeights, ...)
{
    .SL.require("arm")
    if (is.matrix(X)) {
        X = as.data.frame(X)
    }
    if (is.matrix(newX)) {
        newX = as.data.frame(newX)
    }
    fit.glm <- arm::bayesglm(Y ~ factor(t) + ., data = X, family = family,
                             weights = obsWeights)
    pred <- predict(fit.glm, newdata = newX, type = "response")
    fit <- list(object = fit.glm)
    out <- list(pred = pred, fit = fit)
    class(out$fit) <- c("SL.bayesglm")
    return(out)
}
environment(sl_bayesglm) <- asNamespace("SuperLearner")


my_estimate_treatment <- function (dat, adjustVars, glm.trt = NULL, SL.trt = NULL, returnModels = FALSE,
                                   verbose = FALSE, gtol = 0.001, cv.Control = NULL, ...)
{
    if (length(unique(dat$trt)) == 1) {
        eval(parse(text = paste0("dat$g_", unique(dat$trt), "<- 1")))
    }
    else {
        thisY <- as.numeric(dat$trt == max(dat$trt))
        if (!is.null(SL.trt)) {
            if (class(SL.trt) != "SuperLearner") {
                trtMod <- SuperLearner::SuperLearner(Y = thisY, cvControl = cv.Control,
                                                     X = adjustVars, newX = adjustVars,
                                                     SL.library = SL.trt,id = dat$id,
                                                     verbose = verbose, family = "binomial")
            }
            else {
                trtMod <- SL.trt
            }
            dat[[paste0("g_", max(dat$trt))]] <- trtMod$SL.predict
            dat[[paste0("g_", min(dat$trt))]] <- 1 - trtMod$SL.predict
        }
        else if (!is.null(glm.trt) & is.null(SL.trt)) {
            trt_form <- paste("thisY", "~", glm.trt, sep = " ")
            trt_data_in <- as.data.frame(cbind(adjustVars, thisY))
            if (!("glm" %in% class(glm.trt)) & !("speedglm" %in%
                                                 class(glm.trt))) {
                trtMod <- fast_glm(reg_form = stats::as.formula(trt_form),
                                   data = trt_data_in, family = stats::binomial())
            }
            else {
                trtMod <- glm.trt
            }
            suppressWarnings(pred <- predict(trtMod, newdata = trt_data_in,
                                             type = "response"))
            dat[[paste0("g_", max(dat$trt))]] <- pred
            dat[[paste0("g_", min(dat$trt))]] <- 1 - pred
        }
    }
    eval(parse(text = paste0("dat$g_", min(dat$trt), "[dat$g_",
                             min(dat$trt), "< gtol]<- gtol")))
    eval(parse(text = paste0("dat$g_", max(dat$trt), "[dat$g_",
                             max(dat$trt), "< gtol]<- gtol")))
    out <- list()
    out$dat <- dat
    out$trtMod <- NULL
    if (returnModels)
        out$trtMod <- trtMod
    return(out)
}
environment(my_estimate_treatment) <- asNamespace("survtmle")


correct_glm <- function (Y, X, newX, family, obsWeights, model = TRUE, ...)
{
    if (is.matrix(X)) {
        X = as.data.frame(X)
    }
    fit.glm <- glm(Y ~ SMOKER + WA:st1 + AT1 + st2, data = X,
                   family = family, weights = obsWeights, model = model)
    if (is.matrix(newX)) {
        newX = as.data.frame(newX)
    }
    pred <- predict(fit.glm, newdata = newX, type = "response")
    fit <- list(object = fit.glm)
    class(fit) <- "SL.glm"
    out <- list(pred = pred, fit = fit)
    return(out)
}
environment(correct_glm) <- asNamespace('SuperLearner')


misspec_glm <- function (Y, X, newX, family, obsWeights, model = TRUE, ...)
{
    if (is.matrix(X)) {
        X = as.data.frame(X)
    }
    fit.glm <- glm(Y ~ SMOKER + AGE + BMIBL, data = X,
                   family = family, weights = obsWeights, model = model)
    if (is.matrix(newX)) {
        newX = as.data.frame(newX)
    }
    pred <- predict(fit.glm, newdata = newX, type = "response")
    fit <- list(object = fit.glm)
    class(fit) <- "SL.glm"
    out <- list(pred = pred, fit = fit)
    return(out)
}
environment(misspec_glm) <- asNamespace('SuperLearner')


get_sample_hazards <- function(obs, t_max, haz_fn) {
    obs <- dplyr::select(obs, SMOKER) %>% mutate(W = SMOKER) %>%
        dplyr::select(W)

    hazards <- lapply(0:1, function(arm) {
        h <- do.call(rbind, lapply(sort(unique(obs$W)), function(w) {
            sim_haz_fn(t = 1:t_max, A = arm, W = w)
        }))
        rownames(h) <- sort(unique(obs$W))
        return(h)
    })

    obs_haz <- lapply(0:1, function(a) {
        hazards[[a+1]][match(obs$W, sort(unique(obs$W))), ]
    })

    obs_surv <- lapply(0:1, function(a) {
        do.call(rbind, lapply(1:nrow(obs), function(i) {
            hazard_here <- head(c(0, obs_haz[[a+1]][i, ]), -1)
            cumprod(1 - hazard_here)
        }))
    })
    return(list(hazards = obs_haz, survivals = obs_surv))
}

simple_misspec <- function (Y, X, newX, family, obsWeights, model = TRUE, ...)
{
    if (is.matrix(X)) {
        X = as.data.frame(X)
    }
    fit.glm <- glm(Y ~ trt, data = X,
                   family = family, weights = obsWeights, model = model)
    if (is.matrix(newX)) {
        newX = as.data.frame(newX)
    }
    pred <- predict(fit.glm, newdata = newX, type = "response")
    fit <- list(object = fit.glm)
    class(fit) <- "SL.glm"
    out <- list(pred = pred, fit = fit)
    return(out)
}
environment(simple_misspec) <- asNamespace('SuperLearner')

simple_correct <- function (Y, X, newX, family, obsWeights, model = TRUE, ...)
{
    if (is.matrix(X)) {
        X = as.data.frame(X)
    }
    fit.glm <- glm(Y ~ as.logical(trt):OBESEBL + trt:SMOKER + SMOKER, data = X,
                   family = family, weights = obsWeights, model = model)
    if (is.matrix(newX)) {
        newX = as.data.frame(newX)
    }
    pred <- predict(fit.glm, newdata = newX, type = "response")
    fit <- list(object = fit.glm)
    class(fit) <- "SL.glm"
    out <- list(pred = pred, fit = fit)
    return(out)
}
environment(simple_correct) <- asNamespace('SuperLearner')


correct_infcens <- function (Y, X, newX, family, obsWeights, model = TRUE, ...)
{
    if (is.matrix(X)) {
        X = as.data.frame(X)
    }
    fit.glm <- glm(Y ~ trt:SMOKER + SMOKER, data = X,
                   family = family, weights = obsWeights, model = model)
    if (is.matrix(newX)) {
        newX = as.data.frame(newX)
    }
    pred <- predict(fit.glm, newdata = newX, type = "response")
    fit <- list(object = fit.glm)
    class(fit) <- "SL.glm"
    out <- list(pred = pred, fit = fit)
    return(out)
}
environment(correct_infcens) <- asNamespace('SuperLearner')


upscale_haz <- function(surv_curve, timescale) {
    h_init <- surv_curve$hazard
    s_init <- do.call(rbind, lapply(1:nrow(h_init), function(i) {
        cumprod(1 - c(0, h_init[i, ]))
    }))
    s_diff <- t(diff(t(s_init)))
    s_diff_full <- matrix(0, nrow = nrow(s_diff), ncol = ncol(s_diff)*timescale)
    for (i in 1:ncol(s_diff)) {
        s_diff_full[, seq(from = timescale*(i-1)+1, by = 1,
                          length.out = timescale)] <- s_diff[, i] / timescale
    }
    s_init_full <- t(apply(cbind(1, s_diff_full), 1, cumsum))
    dens_full <- survival_curve$new(t = 1:ncol(s_init_full),
                                    survival = s_init_full)
    dens_full$surv2haz()
    dens_full$hazard <- dens_full$hazard[, -ncol(dens_full$hazard)]
    dens_full$survival <- dens_full$survival[, -ncol(dens_full$survival)]
    return(dens_full)
}
environment(upscale_haz) <- asNamespace("MOSS")


cens_glm <- function (Y, X, newX, family, obsWeights, model = TRUE, ...)
{
    if (is.matrix(X)) {
        X = as.data.frame(X)
    }
    fit.glm <- glm(Y ~ factor(t) + I(t<12):as.logical(trt)*SMOKER, data = X,
                   family = family, weights = obsWeights, model = model)
    if (is.matrix(newX)) {
        newX = as.data.frame(newX)
    }
    pred <- predict(fit.glm, newdata = newX, type = "response")
    fit <- list(object = fit.glm)
    class(fit) <- "SL.glm"
    out <- list(pred = pred, fit = fit)
    return(out)
}
environment(cens_glm) <- asNamespace('SuperLearner')


screen_glmnet <- function (Y, X, family, alpha = 1, minscreen = 2, nfolds = 10,
                           nlambda = 100, ...)
{
    .SL.require("glmnet")
    if (!is.matrix(X)) {
        col_names <- colnames(X)
        X <- model.matrix(~-1 + ., X)
    }
    fitCV <- glmnet::cv.glmnet(x = X, y = Y, lambda = NULL, type.measure = "deviance",
                               nfolds = nfolds, family = family$family, alpha = alpha,
                               nlambda = nlambda)

    whichVariable <- (as.numeric(coef(fitCV$glmnet.fit, s = fitCV$lambda.min))[-1] !=
                          0)
    if (sum(whichVariable) < minscreen) {
        warning("fewer than minscreen variables passed the glmnet screen, increased lambda to allow minscreen variables")
        sumCoef <- apply(as.matrix(fitCV$glmnet.fit$beta), 2,
                         function(x) sum((x != 0)))
        newCut <- which.max(sumCoef >= minscreen)
        whichVariable <- (as.matrix(fitCV$glmnet.fit$beta)[,
                                                           newCut] != 0)
    }
    if (exists('col_names')) {
        whichVariable <- as.logical(sapply(col_names, function(col) {
            max(str_detect(names(whichVariable[whichVariable == T]), col))}))
    }
    return(whichVariable)
}
environment(screen_glmnet) <- asNamespace("SuperLearner")

