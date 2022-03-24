
# survtmle ----------------------------------------------------------------

surv_tmle <- function (ftime, ftype, trt, adjustVars, t0 = max(ftime[ftype > 0]),
                       SL.ftime = NULL, SL.ctime = NULL, SL.trt = NULL, glm.ftime = NULL,
                       glm.ctime = NULL, glm.trt = NULL, returnIC = TRUE, maxIter = 10,
                       returnModels = TRUE, ftypeOfInterest = unique(ftype[ftype != 0]),
                       trtOfInterest = unique(trt), method = "hazard", bounds = NULL,
                       verbose = FALSE, tol = 1/(sqrt(length(ftime))),
                       Gcomp = FALSE, gtol = 0.001, targets = max(ftime[ftype > 0]),
                       tvcov_fun = NULL)
{
    call <- match.call(expand.dots = TRUE)
    clean <- checkInputs(ftime = ftime, ftype = ftype, trt = trt, verbose = verbose,
                         t0 = t0, adjustVars = adjustVars, SL.ftime = SL.ftime,
                         SL.ctime = SL.ctime, SL.trt = SL.trt, glm.ftime = glm.ftime,
                         glm.ctime = glm.ctime, glm.trt = glm.trt, returnIC = returnIC,
                         returnModels = returnModels, ftypeOfInterest = ftypeOfInterest,
                         trtOfInterest = trtOfInterest, bounds = bounds,
                         tol = tol, Gcomp = Gcomp, method = method)
    if (method == "hazard") {
        tmle.fit <- haz_tmle(ftime = clean$ftime, ftype = clean$ftype, trt = clean$trt,
                             targets = targets, adjustVars = clean$adjustVars,
                             SL.ftime = clean$SL.ftime, SL.ctime = clean$SL.ctime,
                             SL.trt = clean$SL.trt, glm.ftime = clean$glm.ftime,
                             glm.ctime = clean$glm.ctime, glm.trt = clean$glm.trt,
                             returnIC = returnIC, returnModels = returnModels,
                             ftypeOfInterest = ftypeOfInterest, verbose = verbose,
                             trtOfInterest = trtOfInterest, bounds = bounds,
                             tol = tol, maxIter = maxIter, gtol = gtol,
                             tvcov_fun = tvcov_fun)
    }
    else if (method == "mean") {
        tmle.fit <- surv_ltmle(ftime = clean$ftime, ftype = clean$ftype,
                               trt = clean$trt, t0 = t0, SL.ftime = clean$SL.ftime,
                               adjustVars = clean$adjustVars, SL.trt = clean$SL.trt,
                               SL.ctime = clean$SL.ctime, glm.trt = clean$glm.trt,
                               glm.ftime = clean$glm.ftime, returnIC = returnIC,
                               glm.ctime = clean$glm.ctime, bounds = bounds,
                               returnModels = returnModels, tol = tol, Gcomp = Gcomp,
                               ftypeOfInterest = ftypeOfInterest, verbose = verbose,
                               trtOfInterest = trtOfInterest, gtol = gtol,
                               tvcov_fun = tvcov_fun)
    }

    out <- list(call = call, est = tmle.fit$est, var = tmle.fit$var,
                meanIC = tmle.fit$meanIC, ic = tmle.fit$ic,
                ftimeMod = tmle.fit$ftimeMod, ctimeMod = tmle.fit$ctimeMod,
                trtMod = tmle.fit$trtMod, t0 = t0, ftime = tmle.fit$ftime,
                ftype = tmle.fit$ftype, trt = tmle.fit$trt,
                adjustVars = tmle.fit$adjustVars, init_est = tmle.fit$init_est)
    class(out) <- "survtmle"
    return(out)
}
environment(surv_tmle) <- asNamespace("survtmle")

# hazard_tmle ---------------------------------------------------------


haz_tmle <- function (ftime, ftype, trt, targets = max(ftime[ftype > 0]),
                      adjustVars = NULL, SL.ftime = NULL, SL.ctime = NULL,
                      SL.trt = NULL, glm.ftime = NULL, glm.ctime = NULL,
                      glm.trt = "1", glm.family = "binomial", returnIC = TRUE,
                      returnModels = FALSE, trtOfInterest = unique(trt),
                      ftypeOfInterest = unique(ftype[ftype != 0]),
                      bounds = NULL, verbose = FALSE, maxIter = 100,
                      tol = 1/(sqrt(length(ftime))),
                      gtol = 0.001, tvcov_fun = NULL, ...)
{
    n <- length(ftime)
    id <- seq_len(n)
    t0 <- max(targets)
    dat <- data.frame(id = id, ftime = ftime, ftype = ftype,
                      trt = trt)
    if (!is.null(adjustVars))
        dat <- cbind(dat, adjustVars)
    nJ <- length(ftypeOfInterest)
    allJ <- sort(unique(ftype[ftype != 0]))
    ofInterestJ <- sort(ftypeOfInterest)
    ntrt <- length(trtOfInterest)
    uniqtrt <- sort(trtOfInterest)
    trtOut <- estimateTreatment(dat = dat, ntrt = ntrt, uniqtrt = uniqtrt,
                                adjustVars = adjustVars, SL.trt = SL.trt,
                                glm.trt = glm.trt, returnModels = returnModels,
                                gtol = gtol)
    dat <- trtOut$dat
    trtMod <- trtOut$trtMod
    dataList <- makeDataList(dat = dat, J = allJ, ntrt = ntrt,
                             uniqtrt = uniqtrt, t0 = t0, bounds = bounds)

    # insert time-varying covariates ------------------------------------------
    if (!is.null(tvcov_fun)) {
        dataList <- tryCatch({
            lapply(c("obs", "0", "1"), function(i) {
                tvcov_fun[["update_data"]](dataList[[i]])
            })
        }, error = function(cond) {
            message("tvcov_fun update_data error")
            dataList
        })
        names(dataList) <- c("obs", 0, 1)
        adjustVars <- tryCatch({
            tvcov_fun[["update_adjust"]](adjustVars)
        }, error = function(cond) {
            message("tvcov_fun update_adjust error")
            adjustVars
        })
    }
    censOut <- my_estimate_censoring(dataList = dataList, ntrt = ntrt,
                                     uniqtrt = uniqtrt, t0 = t0, verbose = verbose,
                                     adjustVars = adjustVars, SL.ctime = SL.ctime,
                                     glm.ctime = glm.ctime, glm.family = glm.family,
                                     returnModels = returnModels, gtol = gtol)
    dataList <- censOut$dataList
    ctimeMod <- censOut$ctimeMod
    estOut <- my_estimate_hazards(dataList = dataList, J = allJ, verbose = verbose,
                                  bounds = bounds, adjustVars = adjustVars,
                                  glm.family = glm.family, glm.ftime = glm.ftime,
                                  returnModels = returnModels, SL.ftime = SL.ftime,
                                  ntrt = ntrt)
    dataList <- estOut$dataList
    ftimeMod <- estOut$ftimeMod
    suppressWarnings(if (all(dataList[[1]] == "convergence failure")) {
        return("estimation convergence failure")
    })

    dataList <- update_Vars(dataList = dataList, allJ = allJ, nJ = nJ, ntrt = ntrt,
                            ofInterestJ = ofInterestJ, uniqtrt = uniqtrt, t0 = t0,
                            targets = targets, verbose = verbose)
    dat <- get_haz_IC(dataList = dataList, dat = dat,
                      ofInterestJ = ofInterestJ, allJ = allJ, nJ = nJ,
                      uniqtrt = uniqtrt, ntrt = ntrt, verbose = verbose,
                      targets = targets, t0 = t0)
    infCurves <- dat[, grep("D.j", names(dat))]
    if (!is.numeric(infCurves)) {
        meanIC <- colMeans(infCurves)
        seIC <- sqrt(diag(crossprod(as.matrix(infCurves))/n^2))
    }
    else {
        meanIC <- mean(infCurves)
        seIC <- sqrt(as.vector(crossprod(as.matrix(infCurves))/n^2))
    }
    attr(dataList, "fluc") <- rep(Inf, ntrt * nJ^2)
    ct <- 0

    init_est <- rowNames <- NULL
    for (j in ofInterestJ) {
        for (z in uniqtrt) {
            init_psi_n <- colMeans(dat[paste0("margF", j, ".z", z, ".t", targets)])
            names(init_psi_n) <- paste0("t.", targets)
            init_est <- rbind(init_est, init_psi_n)
            rowNames <- c(rowNames, paste0("A = ", z, ", J = ", j, collapse = ""))
        }
    }
    row.names(init_est) <- rowNames

    while (any(abs(meanIC) > seIC*tol) & ct <= maxIter) {
        ct <- ct + 1
        dataList <- fluctuate_haz(dataList = dataList, ofInterestJ = ofInterestJ,
                                  tol = tol, allJ = allJ, nJ = nJ, uniqtrt = uniqtrt,
                                  ntrt = ntrt, targets = targets, t0 = t0,
                                  verbose = verbose)
        suppressWarnings(if (all(dataList[[1]] == "convergence failure")) {
            return("fluctuation convergence failure")
        })
        dat <- get_haz_IC(dataList = dataList, dat = dat, ofInterestJ = ofInterestJ,
                          allJ = allJ, nJ = nJ, uniqtrt = uniqtrt, ntrt = ntrt,
                          targets = targets, t0 = t0, verbose = verbose)
        infCurves <- dat[, grep("D.j", names(dat))]
        meanIC <- colMeans(infCurves)
        if (verbose) {
            cat("TMLE Iteration ", ct, "  : Pn(IC) = ", round(meanIC, 4), "\n")
        }
    }
    if (ct == maxIter + 1) {
        warning("TMLE fluctuations did not converge. Check that meanIC is adequately \n
                small and proceed with caution.")
    }
    est <- rowNames <- NULL
    for (j in ofInterestJ) {
        for (z in uniqtrt) {
            psi_n <- colMeans(dat[paste0("margF", j, ".z", z, ".t", targets)])
            names(psi_n) <- paste0("t.", targets)
            est <- rbind(est, psi_n)
            rowNames <- c(rowNames, paste0("A = ", z, ", J = ", j, collapse = ""))
        }
    }
    row.names(est) <- rowNames
    var <- crossprod(as.matrix(infCurves))/n^2
    if (!returnModels) {
        del_fitLibrary <- function(x) {
            if ("env" %in% names(x)) {
                x$fitLibrary <- NULL
                x$env <- NULL
            } else {
                for (i in seq_along(x)) {
                    value <- x[[i]]
                    if (class(x[[i]]) %in% c("survtmle" ,"SuperLearner", "list")) {
                        if ("env" %in% names(x[[i]])) {
                            x[[i]]$fitLibrary <- NULL
                            x[[i]]$env <- NULL
                        } else {
                            x[[i]] <- del_fitLibrary(value)
                        }
                    }
                }
            }

            x
        }
        ftimeMod <- del_fitLibrary(ftimeMod)
        ctimeMod <- del_fitLibrary(ctimeMod)
        trtMod <- del_fitLibrary(trtMod)
    }
    out <- list(est = est, var = var, meanIC = meanIC, ic = infCurves,
                trtMod = trtMod, ftimeMod = ftimeMod, ctimeMod = ctimeMod,
                ftime = ftime, ftype = ftype, trt = trt, adjustVars = adjustVars,
                init_est = init_est)
    class(out) <- "survtmle"
    return(out)
}
environment(haz_tmle) <- asNamespace("survtmle")

# update_Variables ----------------------------------------------------



update_Vars <- function (dataList, allJ, ofInterestJ, nJ, uniqtrt, ntrt, t0,
                         targets, verbose, ...)
{
    dataList[2:(ntrt + 1)] <- lapply(dataList[2:(ntrt + 1)],
                                     function(x, allJ)
                                     {
                                         Q_dot <- rowSums(cbind(rep(0, nrow(x)), x[, paste0("Q",allJ, "Haz")]))
                                         x$S.t <- unlist(by((1 - Q_dot), x$id, FUN = cumprod))
                                         x$S.t[x$S.t == 0] <- .Machine$double.neg.eps
                                         S.tminus1 <- c(1, x$S.t[1:(length(x$S.t) - 1)])
                                         S.tminus1[x$t == 1] <- 1
                                         for (j in ofInterestJ) {
                                             x[[paste0("F", j, ".t")]] <- unlist(by(x[, paste0("Q", j, "Haz")]
                                                                                    * S.tminus1, x$id,
                                                                                    FUN = cumsum))
                                         }
                                         x
                                     }, allJ = allJ)
    for (j in ofInterestJ) {
        Fj.t.allZ <- vector(mode = "list", length = ntrt)
        for (i in 1:ntrt) {
            for (targ in targets) {
                Fj.t.allZ[[i]][[as.character(targ)]] <-
                    dataList[[i + 1]][[paste0("F", j, ".t")]][
                        dataList[[i + 1]]$t == targ]
            }
        }
        dataList <- lapply(dataList, function(x, j, uniqtrt, Fj.t.allZ) {
            for (i in seq_along(uniqtrt)) {
                for (k in seq_along(targets)) {
                    ind <- tapply(X = x$id, INDEX = x$id, FUN = NULL)
                    x[[paste0("F", j, ".z", uniqtrt[i], ".t", targets[k])]] <-
                        Fj.t.allZ[[i]][[k]][ind]
                }
            }
            x
        }, j = j, uniqtrt = uniqtrt, Fj.t.allZ = Fj.t.allZ)
    }
    colInd <- which(colnames(dataList[[1]]) %in% c("S.t",
                                                   paste0("F", ofInterestJ, ".t")))
    if (ntrt == 1) {
        merge_vars <- c("id", "t")
    }
    else {
        merge_vars <- c("id", "t", "trt")
    }
    if (length(colInd) == 0) {
        dataList[[1]] <- merge(dataList[[1]],
                               Reduce(rbind, dataList[2:(ntrt + 1)])[,
                                                                     c("id", "t", "trt", "S.t",
                                                                       paste0("F", ofInterestJ, ".t"))],
                               by = merge_vars)
    }
    else {
        dataList[[1]] <-
            merge(x = dataList[[1]][, -colInd],
                  y = Reduce(rbind,
                             dataList[2:(ntrt + 1)])[, c("id", "t", "trt", "S.t",
                                                         paste0("F", ofInterestJ,
                                                                ".t"))],
                  by = merge_vars)
    }
    if (ntrt == 1) {
        dataList[[1]]$trt <- dataList[[1]]$trt.x
        dataList[[1]] <- subset(dataList[[1]], select = -c("trt.x", "trt.y"))
    }
    dataList <- lapply(dataList, function(x, allJ) {
        for (j in allJ) {
            if (length(allJ) > 1) {
                x[[paste0("hazNot", j)]] <- rowSums(cbind(rep(0, nrow(x)),
                                                          x[, paste0("Q",
                                                                     allJ[allJ != j],
                                                                     "Haz")]))
                x[[paste0("hazNot", j)]][x[[paste0("hazNot", j)]] == 1] <-
                    1 - .Machine$double.neg.eps
            }
            else {
                x[[paste0("hazNot", j)]] <- 0
            }
        }
        x
    }, allJ = allJ)
    dataList <- lapply(dataList, function(x, ofInterestJ, uniqtrt) {
        for (z in uniqtrt) {
            for (j in ofInterestJ) {
                for (targ in targets) {
                    x[[paste0("H", j, ".jSelf.z", z, ".t", targ)]] <-
                        (x$ftime >= x$t & x$trt == z) /
                        (x[[paste0("g_", z)]] * x$G_dC) * (1 - x[[paste0("hazNot", j)]]) *
                        ((x$t < targ) * (1 - (x[[paste0("F", j, ".z", z, ".t", targ)]] -
                                                  x[[paste0("F", j, ".t")]])/c(x$S.t)) +
                             as.numeric(x$t == targ))
                    x[[paste0("H", j, ".jNotSelf.z", z, ".t", targ)]] <-
                        -(x$ftime >= x$t & x$trt == z) /
                        (x[[paste0("g_", z)]] * x$G_dC) * (1 - x[[paste0("hazNot", j)]]) *
                        ((x$t < targ) * (x[[paste0("F", j, ".z", z, ".t", targ)]] -
                                             x[[paste0("F", j, ".t")]])/c(x$S.t))
                }
            }
        }
        x
    }, ofInterestJ = ofInterestJ, uniqtrt = uniqtrt)
    dataList
}
environment(haz_tmle) <- asNamespace("survtmle")



# getHazardInfluenceCurve ----------------------------------------------------------

get_haz_IC <- function (dataList, dat, allJ, ofInterestJ, nJ, uniqtrt, t0, targets,
                        verbose, ...)
{
    for (z in uniqtrt) {
        for (j in ofInterestJ) {
            for (targ in targets) {
                dat[[paste0("margF", j, ".z", z, ".t", targ)]] <-
                    mean(dataList[[1]][[paste0("F", j, ".z", z, ".t", targ)]][
                        dataList[[1]]$t == min(dataList[[1]]$t)])
                dat[[paste0("F", j, ".z", z, ".t", targ)]] <-
                    dataList[[1]][[paste0("F", j, ".z", z, ".t", targ)]][
                        dataList[[1]]$t == min(dataList[[1]]$t)]
                thisD <- NULL
                for (jTild in allJ) {
                    H <- paste0("H", j, ".j", ifelse(jTild == j, "Self", "NotSelf"),
                                ".z", z, ".t", targ)
                    thisD <- cbind(thisD, dataList[[1]][[H]] /
                                       (1 - dataList[[1]][[paste0("hazNot", j)]]) *
                                       (dataList[[1]][[paste0("N", jTild)]] -
                                            dataList[[1]][[paste0("Q", jTild, "Haz")]]))
                }
                dat[[paste0("D.j", j, ".z", z, ".t", targ)]] <-
                    unlist(by(rowSums(thisD), dataList[[1]]$id, FUN = sum)) +
                    dat[[paste0("F", j, ".z", z, ".t", targ)]] -
                    dat[[paste0("margF", j, ".z", z, ".t", targ)]]
            }
        }
    }
    return(dat)
}



# fluctuateHazards --------------------------------------------------------


fluctuate_haz <- function (dataList, allJ, ofInterestJ, nJ, uniqtrt, ntrt, t0,
                           targets, verbose, ...)
{
    eps <- NULL
    for (z in uniqtrt) {
        for (j in allJ) {
            cleverCovariatesNotSelf <- NULL
            if (length(ofInterestJ[ofInterestJ != j]) > 0) {
                cleverCovariatesNotSelf <- c(cleverCovariatesNotSelf,
                                             paste0("H", ofInterestJ[ofInterestJ != j],
                                                    ".jNotSelf.z", z, ".t", targets))
            }
            if (j %in% ofInterestJ) {
                cleverCovariatesSelf <- paste0("H", j, ".jSelf.z", z, ".t", targets)
            }
            else {
                cleverCovariatesSelf <- NULL
            }
            dataList <- lapply(dataList, function(x, j, allJ) {
                x$thisScale <- pmin(x[[paste0("u", j)]], 1 - x[[paste0("hazNot", j)]]) -
                    x[[paste0("l", j)]]
                x$thisOffset <- stats::qlogis(pmin((x[[paste0("Q", j, "Haz")]] -
                                                        x[[paste0("l", j)]]) / x$thisScale,
                                                   1 - .Machine$double.neg.eps))
                x$thisOutcome <- (x[[paste0("N",j)]] - x[[paste0("l",j)]]) / x$thisScale
                x
            }, j = j, allJ = allJ)
            fluc.mod <-
                stats::optim(par = rep(0, length(c(cleverCovariatesNotSelf,
                                                   cleverCovariatesSelf))),
                             fn = LogLikelihood_offset, Y = dataList[[1]]$thisOutcome,
                             H = suppressWarnings(as.matrix(
                                 Matrix::Diagonal(x = dataList[[1]]$thisScale) %*%
                                     as.matrix(dataList[[1]][, c(cleverCovariatesNotSelf,
                                                                 cleverCovariatesSelf)])
                             )),
                             offset = dataList[[1]]$thisOffset, method = "BFGS",
                             gr = grad_offset, control = list(reltol = 1e-7, maxit = 5e4))
            if (fluc.mod$convergence != 0) {
                warning("Fluctuation convergence failure. Using initial estimates.\n
                        Proceed with caution")
                beta <- rep(0, length(fluc.mod$par))
            }
            else {
                beta <- fluc.mod$par
            }
            eps <- c(eps, beta)
            dataList <- lapply(dataList, function(x, j) {
                x[[paste0("Q", j, "PseudoHaz")]][x$trt == z] <-
                    plogis(x$thisOffset[x$trt == z] +
                               suppressWarnings(as.matrix(Matrix::Diagonal(
                                   x = x$thisScale[x$trt == z]) %*%
                                       as.matrix(x[x$trt == z,
                                                   c(cleverCovariatesNotSelf,
                                                     cleverCovariatesSelf)])) %*%
                                       as.matrix(beta)))
                x[[paste0("Q", j, "Haz")]][x$trt == z] <-
                    x[[paste0("Q", j, "PseudoHaz")]][x$trt == z] *
                    x$thisScale[x$trt == z] + x[[paste0("l", j)]][x$trt == z]
                x
            }, j = j)
            dataList <- update_Vars(dataList = dataList, allJ = allJ,  nJ = nJ,
                                    ofInterestJ = ofInterestJ, uniqtrt = uniqtrt,
                                    ntrt = ntrt, targets = targets, t0 = t0,
                                    verbose = verbose)
        }
    }
    attr(dataList, "fluc") <- eps
    dataList
}
environment(fluctuate_haz) <- asNamespace("survtmle")




# fluctuate_both_hazards --------------------------------------------------

fluctuate_both_hazards <- function (dataList, allJ, ofInterestJ, nJ, uniqtrt, ntrt, t0,
                                    targets, verbose, ...)
{
    eps <- NULL
    for (z in uniqtrt) {
        for (j in allJ) {
            cleverCovariatesNotSelf <- NULL
            if (length(ofInterestJ[ofInterestJ != j]) > 0) {
                cleverCovariatesNotSelf <- c(cleverCovariatesNotSelf,
                                             paste0("H", ofInterestJ[ofInterestJ != j],
                                                    ".jNotSelf.z", z, ".t", targets))
            }
            if (j %in% ofInterestJ) {
                cleverCovariatesSelf <- paste0("H", j, ".jSelf.z", z, ".t", targets)
            }
            else {
                cleverCovariatesSelf <- NULL
            }
            dataList <- lapply(dataList, function(x, j, allJ) {
                x$thisScale <- pmin(x[[paste0("u", j)]], 1 - x[[paste0("hazNot", j)]]) -
                    x[[paste0("l", j)]]
                x$thisOffset <- stats::qlogis(pmin((x[[paste0("Q", j, "Haz")]] -
                                                        x[[paste0("l", j)]]) /
                                                       x$thisScale,
                                                   1 - .Machine$double.neg.eps))
                x$thisOutcome <- (x[[paste0("N",j)]] - x[[paste0("l",j)]]) /
                    x$thisScale
                x
            }, j = j, allJ = allJ)
            fluc.mod <-
                stats::optim(par = rep(0, length(c(cleverCovariatesNotSelf,
                                                   cleverCovariatesSelf))),
                             fn = LogLikelihood_offset, Y = dataList[[1]]$thisOutcome,
                             H = suppressWarnings(as.matrix(
                                 Matrix::Diagonal(x = dataList[[1]]$thisScale) %*%
                                     as.matrix(dataList[[1]][, c(cleverCovariatesNotSelf,
                                                                 cleverCovariatesSelf)])
                             )),
                             offset = dataList[[1]]$thisOffset, method = "BFGS",
                             gr = grad_offset, control = list(reltol = 1e-7, maxit = 5e4))
            if (fluc.mod$convergence != 0) {
                warning("Fluctuation convergence failure. Using initial estimates.\n
                        Proceed with caution")
                beta <- rep(0, length(fluc.mod$par))
            }
            else {
                beta <- fluc.mod$par
            }
            eps <- c(eps, beta)
            dataList <- lapply(dataList, function(x, j) {
                x[[paste0("Q", j, "PseudoHaz")]][x$trt == z] <-
                    plogis(x$thisOffset[x$trt == z] +
                               suppressWarnings(as.matrix(Matrix::Diagonal(
                                   x = x$thisScale[x$trt == z]) %*%
                                       as.matrix(x[x$trt == z,
                                                   c(cleverCovariatesNotSelf,
                                                     cleverCovariatesSelf)])) %*%
                                       as.matrix(beta)))
                x[[paste0("Q", j, "Haz")]][x$trt == z] <-
                    x[[paste0("Q", j, "PseudoHaz")]][x$trt == z] *
                    x$thisScale[x$trt == z] + x[[paste0("l", j)]][x$trt == z]
                x
            }, j = j)
            dataList <- update_Vars(dataList = dataList, allJ = allJ,  nJ = nJ,
                                    ofInterestJ = ofInterestJ, uniqtrt = uniqtrt,
                                    ntrt = ntrt, targets = targets, t0 = t0,
                                    verbose = verbose)
        }
    }
    attr(dataList, "fluc") <- eps
    dataList
}
environment(fluctuate_haz) <- asNamespace("survtmle")



# EstimateCensoring -------------------------------------------------------


my_estimate_censoring <- function (dataList, adjustVars, t0, SL.ctime = NULL,
                                   glm.ctime = NULL, glm.family,
                                   returnModels = FALSE, verbose = TRUE,
                                   gtol = 0.001, cv.Control = list(), cl = NULL,
                                   ...)
{
    include <- !(dataList[[1]]$t == dataList[[1]]$ftime & dataList[[1]]$C != 1 &
                     dataList[[1]]$t < t0) &
        !(dataList[[1]]$t == dataList[[1]]$ftime &
              dataList[[1]]$C == 1 & dataList[[1]]$t == t0)
    if (!is.null(glm.family)) {
        glm_family <- parse(text = paste0("stats::", glm.family,
                                          "()"))
    }
    if (is.null(SL.ctime)) {
        if (!(any(c("glm", "speedglm") %in% class(glm.ctime)))) {
            if (!all(dataList[[1]]$C == 0)) {
                if (glm.ctime == "Kaplan-Meier") {
                    dat <- dplyr::select(dataList[[1]], id, trt, ftime, ftype)
                    dat <- dplyr::distinct(dat)
                    ctimeMod <- survival::survfit(survival::Surv(
                        time = dat$ftime, event = as.numeric(dat$ftype == 0),
                        type = 'right') ~ 1, type = "kaplan-meier")
                    ctimeMod <- stepfun(ctimeMod$time, c(1, ctimeMod$surv))
                }
                else {
                    ctimeForm <- stats::as.formula(sprintf("%s ~ %s",
                                                           "C", glm.ctime))
                    ctimeMod <- fast_glm(reg_form = ctimeForm,
                                         data = dataList[[1]][include, ],
                                         family = eval(glm_family))
                    if (unique(class(ctimeMod) %in% c("glm", "lm"))) {
                        ctimeMod <- cleanglm(ctimeMod)
                    }
                }
            }
            else {
                dataList <- lapply(dataList, function(x) {
                    x$G_dC <- 1
                    x
                })
                ctimeMod <- "No censoring observed"
                class(ctimeMod) <- "noCens"
            }
        }
        else {
            ctimeMod <- glm.ctime
        }
        if (all(class(ctimeMod) != "noCens")) {
            dataList <- lapply(dataList, function(x) {
                if ("stepfun" %in% class(ctimeMod)) {
                    x$G_dC <- ctimeMod(x$t - 1)
                }
                else {
                    g_dC <- rep(1, length(x[, 1]))
                    if (t0 != 1) {
                        x$t <- x$t - 1
                        suppressWarnings(g_dC <- 1 - predict(ctimeMod, newdata = x,
                                                             type = "response"))
                        x$t <- x$t + 1
                        g_dC[x$t == 1] <- 1
                    }
                    x$G_dC <- as.numeric(unlist(by(g_dC, x$id, FUN = cumprod)))
                }
                x
            })
        }
        else {
            dataList <- lapply(dataList, function(x) {
                x$G_dC <- 1
                x
            })
        }
    }
    else {
        if (sum(names(SL.ctime) == c("obs", "0", "1")) == 3) {
            ctimeMod <- "SL.ctime estimated"
        }
        else if (class(SL.ctime) != "SuperLearner") {
            if (!all(dataList[[1]]$C == 0)) {
                if(inherits(cl, "cluster")) {
                    ctimeMod <- SuperLearner::snowSuperLearner(
                        cluster = cl, Y = dataList[[1]]$C[include],
                        X = dataList[[1]][include, c("t", "trt", names(adjustVars))],
                        id = dataList[[1]]$id[include], family = "binomial",
                        SL.library = SL.ctime, cvControl = cv.Control, verbose = verbose)
                }
                else {
                    ctimeMod <- SuperLearner::SuperLearner(
                        Y = dataList[[1]]$C[include],
                        X = dataList[[1]][include, c("t", "trt", names(adjustVars))],
                        id = dataList[[1]]$id[include], family = "binomial",
                        SL.library = SL.ctime,
                        cvControl = cv.Control, verbose = verbose)
                }
            }
            else {
                dataList <- lapply(dataList, function(x) {
                    x$G_dC <- 1
                    x
                })
                ctimeMod <- "No censoring observed"
                class(ctimeMod) <- "noCens"
            }
        }
        else {
            ctimeMod <- SL.ctime
        }
        if (sum(names(SL.ctime) == c("obs", "0", "1")) == 3) {
            for (name in names(SL.ctime)) {
                dataList[[name]] <- left_join(dataList[[name]],
                                              SL.ctime[[name]],
                                              by = c("id", "t", "trt"))
            }
        }
        else if (class(ctimeMod) != "noCens") {
            dataList <- lapply(dataList, function(x) {
                g_dC <- rep(1, nrow(x))
                if (t0 != 1) {
                    x$t <- x$t - 1
                    g_dC <- rep_along(x$t, 1)
                    t_include <- unique(dataList[[1]]$t[include])
                    g_dC[x$t %in% t_include] <-
                        suppressWarnings(1 - predict(ctimeMod,
                                                     newdata = x[x$t %in% t_include,
                                                                 c("t", "trt", names(adjustVars))],
                                                     onlySL = TRUE)[[1]])
                    x$t <- x$t + 1
                }
                x$G_dC <- as.numeric(unlist(by(g_dC, x$id, FUN = cumprod)))
                x
            })
        }
        else {
            dataList <- lapply(dataList, function(x) {
                x$G_dC <- 1
                x
            })
        }
    }
    dataList <- lapply(dataList, function(x) {
        x$G_dC[x$G_dC < gtol] <- gtol
        x
    })
    out <- list(dataList = dataList, ctimeMod = if (returnModels) {
        ctimeMod
    } else {
        NULL
    })
    return(out)
}
environment(my_estimate_censoring) <- asNamespace("survtmle")


# mean_tmle ---------------------------------------------------------------


surv_ltmle <- function (ftime, ftype, trt, t0 = max(ftime[ftype > 0]),
                        adjustVars = NULL, SL.ftime = NULL, SL.ctime = NULL,
                        SL.trt = NULL, glm.ftime = NULL, glm.ctime = NULL,
                        glm.trt = "1", glm.family = "binomial", returnIC = TRUE,
                        returnModels = FALSE, trtOfInterest = unique(trt),
                        ftypeOfInterest = unique(ftype[ftype != 0]),
                        bounds = NULL, verbose = FALSE, Gcomp = FALSE,
                        gtol = 0.001, tvcov_fun = NULL, ...)
{
    n <- length(ftime)
    id <- seq_len(n)
    dat <- data.frame(id = id, ftime = ftime, ftype = ftype,
                      trt = trt)
    if (!is.null(adjustVars)) {
        dat <- cbind(dat, adjustVars)
    }
    nJ <- length(ftypeOfInterest)
    allJ <- sort(unique(ftype[ftype != 0]))
    ofInterestJ <- sort(ftypeOfInterest)
    ntrt <- length(trtOfInterest)
    uniqtrt <- sort(trtOfInterest)
    trtOut <- estimateTreatment(dat = dat, ntrt = ntrt, uniqtrt = uniqtrt,
                                adjustVars = adjustVars, SL.trt = SL.trt, glm.trt = glm.trt,
                                returnModels = returnModels, gtol = gtol)
    dat <- trtOut$dat
    trtMod <- trtOut$trtMod
    dataList <- makeDataList(dat = dat, J = allJ, ntrt = ntrt,
                             uniqtrt = uniqtrt, t0 = t0, bounds = bounds)
    # insert time-varying covariates ------------------------------------------
    if (!is.null(tvcov_fun)) {
        dataList <- tryCatch({
            lapply(c("obs", "0", "1"), function(i) {
                tvcov_fun[["update_data"]](dataList[[i]])
            })
        }, error = function(cond) {
            message("tvcov_fun update_data error")
            dataList
        })
        names(dataList) <- c("obs", 0, 1)
        adjustVars <- tryCatch({
            tvcov_fun[["update_adjust"]](adjustVars)
        }, error = function(cond) {
            message("tvcov_fun update_adjust error")
            adjustVars
        })
    }
    censOut <- my_estimate_censoring(dataList = dataList, ntrt = ntrt,
                                     uniqtrt = uniqtrt, t0 = t0, verbose = verbose,
                                     adjustVars = adjustVars, SL.ctime = SL.ctime,
                                     glm.ctime = glm.ctime, glm.family = glm.family,
                                     returnModels = returnModels, gtol = gtol)
    dataList <- censOut$dataList
    ctimeMod <- censOut$ctimeMod
    wideDataList <- makeWideDataList(dat = dat, dataList = dataList,
                                     adjustVars = adjustVars, t0 = t0, allJ = allJ,
                                     ntrt = ntrt, uniqtrt = uniqtrt)
    timeAndType <- expand.grid(rev(seq_len(t0)), ofInterestJ)
    ftimeMod <- vector(mode = "list", length = length(ofInterestJ))
    names(ftimeMod) <- paste0("J", ofInterestJ)
    for (j in seq_along(ofInterestJ)) {
        ftimeMod[[j]] <- vector(mode = "list", length = t0)
        names(ftimeMod[[j]]) <- paste0("t", seq_len(t0))
    }
    for (i in seq_len(nrow(timeAndType))) {
        estOut <- estimateIteratedMean(wideDataList = wideDataList,
                                       t = timeAndType[i, 1], whichJ = timeAndType[i, 2],
                                       ntrt = ntrt, uniqtrt = uniqtrt, allJ = allJ, t0 = t0,
                                       SL.ftime = SL.ftime, adjustVars = adjustVars, glm.ftime = glm.ftime,
                                       verbose = verbose, returnModels = returnModels,
                                       bounds = bounds)
        wideDataList <- estOut$wideDataList
        eval(parse(text = paste0("ftimeMod$J", timeAndType[i,
                                                           2], "$t", timeAndType[i, 1], "<-estOut$ftimeMod")))
        wideDataList <- fluctuateIteratedMean(wideDataList = wideDataList,
                                              t = timeAndType[i, 1], whichJ = timeAndType[i, 2],
                                              ntrt = ntrt, uniqtrt = uniqtrt, allJ = allJ, t0 = t0,
                                              SL.ftime = SL.ftime, glm.ftime = glm.ftime, returnModels = returnModels,
                                              bounds = bounds, Gcomp = Gcomp)
    }
    est <- rowNames <- NULL
    for (j in ofInterestJ) {
        for (z in seq_along(uniqtrt)) {
            thisEst <- eval(parse(text = paste("mean(wideDataList[[",
                                               z + 1, "]]$Q", j, "star.1)", sep = "")))
            est <- rbind(est, thisEst)
            rowNames <- c(rowNames, paste(c(uniqtrt[z], j),
                                          collapse = " "))
            eval(parse(text = paste("wideDataList[[1]]$Q", j,
                                    "star.0.Z", uniqtrt[z], " <- rep(thisEst,n)",
                                    sep = "")))
            eval(parse(text = paste("wideDataList[[1]]$Q", j,
                                    "star.1.Z", uniqtrt[z], " <- wideDataList[[(z+1)]]$Q",
                                    j, "star.1", sep = "")))
        }
    }
    row.names(est) <- rowNames
    for (j in ofInterestJ) {
        for (z in seq_along(uniqtrt)) {
            for (t in rev(seq_len(t0))) {
                outcomeName <- ifelse(t == t0, paste("N", j,
                                                     ".", t0, sep = ""), paste("Q", j, "star.",
                                                                               t + 1, sep = ""))
                eval(parse(text = paste("wideDataList[[1]]$D.Z",
                                        uniqtrt[z], ".", j, "star.", t, " <- wideDataList[[1]]$H",
                                        uniqtrt[z], ".", t, "*(wideDataList[[1]][,outcomeName] - wideDataList[[1]]$Q",
                                        j, "star.", t, ")", sep = "")))
            }
            eval(parse(text = paste("wideDataList[[1]]$D.Z",
                                    uniqtrt[z], ".", j, "star.0 <- wideDataList[[1]]$Q",
                                    j, "star.1.Z", uniqtrt[z], " - wideDataList[[1]]$Q",
                                    j, "star.0.Z", uniqtrt[z], sep = "")))
            ind <- eval(parse(text = paste("grep('D.Z", uniqtrt[z],
                                           ".", j, "star', names(wideDataList[[1]]))",
                                           sep = "")))
            eval(parse(text = paste("wideDataList[[1]]$IC",
                                    j, "star.Z", uniqtrt[z], " <- rowSums(cbind(rep(0, nrow(wideDataList[[1]])),wideDataList[[1]][,ind]))",
                                    sep = "")))
        }
    }
    infCurves <- wideDataList[[1]][, grep("IC", names(wideDataList[[1]])),
                                   drop = FALSE]
    meanIC <- apply(infCurves, MARGIN = 2, FUN = mean)
    var <- t(as.matrix(infCurves)) %*% as.matrix(infCurves)/(n^2)
    row.names(var) <- colnames(var) <- rowNames
    out <- list(est = est, var = var, meanIC = meanIC, ic = infCurves,
                trtMod = trtMod, ftimeMod = ftimeMod, ctimeMod = ctimeMod,
                ftime = ftime, ftype = ftype, trt = trt, adjustVars = adjustVars)
    class(out) <- "survtmle"
    return(out)
}
environment(surv_ltmle) <- asNamespace("survtmle")

# my estimate hazards ------------------------------------------------------------
my_estimate_hazards <- function (dataList, J, adjustVars, SL.ftime = NULL,
                                 glm.ftime = NULL, glm.family, returnModels,
                                 bounds, verbose, cv.Control = list(), cl = NULL,
                                 timescale = 1, ...)
{
    ftimeMod <- vector(mode = "list", length = length(J))
    names(ftimeMod) <- paste0("J", J)
    if (!is.null(glm.family)) {
        glm_family <- parse(text = paste0("stats::", glm.family,
                                          "()"))
    }
    if (is.null(SL.ftime)) {
        if (is.null(bounds)) {
            for (j in J) {
                Qj_form <- sprintf("%s ~ %s", paste("N", j, sep = ""),
                                   glm.ftime)
                NlessthanJ <- rep(0, nrow(dataList[[1]]))
                for (i in J[J < j]) {
                    NlessthanJ <- NlessthanJ + dataList[[1]][[paste0("N",
                                                                     i)]]
                }
                if (!("glm" %in% class(glm.ftime[[1]])) & !("speedglm" %in%
                                                            class(glm.ftime[[1]]))) {
                    Qj_mod <- fast_glm(reg_form = stats::as.formula(Qj_form),
                                       data = dataList[[1]][NlessthanJ == 0, ],
                                       family = eval(glm_family))
                    if (unique(class(Qj_mod) %in% c("glm", "lm"))) {
                        Qj_mod <- cleanglm(Qj_mod)
                    }
                }
                else {
                    Qj_mod <- glm.ftime[[paste0("J", j)]]
                }
                ftimeMod[[paste0("J", j)]] <- NULL
                if (returnModels)
                    ftimeMod[[paste0("J", j)]] <- Qj_mod
                dataList <- lapply(dataList, function(x, j) {
                    suppressWarnings(x[[paste0("Q", j, "PseudoHaz")]] <- predict(Qj_mod,
                                                                                 newdata = x,
                                                                                 type = "response"))
                    if (j != min(J)) {
                        x[[paste0("hazLessThan", j)]] <-
                            rowSums(cbind(rep(0, nrow(x)), x[, paste0("Q", J[J < j], "Haz")]))
                        x[[paste0("Q", j, "Haz")]] <-
                            x[[paste0("Q", j, "PseudoHaz")]] * (1 - x[[paste0("hazLessThan", j)]])
                    }
                    else {
                        x[[paste0("hazLessThan", j)]] <- 0
                        x[[paste0("Q", j, "Haz")]] <- x[[paste0("Q", j, "PseudoHaz")]]
                    }
                    x
                }, j = j)
            }
        }
        else {
            for (j in J) {
                Qj_form <- sprintf("%s ~ %s", paste("N", j, sep = ""), glm.ftime)
                X <- stats::model.matrix(stats::as.formula(Qj_form), data = dataList[[1]])
                NlessthanJ <- rep(0, nrow(dataList[[1]]))
                for (i in J[J < j]) {
                    NlessthanJ <- NlessthanJ + dataList[[1]][[paste0("N", i)]]
                }
                dataList <- lapply(dataList, function(x, j) {
                    if (j != min(J)) {
                        x[[paste0("hazLessThan", j)]] <-
                            rowSums(cbind(rep(0, nrow(x)), x[, paste0("Q", J[J < j], "Haz")]))
                    }
                    else {
                        x[[paste0("hazLessThan", j)]] <- 0
                    }
                    x
                }, j = j)
                Ytilde <- (dataList[[1]][[paste0("N", j)]] - dataList[[1]][[paste0("l", j)]]) /
                    (pmin(dataList[[1]][[paste0("u", j)]],
                          1 - dataList[[1]][[paste0("hazLessThan", j)]]) -
                         dataList[[1]][[paste0("l", j)]])
                if (class("glm.ftime") != "list") {
                    Qj_mod <- stats::optim(par = rep(0, ncol(X)),
                                           fn = LogLikelihood, Y = Ytilde, X = X, method = "BFGS",
                                           gr = grad, control = list(reltol = 1e-07, maxit = 50000))
                }
                else {
                    Qj_mod <- glm.ftime[[paste0("J", j)]]
                }
                if (Qj_mod$convergence != 0) {
                    stop("convergence failure")
                }
                else {
                    beta <- Qj_mod$par
                    eval(parse(text = paste0("ftimeMod$J", j, " <- Qj_mod")))
                    dataList <- lapply(dataList, function(x, j) {
                        newX <- stats::model.matrix(stats::as.formula(Qj_form), data = x)
                        x[[paste0("Q", j, "PseudoHaz")]] <- plogis(newX %*% beta)
                        x[[paste0("Q", j, "Haz")]] <- (pmin(x[[paste0("u", j)]],
                                                            1 - x[[paste0("hazLessThan", j)]]) -
                                                           x[[paste0("l", j)]]) *
                            x[[paste0("Q", j, "PseudoHaz")]] + x[[paste0("l", j)]]
                        x
                    }, j = j)
                }
            }
        }
    }
    else if (is.null(glm.ftime)) {
        for (j in J) {
            NlessthanJ <- rep(0, nrow(dataList[[1]]))
            for (i in J[J < j]) {
                NlessthanJ <- NlessthanJ + dataList[[1]][[paste0("N", i)]]
            }
            if (max(c("data.frame", "SuperLearner") %in% class(SL.ftime[[1]])) == 0) {
                if(inherits(cl, "cluster")) {
                    Qj_mod <- SuperLearner::snowSuperLearner(
                        cluster = cl,
                        Y = dataList[[1]][[paste0("N", j)]][NlessthanJ == 0],
                        X = dataList[[1]][NlessthanJ == 0, c("t", "trt", names(adjustVars))],
                        id = dataList[[1]]$id[NlessthanJ ==0],
                        family = stats::binomial(), SL.library = SL.ftime,
                        cvControl = cv.Control, verbose = verbose)
                } else {
                    Qj_mod <- SuperLearner::SuperLearner(
                        Y = dataList[[1]][[paste0("N", j)]][NlessthanJ == 0],
                        X = dataList[[1]][NlessthanJ == 0, c("t", "trt", names(adjustVars))],
                        id = dataList[[1]]$id[NlessthanJ == 0],
                        family = stats::binomial(), SL.library = SL.ftime,
                        cvControl = cv.Control, verbose = verbose)
                }
            }
            else {
                Qj_mod <- SL.ftime[[paste0("J", j)]]
            }
            ftimeMod[[paste0("J", j)]] <- Qj_mod
            dataList <- lapply(dataList, function(x, j) {
                if (max("data.frame" %in% class(SL.ftime[[1]])) > 0) {
                    suppressWarnings(x <- left_join(x, Qj_mod,
                                                    by = c("id", "t", "trt")))
                }
                else {
                    suppressWarnings(x[[paste0("Q", j, "PseudoHaz")]] <-
                                         predict(Qj_mod, onlySL = TRUE,
                                                 newdata = x[, c("t", "trt",
                                                                 names(adjustVars))]
                                         )[[1]])
                }
                if (j != min(J)) {
                    x[[paste0("hazLessThan", j)]] <- rowSums(cbind(rep(0, nrow(x)),
                                                                   x[, paste0("Q", J[J < j], "Haz")]))
                    x[[paste0("Q", j, "Haz")]] <- x[[paste0("Q", j, "PseudoHaz")]] *
                        (1 - x[[paste0("hazLessThan", j)]])
                }
                else {
                    x[[paste0("Q", j, "Haz")]] <- x[[paste0("Q", j, "PseudoHaz")]]
                    x[[paste0("hazLessThan", j)]] <- 0
                }
                x
            }, j = j)
        }
    }
    out <- list()
    out$dataList <- dataList
    out$ftimeMod <- ftimeMod
    return(out)
}
environment(my_estimate_hazards) <- asNamespace('survtmle')
