#### must define in foreach environment?

MOSS_mods <- new.env()


# haz2surv ----------------------------------------------------------------


survival_curve$set(which = "public", name = "haz2surv", function() {
    self$survival <- do.call(rbind, lapply(1:nrow(self$hazard), function(i) {
        hazard_here <- head(c(0, self$hazard[i, ]), -1)
        return(cumprod(1 - hazard_here))
    }))
    return(self)
}, overwrite = TRUE)



# surv2pdf ----------------------------------------------------------------


survival_curve$set(which = "public", name = "surv2pdf", function() {
    self$pdf <- do.call(rbind, lapply(1:self$n(), function(i) {
        c(-diff(self$survival[i, ]), 0)
    }))
    return(self)
}, overwrite = TRUE)



# surv2haz ----------------------------------------------------------------


survival_curve$set(which = "public", name = "surv2haz", function() {
    self$surv2pdf()
    self$hazard <- self$pdf/self$survival
    return(self)
}, overwrite = TRUE)



# my_clever_covariate -----------------------------------------------------


eic$set(which = "public", name = "my_clever_covariate",
        function (t, k_grid)
{
    g <- self$g1W
    if (self$A_intervene == 0)
        g <- 1 - g

    h_matrix <- do.call(cbind, lapply(k_grid, function(k) {
        h <- -as.numeric(self$A == self$A_intervene) * as.numeric(k <= t) /
            (g * self$density_censor$survival[, k]) *
            self$density_failure$survival[, t] / self$density_failure$survival[, k]
    }))

    return(as.vector(t(h_matrix)))
}, overwrite = TRUE)



# my_MOSS_hazard ----------------------------------------------------------


my_MOSS_hazard <- R6Class(classname = "my_MOSS_Hazard", inherit = MOSS_hazard,
                          public = list(t_max = NULL, targets = NULL,
                                        initialize = function(...)
{
    args <- list(...)
    pass <- args[setdiff(names(args), c("t_max", "targets"))]
    do.call(super$initialize, pass)
    self$t_max <- args[['t_max']]
    self$targets <- args[['targets']]
}))



# my_iterate_onestep ------------------------------------------------------


my_MOSS_hazard$set(which = "public", name = "my_iterate_onestep",
                function (method = "l2", epsilon = 1, max_num_iteration = 10,
                          tmle_tolerance = NULL, verbose = FALSE)
{
    self$epsilon <- epsilon
    self$max_num_interation <- max_num_iteration
    if (identical(self$A_intervene, c(1, 0))) {
        if (!is.list(self$density_failure) | !is.list(self$density_censor) |
            length(self$density_failure) != 2 | length(self$density_censor) != 2)
            stop("Treatment and control failure and censoring densities must
                 be supplied to target both curves")
    }
    self$tmle_tolerance <- ifelse(is.null(tmle_tolerance),
                                  1/length(self$A),
                                  tmle_tolerance)
    if (is.null(self$t_max))
        self$t_max <- max(self$T_tilde)
    if (is.null(self$k_grid))
        self$k_grid <- 1:self$t_max
    if (is.null(self$targets))
        self$targets <- self$k_grid

    if (identical(self$A_intervene, c(1, 0))) {
        psi_n <- lapply(1:2, function(i) {colMeans(self$density_failure[[i]]$survival)})
        stats_eic <- lapply(1:2, function(i) {
            self$compute_stats_eic(psi_n = psi_n[[i]],
                                   a_intervene = self$A_intervene[[i]],
                                   density_failure = self$density_failure[[i]],
                                   density_censor = self$density_censor[[i]])
        })
        stats_eic <- do.call(rbind, stats_eic)
        psi_n <- do.call(c, psi_n)
    } else {
        psi_n <- colMeans(self$density_failure$survival)
        stats_eic <- self$compute_stats_eic(psi_n = psi_n,
                                            a_intervene = self$A_intervene,
                                            density_failure = self$density_failure,
                                            density_censor = self$density_censor)
    }
    num_iteration <- 0
    if (verbose) {
        print(cbind("iteration" = num_iteration, stats_eic))
    }
    mean_eic <- as.vector(stats_eic[, "mean"])
    mean_eic_inner_prod_prev <- abs(sqrt(sum(mean_eic^2)))
    mean_eic_inner_prod_current <- mean_eic_inner_prod_prev
    mean_eic_inner_prod_best <- sqrt(sum(mean_eic^2))
    if (is.list(self$density_failure)) {
        self$q_best <- lapply(self$density_failure, function(dens_fail) {
            dens_fail$clone(deep=TRUE)
        })
    } else {
        self$q_best <- self$density_failure$clone(deep = TRUE)
    }
    to_iterate <- TRUE
    if (is.infinite(mean_eic_inner_prod_current) | is.na(mean_eic_inner_prod_current)) {
        to_iterate <- FALSE
    }
    while (mean_eic_inner_prod_current >= self$tmle_tolerance * sqrt(length(self$targets)*2)
           & to_iterate)
    {
        self$my_fit_epsilon(method = method, clipping = self$epsilon)

        if (identical(self$A_intervene, c(1, 0))) {
            psi_n <- lapply(1:2, function(i) {colMeans(self$density_failure[[i]]$survival)})
            stats_eic <- lapply(1:2, function(i) {
                self$compute_stats_eic(psi_n = psi_n[[i]],
                                       a_intervene = self$A_intervene[[i]],
                                       density_failure = self$density_failure[[i]],
                                       density_censor = self$density_censor[[i]])
            })
            stats_eic <- do.call(rbind, stats_eic)
            psi_n <- do.call(c, psi_n)
        } else {
            psi_n <- colMeans(self$density_failure$survival)
            stats_eic <- self$compute_stats_eic(psi_n = psi_n,
                                                a_intervene = self$A_intervene,
                                                density_failure = self$density_failure,
                                                density_censor = self$density_censor)
        }
        num_iteration <- num_iteration + 1
        if (verbose) {
            print(cbind("iteration" = num_iteration, stats_eic))
        }
        mean_eic <- as.vector(stats_eic[, "mean"])
        mean_eic_inner_prod_prev <- mean_eic_inner_prod_current
        mean_eic_inner_prod_current <- abs(sqrt(sum(mean_eic^2)))
        if (is.infinite(mean_eic_inner_prod_current) | is.na(mean_eic_inner_prod_current)) {
            warning("stopping criteria diverged. Reporting best result so far.")
            (break)()
        }
        if (mean_eic_inner_prod_current < mean_eic_inner_prod_best) {
            if (is.list(self$density_failure)) {
                self$q_best <- lapply(self$density_failure, function(dens_fail) {
                    dens_fail$clone(deep=TRUE)
                })
            } else {
                self$q_best <- self$density_failure$clone(deep = TRUE)
            }
            mean_eic_inner_prod_best <- mean_eic_inner_prod_current
        }
        if (num_iteration == self$max_num_interation) {
            (break)()
            warning("Max number of iteration reached, stop TMLE")
        }
        if (all(abs(mean_eic) < as.vector(stats_eic[, "cutoff"]))) {
            (break)()
            warning("mean EIC less than variance based cutoff, stop TMLE updates")
        }
    }
    self$density_failure <- self$q_best
    if (identical(self$A_intervene, c(1, 0))) {
        psi_n <- lapply(1:2, function(i) {
            colMeans(self$density_failure[[i]]$survival)})
    } else {
        psi_n <- colMeans(self$density_failure$survival)
    }
    if (verbose) {
        message(paste("Pn(EIC)=", formatC(mean_eic_inner_prod_best,
                                          format = "e", digits = 2)))
    }
    return(psi_n)
}, overwrite = TRUE)



# compute_stats_eic -------------------------------------------------------


my_MOSS_hazard$set(which = "public", name = "compute_stats_eic",
                function (psi_n, a_intervene, density_failure, density_censor)
{
    eic_fit <- eic$new(A = self$A, T_tilde = self$T_tilde, Delta = self$Delta,
                       density_failure = density_failure,density_censor = density_censor,
                       g1W = self$g1W, psi = psi_n, A_intervene = a_intervene
                       )$all_t(k_grid = self$targets)
    mean_eic <- colMeans(eic_fit)
    var_eic <- as.vector(diag(var(eic_fit)))
    cutoff_eic <- sqrt(var_eic) / (sqrt(length(self$A)) * log(length(self$A)))
    eic_stats <- matrix(data = c(mean_eic, var_eic, cutoff_eic), ncol = 3, byrow = F)
    colnames(eic_stats) <- c("mean", "var", "cutoff")
    eic_stats <- cbind("time" = self$targets, eic_stats)
    return(eic_stats)
}, overwrite = TRUE)



# my_fit_epsilon ----------------------------------------------------------


my_MOSS_hazard$set(which = "public", name = "my_fit_epsilon",
                function (method = "l2", clipping = Inf)
{
    if (is.null(self$t_max))
        self$t_max <- max(self$T_tilde)
    dNt <- self$my_create_dNt()
    if (identical(self$A_intervene, c(1, 0))) {
        h_matrix <- do.call(cbind, lapply(1:2, function(i) {
            self$make_long_h_matrix(A_intervene = self$A_intervene[i],
                                    density_failure = self$density_failure[[i]],
                                    density_censor = self$density_censor[[i]])
        }))
    } else {
        h_matrix <- self$make_long_h_matrix(A_intervene = self$A_intervene,
                                            density_failure = self$density_failure,
                                            density_censor = self$density_censor)
    }
    if (method == "PnD") {
        if (identical(self$A_intervene, c(1, 0))) {
            psi_n <- lapply(1:2, function(i) {colMeans(self$density_failure[[i]]$survival)})
            stats_eic <- lapply(1:2, function(i) {
                self$compute_stats_eic(psi_n = psi_n[[i]],
                                       a_intervene = self$A_intervene[[i]],
                                       density_failure = self$density_failure[[i]],
                                       density_censor = self$density_censor[[i]])
            })
            stats_eic <- do.call(rbind, stats_eic)
            psi_n <- do.call(c, psi_n)
        } else {
            psi_n <- colMeans(self$density_failure$survival)
            stats_eic <- self$compute_stats_eic(psi_n = psi_n,
                                                a_intervene = self$A_intervene,
                                                density_failure = self$density_failure,
                                                density_censor = self$density_censor)
        }
        epsilon_n <- as_tibble(stats_eic)$mean
        l2_norm <- sqrt(sum(epsilon_n^2))
        if (l2_norm >= clipping) {
            epsilon_n <- epsilon_n/l2_norm * clipping
        }
    }
    if (method == "glm") {
        if (identical(self$A_intervene, c(1, 0))) {
            lambda_AW <- self$density_failure[[1]]$hazard
            lambda_AW[self$A == 0, ] <- self$density_failure[[2]]$hazard[self$A == 0, ]
            submodel_fit <- glm(dNt ~ . - 1, family = stats::binomial(),
                                data = as.data.frame(h_matrix),
                                offset = qlogis(as.vector(t(lambda_AW))))
        } else {
            submodel_fit <- glm(dNt ~ . - 1, family = stats::binomial(),
                                data = as.data.frame(h_matrix),
                                offset = qlogis(as.vector(t(lambda_AW))))
        }
        epsilon_n <- unname(submodel_fit$coefficients)
        l2_norm <- sqrt(sum(epsilon_n^2))
        if (l2_norm >= clipping) {
            epsilon_n <- epsilon_n/l2_norm * clipping
        }
        ll <- predict(submodel_fit, type = "response")
        ll[dNt == 0] <- 1 - ll[dNt == 0]
        print(c("loglik" = sum(log(ll)),
                "clipped epsilon" = epsilon_n))
    }
    if (method %in% c("l2", "l1")) {
        epsilon_n <- tryCatch({
            if (method == "l2") {
                alpha <- 0
                norm_func <- MOSS:::norm_l2
                lambda.min.ratio = 0.01
            }
            if (method == "l1") {
                alpha <- 1
                norm_func <- MOSS:::norm_l1
                lambda.min.ratio = 0.9
            }
            ind <- 1
            while (ind == 1) {
                if (identical(self$A_intervene, c(1, 0))) {
                    lambda_AW <- self$density_failure[[1]]$hazard
                    lambda_AW[self$A == 0, ] <-
                        self$density_failure[[2]]$hazard[self$A == 0, ]
                    tmle_offset <- qlogis(as.vector(t(lambda_AW)))
                    enet_fit <- glmnet::glmnet(x = h_matrix, y = dNt,
                        offset = tmle_offset, family = "binomial", alpha = alpha,
                        standardize = FALSE, intercept = FALSE,
                        lambda.min.ratio = lambda.min.ratio, nlambda = 200)
                } else {
                    tmle_offset <- qlogis(as.vector(t(self$density_failure$hazard)))
                    enet_fit <- glmnet::glmnet(x = h_matrix, y = dNt,
                        offset = tmle_offset, family = "binomial", alpha = alpha,
                        standardize = FALSE, intercept = FALSE,
                        lambda.min.ratio = lambda.min.ratio, nlambda = 200)
                }
                norms <- apply(enet_fit$beta, 2, norm_func)
                ind <- max(which(norms <= clipping))
                if (ind > 1) {
                    break
                }
                lambda.min.ratio <- (lambda.min.ratio + 1)/2
            }
            epsilon_n <- enet_fit$beta[, ind]
        }, error = function(e) {
            print("elastic net update failed, returning epsilon = 0")
            return(rep(0, ncol(h_matrix)))
        })
        ll <- predict(enet_fit, newx = h_matrix, newoffset = tmle_offset,
                      s = enet_fit$lambda[ind], type = "response")
        ll[dNt == 0] <- 1 - ll[dNt == 0]
        print(c("loglik" = sum(log(ll)),
                "clipped epsilon" = epsilon_n))

    }
    if (identical(self$A_intervene, c(1, 0))) {
        offset <- h_matrix %*% epsilon_n
        for (i in 1:2) {
            hazard_new <- plogis(qlogis(as.vector(t(self$density_failure[[i]]$hazard))) +
                                     as.vector(offset))
            self$density_failure[[i]]$hazard <- matrix(hazard_new, nrow = length(self$A),
                                                       byrow = TRUE)
            self$density_failure[[i]]$haz2surv()
        }
    } else {
        hazard_new <- plogis(qlogis(as.vector(t(self$density_failure$hazard))) +
                                as.vector(h_matrix %*% epsilon_n))
        # hazard_new <- matrix(hazard_new, nrow = length(self$A), byrow = TRUE)
        # survival_new <- survival_curve$new(t = self$k_grid, hazard = hazard_new)$haz2surv()
        self$density_failure$hazard <- matrix(hazard_new, nrow = length(self$A), byrow = TRUE)
        self$density_failure$haz2surv()
    }
    # return(survival_new)
}, overwrite = TRUE)



# make_long_h_matrix ------------------------------------------------------


my_MOSS_hazard$set(which = "public", name = "make_long_h_matrix",
                function (A_intervene, density_failure, density_censor)
{
    if (is.null(density_failure$survival)) density_failure$haz2surv()
    psi_n <- colMeans(density_failure$survival)
    eic_fit <- eic$new(A = self$A, T_tilde = self$T_tilde, Delta = self$Delta,
                       density_failure = density_failure, density_censor = density_censor,
                       g1W = self$g1W, psi = psi_n, A_intervene = A_intervene)
    if (length(self$targets) > 1) {
    h_matrix <- do.call(cbind, lapply(X = self$targets, function(t)
        eic_fit$my_clever_covariate(t = t, k_grid = self$k_grid)))
    } else {
        h_matrix <- eic_fit$my_clever_covariate(t = self$targets, k_grid = self$k_grid)
    }
    return(h_matrix)
}, overwrite = TRUE)



# my_create_dNt -----------------------------------------------------------


my_MOSS_hazard$set(which = "public", name = "my_create_dNt", function () {
    if (is.null(self$t_max))
        self$t_max <- max(self$T_tilde)
    dNt <- matrix(0, nrow = self$t_max, ncol = length(self$A))
    for (i in 1:length(self$A)) {
        if (self$Delta[i] == 1 & self$T_tilde[i] <= self$t_max) {
            dNt[self$T_tilde[i], i] <- 1
        }
    }
    return(as.vector(dNt))
}, overwrite = TRUE)



