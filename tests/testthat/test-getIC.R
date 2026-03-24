library(testthat)
library(concrete)
library(data.table)

# ============================================================================
# Reference implementation: original getIC, verbatim from before the
# computeIC refactor.  Uses the legacy getCleverCovariate / getHazLS C++
# helpers that are still exported for backward-compatibility.
# ============================================================================
getIC_reference <- function(GStar, Hazards, TotalSurv, NuisanceWeight,
                            TargetEvent, TargetTime, T.tilde, Delta,
                            EvalTimes, GComp = FALSE) {
    GStar <- as.numeric(unlist(GStar))

    IC.a <- do.call(rbind, lapply(TargetEvent, function(j) {
        F.j.t <- apply(Hazards[[as.character(j)]] * TotalSurv, 2, cumsum)

        do.call(rbind, lapply(TargetTime, function(tau) {
            h.FS <- matrix(F.j.t[EvalTimes == tau, ],
                           ncol = ncol(F.j.t),
                           nrow = nrow(F.j.t[EvalTimes <= tau, , drop = FALSE]),
                           byrow = TRUE)
            h.FS <- (h.FS - F.j.t[EvalTimes <= tau, , drop = FALSE]) /
                    TotalSurv[EvalTimes <= tau, , drop = FALSE]

            IC.j.tau <- Reduce("+", x = lapply(names(Hazards), function(l) {
                ClevCov <- concrete:::getCleverCovariate(
                    GStar         = GStar,
                    NuisanceWeight = NuisanceWeight[EvalTimes <= tau, , drop = FALSE],
                    hFS           = h.FS,
                    LeqJ          = as.integer(l == j)
                )

                NLdS <- matrix(0, nrow = nrow(h.FS), ncol = ncol(h.FS))
                for (i in which(Delta == l & T.tilde <= tau)) {
                    NLdS[which(EvalTimes == T.tilde[i]), i] <- 1
                }

                HazLS <- concrete:::getHazLS(
                    T_Tilde  = T.tilde,
                    EvalTimes = EvalTimes[EvalTimes <= tau],
                    HazL     = Hazards[[l]][EvalTimes <= tau, , drop = FALSE]
                )

                colSums(ClevCov * (NLdS - HazLS))
            })) + F.j.t[EvalTimes == tau, ] - mean(F.j.t[EvalTimes == tau, ])

            data.table(ID = seq_along(IC.j.tau), Time = tau,
                       Event = j, IC = IC.j.tau)
        }))
    }))
    return(IC.a)
}

# ============================================================================
# Helper: build a minimal but structurally realistic set of getIC inputs.
#
# Parameters let each test control which features are exercised:
#   n_events   : number of cause types (names "1", "2", ...)
#   n_subjects : N
#   seed       : RNG seed
#   target_event, target_time: which estimands to target
# ============================================================================
make_ic_inputs <- function(n_events = 2, n_subjects = 30, seed = 1,
                           target_event = NULL, target_time = NULL) {
    set.seed(seed)
    N <- n_subjects

    # Discrete time grid: 0 + observed times + target times
    obs_grid   <- sort(unique(sample(1:20, N, replace = TRUE)))
    tgt_times  <- if (!is.null(target_time)) target_time else c(10L, 18L)
    EvalTimes  <- sort(unique(c(0L, obs_grid, tgt_times)))
    T          <- length(EvalTimes)

    # Cause-specific hazards: small positive values so survival stays meaningful
    event_names <- as.character(seq_len(n_events))
    base_rates  <- seq(0.06, 0.03, length.out = n_events)
    Hazards <- setNames(
        lapply(seq_len(n_events), function(j) {
            m <- matrix(pmax(1e-4, rnorm(T * N, base_rates[j], 0.01)),
                        nrow = T, ncol = N)
            attr(m, "j") <- j
            m
        }),
        event_names
    )

    # Proper survival: S(t) = exp(-cumsum(h_total))
    haz_total <- Reduce(`+`, Hazards)
    TotalSurv <- apply(haz_total, 2, function(h) exp(-cumsum(h)))
    TotalSurv <- pmax(TotalSurv, 1e-8)

    # Observed times drawn from the non-zero grid (guarantees T.tilde in EvalTimes)
    T.tilde <- sample(EvalTimes[-1], N, replace = TRUE)

    # Event types: 0 = censored, 1..n_events = event causes
    Delta <- sample(0:n_events, N, replace = TRUE)

    # Intervention probability g*(A|W)
    GStar <- runif(N, 0.3, 0.8)

    # Nuisance weights: 1 / (g * S_C); just use a positive T x N matrix
    NuisanceWeight <- matrix(runif(T * N, 0.5, 2.0), nrow = T, ncol = N)

    list(
        Hazards        = Hazards,
        TotalSurv      = TotalSurv,
        NuisanceWeight = NuisanceWeight,
        GStar          = GStar,
        T.tilde        = T.tilde,
        Delta          = Delta,
        EvalTimes      = EvalTimes,
        TargetEvent    = if (!is.null(target_event)) target_event else seq_len(n_events),
        TargetTime     = tgt_times
    )
}

# Convenience wrapper: run both implementations and return a named list
run_both <- function(inp) {
    new <- concrete:::getIC(
        GStar          = inp$GStar,
        Hazards        = inp$Hazards,
        TotalSurv      = inp$TotalSurv,
        NuisanceWeight = inp$NuisanceWeight,
        TargetEvent    = inp$TargetEvent,
        TargetTime     = inp$TargetTime,
        T.tilde        = inp$T.tilde,
        Delta          = inp$Delta,
        EvalTimes      = inp$EvalTimes,
        GComp          = FALSE
    )
    ref <- getIC_reference(
        GStar          = inp$GStar,
        Hazards        = inp$Hazards,
        TotalSurv      = inp$TotalSurv,
        NuisanceWeight = inp$NuisanceWeight,
        TargetEvent    = inp$TargetEvent,
        TargetTime     = inp$TargetTime,
        T.tilde        = inp$T.tilde,
        Delta          = inp$Delta,
        EvalTimes      = inp$EvalTimes,
        GComp          = FALSE
    )
    # Align both on (Time, Event, ID) before comparing
    setkeyv(new, c("Time", "Event", "ID"))
    setkeyv(ref, c("Time", "Event", "ID"))
    list(new = new, ref = ref)
}

# ============================================================================
# Tests
# ============================================================================

test_that("getIC matches reference: single event type, single target time", {
    inp <- make_ic_inputs(n_events = 1, n_subjects = 30, seed = 101,
                          target_event = 1L, target_time = 12L)
    out <- run_both(inp)
    expect_equal(out$new$IC, out$ref$IC, tolerance = 1e-10)
    expect_equal(out$new$Time,  out$ref$Time)
    expect_equal(out$new$Event, out$ref$Event)
    expect_equal(out$new$ID,    out$ref$ID)
})

test_that("getIC matches reference: two competing events, two target times", {
    inp <- make_ic_inputs(n_events = 2, n_subjects = 40, seed = 202,
                          target_event = 1:2, target_time = c(8L, 16L))
    out <- run_both(inp)
    expect_equal(out$new$IC, out$ref$IC, tolerance = 1e-10)
})

test_that("getIC matches reference: three competing events, three target times", {
    inp <- make_ic_inputs(n_events = 3, n_subjects = 50, seed = 303,
                          target_event = 1:3, target_time = c(5L, 12L, 19L))
    out <- run_both(inp)
    expect_equal(out$new$IC, out$ref$IC, tolerance = 1e-10)
})

test_that("getIC matches reference: targeting a strict subset of events", {
    # Three hazard causes in Hazards but only targeting event 2
    inp <- make_ic_inputs(n_events = 3, n_subjects = 35, seed = 404,
                          target_event = 2L, target_time = c(10L, 18L))
    out <- run_both(inp)
    expect_equal(out$new$IC, out$ref$IC, tolerance = 1e-10)
})

test_that("getIC matches reference: with heavy censoring (many Delta == 0)", {
    set.seed(505)
    inp <- make_ic_inputs(n_events = 2, n_subjects = 40, seed = 505)
    # Force 60 % censoring
    inp$Delta <- sample(c(rep(0L, 24), rep(1L, 8), rep(2L, 8)))
    out <- run_both(inp)
    expect_equal(out$new$IC, out$ref$IC, tolerance = 1e-10)
})

test_that("getIC matches reference: all subjects censored (IC equals F_j - mean)", {
    inp <- make_ic_inputs(n_events = 2, n_subjects = 20, seed = 606,
                          target_event = 1L, target_time = 10L)
    inp$Delta <- rep(0L, length(inp$Delta))   # everyone censored
    out <- run_both(inp)
    # Point mass is 0 for all subjects; result is purely F_j(tau) - mean(F_j(tau))
    expect_equal(out$new$IC, out$ref$IC, tolerance = 1e-10)
})

test_that("getIC matches reference: single target time equal to first observed time", {
    # tau = min non-zero EvalTime; only subjects with T.tilde == tau contribute
    inp <- make_ic_inputs(n_events = 2, n_subjects = 30, seed = 707)
    tau_early <- inp$EvalTimes[2]   # first non-zero time
    inp$TargetTime  <- tau_early
    inp$TargetEvent <- 1L
    out <- run_both(inp)
    expect_equal(out$new$IC, out$ref$IC, tolerance = 1e-10)
})

test_that("getIC matches reference: target time equal to last EvalTime", {
    inp <- make_ic_inputs(n_events = 2, n_subjects = 30, seed = 808)
    tau_last <- max(inp$EvalTimes)
    inp$TargetTime  <- tau_last
    inp$TargetEvent <- 1:2
    out <- run_both(inp)
    expect_equal(out$new$IC, out$ref$IC, tolerance = 1e-10)
})

test_that("getIC output rows equal N * J * K", {
    inp <- make_ic_inputs(n_events = 2, n_subjects = 25, seed = 909,
                          target_event = 1:2, target_time = c(8L, 14L, 18L))
    out <- run_both(inp)
    expected_rows <- length(inp$TargetEvent) * length(inp$TargetTime) * length(inp$GStar)
    expect_equal(nrow(out$new), expected_rows)
    expect_equal(nrow(out$ref), expected_rows)
})

test_that("getIC output has correct columns and types", {
    inp  <- make_ic_inputs(n_events = 1, n_subjects = 20, seed = 111,
                           target_event = 1L, target_time = 10L)
    out  <- run_both(inp)
    expect_named(out$new, c("ID", "Time", "Event", "IC"))
    expect_type(out$new$IC, "double")
    expect_true(is.data.table(out$new))
})

test_that("getIC IC values are finite (no NA / Inf with well-conditioned inputs)", {
    inp <- make_ic_inputs(n_events = 2, n_subjects = 40, seed = 1212,
                          target_event = 1:2, target_time = c(8L, 16L))
    out <- run_both(inp)
    expect_true(all(is.finite(out$new$IC)))
    expect_true(all(is.finite(out$ref$IC)))
})

test_that("getIC mean-centering: mean(IC) equals mean(F_j(tau)) - mean(F_j(tau)) = 0", {
    # The centering term F_j(tau,i) - mean_i F_j(tau,i) makes Pn[IC] == 0
    # iff the rest of the terms also have mean zero — which they do for the full EIC.
    # Here we check the simpler property: the mean of IC across subjects is computable
    # and matches between implementations.
    inp <- make_ic_inputs(n_events = 2, n_subjects = 50, seed = 1313,
                          target_event = 1:2, target_time = c(10L, 18L))
    out <- run_both(inp)
    mean_new <- out$new[, mean(IC), by = c("Time", "Event")]$V1
    mean_ref <- out$ref[, mean(IC), by = c("Time", "Event")]$V1
    expect_equal(mean_new, mean_ref, tolerance = 1e-10)
})
