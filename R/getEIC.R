#' get EICs
#'
#' @param Estimates list
#' @param Data data.table
#' @param Regime list
#' @param TargetEvent numeric vector
#' @param TargetTime numeric vector
#' @param MinNuisance numeric
#' @param GComp boolean
#'
#' @useDynLib concrete
#' @importFrom Rcpp evalCpp
#' @exportPattern "ˆ[[:alpha:]]+"
#' @importFrom stats var

getEIC <- function(Estimates, Data, Regime, TargetEvent, TargetTime, MinNuisance, GComp = FALSE) {
    EvalTimes <- attr(Estimates, "Times")
    # Censored <- 0 %in% Data[[attr(Data, "EventType")]]
    T.tilde <- Data[[attr(Data, "EventTime")]]
    Delta <- Data[[attr(Data, "EventType")]]
    
    for (a in seq_along(Estimates)) {
        NuisanceWeight <- Estimates[[a]][["NuisanceWeight"]]
        GStar <- attr(Estimates[[a]][["PropScore"]], "g.star.obs")
        Hazards <- Estimates[[a]][["Hazards"]]
        TotalSurv <- Estimates[[a]][["EvntFreeSurv"]]
        
        IC.a <- getIC(GStar = GStar, Hazards = Hazards, TotalSurv = TotalSurv,
                      NuisanceWeight = NuisanceWeight, TargetEvent = TargetEvent,
                      TargetTime = TargetTime, T.tilde = T.tilde,
                      Delta = Delta, EvalTimes = EvalTimes, GComp = GComp)
        
        if (GComp)
            Estimates[[a]][["GCompEst"]] <- getGComp(EvalTimes, Hazards, TotalSurv, TargetTime)
        
        Estimates[[a]][["SummEIC"]] <- summarizeIC(IC.a)
        Estimates[[a]][["IC"]] <- IC.a
    }
    return(Estimates)
}

getIC <- function(GStar, Hazards, TotalSurv, NuisanceWeight, TargetEvent, TargetTime, 
                  T.tilde, Delta, EvalTimes, GComp) {
    Target <- expand.grid("Time" = TargetTime, "Event" = TargetEvent)
    UniqueEvents <- setdiff(sort(unique(Delta)), 0)
    GStar <- as.numeric(unlist(GStar))
    
    IC.a <- do.call(rbind, lapply(TargetEvent, function(j) {
        # Cumulative Incidence Function
        F.j.t <- apply(Hazards[[as.character(j)]] * TotalSurv, 2, cumsum)
        do.call(rbind, lapply(TargetTime, function(tau) {
            # The event-related (F(t) and S(t)) contributions to the clever covariate (h)  
            h.FS <- matrix(F.j.t[EvalTimes == tau, ], 
                           ncol = ncol(F.j.t), 
                           nrow = nrow(F.j.t[EvalTimes <= tau, ]), 
                           byrow = TRUE)
            h.FS <- (h.FS - F.j.t[EvalTimes <= tau, ]) / TotalSurv[EvalTimes <= tau, ]
            
            IC.j.tau <- Reduce("+", x = lapply(names(Hazards), function(l) {
                ClevCov <- getCleverCovariate(GStar = GStar, 
                                              NuisanceWeight = NuisanceWeight[EvalTimes <= tau, ], 
                                              hFS = h.FS, 
                                              LeqJ = as.integer(l == j))
                
                NLdS <- matrix(data = 0, nrow = nrow(h.FS), ncol = ncol(h.FS))
                for (i in which(Delta == l & T.tilde <= tau)) {
                    NLdS[which(EvalTimes == T.tilde[i]), i] <- 1
                }
                
                HazLS <- getHazLS(T_Tilde = T.tilde, 
                                  EvalTimes = EvalTimes[EvalTimes <= tau], 
                                  HazL = Hazards[[l]][EvalTimes <= tau, ])
                
                return(colSums(ClevCov * (NLdS - HazLS)))
            })) + F.j.t[EvalTimes == tau, ] - mean(F.j.t[EvalTimes == tau, ])
            if (anyNA(IC.j.tau))
                stop("IC overflow: either increase MinNuisance or specify a target estimand ", 
                     " (Target Event, Target Time, & Intervention) with more support in the data.")
            return(data.table("ID" = seq_along(IC.j.tau), "Time" = tau, "Event" = j, "IC" = IC.j.tau))
        }))
    }))
    return(IC.a)
}

getGComp <- function(EvalTimes, Hazards, TotalSurv, TargetTime) {
    F.j.tau <- Event <- Time <- NULL
    Risks <- do.call(rbind, lapply(Hazards, function(haz.j) {
        Risk.a <- sapply(1:ncol(haz.j), function(i) {
            cumsum(TotalSurv[, i] * haz.j[, i])
        })
        Risk.a <- cbind("Event" = as.numeric(attr(haz.j, "j")),
                        "Time" = EvalTimes[EvalTimes %in% TargetTime],
                        "F.j.tau" = rowMeans(subset(Risk.a, EvalTimes %in% TargetTime)))
        return(Risk.a)
    }))
    Risks <- as.data.table(Risks)
    Risks <- rbind(Risks, Risks[, list("Event" = -1, "F.j.tau" = 1 - sum(F.j.tau)), by = "Time"])
    return(Risks[, list("Event" = Event, "Time" = Time, "Risk" = F.j.tau)])
}

summarizeIC <- function(IC.a) {
    IC <- NULL
    IC.a <- rbind(IC.a,
                  IC.a[, list("Event" = -1, "IC" = -sum(IC)), by = c("ID", "Time")])
    summIC <- IC.a[, list("PnEIC" = mean(IC), 
                          "seEIC" = sqrt(stats::var(IC)),
                          "seEIC/(sqrt(n)log(n))" = 
                            sqrt(stats::var(IC)) / (3 * sqrt(.N) * log(.N))),
                   by = c("Time", "Event")]
    return(summIC)
}
