#' Intervention
#' a function that takes the observed interventions and transforms them into the desired interventions
#' @param a the numeric vector of observed interventions, or a data.table? for longitudinal? or long format?
#'
#' @return a.intervened, a numeric vector equal in size to a, representing the desired intervention
#'
#' @export
#'
#' @examples
Intervention <- list("A = a" = function(A, CovDataTable) {
    regime <- #getDesiredTrtRegime(A, CovDataTable)
    attr(regime, "g.star") <- #getPropScoreOfDesiredTrtRegime()
        # maybe like predict(TrtFit, A = regime, CovDataTable) or
        # as.numeric(A == regime)
return(regime)
})
