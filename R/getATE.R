getATE <- function(ConcreteOutput) {
    TmleIC <- ConcreteOutput$ic

    TmleEst <- ConcreteOutput$estimate$tmle
    TmleEst <- melt(TmleEst, id.vars = c("A", "time"),
                variable.name = "J", value.name = "Psi")
    TmleEst <- dcast(TmleEst, time + J ~ A, value.var = "Psi")
    TmleEst[, ATE := get('1') - get('0')]

    TmleSE <- ConcreteOutput$se
    TmleSE <- data.table("name" = names(TmleSE), "val" = TmleSE)
    TmleSE <- TmleSE[, c("J", "A", "time") := tstrsplit(name, "\\.(t|a)")]
    TmleSE[, time := as.numeric(time)]
    TmleSE[, J := gsub("j", "", J)]
    TmleSE <- dcast(TmleSE, J + time ~ A, value.var = "val")
    TmleSE[, "se" := sqrt(get('0')^2 + get('1')^2)]

    ATE_out <- merge(TmleEst[, -c("0", "1")], TmleSE[, -c("0", "1")],
                     by = c("time", "J"))

    return(ATE_out)
}
