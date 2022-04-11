#' @title Performance metrics - total biomass
#
#' @description Calculate total biomass performance metrics
#'
#' @param mse HERE:
#'
#' @export
perform.tb <- function(mse){
    tmp <- mse$est[,,,"tb",]
    if(length(dim(tmp)) == 4) dimen <- 1:3 else dimen <- 1:2
    res <- apply(tmp, dimen, function(x) c(mean(x, na.rm=TRUE),
                                            quantile(x, probs = c(0.025,0.5,0.975), na.rm=TRUE),
                                            sd(x, na.rm=TRUE)))
    rownames(res) <- c("mean","lo","med","up","sd")
    return(res)
}



#' @title Performance metrics - spawning stock biomass
#
#' @description Calculate spawning stock biomass performance metrics
#'
#' @param mse HERE:
#'
#' @export
perform.ssb <- function(mse){
    tmp <- mse$est[,,,"ssb",]
    if(length(dim(tmp)) == 4) dimen <- 1:3 else dimen <- 1:2
    res <- apply(tmp, dimen, function(x) c(mean(x, na.rm=TRUE),
                                              quantile(x, probs = c(0.025,0.5,0.975), na.rm=TRUE),
                                              sd(x, na.rm=TRUE)))
    rownames(res) <- c("mean","lo","med","up","sd")
    return(res)
}



#' @title Performance metrics - yield
#
#' @description Calculate yield performance metrics
#'
#' @param mse HERE:
#'
#' @export
perform.yield <- function(mse){
    tmp <- mse$est[,,,"y",]
    if(length(dim(tmp)) == 4) dimen <- 1:3 else dimen <- 1:2
    res <- apply(tmp, dimen, function(x) c(mean(x, na.rm=TRUE),
                                              quantile(x, probs = c(0.025,0.5,0.975), na.rm=TRUE),
                                              sd(x, na.rm=TRUE)))
    rownames(res) <- c("mean","lo","med","up","sd")
    return(res)
}


#' @title Performance metrics - risk
#
#' @description Calculate risk performance metrics
#'
#' @param mse HERE:
#'
#' @export
perform.risk <- function(mse){
    tmp <- mse$est[,,,"ssb",]
    if(length(dim(tmp)) == 4) dimen <- 1:3 else dimen <- 1:2
    res <- apply(tmp, dimen, function(x) mean(x < mse$pars$refs[3], na.rm=TRUE)) ## make refs list with names => easier to check
    return(res)
}



#' @title Performance metrics - yield variability
#
#' @description Calculate yield variability performance metrics
#'
#' @param mse HERE:
#'
#' @export
perform.yielddiff <- function(mse){
    tmp <- mse$est[,,,"y",]
    if(length(dim(tmp)) == 4) dimen <- 2:4 else dimen <- 2:3
    tmp2 <- apply(tmp, dimen, diff)
    if(length(dim(tmp)) == 4) dimen <- 2:3 else dimen <- 2
    res <- apply(tmp2, dimen, function(x) c(mean(x, na.rm=TRUE),
                                             quantile(x, probs = c(0.025,0.5,0.975), na.rm=TRUE),
                                             sd(x, na.rm=TRUE)))
    rownames(res) <- c("mean","lo","med","up","sd")
    return(res)
}


#' @title Performance metrics
#
#' @description Calculate all performance metrics
#'
#' @param mse HERE:
#'
#' @export
perform <- function(mse){

    ## Checks
    if(!any(names(mse) == "est")) stop("Estimates missing. Did you run 'run.mse'?")

    ## TB
    tb <- perform.tb(mse)

    ## SSB
    ssb <- perform.ssb(mse)

    ## Yield
    yield <- perform.yield(mse)

    ## Risk
    risk <- perform.risk(mse)

    ## Yield diff
    yielddiff <- perform.yielddiff(mse)

    ## Return results
    res <- list(tb = tb, ssb = ssb, yield = yield, risk = risk, yielddiff = yielddiff)
    return(res)

}
