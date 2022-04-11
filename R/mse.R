#' @title Run one MSE rep
#
#' @description Run one MSE rep
#'
#' @param x HERE:
#' @param pars HERE:
#' @param proyears HERE:
#' @param first.tb.is.obs logical, HERE:
#'
#' @export
run.mse.one.rep <- function(x, pars, proyears = 5, first.tb.is.obs = TRUE){

    ## Scenarios (r-K combinations)
    ## ---------------------------------
    scenarios <- expand.grid(r = pars$r, K = pars$K, n = pars$n, Fmsy = pars$Fmsy)
    nms <- length(pars$MStypes)
    nscen <- nrow(scenarios)

    ## Make containers
    ## ---------------------------------
    tmp <- array(NA, c(proyears, nms, nscen))
    tb <- tbObs <- sp <- spPro <- y <- yImp <-
        ssb <- ssbObs <- tacVar <-
            Fint <- Frea <- tmp
    tb[1,,] <- tbObs[1,,] <- pars$Bcur
    Frea[1,,] <- pars$Fcur

    ## Noise (same for each scenario and MS)
    ## ---------------------------------
    noisePro <- sim.noise(pars$sigmas[1], proyears)
    noiseObs <- sim.log.noise(pars$sigmas[2], proyears)
    noiseImp <- sim.log.noise(pars$sigmas[3], proyears)

    ## Loop over scenarios, MS, and years
    ## ---------------------------------
    for(scen in 1:nscen){
        for(ms in 1:nms){
            mstype <- pars$MStypes[ms]
            if(mstype == "noF") Frea[1,ms,] <- 0.0
            for(year in 1:proyears){

                ## use current TB as perceived TB
                if(first.tb.is.obs && year == 1){
                    tb[year,ms,scen] <- add.log.noise(tbObs[year,ms,scen], noiseObs[year])
                }

                ## ---------------------------------
                ## SSB as proportion of TB dependent on F
                ssb[year,ms,scen] <- calc.ssb(tb[year,ms,scen], Frea[year,ms,scen], pars$matpars)

                ## ---------------------------------
                ## Surplus production (sp)
                sp[year,ms,scen] <- calc.prod.thorson(tb[year,ms,scen], scenarios[scen,1], scenarios[scen,2],
                                                      scenarios[scen,3], scenarios[scen,4])
                ## White process noise on surplus production (environment, interspecific interactions, etc.)
                ## spPro[year,ms,scen] <- add.noise(sp[year,ms,scen], noisePro[year])
                spPro[year,ms,scen] <- add.noise(sp[year,ms,scen], sim.noise(pars$sigmas[1] * tb[year,ms,scen], n = 1)) ## not reproducible!

                ## ---------------------------------
                ## White observation noise on TB
                if(!first.tb.is.obs || year > 1) tbObs[year,ms,scen] <- add.log.noise(tb[year,ms,scen], noiseObs[year])
                ## SSB as proportion of TB dependent on F
                ssbObs[year,ms,scen] <- calc.ssb(tbObs[year,ms,scen], Frea[year,ms,scen], pars$matpars)

                ## Simulated management based on HCR
                ## ---------------------------------
                ## F according to MSE
                Fint[year,ms,scen] <- apply.ms(pars$refs[1], mstype,
                                               SSB = ssbObs[year,ms,scen], Bpa = pars$refs[2])
                ## Yield based on target F
                y[year,ms,scen] <- tbObs[year,ms,scen] * Fint[year,ms,scen]
                ## White implementation noise on yield
                yImp[year,ms,scen] <- add.log.noise(y[year,ms,scen], noiseImp[year])

                ## Next years biomass
                ## ---------------------------------
                if(year != proyears){
                    tmp <- tb[year,ms,scen] + spPro[year,ms,scen] - yImp[year,ms,scen]
                    tb[year+1,ms,scen] <- ifelse(tmp < 0, 0.001, tmp)
                    Frea[year+1,ms,scen] <- yImp[year,ms,scen] / tb[year,ms,scen]
                }
            }
            ## TAC difference
            ## ---------------------------------
            tacVar[,ms,scen] <- c(y[-1,ms,scen],NaN) - y[,ms,scen]
        }
    }

    ## Combine results
    ## ---------------------------------
    res <- abind::abind(tb = tb, tbObs = tbObs,
                y = y, yImp = yImp, sp = sp, spPro = spPro,
                ssb = ssb, ssbObs = ssbObs,
                tacVar = tacVar, Fint = Fint,
                Frea = Frea, along = 4)
    return(res)
}


#' @title Run MSE
#
#' @description Run MSE
#'
#' @param pars HERE:
#' @param proyears HERE:
#' @param nsim HERE:
#' @param ncores HERE:
#'
#' @export
run.mse <- function(pars, proyears = 5, nsim = 10, ncores = 1){

    ## Check input
    ## ---------------------------------
    if(!any(names(pars) == "MStypes")) pars$MStypes <- c("noF","FMSY","ICES")
    if(!any(names(pars) == "sigmas")) pars$sigmas <- c(0.15,0.15,0.001)
    if(length(pars$sigmas) < 3) stop(paste0("Only ",length(pars$sigmas)," sigmas provided, but 3 sigmas required: biomass process error, biomass observation error, and TAC implementation error!"))
    if(!any(names(pars) == "Bcur")) stop("Bcur missing in pars!")
    if(!any(names(pars) == "Fcur")) stop("Fcur missing in pars!")
    if(!any(names(pars) == "refs")) stop("refs missing in pars!")
    if(!any(names(pars) == "matpars")) stop("matpars missing in pars!")
    if(!any(names(pars) == "r")) stop("r missing in pars!")
    if(!any(names(pars) == "K")) stop("K missing in pars!")

    ## Scenarios (r-K combinations)
    ## ---------------------------------
    scenarios <- expand.grid(r = pars$r, K = pars$K, n = pars$n, Fmsy = pars$Fmsy)

    ## Parallel setup (needed for windows)
    ## ---------------------------------
    if(ncores > 1){
        if(snowfall::sfIsRunning()) snowfall::sfStop()

        snowfall::sfInit(parallel=TRUE, cpus = ncores)
        ex <- Filter(function(x) is.function(get(x, .GlobalEnv)), ls(.GlobalEnv))
        snowfall::sfExport(list=ex)

##         snowfall::sfLibrary(pk, character.only = TRUE, verbose=FALSE)

        res <- snowfall::sfClusterApplyLB(1:nsim, run.mse.one.rep, pars = pars,
                                          proyears = proyears)

        sfStop()

        ## older (only works on UNIX)
        ## res <- parallel::mclapply(1:nsim,
        ##                           function(x) run.mse.one.rep(pars, proyears),
        ##                           mc.cores = ncores)


    }else{
        res <- lapply(1:nsim,
                      function(x) run.mse.one.rep(x, pars, proyears))
    }

    ## Rearrange results to array with dim: year, MS, scenario, quantity, nsim
    est <- do.call(abind::abind, c(res, along = 5))

    ## Return results
    ret <- list(pars = pars, scenarios = scenarios,
                proyears = proyears, nsim = nsim,
                est = est, years = pars$yearcur:(pars$yearcur+proyears-1))

    return(ret)
}
