#' @title Run one MSE rep
#
#' @description Run one MSE rep
#'
#' @param x HERE:
#' @param pars HERE:
#' @param first.esb.is.obs logical, HERE:
#'
#' @export
run.mse.one.rep <- function(x, pars, first.esb.is.obs = TRUE){

    proyears <- pars$proyears

    ## Scenarios (r-K combinations)
    ## ---------------------------------
    scenarios <- define.scenarios(pars)
    nms <- length(pars$HCR)
    nscen <- nrow(scenarios)
    ssbFlag <- ifelse(any(names(pars) == "ssb.pars"), TRUE, FALSE)

    ## Make containers
    ## ---------------------------------
    tmp <- array(NA, c(proyears, nms, nscen))
    esb <- esbObs <- sp <- spPro <- y <- yImp <-
        ssb <- ssbObs <- tacVar <-
            Fint <- Frea <- tmp
    esb[1,,] <- esbObs[1,,] <- pars$Bcur
    Frea[1,,] <- pars$Fcur

    ## Noise (same for each scenario and MS)
    ## ---------------------------------
    noiseProList <- lapply(pars$sdPro, sim.log.noise, n = proyears) ## CHECK: why no log?
    noiseObsList <- lapply(pars$sdObs, sim.log.noise, n = proyears)
    noiseImpList <- lapply(pars$sdImp, sim.log.noise, n = proyears)

    ## Loop over scenarios, MS, and years
    ## ---------------------------------
    for(scen in 1:nscen){

        indPro <- match(scenarios[scen, "sdPro"], pars$sdPro)[1]
        indObs <- match(scenarios[scen, "sdObs"], pars$sdObs)[1]
        indImp <- match(scenarios[scen, "sdImp"], pars$sdImp)[1]


        for(ms in 1:nms){
            HCR <- pars$HCR[ms]
            if(HCR == "noF") Frea[1,ms,] <- 0.0
            for(year in 1:proyears){

                ## use current ESB as perceived ESB
                if(first.esb.is.obs && year == 1){
                    esb[year,ms,scen] <- add.log.noise(esbObs[year,ms,scen], noiseObsList[[indObs]][year])
                }

                ## ---------------------------------
                ## SSB as proportion of ESB dependent on F
                if(ssbFlag){
                    ssb[year,ms,scen] <- esb2ssb(esb[year,ms,scen], Frea[year,ms,scen], pars$ssb.pars)
                }

                ## ---------------------------------
                ## Surplus production (sp)
                sp[year,ms,scen] <- calc.prod(esb[year,ms,scen], as.list(scenarios[scen,]), type = pars$prod.type)
                ## White process noise on surplus production (environment, interspecific interactions, etc.)
                spPro[year,ms,scen] <- add.log.noise(sp[year,ms,scen], noiseProList[[indPro]][year])
                ## TODO: if sdPro == CV:
                ## spPro[year,ms,scen] <- add.noise(sp[year,ms,scen],
                ##                                  sim.noise(scenarios[scen,5] * esb[year,ms,scen], n = 1))

                ## ---------------------------------
                ## White observation noise on ESB
                if(!first.esb.is.obs || year > 1){
                    esbObs[year,ms,scen] <- add.log.noise(esb[year,ms,scen], noiseObsList[[indObs]][year])
                }
                ## SSB as proportion of ESB dependent on F
                if(ssbFlag){
                    ssbObs[year,ms,scen] <- esb2ssb(esbObs[year,ms,scen], Frea[year,ms,scen], pars$ssb.pars)
                }


                ## Simulated management based on HCR
                ## ---------------------------------
                ## F according to MSE
                if(pars$biomass.type.hcr[ms] == "SSB"){
                    if(!ssbFlag) stop("'pars$biomass.type.hcr' chosen, but no conversion parameters ('pars$ssb.pars') provided!")

                    ## Find biomass next year using HCR
                    yimp.cur <- add.log.noise(esbObs[year,ms,scen] *
                                      apply.hcr(HCR, biomass = ssbObs[year,ms,scen], refs = pars$refs),
                                      noiseImpList[[indImp]][year])
                    esb.next <- esb[year,ms,scen] + spPro[year,ms,scen] - yimp.cur
                    esb.next <- ifelse(esb.next < 0, 0.001, esb.next)
                    if(!first.esb.is.obs || year > 1){
                        esb.next <- add.log.noise(esb.next, noiseObsList[[indObs]][year])
                    }
                    ssb.pred <- esb2ssb(esb.next, yimp.cur/esb[year,ms,scen], pars$ssb.pars)

                    ## Account for blim and predicted biomass
                    if(HCR == "ICES3" && ssb.pred < pars$refs$Blim){

                        ## which biomass if F = 0, if still below blim -> F = 0
                        yimp.cur <- add.log.noise(esbObs[year,ms,scen] *
                                                  apply.hcr("noF", biomass = ssbObs[year,ms,scen], refs = pars$refs),
                                                  noiseImpList[[indImp]][year])
                        esb.next <- esb[year,ms,scen] + spPro[year,ms,scen] - yimp.cur
                        esb.next <- ifelse(esb.next < 0, 0.001, esb.next)
                        if(!first.esb.is.obs || year > 1){
                            esb.next <- add.log.noise(esb.next, noiseObsList[[indObs]][year])
                        }
                        ssb.pred <- esb2ssb(esb.next, yimp.cur/esb[year,ms,scen], pars$ssb.pars)

                        if(ssb.pred < pars$refs$Blim){
                            Fint[year,ms,scen] <- 0
                        }else{
                            ## otherwise highest F that results in b > blim
                            func <- function(f){
                                ## TODO: remove implementation noise?
                                yimp.cur <- add.log.noise(esbObs[year,ms,scen] * f,
                                                          noiseImpList[[indImp]][year])
                                esb.next <- esb[year,ms,scen] + spPro[year,ms,scen] - yimp.cur
                                esb.next <- ifelse(esb.next < 0, 0.001, esb.next)
                                if(!first.esb.is.obs || year > 1){
                                    esb.next <- add.log.noise(esb.next, noiseObsList[[indObs]][year])
                                }
                                ssb.pred <- esb2ssb(esb.next, yimp.cur/esb[year,ms,scen], pars$ssb.pars)
                                ## ifelse inside optimise not ideal
                                return(f * ifelse(pars$refs$Blim - ssb.pred < 0, 1, -1))
                            }
                            opt <- optimise(func, c(1e-6, 10), maximum = TRUE, tol = 0.001)
                            Fint[year,ms,scen] <- opt$maximum
                        }
                    }else{
                        Fint[year,ms,scen] <- apply.hcr(HCR, biomass = ssbObs[year,ms,scen],
                                                        biomass.blim = ssb.pred,
                                                        refs = pars$refs)
                    }

                }else if(pars$biomass.type.hcr[ms] == "ESB"){
                    ## TODO: implement ICES3 here!
                    ## Find biomass next year using HCR
                    yimp.cur <- add.log.noise(esbObs[year,ms,scen] *
                                      apply.hcr(HCR, biomass = ssbObs[year,ms,scen], refs = pars$refs),
                                      noiseImpList[[indImp]][year])
                    esb.next <- esb[year,ms,scen] + spPro[year,ms,scen] - yimp.cur
                    esb.next <- ifelse(esb.next < 0, 0.001, esb.next)
                    if(!first.esb.is.obs || year > 1){
                        esb.next <- add.log.noise(esb.next, noiseObsList[[indObs]][year])
                    }
                    Fint[year,ms,scen] <- apply.hcr(HCR,
                                                    biomass = esbObs[year,ms,scen],
                                                    biomass.blim = esb.next,
                                                    refs = pars$refs)
                }
                ## Yield based on target F
                y[year,ms,scen] <- esbObs[year,ms,scen] * Fint[year,ms,scen]
                ## White implementation noise on yield
                yImp[year,ms,scen] <- add.log.noise(y[year,ms,scen], noiseImpList[[indImp]][year])

                ## Next years biomass
                ## ---------------------------------
                if(year != proyears){
                    tmp <- esb[year,ms,scen] + spPro[year,ms,scen] - yImp[year,ms,scen]
                    esb[year+1,ms,scen] <- ifelse(tmp < 0, 0.001, tmp)
                    Frea[year+1,ms,scen] <- yImp[year,ms,scen] / esb[year,ms,scen]
                }
            }
            ## TAC difference
            ## ---------------------------------
            tacVar[,ms,scen] <- c(y[-1,ms,scen],NaN) - y[,ms,scen]
        }
    }

    ## Combine results
    ## ---------------------------------
    if(ssbFlag){
        res <- abind::abind(esb = esb, esbObs = esbObs,
                            y = y, yImp = yImp, sp = sp, spPro = spPro,
                            ssb = ssb, ssbObs = ssbObs,
                            tacVar = tacVar, Fint = Fint,
                            Frea = Frea, along = 4)
    }else{
        res <- abind::abind(esb = esb, esbObs = esbObs,
                            y = y, yImp = yImp, sp = sp, spPro = spPro,
                            tacVar = tacVar, Fint = Fint,
                            Frea = Frea, along = 4)
    }
    return(res)
}


#' @title Run MSE
#
#' @description Run MSE
#'
#' @param pars HERE:
#' @param mc.cores HERE:
#' @param verbose Print information
#'
#' @export
run.mse <- function(pars, mc.cores = 1, verbose = TRUE){

    ## Check input
    ## ---------------------------------
    if(!any(names(pars) == "HCR")) pars$HCR <- c("noF","FMSY","ICES")
    ## if(!any(names(pars) == "sigmas")) pars$sigmas <- c(0.15,0.15,0.001)
    ## if(length(pars$sigmas) < 3) stop(paste0("Only ",length(pars$sigmas)," sigmas provided, but 3 sigmas required: biomass process error, biomass observation error, and TAC implementation error!"))
    if(!any(names(pars) == "Bcur")) stop("Bcur missing in pars!")
    if(!any(names(pars) == "Fcur")) stop("Fcur missing in pars!")
    if(!any(names(pars) == "refs")) stop("refs missing in pars!")
    ## if(!any(names(pars) == "matpars")) stop("matpars missing in pars!")
    if(!any(names(pars) == "r")) stop("r missing in pars!")
    if(!any(names(pars) == "K")) stop("K missing in pars!")

    if(!any(names(pars) == "biomass.type.hcr")) stop("Please provide information which biomass type should be used for the HCR ('pars$biomass.type.hcr'), choose either 'SSB' or 'ESB'. Note, that this should correspond to the biomass type of the biomass reference points in 'pars$refs'!")
    if(length(pars$biomass.type.hcr) != length(pars$HCR) && length(pars$biomass.type.hcr) == 1){
        if(verbose) writeLines(paste0("Only one biomass.type.hcr specified, but ",length(pars$HCR), " HCRs specified. Assuming the same biomass.type for all HCRs!"))
        pars$biomass.type.hcr <- rep(pars$biomass.type.hcr, length(pars$HCR))
    }else if(length(pars$biomass.type.hcr) != length(pars$HCR)){
        stop("Length of biomass types (pars$biomass.type.hcr) does not correspond to the numbers of HCRs (pars$HCR). Please check!")
    }

    if(!any(names(pars) == "proyears")){
        if(verbose) writeLines("Number of projection years ('pars$proyears') not specified. Setting to 5!")
        pars$proyears <- 5
    }
    if(!any(names(pars) == "nrep")){
        if(verbose) writeLines("Number of replicates ('pars$nrep') not specified. Setting to 10!")
        pars$nrep <- 10
    }

    proyears <- pars$proyears
    nrep <- pars$nrep

    ## Scenarios (r-K combinations)
    ## ---------------------------------
    scenarios <- define.scenarios(pars)


    ## Run mse
    ## ---------------------------------
    if(mc.cores > 1){
        res <- mclapply.all.os(as.list(1:nrep), function(x){
            run.mse.one.rep(x, pars, proyears)
        }, mc.cores = mc.cores)
    }else{
        res <- lapply(1:nrep,
                      function(x) run.mse.one.rep(x, pars, proyears))
    }

    ## Rearrange results to array with dim: year, MS, scenario, quantity, nrep
    est <- do.call(abind::abind, c(res, along = 5))

    ## Return results
    ret <- list(pars = pars, scenarios = scenarios,
                est = est, years = pars$yearcur:(pars$yearcur+proyears-1))

    return(ret)
}
