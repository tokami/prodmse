
#' @title Calculate Standard Error (SE)
#
#' @description Calculates Standard Error (SE)
#'
#' @param x numeric
#' @param na.rm logical; delete NAs?
#'
#' @export
se <- function(x, na.rm = FALSE){
    n <- ifelse(na.rm, length(na.omit(x)), length(x))
    sqrt(var(x, na.rm=na.rm) / n)
}


#' @title Calculate production based on Schaefer
#
#' @description Calculates production based on the Schaefer model
#'
#' @param B biomass
#' @param r intrinsic growth rate
#' @param K carrying capacity
#'
#' @export
calc.prod <- function(B, r, K){
    B + B * r * (1 - B/K)
}

#' @title Calculate production based on Thorson
#
#' @description Calculates production based on the Thorson model
#'
#' @param B biomass
#' @param r intrinsic growth rate
#' @param K carrying capacity
#'
#' @export
calc.prod.thorson <- function(B, r, K, n, fmsy){
    gamma <- n^(n/(n-1))/(n-1)
    bmsy_K <- (1/n)^(1/(n-1))
    gamma * bmsy_K * fmsy * B - gamma * bmsy_K * fmsy * K * (B/K)^n
}


#' @title Simulate log normal noise
#
#' @description Simulate log normal noise
#'
#' @param sigma standard deviation
#' @param n Number of observations
#'
#' @export
sim.log.noise <- function(sigma, n = 1){
    rnorm(n,0,sigma) - sigma^2/2
    ## CHECK: no exp here?
}


#' @title Simulate normal noise
#
#' @description Simulate normal noise
#'
#' @param sigma standard deviation
#' @param n Number of observations
#'
#' @export
sim.noise <- function(sigma, n = 1){
    rnorm(n,0,sigma)
}


#' @title Add log normal noise to value
#
#' @description Add log normal noise to value
#'
#' @param val value
#' @param noise noise
#'
#' @export
add.log.noise <- function(val, noise){
    exp(log(val) + noise)
}


#' @title Add normal noise to value
#
#' @description Add normal noise to value
#'
#' @param val value
#' @param noise noise
#'
#' @export
add.noise <- function(val, noise){
    val + noise
}


#' @title Calculate SSB
#
#' @description Calculate SSB
#'
#' @param B biomass
#' @param F fishing mortality
#' @param matpars maturity parameters
#'
#' @export
calc.ssb <- function(B, F, matpars){
    ## B / (( matpars[1] * F * matpars[2]/matpars[3]) + matpars[4]) ## version 1
    B * (matpars[1] * F + matpars[2])
}



#' @title Run fmsy
#
#' @description Run fmsy
#'
#' @param mse HERE:
#' @param yr HERE:
#' @param fmsy HERE:
#' @param nsim HERE:
#' @param ncores HERE:
#'
#' @export
run.fmsy <- function(mse, yr = NULL, fmsy = seq(0,0.6,0.02), nsim = NULL, ncores = 1){
    dims <- dim(mse$est)
    if(is.null(yr)){
        years <- 1:dims[1]
    }else{
        years <- which(mse$years %in% yr)
    }
    if(is.null(nsim)) nsim <- mse$nsim
    pars <- mse$pars


    if(ncores > 1){
        if(snowfall::sfIsRunning()) snowfall::sfStop()

        snowfall::sfInit(parallel=TRUE, cpus = ncores)
        ex <- Filter(function(x) is.function(get(x, .GlobalEnv)), ls(.GlobalEnv))
        snowfall::sfExport(list=ex)
##         snowfall::sfLibrary(pk, character.only = TRUE, verbose=FALSE)

        tmp0 <- snowfall::sfClusterApplyLB(1:length(fmsy),
                                          function(x, pars, proyears, nsim){
                                              pars$refs[1] <- fmsy[x]
                                              run.mse(pars, proyears, nsim, ncores = 1)
                                          }, pars = pars, proyears = mse$proyears, nsim = nsim)

        sfStop()

    }else{

        tmp0 <- lapply(1:length(fmsy),
                      function(x, pars, proyears, nsim){
                          pars$refs[1] <- fmsy[x]
                          run.mse(pars, proyears, nsim, ncores = 1)
                      }, pars = pars, proyears = mse$proyears, nsim = nsim)

    }

    res <- array(NA, c(length(fmsy),dims[4]+1,dims[2]))
    for(fs in 1:length(fmsy)){
        tmp <- tmp0[[fs]]
        res[fs,1:11,] <- t(apply(tmp$est[years,,,,,drop=FALSE], c(2,4), mean, na.rm = TRUE))
        tmp2 <- t(apply(tmp$est[years,,,,,drop=FALSE], c(2,4), sd, na.rm = TRUE))
        res[fs,9,] <- tmp2[9,] * 100 / res[fs,4,]  ## SE of tacVar * 100 / mean of yImp
        tmp2 <- t(apply(tmp$est[years,,,,,drop=FALSE], c(2,4), quantile, probs = 0.05 , na.rm = TRUE))
        res[fs,12,] <- tmp2[8,]
    }

    ret <- list()
    for(i in 1:dim(res)[3]){
        ret[[i]] <- cbind(fmsy,res[,,i])
        colnames(ret[[i]]) <- c("Fmsy",names(mse$est[1,1,1,,1]),"ssbObsQuant5")
    }

    return(ret)

}



#' @title Run explo
#
#' @description Run explo
#'
#' @param mse HERE:
#' @param yr HERE:
#' @param fmsy HERE:
#' @param btrigger HERE:
#' @param nsim HERE:
#' @param ncores HERE:
#'
#' @export
run.explo <- function(mse, yr = NULL, fmsy = seq(0.1,1,0.05),
                      btrigger = seq(1000,5000,200),
                      nsim = NULL, ncores = 1){
    dims <- dim(mse$est)
    if(is.null(yr)){
        years <- 1:dims[1]
    }else{
        years <- which(mse$years %in% yr)
    }
    if(is.null(nsim)) nsim <- mse$nsim
    pars <- mse$pars

    ## all scenarios (f & btrigger combis)
    scens <- expand.grid(fmsy, btrigger)
    nscen <- nrow(scens)


    if(ncores > 1){
        if(snowfall::sfIsRunning()) snowfall::sfStop()

        snowfall::sfInit(parallel=TRUE, cpus = ncores)
        ex <- Filter(function(x) is.function(get(x, .GlobalEnv)), ls(.GlobalEnv))
        snowfall::sfExport(list=ex)
##         snowfall::sfLibrary(pk, character.only = TRUE, verbose=FALSE)

        tmp0 <- snowfall::sfClusterApplyLB(1:nscen,
                                          function(x, pars, proyears, nsim){
                                              pars$refs[1] <- scens[x,1]
                                              pars$refs[2] <- scens[x,2]
                                              run.mse(pars, proyears, nsim, ncores = 1)
                                          }, pars = pars, proyears = mse$proyears,
                                          nsim = nsim)

        sfStop()

    }else{

        tmp0 <- lapply(1:length(fmsy),
                      function(x, pars, proyears, nsim){
                          pars$refs[1] <- fmsy[x]
                          run.mse(pars, proyears, nsim, ncores = 1)
                      }, pars = pars, proyears = mse$proyears, nsim = nsim)

    }

    res <- array(NA, c(nscen,dims[4]+2,dims[2]))
    for(scen in 1:nscen){
        tmp <- tmp0[[scen]]
        res[scen,1:11,] <- t(apply(tmp$est[years,,,,,drop=FALSE], c(2,4), mean, na.rm = TRUE))
        tmp2 <- t(apply(tmp$est[years,,,,,drop=FALSE], c(2,4), sd, na.rm = TRUE))
        res[scen,9,] <- tmp2[9,] * 100 / res[scen,4,]  ## SE of tacVar * 100 / mean of yImp
        tmp2 <- t(apply(tmp$est[years,,,,,drop=FALSE], c(2,4), quantile, probs = 0.05 , na.rm = TRUE))
        res[scen,12,] <- tmp2[8,]
        tmp2 <- t(apply(tmp$est[years,,,,,drop=FALSE], c(2,4), function(x){
            mean(x < pars$refs[3], na.rm = TRUE)

        }))
        res[scen,13,] <- tmp2[7,]
    }

    ret <- list()
    for(i in 1:dim(res)[3]){
        ret[[i]] <- cbind(scens,res[,,i])
        colnames(ret[[i]]) <- c("Fmsy","Btrigger",names(mse$est[1,1,1,,1]),"ssbObsQuant5","risk")
    }

    return(ret)

}
