
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


#' @title Calculate production
#
#' @description Calculates production
#'
#' @param B biomass
#' @param pars List with parameters (r, K, n) plus m for Fletcher and fmsy for
#' @param type Either 1 for "Schaefer", 2 for  "Fletcher" or 3 for "Thorson". Default is 2 ("Fletcher").
#'
#' @details
#'
#' @export
calc.prod <- function(B, pars, type = 2){

    if(type == "Schaefer" || type == 1){
        r <- pars[["r"]]
        K <- pars[["K"]]
        sp <- r * B * (1 - B/K)
    }else if(type == "Fletcher" || type == 2){
        r <- pars[["r"]]
        K <- pars[["K"]]
        n <- pars[["n"]]
        m <- pars[["m"]]
        gamma <- n^(n/(n-1))/(n-1)
        sp <- gamma * m/K * B * (1 - (B/K)^(n-1))
    }else if(type == "Thorson" || type == 3){
        r <- pars[["r"]]
        K <- pars[["K"]]
        n <- pars[["n"]]
        fmsy <- pars[["fmsy"]]
        gamma <- n^(n/(n-1))/(n-1)
        bmsy_K <- (1/n)^(1/(n-1))
        sp <- gamma * bmsy_K * fmsy * B - gamma * bmsy_K * fmsy * K * (B/K)^n
    }else{
        stop("Only type = 1 (Fletcher) and type = 2 (Thorson) implemented!")
    }

    return(sp)
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


#' @title Convert ESB to SSB
#
#' @description Convert ESB to SSB
#'
#' @param ESB Exploitable stock biomass
#' @param F fishing mortality
#' @param pars conversion parameters
#'
#' @return Spawning stock biomass (SSB)
#'
#' @export
esb2ssb <- function(ESB, F, pars, exp = TRUE){

    if(!pars$log){
        esb <- pars[[1]][1] + ESB * pars[[1]][2] + F * pars[[1]][3]
    }else{
        esb <- pars[[1]][1] + log(ESB + 1e-10) * pars[[1]][2] + F * pars[[1]][3]
        if(exp){
            esb <- exp(esb)
        }
    }

    return(esb)
}


#' @title Estimate parameters for ESB to SSB conversion
#
#' @description Estimate parameters for ESB to SSB conversion
#'
#' @param data data frame with a column, ESB exploitable stock biomass, SSB spawning stock biomass, and F fishing mortality
#' @param plot Plot the results? Default: TRUE
#' @param log Fit ont log scale?
#' @param verbose Print info?
#'
#' @export
est.pars.esb2ssb <- function(data, plot = TRUE, log = FALSE, verbose = TRUE){

    if(!any(colnames(data) == "SSB")){
        data$SSB <- data$ssb
    }
    if(!any(colnames(data) == "ESB")){
        data$ESB <- data$esb
    }
    if(!any(colnames(data) == "F")){
        data$F <- data$fbar
    }

    if(log){
        func <- function(pars, ESB, SSB, F, log){
            pars <- c(list(pars), list(log = log))
            return(sqrt(mean((log(SSB + 1e-10) - esb2ssb(ESB, F, pars, exp = FALSE))^2)))
        }
    }else{
        func <- function(pars, ESB, SSB, F, log){
            pars <- c(list(pars), list(log = log))
            return(sqrt(mean((SSB - esb2ssb(ESB, F, pars))^2)))
        }
    }

    opt <- nlminb(c(0,0,0), func, SSB = data[,"SSB"], ESB = data[,"ESB"], F = data[,"F"], log = log)
    res <- c(pars = list(opt$par), list(log = log))

    if(verbose) writeLines(paste0(opt$message))

    if(plot){
        esb.plot <- seq(0, 1.2 * max(data$ESB, na.rm = TRUE), length.out = 1e3)
        ssb.plot <- seq(0, 1.2 * max(data$SSB, na.rm = TRUE), length.out = 1e3)
        fs <- quantile(data$F, probs = c(0.01,0.5,0.99), na.rm = TRUE)
        plot(esb.plot, ssb.plot, ty = "n",
             xlab = "ESB", ylab = "SSB")
        lines(esb.plot, esb.plot, col = "grey70", lty = 2, lwd = 2)
        for(i in 1:length(fs)){
            lines(esb.plot, esb2ssb(esb.plot, fs[i], res), lwd = 2, col = i+1)
        }
        points(data$ESB, data$SSB, col = ifelse(data$F >= fs[2], 4, 2))
        legend("bottomright", legend = c("F < median", "F >= median", "1-1",paste0("F = ", signif(fs,2))),
               col = c(2,4,"grey70", 1:length(fs) + 1), lty = c(NA,NA,2, rep(1,length(fs))),
               lwd = c(NA,NA,rep(2,length(fs)+1)), pch = c(1,1,rep(NA,length(fs)+1)))
    }

    return(res)
}



#' @title Run fmsy
#
#' @description Run fmsy
#'
#' @param mse HERE:
#' @param yr HERE:
#' @param fmsy HERE:
#' @param nsim HERE:
#' @param mc.cores HERE:
#'
#' @export
run.fmsy <- function(mse, yr = NULL, fmsy = seq(0,0.6,0.02), nsim = NULL, mc.cores = 1){

    dims <- dim(mse$est)
    if(is.null(yr)){
        years <- 1:dims[1]
    }else{
        years <- which(mse$years %in% yr)
    }
    if(is.null(nsim)) nsim <- mse$pars$nsim
    pars <- mse$pars

    pars$nsim <- nsim

    if(mc.cores > 1){

        tmp0 <- mclapply.all.os(as.list(1:length(fmsy)),
                      function(x, pars){
                          pars$refs[1] <- fmsy[x]
                          run.mse(pars, mc.cores = 1)
                      }, pars = pars, mc.cores = mc.cores)

    }else{

        tmp0 <- lapply(1:length(fmsy),
                      function(x, pars){
                          pars$refs[1] <- fmsy[x]
                          run.mse(pars, mc.cores = 1)
                      }, pars = pars)

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
#' @param mc.cores HERE:
#'
#' @export
run.explo <- function(mse, yr = NULL, fmsy = seq(0.1,1,0.05),
                      btrigger = seq(1000,5000,200),
                      nsim = NULL, mc.cores = 1){
    dims <- dim(mse$est)
    if(is.null(yr)){
        years <- 1:dims[1]
    }else{
        years <- which(mse$years %in% yr)
    }
    if(is.null(nsim)) nsim <- mse$pars$nsim
    pars <- mse$pars

    ## all scenarios (f & btrigger combis)
    scens <- expand.grid(fmsy, btrigger)
    nscen <- nrow(scens)

    pars$nsim <- nsim

    scenarios.par <- define.scenarios(pars)
    nscen.par <- nrow(scenarios.par)
    nms <- length(pars$HCR)

    ## Run
    ## if(mc.cores > 1){
    ##     tmp0 <- mclapply.all.os(as.list(1:nscen),
    ##                   function(x, pars){
    ##                       pars$refs[["Fmsy"]] <- scens[x,1]
    ##                       pars$refs[["Btrigger"]] <- scens[x,2]
    ##                       run.mse(pars, mc.cores = 1)
    ##                   }, pars = pars, mc.cores = mc.cores)
    ## }else{
        tmp0 <- lapply(1:nscen,
                      function(x, pars){
                          pars$refs[["Fmsy"]] <- scens[x,1]
                          pars$refs[["Btrigger"]] <- scens[x,2]
                          run.mse(pars, mc.cores = mc.cores)
                      }, pars = pars)
    ## }


    res <- array(NA, c(nms,nscen.par,nscen,dims[4]+2))
    for(ms in 1:nms){
        for(scen.par in 1:nscen.par){
            for(scen in 1:nscen){
                tmp <- tmp0[[scen]]
                res[ms,scen.par,scen,1:11] <- t(
                    apply(tmp$est[years,ms,scen.par,,,drop=FALSE],
                          c(2,4), mean, na.rm = TRUE))
                tmp2 <- t(apply(tmp$est[years,ms,scen.par,,,drop=FALSE],
                                c(2,4), sd, na.rm = TRUE))
                res[ms,scen.par,scen,9] <- tmp2[9,] * 100 / res[ms,scen.par,scen,4]  ## SE of tacVar * 100 / mean of yImp
                tmp2 <- t(apply(tmp$est[years,ms,scen.par,,,drop=FALSE],
                                c(2,4), quantile, probs = 0.05 , na.rm = TRUE))
                res[ms,scen.par,scen,12] <- tmp2[8,]
                tmp2 <- t(apply(tmp$est[years,ms,scen.par,,,drop=FALSE],
                                c(2,4), function(x){
                                    mean(x < pars$refs[3], na.rm = TRUE)
                }))
                res[ms,scen.par,scen,13] <- tmp2[7,]
            }
        }
    }

    ret <- vector("list",nms)
    names(ret) <- pars$HCR
    for(ms in 1:nms){
        ret[[ms]] <- vector("list",nscen.par)
        names(ret[[ms]]) <- paste0("scen_",seq(nscen.par))
        for(scen.par in 1:nscen.par){
            ret[[ms]][[scen.par]] <- cbind(scens,res[ms,scen.par,,])
            colnames(ret[[ms]][[scen.par]]) <- c("Fmsy","Btrigger",
                                                 names(mse$est[1,1,1,,1]),
                                                 "ssbObsQuant5","risk")
        }
    }

    ## for(i in 1:dim(res)[3]){
    ##     ret[[i]] <- cbind(scens,res[,,i])
    ##     colnames(ret[[i]]) <- c("Fmsy","Btrigger",names(mse$est[1,1,1,,1]),"ssbObsQuant5","risk")
    ## }

    return(ret)

}


#' @title Define scenarios
#
#' @description Define scenarios
#'
#' @param pars HERE:
#'
#' @export
define.scenarios <- function(pars){
    if(pars$prod.type == 1 || pars$prod.type == "Schaefer"){
        if(!any(names(pars) == "r")) stop("r not provided in pars!")
        if(!any(names(pars) == "K")) stop("K not provided in pars!")
        if(!any(names(pars) == "sdPro")) pars$sdPro <- 0.0
        if(!any(names(pars) == "sdObs")) pars$sdObs <- 0.0
        if(!any(names(pars) == "sdImp")) pars$sdImp <- 0.0
        scenarios <- expand.grid(r = pars$r, K = pars$K,
                                 sdPro = pars$sdPro, sdObs = pars$sdObs, sdImp = pars$sdImp)
    }else if(pars$prod.type == 2 || pars$prod.type == "Fletcher"){
        if(!any(names(pars) == "r")) stop("r not provided in pars!")
        if(!any(names(pars) == "K")) stop("K not provided in pars!")
        if(!any(names(pars) == "n")) stop("n not provided in pars!")
        if(!any(names(pars) == "m")) stop("m not provided in pars!")
        if(!any(names(pars) == "sdPro")) pars$sdPro <- 0.0
        if(!any(names(pars) == "sdObs")) pars$sdObs <- 0.0
        if(!any(names(pars) == "sdImp")) pars$sdImp <- 0.0
        scenarios <- expand.grid(r = pars$r, K = pars$K, n = pars$n, m = pars$m,
                                 sdPro = pars$sdPro, sdObs = pars$sdObs, sdImp = pars$sdImp)
    }else if(pars$prod.type == 3 || pars$prod.type == "Thorson"){
        if(!any(names(pars) == "r")) stop("r not provided in pars!")
        if(!any(names(pars) == "K")) stop("K not provided in pars!")
        if(!any(names(pars) == "n")) stop("n not provided in pars!")
        if(!any(names(pars) == "fmsy")) stop("fmsy not provided in pars!")
        if(!any(names(pars) == "sdPro")) pars$sdPro <- 0.0
        if(!any(names(pars) == "sdObs")) pars$sdObs <- 0.0
        if(!any(names(pars) == "sdImp")) pars$sdImp <- 0.0
        scenarios <- expand.grid(r = pars$r, K = pars$K, n = pars$n, Fmsy = pars$Fmsy,
                                 sdPro = pars$sdPro, sdObs = pars$sdObs, sdImp = pars$sdImp)
    }
    return(scenarios)
}


#' @title mclapply.windows
#' @description Alternative parallelisation for windows
#'
#' @importFrom parallel detectCores makeCluster clusterExport parLapply stopCluster
#'
#' @details Reference: https://www.r-bloggers.com/2014/07/implementing-mclapply-on-windows-a-primer-on-embarrassingly-parallel-computation-on-multicore-systems-with-r/
#'
#' @export
mclapply.windows <- function(...,mc.cores = parallel::detectCores()-1) {
    ## Create a cluster
    size.of.list <- length(list(...)[[1]])
    cl <- parallel::makeCluster(spec = min(size.of.list, mc.cores) )

    ## Find out the names of the loaded packages
    loaded.package.names <- c(
        ## Base packages
        sessionInfo()$basePkgs,
        ## Additional packages
        names( sessionInfo()$otherPkgs ))

    tryCatch( {
       ## Copy over all of the objects within scope to all clusters
       this.env <- environment()
       while( identical( this.env, globalenv() ) == FALSE ) {
           parallel::clusterExport(cl,
                         ls(all.names=TRUE, env=this.env),
                         envir=this.env)
           this.env <- parent.env(environment())
       }
       parallel::clusterExport(cl,
                     ls(all.names=TRUE, env=globalenv()),
                     envir=globalenv())

       ## Load the libraries on all the clusters
       ## N.B. length(cl) returns the number of clusters
       parallel::parLapply(
                     cl,
                     1:length(cl),
                     function(xx){
           lapply(loaded.package.names, function(yy) {
               require(yy , character.only=TRUE)})
       })

       ## Run the lapply in parallel
       return(parallel::parLapply( cl, ...))

    }, finally = {
       ## Stop the cluster
       parallel::stopCluster(cl)
    })
}


#' @title mclapply.all.os
#' @description mclapply comaptible with all OS
#'
#' @importFrom parallel mclapply
#'
#' @export
mclapply.all.os <- switch(
    Sys.info()[['sysname']],
   Windows = {mclapply.windows},
   Linux   = {parallel::mclapply},
   Darwin  = {parallel::mclapply}
)



#' @title Extract relevant parameters from fitted spict object
#'
#' @description Extract relevant parameters from fitted spict object
#'
#' @param fit Fitted spict object
#' @param verbose Print information
#'
#' @importFrom spict get.par
#'
#' @export
get.pars.spict <- function(fit, verbose = TRUE){
    pars <- list()

    ## Production parameters
    pars$r <- get.par("logr", fit, exp = TRUE)[,2]
    pars$K <- get.par("logK", fit, exp = TRUE)[,2]
    pars$n <- get.par("logn", fit, exp = TRUE)[,2]
    pars$m <- get.par("logm", fit, exp = TRUE)[,2]

    ## Noise parameters
    if(fit$inp$phases$logsdb == -1){
        if(verbose) writeLines("Process uncertainty was fixed! Setting to NA")
        pars$sdPro <- NA
    }else{
        pars$sdPro <- get.par("logsdb", fit, exp = TRUE)[,2]
    }
    if(fit$inp$phases$logsdi == -1){
        if(verbose) writeLines("Observation uncertainty was fixed! Setting to NA")
        pars$sdObs <- NA
    }else{
        pars$sdObs <- get.par("logsdi", fit, exp = TRUE)[,2]
    }


    ## Current states
    pars$Bcur <- get.par("logB", fit, exp = TRUE)[
        which.min(abs(fit$inp$time - ceiling(fit$inp$timerangeObs[2]))),2]
    pars$Fcur <- get.par("logFnotS", fit, exp = TRUE)[
        which.min(abs(fit$inp$time - ceiling(fit$inp$timerangeObs[2]))),2]

    ## Referece points
    pars$refs <- list()
    Bmsy <- get.par("Bmsy", fit)[,2]
    pars$refs$Fmsy <- get.par("Fmsy", fit)[,2]
    pars$refs$Bpa <- Bmsy/2
    pars$refs$Blim <- Bmsy/3
    pars$refs$Btrigger <- Bmsy/2

    return(pars)
}
