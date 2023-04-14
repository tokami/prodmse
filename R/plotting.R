
#' @title Plot tradeoff
#
#' @description Plot tradeoff
#'
#' @param mse HERE:
#'
#' @export
plotprodmse.td <- function(mse){

    mse <- res

    ##
    nscen <- nrow(mse$scenarios)
    nms <- length(mse$pars$HCR)
    nall <- nscen * nms
    Bpa <- mse$pars$ref[2]
    Blim <- mse$pars$ref[3]

    ##
    cols <- c("darkgreen","dodgerblue","darkorange","purple","red")
    tmp <- nscen / 3
    if(tmp %% 1 == 0){
        parrow <- 3
        parcol <-tmp
    }else{
        tmp <- nscen / 4
        if(tmp %% 1 == 0){
            parrow <- 4
            parcol <- tmp
        }else{
            tmp <- nscen / 5
            if(tmp %% 1 == 0){
                parrow <- 5
                parcol <- tmp
            }else{
                tmp <- nscen / 2
                if(tmp %% 1 == 0){
                    parrow <- 2
                    parcol <- tmp
                }else{
                    parrow <- 1
                    parcol <- nscen / 1
                }
            }
        }
    }

    ## plot
    opar <- par()
    par(mfrow = c(parrow,parcol), mar = c(2,1,2,0), oma=c(3,5,2,3))
    tmp <- mse$est[,,,"y",]
    if(length(dim(tmp)) == 4) dimen <- 2:3 else dimen <- 2
    yield <- apply(tmp, dimen, function(x) mean(x, na.rm=TRUE))
    tmp <- mse$est[,,,"ssb",]
    if(length(dim(tmp)) == 4) dimen <- 2:3 else dimen <- 2
    risk <- apply(tmp, dimen, function(x) mean(x < Blim, na.rm = TRUE))

    ylim <- range(yield)
    xlim <- c(0,1.4) * range(risk, 0.051)
    for(scen in 1:nscen){
        if(length(dim(tmp)) == 4){
            xplot <- risk[,scen]
            yplot <- yield[,scen]
            pchs <- seq(nrow(yield))
        } else {
            xplot <- risk
            yplot <- yield
            pchs <- seq(length(yield))
        }
        plot(xplot, yplot,
             ty= "n", ylim = ylim, xlim = xlim,
             xaxt = "n", xlab = "",
             yaxt = "n", ylab = "")
        points(xplot, yplot, cex = 1.8,
               col = cols[1:nms],
               pch = pchs, lwd = 1.5)
        abline(v = 0.05, lty = 2, lwd = 1.5)
        if(scen %in% (nscen-parcol+1):nscen) axis(1)
        if(scen %in% seq(1,nscen,by=parcol)) axis(2)
        if(scen == 1) mtext("Yield", 2, 2.5, outer = TRUE)
        mtext(paste0("Scenario ", scen), 3, 1, font = 2)
        if(scen == nscen) legend("topright", legend = mse$pars$HCR,
                                 col = cols[1:nms], bg="white",
                                 pch = pchs)
    }
    mtext("Risk", 1, 1
        , outer=TRUE)

}


#' @title Plot time series
#
#' @description Plot time series
#'
#' @param mse HERE:
#'
#' @export
plotprodmse.ts <- function(mse){

    ##
    nscen <- nrow(mse$scenarios)
    nms <- length(mse$pars$HCR)
    nall <- nscen * nms
    Bpa <- mse$pars$ref[2]
    Blim <- mse$pars$ref[3]


    ##
    years <- 1:mse$pars$proyears
    cols <- c("darkgreen","dodgerblue","darkorange","purple","red")

    ## plot
    opar <- par()
    par(mfrow = c(4,nscen), mar = c(1,1,0,0), oma=c(5,5,4,3))
    ## biomass
    tmp <- mse$est[,,,"ssb",]
    if(length(dim(tmp)) == 4) dimen <- 1:3 else dimen <- 1:2
    tmp2 <- apply(tmp, dimen, mean, na.rm = TRUE)
    ylim <- c(0.95,1.3) * range(tmp2,Bpa,Blim, na.rm = TRUE)
    for(scen in 1:nscen){
        if(length(dim(tmp)) == 4) yplot <- tmp2[,1,scen] else yplot <- tmp2[,1]
        plot(years, yplot, ylim = ylim,
             ty="n",
             xaxt = "n", xlab = "",
             yaxt = "n", ylab = "")
        for(ms in 1:nms){
            if(length(dim(tmp)) == 4) yplot <- tmp2[,ms,scen] else yplot <- tmp2[,ms]
            lines(years, yplot, col = cols[ms], lwd = 1.5)
        }
        if(scen == 1) axis(2)
        if(scen == 1) mtext("SSB", 2, 3.5)
        abline(h = Blim, lty=1, col = "grey70")
        mtext(paste0("Scenario ", scen), 3, 1, font = 2)
        if(scen == nscen) legend("topright", legend = mse$pars$HCR,
                                 col = cols[1:nms], bg="white", lwd=1.5)
    }
    ## yield
    tmp <- mse$est[,,,"y",]
    if(length(dim(tmp)) == 4) dimen <- 1:3 else dimen <- 1:2
    tmp2 <- apply(tmp, dimen, mean, na.rm = TRUE)
    ylim <- c(0.95,1.3) * range(tmp2, na.rm = TRUE)
    for(scen in 1:nscen){
        if(length(dim(tmp)) == 4) yplot <- tmp2[,1,scen] else yplot <- tmp2[,1]
        plot(years, yplot, ylim = ylim,
             ty="n",
             xaxt = "n", xlab = "",
             yaxt = "n", ylab = "")
        for(ms in 1:nms){
            if(length(dim(tmp)) == 4) yplot <- tmp2[,ms,scen] else yplot <- tmp2[,ms]
            lines(years, yplot, col = cols[ms], lwd = 1.5)
        }
        if(scen == 1) axis(2)
        if(scen == 1) mtext("Yield", 2, 3.5)
    }
    ## risk
    tmp <- mse$est[,,,"ssb",]
    if(length(dim(tmp)) == 4) dimen <- 1:3 else dimen <- 1:2
    tmp2 <- apply(tmp, dimen, function(x) mean(x < Blim , na.rm = TRUE))
    ylim <- c(0.99,1.01) * range(tmp2,0,0.05,1, na.rm = TRUE)
    for(scen in 1:nscen){
        if(length(dim(tmp)) == 4) yplot <- tmp2[,1,scen] else yplot <- tmp2[,1]
        plot(years, yplot, ylim = ylim,
             ty="n",
             xaxt = "n", xlab = "",
             yaxt = "n", ylab = "")
        for(ms in 1:nms){
            if(length(dim(tmp)) == 4) yplot <- tmp2[,ms,scen] else yplot <- tmp2[,ms]
            lines(years, yplot, col = cols[ms], lwd = 1.5)
        }
        if(scen == 1) axis(2)
        if(scen == 1) mtext("Risk", 2, 3.5)
        abline(h = 0.05, lty=1, col = "grey70")
    }
    ## yield diff
    tmp <- mse$est[,,,"y",]
    if(length(dim(tmp)) == 4) dimen <- 2:4 else dimen <- 2:3
    tmp1 <- apply(tmp, dimen, diff)
    if(length(dim(tmp)) == 4) dimen <- 1:3 else dimen <- 1:2
    tmp2 <- apply(tmp1, dimen, mean , na.rm = TRUE)
    ylim <- c(0.95,1.3) * range(tmp2,0.0,na.rm = TRUE)
    for(scen in 1:nscen){
        if(length(dim(tmp)) == 4) yplot <- tmp2[,1,scen] else yplot <- tmp2[,1]
        plot(years, c(yplot,NA), ylim = ylim,
             ty="n",
             xaxt = "n", xlab = "",
             yaxt = "n", ylab = "")
        for(ms in 1:nms){
            if(length(dim(tmp)) == 4) yplot <- tmp2[,ms,scen] else yplot <- tmp2[,ms]
            lines(years, c(yplot,NA), col = cols[ms], lwd = 1.5)
        }
        if(scen == 1) axis(2)
        if(scen == 1) mtext("Difference in yield", 2, 3.5)
        abline(h = 0.0, lty=1, col = "grey70")
        axis(1)
    }
    mtext("Projection years", 1, 3, outer = TRUE)

}



#' @title Plot distribution
#
#' @description Plot distribution
#'
#' @param mse HERE:
#'
#' @export
plotprodmse.dist <- function(mse, outline = TRUE){

    ##
    nscen <- nrow(mse$scenarios)
    nms <- length(mse$pars$HCR)
    nall <- nscen * nms
    Bpa <- mse$pars$ref[2]
    Blim <- mse$pars$ref[3]

    ##
    cols <- c("darkgreen","dodgerblue","darkorange","purple","red")

    ## plot
    opar <- par()
    par(mfrow = c(4,nscen), mar = c(1,1,0,0), oma=c(5,5,4,3))
    ## biomass
    tmp <- mse$est[,,,"ssb",]
    ylim <- range(tmp, na.rm = TRUE)
    for(scen in 1:nscen){
        if(length(dim(tmp)) == 4) yplot <- tmp[,,scen,] else yplot <- tmp
        tmp2 <- apply(yplot,2,rbind)
        boxplot(tmp2, ylim = ylim,
                col = cols[1:nms],
                outline = outline,
                xaxt = "n", xlab = "",
                yaxt = "n", ylab = "")
        if(scen == 1) axis(2)
        if(scen == 1) mtext("SSB", 2, 3.5)
        mtext(paste0("Scenario ", scen), 3, 1, font = 2)
    }
    ## yield
    tmp <- mse$est[,,,"y",]
    ylim <- range(tmp, na.rm = TRUE)
    for(scen in 1:nscen){
        if(length(dim(tmp)) == 4) yplot <- tmp[,,scen,] else yplot <- tmp
        tmp2 <- apply(yplot,2,rbind)
        boxplot(tmp2, ylim = ylim,
                col = cols[1:nms],
                outline = outline,
                xaxt = "n", xlab = "",
                yaxt = "n", ylab = "")
        if(scen == 1) axis(2)
        if(scen == 1) mtext("Yield", 2, 3.5)
    }
    ## risk
    tmp <- mse$est[,,,"ssb",]
    if(length(dim(tmp)) == 4) dimen <- 2:4 else dimen <-  2:3
    tmp2 <- apply(tmp, dimen, function(x) mean(x < Blim, na.rm = TRUE))
    ylim <- range(tmp2,0, na.rm = TRUE)
    for(scen in 1:nscen){
        if(length(dim(tmp)) == 4) yplot <- tmp2[,scen,] else yplot <- tmp2
        tmp3 <- apply(yplot,1,rbind)
        boxplot(tmp3, ylim = ylim,
                col = cols[1:nms],
                outline = outline,
                xaxt = "n", xlab = "",
                yaxt = "n", ylab = "")
        if(scen == 1) axis(2)
        if(scen == 1) mtext("Risk", 2, 3.5)
    }
    ## yield diff
    tmp <- mse$est[,,,"y",]
    if(length(dim(tmp)) == 4) dimen <- 2:4 else dimen <- 2:3
    tmp2 <- apply(tmp, dimen, diff)
    ylim <- range(tmp2, na.rm = TRUE)
    for(scen in 1:nscen){
        if(length(dim(tmp)) == 4) yplot <- tmp2[,,scen,] else yplot <- tmp2
        tmp3 <- apply(yplot,2,rbind)
        boxplot(tmp3, ylim = ylim,
                col = cols[1:nms],
                outline = outline,
                xaxt = "n", xlab = "",
                yaxt = "n", ylab = "")
        if(scen == 1) axis(2)
        if(scen == 1) mtext("Yield diff", 2, 3.5)
        axis(1, at = seq(mse$pars$HCR), labels = mse$pars$HCR)
    }
    mtext("Management strategies", 1, 3, outer = TRUE)
}


#' @title Plot single time series
#
#' @description Plot single time series
#'
#' @param mse HERE:
#' @param quant HERE:
#' @param col HERE:
#' @param scenario HERE:
#' @param ms HERE:
#' @param ylab HERE:
#' @param yr HERE:
#' @param rep HERE:
#'
#' @export
plotprodmse.ts.single <- function(mse, quant = c("ssbObs","yImp"), col = c("darkorange","goldenrod2"),
                           scenario = 1, ms = 1, ylab = "Biomass '000t",
                           yr = NULL, rep = NULL){

    ##
    nscen <- nrow(scenario)
    dims <- dim(mse$est)
    if(is.null(yr)){
        years <- 1:dims[1]
    }else{
        years <- which(mse$years %in% yr)
    }
    if(is.null(rep)) rep <- 1:dims[5]

    if(length(rep) > 1){
        tmp <- apply(mse$est, c(1:4), mean, na.rm = TRUE)
    }else{
        tmp <- apply(mse$est[,,,,rep, drop=FALSE], c(1:4), mean, na.rm = TRUE)
    }
    ylim <- c(0.9,1.1) * range(tmp[years,ms,scenario,quant], na.rm = TRUE)
    quantL <- plyr::revalue(quant, c("ssbObs" = "SSB perceived",
                                    "yImp" = "Yield after Impl error"))
    ## plot
    opar <- par()
    par(mfrow = c(1,length(scenario)), mar = c(1,1,0,0), oma=c(5,5,4,3))
    for(scen in scenario){
        plot(yr, tmp[years,ms[1],scen,quant[1]],
             ylim = ylim,
             ty="n",
             xaxt = "n", xlab = "",
             yaxt = "n", ylab = "")
        for(q in 1:length(quant)){
            for(m in 1:length(ms)){
                lines(yr, tmp[years, ms[m], scen, quant[q]],
                      col = col[q], lty = m, lwd = 1.5, ty='b', pch=16)
            }
        }
        axis(2)
        axis(1)
        mtext(ylab, 2, 3.5)
        legend("topright", legend = quantL,
               col = col, bg="white",
               lwd=2, lty = 1:length(ms))
    }
    mtext("Projection years", 1, 3, outer = TRUE)

}



#' @title Plot fmsy
#
#' @description Plot fmsy
#'
#' @param mse HERE:
#' @param ms HERE:
#'
#' @export
plotprodmse.fmsy <- function(x, ms = 1, scen = 1,
                             b.unit = 1
                             ){

    tmp <- x[[ms]][[scen]]
    nami <- colnames(tmp)
    esbind <- which(nami == "esb")
    ssbind <- which(nami == "ssb")
    yind <- which(nami == "y")
    tacVarind <- which(nami == "tacVar")

    ##
    cex <- 1.3

    ylim1 <- c(0.0,1.1) * range(tmp[,c(esbind,ssbind,yind)],na.rm=TRUE)
    ylim2 <- c(0.0,1.1) * range(tmp[,c(tacVarind)],na.rm=TRUE)

    par(mfrow=c(1,1), mar = c(5,5,2,5), oma=c(0,0,0,0))
    plot(tmp[,1], tmp[,esbind],
         ty = "n",
         xlab = "", ylab = "",
         xaxt = "n", yaxt = "n",
         ylim = ylim1)
    lines(tmp[,1], tmp[,esbind], ty="l", col="dodgerblue4",pch=16,lwd=2,cex=cex)
    lines(tmp[,1], tmp[,ssbind], ty="l", col="dodgerblue2",pch=16,lwd=2,cex=cex)
    lines(tmp[,1], tmp[,yind], ty="l", col="darkorange",pch=16,lwd=2,cex=cex)
    ## points(tmp[,1], tmp[,esbind], ty="b", col="dodgerblue4",pch=16,lwd=1.5,cex=cex)
    ## points(tmp[,1], tmp[,ssbind], ty="b", col="dodgerblue2",pch=16,lwd=1.5,cex=cex)
    ## points(tmp[,1], tmp[,yind], ty="b", col="darkorange",pch=16,lwd=1.5,cex=cex)
    axis(1)
    axis(2)
    mtext(expression("F ["*yr^{-1}*"]"),1,3)
    if(b.unit == 1){
        addi <- "[t]"
    }else if(b.unit == 1e3){
        addi <- "['000t]"
    }else if(b.unit == 1e6){
        addi <- "[Mt]"
    }
    mtext(paste0("Biomass ",addi), 2, 3)
    par(new=TRUE)
    plot(tmp[,1], tmp[,tacVarind],
         ty = "n",axes = FALSE,
         xlab = "", ylab = "",
         xaxt = "n", yaxt = "n",
         ylim = ylim2)
    lines(tmp[,1], tmp[,tacVarind], ty="l",
          col="darkolivegreen4",pch=16,lwd=2,cex=cex)
    ## points(tmp[,1], tmp[,tacVarind], ty="b", col="darkolivegreen4",pch=16,lwd=1.5,cex=cex)
    axis(4)
    mtext("Interannual catch variability [%]",4,3)
    box()
    legend("topright", legend = c("ESB", "SSB",
                                  "Catch", "Catch variability"),
           col = c("dodgerblue4","dodgerblue2","darkorange","darkolivegreen4"),
           lty=1, lwd=2, bg = "white")
}



#' @title Plot fmsy2
#
#' @description Plot fmsy2
#'
#' @param mse HERE:
#' @param ms HERE:
#'
#' @export
plotprodmse.fmsy2 <- function(x, ms = 1, scen = 1,
                              b.unit = 1, c.unit = 1
                              ){

    tmp <- x[[ms]][[scen]]
    nami <- colnames(tmp)
    ssbind <- which(nami == "ssb")
    ssb5ind <- which(nami == "ssbObsQuant5")
    yind <- which(nami == "y")


    ##
    cex <- 1.3

    ylim1 <- c(0.0,1.1) * range(tmp[,c(ssbind,ssb5ind)] / b.unit,na.rm=TRUE)
    ylim2 <- c(0.0,1.1) * range(tmp[,c(yind)] / c.unit,na.rm=TRUE)

    par(mfrow=c(1,1), mar = c(5,5,2,5), oma=c(0,0,0,0))
    plot(tmp[,1], tmp[,ssbind],
         ty = "n",
         xlab = "", ylab = "",
         xaxt = "n", yaxt = "n",
         ylim = ylim1)
    abline(h = pars$refs[[3]] / b.unit, lwd=1.5, col = "darkred")
    lines(tmp[,1], tmp[,ssbind] / b.unit, ty="l",
          col="dodgerblue2",pch=16,lwd=2,cex=cex)
    lines(tmp[,1], tmp[,ssb5ind] / b.unit, ty="l",
          col="dodgerblue4",pch=16,lwd=2,cex=cex)
    ## points(tmp[,1], tmp[,ssbind] / b.unit, ty="b",
    ##        col="dodgerblue2",pch=16,lwd=1.5,cex=cex)
    ## points(tmp[,1], tmp[,ssb5ind] / b.unit, ty="b",
    ##        col="dodgerblue4",pch=16,lwd=1.5,cex=cex)
    axis(1)
    axis(2)
    mtext(expression("F ["*yr^{-1}*"]"),1,3)
    if(b.unit == 1){
        addi <- "[t]"
    }else if(b.unit == 1e3){
        addi <- "['000t]"
    }else if(b.unit == 1e6){
        addi <- "[Mt]"
    }
    mtext(paste0("Biomass ",addi), 2, 3)
    par(new=TRUE)
    plot(tmp[,1], tmp[,yind],
         ty = "n",axes = FALSE,
         xlab = "", ylab = "",
         xaxt = "n", yaxt = "n",
         ylim = ylim2)
    lines(tmp[,1], tmp[,yind] / c.unit, ty="l",
          col="darkorange",pch=16,lwd=2,cex=cex)
    ## points(tmp[,1], tmp[,yind] / c.unit, ty="b",
    ##        col="darkorange",pch=16,lwd=1.5,cex=cex)
    axis(4)
    if(c.unit == 1){
        addi <- "[t]"
    }else if(c.unit == 1e3){
        addi <- "['000t]"
    }else if(c.unit == 1e6){
        addi <- "[Mt]"
    }
    mtext(paste0("Catch ",addi),4,3)
    box()
    legend("topright", legend = list("SSB",
                                     "5% Perc. of SSB",
                                     expression(B[lim]),
                                     "Catch"),
           col = c("dodgerblue2","dodgerblue4","darkred","darkorange"),
           lty=1, lwd=2, bg = "white")

}



#' @title Plot explo
#
#' @description Plot explo
#'
#' @param mse HERE:
#' @param pars HERE:
#' @param ms HERE:
#'
#' @importFrom RColorBrewer brewer.pal
#'
#' @export
plotprodmse.explo <- function(x, pars, ms = 1, scen = 1,
                              yield.unit = 1, ssb.unit = 1,
                              btrigger.unit = 1, digits = 3){

    x <- x[[ms]][[scen]]
    x <- x[order(x$Fmsy),]
    x <- x[order(x$Btrigger),]
    x$Btrigger <- x$Btrigger / btrigger.unit
    n <- nrow(x)
    uniF <- unique(x$Fmsy)
    nF <- length(uniF)
    uniB <- unique(x$Btrigger)
    nB <- length(uniB)
    fmsy <- pars$refs[[1]]
    btrigger <- pars$refs[[4]] / btrigger.unit

    find <- which.min(abs(uniF-fmsy))
    polyF <- uniF[find] + c(-1,1) * diff(uniF)[c(find,find+1)]/2
    bind <- which.min(abs(uniB-btrigger))
    polyB <- uniB[bind] + c(-1,1) * diff(uniB)[c(bind,bind+1)]/2


    ##
    cexText <- 0.95
    col <- RColorBrewer::brewer.pal(10, "RdYlGn")[-c(1:2,9:10)]
    ncols <- length(col)

    opar <- par(no.readonly = TRUE)
    layout(matrix(1:4,ncol=2,nrow=2,byrow=TRUE),
           widths = c(7,7), heights = c(7,7))
    on.exit(par(opar))
    par(mar = c(3.5,3,3,1), oma = c(3,3,1,1))

    ## yImp
    quant <- "yImp"
    subx <- as.data.frame(cbind(Fmsy = x$Fmsy, Btrigger = x$Btrigger, z = x[quant]))
    tab <- as.matrix(reshape2::dcast(subx, Btrigger ~ Fmsy)[,-1])
    zlim <- range(tab)
    breaks <- seq(zlim[1], zlim[2], length.out = ncols+1)
    ## par(mar = c(6,6,5,0.4))
    plot(x$Btrigger, x$Fmsy, ty='n',
         xaxt="n", yaxt="n",
         xlab="",ylab="")
    image(uniB, uniF, tab, zlim=zlim, breaks=breaks,add = TRUE, col = col)
    abline(v = uniB[-nB] + diff(uniB)/2, col = "grey80", lty=2)
    abline(h = uniF[-nF] + diff(uniF)/2, col = "grey80", lty=2)
    polygon(x = c(polyB,rev(polyB)),
            y = rep(polyF,each=2),
             lwd = 1.3)
    text(x$Btrigger, x$Fmsy, labels=signif(x[[quant]] / yield.unit,digits),
         col = ifelse(x$risk > 0.05, "darkred", "black"), cex = cexText)
    axis(1, cex.axis=1.1)
    axis(2, cex.axis=1.1)
    ## mtext(bquote(B[triggger]),1,3,cex=1.3)
    ## mtext(bquote(F[target]),2,3,cex=1.3)
    if(yield.unit == 1){
        addi <- "[t]"
    }else if(yield.unit == 1e3){
        addi <- "['000t]"
    }else if(yield.unit == 1e6){
        addi <- "[Mt]"
    }
    mtext(paste0("Catch after Implementation error ",addi),
          3, 1, font=2,cex=1.2)
    ## par(mar = c(6,0.2,5,4))
    ## imageScale(z = tab, ylim = zlim, zlim = zlim, col = col, breaks = breaks,
    ##            axis.pos = 4)

    ## risk
    quant <- "risk"
    subx <- as.data.frame(cbind(Fmsy = x$Fmsy, Btrigger = x$Btrigger, z = x[quant]))
    tab <- as.matrix(reshape2::dcast(subx, Btrigger ~ Fmsy)[,-1])
    tab[tab > 0.05] <- 0.051
    zlim <- c(0,0.05)
    breaks <- c(seq(zlim[1], zlim[2], length.out = ncols),1e6)
    ## par(mar = c(6,6,5,0.4))
    plot(x$Btrigger, x$Fmsy, ty='n',
         xaxt="n", yaxt="n",
         xlab="",ylab="")
    image(uniB, uniF, tab, zlim = zlim, breaks = breaks, add = TRUE, col = rev(col))
    abline(v = uniB[-nB] + diff(uniB)/2, col = "grey80", lty=2)
    abline(h = uniF[-nF] + diff(uniF)/2, col = "grey80", lty=2)
    polygon(x = c(polyB,rev(polyB)),
            y = rep(polyF,each=2),
             lwd = 1.3)
    text(x$Btrigger, x$Fmsy, labels=round(x[[quant]],2),
         col = ifelse(x$risk > 0.05, "darkred", "black"), cex = cexText)
    axis(1, cex.axis=1.1)
    axis(2, cex.axis=1.1)
    ## mtext(bquote(B[triggger]),1,3,cex=1.3)
    ## mtext(bquote(F[target]),2,3,cex=1.3)
    mtext("Risk", 3, 1, font=2,cex=1.2)
    ## par(mar = c(6,0.2,5,4))
    ## imageScale(z = tab, ylim = zlim, zlim = zlim, col = rev(col), breaks = breaks,
    ##             axis.pos = 4)


    ## TAC var
    quant <- "tacVar"
    subx <- as.data.frame(cbind(Fmsy = x$Fmsy, Btrigger = x$Btrigger, z = x[quant]))
    tab <- as.matrix(reshape2::dcast(subx, Btrigger ~ Fmsy)[,-1])
    zlim <- range(tab)
    breaks <- seq(zlim[1], zlim[2], length.out = ncols+1)
    ## par(mar = c(6,6,5,0.4))
    plot(x$Btrigger, x$Fmsy, ty='n',
         xaxt="n", yaxt="n",
         xlab="",ylab="")
    image(uniB, uniF, tab, zlim=zlim, breaks=breaks, add = TRUE, col = rev(col))
    abline(v = uniB[-nB] + diff(uniB)/2, col = "grey80", lty=2)
    abline(h = uniF[-nF] + diff(uniF)/2, col = "grey80", lty=2)
    polygon(x = c(polyB,rev(polyB)),
            y = rep(polyF,each=2),
             lwd = 1.3)
    text(x$Btrigger, x$Fmsy, labels=signif(x[[quant]],digits),
         col = ifelse(x$risk > 0.05, "darkred", "black"), cex = cexText)
    axis(1, cex.axis=1.1)
    axis(2, cex.axis=1.1)
    ## mtext(bquote(B[triggger]),1,3,cex=1.3)
    ## mtext(bquote(F[target]),2,3,cex=1.3)
    mtext("Inter-annual catch variability [%]", 3, 1, font=2,cex=1.2)
    ## par(mar = c(6,0.2,5,4))
    ## imageScale(z = tab, ylim = zlim, zlim = zlim, col = rev(col), breaks = breaks,
    ##             axis.pos = 4)

    ## SSB
    quant <- "ssb"
    subx <- as.data.frame(cbind(Fmsy = x$Fmsy, Btrigger = x$Btrigger, z = x[quant]))
    tab <- as.matrix(reshape2::dcast(subx, Btrigger ~ Fmsy)[,-1])
    tab[tab < pars$refs[[3]]] <- pars$refs[[3]]
    tab[tab > pars$refs[[2]]] <- pars$refs[[2]]
    zlim <- c(pars$refs[[3]], pars$refs[[2]])
    breaks <- seq(pars$refs[[3]], pars$refs[[2]], length.out = ncols+1)
    ## par(mar = c(6,6,5,0.4))
    plot(x$Btrigger, x$Fmsy, ty='n',
         xaxt="n", yaxt="n",
         xlab="",ylab="")
    image(uniB, uniF, tab, zlim = zlim, breaks = breaks, add = TRUE, col = col)
    abline(v = uniB[-nB] + diff(uniB)/2, col = "grey80", lty=2)
    abline(h = uniF[-nF] + diff(uniF)/2, col = "grey80", lty=2)
    polygon(x = c(polyB,rev(polyB)),
            y = rep(polyF,each=2),
             lwd = 1.3)
    text(x$Btrigger, x$Fmsy, labels=signif(x[[quant]] / ssb.unit, digits),
         col = ifelse(x$risk > 0.05, "darkred", "black"), cex = cexText)
    axis(1, cex.axis=1.1)
    axis(2, cex.axis=1.1)
    ## mtext(bquote(B[triggger]),1,3,cex=1.3)
    ## mtext(bquote(F[target]),2,3,cex=1.3)
    if(ssb.unit == 1){
        addi <- "[t]"
    }else if(ssb.unit == 1e3){
        addi <- "['000t]"
    }else if(ssb.unit == 1e6){
        addi <- "[Mt]"
    }
    mtext(paste0("SSB ",addi), 3, 1, font=2,cex=1.2)
    ## par(mar = c(6,0.2,5,4))
    ## imageScale(z = tab, ylim = zlim, zlim = zlim, col = col, breaks = breaks,
    ##            axis.pos = 4)

    if(btrigger.unit == 1){
        addi <- "[t]"
    }else if(btrigger.unit == 1e3){
        addi <- "['000t]"
    }else if(btrigger.unit == 1e6){
        addi <- "[Mt]"
    }
    mtext(bquote(B[triggger]~.(addi)),1,1,cex=1.3, outer = TRUE)
    mtext(bquote(F[target] ~ "["*yr^{-1}*"]"),2,0,cex=1.3, outer = TRUE)

}
