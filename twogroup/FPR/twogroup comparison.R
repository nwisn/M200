

runalltests <- function(X, Y, nboot=300){
    
    require(DescTools)
    
    # function that takes the matrix of bootstraps with the original data, and computes p-values for different statistics
    # returns a bunch of good stuff, including all the bootstrap tests statistics
    boot_stats <- function(xy, X, Y, middle, spread, type="equal.var", rank=F){
        nx <- length(X)
        ny <- length(Y)
        x.mean <- apply(xy[, (1:nx)], 1, middle)
        y.mean <- apply(xy[,-(1:nx)], 1, middle)
        x.sd <- apply(xy[, (1:nx)], 1, spread)
        y.sd <- apply(xy[,-(1:nx)], 1, spread)
        t <- tstat(x.mean, y.mean, x.sd, y.sd, nx, ny, type=type)
        if(rank==T){xyr <- rank(c(X,Y)); X <- xyr[1:nx];  Y <- xyr[-(1:ny)] }
        t0 <- tstat(middle(X), middle(Y), spread(X), spread(Y), nx, ny, type=type)
        p <- sum(abs(t) >= abs(t0))/dim(xy)[1]
        return( list(p = p, t0 = t0, t = t, X.mu = x.mean, Y.mu = y.mean, X.s = x.sd, Y.s = y.sd))
    }
    
    # function that computes t-like statistics
    tstat <- function(Xbar, Ybar, Xsd, Ysd, Xn, Yn, type="equal.var"){
        d <- Xbar - Ybar
        if(type == "equal.var"){
            s <- sqrt( ((Xn-1)*Xsd^2 + (Yn-1)*Ysd^2) / (Xn + Yn - 2) ) * sqrt((1/Xn) + (1/Yn))
        } else if(type == "unequal.var"){
            s <- sqrt( ((Xsd^2)/Xn) + ((Ysd^2)/Yn) )  # Welch's correction
        } else if(type == "Cohen") {
            s <- sqrt( ((Xn-1)*Xsd^2 + (Yn-1)*Ysd^2) / (Xn + Yn - 2) ) 
        } else if(type == "delta") {
            s <- 1
        }
        t <- d/s
        return(t)
    }
    
    
    
    # Get the sample sizes
    nX=length(X)
    nY=length(Y)
    
    # dither noise for medians
    noise <- matrix(rnorm(nboot*(nX+nY),0, min(sd(X),sd(Y))/(nX+nY)^(3/2)), nboot, nX+nY)
    
    # --------------------------------------------------------------------------------------------------------------
    # Bigbox tests
    # --------------------------------------------------------------------------------------------------------------
    # bigbox 
    xy <- matrix(sample(c(X,Y), size = nboot * length(c(X,Y)), replace = TRUE), nboot, length(c(X,Y)))  # bigbox
    
    # for means
    onebox.mean <- list()
    onebox.mean[[1]] <- boot_stats(xy, X, Y, middle=mean, spread=sd, type="delta")  # effect size
    onebox.mean[[2]] <- boot_stats(xy + noise, X, Y, middle=mean, spread=sd, type="equal.var")  # standard deviation
    onebox.mean[[3]] <- boot_stats(xy + noise, X, Y, middle=mean, spread=function(x) MeanAD(x, FUN=mean), type="equal.var")  # mean absolute deviation about the mean
    onebox.mean[[4]] <- boot_stats(xy + noise, X, Y, middle=mean, spread=function(x) MeanAD(x, FUN=median), type="equal.var")  # mean absolute deviation about the median
    onebox.mean[[5]] <- boot_stats(xy + noise, X, Y, middle=mean, spread=function(x) mad(x, center=mean(x)), type="equal.var" ) # median absolute deviation about the mean
    onebox.mean[[6]] <- boot_stats(xy + noise, X, Y, middle=mean, spread=function(x) mad(x, center=median(x)), type="equal.var" ) # median absolute deviation about the median
    
    
    # for medians
    onebox.median <- list()
    onebox.median[[1]] <- boot_stats(xy + noise, X, Y, middle=median, spread=sd, type="delta") # effect size
    onebox.median[[2]] <- boot_stats(xy + noise, X, Y, middle=median, spread=sd, type="equal.var") # standard deviation
    onebox.median[[3]] <- boot_stats(xy + noise, X, Y, middle=median, spread=function(x) MeanAD(x, FUN=mean), type="equal.var")  # mean absolute deviation about the mean
    onebox.median[[4]] <- boot_stats(xy + noise, X, Y, middle=median, spread=function(x) MeanAD(x, FUN=median), type="equal.var")  # mean absolute deviation about the median
    onebox.median[[5]] <- boot_stats(xy + noise, X, Y, middle=median, spread=function(x) mad(x, center=mean(x)), type="equal.var" ) # median absolute deviation about the mean
    onebox.median[[6]] <- boot_stats(xy + noise, X, Y, middle=median, spread=function(x) mad(x, center=median(x)), type="equal.var" ) # median absolute deviation about the median
    
    
    # --------------------------------------------------------------------------------------------------------------
    # bigbox residuals
    xy.res <- matrix(sample(c(X-center(X), Y-center(Y)), size = nboot * length(c(X,Y)), replace = TRUE), nboot, length(c(X,Y)))  # bigbox residuals
    
    # for means
    onebox.mean.res <- list()
    onebox.mean.res[[1]] <- boot_stats(xy.res, X, Y, middle=mean, spread=sd, type="delta")  # effect size
    onebox.mean.res[[2]] <- boot_stats(xy.res + noise, X, Y, middle=mean, spread=sd, type="equal.var")  # standard deviation
    onebox.mean.res[[3]] <- boot_stats(xy.res + noise, X, Y, middle=mean, spread=function(x) MeanAD(x, FUN=mean), type="equal.var")  # mean absolute deviation about the mean
    onebox.mean.res[[4]] <- boot_stats(xy.res + noise, X, Y, middle=mean, spread=function(x) MeanAD(x, FUN=median), type="equal.var")  # mean absolute deviation about the median
    onebox.mean.res[[5]] <- boot_stats(xy.res + noise, X, Y, middle=mean, spread=function(x) mad(x, center=mean(x)), type="equal.var" ) # median absolute deviation about the mean
    onebox.mean.res[[6]] <- boot_stats(xy.res + noise, X, Y, middle=mean, spread=function(x) mad(x, center=median(x)), type="equal.var" ) # median absolute deviation about the median
    
    
    # for medians
    onebox.median.res <- list()
    onebox.median.res[[1]] <- boot_stats(xy.res + noise, X, Y, middle=median, spread=sd, type="delta") # effect size
    onebox.median.res[[2]] <- boot_stats(xy.res + noise, X, Y, middle=median, spread=sd, type="equal.var") # standard deviation
    onebox.median.res[[3]] <- boot_stats(xy.res + noise, X, Y, middle=median, spread=function(x) MeanAD(x, FUN=mean), type="equal.var")  # mean absolute deviation about the mean
    onebox.median.res[[4]] <- boot_stats(xy.res + noise, X, Y, middle=median, spread=function(x) MeanAD(x, FUN=median), type="equal.var")  # mean absolute deviation about the median
    onebox.median.res[[5]] <- boot_stats(xy.res + noise, X, Y, middle=median, spread=function(x) mad(x, center=mean(x)), type="equal.var" ) # median absolute deviation about the mean
    onebox.median.res[[6]] <- boot_stats(xy.res + noise, X, Y, middle=median, spread=function(x) mad(x, center=median(x)), type="equal.var" ) # median absolute deviation about the median
    
    

    # --------------------------------------------------------------------------------------------------------------
    # bigbox rank
    
    xy.rank <- t(apply(xy, 1, rank)) # bigbox rank
    
    # for means
    onebox.rank <- list()
    onebox.rank[[1]] <- boot_stats(xy.rank, X, Y, middle=mean, spread=sd, type="delta", rank=T)  # effect size
    onebox.rank[[2]] <- boot_stats(xy.rank + noise, X, Y, middle=mean, spread=sd, type="equal.var", rank=T)  # standardized rank difference
    
    
    
    # --------------------------------------------------------------------------------------------------------------
    # twobox tests
    # --------------------------------------------------------------------------------------------------------------
    x <- matrix(sample(X-center(X), size = nboot * nX, replace = TRUE), nboot, nX)  # bigbox residuals
    y <- matrix(sample(Y-center(Y), size = nboot * nY, replace = TRUE), nboot, nY)  # bigbox residuals
    xy.two <- cbind(x,y)
    
    
    # for means
    twobox.mean <- list()
    twobox.mean[[1]] <- boot_stats(xy.two, X, Y, middle=mean, spread=sd, type="delta")  # effect size
    twobox.mean[[2]] <- boot_stats(xy.two + noise, X, Y, middle=mean, spread=sd, type="unequal.var")  # standard deviation
    twobox.mean[[3]] <- boot_stats(xy.two + noise, X, Y, middle=mean, spread=function(x) MeanAD(x, FUN=mean), type="unequal.var")  # mean absolute deviation about the mean
    twobox.mean[[4]] <- boot_stats(xy.two + noise, X, Y, middle=mean, spread=function(x) MeanAD(x, FUN=median), type="unequal.var")  # mean absolute deviation about the median
    twobox.mean[[5]] <- boot_stats(xy.two + noise, X, Y, middle=mean, spread=function(x) mad(x, center=mean(x)), type="unequal.var" ) # median absolute deviation about the mean
    twobox.mean[[6]] <- boot_stats(xy.two + noise, X, Y, middle=mean, spread=function(x) mad(x, center=median(x)), type="unequal.var" ) # median absolute deviation about the median
    
    
    # for medians
    twobox.median <- list()
    twobox.median[[1]] <- boot_stats(xy.two + noise, X, Y, middle=median, spread=sd, type="delta") # effect size
    twobox.median[[2]] <- boot_stats(xy.two + noise, X, Y, middle=median, spread=sd, type="unequal.var") # standard deviation
    twobox.median[[3]] <- boot_stats(xy.two + noise, X, Y, middle=median, spread=function(x) MeanAD(x, FUN=mean), type="unequal.var")  # mean absolute deviation about the mean
    twobox.median[[4]] <- boot_stats(xy.two + noise, X, Y, middle=median, spread=function(x) MeanAD(x, FUN=median), type="unequal.var")  # mean absolute deviation about the median
    twobox.median[[5]] <- boot_stats(xy.two + noise, X, Y, middle=median, spread=function(x) mad(x, center=mean(x)), type="unequal.var" ) # median absolute deviation about the mean
    twobox.median[[6]] <- boot_stats(xy.two + noise, X, Y, middle=median, spread=function(x) mad(x, center=median(x)), type="unequal.var" ) # median absolute deviation about the median
    
    
    # --------------------------------------------------------------------------------------------------------------
    # formula tests
    # --------------------------------------------------------------------------------------------------------------
    ttest <- t.test(X,Y, var.equal=T)
    welchtest <- t.test(X,Y, var.equal=F)
    utest <- wilcox.test(X,Y)
    
    # --------------------------------------------------------------------------------------------------------------
    # Collect p-values
    # --------------------------------------------------------------------------------------------------------------
    
    pval <- numeric()
    nlist <- length(onebox.mean)
    for (i in 1:nlist){
        pval[i] <- onebox.mean[[i]][[1]]
        pval[i+nlist] <- onebox.median[[i]][[1]]
        pval[i+2*nlist] <- onebox.mean.res[[i]][[1]]
        pval[i+3*nlist] <- onebox.median.res[[i]][[1]]
        pval[i+4*nlist] <- twobox.mean[[i]][[1]]
        pval[i+5*nlist] <- twobox.median[[i]][[1]]
    }
    
    current <- length(pval)
    for (i in 1:length(onebox.rank)){
        pval[i+current] <- onebox.rank[[i]][[1]]
    }
    
    current <- length(pval)
    pval[1+current] <- ttest$p.value
    pval[2+current] <- welchtest$p.value
    pval[3+current] <- utest$p.value
    
    testnames <- c(
        "onebox mean difference",
        "onebox mean t (stddev)",
        "onebox mean t (mean AD mean)",
        "onebox mean t (mean AD median)",
        "onebox mean t (median AD mean)",
        "onebox mean t (median AD median)",
        
        "onebox median difference",
        "onebox median t (stddev)",
        "onebox median t (mean AD mean)",
        "onebox median t (mean AD median)",
        "onebox median t (median AD mean)",
        "onebox median t (median AD median)",
        
        "onebox residuals mean difference",
        "onebox residuals mean t (stddev)",
        "onebox residuals mean t (mean AD mean)",
        "onebox residuals mean t (mean AD median)",
        "onebox residuals mean t (median AD mean)",
        "onebox residuals mean t (median AD median)",
        
        "onebox residuals median difference",
        "onebox residuals median t (stddev)",
        "onebox residuals median t (mean AD mean)",
        "onebox residuals median t (mean AD median)",
        "onebox residuals median t (median AD mean)",
        "onebox residuals median t (median AD median)",
        
        "twobox residuals mean difference",
        "twobox residuals mean t (stddev)",
        "twobox residuals mean t (mean AD mean)",
        "twobox residuals mean t (mean AD median)",
        "twobox residuals mean t (median AD mean)",
        "twobox residuals mean t (median AD median)",
        
        "twobox residuals median difference",
        "twobox residuals median t (stddev)",
        "twobox residuals median t (mean AD mean)",
        "twobox residuals median t (mean AD median)",
        "twobox residuals median t (median AD mean)",
        "twobox residuals median t (median AD median)",
        
        "onebox rank difference",
        "onebox rank t (stddev)",
        
        "Student's t-test",
        "Welch's t-test",
        "Wilcoxon rank-sum test"
    )
    
    names(pval) <- testnames
    pcolumn <- cbind(pval)
    
    return(pcolumn)
}


# function to run the sims at a single lattice point
require(fGarch)
runsims <- function(nsim, nboot, n1, n2, delta, sigma, skew=0){
    pvals <- matrix(NA, nrow=41, ncol=nsim)
    for (i in 1:nsim){
        print(paste("running sim",i, "of", nsim))
        if (skew==0){
            X <- rnorm(n1, 0, 1)
            Y <- rnorm(n2, 0+delta, sigma)
        } else {
            X <- rsnorm(n1, 0, 1, -skew)
            Y <- rsnorm(n2, 0+delta, sigma, skew)
        }
        results <- runalltests(X,Y, nboot)
        pvals[,i] <- results
        rownames(pvals) <- rownames(results)
    }
    FPR <- rowSums(pvals < 0.05, na.rm=T)/nsim
    return( list(pvals = pvals, fpr = FPR))
}







# confidence interval on .05 false positive rate
f.l <- function(x) 1.96 +(x-.05)/sqrt((.05*(1-.05))/nsim)
f.u <- function(x) -1.96 +(x-.05)/sqrt((.05*(1-.05))/nsim)
uniroot(f, interval=c(0, 1))$root
fpr.l <- uniroot(f.l, interval=c(0, 1))$root
fpr.u <- uniroot(f.u, interval=c(0, 1))$root

# order rows into better blocks
roworder <- c(  39, # Student
                seq(1,6,1), # onebox mean
                seq(1+6*1, 6+6*1, 1), # onebox median
                c(41, 37,38), # onebox rank, Wilcox
                seq(1+6*2, 6+6*2, 1), # onebox residuals mean
                seq(1+6*3, 6+6*3, 1), # onebox residuals median
                40, # Welch
                seq(1+6*4, 6+6*4, 1), # twobox mean
                seq(1+6*5, 6+6*5, 1) # twobox median
)








# Run the tests

n1 <- c(2, 3, 4, 5, 6, 7, 8, 9, 10)*3
n2 <- c(2, 3, 4, 5, 6, 7, 8, 9, 10)

#n1 <- c(2,3)
#n2 <- c(2,3)
delta <- seq(0,2,.25)
sigma <- c(1, 6)


nsim <- 300; nboot <- 300
FPR <- matrix(NA, nrow=41, ncol=length(n1))
simresult <- list()
for (n.i in 1:length(n1)){
    print(paste("*** evaluating at sample size", n1[n.i], ", ",n2[n.i] ))
    simresult[[n.i]] <- runsims(nsim, nboot, n1[n.i], n2[n.i], delta[1], sigma[2], skew=0)
    FPR[,n.i] <- simresult[[n.i]]$fpr
    rownames(FPR) <- names(simresult[[1]]$fpr)
    colnames(FPR) <- n1
} 




# heatmap

boxlabels <- FPR 
boxlabels[FPR <fpr.u & FPR>fpr.l] <- NA

require(gplots)
#my_palette <- colorRampPalette(c("white", "yellow", "red"))(n = 299)
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)
col_breaks = c(seq(0,fpr.l, length=100),  # for red
               seq(fpr.l,fpr.u, length=100),              # for yellow
               seq(fpr.u, max(FPR), length=100))              # for green


#lmat <- rbind( c(1,3,4), c(2,1,4) )
lhei <- c(1,6)
lwid <- c(.1, 5, 0)
heatmap.2(FPR[roworder,], Rowv=F, Colv=F, col=my_palette, breaks=col_breaks, trace="none", 
          lwid=lwid, 
          #lmat=lmat, 
          lhei=lhei,
          cellnote=round(boxlabels[roworder,],2), notecol="black", cexRow=.7, cexCol=.7, notecex=.4, main="Unequal n, Unequal Variance", xlab="Group Sample Size",
          rowsep =c(7,7+6, 7+6+3, 7+6+3+6, 7+6+3+6+6, 7+6+3+6+6+7, 7+6+3+6+6+7+6),
          sepcolor="black",
          sepwidth=c(.1,.1),
          margins = c(5, 20),
          RowSideColors = c(rep("red",7), rep("orange",6), rep("yellow",3), rep("green",6), rep("blue",6), rep("darkblue",7), rep("violet",6))
)


#par(lend = 1)           # square line ends for the color legend
#legend("topleft",      # location of the legend on the heatmap plot
#       legend = c("Onebox Mean", "Onebox Median", "Onebox Rank", "Onebox Mean Residuals", "Onebox Median Residuals", "Twobox Mean", "Twobox Median"), # category labels
#       col = c("red", "orange", "yellow", "green", "blue", "darkblue", "violet"),  # color key
#       lty= 1,             # line style
#       lwd =5,            # line width
#       cex=.4
#)


save(simresult, FPR, file="~/Desktop/n1 gt n2 unequalvar.RData")
#load(file="~/Desktop/Twogroup Simulations (May 14)/equalvar.RData")



