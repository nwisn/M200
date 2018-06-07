# Monte Carlo simulation 
# Nicholas Wisniewski

splitAt <- function(x, pos) unname(split(x, cumsum(seq_along(x) %in% (pos+1)))) # splitting function

difference <- function(X, Y, estimator){
    # difference in means
    d <- estimator(X) - estimator(Y)
    return(d)
}

rankt <- function(X, Y, t=difference, estimator=mean){
    # difference in ranks
    XY <- c(X,Y)
    rXY <- rank(XY)
    groups.rank <- splitAt(rXY, length(X))
    rd <- t(groups.rank[[1]], groups.rank[[2]], estimator)
    return(rd)
}

ratio <- function(X, Y, estimator){
    # ratio statistic
    r <- estimator(X) / estimator(Y)
    return(r)
}

formula.tests <- function(X, Y){
    # formula based tests
    Student <- t.test(X, Y, var.equal=T)
    Welch <- t.test(X, Y, var.equal=F)
    Wilcoxon <- wilcox.test(X, Y, correct=F)
    KS <- ks.test(X, Y)
    pvalues <- c(Student$p.value, Welch$p.value, Wilcoxon$p.value, KS$p.value)
    names(pvalues) <- c("Student", "Welch", "Wilcox", "KS")
    return(pvalues)
}

test.stats <- function(X, Y, estimator=mean, t1=difference, t2=rankt) {
    d <- t1(X, Y, estimator) # mean difference
    rd <- t2(X, Y, t1, estimator) # rank difference
    stats <- c(d, rd)
    names(stats) <- c("d", "rankd")
    return(stats)
}

# bigbox; compute pvalues for all tests except twobox
bigbox.all <- function(X, Y, nboot, t=difference, estimator=mean){
    n <- c(length(X), length(Y)) # sample size
    pvalue.formulas <- formula.tests(X, Y) # formula-based tests
    tstats <- test.stats(X, Y, estimator, t, rankt)  # test statistics
    d <- matrix(NA, nrow=nboot, ncol=length(tstats), dimnames=list(1:nboot,c("bigbox","rank")))  # preallocate
    perm <- matrix(NA, nrow=nboot, ncol=length(tstats), dimnames=list(1:nboot,c("permutation","rankperm")))  # preallocate
    
    XY <- c(X,Y) 
    for(i in 1:nboot){
        xy <- sample(XY, size=sum(n), replace=T)
        groups <- splitAt(xy, n[1])
        d[i,] <- test.stats(groups[[1]], groups[[2]], estimator, t, rankt)
        
        xy.perm <- sample(XY, size=sum(n), replace=F)
        groups.perm <- splitAt(xy.perm, n[1])
        perm[i,] <- test.stats(groups.perm[[1]], groups.perm[[2]], estimator, t, rankt)
    }
    #pvalue <- colSums(abs(d) >= abs(tstats)) / nboot
    #pvalue.perm <- colSums(abs(perm) >= abs(tstats)) / nboot
    pvalue <- colSums(sweep(abs(d),2,abs(tstats),"-") >0) / nboot
    pvalue.perm <- colSums(sweep(abs(perm),2,abs(tstats),"-") >0) / nboot
    return(c(pvalue, pvalue.perm, pvalue.formulas))
}

# twobox; compute pvalue for twobox
twobox <- function(X, Y, nboot, t=difference, estimator=mean){
    n <- c(length(X), length(Y)) # sample size
    tstats <- test.stats(X, Y, estimator, t, rankt)[1]  # test statistics
    d <- matrix(NA, nrow=nboot, ncol=length(tstats), dimnames=list(1:nboot,c("twobox")))  # preallocate
    
    X0 <- X - estimator(X)
    Y0 <- Y - estimator(Y)
    for(i in 1:nboot){
        x <- sample(X0, size=n[1], replace=T)
        y <- sample(Y0, size=n[2], replace=T)
        d[i] <- test.stats(x, y, estimator, t)[1]
    }
    #pvalue <- colSums(abs(d) >= abs(tstats)) / nboot
    pvalue <- colSums(sweep(abs(d),2,abs(tstats),"-") >0) / nboot
    return(pvalue)
}

# run both bigbox and twobox, as well as formulas
twogroup <- function(X, Y, nboot, t=difference, estimator=mean){
    one <- bigbox.all(X, Y, nboot, difference, mean)
    two <- twobox(X, Y, nboot, difference, mean)
    return(c(two, one)) #pvalues
}

simulate.base <- function(shift, sigma, xi1, xi2, n1, n2, nsim, nboot){
    require(fGarch)
    testnames <- c("twobox", "bigbox", "rank", "permutation", "rankperm", "Student","Welch", "Wilcox", "KS") 
    pvalues <- matrix(NA, nrow=nsim, ncol=9, dimnames=list(1:nsim, testnames))
    for(i in 1:nsim){
        if(xi1==0) {X <- rnorm(n1, mean=0, sd=1)} else {X <- rsnorm(n1, mean=0, sd=1, xi=xi1)}
        if(xi2==0) {Y <- rnorm(n2, mean=shift, sd=sigma)} else {Y <- rsnorm(n2, mean=shift, sd=sigma, xi=xi2)}
        pvalues[i,] <- twogroup(X, Y, nboot)
    }
    return(pvalues)
}

simulate.curve <- function(shift, sigma, xi1, xi2, samplesize, nsim, nboot) {
    testnames <- c("twobox", "bigbox", "rank", "permutation", "rankperm", "Student","Welch", "Wilcox", "KS") 
    rate <- matrix(NA, nrow=length(shift), ncol=9, dimnames=list(shift, testnames))
    for (i in 1:length(shift)){
        print(paste("     d =",shift[i]))
        sims <- simulate.base(shift[i], sigma, xi1, xi2, samplesize/2, samplesize/2, nsim, nboot)
        rate[i,] <- colSums(sims < 0.05)/nsim
    }
    return(rate)
}

simulate <- function(shift, sigma, xi, samplesize, nsim, nboot){
    mdata <- list()
    for(i in 1:length(samplesize)){
        print(paste("n =",samplesize[i]))
        sims <- simulate.curve(shift, sigma, xi1=xi, xi2=-xi, samplesize=samplesize[i], nsim, nboot)
        sims <- as.data.frame(sims)
        sims$delta <- rownames(sims)
        sims$n <- samplesize[i]
        sims$sigma <- sigma
        sims$xi <- xi
        mdata[[i]] <- melt(sims, id=c("n", "delta","sigma", "xi"))
    }
    merged.data.frame = Reduce(function(...) merge(..., all=T), mdata)
    return(merged.data.frame)
}

FPRplot <- function(data){
    data2 <- subset(data, data$delta==0)
    powerplot <- ggplot(data = data2, aes(x=n, y=value, group=variable)) + geom_line(aes(colour=variable, linetype=variable), size=1) + 
        geom_hline(aes(yintercept=.05), colour="black", linetype="dashed") +
        ggtitle(paste("FPR (sd=", sigma, ", n=",n1+n2, ", nboot=",nboot,", nsim=",nsim, ")", sep="")) +
        geom_point() + 
        ylim(0, max(data2$value)) +
        labs(x="sample size",y="false positive rate")  + theme_linedraw() + 
        theme(legend.key.size = unit(1.5, "cm")) +
        theme(legend.text = element_text(size = 12, colour = "black", angle = 0)) +
        theme(axis.title.y = element_text(size = rel(1.5), angle = 90)) +
        theme(axis.title.x = element_text(size = rel(1.5), angle = 0)) +
        theme(plot.title = element_text(size = rel(1.1), angle = 0)) 
    return(powerplot)
}

powerplot <- function(data, samplesize, sigma, xi){
    powerplot <- ggplot(data = subset(data, data$n==samplesize & data$sigma==sigma & data$xi==xi), aes(x=delta, y=value, group=variable)) + geom_line(aes(colour=variable, linetype=variable), size=1) + 
        geom_hline(aes(yintercept=.8), colour="#990000", linetype="dashed") +
        geom_hline(aes(yintercept=.05), colour="#990000", linetype="dashed") +
        ggtitle(paste("Power Curves\n n=", samplesize, ", s=", sigma, ", xi=",xi, "\nnboot=",nboot,", nsim=",nsim, sep="")) +
        geom_point() + 
        labs(x="difference in means",y="positive rate")  + theme_linedraw() + 
        theme(legend.key.size = unit(1.5, "cm")) +
        theme(legend.text = element_text(size = 12, colour = "black", angle = 0)) +
        theme(axis.title.y = element_text(size = rel(1.5), angle = 90)) +
        theme(axis.title.x = element_text(size = rel(1.5), angle = 0)) +
        theme(plot.title = element_text(size = rel(1.1), angle = 0)) 
    return(powerplot)
}


# --------------------------------------------------------------------------------
# Analysis

#-------------------------------
# Equal variance
require(gridExtra)
nsim <- 400
nboot <- nsim
samplesize <-  c(4,6,8,10,12,14,16,18,20,30,40,60,80,100)
#samplesize <-  c(4,6)
deltas <- seq(0,2,2/8)

#equal variance
result <- simulate(deltas, sigma=1, xi=0, samplesize = samplesize, nsim=nsim, nboot=nboot)
FPRplot(result)
grid.arrange(powerplot(result, 4,1,0), powerplot(result, 6,1,0), 
             powerplot(result, 8,1,0), powerplot(result, 10,1,0),
             powerplot(result, 20,1,0), powerplot(result, 30,1,0),
             powerplot(result, 60,1,0), powerplot(result, 100,1,0),
             ncol=2)

save(result, file="/Users/drwho/Dropbox/Encrypted/M200 Advanced Experimental Statistics/M200.4/twogroup/equalvar.RData")

#unequal variance
nsim <- 400; nboot <- nsim
deltas <- seq(0,5,5/8)
result.unequal <- simulate(deltas, sigma=4, xi=0, samplesize = samplesize, nsim=nsim, nboot=nboot)
save(result.unequal, file="/Users/drwho/Dropbox/Encrypted/M200 Advanced Experimental Statistics/M200.4/twogroup/unequalvar.RData")
result <- result.unequal
FPRplot(result)
grid.arrange(powerplot(result, 4,1,0), powerplot(result, 6,1,0), 
             powerplot(result, 8,1,0), powerplot(result, 10,1,0),
             powerplot(result, 20,1,0), powerplot(result, 30,1,0),
             powerplot(result, 60,1,0), powerplot(result, 100,1,0),
             ncol=2)

#skew
nsim <- 400; nboot <- nsim
deltas <- seq(0,5,5/8)
result.skew <- simulate(deltas, sigma=1, xi=1000, samplesize = samplesize, nsim=nsim, nboot=nboot)
save(result.skew, file="/Users/drwho/Dropbox/Encrypted/M200 Advanced Experimental Statistics/M200.4/twogroup/skew.RData")
result <- result.skew
FPRplot(result)
grid.arrange(powerplot(result, 4,1,0), powerplot(result, 6,1,0), 
             powerplot(result, 8,1,0), powerplot(result, 10,1,0),
             powerplot(result, 20,1,0), powerplot(result, 30,1,0),
             powerplot(result, 60,1,0), powerplot(result, 100,1,0),
             ncol=2)







