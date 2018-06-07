library(knitr)
knit_hooks$set(purl = function(before, options) {
  if (before) return()
  input  = current_input()  # filename of input document
  output = paste(tools::file_path_sans_ext(input), 'R', sep = '.')
  if (knitr:::isFALSE(knitr:::.knitEnv$tangle.start)) {
    assign('tangle.start', TRUE, knitr:::.knitEnv)
    unlink(output)
  }
  cat(options$code, file = output, sep = '\n', append = TRUE)
})
df <- CO2[,c(3,5)] # import data
summary(df)
head(df,10)

# split the dataframe into two groups
X <- df[df$Treatment=="nonchilled","uptake"]  # nonchilled group
Y <- df[df$Treatment=="chilled","uptake"]   # chilled group

library(vioplot)
par(mfrow=c(1,2)) # 2 images

# boxplot with stripchart overlayed
boxplot(uptake ~ Treatment, data=df, main="CO2 uptake", ylab="uptake rate (umol/m^2 sec)")
stripchart(uptake ~ Treatment, method="jitter", add=T, data=df, vertical=T, pch=16, col="red")

# violin plot

vioplot(X,Y, names=c("nonchilled","chilled"), col="red", wex=.7, h=3)

par(mfrow=c(1,2))
qqnorm(X, main="Nonchilled QQ Plot")
qqline(X, col="red")

qqnorm(Y, main="Chilled QQ Plot")
qqline(Y, col="red")
shapiro.test(X)
shapiro.test(Y)
library(car)
bartlett.test(uptake ~ Treatment, data=df)
leveneTest(uptake ~ Treatment, data=df)
fligner.test(uptake ~ Treatment, data=df)
# Student's t-test
ttest <- t.test(X, Y, var.equal = T)
pval <- ttest$p.value          # p-value
# get effect-size and confidence interval
delta <- -diff(ttest$estimate)  # effect size
delta.CI <- ttest$conf.int    # 95% confidence interval 
# format effect size with confidence interval for printing
label.d <- paste("delta = ", round(delta,2), "(", paste(round(delta.CI,2), collapse=","), ")", sep="")
label.p <- paste("pval = ", signif(pval,2), sep="")
noquote(label.d)
noquote(label.p)
wilcox <- wilcox.test(X, Y, conf.int=T, correct=F, exact=F)
pval.U <- wilcox$p.value
dU <- wilcox$estimate     # effect size
dU.CI <- wilcox$conf.int  # confidence interval
# format effect size with confidence interval for printing
label.d <- paste("d= ", round(dU,2), "(", paste(round(dU.CI,2), collapse=","), ")", sep="")
label.p <- paste("p= ", signif(pval.U,2), sep="")
noquote(label.d)
noquote(label.p)
# test statistic: difference in means
dstat <- function(x1, x2, center=mean){
    d <- center(x1) - center(x2) 
    return(d)
}

# test statistic: Student's t
tstat <- function(x1, x2, center=mean, spread=sd){ 
    d <- center(x1) - center(x2) # difference in means
    n1 <- length(x1)
    n2 <- length(x2)
    s <- ((n1-1)*spread(x1)^2 + (n2-1)*spread(x2)^2) / (n1+n2-2) / sqrt((1/n1)+(1/n2))
    t <- d/s
    return(t)
} 

nboot <- 10000  # number of bootstraps
d <- numeric()
d0 <- tstat(X,Y)  # actual test statistic
for(i in 1:nboot){
    x <- sample( c(X,Y), length(X), replace=T)  # resample new x group
    y <- sample( c(X,Y), length(Y), replace=T)  # resample new y group
    d[i] <- tstat(x, y)  # bootstrap test statistic
}
p.onebox <- sum(abs(d)>=abs(d0))/nboot  # compute two-sided p-value
par(mfrow=c(1,2)) # 2 images
hist(d, prob=T, breaks="FD", col="white", main=paste("Null Distribution \n p=",p.onebox), xlab="Mean difference d", xlim=c(min(min(d),-abs(d0)),max(max(d),abs(d0))))  # histogram
abline(v=c(d0,-d0), col="red")
lines(seq(min(d), max(d), length.out = 300) -> y, dnorm(y, mean=mean(d), sd=sd(d)), col = "black", lty="dotted", lwd=5) # add dotted best-fit gaussian

qqnorm(d, main="QQ Plot of Bootstraps"); qqline(d, col="red") # Q-Q plot
# bootstrap using vectorization
boot.1 <- function(X, Y, nboot, center = function(x) sum(x) / length(x)){
    x <- matrix(sample(c(X,Y), size = nboot * length(X), replace = TRUE), nboot, length(X))
    y <- matrix(sample(c(X,Y), size = nboot * length(Y), replace = TRUE), nboot, length(Y))
    x.mean <- apply(x, 1, center)
    y.mean <- apply(y, 1, center)
    d <- x.mean - y.mean             # null distribution
    d0 <- dstat(X,Y, center)         # actual test statistic
    p <- sum(abs(d)>=abs(d0))/nboot  # two-sided p-value
    output <- list(pval = p, statistic = d0, bootstraps = d)
    return(output)
}
tstat <- function(x1, x2, center=mean, spread=sd){
    d <- center(x1) - center(x2) # difference in means
    s <- sqrt( ((spread(x1)^2)/length(x1)) + ((spread(x2)^2)/length(x2)) )  # Welch's correction
    t <- d/s
    return(t)
} 
t0 <- tstat(X,Y)  # actual test statistic
X0 <- X - mean(X)  # recenter X
Y0 <- Y - mean(Y)  # recenter Y
t <- numeric()  # initialize
for(i in 1:nboot){
    x <- sample( X0, length(X), replace=T)  # resample new x group
    y <- sample( Y0, length(Y), replace=T)  # resample new y group
    t[i] <- tstat(x, y)  # bootstrap test statistic
}
p.twobox <- sum(abs(t)>=abs(t0))/nboot  # compute two-sided p-value
par(mfrow=c(1,2)) # 2 images
hist(d, prob=T, breaks="FD", col="white", main=paste("Null Distribution \n p=",p.twobox), xlab="Mean difference d", xlim=c(min(min(d),-abs(d0)),max(max(d),abs(d0))))  # histogram
abline(v=c(d0,-d0), col="red")
lines(seq(min(d), max(d), length.out = 300) -> y, dnorm(y, mean=mean(d), sd=sd(d)), col = "black", lty="dotted", lwd=5) # add dotted best-fit gaussian

qqnorm(d, main="QQ Plot of Bootstraps"); qqline(d, col="red") # Q-Q plot

rstat <- function(XY, nX, center=mean){
    iX <- 1:nX  # indices of X in XY
    rXY <- rank(XY)  # rank all
    rX <- rXY[iX]  # split ranks into X
    rY <- rXY[-iX] # split ranks into Y
    rd <- dstat(rX, rY, center)
    return(rd)
}
XY <- c(X,Y)  # combine data into one-box
d0 <- rstat(XY, length(X))  # actual test statistic
d <- numeric()
for(i in 1:nboot){
    xy <- sample( XY, replace=T)  # resample new x group
    d[i] <- rstat(xy, length(X))  # bootstrap test statistic
}
p.rank <- sum(abs(d)>=abs(d0))/nboot  # compute two-sided p-value
par(mfrow=c(1,2)) # 2 images
hist(d, prob=T, breaks="FD", col="white", main=paste("Null Distribution \n p=",p.rank), xlab="Mean difference d", xlim=c(min(min(d),-abs(d0)),max(max(d),abs(d0))))  # histogram
abline(v=c(d0,-d0), col="red")
lines(seq(min(d), max(d), length.out = 300) -> y, dnorm(y, mean=mean(d), sd=sd(d)), col = "black", lty="dotted", lwd=5) # add dotted best-fit gaussian

qqnorm(d, main="QQ Plot of Bootstraps"); qqline(d, col="red") # Q-Q plot
d0 <- dstat(X,Y)  # actual test statistic
d <- numeric()  # initialize
for(i in 1:nboot){
    x <- sample( X, length(X), replace=T)  # resample new x group
    y <- sample( Y, length(Y), replace=T)  # resample new y group
    d[i] <- dstat(x, y)  # bootstrap test statistic
}
boot.CI <- quantile(d,c(.025,.975))  # find 95% confidence interval

par(mfrow=c(1,2)) # 2 images
hist(d, prob=T, breaks="FD", col="white", main=paste("Bootstrap Distribution \n d=",round(d0,2),",(",paste(round(boot.CI,2),collapse=","),")"), 
     xlab="Mean difference d", xlim=c(min(min(d),-abs(d0)),max(max(d),abs(d0))))  # histogram
abline(v=c(d0,boot.CI), col=c("red","blue","blue"))
lines(seq(min(d), max(d), length.out = 300) -> y, dnorm(y, mean=mean(d), sd=sd(d)), col = "black", lty="dotted", lwd=5) # add dotted best-fit gaussian

qqnorm(d, main="QQ Plot of Bootstraps"); qqline(d, col="red") # Q-Q plot
library(boot)

# test function
diff.means <- function(dat, w) {   
     ix1 <- which(dat[,"Treatment"]=="nonchilled")
     m1 <- sum(dat[ix1,"uptake"] * w[ix1])
     m2 <- sum(dat[-ix1,"uptake"] * w[-ix1])
     m1 - m2
}

result <- boot(df, diff.means, R = nboot, stype = "w", strata = df$Treatment)
bca <- boot.ci(result, type="all")



# function that computes t-like statistics
tstat <- function(X.mu, Y.mu, X.sd, Y.sd, X.n, Y.n){
    d <- X.mu - Y.mu
    s <- sqrt( ((X.sd^2)/X.n) + ((Y.sd^2)/Y.n) )  # Welch's correction
    t <- d/s
    return(t)
}

# two-box bootstrap using vectorization
boot.2 <- function(X, Y, nboot, center = function(x) sum(x) / length(x), spread = sd){
    x <- matrix(sample( X-center(X), size = nboot * length(X), replace = TRUE), nboot, length(X))
    y <- matrix(sample( Y-center(Y), size = nboot * length(Y), replace = TRUE), nboot, length(Y))
    x.mean <- apply(x, 1, center)
    y.mean <- apply(y, 1, center)
    x.sd <- apply(x, 1, spread)
    y.sd <- apply(y, 1, spread)
    t <- tstat(x.mean, y.mean, x.sd, y.sd, length(X), length(Y))  # vectorized
    t0 <- tstat(center(X), center(Y), spread(X), spread(Y), length(X), length(Y))
    p <- sum(abs(t) >= abs(t0))/nboot
    output <- list(pval = p, statistic = t0, bootstraps = t)
    return(output)
}

rstat <- function(XY, nX, center=mean){
    iX <- 1:nX  # indices of X in XY
    rXY <- rank(XY)  # rank all
    rX <- rXY[iX]  # split ranks into X
    rY <- rXY[-iX] # split ranks into Y
    rd <- dstat(rX, rY, center)
    return(rd)
}

XY <- c(X,Y)  # combine data into one-box
r0 <- rstat(XY, length(X))  # actual test statistic
d <- numeric()
for(i in 1:nboot){
    xy <- sample( XY, replace=T)  # resample new x group
    d[i] <- rstat(xy, length(X))  # bootstrap test statistic
}
p.rank <- sum(abs(d)>=abs(r0))/nboot  # compute two-sided p-value

require(DescTools)

onebox.med <- boot.1(X,Y,nboot, center=median)
twobox.med <- boot.2(X,Y,nboot, center=median, spread=function(x) MeanAD(x, FUN=median))

# Bootstrap effect size and confidence interval using vectorization
boot.effect <- function(X, Y, nboot, alpha=0.05, center = function(x) sum(x) / length(x)){
    x <- matrix(sample( X, size = nboot * length(X), replace = TRUE), nboot, length(X))
    y <- matrix(sample( Y, size = nboot * length(Y), replace = TRUE), nboot, length(Y))
    x.mean <- apply(x, 1, center)
    y.mean <- apply(y, 1, center)
    d <- x.mean - y.mean             # null distribution
    d0 <- dstat(X,Y, center)         # actual test statistic
    output <- list(effect = d0, ci = quantile(d, c(alpha/2, 1-(alpha/2))), bootstraps = d)
    return(invisible(output))
}

effect.med <- boot.effect(X,Y,nboot, center=median)

p.comparison <- data.frame(t=round(pval,4),  
                         #LMboot=pval.lmboot[-4], 
                         onebox=p.onebox, 
                         twobox=p.twobox, 
                         rank=round(pval.U[-4],4), 
                         rankboot=p.rank,
                         onebox.med=onebox.med$pval,
                         twobox.med=twobox.med$pval)
rownames(p.comparison) <- c("p-value")
colnames(p.comparison) <- c("Student's t-test", "Bootstrap (one-box)", "Bootstrap (two-box)", "Mann-Whitney U test", "Bootstrap (rank)", "Bootstrap (one-box, median)", "Bootstrap (two-box, median)")
table <- t(p.comparison)
library(knitr)
kable(table, digits=4)
d.comparison <- data.frame(d=c(round(delta,4),round(delta.CI,4)),  
                         #LMboot=pval.lmboot[-4], 
                         u=c(round(dU,4), round(dU.CI,4)),
                         med1=c(effect.med$effect, effect.med$ci),
                         boot=c(d0, boot.CI),
                         bca=c(bca$t0, bca$bca[4:5]),
                         normal = c(bca$t0, bca$normal[2:3]),
                         basic = c(bca$t0, bca$basic[4:5]),
                         pct= c(bca$t0, bca$percent[4:5]))

                         
rownames(d.comparison) <- c("d", "Lower 95% CL", "Upper 95% CL")
colnames(d.comparison) <- c("Student", "Mann-Whitney","Bootstrap (median)", "Bootstrap (naive percentile)", "bca", "normal", "basic", "percentile")
table <- t(d.comparison)
library(knitr)
kable(table, digits=4)
d.onebox <- function(xy, ix, center=mean){
    x.mean <- apply(xy[,ix], 1, center)
    y.mean <- apply(xy[,-ix], 1, center)
    x.mean - y.mean
}

d.twobox <- function(x, y, center=mean){
    x.mean <- apply(x, 1, center)
    y.mean <- apply(y, 1, center)
    x.mean - y.mean
}

rank.onebox <- function(xy, ix, center=mean){
    rxy <- t(apply(xy, 1, rank))
    rx.mean <- apply(rxy[,ix], 1, center)
    ry.mean <- apply(rxy[,-ix], 1, center)
    rx.mean - ry.mean
}

d.effect <- function(x1, x2, center=mean) center(x1) - center(x2) # test statistic

rank.effect <- function(XY, nX, center){
    iX <- 1:nX  # indices of X in XY
    rXY <- rank(XY)  # rank all
    rX <- rXY[iX]  # split ranks into X
    rY <- rXY[-iX] # split ranks into Y
    rd <- d.effect(rX, rY, center)
    return(rd)
}

boot.interval <- function(x, y, alpha=0.05, center=mean){
    x.mean <- apply(x, 1, center)
    y.mean <- apply(y, 1, center)
    d <- x.mean - y.mean
    quantile(d, c(alpha/2, 1-alpha/2))
}

formula.tests <- function(X, Y, alpha=0.05){
    ttest <- t.test(X, Y, var.equal=T, conf.level=(1-alpha)) # Student
    ttest2 <- t.test(X, Y, var.equal=F, conf.level=(1-alpha)) # Welch
    Utest <- wilcox.test(X, Y, exact=F, correct=F, conf.int=T, conf.level=(1-alpha))
    pval <- data.frame(ttest = ttest$p.value, welch = ttest2$p.value, mw = Utest$p.value)
    ci <- data.frame(ttest = ttest$conf.int, welch = ttest2$conf.int, mw = Utest$conf.int)
    effect <- data.frame(ttest=-diff(ttest$estimate), utest=Utest$estimate)
    return(list(effect=effect, ci=ci, p=pval))
}
# the main function
twogroup <- function(X, Y, R=1000, alpha=0.05, center = function(x) sum(x) / length(x)){
    nx <- length(X)
    ny <- length(Y)
    ntot <- nx + ny
    ix <- 1:nx 
    
    # formula-based tests
    results <- formula.tests(X, Y, alpha)  # formula tests
    
    # effect sizes
    d.0 <- d.effect(X, Y, center)  # effect size for mean difference
    d.0.med <- d.effect(X, Y, center=median)  # effect size for median difference
    d.r0 <- rank.effect(c(X,Y), nx, center)  # effect size for rank difference
    
    # one-box statistics
    xy <- matrix(sample(c(X,Y), R*ntot, replace = T), R, ntot)
    d.1 <- d.onebox(xy, ix, center)  # one-box mean null
    d.1.med <- d.onebox(xy, ix, median)  # one-box median null
    d.r <- rank.onebox(xy, ix, center)  # rank null
    
    # two-box sampling
    x <- matrix(sample(X, R*nx, replace = T), R, nx)
    y <- matrix(sample(Y, R*ny, replace = T), R, ny)
    d.2 <- d.twobox(x-center(X), y-center(Y), center)  # two-box null
    d.2.med <- d.twobox(x-median(X), y-median(Y), center=median)  # two-box null
    
    # confidence intervals
    d.ci <- boot.interval(x, y, alpha, center)  # confidence interval on mean difference
    d.ci.med <- boot.interval(x, y, alpha, median)  # confidence interval on median difference
    
    # output tables
    pvalue <- data.frame(Student = results$p$ttest,
                         Welch = results$p$welch,
                         onebox.mean = sum(abs(d.1) >= abs(d.0))/R, 
                         twobox.mean = sum(abs(d.2) >= abs(d.0))/R,
                         onebox.median = sum(abs(d.1.med) >= abs(d.0.med))/R, 
                         twobox.median = sum(abs(d.2.med) >= abs(d.0.med))/R,
                         rankboot = sum(abs(d.r) >= abs(d.r0))/R,
                         Wilcox = results$p$mw)
    ci <- data.frame(d.mean = d.ci, 
                     d.median = d.ci.med,
                     d.ttest = results$ci$ttest,
                     d.welch = results$ci$welch,
                     d.wilcox = results$ci$mw)
    effect <- data.frame(d.mean = d.0, d.median = d.0.med, u = results$effect$utest)
    list(effect=effect, ci=ci, pvalue=pvalue)
}
alltest <- twogroup(X,Y, 10000)
# fast bootstrap
x <- matrix(sample(c(X,Y), size = nboot * length(X), replace = TRUE), nboot, length(X))
y <- matrix(sample(c(X,Y), size = nboot * length(Y), replace = TRUE), nboot, length(Y))
x.mean <- apply(x, 1, mean)
y.mean <- apply(y, 1, mean)
d <- x.mean - y.mean        # null distribution
d0 <- dstat(X,Y)            # actual test statistic
sum(abs(d)>=abs(d0))/nboot  # two-sided p-value
# bootstrap using for-loop implementation
boot.loop <- function(X, Y, nboot){
    d <- numeric()
    d0 <- dstat(X,Y)  # actual test statistic
    for(i in 1:nboot){
        x <- sample( c(X,Y), length(X), replace=T)  # resample new x group
        y <- sample( c(X,Y), length(Y), replace=T)  # resample new y group
        d[i] <- dstat(x, y)  # bootstrap test statistic
    }
    p <- sum(abs(d)>=abs(d0))/nboot  # compute two-sided p-value
    return(list(pval = p, statistic = d0, bootstraps = d))
}

# bootstrap using vectorization
boot.1 <- function(X, Y, nboot, center = function(x) sum(x) / length(x)){
    x <- matrix(sample(c(X,Y), size = nboot * length(X), replace = TRUE), nboot, length(X))
    y <- matrix(sample(c(X,Y), size = nboot * length(Y), replace = TRUE), nboot, length(Y))
    x.mean <- apply(x, 1, center)
    y.mean <- apply(y, 1, center)
    d <- x.mean - y.mean             # null distribution
    d0 <- dstat(X,Y, center)         # actual test statistic
    p <- sum(abs(d)>=abs(d0))/nboot  # two-sided p-value
    output <- list(pval = p, statistic = d0, bootstraps = d)
    return(output)
}

library(microbenchmark)
timings <- microbenchmark(boot.loop(X,Y,100), boot.1(X,Y,100,mean), boot.1(X,Y,100) )
library(knitr)
kable(summary(timings))
# faster mean and standard deviation
mean.fast <- function(x) sum(x) / length(x)
sd.fast <- function(x){
    mu <- mean.fast(x)
    sqrt(sum((x - mu)^2) / (length(x)-1))
}


# function that computes test statistics and effect sizes
tstat <- function(Xbar, Ybar, Xsd, Ysd, Xn, Yn, type="unequal.var"){
    d <- Xbar - Ybar
    if(type == "equal.var"){ # Student's statistic
        s <- sqrt( ((Xn-1)*Xsd^2 + (Yn-1)*Ysd^2) / (Xn + Yn - 2) ) * sqrt((1/Xn) + (1/Yn))
    } else if(type == "unequal.var"){  # Welch's statistic
        s <- sqrt( ((Xsd^2)/Xn) + ((Ysd^2)/Yn) )  
    } else if(type == "Cohen") {  # Cohen's d
        s <- sqrt( ((Xn-1)*Xsd^2 + (Yn-1)*Ysd^2) / (Xn + Yn - 2) ) 
    } else if(type == "delta") {  # simple difference in means
        s <- 1
    }
    t <- d/s
    return(t)
}

# function to bootstrap pvalues. can do any bootstrap test
boot_stats <- function(xy, X, Y, middle=mean.fast, spread=sd.fast, type="unequal.var", rank=F){
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

# use common names for the center (mean, median) and spread (stddev, abs dev)
center <- mean.fast
spread <- sd.fast
nboot <- 10000


# one-box test
xy <- matrix(sample(c(X,Y), size = nboot * length(c(X,Y)), replace = TRUE), nboot, length(c(X,Y))) 
onebox <- boot_stats(xy, X, Y, middle=center, spread=spread, type="equal.var")

# rank test
xy <- matrix(sample(c(X,Y), size = nboot * length(c(X,Y)), replace = TRUE), nboot, length(c(X,Y)))  
xy.rank <- t(apply(xy, 1, rank)) # rank
onebox.rank <- boot_stats(xy.rank, X, Y, middle=mean.fast, type="delta", rank=T)  

# two-box
x <- matrix(sample(X-center(X), size = nboot * length(X), replace = TRUE), nboot, length(X))  # x residuals
y <- matrix(sample(Y-center(Y), size = nboot * length(Y), replace = TRUE), nboot, length(Y))  # y residuals
twobox <- boot_stats(cbind(x,y), X, Y, middle=center, spread=spread, type="unequal.var")

# effect size with confidence interval
x <- matrix(sample(X, size = nboot * length(X), replace = TRUE), nboot, length(X))  # x
y <- matrix(sample(Y, size = nboot * length(Y), replace = TRUE), nboot, length(Y))  # y
effectsize <- boot_stats(cbind(x,y), X, Y, middle=center, type="delta")
effectsize.d <- boot_stats(cbind(x,y), X, Y, middle=center, type="Cohen")

# two-box
x <- matrix(sample(X-center(X), size = nboot * length(X), replace = TRUE), nboot, length(X))  # residuals
y <- matrix(sample(Y-center(Y), size = nboot * length(Y), replace = TRUE), nboot, length(Y))  # residuals
    
# dither noise for small samples
noise <- matrix(rnorm(nboot*length(c(X,Y)), 0, min(sd(X),sd(Y))/length(c(X,Y))^(3/2)), nboot, length(c(X,Y)))

# compute p-values
twobox <- boot_stats(cbind(x,y) + noise, X, Y, middle=center, spread=spread, type="unequal.var")
one <- boot.1(X,Y,nboot)                               # one-box NHST on means
one.med <- boot.1(X,Y,nboot, center=median)            # one-box NHST on medians
one.rank <- boot.1(X,Y,nboot, rank=T)                  # one-box NHST on ranks
two <- boot.2(X,Y,nboot)                               # two-box NHST on means
two.med <- boot.2(X,Y,nboot, center=median)            # two-box NHST on medians
effect <- boot.effect(X,Y,nboot)                       # mean difference effect size and confidence interval
effect.med <- boot.effect(X,Y,nboot, center=median)    # median difference effect size and confidence interval
one$pval
one.med$pval
one.rank$pval
two$pval
two.med$pval
