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
df <- PlantGrowth
df$group <- factor(df$group, labels = c("Control", "Treatment 1", "Treatment 2"))
df
#df.rm <- read.csv("/Users/drwho/Dropbox/Encrypted/M200 Advanced Experimental Statistics/M200.4/repeated measures dataset.csv")
#df.aov <- read.csv("/Users/drwho/Dropbox/Encrypted/M200 Advanced Experimental Statistics/M200.4/anova dataset.csv")
#df <- df.aov
#df
boxplot(weight ~ group, data=df) 
stripchart(weight ~ group, add=T, method='jitter', vertical=T, pch=16, col="red", data = df)
abline(h=mean(df$weight), lty='dotted')  # draw dotted line at the grand mean
model <- lm(weight ~ group, data=df)
result <- anova(model)
result
p.standard <- result$`Pr(>F)`[1]
# get the F-statistic from lm()
F0 = result$`F value`[1]

# Bootstrap the F-statistic to compute the p-value
deals <- 10000
bootF <- rep(NA,deals)
for(i in 1:deals){
    df$boot <- sample(df$weight, replace=T)
    bootresult <- anova(lm(boot ~ group, data=df))
    bootF[i] <- bootresult$`F value`[1]
}
pval <- sum(bootF > F0)/deals
hist(bootF, breaks='FD', main="Bootstrapped null distribution")
abline(v=F0, col='red')
text(x=F0, y=100, labels=paste("p=",pval,sep=""))
dfB <- result$Df[1]
dfW <- result$Df[2]
qqplot( qf( seq(0,1,length=deals), dfB, dfW), bootF, 
        main="QQ plot for F-distribution", 
        xlab="Theoretical F Quantiles", 
        ylab="Bootstrap F-like Quantiles")
abline(0, 1, col='red')
pairwise.t.test(df$weight, df$group, p.adjust="bonferroni")
bartlett.test(weight ~ group, data=df)
fligner.test(weight ~ group, data=df)
library(car)
leveneTest(weight ~ group, data=df)
par(mfrow=c(2,2))
plot(model)
kruskal.test(weight~group, data=df)
kw <- kruskal.test(weight~group, data=df)
p.kw <- kw$p.value
pairwise.wilcox.test(df$weight, df$group, p.adjust.method="bonferroni", exact=F, correct=F)
# Define functions to compute the Sum of Squares, and Sum of Abs
ssq <- function(x) sum(x^2)
sab <- function(x) sum(abs(x))

SSstat <- function(Y, group, center=mean, agg=ssq){
    groupMeans = ave( Y, group, FUN = center)
    SSB = agg(groupMeans - center(Y))
    return(SSB)
}

# A more general function to compute the F-statistic
Fstat <- function(Y, group, center=mean, agg=ssq) {   
    groupMeans = ave( Y, group, FUN = center)
    SSB = agg(groupMeans - center(Y))
    SSW = agg(Y - groupMeans)
    return( SSB / SSW)
}
# Run it on the dataframe to get the test statistic:
SS0 <- SSstat(df$weight,df$group, center=mean, agg=ssq)
F0 <- Fstat(df$weight,df$group, center=mean, agg=ssq)

# Bootstrap the SS-statistic to compute the p-value
bootSS <- rep(NA,deals)
for(i in 1:deals){
    df$boot <- sample(df$weight, replace=T)
    bootSS[i] <- SSstat(df$boot, df$group, center=mean, agg=ssq)
}
pval.SS <- sum(bootSS > SS0)/deals

# Bootstrap the F-statistic to compute the p-value
bootF <- rep(NA,deals)
for(i in 1:deals){
    df$boot <- sample(df$weight, replace=T)
    bootF[i] <- Fstat(df$boot, df$group, center=mean, agg=ssq)
}
pval.F <- sum(bootF > F0)/deals
hist(bootSS, breaks='FD', main="Bootstrapped null distribution")
abline(v=SS0, col='red')
text(x=SS0, y=100, labels=paste("p=",pval.SS,sep=""))
hist(bootF, breaks='FD', main="Bootstrapped null distribution")
abline(v=F0, col='red')
text(x=F0, y=100, labels=paste("p=",pval,sep=""))
#dfB <- nlevels(df$group) - 1
#dfW <- nrow(df) - nlevels(df$group)
qqplot( qf( seq(0,1,length=deals), dfB, dfW), bootF, 
        main="QQ plot for F-distribution", 
        xlab="Theoretical F Quantiles", 
        ylab="Bootstrap F-like Quantiles")
abline(0, dfB/dfW, col='red')
# Run it on the dataframe to get the test statistic:
SS0 <- SSstat(df$weight,df$group, center=mean, agg=sab)
F0 <- Fstat(df$weight,df$group, center=mean, agg=sab)

# Bootstrap the S-abs statistic to compute the p-value
bootSS <- rep(NA,deals)
for(i in 1:deals){
    df$boot <- sample(df$weight, replace=T)
    bootSS[i] <- SSstat(df$boot, df$group, center=mean, agg=sab)
}
pval.SS <- sum(bootSS > SS0)/deals

# Bootstrap the F-abs statistic to compute the p-value
bootF <- rep(NA,deals)
for(i in 1:deals){
    df$boot <- sample(df$weight, replace=T)
    bootF[i] <- Fstat(df$boot, df$group, center=mean, agg=sab)
}
pval.F <- sum(bootF > F0)/deals
hist(bootSS, breaks='FD', main="Bootstrapped null distribution")
abline(v=SS0, col='red')
text(x=SS0, y=100, labels=paste("p=",pval.SS,sep=""))
hist(bootF, breaks='FD', main="Bootstrapped null distribution")
abline(v=F0, col='red')
text(x=F0, y=100, labels=paste("p=",pval.F,sep=""))
qqplot( qchisq( seq(0,1,length=deals), 8), bootF, 
        main="QQ plot for F-abs distribution", 
        xlab="Theoretical Chi-Squared Quantiles", 
        ylab="Bootstrap F-abs Quantiles")
abline(0,1/dfW, col='red')
# Run it on the dataframe to get the test statistic:
F0 = Fstat(rank(df$weight),df$group, center=mean, agg=ssq)

# Bootstrap the F-statistic to compute the p-value
bootF <- rep(NA,deals)
for(i in 1:deals){
    df$boot <- sample(df$weight, replace=T)
    bootF[i] <- Fstat(rank(df$boot), df$group, center=mean, agg=ssq)
}
pval <- sum(bootF > F0)/deals
hist(bootF, breaks='FD', main="Bootstrapped null distribution")
abline(v=F0, col='red')
text(x=F0, y=100, labels=paste("p=",pval,sep=""))
qqplot( qf( seq(0,1,length=deals), dfB, dfW), bootF, 
        main="QQ plot for F-distribution", 
        xlab="Theoretical F Quantiles", 
        ylab="Bootstrap F-like Quantiles")
# Run it on the dataframe to get the test statistic:
SS0 <- SSstat(df$weight,df$group, center=mean, agg=ssq)
F0 <- Fstat(df$weight,df$group, center=mean, agg=ssq)

# Bootstrap the SS-statistic to compute the p-value
bootSS <- rep(NA,deals)
residuals <- df$weight - ave(df$weight, df$group)
for(i in 1:deals){
    df$boot <- ave(residuals, df$group, FUN=function(x) sample(x, replace=TRUE))
    bootSS[i] <- SSstat(df$boot, df$group, center=mean, agg=ssq)
}
pval.SS <- sum(bootSS > SS0)/deals

# Bootstrap the F-statistic to compute the p-value
bootF <- rep(NA,deals)
residuals <- df$weight - ave(df$weight, df$group)
for(i in 1:deals){
    df$boot <- ave(residuals, df$group, FUN=function(x) sample(x, replace=TRUE))
    bootF[i] <- Fstat(df$boot, df$group, center=mean, agg=ssq)
}
pval.F <- sum(bootF > F0)/deals
hist(bootSS, breaks='FD', main="Bootstrapped null distribution")
abline(v=SS0, col='red')
text(x=SS0, y=100, labels=paste("p=",pval.SS,sep=""))
hist(bootF, breaks='FD', main="Bootstrapped null distribution")
abline(v=F0, col='red')
text(x=F0, y=100, labels=paste("p=",pval,sep=""))
#dfB <- nlevels(df$group) - 1
#dfW <- nrow(df) - nlevels(df$group)
qqplot( qf( seq(0,1,length=deals), dfB, dfW), bootF, 
        main="QQ plot for F-distribution", 
        xlab="Theoretical F Quantiles", 
        ylab="Bootstrap F-like Quantiles")
abline(0, dfB/dfW, col='red')
# Run it on the dataframe to get the test statistic:
SS0 <- SSstat(df$weight,df$group, center=mean, agg=sab)
F0 <- Fstat(df$weight,df$group, center=mean, agg=sab)

# Bootstrap the S-abs statistic to compute the p-value
bootSS <- rep(NA,deals)
for(i in 1:deals){
    df$boot <- ave(residuals, df$group, FUN=function(x) sample(x, replace=TRUE))
    bootSS[i] <- SSstat(df$boot, df$group, center=mean, agg=sab)
}
pval.SS <- sum(bootSS > SS0)/deals

# Bootstrap the F-abs statistic to compute the p-value
bootF <- rep(NA,deals)
for(i in 1:deals){
    df$boot <- ave(residuals, df$group, FUN=function(x) sample(x, replace=TRUE))
    bootF[i] <- Fstat(df$boot, df$group, center=mean, agg=sab)
}
pval.F <- sum(bootF > F0)/deals
hist(bootSS, breaks='FD', main="Bootstrapped null distribution")
abline(v=SS0, col='red')
text(x=SS0, y=100, labels=paste("p=",pval.SS,sep=""))
hist(bootF, breaks='FD', main="Bootstrapped null distribution")
abline(v=F0, col='red')
text(x=F0, y=100, labels=paste("p=",pval.F,sep=""))
qqplot( qchisq( seq(0,1,length=deals), 8), bootF, 
        main="QQ plot for F-abs distribution", 
        xlab="Theoretical Chi-Squared Quantiles", 
        ylab="Bootstrap F-abs Quantiles")
abline(0,1/dfW, col='red')
# faster mean and standard deviation
mean.fast <- function(x) sum(x) / length(x)
sd.fast <- function(x){
    mu <- mean.fast(x)
    sqrt(sum((x - mu)^2) / (length(x)-1))
}

# Define functions to compute the Sum of Squares, and Sum of Abs
ssq <- function(x) sum(x^2)
sab <- function(x) sum(abs(x))

# A more general function to compute the F-statistic
Fstat <- function(Y, group, center=mean.fast, agg=ssq) {   
    groupMeans = ave( Y, group, FUN = center)
    SSB = agg(groupMeans - center(Y))
    SSW = agg(Y - groupMeans)
    return( SSB / SSW)
}

Y <- df$weight
G <- df$group
center <- mean.fast
nboot <- 10000
# one-box
y <- matrix(sample(Y, size = nboot * length(Y), replace = TRUE), nboot, length(Y) )
Fboot <- apply(y, 1, Fstat, group=G)
F0 <- Fstat(Y, G)
p1 <- sum(Fboot >= F0) / nboot

# one-box residuals
y <- matrix(sample( Y-ave( Y, G, FUN=center), size = nboot * length(Y), replace = TRUE), nboot, length(Y) )  # residuals
Fboot <- apply(y, 1, Fstat, group=G)
F0 <- Fstat(Y, G)
p1r <- sum(Fboot >= F0) / nboot

# one-box ranks
y <- matrix(sample( Y, size = nboot * length(Y), replace = TRUE), nboot, length(Y) )  # residuals
y.rank <- t(apply(y, 1, rank))
Fboot <- apply(y.rank, 1, Fstat, group=G)
F0 <- Fstat(rank(Y), G)
p.rank <- sum(Fboot >= F0) / nboot

# many-box residuals
residuals <- Y-ave( Y, G, FUN=center)
y <- replicate(nboot, ave(residuals, G, FUN=function(x) sample(x, replace=TRUE)), simplify="array")
Fboot <- apply(y, 2, Fstat, group=G)
F0 <- Fstat(Y, G)
p2r <- sum(Fboot >= F0) / nboot

# one-box
y <- matrix(sample(Y, size = nboot * length(Y), replace = TRUE), nboot, length(Y) )
Fboot <- apply(y, 1, Fstat, group=G, agg=sab)
F0 <- Fstat(Y, G, agg=sab)
p1ab <- sum(Fboot >= F0) / nboot

# one-box residuals
y <- matrix(sample( Y-ave( Y, G, FUN=center), size = nboot * length(Y), replace = TRUE), nboot, length(Y) )  # residuals
Fboot <- apply(y, 1, Fstat, group=G, agg=sab)
F0 <- Fstat(Y, G, agg=sab)
p1rab <- sum(Fboot >= F0) / nboot

# many-box residuals
residuals <- Y-ave( Y, G, FUN=center)
y <- replicate(nboot, ave(residuals, G, FUN=function(x) sample(x, replace=TRUE)), simplify="array")
Fboot <- apply(y, 2, Fstat, group=G, agg=sab)
F0 <- Fstat(Y, G, agg=sab)
p2rab <- sum(Fboot >= F0) / nboot



# one-box
y <- matrix(sample(Y, size = nboot * length(Y), replace = TRUE), nboot, length(Y) )
Fboot <- apply(y, 1, Fstat, group=G, center=median, agg=sab)
F0 <- Fstat(Y, G, center=median, agg=sab)
p1mab <- sum(Fboot >= F0) / nboot

# one-box residuals
y <- matrix(sample( Y-ave( Y, G, FUN=center), size = nboot * length(Y), replace = TRUE), nboot, length(Y) )  # residuals
Fboot <- apply(y, 1, Fstat, group=G, center=median, agg=sab)
F0 <- Fstat(Y, G, agg=sab, center=median)
p1rmab <- sum(Fboot >= F0) / nboot


# many-box residuals
residuals <- Y-ave( Y, G, FUN=center)
y <- replicate(nboot, ave(residuals, G, FUN=function(x) sample(x, replace=TRUE)), simplify="array")
Fboot <- apply(y, 2, Fstat, group=G, center=median, agg=sab)
F0 <- Fstat(Y, G, center=median, agg=sab)
p2rmab <- sum(Fboot >= F0) / nboot

tests <- c(p.standard, p.kw, p1, p1r, p.rank, p2r,
  p1ab, p1rab, p2rab,
  p1mab, p1rmab, p2rmab)

names(tests) <- c("Standard ANOVA", "Kruskal-Wallis",  "one-box L2", "one-box residuals L2", "rank", "many-box residuals L2", 
                  "one-box L1", "one-box residuals L1", "many-box residuals L1",
                  "one-box, median L1", "one-box residuals, median L1", "many-box residuals, median L1")

ptab <- cbind(tests[c(1,3,4,7,8,10,11,6,9,12,2,5)])
colnames(ptab) <- "p-value"

library(knitr)
kable(round(ptab,4))
# some test stuff
F0 <- Fstat(df$weight,df$group, center=mean.fast, agg=ssq)

# Bootstrap the F-statistic to compute the p-value
bootF <- rep(NA,deals)
residuals <- df$weight - ave(df$weight, df$group)
for(i in 1:deals){
    df$boot <- ave(residuals, df$group, FUN=function(x) sample(x, replace=TRUE))
    bootF[i] <- Fstat(df$boot, df$group, center=mean, agg=ssq)
}
pval.F <- sum(bootF > F0)/deals


# many-box residuals
residuals <- Y-ave( Y, G, FUN=center)
y <- replicate(nboot, ave(residuals, G, FUN=function(x) sample(x, replace=TRUE)), simplify="array")
Fboot <- apply(y, 2, Fstat, group=G)
F0 <- Fstat(Y, G)
p2r <- sum(Fboot >= F0) / nboot


