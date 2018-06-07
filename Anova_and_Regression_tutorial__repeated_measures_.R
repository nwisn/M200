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
# Initialization
set.seed(100)
n <- 15
sigma <- 100

# Monte Carlo 
t1normal <- rnorm(n, mean=600, sd=sigma)
t1super <- rnorm(n, mean=600, sd=sigma)
t2normal <- rnorm(n, mean=600, sd=sigma)
t2super <- rnorm(n, mean=650, sd=sigma)
t3normal <- rnorm(n, mean=600, sd=sigma)
t3super <- rnorm(n, mean=800, sd=sigma)

# Put together into a dataframe
df <- data.frame(distance = c(t1normal, t1super, t2normal, t2super, t3normal, t3super),
                 time = c(rep("time1", 2*n), rep("time2", 2*n), rep("time3", 2*n)),
                 diet =  rep(c(rep("normal", n), rep("superfood", n)), 3),
                 frog = rep(1:30,3))
head(df,20)
dev.off()
par(mar=c(10,3,2,2))
boxplot(distance ~ time + diet, horizontal=F, las=2, data=df, col=c(rep("red",3),rep("blue",3)))
stripchart(distance ~ time + diet, data=df, add=T, method='jitter', vertical=T, pch=16, col="grey")
abline(h=mean(df$distance), lty='dotted')  # draw dotted line at the grand mean
par(mar=c(7,5,2,2))
par(mfrow=c(1,2))
boxplot(distance ~ time, horizontal=F, las=2, data=subset(df, df$diet=="normal"), main="normal diet", ylab="distance jumped (cm)")
stripchart(distance ~ time, data=subset(df, df$diet=="normal"), add=T, method='jitter', vertical=T, pch=16, col="red")

boxplot(distance ~ time, horizontal=F, las=2, data=subset(df, df$diet=="superfood"), main="superfood diet", ylab="distance jumped (cm)")
stripchart(distance ~ time, data=subset(df, df$diet=="superfood"), add=T, method='jitter', vertical=T, pch=16, col="blue")
library(lme4)
library(lmerTest)
model <- lmer(distance ~ diet * time + (1|frog), data=df)
result <- anova(model)
pval.lm <- result$`Pr(>F)`
result
interaction.plot(df$time, df$diet, df$distance, 
                 fun = mean,
                 trace.label = "Diet",
                 xlab="Frog time", 
                 ylab="distance jumped",
                 col=c("red","blue"),
                 lwd = 3, lty = 1)

par(mar=c(7,5,2,2))
par(mfrow=c(1,2))
plot(distance ~ time + diet, data=df, las=2, xlab="",  ylab="distance jumped (cm)")
# get the F-statistic from lm()
F0 = result$`Pr(>F)`

# Bootstrap the F-statistics to compute the p-values
deals <- 1000
bootF <- matrix(NA, nrow=deals, ncol=length(result$`Pr(>F)`))
for(i in 1:deals){
    df$boot <- sample(df$distance, replace=T)
    bootresult <- anova(lmer(boot ~ diet * time + (1|frog), data=df))
    bootF[i,] <- bootresult$`Pr(>F)`
}
pval.lmboot <- colSums(sweep(bootF, 2, F0,"-") > 0) / deals
par(mfrow=c(1,3))
xlabel <- c("F diet", "F time", "F interaction")
for(i in 1:3){
    hist(bootF[,i], breaks='FD', main= paste("Bootstrapped p=",pval.lmboot[i],sep=""), xlab=xlabel[i])
    abline(v=F0[i], col='red')
    text(x=F0[i], y=100, labels=paste("p=",pval.lmboot[i],sep=""))
}
df.time <- split(df, df$time)

anova(lm(distance ~ diet, data=df.time$time1)) # time 1
anova(lm(distance ~ diet, data=df.time$time2)) # time 2
anova(lm(distance ~ diet, data=df.time$time3)) # time 3
par(mfrow=c(1,3))
a1 <- anova(lm(distance ~ diet, data=df.time$time1)) # time 1 
a2 <- anova(lm(distance ~ diet, data=df.time$time2)) # time 2
a3 <- anova(lm(distance ~ diet, data=df.time$time3)) # time 3
boxplot(distance ~ diet, data=df.time$time1, ylab="distance jumped (cm)", main=c("time 1",paste("p=",round(a1$`Pr(>F)`[1],4))))
boxplot(distance ~ diet, data=df.time$time2, ylab="distance jumped (cm)", main=c("time 2",paste("p=",round(a2$`Pr(>F)`[1],4))))
boxplot(distance ~ diet, data=df.time$time3, ylab="distance jumped (cm)", main=c("time 3",paste("p=",signif(a3$`Pr(>F)`[1],2))))

df.diet <- split(df, df$diet)

anova(lmer(distance ~ time + (1|frog), data=df.diet$normal)) #normal
anova(lmer(distance ~ time + (1|frog), data=df.diet$superfood)) #superfood
a1 <- anova(lmer(distance ~ time + (1|frog), data=df.diet$normal)) #normal
a2 <- anova(lmer(distance ~ time + (1|frog), data=df.diet$superfood)) #superfood
par(mfrow=c(1,2))
boxplot(distance ~ time, data=df.diet$normal, ylab="distance jumped (cm)", main=c("Normal food", paste("p=", round(a1$`Pr(>F)`[1],4))))
boxplot(distance ~ time, data=df.diet$superfood, ylab="distance jumped (cm)", main=c("Superfood", paste("p=", round(a2$`Pr(>F)`[1],4))))
pairwise.t.test(df.diet$superfood$distance, df.diet$superfood$time, paired=T, p.adjust="bonferroni") # Super diet
bartlett.test(distance ~ interaction(diet,time), data=df)
fligner.test(distance ~ interaction(diet,time), data=df)
library(car)
leveneTest(distance ~ interaction(diet,time), data=df)
plot(model)
library(lattice)
lattice::xyplot(distance~time| frog, groups=diet, data=df, type=c('p','r'), auto.key=F)
lattice::xyplot(fitted(model)~time| frog, groups=diet, data=df, type=c('p','r'), auto.key=F)
model <-  lmer(rank(distance) ~ diet * time + (1|frog), data=df)
result <- anova(model)
pval.rank <- result$`Pr(>F)`
result
kruskal.test(distance~diet, data=df.time$time1) # diet at t1
kruskal.test(distance~diet, data=df.time$time2) # diet at t2
kruskal.test(distance~diet, data=df.time$time3) # diet at t3
par(mfrow=c(1,3))
a1 <- kruskal.test(distance~diet, data=df.time$time1) # diet at t1
a2 <- kruskal.test(distance~diet, data=df.time$time2) # diet at t2
a3 <- kruskal.test(distance~diet, data=df.time$time3) # diet at t3
boxplot(distance ~ diet, data=df.time$time1, ylab="distance jumped (cm)", main=c("time 1",paste("p=",round(a1$p.value,4))))
boxplot(distance ~ diet, data=df.time$time2, ylab="distance jumped (cm)", main=c("time 2",paste("p=",round(a2$p.value,4))))
boxplot(distance ~ diet, data=df.time$time3, ylab="distance jumped (cm)", main=c("time 3",paste("p=",signif(a3$p.value,2))))

friedman.test(distance~time|frog, data=df.diet$normal) # time differences in normal diet subgroup
friedman.test(distance~time|frog, data=df.diet$superfood) # time differences in superfood diet subgroup
a1 <- friedman.test(distance~time|frog, data=df.diet$normal) # time differences in normal diet subgroup
a2 <- friedman.test(distance~time|frog, data=df.diet$superfood) # time differences in superfood diet subgroup
par(mfrow=c(1,2))
boxplot(distance ~ time, data=df.diet$normal,  ylab="distance jumped (cm)",main=c("Normal food", paste("p=", round(a1$p.value,4))))
boxplot(distance ~ time, data=df.diet$superfood, ylab="distance jumped (cm)", main=c("Superfood", paste("p=", round(a2$p.value,4))))
pairwise.wilcox.test(df.diet$superfood$distance, df.diet$superfood$time, paired=T, p.adjust="bonferroni") # Superfood
par(mfrow=c(1,1))
boxplot(distance ~ time, data=df.diet$superfood, ylab="distance jumped (cm)", main=c("Superfood"))
# Define functions to compute the Sum of Squares, and Sum of Abs
ssq <- function(x) sum(x^2)
sab <- function(x) sum(abs(x))

# Compute the SS_subject term
SSsubject <- function( X, subject, f2, center=mean, agg=ssq){
    long <- data.frame(cbind(subject, f2, X))
    wide <- reshape(long, v.names="X", idvar = "subject", timevar = "f2",direction = "wide")
    k <- ncol(wide[,-1])
    subjectMeans <- rowMeans(wide[,-1])
    SSsubject <- k*agg(subjectMeans - center(X))
    return( SSsubject)
}

# A more general function to compute the F-statistic
Fstat.2wayRM <- function( X, f1, f2, subject, center=mean, agg=ssq, rank=FALSE) { 
    if (rank==FALSE){Y = X} else{Y = rank(-X)}
    mu <- center(Y)  #Grand Mean
    groupMeans.f1 <- ave( Y, f1, FUN=center ) #Pool f1
    groupMeans.f2 <- ave( Y, f2, FUN=center ) #Pool f2
    groupMeans.int <- ave( Y, f1, f2, FUN=center ) #Cell means
    SSB.f1 <- agg(groupMeans.f1 - mu)
    SSB.f2 <- agg(groupMeans.f2 - mu)
    SSB.int <- agg( groupMeans.int - (groupMeans.f1 + groupMeans.f2 - mu) )
    SSW <- agg(Y - groupMeans.int)
    SSsub <- SSsubject(Y, subject, f2, center=center, agg=agg)
    SSe <- SSW - SSsub
    return( c(f1=SSB.f1 , f2=SSB.f2, int=SSB.int) /  SSe )
}

# Recenter each of the timepoints before resampling
center.time <- function( X, f2, center=mean){
    timeMeans <- ave( X, f2, FUN=center )
    return(X - timeMeans)
}

# Convert wide dataframe to long
wideToLong <- function(wide, timecols, timenames){
    long <- reshape(wide, varying = timecols, v.names = "X", timevar = "f2", times = timenames, direction = "long")
    return(long)
}


X <- df$distance
f1 <- df$diet
f2 <- df$time  # f2 must be time in order for the code to work
subject <- df$frog
deals <- 10000

# initialize some new dataframes
Y <- center.time(X, f2, center=mean) # recenter timepoints
long <- data.frame(subject, f1, f2, Y)  # put into a long dataframe
wide <- reshape(long, v.names="Y", timevar="f2", idvar="subject", direction="wide") # convert to short dataframe to resample subjects

# initialize bootstrap holders
boot.wide <- wide # initialize
bootF <- matrix(NA, nrow=deals, ncol=3) # preallocate
indices <- 1:nrow(wide) # initialize

# Bootstrap
F0 <- Fstat.2wayRM(X, f1, f2, subject, center=mean, agg=ssq, rank=F)
for(i in 1:deals){
    boot.i <- sample(indices, replace=T)  # resample subjects
    boot.wide[,3:5] <- wide[boot.i,3:5]  # get resampled data, which keeps the temporal samples linked
    boot.long <- wideToLong(boot.wide, names(boot.wide)[3:5])  # convert to long format
    bootF[i,] <- Fstat.2wayRM(boot.long$X, boot.long$f1, boot.long$f2, boot.long$subject, center=mean, agg=ssq, rank=F)
}
pval.L2 <- colSums(sweep(bootF, 2, F0,"-") > 0) / deals
   
par(mfrow=c(1,3))
xlabel <- c("F diet", "F time", "F interaction")
for(i in 1:3){
    hist(bootF[,i], breaks='FD',  main= paste("Bootstrapped p=",pval.L2[i],sep=""), xlab=xlabel[i])
    abline(v=F0[i], col='red')
    text(x=F0[i], y=100, labels=paste("p=",pval.L2[i],sep=""))
}
dof <- result$Df
titlevec <- c("QQ plot for F diet", "QQ plot for F time", "QQ plot for F interaction")
par(mfrow=c(1,3))
for (i in 1:3){
    qqplot( qf( seq(0,1,length=deals),dof[i], dof[4]), bootF, 
        main=titlevec[i], 
        xlab="Theoretical F Quantiles", 
        ylab="Bootstrap F-like Quantiles")
abline(0, dof[i]/dof[4], col='red')
}

F0 <- Fstat.2wayRM(X, f1, f2, subject, center=mean, agg=sab, rank=F)
for(i in 1:deals){
    boot.i <- sample(indices, replace=T)  # resample subjects
    boot.wide[,3:5] <- wide[boot.i,3:5]  # get resampled data, which keeps the temporal samples linked
    boot.long <- wideToLong(boot.wide, names(boot.wide)[3:5])  # convert to long format
    bootF[i,] <- Fstat.2wayRM(boot.long$X, boot.long$f1, boot.long$f2, boot.long$subject, center=mean, agg=sab, rank=F)
}
pval.L1 <- colSums(sweep(bootF, 2, F0,"-") > 0) / deals

par(mfrow=c(1,3))
xlabel <- c("F diet", "F time", "F interaction")
for(i in 1:3){
    hist(bootF[,i], breaks='FD',  main= paste("Bootstrapped p=",pval.L1[i],sep=""), xlab=xlabel[i])
    abline(v=F0[i], col='red')
    text(x=F0[i], y=100, labels=paste("p=",pval.L1[i],sep=""))
}
F0 <- Fstat.2wayRM(X, f1, f2, subject, center=mean, agg=ssq, rank=T)
for(i in 1:deals){
    boot.i <- sample(indices, replace=T)  # resample subjects
    boot.wide[,3:5] <- wide[boot.i,3:5]  # get resampled data, which keeps the temporal samples linked
    boot.long <- wideToLong(boot.wide, names(boot.wide)[3:5])  # convert to long format
    bootF[i,] <- Fstat.2wayRM(boot.long$X, boot.long$f1, boot.long$f2, boot.long$subject, center=mean, agg=ssq, rank=T)
}
pval.rankboot <- colSums(sweep(bootF, 2, F0,"-") > 0) / deals

par(mfrow=c(1,3))
xlabel <- c("F diet", "F time", "F interaction")
for(i in 1:3){
    hist(bootF[,i], breaks='FD',  main= paste("Bootstrapped p=",pval.rankboot[i],sep=""), xlab=xlabel[i])
    abline(v=F0[i], col='red')
    text(x=F0[i], y=100, labels=paste("p=",pval.rankboot[i],sep=""))
}
par(mfrow=c(1,3))
for (i in 1:3){
    qqplot( qf( seq(0,1,length=deals),dof[i], dof[4]), bootF, 
        main=titlevec[i], 
        xlab="Theoretical F Quantiles", 
        ylab="Bootstrap F-like Quantiles")
abline(0, dof[i]/dof[4], col='red')
}
comparison <- data.frame(LM=round(pval.lm[-4],4),  
                         #LMboot=pval.lmboot[-4], 
                         L1=pval.L1, 
                         L2=pval.L2, 
                         rank=round(pval.rank[-4],4), 
                         rankboot=pval.rankboot)
rownames(comparison) <- c("diet", "time", "interaction")
colnames(comparison) <- c("built-in ANOVA", "bootstrap with L1-norm", "bootstrap with L2-norm", "rank", "bootstrap rank")
table <- t(comparison)
library(knitr)
kable(table, digits=4)
