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
effectsizes <- seq(0,2,.25)
samplesizes <- c(4,6,8,10,20,30,40)
library(pwr)
power.theory <- array(NA, c(length(effectsizes), length(samplesizes)) )
for (i in 1:length(effectsizes)){
    power.out <- pwr.t.test(n = samplesizes/2 , d = effectsizes[i], sig.level = 0.05, type = "two.sample")
    power.theory[i,] <- power.out$power
}
# Plot power curve
matplot(power.theory, type='l', col=1:length(samplesizes), xlab="effect size", ylab="power", axes=F, main="Theoretical Power") 
legend("bottomright", inset=.05, legend=samplesizes, pch=16, horiz=TRUE, col=1:length(samplesizes))
axis(side=1,at=1:length(effectsizes),labels=effectsizes)
axis(2)
abline(h=.8)

nboot=1000
test <- array(NA, c(length(effectsizes), length(samplesizes), nboot))
for (i in 1:length(effectsizes)){
    for (j in 1:length(samplesizes)){
        for (k in 1:nboot){
            x <- rnorm(samplesizes[j]/2, 0 + effectsizes[i], 1)
            y <- rnorm(samplesizes[j]/2, 0, 1)
            test.out <- t.test(x,y)  # bootstrap test statistic
            test[i,j,k] <- test.out$p.value < 0.05
        }
    }
}
power.mc <- array(NA, c(length(effectsizes), length(samplesizes)) )
for (i in 1:length(effectsizes)){
    for (j in 1:length(samplesizes)){
        power.mc[i,j] <- sum(test[i,j,])/nboot
    }
}
matplot(power.mc, type='l', col=1:length(samplesizes), xlab="effect size", ylab="power", axes=F, main="Monte Carlo Power") 
legend("bottomright", inset=.05, legend=samplesizes, pch=16, horiz=TRUE, col=1:length(samplesizes))
axis(side=1,at=1:length(effectsizes),labels=effectsizes)
axis(2)
abline(h=.8)
bigN <- 20  # total sample size (both groups combined)
X <- rnorm(bigN/2, 0,1)
Y <- rnorm(bigN/2, 1,1)
dstat <- function(x1, x2, center=mean) center(x1) - center(x2) # test statistic
d0 <- dstat(X,Y)  # actual test statistic
dnull <- array(NA, c(length(samplesizes), nboot)); 
for (j in 1:length(samplesizes)){
    for (k in 1:nboot){
        x <- sample( c(X,Y), samplesizes[j]/2, replace=T)  # resample new x group
        y <- sample( c(X,Y), samplesizes[j]/2, replace=T)  # resample new y group
        dnull[j,k] <- dstat(x, y)  # bootstrap test statistic
    }
}
crit <- apply(dnull,1,quantile, c(0.025,0.975))
deffect <- array(nboot); 
for (k in 1:nboot){
    x <- sample( c(X,Y) + effectsizes[4], bigN/2, replace=T)  # resample new x group
    y <- sample( c(X,Y), bigN/2, replace=T)  # resample new y group
    deffect[k] <- dstat(x, y)  # bootstrap test statistic
}
hist(dnull[5,], breaks="FD", xlim=c(-2,3), col='blue', xlab="Effect Size (difference in means)", main=paste("Power at n=",bigN/2, ", d=",effectsizes[4], sep=""))
hist(deffect, breaks="FD", add=T, col='red')
abline(v=crit[,4], lwd=4)
legend("bottomright", inset=.05, legend=c("Null","Effect"), pch=16, horiz=F, col=c("blue","red"))
deffect <- array(NA, c(length(effectsizes), length(samplesizes), nboot)); 
for (i in 1:length(effectsizes)){
    for (j in 1:length(samplesizes)){
        for (k in 1:nboot){
            x <- sample( c(X,Y) + effectsizes[i], samplesizes[j]/2, replace=T)  # resample new x group
            y <- sample( c(X,Y), samplesizes[j]/2, replace=T)  # resample new y group
            deffect[i,j,k] <- dstat(x, y)  # bootstrap test statistic
        }
    }
}
power <- array(NA, c(length(effectsizes), length(samplesizes))); 
for (i in 1:length(effectsizes)){
    for (j in 1:length(samplesizes)){
        power[i,j] <- (sum(deffect[i,j,] <= crit[1,j]) + sum(deffect[i,j,] >= crit[2,j]) ) / nboot  # two sided
    }
}
matplot(power, type='l', col=1:length(samplesizes), xlab="effect size", ylab="power", axes=F, main=paste("Bootrapped from n=",bigN)) 
legend("bottomright", inset=.05, legend=samplesizes, pch=16, horiz=TRUE, col=1:length(samplesizes))
axis(side=1,at=1:length(effectsizes),labels=effectsizes)
axis(2)
abline(h=.8)
bigN <- 2000 # total sample size (both groups combined)
X <- rnorm(bigN/2, 0,1)
Y <- rnorm(bigN/2, 1,1)

dstat <- function(x1, x2, center=mean) center(x1) - center(x2) # test statistic
d0 <- dstat(X,Y)  # actual test statistic
dnull <- array(NA, c(length(samplesizes), nboot)); 
for (j in 1:length(samplesizes)){
    for (k in 1:nboot){
        x <- sample( c(X,Y), samplesizes[j]/2, replace=T)  # resample new x group
        y <- sample( c(X,Y), samplesizes[j]/2, replace=T)  # resample new y group
        dnull[j,k] <- dstat(x, y)  # bootstrap test statistic
    }
}
crit <- apply(dnull,1,quantile, c(0.025,0.975))

deffect <- array(NA, c(length(effectsizes), length(samplesizes), nboot)); 
for (i in 1:length(effectsizes)){
    for (j in 1:length(samplesizes)){
        for (k in 1:nboot){
            x <- sample( c(X,Y) + effectsizes[i], samplesizes[j]/2, replace=T)  # resample new x group
            y <- sample( c(X,Y), samplesizes[j]/2, replace=T)  # resample new y group
            deffect[i,j,k] <- dstat(x, y)  # bootstrap test statistic
        }
    }
}

power <- array(NA, c(length(effectsizes), length(samplesizes))); 
for (i in 1:length(effectsizes)){
    for (j in 1:length(samplesizes)){
        power[i,j] <- (sum(deffect[i,j,] <= crit[1,j]) + sum(deffect[i,j,] >= crit[2,j]) ) / nboot  # two sided
    }
}

matplot(power, type='l', col=1:length(samplesizes), xlab="effect size", ylab="power", axes=F, main=paste("Bootrapped from n=",bigN)) 
legend("bottomright", inset=.05, legend=samplesizes, pch=16, horiz=TRUE, col=1:length(samplesizes))
axis(side=1,at=1:length(effectsizes),labels=effectsizes)
axis(2)
abline(h=.8)
bigN <- 6 # total sample size (both groups combined)
X <- rnorm(bigN/2, 0,1)
Y <- rnorm(bigN/2, 1,1)

dstat <- function(x1, x2, center=mean) center(x1) - center(x2) # test statistic
d0 <- dstat(X,Y)  # actual test statistic
dnull <- array(NA, c(length(samplesizes), nboot)); 
for (j in 1:length(samplesizes)){
    for (k in 1:nboot){
        x <- sample( c(X,Y), samplesizes[j]/2, replace=T)  # resample new x group
        y <- sample( c(X,Y), samplesizes[j]/2, replace=T)  # resample new y group
        dnull[j,k] <- dstat(x, y)  # bootstrap test statistic
    }
}
crit <- apply(dnull,1,quantile, c(0.025,0.975))

deffect <- array(NA, c(length(effectsizes), length(samplesizes), nboot)); 
for (i in 1:length(effectsizes)){
    for (j in 1:length(samplesizes)){
        for (k in 1:nboot){
            x <- sample( c(X,Y) + effectsizes[i], samplesizes[j]/2, replace=T)  # resample new x group
            y <- sample( c(X,Y), samplesizes[j]/2, replace=T)  # resample new y group
            deffect[i,j,k] <- dstat(x, y)  # bootstrap test statistic
        }
    }
}

power <- array(NA, c(length(effectsizes), length(samplesizes))); 
for (i in 1:length(effectsizes)){
    for (j in 1:length(samplesizes)){
        power[i,j] <- (sum(deffect[i,j,] <= crit[1,j]) + sum(deffect[i,j,] >= crit[2,j]) ) / nboot  # two sided
    }
}

matplot(power, type='l', col=1:length(samplesizes), xlab="effect size", ylab="power", axes=F, main=paste("Bootrapped from n=",bigN)) 
legend("bottomright", inset=.05, legend=samplesizes, pch=16, horiz=TRUE, col=1:length(samplesizes))
axis(side=1,at=1:length(effectsizes),labels=effectsizes)
axis(2)
abline(h=.8)
bigN <- 20 # total sample size (both groups combined)
X <- rnorm(bigN/2, 0,1)
Y <- rnorm(bigN/2, 1,1)

d0 <- dstat(X,Y)  # actual test statistic
dnull <- array(NA, c(length(samplesizes), nboot)); 
for (j in 1:length(samplesizes)){
    for (k in 1:nboot){
        x <- sample( X-d0, samplesizes[j]/2, replace=T)  # resample new x group
        y <- sample( Y, samplesizes[j]/2, replace=T)  # resample new y group
        dnull[j,k] <- dstat(x, y)  # bootstrap test statistic
    }
}
crit <- apply(dnull,1,quantile, c(0.025,0.975))
deffect <- array(NA, c(length(effectsizes), length(samplesizes), nboot)); 
for (i in 1:length(effectsizes)){
    for (j in 1:length(samplesizes)){
        for (k in 1:nboot){
            x <- sample( X-d0+effectsizes[i], samplesizes[j]/2, replace=T)  # resample new x group
            y <- sample( Y, samplesizes[j]/2, replace=T)  # resample new y group
            deffect[i,j,k] <- dstat(x, y)  # bootstrap test statistic
        }
    }
}
power <- array(NA, c(length(effectsizes), length(samplesizes))); 
for (i in 1:length(effectsizes)){
    for (j in 1:length(samplesizes)){
        power[i,j] <- (sum(deffect[i,j,] <= crit[1,j]) + sum(deffect[i,j,] >= crit[2,j]) ) / nboot  # two sided
    }
}
matplot(power, type='l', col=1:length(samplesizes), xlab="effect size", ylab="power", axes=F,  main=paste("n=",bigN)) 
legend("bottomright", inset=.05, legend=samplesizes, pch=16, horiz=TRUE, col=1:length(samplesizes))
axis(side=1,at=1:length(effectsizes),labels=effectsizes)
axis(2)
abline(h=.8)
