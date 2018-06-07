load("/Users/drwho/equalvar.RData")
load("/Users/drwho/unequalvar.RData")
load("/Users/drwho/skew.RData")

FPRplot <- function(data, sigma, xi, nboot, nsim){
    data2 <- subset(data, data$delta==0)
    powerplot <- ggplot(data = data2, aes(x=n, y=value, group=variable)) + geom_line(aes(colour=variable, linetype=variable), size=1) + 
        geom_hline(aes(yintercept=.05), colour="black", linetype="dashed") +
        ggtitle(paste("FPR (sd=", sigma, ", xi=",xi, ", nboot=",nboot,", nsim=",nsim, ")", sep="")) +
        #geom_point() + 
        ylim(0, max(data2$value)) +
        labs(x="sample size",y="false positive rate")  + theme_linedraw() + 
        theme(legend.key.size = unit(1.5, "cm")) +
        theme(legend.text = element_text(size = 12, colour = "black", angle = 0)) +
        theme(axis.title.y = element_text(size = rel(1.5), angle = 90)) +
        theme(axis.title.x = element_text(size = rel(1.5), angle = 0)) +
        theme(plot.title = element_text(size = rel(1.1), angle = 0)) 
    return(powerplot)
}

powerplot <- function(data, samplesize, sigma, xi, nboot, nsim){
    powerplot <- ggplot(data = subset(data, data$n==samplesize & data$sigma==sigma & data$xi==xi), aes(x=delta, y=value, group=variable)) + geom_line(aes(colour=variable, linetype=variable), size=1) + 
        geom_hline(aes(yintercept=.8), colour="#990000", linetype="dashed") +
        geom_hline(aes(yintercept=.05), colour="#990000", linetype="dashed") +
        ggtitle(paste("Power Curves\n n=", samplesize, ", s=", sigma, ", xi=",xi, "\nnboot=",nboot,", nsim=",nsim, sep="")) +
        #geom_point() + 
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
require(gridExtra)


#-------------------------------
# Equal variance
result.noKS <- subset(result, result$variable !="KS")
FPRplot(result.noKS, sigma=1, xi=0, nboot=1000, nsim=1000)

result.0 <- subset(result.noKS, result.noKS$delta==0)
result.0

grid.arrange(powerplot(result.noKS, 4,1,0 ,1000,1000), powerplot(result.noKS, 6,1,0,1000,1000), 
             powerplot(result.noKS, 8,1,0,1000,1000), powerplot(result.noKS, 10,1,0,1000,1000),
             powerplot(result.noKS, 20,1,0,1000,1000), powerplot(result.noKS, 30,1,0,1000,1000),
             #powerplot(result.noKS, 60,1,0), powerplot(result.noKS, 100,1,0),
             ncol=2)


# Unequal variance
result.uq.noKS <- subset(result.unequal, result.unequal$variable !="KS")
FPRplot(result.uq.noKS, sigma=4, xi=0, nboot=1000, nsim=1000)

result.uq.0 <- subset(result.uq.noKS, result.uq.noKS$delta==0)
result.uq.0

grid.arrange(powerplot(result.uq.noKS, 4,4,0 ,1000,1000), powerplot(result.uq.noKS, 6,4,0,1000,1000), 
             powerplot(result.uq.noKS, 8,4,0,1000,1000), powerplot(result.uq.noKS, 10,4,0,1000,1000),
             powerplot(result.uq.noKS, 20,4,0,1000,1000), powerplot(result.uq.noKS, 30,4,0,1000,1000),
             #powerplot(result, 60,1,0), powerplot(result, 100,1,0),
             ncol=2)


# Skew
result.skew.noKS <- subset(result.skew, result.skew$variable !="KS")
FPRplot(result.skew.noKS, sigma=1, xi=1000, nboot=1000, nsim=1000)

result.skew.0 <- subset(result.skew.noKS, result.skew.noKS$delta==0)
result.skew.0

grid.arrange(powerplot(result.skew.noKS, 4,1,1000 ,1000,1000), powerplot(result.skew.noKS, 6,1,1000,1000,1000), 
             powerplot(result.skew.noKS, 8,1,1000,1000,1000), powerplot(result.skew.noKS, 10,1,1000,1000,1000),
             powerplot(result.skew.noKS, 20,1,1000,1000,1000), powerplot(result.skew.noKS, 30,1,1000,1000,1000),
             #powerplot(result, 60,1,0), powerplot(result, 100,1,0),
             ncol=2)




set.seed(100)
nboot <-300
nsim= 300
estimator=median

samplesize <- 50
p2r=numeric()
p1 = numeric()
p1r = numeric()
for (j in 1:nsim){
    x <- rnorm(samplesize/2, 0, 1)
    y <- rnorm(samplesize/2, 0, 4)
    d0 <- estimator(x) - estimator(y)
    d1 <- numeric()
    d1r <- numeric()
    d2r <- numeric()
    for (i in 1:nboot){
        print(paste(i,j, sep=","))
        xy1 <- sample(c(x,y), replace=T)
        xy1r <- sample(c(x-estimator(x),y-estimator(y)), replace=T)
        x2 <- sample(x-estimator(x), replace=T)
        y2 <- sample(y-estimator(y), replace=T)
        d1[i] <- estimator(xy1[1:(samplesize/2)]) - estimator(xy1[-(1:(samplesize/2))])
        d1r[i] <- estimator(xy1r[1:(samplesize/2)]) - estimator(xy1r[-(1:(samplesize/2))])
        d2r[i] <- estimator(x2) - estimator(y2)
    }
    p2r[j] <- sum(abs(d2r)>=d0)/nboot
    p1[j] <- sum(abs(d1)>=d0)/nboot
    p1r[j] <- sum(abs(d1r)>=d0)/nboot
}

sum(p2r < 0.05)/nsim
sum(p1 < 0.05)/nsim
sum(p1r < 0.05)/nsim

hist(d2r, breaks="FD", col="red")
hist(d1, breaks="FD", add=T, col="blue")
hist(d1r, breaks="FD", add=T, col="green")
abline(v=d0)


boxplot(x,y)
boxplot(x-estimator(x),y-estimator(y))



boxplot(xy1, xy1r, x2,y2)











load("/Users/drwho/Dropbox/Encrypted/M200 Advanced Experimental Statistics/M200.4/twogroup/equalvar.RData")
result <- subset(result, result$variable!="KS")
result.unequal <- subset(result.unequal, result.unequal$variable!="KS")
result.skew <- subset(result.skew, result.skew$variable!="KS")


grid.arrange(powerplot(result, 4,1,0), powerplot(result, 6,1,0), 
             powerplot(result, 8,1,0), powerplot(result, 10,1,0),
             powerplot(result, 20,1,0), powerplot(result, 30,1,0),
             #powerplot(result, 60,1,0), powerplot(result, 100,1,0),
             ncol=2)

grid.arrange(powerplot(result.skew, 4,1,0), powerplot(result.skew, 6,1,0), 
             powerplot(result.skew, 8,1,0), powerplot(result.skew, 10,1,0),
             powerplot(result.skew, 80,1,0), powerplot(result.skew, 100,1,0),
             #powerplot(result, 60,1,0), powerplot(result, 100,1,0),
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

ntemp = 30; dtemp = 20/8
temp.un <- subset(result.unequal, result.unequal$n==ntemp & result.unequal$delta==dtemp)
temp <- subset(result, result$n==ntemp & result$delta==dtemp)
temp.sk <- subset(result.skew, result.skew$n==ntemp & result.skew$delta==dtemp)
temp[order(temp$value, decreasing=T), ]
temp.un[order(temp.un$value, decreasing=T), ]
temp.sk[order(temp.sk$value, decreasing=T), ]

ntemp = 30; dtemp = 15/8
temp.un2 <- subset(result.unequal, result.unequal$n==ntemp & result.unequal$delta %in% c(0,dtemp) )
temp.un2.w <- reshape(temp.un2, timevar="delta", idvar=c("n","sigma","xi","variable"), direction = "wide")
temp.un2.w$corrected <- temp.un2.w[,6] - temp.un2.w[,5]
temp.un2.w[order(temp.un2.w$corrected, decreasing=T),]


ntemp = 30; dtemp = 15/8
temp.un3 <- subset(result.skew, result.skew$n==ntemp & result.skew$delta %in% c(0,dtemp) )
temp.un3.w <- reshape(temp.un3, timevar="delta", idvar=c("n","sigma","xi","variable"), direction = "wide")
temp.un3.w$corrected <- temp.un3.w[,6] - temp.un3.w[,5]
temp.un3.w[order(temp.un3.w$corrected, decreasing=T),]
