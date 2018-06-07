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
BFnormal <- rnorm(n, mean=600, sd=sigma)
BFsuper <- rnorm(n, mean=600, sd=sigma)
LFnormal <- rnorm(n, mean=600, sd=sigma)
LFsuper <- rnorm(n, mean=650, sd=sigma)
TFnormal <- rnorm(n, mean=600, sd=sigma)
TFsuper <- rnorm(n, mean=800, sd=sigma)

# Put together into a dataframe
df <- data.frame(distance = c(BFnormal, BFsuper, LFnormal, LFsuper, TFnormal, TFsuper),
                 species = c(rep("Bullfrog", 2*n), rep("Leopardfrog", 2*n), rep("Treefrog", 2*n)),
                 diet =  rep(c(rep("normal", n), rep("superfood", n)), 3))
head(df,20)
dev.off()
par(mar=c(10,3,2,2))
boxplot(distance ~ species + diet, horizontal=F, las=2, data=df, col=c(rep("red",3),rep("blue",3)))
stripchart(distance ~ species + diet, data=df, add=T, method='jitter', vertical=T, pch=16, col="grey")
abline(h=mean(df$distance), lty='dotted')  # draw dotted line at the grand mean
par(mar=c(7,5,2,2))
par(mfrow=c(1,2))
boxplot(distance ~ species, horizontal=F, las=2, data=subset(df, df$diet=="normal"), main="normal diet", ylab="distance jumped (cm)")
stripchart(distance ~ species, data=subset(df, df$diet=="normal"), add=T, method='jitter', vertical=T, pch=16, col="red")

boxplot(distance ~ species, horizontal=F, las=2, data=subset(df, df$diet=="superfood"), main="superfood diet", ylab="distance jumped (cm)")
stripchart(distance ~ species, data=subset(df, df$diet=="superfood"), add=T, method='jitter', vertical=T, pch=16, col="blue")
model <- lm(distance ~ diet * species, data=df)
result <- anova(model)
pval.lm <- result$`Pr(>F)`
result

interaction.plot(df$species, df$diet, df$distance, 
                 fun = mean,
                 trace.label = "Diet",
                 xlab="Frog species", 
                 ylab="distance jumped",
                 col=c("red","blue"),
                 lwd = 3, lty = 1)

par(mar=c(7,5,2,2))
par(mfrow=c(1,2))
plot(distance ~ species + diet, data=df, las=2, xlab="",  ylab="distance jumped (cm)")
# get the F-statistic from lm()
F0 = result$`F value`

# Bootstrap the F-statistics to compute the p-values
deals <- 10000
bootF <- matrix(NA, nrow=deals, ncol=length(result$`F value`))
for(i in 1:deals){
    df$boot <- sample(df$distance, replace=T)
    bootresult <- anova(lm(boot ~ diet * species, data=df))
    bootF[i,] <- bootresult$`F value`
}
pval.lmboot <- colSums(sweep(bootF, 2, F0,"-") > 0) / deals
par(mfrow=c(1,3))
xlabel <- c("F diet", "F frog-species", "F interaction")
for(i in 1:3){
    hist(bootF[,i], breaks='FD', main= paste("Bootstrapped p=",pval.lmboot[i],sep=""), xlab=xlabel[i])
    abline(v=F0[i], col='red')
    text(x=F0[i], y=100, labels=paste("p=",pval.lmboot[i],sep=""))
}
df.species <- split(df, df$species)

anova(lm(distance ~ diet, data=df.species$Bullfrog)) # Bullfrogs
anova(lm(distance ~ diet, data=df.species$Leopardfrog)) # Leopardfrogs
anova(lm(distance ~ diet, data=df.species$Treefrog)) # Treefrogs
par(mfrow=c(1,3))
a1 <- anova(lm(distance ~ diet, data=df.species$Bullfrog)) # Bullfrogs
a2 <- anova(lm(distance ~ diet, data=df.species$Leopardfrog)) # Leopardfrogs
a3 <- anova(lm(distance ~ diet, data=df.species$Treefrog)) # Treefrogs
boxplot(distance ~ diet, data=df.species$Bullfrog,  ylab="distance jumped (cm)", main=c("Bull frogs",paste("p=",round(a1$`Pr(>F)`[1],4))))
boxplot(distance ~ diet, data=df.species$Leopardfrog, ylab="distance jumped (cm)", main=c("Leopard frogs",paste("p=",round(a2$`Pr(>F)`[1],4))))
boxplot(distance ~ diet, data=df.species$Treefrog, ylab="distance jumped (cm)", main=c("Tree frogs",paste("p=",signif(a3$`Pr(>F)`[1],2))))

df.diet <- split(df, df$diet)

anova(lm(distance ~ species, data=df.diet$normal)) # Normal diet
anova(lm(distance ~ species, data=df.diet$superfood)) # Super diet
a1 <- anova(lm(distance ~ species, data=df.diet$normal)) # Normal diet
a2 <- anova(lm(distance ~ species, data=df.diet$superfood)) # Super diet
par(mfrow=c(1,2))
boxplot(distance ~ species, data=df.diet$normal, ylab="distance jumped (cm)", main=c("Normal food", paste("p=", round(a1$`Pr(>F)`[1],4))))
boxplot(distance ~ species, data=df.diet$superfood, ylab="distance jumped (cm)", main=c("Superfood", paste("p=", round(a2$`Pr(>F)`[1],4))))
pairwise.t.test(df.diet$superfood$distance, df.diet$superfood$species, p.adjust="bonferroni") # Super diet
bartlett.test(distance ~ interaction(diet,species), data=df)
fligner.test(distance ~ interaction(diet,species), data=df)
library(car)
leveneTest(distance ~ interaction(diet,species), data=df)
par(mfrow=c(2,2))
plot(model)
model <- lm(rank(df$distance) ~ diet * species, data=df)
result <- anova(model)
pval.rank <- result$`Pr(>F)`
result
kruskal.test(distance~diet, data=df.species$Bullfrog) # diet in bullfrogs
kruskal.test(distance~diet, data=df.species$Leopardfrog) # diet in leopardfrogs
kruskal.test(distance~diet, data=df.species$Treefrog) # diet in treefrogs
par(mfrow=c(1,3))
a1 <- kruskal.test(distance~diet, data=df.species$Bullfrog) # diet in bullfrogs
a2 <- kruskal.test(distance~diet, data=df.species$Leopardfrog) # diet in leopardfrogs
a3 <- kruskal.test(distance~diet, data=df.species$Treefrog) # diet in treefrogs
boxplot(distance ~ diet, data=df.species$Bullfrog, ylab="distance jumped (cm)", main=c("Bull frogs",paste("p=",round(a1$p.value,4))))
boxplot(distance ~ diet, data=df.species$Leopardfrog, ylab="distance jumped (cm)", main=c("Leopard frogs",paste("p=",round(a2$p.value,4))))
boxplot(distance ~ diet, data=df.species$Treefrog, ylab="distance jumped (cm)", main=c("Tree frogs",paste("p=",signif(a3$p.value,2))))

kruskal.test(distance~species, data=df.diet$normal) # species differences in normal diet subgroup
kruskal.test(distance~species, data=df.diet$superfood) # species differences in superfood diet subgroup
a1 <- kruskal.test(distance ~ species, data=df.diet$normal) # Normal diet
a2 <- kruskal.test(distance ~ species, data=df.diet$superfood) # Super diet
par(mfrow=c(1,2))
boxplot(distance ~ species, data=df.diet$normal,  ylab="distance jumped (cm)",main=c("Normal food", paste("p=", round(a1$p.value,4))))
boxplot(distance ~ species, data=df.diet$superfood, ylab="distance jumped (cm)", main=c("Superfood", paste("p=", round(a2$p.value,4))))
pairwise.wilcox.test(df.diet$superfood$distance, df.diet$superfood$species, p.adjust="bonferroni") # Superfood
par(mfrow=c(1,1))
boxplot(distance ~ species, data=df.diet$superfood, ylab="distance jumped (cm)", main=c("Superfood"))
# Define functions to compute the Sum of Squares, and Sum of Abs
ssq <- function(x) sum(x^2)
sab <- function(x) sum(abs(x))

# A more general function to compute the F-statistic
Fstat.2way <- function( X, f1, f2, center=mean, agg=ssq, rank=FALSE) { 
    if (rank==FALSE){Y = X} else{Y = rank(-X)}
    mu = center(Y)  #Grand Mean
    groupMeans.f1 = ave( Y, f1, FUN=center ) #Pool f1
    groupMeans.f2 = ave( Y, f2, FUN=center ) #Pool f2
    groupMeans.int = ave( Y, f1, f2, FUN=center ) #Cell means
    SSB.f1 = agg(groupMeans.f1 - mu)
    SSB.f2 = agg(groupMeans.f2 - mu)
    SSB.int = agg( groupMeans.int - (groupMeans.f1 + groupMeans.f2 - mu) )
    SSW = agg(Y - groupMeans.int)
    return( c(f1=SSB.f1 , f2=SSB.f2, int=SSB.int) /  SSW )
}

X <- df$distance
f1 <- df$diet
f2 <- df$species
deals <- 10000

F0 <- Fstat.2way(X, f1, f2, center=mean, agg=ssq, rank=F)
bootF <- matrix(NA, nrow=deals, ncol=3) #preallocate
for(i in 1:deals){
    bootX <- sample(X, replace=T)
    bootF[i,] <- Fstat.2way(bootX, f1, f2, center=mean, agg=ssq, rank=F)
}
pval.L2 <- colSums(sweep(bootF, 2, F0,"-") > 0) / deals
   
par(mfrow=c(1,3))
xlabel <- c("F diet", "F frog-species", "F interaction")
for(i in 1:3){
    hist(bootF[,i], breaks='FD',  main= paste("Bootstrapped p=",pval.L2[i],sep=""), xlab=xlabel[i])
    abline(v=F0[i], col='red')
    text(x=F0[i], y=100, labels=paste("p=",pval.L2[i],sep=""))
}
dof <- result$Df
titlevec <- c("QQ plot for F diet", "QQ plot for F frog-species", "QQ plot for F interaction")
par(mfrow=c(1,3))
for (i in 1:3){
    qqplot( qf( seq(0,1,length=deals),dof[i], dof[4]), bootF, 
        main=titlevec[i], 
        xlab="Theoretical F Quantiles", 
        ylab="Bootstrap F-like Quantiles")
abline(0, dof[i]/dof[4], col='red')
}

F0 <- Fstat.2way(X, f1, f2, center=mean, agg=sab, rank=F)
bootF <- matrix(NA, nrow=deals, ncol=3) #preallocate
for(i in 1:deals){
    bootX <- sample(X, replace=T)
    bootF[i,] <- Fstat.2way(bootX, f1, f2, center=mean, agg=sab, rank=F)
}
pval.L1 <- colSums(sweep(bootF, 2, F0,"-") > 0) / deals
par(mfrow=c(1,3))
xlabel <- c("F diet", "F frog-species", "F interaction")
for(i in 1:3){
    hist(bootF[,i], breaks='FD',  main= paste("Bootstrapped p=",pval.L1[i],sep=""), xlab=xlabel[i])
    abline(v=F0[i], col='red')
    text(x=F0[i], y=100, labels=paste("p=",pval.L1[i],sep=""))
}
# Run it on the dataframe to get the test statistic:
F0 <- Fstat.2way(X, f1, f2, center=mean, agg=ssq, rank=T)
bootF <- matrix(NA, nrow=deals, ncol=3) #preallocate
for(i in 1:deals){
    bootX <- sample(X, replace=T)
    bootF[i,] <- Fstat.2way(bootX, f1, f2, center=mean, agg=ssq, rank=T)
}
pval.rankboot <- colSums(sweep(bootF, 2, F0,"-") > 0) / deals
par(mfrow=c(1,3))
xlabel <- c("F diet", "F frog-species", "F interaction")
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
rownames(comparison) <- c("diet", "species", "interaction")
colnames(comparison) <- c("built-in ANOVA", "bootstrap with L1-norm", "bootstrap with L2-norm", "rank", "bootstrap rank")
table <- t(comparison)
library(knitr)
kable(table, digits=4)
