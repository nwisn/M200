# make distribution figures

require(fGarch)
require(ggplot2)

temp1.1 <- rnorm(10000, mean = 0, sd = 1)
temp1.2 <- rnorm(10000, mean = 0, sd = 1)
temp2.1 <- rnorm(10000, mean = 0, sd = 1)
temp2.2 <- rnorm(10000, mean = 0, sd = 4)
temp3.1 <- rsnorm(10000, mean = 0, sd = 1, xi = -1000)
temp3.2 <- rsnorm(10000, mean = 0, sd = 1, xi = 1000)

textsize <- 20
tempdata <- data.frame(values = c(temp1.1, temp1.2), group = c(rep("X", 10000), rep("Y", 10000)))
p1 <- ggplot(tempdata, aes(group, values)) + geom_violin(aes(fill = group), adjust=3) + geom_boxplot(width=0.1) + theme_classic()
p1 <- p1 + theme(axis.text=element_text(size=textsize), 
         axis.title=element_text(size=textsize,face="bold")) + xlab("a") + theme(legend.position="none")

tempdata2 <- data.frame(values = c(temp2.1, temp2.2), group = c(rep("X", 10000), rep("Y", 10000)))
p2 <- ggplot(tempdata2, aes(group, values)) + geom_violin(aes(fill = group), adjust=3) + geom_boxplot(width=0.1) + theme_classic()
p2 <- p2 + theme(axis.text=element_text(size=textsize),axis.title.y = element_blank(),
          axis.title=element_text(size=textsize,face="bold")) + xlab("b") + theme(legend.position="none")

tempdata3 <- data.frame(values = c(temp3.1, temp3.2), group = c(rep("X", 10000), rep("Y", 10000)))
p3 <- ggplot(tempdata3, aes(group, values)) + geom_violin(aes(fill = group), adjust=3) + geom_boxplot(width=0.1) + theme_classic()
p3 <- p3 + theme(axis.text=element_text(size=textsize), axis.title.y = element_blank(),
           axis.title=element_text(size=textsize,face="bold")) + xlab("c") + theme(legend.position="none")

grid.arrange(p1,p2,p3,ncol=3)



boot.samples = matrix(sample(CommuteAtlanta$Time, size = B * n, replace = TRUE), B, n)

set <- matrix(c(150,10))
fisher.test(set)
