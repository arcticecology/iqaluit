library(riojaPlot)
library(readxl)
library(dplyr)
library(rioja)
library(analogue)
library(vegan)
library(ggpalaeo)
library(tidyverse)
library(mgcv)

#
#Environmental change observed from two urban Arctic lakes in Iqaluit, Nunavut
#Nishikawa C., Medeiros A.S., Eamer, J., and Quinlan R.
#R Script for production of stratigraphic plots
#See Fortin et al. 2015 for reconstruction model
#
#
#setwd("C:\\Users\\andre\\OneDrive\\Desktop\\2024\\Connor\\IQ04 submission\\")

coreinput<-read.csv(file="IQ04.csv", row.names=1)

BCI<-coreinput

chron <- cbind(BCI[,1] , BCI[,2])
colnames(chron) <- c("Depth","Year")
chron=as.data.frame(chron)

BCI<-BCI[,-cbind(1:2)] #first two columns
core <- BCI / rowSums(BCI) * 100

# Fortin et al. 2015 used filtered >=2%-in-2-lakes spp criteria in their analyses
# this code is to retain only taxa that are >=2% in 2 intervals in the training set
N <- 2
M <- 2
i <- colSums(core >= M) >= N
spp_red <- core[, i, drop = FALSE]
# 61 chironomid taxa remain in training set
ncol(spp_red)

fos<-spp_red[,-cbind(22:23)] #remove undifferentiated from plots

nms <- colnames(fos)
nms3<-gsub("\\.", " ", nms)




depth<-as.numeric(chron$Depth)



#diversity

H<-diversity(core)
H<-as.data.frame(H)
J <- H/log(specnumber(core))
J<-as.data.frame(J)
#DCA and PCA
dca <- decorana(fos, iweigh=1)
sc1 <- scores(dca, display="sites", choices=1:2)
DC1<-as.data.frame(sc1)

H1.pca <- rda(fos)
sc.pca <-scores(H1.pca,display="sites", choices=1:1)
PC1<-as.data.frame(sc.pca)


recon_IQ04 <-read.csv(file="reconstruction_IQ04_2024.csv", row.names=1)
recon_IQ01 <-read.csv(file="reconstruction_IQ01_rev.csv", row.names=1)
metstation <-read.csv(file="en_climate_monthly_2402590.csv", row.names=1)
library(ggplot2)

#ESM5

p4 <- ggplot()+
  ggtitle("Iqaluit Temperature 1930-2010")+ # for the main title
  xlab("Year")+ # for the x axis label
  ylab(expression("Temperature Anomaly ("*~degree*C*")"))+ # for the y axis label
  geom_point(data=recon_IQ04, aes(x=Year, y=Anom, colour="IQ04"), shape=15) + xlim(1930,2010) +
  geom_smooth(data=recon_IQ04, aes(x=Year, y=Anom), method=lm , color="darkgreen", fill="#69b3a2", se=TRUE) +
  geom_point(data=recon_IQ01, aes(x=Year, y=Anom, colour="IQ01"), shape=17) +
  geom_smooth(data=recon_IQ01, aes(x=Year, y=Anom), method=lm , color="red", fill="#FF9999", se=TRUE) +
  geom_point(data=metstation, aes(x=Year, y=Tanon, colour="MET")) +
  geom_smooth(data=metstation, aes(x=Year, y=Tanon), method=lm , color="blue", fill="#99CCFF", se=TRUE) +
  labs(colour="Legend")
  
p4 + theme(legend.position = c(0.1, 0.8))


p5 <- ggplot()+
  ggtitle("Iqaluit Temperature 1900-2010")+ # for the main title
  xlab("Year")+ # for the x axis label
  ylab(expression("Temperature Anomaly ("*~degree*C*")"))+ # for the y axis label
  geom_point(data=recon_IQ04, aes(x=Year, y=Anom, colour="IQ04"), shape=15) + xlim(1930,2010) +
  geom_smooth(data=recon_IQ04, aes(x=Year, y=Anom), method = "gam", formula = y ~ s(x, k = 3), size = 1 , color="darkgreen", fill="#69b3a2", se=TRUE) +
  geom_point(data=recon_IQ01, aes(x=Year, y=Anom, colour="IQ01"), shape=17) +
  geom_smooth(data=recon_IQ01, aes(x=Year, y=Anom), method = "gam", formula = y ~ s(x, k = 3), size = 1 , color="red", fill="#FF9999", se=TRUE) +
  geom_point(data=metstation, aes(x=Year, y=Tanon, colour="MET")) +
  geom_smooth(data=metstation, aes(x=Year, y=Tanon), method = "gam", formula = y ~ s(x, k = 3), size = 1 , color="blue", fill="#99CCFF", se=TRUE) +
  labs(colour="Legend")

p5 + theme(legend.position = c(0.1, 0.8))


re04<-as.data.frame(recon_IQ04$Comp02)
colnames(re04) <- c("Comp02")


err2<-as.data.frame(recon_IQ04$err2)
colnames(err2) <- c("err2")
err1<-as.data.frame(recon_IQ04$err1)
colnames(err1) <- c("err1")


diss <- dist((core/100)^2)

clust <- chclust(diss, method="coniss")
bstick(clust, 5)
plot(clust, hang=-1)

riojaPlot(fos, chron, 
          yvar.name="Depth",
          sec.yvar.name="Year",
          plot.sec.axis=TRUE,
          scale.percent=TRUE, 
          plot.exag=TRUE)

#IQ04 stratigraphy 

tiff("IQ04.tiff", width = 6, height = 4, units = 'in', res = 300)


rp1 <- riojaPlot(fos, chron, 
                 ymin = 0, ymax=20, yinterval = 2, 
                 yvar.name="Depth",
                 sec.yvar.name="Year",
                 plot.sec.axis=TRUE,
                 ylabPos=0.3,
                 y.rev=TRUE,
                 yTop=0.6,
                 yBottom=0.20,
                 scale.percent=TRUE, plot.groups=FALSE, do.clust = TRUE,
                 plot.poly=FALSE, plot.line=FALSE, lwd.bar=3,
                 plot.clust=TRUE,  
                 plot.cumul=TRUE, cex.cumul=0.6,
                 srt.xlabel=45, plot.bar=TRUE, xlabels=nms3,labels.break.long=FALSE,
                 tcl=-0.1, cex.yaxis=0.5, cex.xlabel=0.5, col.bar="black",
                 cex.xaxis=0.5, xRight = 0.7, plot.line=FALSE)

#par(fig=c(0.8,1,0,1), new=TRUE)


tmp2 <- data.frame(x=re04$Comp02, y=chron$Depth)
sdev<-err2$err2
v1<-err1$err1
fun.gam3 <- function(x, y, i, nm, style) {
  lo <- loess(re04$Comp02 ~ chron$Depth)
  lines(lo$fitted, chron$Depth, col = "darkgray", lwd=2)
  arrows(re04$Comp0-sdev, chron$Year, re04$Comp0+sdev, chron$Depth, length=0.05, angle=180, code=3)
  
}
fun.gam2 <- function(x, y, i, nm, style) {
  lo <- loess(re04$Comp02 ~ chron$Depth)
  lines(lo$fitted, chron$Depth, col = "red", lwd=2)
  polygon(c(re04$Comp0+sdev, rev(re04$Comp0-sdev)), c(chron$Depth, rev(chron$Depth)), 
          col = rgb(0.7, 0.7, 0.7,0.5))
  arrows(re04$Comp0-v1, chron$Depth, re04$Comp0+v1, chron$Depth, length=0.05, angle=180, code=3)
}

min2<- data.frame(min=1.7,max=2.7) #scale the x
min3<- data.frame(min=c(-0.5, -0.5),max=c(0.8, 0.8)) #scale the x
rp2 <- riojaPlot(H, chron, yvar.name="Depth",
                 riojaPlot=rp1,xGap = 0.01,
                 xRight=0.76, minmax = min2,
                 scale.minmax=FALSE, plot.bar=F, 
                 plot.line=T, plot.symb=TRUE, symb.cex=0.5)

rp3 <- riojaPlot(DC1, chron, yvar.name="Depth",
                 riojaPlot=rp2,xGap = 0.01, xSpace = 0.01,
                 xRight=0.88, minmax = min3,
                 scale.minmax=FALSE, plot.bar=F, 
                 plot.line=T, plot.symb=TRUE, symb.cex=0.5)

riojaPlot(re04, chron, 
                 riojaPlot=rp3, xGap = 0.001, yvar.name="Depth",
                 xRight=0.955, scale.minmax=FALSE, 
                 plot.bar=FALSE, plot.line=FALSE, 
                 plot.symb=TRUE, symb.cex=0.5, fun.xfront=fun.gam2)

#rp3 <- addRPClustZone(rp3, clust, nZone=2, lwd=2, lty=1, col=rgb(0,0,0,0.5))
rp4<-addRPZone(rp3, 2.8, 1.8, xLeft=NULL, xRight=NULL, col=rgb(0,0,0,0.5), 
               alpha=0.1, border=NA, verbose=TRUE, lty=2, lwd=2)
#rp5 <- addRPClustZone(rp2, clust, nZone=3,  lwd=2, lty=1, col=rgb(0,0,0,0.5))


dev.off()



#
#
# script for IQ01
#
#
coreinput<-read.csv(file="IQ01.csv", row.names=1)

BCI<-coreinput

chron <- cbind(BCI[,1] , BCI[,2])
colnames(chron) <- c("Depth","Year")
chron=as.data.frame(chron)

BCI<-BCI[,-cbind(1:2)] #first two columns
core <- BCI / rowSums(BCI) * 100

# Fortin et al. 2015 used filtered >=2%-in-2-lakes spp criteria in their analyses
# this code is to retain only taxa that are >=2% in 2 intervals in the training set
N <- 2
M <- 2
i <- colSums(core >= M) >= N
spp_red <- core[, i, drop = FALSE]
# 61 chironomid taxa remain in training set
ncol(spp_red)

fos<-spp_red[,-cbind(24:25)] #remove undifferentiated from plots


nms <- colnames(fos)
nms3<-gsub("\\.", " ", nms)




depth<-as.numeric(chron$Depth)


diss <- dist((core/100)^2)

clust <- chclust(diss, method="coniss")
bstick(clust, 5)
plot(clust, hang=-1)
#diversity

H<-diversity(core)
H<-as.data.frame(H)
J <- H/log(specnumber(core))
J<-as.data.frame(J)
#DCA and PCA
dca <- decorana(fos, iweigh=1)
sc1.dca <- scores(dca, display="sites", choices=1:2)
DC1<-as.data.frame(sc1.dca)

H1.pca <- rda(fos)
sc.pca <-scores(H1.pca,display="sites", choices=1:1)
PC1<-as.data.frame(sc.pca)

recon_IQ01 <-read.csv(file="reconstruction_IQ01_rev.csv", row.names=1)
re04<-as.data.frame(recon_IQ01$Comp02)
colnames(re04) <- c("Comp02")

err1<-as.data.frame(recon_IQ01$err1)
colnames(err1) <- c("err1")

err2<-as.data.frame(recon_IQ01$err2)
colnames(err2) <- c("err2")


library(mgcv)



fun.gam3a <- function(x, y, i, nm, style) {
  tmp <- data.frame(x=y, y=x)
  lw1 <- loess(y ~ x,data=tmp)
  j <- order(data$y)
  lines(tmp$x[j],lw1$fitted[j],col="red",lwd=3)
}


riojaPlot(fos, chron, 
          yvar.name="Depth",
          sec.yvar.name="Year",
          plot.sec.axis=TRUE,
          scale.percent=TRUE, 
          plot.exag=TRUE)


min1<- data.frame(min=7,max=9) #scle the chironomid reconstruction
min5<- data.frame(min=2,max=3) #scale the x

#IQ01  stratigraphy

tiff("IQ01.tiff", width = 6, height = 4, units = 'in', res = 300)
rp1b <- riojaPlot(fos, chron, 
                 ymax = 25, 
                 ymin=0, yinterval = 2, 
                 yvar.name = "Depth",
                 sec.yvar.name="Year",
                 plot.sec.axis=TRUE,
                 y.rev=TRUE,
                 yTop=0.6,
                 yBottom=0.20, xlabels=nms3,
                 scale.percent=TRUE, plot.groups=FALSE, do.clust = TRUE,
                 plot.poly=FALSE, plot.line=FALSE, lwd.bar=3,
                 plot.clust=TRUE,   ylabPos=0.3,
                 plot.cumul=TRUE, cex.cumul=0.6,
                 srt.xlabel=45, plot.bar=TRUE,labels.break.long=FALSE,
                 tcl=-0.1, cex.yaxis=0.5, cex.xlabel=0.5, col.bar="black",
                 cex.xaxis=0.5, xRight = 0.7, plot.line=FALSE)


tmp2 <- data.frame(x=re04$Comp02, y=chron$Year)


sdev<-err2$err2
v1<-err1$err1

ave<-re04$Comp02
fun.gam2 <- function(x, y, i, nm, style) {
  tmp <- data.frame(x=y, y=x)
  x <- sort(tmp2$x)
  y<-tmp2$y[order(tmp2$y)]
  loess_fit <- loess(y~x)
   lines(x, predict(loess_fit), col = "red", lwd=2)
  arrows(ave-sdev, derp, ave+sdev, derp, length=0.05, angle=180, code=3)
  
}

fun.gam3 <- function(x, y, i, nm, style) {
  lo <- loess(re04$Comp02 ~ chron$Depth)
  lines(lo$fitted, chron$Depth, col = "darkgray", lwd=2)
  arrows(re04$Comp0-sdev, chron$Depth, re04$Comp0+sdev, chron$Depth, length=0.05, angle=180, code=3)
}

fun.gam3b <- function(x, y, i, nm, style) {
  lo <- loess(re04$Comp02 ~ chron$Depth)
  lines(lo$fitted, chron$Depth, col = "red", lwd=2)
  polygon(c(re04$Comp0+sdev, rev(re04$Comp0-sdev)), c(chron$Depth, rev(chron$Depth)), 
          col = rgb(0.7, 0.7, 0.7,0.5))
  arrows(re04$Comp0-v1, chron$Depth, re04$Comp0+v1, chron$Depth, length=0.05, angle=180, code=3)
}

rp2b <- riojaPlot(H, chron, yvar.name="Depth",
                 riojaPlot=rp1b,xGap = 0.009,minmax=min5,
                 xRight=0.76, 
                 scale.minmax=FALSE, plot.bar=F, 
                 plot.line=T, plot.symb=TRUE, symb.cex=0.6)

rp3b <- riojaPlot(DC1, chron, yvar.name="Depth",
                 riojaPlot=rp2b,xGap = 0.01, xSpace=0.01,
                 xRight=0.89, minmax=min3,
                 scale.minmax=FALSE, plot.bar=F, 
                 plot.line=T, plot.symb=TRUE, symb.cex=0.6)

riojaPlot(re04, chron, 
          riojaPlot=rp3b, xGap = 0.01, yvar.name="Depth",
          xRight=0.955, scale.minmax=FALSE, 
          plot.bar=FALSE, plot.line=FALSE, 
          plot.symb=TRUE, symb.cex=0.6,minmax=min1, fun.xfront=fun.gam3b)

rp3 <- addRPClustZone(rp3b, clust, nZone=2,  lwd=2, lty=1, col=rgb(0,0,0,0.5))


dev.off()

