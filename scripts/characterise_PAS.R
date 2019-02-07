source("http://bioconductor.org/biocLite.R")
require("GO.db")
require("limma")
require("topGO")
require("biomaRt")
require("org.Rn.eg.db")
library(grDevices)
library(Rsamtools)
library(IRanges)
library(GenomicRanges)
library(rtracklayer)
library(geneplotter)
require("multtest")
require("mclust")
rm(list = ls())

source("./scripts/nested_functions.R")


load("./annotation/rn5/GRanges_comprehensive_transcriptome_rat_24_nov_2015.RData")
myUTR        <- import.gff("./annotation/rn5/Lngf_sub.gtf",format="gtf")
anno_ngf     <- read.table("./data/anno_ngf.tab",header=T,sep="\t")
names(myUTR) <- as.character(myUTR$ID)
ix1          <- match(names(myUTR),anno_ngf$uniqueID) 
myUTR        <- myUTR[!is.na(ix1),]
ix1          <- match(names(myUTR),anno_ngf$uniqueID) 
anno_ngf     <- anno_ngf[ix1,]
myUTR           <- myUTR[anno_ngf$newL>=10,]
anno_ngf        <- anno_ngf[anno_ngf$newL>=10,]

myGRpaf                 <- import.gff(con="./annotation/rn5/myPAS.gtf",format="gtf")
motifs                  <- unique(as.character(myGRpaf$rep.motifID..length.end..))
myPA                    <- lapply(motifs,function(x)return(myGRpaf[myGRpaf$rep.motifID..length.end..==x,]))
names(myPA)             <- motifs
myPA[[length(myPA)+1]]  <- myGRpaf
names(myPA)[13]         <-"all"


# A. Create fasta sequence of 3' UTR w/o introns BUT with additional 100 nt at the end
Code in   : /dependencies/extractFASTA.R
Stored in : /home/rluisier/data/Riccio/Exp_1/Dec2016/utrid/APA/L3_extended_no_introns.fasta

# B. Get motif in sequence
Code in   : /dependencies/findPASinUTR.R
Stored in : /home/rluisier/data/Riccio/Exp_1/Dec2016/utrid/myRes_APA_12032016.RData


# C. Characterise their content
load("/home/rluisier/data/Riccio/Exp_1/Dec2016/utrid/myRes_APA_12032016.RData")
load("/Users/luisier/Desktop/DataAnalsyisRiccio/Dec2016/myRes_APA_12032016.RData")
load("/Users/luisier/Desktop/DataAnalysisRiccio/Dec2016/myRes_APA_12032016.RData")
#Be careful as not necessary the same order


# F1. Focus [-100,+50] -- for paper; distribution no. of PAS
require(gplots)
GetNoPAS <- function(target=myPA[[1]],utrGR=focusGR){
  # o) Remove redundant targets
  target                <- target[!duplicated(target)]
  # a) Find overlap
  myHits               <- subjectHits(findOverlaps(query=target,subject=utrGR,ignore.strand=F))
  # b) Find how many Hits in region of interest
  no.hits              <- as.data.frame(table(myHits))#get number of sites per isoform
  out                  <- rep(0,length(utrGR))
  names(out)           <- utrGR$ID
  out[as.numeric(as.character(no.hits$myHits))] <- as.numeric(as.character(no.hits$Freq))
  return(out)
}

GetFracWithSites <- function(my.nos=my.no.sites[,1],sel=rep(TRUE,length(anno_ngf$newL))){
  cons.frac <- sum(my.nos[sel&anno_ngf$is.conservative]>0)/sum(anno_ngf$is.conservative[sel])
  cons.new  <- sum(my.nos[sel&!anno_ngf$is.conservative]>0)/sum(!anno_ngf$is.conservative[sel])
  return(c(cons.frac,cons.new))
}

IsPAdb <- function(target=myPAdb,utrGR=focusGR){
  # o) Remove redundant targets
  target                <- target[!duplicated(target)]
  # a) Find overlap
  myHits               <- subjectHits(findOverlaps(query=target,subject=utrGR,ignore.strand=F))
  # b) Find how many Hits in region of interest
  no.hits              <- as.data.frame(table(myHits))#get number of sites per isoform
  out                  <- rep(0,length(utrGR))
  names(out)           <- utrGR$ID
  out[as.numeric(as.character(no.hits$myHits))] <- as.numeric(as.character(no.hits$Freq))
  return(out)
}

GetFracWithSitesLength <- function(my.nos=my.no.sites[,1],BINS=BINS1,sel=rep(TRUE,length(anno_ngf$newL))){
  my.nos <- my.nos[sel]
  myidx  <- factor(as.character(BINS[sel]))
  return(tapply(my.nos,INDEX=myidx,FUN=function(Z)return(sum(Z>0)/length(Z))))
}

TestFrac1 <- function(n1,n2,p1,p2){
  p  <- (p1*n1+p2*n2)/(n1+n2)
  SE <- sqrt(p*(1-p)*(1/n1+1/n2))
  Z  <- (p1-p2)/SE
  PVal<- 2*pnorm(-abs(Z))#2sided
  return(list(Z,PVal))
}

TestFrac2 <- function(n1,n2,p1,p2){
  SE <- sqrt((p1*(1-p1)/n1)+(p2*(1-p2)/n2))
  Z  <- (p1-p2)/SE
  PVal<- 2*pnorm(-abs(Z))#2sided
  return(list(Z,PVal))
}

myno1 <- 100*c(80,120,170)/689
myno2 <- 100*c(50,70,100)/737
pdf("~/Desktop/tesmp.pdf")
par(mfrow=c(1,2))
mp<-barplot(myno1,las=1,frame=F)
mtext(side=3,line=0,text=round(myno1,digit=1),at=mp)
mp<-barplot(myno2,las=1,frame=F)
mtext(side=3,line=0,text=round(myno2,digit=1),at=mp)
dev.off()

#BINS0               <- cut(log10(anno_ngf$newL),breaks=seq(1,4.5,0.25),include.lowest=T)
#BINS1               <- cut(log10(anno_ngf$newL),breaks=c(1,1.5,2,2.25,2.75,3.25,3.75,4.25),include.lowest=T)
#BINS2                <- cut(log10(anno_ngf$newL),breaks=c(1.5,2,2.25,2.5,2.75,3.0,3.25,3.5,3.75,4.0,4.25),include.lowest=T)
BINS2                <- cut(log10(anno_ngf$newL),breaks=c(1.5,2,2.25,2.75,3.25,3.75,4.25),include.lowest=T)
BINS2                <- cut(log10(anno_ngf$newL),breaks=seq(from=1.75,t=4.25,by=0.25),include.lowest=T)

focusGR            <- myUTR
focusGR            <- resize(x=focusGR,fix="start",width=width(focusGR)+25,use.names=T)
focusGR            <- resize(x=focusGR,fix="end",width=100,use.names=T)

focusGR            <- myUTR
focusGR            <- resize(x=focusGR,fix="start",width=width(focusGR)+50,use.names=T)
focusGR            <- resize(x=focusGR,fix="end",width=150,use.names=T)

focusAnno          <- anno_ngf[match(names(focusGR),anno_ngf$uniqueID),]
my.no.sites        <- do.call(what=cbind,args=lapply(myPA,function(X)return(GetNoPAS(target=X,utrGR=focusGR))))

my.no.sites.sub    <- cbind(my.no.sites[,c("AATAAA","ATTAAA")],apply(my.no.sites[,-match(c("AATAAA","ATTAAA","all"),colnames(my.no.sites))],1,sum))
colnames(my.no.sites.sub)[3] <- "remaining" 
my.no.padb      <- IsPAdb(target=myPAdb,utrGR=focusGR)
no.padb       <- c(GetFracWithSites(my.nos=my.no.padb,sel=rep(TRUE,length(my.no.padb))),
                   GetFracWithSitesLength(my.nos=my.no.padb,BINS=BINS2,sel=!focusAnno$is.conservative))


my.frac.sites.1    <- apply(my.no.sites,2,function(Z)return(GetFracWithSites(my.nos=Z,sel=rep(TRUE,length(Z)))))
#my.frac.sites.2   <- apply(my.no.sites,2,function(Z)return(GetFracWithSites(my.nos=Z,sel=focusAnno$is.conservative|myUTR$is.pas=="TRUE")))
my.frac.sites.2    <- apply(my.no.sites,2,function(Z)return(GetFracWithSites(my.nos=Z,sel=myUTR$is.pas=="TRUE")))#this one is more logic given that we want to understand for all the sites with a PAS, which one are present most of the time; we evidence on those we see and simply want to understand whether the frequency is the same
my.frac.sites.3    <- apply(my.no.sites.sub,2,function(Z)return(GetFracWithSites(my.nos=Z,sel=myUTR$is.pas=="TRUE")))

my.frac.sites.3    <- cbind(my.frac.sites.3,1-my.frac.sites.1[,13])

pdf("/Users/luisier/Desktop/DataAnalsyisRiccio/Dec2016/utrid/pas_stringent/characterise_PAS.pdf")
pdf("/Users/luisier/Desktop/DataAnalysisRiccio/Dec2016/utrid/pas_stringent/characterise_PAS_150_sub.pdf")

pdf("/Users/luisier/Desktop/DataAnalysisRiccio/Dec2016/utrid/pas_stringent/characterise_PAS_75.pdf")

par(mfrow=c(2,3))
mp=barplot(my.frac.sites.1[,13]*100,col="white",las=1,cex.axis=0.5,ylim=c(0,100))
mtext(side=2,line=3,text="% of 3' UTR with PAS in ]-100:+50[ nt around 3' end",cex=0.5)
grid()
mtext(side=3,line=0,text=round(my.frac.sites.1[,13]*100,digits=2),at=mp,cex=0.5)
barplot(my.frac.sites.2[,-13]*100,beside=T,col=c("white","black"),las=1,cex.axis=0.5,cex=0.5,ylim=c(0,100))
mtext(side=2,line=3,text="% of 3' UTR with PAS in ]-100:+50[ nt around 3' end with at least one PAS",cex=0.5)
barplot(my.frac.sites.3*100,beside=T,col=c("white","black"),las=1,cex.axis=0.5,cex=0.5,ylim=c(0,100))
mtext(side=2,line=3,text="% of 3' UTR with PAS in ]-100:+50[ nt around 3' end with at least one PAS",cex=0.5)
grid()
#barplot(no.padb[-2]*100,beside=T,col=c("white",rep("black",7)),las=1,cex.axis=0.5)
#mtext(side=2,line=3,text="% of 3' UTR with PAdb in ]-100:+50[ nt around 3' end",cex=0.5)

#Compare;perform statistical test to show significant difference in PAS usage between conservative and novel
#https://onlinecourses.science.psu.edu/stat414/node/268
#The test statistic for testing the difference in two population proportions, that is, for testing the null hypothesis p1-p2=0:
TestFrac1(n1=sum(myUTR$is.pas=="TRUE"&anno_ngf$is.conservative),n2=sum(myUTR$is.pas=="TRUE"&!anno_ngf$is.conservative),p1=my.frac.sites.3[1,1],p2=my.frac.sites.3[2,1])
TestFrac2(n1=sum(myUTR$is.pas=="TRUE"&anno_ngf$is.conservative),n2=sum(myUTR$is.pas=="TRUE"&!anno_ngf$is.conservative),p1=my.frac.sites.3[1,1],p2=my.frac.sites.3[2,1])

TestFrac1(n1=sum(myUTR$is.pas=="TRUE"&anno_ngf$is.conservative),n2=sum(myUTR$is.pas=="TRUE"&!anno_ngf$is.conservative),p1=my.frac.sites.3[1,2],p2=my.frac.sites.3[2,2])
TestFrac2(n1=sum(myUTR$is.pas=="TRUE"&anno_ngf$is.conservative),n2=sum(myUTR$is.pas=="TRUE"&!anno_ngf$is.conservative),p1=my.frac.sites.3[1,2],p2=my.frac.sites.3[2,2])

TestFrac1(n1=sum(myUTR$is.pas=="TRUE"&anno_ngf$is.conservative),n2=sum(myUTR$is.pas=="TRUE"&!anno_ngf$is.conservative),p1=my.frac.sites.3[1,3],p2=my.frac.sites.3[2,3])
TestFrac2(n1=sum(myUTR$is.pas=="TRUE"&anno_ngf$is.conservative),n2=sum(myUTR$is.pas=="TRUE"&!anno_ngf$is.conservative),p1=my.frac.sites.3[1,3],p2=my.frac.sites.3[2,3])


myLev<-c("[1.75,2]","(2,2.25]","(2.25,2.5]","(2.5,2.75]","(2.75,3]","(3,3.25]","(3.25,3.5]","(3.5,3.75]","(3.75,4]","(4,4.25]")

avg.length <- apply(my.no.sites[myUTR$is.pas=="TRUE",],2,function(Z)return(tapply(Z,INDEX=BINS2[myUTR$is.pas=="TRUE"],FUN=mean)))
avg.length.cons<-apply(my.no.sites[anno_ngf$is.conservative&myUTR$is.pas=="TRUE",],2,function(Z)return(tapply(Z,INDEX=factor(as.character(BINS2[anno_ngf$is.conservative&myUTR$is.pas=="TRUE"]),levels=myLev),FUN=mean)))
avg.length.new<-apply(my.no.sites[myUTR$is.pas=="TRUE"&!anno_ngf$is.conservative,],2,function(Z)return(tapply(Z,INDEX=factor(as.character(BINS2[myUTR$is.pas=="TRUE"&!anno_ngf$is.conservative]),levels=myLev),FUN=mean)))

par(mfrow=c(2,2))
barplot(avg.length[,13],ylim=c(0,2))
grid()
abline(h=mean(my.no.sites[myUTR$is.pas=="TRUE",13]),col="red")
barplot(avg.length.cons[,13],ylim=c(0,2))
barplot(avg.length.new[,13],ylim=c(0,2))

dev.off()

#This version is potentially limited by the size of the group of interest
EnrichmentAnalysisNumber <- function(dat=my.no.sites[,1],no.boot=500,myFocus=BINS1==levels(BINS1)[3],sel=anno_ngf$is.conservative){
  
  dat         <- dat[sel]
  myFocus     <- myFocus[sel] 
  dat         <- dat[!is.na(myFocus)]
  myFocus     <- myFocus[!is.na(myFocus)]
  
  #Foreground: how many sites in this range
  Ntot    <- sum(dat[myFocus])
  #Background
  myBG     <- unlist(lapply(c(1:no.boot),function(Z)return(sum(sample(dat,size=sum(myFocus),replace=TRUE)))))
  #Compute Z-score of enrichment
  myZ      <- (Ntot-mean(myBG))/sd(myBG)
  p.greater<- sum(myBG<Ntot)/no.boot
  p.smaller<- sum(myBG>Ntot)/no.boot
  return(c(myZ,p.greater,p.smaller))
}


EnrichmentAnalysisAvg <- function(dat=my.no.sites[,1],no.boot=500,myFocus=BINS1==levels(BINS1)[3],sel=anno_ngf$is.conservative){
  
  dat     <- dat[sel]
  myFocus <- myFocus[sel] 
  
  #Foreground: how many sites in this range
  Avg    <- mean(dat[myFocus])
  Sd     <- sd(dat[myFocus])
  
  #Background
  myBG     <- unlist(lapply(c(1:no.boot),function(Z)return(mean(sample(dat,size=sum(myFocus),replace=TRUE)))))
  #Compute Z-score of enrichment
  myZ      <- (Avg-mean(myBG))/sqrt(sd(myBG)^2+Sd^2)
  p.greater<- sum(myBG<Avg)/no.boot
  p.smaller<- sum(myBG>Avg)/no.boot
  return(c(myZ,p.greater,p.smaller))
}



EnrichmentAnalysisHT <- function(dat=my.no.sites[,1],myFocus=BINS2==levels(BINS2)[3],sel=myUTR$is.pas=="TRUE"){
  
  dat         <- dat[sel]>0
  myFocus     <- myFocus[sel] 
  dat         <- dat[!is.na(myFocus)]
  myFocus     <- myFocus[!is.na(myFocus)]
  
  A       <- sum(dat&myFocus)
  B       <- sum(dat)-A
  C       <- sum(myFocus)-A
  D       <- length(dat)-A-B-C
  
  return(c(fisher.test(matrix(c(A,B,C,D),nrow=2,ncol=2),alternative="greater")$p.value,
           fisher.test(matrix(c(A,B,C,D),nrow=2,ncol=2),alternative="less")$p.value))
  
}



plotHeatmap <- function(mat=t(myZval.per.pos[,-13]),mytitle){
  require(colorRamps)
  require(gplots)
  
  mycols          <- c(colorpanel(n=50, low="cyan",high="black"),
                       colorpanel(n=50, low="black",high="magenta"))
  
  mycols          <- c(colorpanel(n=50, low="cyan",high="white"),
                       colorpanel(n=50, low="white",high="magenta"))
  
  b               <- seq(-max(abs(mat)),max(abs(mat)),length=101)
  hc2             <- hclust(dist(mat,"man"), method="ward.D")
  heatmap.2(mat,keysize=1,mar=c(10,10),col=mycols,breaks=b,
            scale="none",
            dendro="row",Rowv=as.dendrogram(hc2),Colv=FALSE,
            key=TRUE,symkey=TRUE, density.info="none", trace="none",
            cexCol=0.8, cexRow=0.8, font.lab=1,main=mytitle)
}


#BINS2                <- cut(log10(anno_ngf$newL),breaks=c(1.5,2,2.25,2.75,3.25,3.75,4.25),include.lowest=T)
myBINS=BINS2
myZval.per.pos <- do.call(lapply(c(1:13),function(W)return(do.call(lapply(levels(myBINS),function(Z)return(EnrichmentAnalysisNumber(dat=my.no.sites[,W],no.boot=500,myFocus=myBINS==Z,sel=myUTR$is.pas=="TRUE")[1])),what=c))),what=cbind)
colnames(myZval.per.pos)<-colnames(my.no.sites)
rownames(myZval.per.pos)<-levels(myBINS)


#Here we want to characterise the strenght of 3' UTR in function of the number of PAS present in the last 150 nt; you can either do this using the sum or using the average. With the average I am not sure about the stat.
#In this case I do not really care about which PAS, but I am rather intersted in the total number of PAS present irrespective of PAS sequence.
mySum    <- tapply(my.no.sites[my.no.sites[,13]>0,13],INDEX=myBINS[my.no.sites[,13]>0],FUN=sum)
myMean   <- tapply(my.no.sites[my.no.sites[,13]>0,13],INDEX=myBINS[my.no.sites[,13]>0],FUN=mean)
mySD     <- tapply(my.no.sites[my.no.sites[,13]>0,13],INDEX=myBINS[my.no.sites[,13]>0],FUN=sd)
#Z-score must be computed as if this was Poisson distributed and not like Gaussian
myZscore <- (myMean-mean(my.no.sites[my.no.sites[,13]>0,13]))/sqrt(mySD^2+sd(my.no.sites[my.no.sites[,13]>0,13])^2) 

myMed   <- tapply(my.no.sites[my.no.sites[,13]>0,13],INDEX=myBINS[my.no.sites[,13]>0],FUN=median)
myMad   <- tapply(my.no.sites[my.no.sites[,13]>0,13],INDEX=myBINS[my.no.sites[,13]>0],FUN=mad)
myZ2   <- (myMad-median(my.no.sites[my.no.sites[,13]>0,13]))/sqrt(myMad^2+mad(my.no.sites[my.no.sites[,13]>0,13])^2) 

#Here we study the frequency of site with specific PAS in order to get an idead of the preferential usage of different PAS according to 3' UTR length (please note this is not necessarily associated with the strenght of 3' UTR)
#Fiisher count test
myBINS=BINS2
myZPVal.per.pos.more <- -log10(do.call(lapply(c(1:13),function(W)return(do.call(lapply(levels(myBINS),function(Z)return(EnrichmentAnalysisHT(dat=my.no.sites[,W],myFocus=myBINS==Z,sel=myUTR$is.pas=="TRUE")[1])),what=c))),what=cbind))
myZPVal.per.pos.less <- -log10(do.call(lapply(c(1:13),function(W)return(do.call(lapply(levels(myBINS),function(Z)return(EnrichmentAnalysisHT(dat=my.no.sites[,W],myFocus=myBINS==Z,sel=myUTR$is.pas=="TRUE")[2])),what=c))),what=cbind))
myPval <- myZPVal.per.pos.more
myPval[abs(myZPVal.per.pos.less)>abs(myZPVal.per.pos.more)] <- -myZPVal.per.pos.less[abs(myZPVal.per.pos.less)>abs(myZPVal.per.pos.more)]
colnames(myPval)<-colnames(my.no.sites)
rownames(myPval)<-levels(myBINS)


pdf("~/Desktop/DataAnalysisRiccio/Dec2016/utrid/pas_stringent/prefential_positioning.pdf")

plotHeatmap(mat=t(myPval[,-13]),mytitle="P-Val enrichment")
par(mfrow=c(3,2))
for(i in c(1:13)){
  barplot(tapply(my.no.sites[my.no.sites[,13]>0,i],INDEX=myBINS[my.no.sites[,13]>0],FUN=function(X)return(sum(X>0)/length(X))),main=colnames(my.no.sites)[i])
}

plotHeatmap(mat=t(myZval.per.pos[,-13]),mytitle="Z-score enrichment")



par(mfrow=c(2,2))
barplot(myZval.per.pos[,13],las=1,ylab="Z-score for number sites")
barplot(myMean,las=1,ylab="average number sites")
barplot(myZscore,las=1,ylab="Z-score for increase number of sites")
barplot(table(myBINS))

myMean.cons <- tapply(my.no.sites[my.no.sites[,13]>0&anno_ngf$is.conservative,13],INDEX=myBINS[my.no.sites[,13]>0&anno_ngf$is.conservative],FUN=mean)
myMean.new  <- tapply(my.no.sites[my.no.sites[,13]>0&!anno_ngf$is.conservative,13],INDEX=myBINS[my.no.sites[,13]>0&!anno_ngf$is.conservative],FUN=mean)
mySD.cons   <- tapply(my.no.sites[my.no.sites[,13]>0&anno_ngf$is.conservative,13],INDEX=myBINS[my.no.sites[,13]>0&anno_ngf$is.conservative],FUN=sd)
mySD.new    <- tapply(my.no.sites[my.no.sites[,13]>0&!anno_ngf$is.conservative,13],INDEX=myBINS[my.no.sites[,13]>0&!anno_ngf$is.conservative],FUN=sd)
myZ.cons    <- myMean.cons/mySD.cons
myZ.new    <- myMean.new/mySD.new

myZscore.cons <- (myMean.cons-mean(my.no.sites[my.no.sites[,13]>0&anno_ngf$is.conservative,13]))/sqrt(mySD.cons^2+sd(my.no.sites[my.no.sites[,13]>0&anno_ngf$is.conservative,13])^2) 
myZscore.new <- (myMean.new-mean(my.no.sites[my.no.sites[,13]>0&!anno_ngf$is.conservative,13]))/sqrt(mySD.new^2+sd(my.no.sites[my.no.sites[,13]>0&!anno_ngf$is.conservative,13])^2) 

barplot(myMean.cons,ylab="avg no. sites cons")
barplot(myMean.new,ylab="avg no. sites new")

barplot(myZscore.cons,ylab="Z-score cons")
barplot(myZscore.new,ylab="Z-score new")
dev.off()

#Test significance of the enrichmetn using either Poisson or Negative Binomial
dat         <-  data.frame(no.pas=my.no.sites[my.no.sites[,13]>0,13],length=BINS2[my.no.sites[,13]>0],tL=log10(anno_ngf$newL)[my.no.sites[,13]>0])
require(ggplot2)
ggplot(dat, aes(no.pas, fill = length)) +
  geom_histogram(binwidth=0.5) +
  facet_grid(length ~ ., margins=TRUE, scales="free")

with(dat, tapply(no.pas, length, function(x) {
  sprintf("M (SD) = %1.2f (%1.2f)", mean(x), sd(x))
}))

#Variance equals means
require(MASS)
summary(m1 <- glm.nb(no.pas ~length, data = dat))
summary(m1 <- glm.nb(no.pas ~tL, data = dat))

barplot(summary(m1 <- glm.nb(no.pas ~0+length, data = dat))$coeff[,3])

m2 <- glm(no.pas ~0+length, family = "poisson", data = dat)
barplot(summary(m2)$coeff[,3])


Coefficients:
  Estimate Std. Error z value Pr(>|z|)    
(Intercept)       -0.38349    0.07732  -4.960 7.05e-07 ***
  length(1.5,2]      0.50391    0.08094   6.226 4.78e-10 ***
  length(2,2.25]     0.66723    0.07957   8.385  < 2e-16 ***
  length(2.25,2.75]  0.72018    0.07785   9.250  < 2e-16 ***
  length(2.75,3.25]  0.75126    0.07768   9.671  < 2e-16 ***
  length(3.25,3.75]  0.83407    0.07786  10.713  < 2e-16 ***
  length(3.75,4.25]  0.79916    0.07987  10.006  < 2e-16 ***
  
  
  #Check that the model is correct and not Poisson  
  #negative binomial models assume the conditional means are not equal to the conditional variances. This inequality is captured by estimating a dispersion parameter (not shown in the output) that is held constant in a Poisson model. Thus, the Poisson model is actually nested in the negative binomial model. We can then use a likelihood ratio test to compare these two and test this model assumption.
  m1 <- glm.nb(no.pas ~length, data = dat)
  m2 <- glm(no.pas ~length, family = "poisson", data = dat)
  m2 <- glm(no.pas ~tL, family = "poisson", data = dat)
  
  
  X2 <- 2 * (logLik(m1) - logLik(m2))
  X2
  'log Lik.' 907.6986 (df=8)
  pchisq(X2, df = 1, lower.tail=FALSE)
  'log Lik.' 2.16e-203 (df=5)
  #This very large chi-square strongly suggests the negative binomial model, which estimates the dispersion parameter, is more appropriate than the Poisson model.
  (est <- cbind(Estimate = coef(m1), confint(m1)))
  exp(est)
  
  newdata2 <- data.frame(tL = seq(from = min(dat$tL), to = max(dat$tL), length.out = 100))
  newdata2 <- cbind(newdata2, predict(m2, newdata2, type = "link", se.fit=TRUE))
  newdata2 <- within(newdata2, {
    no.pas <- exp(fit)
    LL <- exp(fit - 1.96 * se.fit)
    UL <- exp(fit + 1.96 * se.fit)
  })
  
  ggplot(newdata2,aes(tL, no.pas))+
    geom_ribbon(aes(ymin = LL, ymax = UL), alpha = .25) +
    geom_line( size = 2) +
    labs(x = "3' UTR length", y = "Predicted # PAS")
  
  
  
  
  #Version before July
  BINS0               <- cut(log10(anno_ngf$newL),breaks=seq(1,4.5,0.25),include.lowest=T)
  BINS1               <- cut(log10(anno_ngf$newL),breaks=c(1,1.5,2,2.25,2.75,3.25,3.75,4.25),include.lowest=T)
  BINS2               <- cut(log10(anno_ngf$newL),breaks=c(1,1.5,2,2.25,2.5,2.75,3.0,3.25,3.5,3.75,4.0,4.25),include.lowest=T)
  
  pdf("/Users/luisier/Desktop/DataAnalsyisRiccio/Dec2016/utrid/pas_stringent/distribution_pas_along_sites.pdf")
  par(mfrow=c(2,2))
  barplot(my.frac.sites.1[,13]*100,col="white",las=1,cex.axis=0.5)
  mtext(side=2,line=3,text="% of 3' UTR with PAS in ]-100:+50[ nt around 3' end",cex=0.5)
  barplot(my.frac.sites.2[,-13]*100,beside=T,col=c("white","black"),las=1,cex.axis=0.5,cex=0.5)
  barplot(no.padb.2[-2]*100,beside=T,col=c("white",rep("black",6)),las=1,cex.axis=0.5)
  mtext(side=2,line=3,text="% of 3' UTR with PAdb in ]-100:+50[ nt around 3' end",cex=0.5)
  
  #Add 3'UTR length extended versus normal; extended only
  is.extended                   <- unique(anno_ngf$txID[(anno_ngf$newL-anno_ngf$initL)>10])
  myL                           <- list(
    Lconservative=anno_ngf$initL[match(is.extended,anno_ngf$txID)],
    Lnew=anno_ngf$maxL[match(is.extended,anno_ngf$txID)]
  )
  
  boxplot(myL,frame=F,outline=F,las=1,cex.axis=0.7)
  mtext(side=2,line=3,text="3'UTR length",cex=0.7)
  mtext(side=3,line=0,text=paste("# extended 3'UTR=",length(is.extended),sep=""),cex=0.6)
  
  par(mfrow=c(2,2))
  barplot(my.frac.sites.5[,13]*100,col="white",las=1,cex.axis=0.5)
  mtext(side=2,line=3,text="% of 3' UTR with PAS in ]-75:+50[ nt around 3' end",cex=0.5)
  barplot(my.frac.sites.6[,-13]*100,beside=T,col=c("white","black"),las=1,cex.axis=0.5)
  barplot(no.padb.4[-2]*100,beside=T,col=c("white",rep("black",6)),las=1,cex.axis=0.5)
  mtext(side=2,line=3,text="% of 3' UTR with PAdb in ]-75:+50[ nt around 3' end",cex=0.5)
  
  #Add 3'UTR length extended versus normal; extended only
  is.extended                   <- unique(anno_ngf$txID[(anno_ngf$newL-anno_ngf$initL)>10])
  myL                           <- list(
    Lconservative=anno_ngf$initL[match(is.extended,anno_ngf$txID)],
    Lnew=anno_ngf$maxL[match(is.extended,anno_ngf$txID)]
  )
  
  boxplot(myL,frame=F,outline=F,las=1,cex.axis=0.7)
  mtext(side=2,line=3,text="3'UTR length",cex=0.7)
  mtext(side=3,line=0,text=paste("# extended 3'UTR=",length(is.extended),sep=""),cex=0.6)
  
  par(mfrow=c(2,2))
  barplot(my.frac.sites.7[,13]*100,col="white",las=1,cex.axis=0.5)
  mtext(side=2,line=3,text="% of 3' UTR with PAS in ]-50:+50[ nt around 3' end",cex=0.5)
  barplot(my.frac.sites.8[,-13]*100,beside=T,col=c("white","black"),las=1,cex.axis=0.5)
  barplot(no.padb.6[-2]*100,beside=T,col=c("white",rep("black",6)),las=1,cex.axis=0.5)
  mtext(side=2,line=3,text="% of 3' UTR with PAdb in ]-50:+50[ nt around 3' end",cex=0.5)
  
  #Add 3'UTR length extended versus normal; extended only
  is.extended                   <- unique(anno_ngf$txID[(anno_ngf$newL-anno_ngf$initL)>10])
  myL                           <- list(
    Lconservative=anno_ngf$initL[match(is.extended,anno_ngf$txID)],
    Lnew=anno_ngf$maxL[match(is.extended,anno_ngf$txID)]
  )
  
  boxplot(myL,frame=F,outline=F,las=1,cex.axis=0.7)
  mtext(side=2,line=3,text="3'UTR length",cex=0.7)
  mtext(side=3,line=0,text=paste("# extended 3'UTR=",length(is.extended),sep=""),cex=0.6)
  
  dev.off()
  
  
  sel.RUD1 <- read.csv("/Users/luisier/Desktop/DataAnalsyisRiccio/Dec2016/axonal_remodelling/final_analysis_stringent/myselRUD_13062016.csv")[,-1]
  
  
  GetNoPAS <- function(target=myPA[[13]],utrGR=focusGR){
    # o) Remove redundant targets
    target                <- target[!duplicated(target)]
    # a) Find overlap
    myHits               <- subjectHits(findOverlaps(query=target,subject=utrGR,ignore.strand=F))
    # b) Find how many Hits in region of interest
    no.hits              <- as.data.frame(table(myHits))#get number of sites per isoform
    out                  <- rep(0,length(utrGR))
    names(out)           <- utrGR$ID
    out[as.numeric(as.character(no.hits$myHits))] <- as.numeric(as.character(no.hits$Freq))
    return(out)
  }
  
  
  
  
  
  
  
  GetFracWithSitesRemod <- function(my.nos=my.no.sites[,1],proxi=as.character(sel.RUD1[,3]),dista=as.character(sel.RUD1[,4])){
    
    
    proxi <- proxi[!is.na(proxi)]
    dista <- dista[!is.na(dista)]
    
    frac.prox  <- sum(my.nos[match(proxi,rownames(my.no.sites))]>0)/length(proxi)
    frac.dist  <- sum(my.nos[match(dista,rownames(my.no.sites))]>0)/length(dista)
    
    return(c(frac.prox,frac.dist))
  }
  
  
  
  GetAverageNumberWithLength <- function(colix=1,BINS=BINS1,sel=rep(TRUE,length(anno_ngf$newL))){
    my.nos <- my.no.sites[sel,colix]
    myidx  <- factor(as.character(BINS[sel]))
    return(tapply(my.nos,INDEX=myidx,FUN=mean))
  }
  
  GetSDNumberWithLength <- function(colix=1,BINS=BINS1,sel=rep(TRUE,length(anno_ngf$newL))){
    my.nos <- my.no.sites[sel,colix]
    myidx  <- factor(as.character(BINS[sel]))
    return(tapply(my.nos,INDEX=myidx,FUN=sd))
  }
  
  avg.per.length.cons <- do.call(what=rbind,args=lapply(c(1:ncol(my.no.sites)),function(Z)return(GetAverageNumberWithLength(colix=Z,BINS=BINS1,sel=anno_ngf$is.conservative))))
  avg.per.length.new <- do.call(what=rbind,args=lapply(c(1:ncol(my.no.sites)),function(Z)return(GetAverageNumberWithLength(colix=Z,BINS=BINS1,sel=!anno_ngf$is.conservative))))
  avg.per.length.tot <- do.call(what=rbind,args=lapply(c(1:ncol(my.no.sites)),function(Z)return(GetAverageNumberWithLength(colix=Z,BINS=BINS1,sel=rep(TRUE,length(my.no.sites[,1]))))))
  
  sd.per.length.cons <- do.call(what=rbind,args=lapply(c(1:ncol(my.no.sites)),function(Z)return(GetSDNumberWithLength(colix=Z,BINS=BINS1,sel=anno_ngf$is.conservative))))
  sd.per.length.new <- do.call(what=rbind,args=lapply(c(1:ncol(my.no.sites)),function(Z)return(GetSDNumberWithLength(colix=Z,BINS=BINS1,sel=!anno_ngf$is.conservative))))
  sd.per.length.tot <- do.call(what=rbind,args=lapply(c(1:ncol(my.no.sites)),function(Z)return(GetSDNumberWithLength(colix=Z,BINS=BINS1,sel=rep(TRUE,length(my.no.sites[,1]))))))
  
  
  
  
  
  
  avg.per.length.cons <- do.call(what=rbind,args=lapply(c(1:ncol(my.no.sites)),function(Z)return(GetAverageNumberWithLength(colix=Z,BINS=BINS2,sel=anno_ngf$is.conservative))))
  avg.per.length.new <- do.call(what=rbind,args=lapply(c(1:ncol(my.no.sites)),function(Z)return(GetAverageNumberWithLength(colix=Z,BINS=BINS2,sel=!anno_ngf$is.conservative))))
  avg.per.length.tot <- do.call(what=rbind,args=lapply(c(1:ncol(my.no.sites)),function(Z)return(GetAverageNumberWithLength(colix=Z,BINS=BINS2,sel=rep(TRUE,length(my.no.sites[,1]))))))
  
  sd.per.length.cons <- do.call(what=rbind,args=lapply(c(1:ncol(my.no.sites)),function(Z)return(GetSDNumberWithLength(colix=Z,BINS=BINS2,sel=anno_ngf$is.conservative))))
  sd.per.length.new <- do.call(what=rbind,args=lapply(c(1:ncol(my.no.sites)),function(Z)return(GetSDNumberWithLength(colix=Z,BINS=BINS2,sel=!anno_ngf$is.conservative))))
  sd.per.length.tot <- do.call(what=rbind,args=lapply(c(1:ncol(my.no.sites)),function(Z)return(GetSDNumberWithLength(colix=Z,BINS=BINS2,sel=rep(TRUE,length(my.no.sites[,1]))))))
  
  
  mp<-barplot(avg.per.length.tot[13,c(7,1:6)])
  error.bar(x=mp, y=avg.per.length.tot[13,c(7,1:6)], upper=sd.per.length.tot[13,c(7,1:6)], length=0.01,col="black")
  
  
  #Test significance of the enrichmetn using either Poisson or Negative Binomial
  dat         <-  data.frame(no.pas=my.no.sites[,13],length=BINS1,tL=log10(anno_ngf$newL))
  require(ggplot2)
  ggplot(dat, aes(no.pas, fill = length)) +
    geom_histogram(binwidth=0.5) +
    facet_grid(length ~ ., margins=TRUE, scales="free")
  
  with(dat, tapply(no.pas, length, function(x) {
    sprintf("M (SD) = %1.2f (%1.2f)", mean(x), sd(x))
  }))
  
  #Variance equals means
  require(MASS)
  summary(m1 <- glm.nb(no.pas ~length, data = dat))
  summary(m1 <- glm.nb(no.pas ~tL, data = dat))
  
  Coefficients:
    Estimate Std. Error z value Pr(>|z|)    
  (Intercept)       -0.38349    0.07732  -4.960 7.05e-07 ***
    length(1.5,2]      0.50391    0.08094   6.226 4.78e-10 ***
  length(2,2.25]     0.66723    0.07957   8.385  < 2e-16 ***
  length(2.25,2.75]  0.72018    0.07785   9.250  < 2e-16 ***
  length(2.75,3.25]  0.75126    0.07768   9.671  < 2e-16 ***
  length(3.25,3.75]  0.83407    0.07786  10.713  < 2e-16 ***
  length(3.75,4.25]  0.79916    0.07987  10.006  < 2e-16 ***
  
  
  #Check that the model is correct and not Poisson  
  #negative binomial models assume the conditional means are not equal to the conditional variances. This inequality is captured by estimating a dispersion parameter (not shown in the output) that is held constant in a Poisson model. Thus, the Poisson model is actually nested in the negative binomial model. We can then use a likelihood ratio test to compare these two and test this model assumption.
  m1 <- glm.nb(no.pas ~length, data = dat)
  m2 <- glm(no.pas ~length, family = "poisson", data = dat)
  m2 <- glm(no.pas ~tL, family = "poisson", data = dat)
  
  
  X2 <- 2 * (logLik(m1) - logLik(m2))
  X2
  'log Lik.' 907.6986 (df=8)
  pchisq(X2, df = 1, lower.tail=FALSE)
  'log Lik.' 2.16e-203 (df=5)
  #This very large chi-square strongly suggests the negative binomial model, which estimates the dispersion parameter, is more appropriate than the Poisson model.
  (est <- cbind(Estimate = coef(m1), confint(m1)))
  exp(est)
  
  newdata2 <- data.frame(tL = seq(from = min(dat$tL), to = max(dat$tL), length.out = 100))
  newdata2 <- cbind(newdata2, predict(m2, newdata2, type = "link", se.fit=TRUE))
  newdata2 <- within(newdata2, {
    no.pas <- exp(fit)
    LL <- exp(fit - 1.96 * se.fit)
    UL <- exp(fit + 1.96 * se.fit)
  })
  
  ggplot(newdata2,aes(tL, no.pas))+
    geom_ribbon(aes(ymin = LL, ymax = UL), alpha = .25) +
    geom_line( size = 2) +
    labs(x = "3' UTR length", y = "Predicted # PAS")
  
  StudyNumberPAS_v2 <- function(utrGR=focusGR,utrAnno=focusAnno,target=myPA[[1]],version="all"){
    
    if(version=="all"){
      print("no modif required")
    }
    if(version=="conservative"){
      utrGR   <- utrGR[utrAnno$is.conservative,]
      utrAnno <- utrAnno[utrAnno$is.conservative,]
    }
    if(version=="new"){
      utrGR   <- utrGR[!utrAnno$is.conservative,]
      utrAnno <- utrAnno[!utrAnno$is.conservative,]
    }
    
    BINS1                <- cut(log10(utrAnno$newL),breaks=c(1,1.5,2,2.25,2.75,3.25,3.75,4.25),include.lowest=T)
    print("debug1")
    Lev                  <- levels(BINS1)
    print("debug2")
    myListUTR            <- tapply(c(1:length(utrGR)),INDEX=BINS1,FUN=function(x)return(utrGR[x]))
    print("debug3")
    print(length(myListUTR))
    
    # Compute number of sites per category
    GetNumberSites <- function(subUTR=myListUTR[[1]]){
      if(length(subUTR)>10){
        myHits      <- subjectHits(findOverlaps(query=target,subject=subUTR,ignore.strand=F))
        print("debug4")
        no.missing  <- sum(!c(1:length(subUTR))%in%unique(myHits))
        print("debug5")
        if(length(myHits)==0){
          out         <- rep(0,length(subUTR))
          Out         <- rep(0,no.missing)
        }
        
        
        
        if(length(myHits>0)){
          out         <- as.data.frame(table(table(myHits)))#get number of sites per isoform
          print("debug6")
          temp       <- mean(c(as.data.frame(table(myHits))$Freq,rep(0,no.missing)))#Average number of PAS per isofor
          print("debug7")
          Out        <- c(rep(0,no.missing),do.call(lapply(c(1:nrow(out)),function(x)return(rep(out$Var1[x],out$Freq[x]))),what=c))#Total number of motifs
        }
      }
      else{
        temp=NA
        Out=NA
      }
      return(list(temp,Out))
    }
    
    no.sites                                <-  lapply(myListUTR,GetNumberSites)
    
    
    my.all.sites                            <- lapply(c(1:length(no.sites)),function(x)return(no.sites[[x]][[2]]))
    my.all.sites[[length(my.all.sites)+1]]  <- GetNumberSites(utrGR[utrAnno$is.conservative,])[[2]]
    my.all.sites                            <- my.all.sites[c((length(Lev)+1),c(1:length(Lev)))]
    names(my.all.sites)                     <- c("conservative",Lev)
    
    my.no.with.sites                       <- rbind(unlist(lapply(my.all.sites,function(x)return(sum(x!=0)))),
                                                    unlist(lapply(my.all.sites,function(x)return(length(x))))
    )
    
    rownames(my.no.with.sites)<- c("with sites", "tot")
    
    out          <- unlist(lapply(c(1:length(no.sites)),function(x)return(no.sites[[x]][[1]])))
    out          <- c(GetNumberSites(utrGR[utrAnno$is.conservative,])[[1]],out)#Average number of PAS per isoform
    names(out)   <-c("conservative",Lev)
    
    return(list(out,my.all.sites,my.no.with.sites))
  }
  
  
  
  
  freqwithsites      <- do.call(lapply(c(1:length(my.all.sites)),function(x)return(my.all.sites[[x]][[3]][1,]/my.all.sites[[x]][[3]][2,])),what=rbind)
  tot                 <- do.call(lapply(c(1:length(my.all.sites)),function(x)return(my.all.sites[[x]][[1]])),what=rbind)
  rownames(tot)<-rownames(freqwithsites)       <- names(myPA)
  
  mat                 <- t(scale(t(tot[-13,-1]),center = TRUE, scale = TRUE))
  mycols              <- colorRampPalette( c("cyan", "black", "magenta"), space="rgb")(31)
  b                   <- c(seq((-max(abs(mat))),0,length=16),seq(0,max(abs(mat)),length=16))
  
  pdf(paste(outdir,"no.sites.PAS.150.pdf",sep=""))
  heatmap.2(as.matrix(mat),keysize=1,col=mycols,breaks=b, scale="none",
            Rowv=TRUE,Colv=FALSE,dendr="row",key=TRUE,symkey=FALSE, density.info="none", trace="none",cexCol=1.0, cexRow=1.0, font.lab=1,las=2)
  mtext(side=3,line=0,text="[-100,+50]")
  barplot(tot[13,],col=c("white","magenta","blue","green","orange","grey","black"),ylab="average no.sites")
  dev.off()
  
  mat                 <- freqwithsites[-13,-1]
  mycols              <- colorRampPalette( c("cyan", "black", "magenta"), space="rgb")(31)
  b                   <- c(seq((-max(abs(mat))),0,length=16),seq(0,max(abs(mat)),length=16))
  
  pdf(paste(outdir,"freq.sites.PAS.150.pdf",sep=""))
  heatmap.2(as.matrix(mat),keysize=1,col=mycols,breaks=b, scale="none",
            Rowv=TRUE,Colv=FALSE,dendr="row",key=TRUE,symkey=FALSE, density.info="none", trace="none",cexCol=1.0, cexRow=1.0, font.lab=1,las=2)
  mtext(side=3,line=0,text="[-100,+50]")
  barplot(freqwithsites[13,],col=c("white","magenta","blue","green","orange","grey","black"),ylab="average freq sites")
  dev.off()
  
  
  # F2. Focus [-75,+50]
  require(gplots)
  focusGR            <- myUTR
  focusGR            <- resize(x=focusGR,fix="start",width=width(focusGR)+50,use.names=T)
  focusGR            <- resize(x=focusGR,fix="end",width=125,use.names=T)
  focusAnno          <- anno_ngf[match(names(focusGR),anno_ngf$uniqueID),]
  
  my.all.sites        <- lapply(myPA, function(x)return(StudyNumberPAS_v2(utrGR=focusGR,utrAnno=focusAnno,target=x,version="all")))#takes some time
  
  
  freqwithsites      <- do.call(lapply(c(1:length(my.all.sites)),function(x)return(my.all.sites[[x]][[3]][1,]/my.all.sites[[x]][[3]][2,])),what=rbind)
  tot                 <- do.call(lapply(c(1:length(my.all.sites)),function(x)return(my.all.sites[[x]][[1]])),what=rbind)
  rownames(tot)<-rownames(freqwithsites)       <- names(myPA)
  
  mat                 <- t(scale(t(tot[-13,-1]),center = TRUE, scale = TRUE))
  mycols              <- colorRampPalette( c("cyan", "black", "magenta"), space="rgb")(31)
  b                   <- c(seq((-max(abs(mat))),0,length=16),seq(0,max(abs(mat)),length=16))
  
  pdf(paste(outdir,"no.sites.PAS.125.pdf",sep=""))
  heatmap.2(as.matrix(mat),keysize=1,col=mycols,breaks=b, scale="none",
            Rowv=TRUE,Colv=FALSE,dendr="row",key=TRUE,symkey=FALSE, density.info="none", trace="none",cexCol=1.0, cexRow=1.0, font.lab=1,las=2)
  mtext(side=3,line=0,text="[-75,+50]")
  barplot(tot[13,],col=c("white","magenta","blue","green","orange","grey","black"),ylab="average no.sites")
  dev.off()
  
  mat                 <- freqwithsites[-13,-1]
  mycols              <- colorRampPalette( c("cyan", "black", "magenta"), space="rgb")(31)
  b                   <- c(seq((-max(abs(mat))),0,length=16),seq(0,max(abs(mat)),length=16))
  
  pdf(paste(outdir,"freq.sites.PAS.125.pdf",sep=""))
  heatmap.2(as.matrix(mat),keysize=1,col=mycols,breaks=b, scale="none",
            Rowv=TRUE,Colv=FALSE,dendr="row",key=TRUE,symkey=FALSE, density.info="none", trace="none",cexCol=1.0, cexRow=1.0, font.lab=1,las=2)
  mtext(side=3,line=0,text="[-75,+50]")
  barplot(freqwithsites[13,],col=c("white","magenta","blue","green","orange","grey","black"),ylab="average freq sites")
  dev.off()
  
  #Only Conservative
  my.all.sites        <- lapply(myPA, function(x)return(StudyNumberPAS_v2(utrGR=focusGR,utrAnno=focusAnno,target=x,version="conservative")))#takes some time
  freqwithsites      <- do.call(lapply(c(1:length(my.all.sites)),function(x)return(my.all.sites[[x]][[3]][1,]/my.all.sites[[x]][[3]][2,])),what=rbind)
  tot                 <- do.call(lapply(c(1:length(my.all.sites)),function(x)return(my.all.sites[[x]][[1]])),what=rbind)
  rownames(tot)<-rownames(freqwithsites)       <- names(myPA)
  mat                 <- t(scale(t(tot[-13,-1]),center = TRUE, scale = TRUE))
  mycols              <- colorRampPalette( c("cyan", "black", "magenta"), space="rgb")(31)
  b                   <- c(seq((-max(abs(mat))),0,length=16),seq(0,max(abs(mat)),length=16))
  pdf(paste(outdir,"no.sites.PAS.125.conservative.pdf",sep=""))
  heatmap.2(as.matrix(mat),keysize=1,col=mycols,breaks=b, scale="none",
            Rowv=TRUE,Colv=FALSE,dendr="row",key=TRUE,symkey=FALSE, density.info="none", trace="none",cexCol=1.0, cexRow=1.0, font.lab=1,las=2)
  mtext(side=3,line=0,text="[-75,+50]")
  barplot(tot[13,],col=c("white","magenta","blue","green","orange","grey","black"),ylab="average no.sites")
  dev.off()
  mat                 <- freqwithsites[-13,-1]
  mycols              <- colorRampPalette( c("cyan", "black", "magenta"), space="rgb")(31)
  b                   <- c(seq((-max(abs(mat))),0,length=16),seq(0,max(abs(mat)),length=16))
  pdf(paste(outdir,"freq.sites.PAS.125.conservative.pdf",sep=""))
  heatmap.2(as.matrix(mat),keysize=1,col=mycols,breaks=b, scale="none",
            Rowv=TRUE,Colv=FALSE,dendr="row",key=TRUE,symkey=FALSE, density.info="none", trace="none",cexCol=1.0, cexRow=1.0, font.lab=1,las=2)
  mtext(side=3,line=0,text="[-75,+50]")
  barplot(freqwithsites[13,],col=c("white","magenta","blue","green","orange","grey","black"),ylab="average freq sites")
  dev.off()
  
  
  #Only new
  my.all.sites        <- lapply(myPA, function(x)return(StudyNumberPAS_v2(utrGR=focusGR,utrAnno=focusAnno,target=x,version="new")))#takes some time
  freqwithsites      <- do.call(lapply(c(1:length(my.all.sites)),function(x)return(my.all.sites[[x]][[3]][1,]/my.all.sites[[x]][[3]][2,])),what=rbind)
  tot                 <- do.call(lapply(c(1:length(my.all.sites)),function(x)return(my.all.sites[[x]][[1]])),what=rbind)
  rownames(tot)<-rownames(freqwithsites)       <- names(myPA)
  mat                 <- t(scale(t(tot[-13,-1]),center = TRUE, scale = TRUE))
  mycols              <- colorRampPalette( c("cyan", "black", "magenta"), space="rgb")(31)
  b                   <- c(seq((-max(abs(mat))),0,length=16),seq(0,max(abs(mat)),length=16))
  pdf(paste(outdir,"no.sites.PAS.125.new.pdf",sep=""))
  heatmap.2(as.matrix(mat),keysize=1,col=mycols,breaks=b, scale="none",
            Rowv=TRUE,Colv=FALSE,dendr="row",key=TRUE,symkey=FALSE, density.info="none", trace="none",cexCol=1.0, cexRow=1.0, font.lab=1,las=2)
  mtext(side=3,line=0,text="[-75,+50]")
  barplot(tot[13,],col=c("white","magenta","blue","green","orange","grey","black"),ylab="average no.sites")
  dev.off()
  mat                 <- freqwithsites[-13,-1]
  mycols              <- colorRampPalette( c("cyan", "black", "magenta"), space="rgb")(31)
  b                   <- c(seq((-max(abs(mat))),0,length=16),seq(0,max(abs(mat)),length=16))
  pdf(paste(outdir,"freq.sites.PAS.125.new.pdf",sep=""))
  heatmap.2(as.matrix(mat),keysize=1,col=mycols,breaks=b, scale="none",
            Rowv=TRUE,Colv=FALSE,dendr="row",key=TRUE,symkey=FALSE, density.info="none", trace="none",cexCol=1.0, cexRow=1.0, font.lab=1,las=2)
  mtext(side=3,line=0,text="[-75,+50]")
  barplot(freqwithsites[13,],col=c("white","magenta","blue","green","orange","grey","black"),ylab="average freq sites")
  dev.off()
  
  
  #### From here seems better for the paper
  
  # C.1 Compute distance to next from closest to the end until the furthermore in conservative
  IX               <- match(anno_ngf$uniqueID,names(myRes[[1]]$no.motifs))[anno_ngf$is.conservative]
  my.no.motifs     <- myRes[[1]]$no.motifs[IX]
  my.dist.to.start <- myRes[[1]]$dist_from_end[IX]
  
  
  mysel            <- which(my.no.motifs>1)
  dist1            <- unlist(lapply(my.dist.to.start[mysel],function(x)return(diff(sort(unlist(x),decreasing=F))[1])))#It seems that I always select the first
  dist1            <- dist1[dist1>=25]
  mysel            <- which(my.no.motifs>2)
  dist2            <- unlist(lapply(my.dist.to.start[mysel],function(x)return(diff(sort(unlist(x),decreasing=F))[2])))
  dist2            <- dist2[dist2>=25]
  
  mysel            <- which(my.no.motifs>3)
  dist3            <- unlist(lapply(my.dist.to.start[mysel],function(x)return(diff(sort(unlist(x),decreasing=F))[3])))
  mysel            <- which(my.no.motifs>4)
  dist4            <- unlist(lapply(my.dist.to.start[mysel],function(x)return(diff(sort(unlist(x),decreasing=F))[4])))
  distall          <- data.frame(distTonext=c(dist1,dist2,dist3,dist4),category=c(rep(1,length(dist1)),rep(2,length(dist2)),rep(3,length(dist3)),rep(4,length(dist4))))
  
  
  # D. Analyse positional preference of the PAS from 3' end using the myPA found on GrossSegments (should result in identical results as previous section)
  subUTR  <- myUTR[anno_ngf$is.conservative,]
  subAnno <- anno_ngf[anno_ngf$is.conservative,]
  
  GetMyDist <- function(IX,subUTR=subUTR,subAnno=subAnno){
    gOver   <- findOverlaps(query=subUTR,subject=myPA[[IX]],ignore.strand=F)
    IX1     <- queryHits(gOver)
    IX2     <- subjectHits(gOver)
    POS     <- as.character(strand(subUTR[IX1,]))=="+"
    
    dist1   <- end(subUTR[IX1,])[POS] - start(myPA[[IX]][IX2,])[POS]
    dist2   <- -start(subUTR[IX1,])[!POS] + end(myPA[[IX]][IX2,])[!POS]
    
    mydist <- c(dist1,dist2)
    myL    <- c((subAnno$newL[IX1])[POS],(subAnno$newL[IX1])[!POS])
    myIX   <- c(IX1[POS],IX1[!POS])
    
    
    testF <- function(myRang){return(tapply(mydist,INDEX=factor(myIX),FUN=function(x){
      if(length(x>=myRang)){
        return(x[(sort(x,decreasing=FALSE,index.return=T)$ix)[myRang]])
        
      }
      else{return(NA)}
    }
    ))}
    myDist <- lapply(X=c(1,2,3,4,5),FUN=testF)
    return(myDist)
    
  }
  
  myDist_I0 <- lapply(c(1:length(myPA)),FUN=function(X)return(GetMyDist(IX=X,subUTR=myUTR[anno_ngf$is.conservative,],subAnno=anno_ngf[anno_ngf$is.conservative,])))
  myDist_IN <- lapply(c(1:length(myPA)),FUN=function(X)return(GetMyDist(IX=X,subUTR=myUTR[!anno_ngf$is.conservative,],subAnno=anno_ngf[!anno_ngf$is.conservative,])))
  
  
  
  pdf(paste(outdir,"distance_from_3end.pdf",sep=""))
  par(mar=c(1,2,2,1))
  layout(matrix(c(1:15),3,5,byrow = TRUE))
  for(I in c(1:length(myPA))){
    for(J in c(1:5)){
      plot(density(log10(myDist_IN[[I]][[J]][!is.na(myDist_IN[[I]][[J]])])),main="",frame=FALSE)
      lines(density(log10(myDist_I0[[I]][[J]][!is.na(myDist_I0[[I]][[J]])])),col="red")
      lines(density(log10(anno_ngf$newL)),col="grey",lty=2)
      mtext(side=3,line=0,text=names(myPA)[I])
    }
  }
  dev.off()
  
  
  GetMaxDensity <- function(myX=myDist_IN[[2]][[1]]){
    mydens <- density(log10(myX))
    return(10^min(mydens$x[which(mydens$y==max(mydens$y))]))
  }
  
  myPos <- do.call(args=lapply(X=myDist_I0,
                               FUN=function(nDist)return(unlist(
                                 lapply(X=nDist,FUN=function(Z)return(GetMaxDensity(Z[!is.na(Z)])))
                               )
                               )
  )
  ,what=rbind
  )
  
  
  localMaxima <- function(x) {
    # Use -Inf instead if x is numeric (non-integer)
    y <- diff(c(-.Machine$integer.max, x)) > 0L
    rle(y)$lengths
    y <- cumsum(rle(y)$lengths)
    y <- y[seq.int(1L, length(y), 2L)]
    if (x[[1]] == x[[2]]) {
      y <- y[-1]
    }
    y
  }
  
  
  peakfinder <- function(d){
    dh <- hist(d,plot=FALSE)
    ins <- dh[["density"]]
    nbins <- length(ins)
    ss <- which(rank(ins)%in%seq(from=nbins-2,to=nbins)) ## pick the top 3 intensities
    return(10^dh[["mids"]][ss])
  }
  
  myTest <- function(mytemp=myDist_I0[[2]][[1]]){
    mytemp <- log10(mytemp[!is.na(mytemp)])
    return(peakfinder(mytemp))
  }
  
  myTest(myDist_I0[[1]][[1]])
  myTest(myDist_I0[[2]][[1]])
  myTest(myDist_I0[[3]][[1]])
  
  #calculate turning points (extrema)
  require(pastecs)
  tp  <-turnpoints(ts_y)
  #plot
  pdf(paste(outdir,"test.pdf",sep=""))
  plot(mydens)
  points(mydens$x[tp$tppos],mydens$y[tp$tppos],col="red")
  dev.off()
  
  rownames(myPos)<- names(myPA)
  # B. Analyse distance between 2 PAS --> please note that in this analysis I simply do not consider introns which will bias the signal
  GetDistToNext <- function(IX,subUTR=subUTR,subAnno=subAnno){
    gOver   <- findOverlaps(query=subUTR,subject=myPA[[IX]],ignore.strand=F)
    IX1     <- queryHits(gOver)
    IX2     <- subjectHits(gOver)
    POS     <- as.character(strand(subUTR[IX1,]))=="+"
    
    dist1   <- end(subUTR[IX1,])[POS] - start(myPA[[IX]][IX2,])[POS]
    dist2   <- -start(subUTR[IX1,])[!POS] + end(myPA[[IX]][IX2,])[!POS]
    
    mydist <- c(dist1,dist2)
    myL    <- c((subAnno$newL[IX1])[POS],(subAnno$newL[IX1])[!POS])
    myIX   <- c(IX1[POS],IX1[!POS])
    
    
    testF <- function(myRang){return(tapply(mydist,INDEX=factor(myIX),FUN=function(x){
      if(length(x>=myRang)){
        X1 <- x[(sort(x,decreasing=FALSE,index.return=T)$ix)[myRang]]
        X2 <- x[(sort(x,decreasing=FALSE,index.return=T)$ix)[1]]
        return(X1-X2)
        
      }
      else{return(NA)}
    }
    ))}
    myDist <- lapply(X=c(2,3,4,5),FUN=testF)
    return(myDist)
    
  }
  
  
  myDistp_I0 <- lapply(c(1:length(myPA)),FUN=function(X)return(GetDistToNext(IX=X,subUTR=myUTR[anno_ngf$is.conservative,],subAnno=anno_ngf[anno_ngf$is.conservative,])))
  myDistp_IN <- lapply(c(1:length(myPA)),FUN=function(X)return(GetDistToNext(IX=X,subUTR=myUTR[!anno_ngf$is.conservative,],subAnno=anno_ngf[!anno_ngf$is.conservative,])))
  
  
  pdf(paste(outdir,"distance_to_next.pdf",sep=""))
  par(mar=c(1,2,2,1))
  layout(matrix(c(1:12),3,4,byrow = TRUE))
  for(I in c(1:length(myPA))){
    for(J in c(1:4)){
      plot(density(log10(myDistp_IN[[I]][[J]][!is.na(myDistp_IN[[I]][[J]])])),main="",frame=FALSE,xlim=c(0,5))
      lines(density(log10(myDistp_I0[[I]][[J]][!is.na(myDistp_I0[[I]][[J]])])),col="red")
      lines(density(log10(anno_ngf$newL)),col="grey",lty=2)
      lines(density(log10(c(myDistp_IN[[I]][[J]][!is.na(myDistp_IN[[I]][[J]])],myDistp_I0[[I]][[J]][!is.na(myDistp_I0[[I]][[J]])]))),col="green")
      mtext(side=3,line=0,text=names(myPA)[I])
    }
  }
  dev.off()
  
  
  
  --> I need to understand why this result is not consistent with previous finding, especially for the I0.
  
  require(gplots)
  focusGR            <- myUTR
  focusGR            <- resize(x=focusGR,fix="start",width=width(focusGR)+50,use.names=T)
  focusGR            <- resize(x=focusGR,fix="end",width=125,use.names=T)
  focusAnno          <- anno_ngf[match(names(focusGR),anno_ngf$uniqueID),]
  
  #Get no. sites per motif ID per isoforms
  my.no.sites <- do.call(what=cbind,args=lapply(myPA,function(x)return(GetNoPAS(target=x,utrGR=focusGR))))
  no.tot.sites<- apply(my.no.sites[,-13],1,sum)
  
  
  #
  # A.
  #
  BINS                <- cut(log10(focusAnno$newL),breaks=c(1,1.5,2,2.25,2.75,3.25,3.75,4.25),include.lowest=T)
  temp                <- tapply(no.tot.sites, INDEX=BINS,function(x)return(x))
  out                 <- matrix(nrow=length(temp),ncol=length(temp))
  for(i in c(1:length(temp))){
    for(j in c(1:length(temp))){
      out[i,j]<- t.test(temp[[i]],temp[[j]])$p.value
    }
  }
  colnames(out)<-rownames(out)<-names(temp)
  
  
  
  #Focus on those which have at least 2 isoforms (to have a good idea on the synchrony)
  
  mysel          <- rep(TRUE,length(focusAnno$no.iso))
  outdir <- "/home/rluisier/data/Riccio/Exp_1/Dec2016/utrid/pas_stringent/all/"
  mysel          <- focusAnno$no.iso>1#24'275; new 26'368
  outdir <- "/home/rluisier/data/Riccio/Exp_1/Dec2016/utrid/pas_stringent/with_atleast_2_iso/"
  mysel          <- focusAnno$no.iso>2#new =  19578
  outdir <- "/home/rluisier/data/Riccio/Exp_1/Dec2016/utrid/pas_stringent/with_atleast_3_iso/"
  mysel          <- focusAnno$no.iso>3#13'501; new =12'645
  outdir <- "/home/rluisier/data/Riccio/Exp_1/Dec2016/utrid/pas_stringent/with_atleast_4_iso/"
  
  mysel          <- focusAnno$no.iso>4#7'185
  outdir <- "/home/rluisier/data/Riccio/Exp_1/Dec2016/utrid/pas_stringent/with_atleast_5_iso/"
  
  
  #Compute average no and average frewuency of sites in FG
  avg.no.sites   <- do.call(what=rbind,lapply(c(1:ncol(my.no.sites)), FUN=function(Z)return(tapply(my.no.sites[mysel,Z],INDEX=BINS[mysel],FUN=mean))))
  avg.frac.with.sites<- do.call(what=rbind,lapply(c(1:ncol(my.no.sites)), FUN=function(Z)return(tapply(my.no.sites[mysel,Z],INDEX=BINS[mysel],FUN=function(Y)return(sum(Y!=0)/length(Y))))))
  rownames(avg.no.sites)<-rownames(avg.frac.with.sites)<- colnames(my.no.sites)
  
  
  #Compute average no and average frewuency of sites in BG == shuffled background
  Nboot=1000
  myIdx <- matrix(NA,nrow=length(BINS[mysel]),ncol=Nboot)
  for(k in c(1:Nboot)){
    myIdx[,k]<- sample(as.character(BINS[mysel]))
  }
  
  avg.no.sites.bg         <- lapply(c(1:ncol(my.no.sites)),FUN=function(Z)return(do.call(what=cbind,lapply(c(1:Nboot), FUN=function(W)return(tapply(my.no.sites[mysel,Z],INDEX=factor(myIdx[,W]),FUN=mean))))))
  avg.frac.with.sites.bg <- lapply(c(1:ncol(my.no.sites)),FUN=function(Z)return(do.call(what=cbind,lapply(c(1:Nboot), FUN=function(W)return(tapply(my.no.sites[mysel,Z],INDEX=factor(myIdx[,W]),FUN=function(Q)return(sum(Q!=0)/length(Q))))))))
  
  
  
  
  # Z-value
  mymean.no    <- do.call(lapply(c(1:length(avg.no.sites.bg)), FUN=function(X)return(apply(avg.no.sites.bg[[X]],1,mean))),what=rbind)
  mysd.no      <- do.call(lapply(c(1:length(avg.no.sites.bg)), FUN=function(X)return(apply(avg.no.sites.bg[[X]],1,sd))),what=rbind)
  mymean.frac  <- do.call(lapply(c(1:length(avg.frac.with.sites.bg)), FUN=function(X)return(apply(avg.frac.with.sites.bg[[X]],1,mean))),what=rbind)
  mysd.frac    <- do.call(lapply(c(1:length(avg.frac.with.sites.bg)), FUN=function(X)return(apply(avg.frac.with.sites.bg[[X]],1,sd))),what=rbind)
  myzval.frac  <- (avg.frac.with.sites-mymean.frac)/mysd.frac
  myzval.no    <- (avg.no.sites-mymean.no)/mysd.no
  
  
  #Fisher exact test on each category and each isoforms
  mytest.greater <-   do.call(args=lapply(c(1:ncol(my.no.sites)),
                                          FUN=function(Z)return(do.call(args=lapply(levels(BINS[mysel]),function(X)return(fisher.test(BINS[mysel]==X,my.no.sites[mysel,Z]>0,alternative="greater")$p.value)),what=c))),what=rbind)
  mytest.less <-   do.call(args=lapply(c(1:ncol(my.no.sites)),
                                       FUN=function(Z)return(do.call(args=lapply(levels(BINS[mysel]),function(X)return(fisher.test(BINS[mysel]==X,my.no.sites[mysel,Z]>0,alternative="less")$p.value)),what=c))),what=rbind)
  mytest.greater <- -log10(mytest.greater)
  mytest.less    <- -log10(mytest.less)
  myfisher.test  <- mytest.greater
  myfisher.test[mytest.greater<mytest.less]<- -(mytest.less[mytest.greater<mytest.less])
  colnames(mytest.greater)<-colnames(mytest.less)<-colnames(myfisher.test)<- levels(BINS)
  rownames(mytest.greater)<-rownames(mytest.less)<-rownames(myfisher.test)<- colnames(my.no.sites)
  
  #>cor(as.vector(myzval.no ),as.vector(myfisher.test))
  #[1] 0.8874058
  #> cor(as.vector(myzval.frac),as.vector(myfisher.test))
  #[1] 0.9709798
  #> cor(as.vector(myzval.frac),as.vector(myzval.no))
  #[1] 0.9392339
  
  #Fiserh exact test between conservative and new
  mytest           <-cbind(
    unlist(lapply(c(1:ncol(my.no.sites[mysel,])),function(X)return(fisher.test(focusAnno$is.conservative[mysel],my.no.sites[mysel,X]>0,alternative="greater")$p.value))),
    unlist(lapply(c(1:ncol(my.no.sites[mysel,])),function(X)return(fisher.test(focusAnno$is.conservative[mysel],my.no.sites[mysel,X]>0,alternative="less")$p.value)))
  )
  rownames(mytest)<- colnames(my.no.sites)
  colnames(mytest) <- c("more in rn5","less in rn5")
  
  
  rel.frac.rn5         <- apply(my.no.sites[anno_ngf$is.conservative,],2,function(x)return(sum(x>0)/length(x)))[-13]
  rel.frac.novel       <- apply(my.no.sites[(!anno_ngf$is.conservative),],2,function(x)return(sum(x>0)/length(x)))[-13]
  names(rel.frac.rn5)<-names(rel.frac.novel)<- colnames(my.no.sites)[-13]
  
  
  pdf(paste(outdir,"PAS_zval_no_iso_all.pdf",sep=""))
  par(mfrow=c(2,2))
  pie(rel.frac.rn5)
  mtext(side=3,text=paste("Rn5 (n=",sum(anno_ngf$is.conservative&mysel),")",sep=""),line=0)
  pie(rel.frac.novel)
  mtext(side=3,text=paste("Novel (n=",sum(!anno_ngf$is.conservative),")",sep=""),line=0)
  par(mfrow=c(1,1))
  
  
  tot                 <- mytest.greater[-13,]
  mycols              <- colorRampPalette( c("white",rgb(221/255,119/255,136/255),rgb(170/255,68/255,85/255),rgb(119/255,17/255,34/255)), space="rgb")(201)
  b                   <- seq(0,max(tot),length.out=202)
  heatmap.2(tot,keysize=1,col=mycols,breaks=b, scale="none",Rowv=TRUE,Colv=FALSE,dendr="row",key=TRUE,symkey=FALSE,distfun=function(x)return(as.dist(1-cor(t(x)))),hclustfun=function(x)return(hclust(x,method="average")), density.info="none", trace="none",cexCol=1.0, cexRow=1.0, font.lab=1,las=2)
  mtext(side=3,line=0,text="fisher exact test P-val average fraction with sites")
  
  temp<-rbind(
    apply(my.no.sites[anno_ngf$is.conservative,],2,function(x)return(sum(x>0)/sum(anno_ngf$is.conservative))),
    apply(my.no.sites[!anno_ngf$is.conservative,],2,function(x)return(sum(x>0)/sum(!anno_ngf$is.conservative)))
  )
  
  barplot(temp,col=c(rgb(170/255,68/255,153/255),rgb(51/255,34/255,136/255)),beside=T,las=1,cex.names=0.3)
  
  tot                 <- myzval.frac[-13,]
  mycols              <- colorRampPalette( c(rgb(51/255,34/255,136/255),rgb(136/255,204/255,238/255),"white",rgb(204/255,102/255,119/255),rgb(136/255,34/255,85/255)), space="rgb")(101)
  b                   <- c(seq((-max(abs(tot))),0,length=51),seq(0,max(abs(tot)),length=51))
  heatmap.2(tot,keysize=1,col=mycols,breaks=b, scale="none",Rowv=TRUE,Colv=FALSE,dendr="row",key=TRUE,symkey=FALSE, density.info="none", trace="none",cexCol=1.0, cexRow=1.0, font.lab=1,las=2)
  mtext(side=3,line=0,text="Z-valwith average fraction with sites")
  
  tot                 <- avg.frac.with.sites[-13,]
  tot                 <- t(scale(t(tot[-13,]),center = TRUE, scale = TRUE))
  b                   <- c(seq((-max(abs(tot))),0,length=101),seq(0,max(abs(tot)),length=101))
  mycols              <- colorRampPalette( c(rgb(51/255,34/255,136/255),rgb(136/255,204/255,238/255),"white",rgb(204/255,102/255,119/255),rgb(136/255,34/255,85/255)), space="rgb")(201)
  heatmap.2(tot,keysize=1,col=mycols,breaks=b, scale="none",Rowv=TRUE,Colv=FALSE,dendr="row",key=TRUE,symkey=FALSE, density.info="none", trace="none",cexCol=1.0, cexRow=1.0, font.lab=1,las=2)
  mtext(side=3,line=0,text="fraction with sites")
  
  tot                 <- myzval.no [-13,]
  b                   <- c(seq((-max(abs(tot))),0,length=101),seq(0,max(abs(tot)),length=101))
  heatmap.2(tot,keysize=1,col=mycols,breaks=b, scale="none",Rowv=TRUE,Colv=FALSE,dendr="row",key=TRUE,symkey=FALSE, density.info="none", trace="none",cexCol=1.0, cexRow=1.0, font.lab=1,las=2)
  mtext(side=3,line=0,text="Z-valwith no. sites")
  
  
  tot                 <- avg.no.sites[-13,]
  tot                 <- t(scale(t(tot[-13,]),center = TRUE, scale = TRUE))
  b                   <- c(seq((-max(abs(tot))),0,length=101),seq(0,max(abs(tot)),length=101))
  mycols              <- colorRampPalette( c(rgb(51/255,34/255,136/255),rgb(136/255,204/255,238/255),"white",rgb(204/255,102/255,119/255),rgb(136/255,34/255,85/255)), space="rgb")(201)
  heatmap.2(tot,keysize=1,col=mycols,breaks=b, scale="none",Rowv=TRUE,Colv=FALSE,dendr="row",key=TRUE,symkey=FALSE, density.info="none", trace="none",cexCol=1.0, cexRow=1.0, font.lab=1,las=2)
  mtext(side=3,line=0,text="no. sites")
  
  #avg.frac.with.sites  <- cbind(avg.frac.with.sites,do.call(what=c,lapply(c(1:ncol(my.no.sites)), FUN=function(Z)return(sum(my.no.sites[anno_ngf$is.conservative,Z]>0)/sum(anno_ngf$is.conservative))))
  
  #avg.no.sites<- cbind(avg.no.sites,do.call(what=c,lapply(c(1:ncol(my.no.sites)), FUN=function(Z)return(mean(my.no.sites[anno_ngf$is.conservative,Z]))))
  
  barplot(avg.no.sites[13,],col=c("white","magenta","blue","green","orange","grey","black"),ylab="average no.sites")
  
  
  dev.off()
  
  #Test of significance between fraction of sites with motifs between the different categories using GLM
  data        <-  data.frame(count=my.no.sites[,13],length=BINS)
  data        <-  data.frame(count=my.no.sites[,13],length.1=BINS==levels(BINS)[1],
                             length.2=BINS==levels(BINS)[2],
                             length.3=BINS==levels(BINS)[3],
                             length.4=BINS==levels(BINS)[4],
                             length.5=BINS==levels(BINS)[5],
                             length.6=BINS==levels(BINS)[6],
                             length.7=BINS==levels(BINS)[7],
                             length=BINS)
  
  test0        <-  glm(count~1,family=poisson,data=data)
  mytest <- list(
    test1        <-  glm(count~1+length.1,family=poisson,data=data),
    test2        <-  glm(count~1+length.2,family=poisson,data=data),
    test3        <-  glm(count~1+length.3,family=poisson,data=data),
    test4        <-  glm(count~1+length.4,family=poisson,data=data),
    test5        <-  glm(count~1+length.5,family=poisson,data=data),
    test6        <-  glm(count~1+length.6,family=poisson,data=data),
    test7        <-  glm(count~1+length.7,family=poisson,data=data)
  )
  test0        <-  glm(count~1+length,family=poisson,data=data)
  summary(test0)
  zvals        <- unlist(lapply(mytest,function(x)return(summary(x)$coefficients[2,3])))
  
  
  
  testall        <-  glm(count~1+length,family=poisson,data=data)
  
  
  test.fitted <-  predict(test0, type = "response")
  data        <- data[match(names(test.fitted),rownames(data)),]
  confint(test0)
  test0$fitted.values
  
  PoisonFun <- function(lambda,x){return(exp(-lambda)*lambda^x/factorial(x))}
  mylambda <- test0$coef
  myX      <- seq(from=0,b=0.1,to=10)
  
  pdf(paste(outdir,"test_poisso.pdf",sep=""))
  par(mfrow=c(2,4))
  hist(my.no.sites[,3],freq=FALSE)
  for(i in c(1:length(levels(BINS3)))){
    hist(mylist[[i]],freq=F,ylim=c(0,1))
    lines(myX,PoisonFun(mylambda[i],myX),col=mycols[i],lwd=2)
  }
  par(mfrow=c(1,1))
  i=1
  plot(myX,PoisonFun(mylambda[i],myX),col=mycols[i],lwd=2,type="l",frame=F,las=1,ylim=c(0,0.6))
  for(i in c(2:length(levels(BINS3)))){
    lines(myX,PoisonFun(mylambda[i],myX),col=mycols[i],lwd=2)
  }
  dev.off()
  
  
  
  
  
  ##### Study preferential positioning of the motif along 3' UTR: both distance from 3' end and distance to next motif
  
  GetNoPAS <- function(target=myPA[[13]]){
    # o) Remove redundant targets
    target <- target[!duplicated(target)]
    extUTR <- myUTR
    extUTR <- resize(x=extUTR,fix="start",width=width(extUTR)+100,use.names=T)
    
    # a) Find overlap
    gOver              <- findOverlaps(query=target,subject=extUTR,ignore.strand=F)
    quj                <- queryHits(gOver)
    subj               <- subjectHits(gOver)
    
    
    # b) Find all distance from 3' end
    g1                 <- as.numeric(start(target))[quj]
    g2                 <- as.numeric(end(extUTR))[subj]
    distoend           <- g1-g2+100
    POS                <- which(as.character(strand(extUTR))=="+")
    ix                 <- subj%in%POS
    g1[!ix]            <- as.numeric(end(target))[quj[!ix]]
    g2[!ix]            <- as.numeric(start(extUTR))[subj[!ix]]
    distoend[!ix]      <- (-g1[!ix]+g2[!ix]+100)
    
    # c) Extract closest by uniqueID
    myClosest<- tapply(distoend,INDEX=factor(anno_ngf$uniqueID[subj]),FUN=function(X)return(X[which(abs(X)==min(abs(X)))[1]]))
    
    # d) Extract all distance by uniqueID
    myDist  <- tapply(distoend,INDEX=factor(anno_ngf$uniqueID[subj]),FUN=function(X)return(X))
    
    # e) Create vector of position and frequency
    test1 <- as.data.frame(table(myClosest))
    test2 <- as.data.frame(table(unlist(myDist)))
    
    mypos2                                                    <- c(min(as.numeric(as.character(test2[,1]))):max(as.numeric(as.character(test2[,1]))))
    ntest2                                                    <- rep(0,length(mypos2))
    ntest2[match(as.numeric(as.character(test2[,1])),mypos2)] <- as.numeric(as.character(test2[,2]))
    
    ntest1                                                    <- rep(0,length(mypos2))
    ntest1[match(as.numeric(as.character(test1[,1])),mypos2)] <- as.numeric(as.character(test1[,2]))
    
    myFreq.all                                                <- cbind(mypos2,ntest2)
    myFreq.first                                              <- cbind(mypos2,ntest1)
    colnames(myFreq.all)<-colnames(myFreq.first)<- c("distance_from_3","all")
    
    
    for(i in c("-1","0","1","2","3")){
      sel1   <- names(myClosest)%in%as.character(anno_ngf$uniqueID)[anno_ngf$iso_cons_merged==i]
      test1  <- as.data.frame(table(myClosest[sel1]))
      sel2   <- names(myDist)%in%as.character(anno_ngf$uniqueID)[anno_ngf$iso_cons_merged==i]
      test2  <- as.data.frame(table(unlist(myDist[sel2])))
      
      ntest2                                                    <- rep(0,length(mypos2))
      ntest2[match(as.numeric(as.character(test2[,1])),mypos2)] <- as.numeric(as.character(test2[,2]))
      
      ntest1                                                    <- rep(0,length(mypos2))
      ntest1[match(as.numeric(as.character(test1[,1])),mypos2)] <- as.numeric(as.character(test1[,2]))
      
      myFreq.all                                                <- cbind(myFreq.all,ntest2)
      myFreq.first                                              <- cbind(myFreq.first,ntest1)
      
    }
    colnames( myFreq.all)[c(3:7)]<- colnames(myFreq.first)[c(3:7)]<-c("-1","0","1","2","3")
    out <- list(myClosest,myDist,myFreq.all, myFreq.first)
    names(out)<- c("myClosest","myAll","myFreq.all","myFreq.first")
    return(out)
  }
  
  GetNoPAS_v2 <- function(target=myPA[[13]]){
    # o) Remove redundant targets
    target <- target[!duplicated(target)]
    extUTR <- myUTR
    extUTR <- resize(x=extUTR,fix="start",width=width(extUTR)+100,use.names=T)#add 100 nt the end of each 3' UTR
    
    # a) Find overlap
    gOver              <- findOverlaps(query=target,subject=extUTR,ignore.strand=F)
    quj                <- queryHits(gOver)
    subj               <- subjectHits(gOver)
    
    
    # b) Find all distance from 3' end
    g1                 <- as.numeric(start(target))[quj]
    g2                 <- as.numeric(end(extUTR))[subj]
    distoend           <- g1-g2+100
    POS                <- which(as.character(strand(extUTR))=="+")
    ix                 <- subj%in%POS
    g1[!ix]            <- as.numeric(end(target))[quj[!ix]]
    g2[!ix]            <- as.numeric(start(extUTR))[subj[!ix]]
    distoend[!ix]      <- (-g1[!ix]+g2[!ix]+100)
    
    # c) Extract closest by uniqueID
    myClosest<- tapply(distoend,INDEX=factor(anno_ngf$uniqueID[subj]),FUN=function(X)return(X[which(abs(X)==min(abs(X)))[1]]))
    
    # d) Extract all distance by uniqueID
    myDist  <- tapply(distoend,INDEX=factor(anno_ngf$uniqueID[subj]),FUN=function(X)return(X))
    
    # e) Create vector of position and frequency
    test1 <- as.data.frame(table(myClosest))
    test2 <- as.data.frame(table(unlist(myDist)))
    
    mypos2                                                    <- c(min(as.numeric(as.character(test2[,1]))):max(as.numeric(as.character(test2[,1]))))
    ntest2                                                    <- rep(0,length(mypos2))
    ntest2[match(as.numeric(as.character(test2[,1])),mypos2)] <- as.numeric(as.character(test2[,2]))
    
    ntest1                                                    <- rep(0,length(mypos2))
    ntest1[match(as.numeric(as.character(test1[,1])),mypos2)] <- as.numeric(as.character(test1[,2]))
    
    myFreq.all                                                <- cbind(mypos2,ntest2)
    myFreq.first                                              <- cbind(mypos2,ntest1)
    colnames(myFreq.all)<-colnames(myFreq.first)<- c("distance_from_3","all")
    
    
    for(i in c("-1","0","1","2","3")){
      sel1   <- names(myClosest)%in%as.character(anno_ngf$uniqueID)[anno_ngf$iso_cons_merged==i]
      test1  <- as.data.frame(table(myClosest[sel1]))
      sel2   <- names(myDist)%in%as.character(anno_ngf$uniqueID)[anno_ngf$iso_cons_merged==i]
      test2  <- as.data.frame(table(unlist(myDist[sel2])))
      
      ntest2                                                    <- rep(0,length(mypos2))
      ntest2[match(as.numeric(as.character(test2[,1])),mypos2)] <- as.numeric(as.character(test2[,2]))
      
      ntest1                                                    <- rep(0,length(mypos2))
      ntest1[match(as.numeric(as.character(test1[,1])),mypos2)] <- as.numeric(as.character(test1[,2]))
      
      myFreq.all                                                <- cbind(myFreq.all,ntest2)
      myFreq.first                                              <- cbind(myFreq.first,ntest1)
      
    }
    colnames( myFreq.all)[c(3:7)]<- colnames(myFreq.first)[c(3:7)]<-c("-1","0","1","2","3")
    out <- list(myClosest,myDist,myFreq.all, myFreq.first)
    names(out)<- c("myClosest","myAll","myFreq.all","myFreq.first")
    return(out)
  }
  
  
  
  pos.PAS<- lapply(myPA,FUN=function(X)return(GetNoPAS(target=X)))
  
  outdir <- "/home/rluisier/data/Riccio/Exp_1/Dec2016/utrid/pas_stringent/preferential_positioning/"
  
  mycols <- c("magenta","blue","orange","green","grey")
  pdf(paste(outdir,"position_first_peaks_v2.pdf",sep=""))
  par(mfrow=c(5,1),mar=c(3,3,3,3))
  for(i in c(1:12)){
    for(k in c(1:5)){
      myPos<- pos.PAS[[i]][[4]][,1]
      myNum<- pos.PAS[[i]][[4]][,k+2]
      myNum<- myNum/max(myNum)
      plot(x=rev(rev(myPos)[c(50:300)]),y=rev(rev(myNum)[c(50:300)]),type="l",col=mycols[k],las=1,frame=FALSE,lwd=2)
      if(k==1){mtext(side=3,line=0,text=names(pos.PAS)[i])}
    }
  }
  dev.off()
  
  
  pdf(paste(outdir,"position_all_peaks.pdf",sep=""))
  par(mfrow=c(5,1),mar=c(3,3,3,3))
  for(i in c(1:12)){
    for(k in c(1:5)){
      myPos<- pos.PAS[[i]][[3]][,1]
      myNum<- pos.PAS[[i]][[3]][,k+2]
      myNum<- myNum/max(myNum)
      plot(x=rev(rev(myPos)[c(50:300)]),y=rev(rev(myNum)[c(50:300)]),type="l",col=mycols[k],las=1,frame=FALSE,lwd=2)
      if(k==1){mtext(side=3,line=0,text=names(pos.PAS)[i])}
    }
  }
  dev.off()
  
  
  
  #Study distance to next
  
  temp1 <- lapply(c(1:length(pos.PAS)),function(x)return(pos.PAS[[x]][[2]]))
  #Focus on Conservative initially
  temp2 <- lapply(c(1:length(temp1)),function(x)return(temp1[[x]][grep(names(temp1[[x]]),pattern="\\.0")]))
  temp2 <- temp1
  #Select only those which have more than one site
  temp3 <- lapply(c(1:length(temp2)),function(x)return(temp2[[x]][unlist(lapply(temp2[[x]],length))>1]))
  #Get distance to next
  GetDistToNext  <- function(x)return(diff(sort(x,decreasing=F)))
  GetDistToFirst <- function(x){
    x <- x[x<=0]
    return(diff(sort(x,decreasing=F))[1])
  }
  distonext  <- lapply(c(1:length(temp3)),function(x)return(log10(unlist(lapply(temp3[[x]],GetDistToNext)))))
  distofirst <- lapply(c(1:length(temp3)),function(x)return(log10(unlist(lapply(temp3[[x]],GetDistToFirst)))))
  motifs     <- names(myPA)
  pdf(paste(outdir,"distonext.pdf",sep=""))
  par(mfrow=c(5,2))
  lapply(c(1:length(distonext)),function(x)return(plot(density(distonext[[x]]),main=motifs[x])))
  par(mfrow=c(5,2))
  lapply(c(1:length(distofirst)),function(x)return(plot(density(distofirst[[x]],na.rm=T),main=motifs[x])))
  dev.off()
  #                                          #
  #          ----- To be continued ------    #
  #                                          #
  
  
  dat <- my.all.sites[[1]][[3]][,-1]
  myPval <-myZval<- matrix(0,nrow=6,ncol=6)
  rownames(myPval)<-colnames(myPval)<-rownames(myZval)<-colnames(myZval)<-colnames(dat)
  for(i in c(1:ncol(dat))){
    for(j in c(1:ncol(dat))){
      myZval[i,j] <- MyTestFrac(n1=dat[2,i],n2=dat[2,j],p1=dat[1,i]/dat[2,i],p2=dat[1,j]/dat[2,j])[[1]]
      myPval[i,j] <- -log10(MyTestFrac(n1=dat[2,i],n2=dat[2,j],p1=dat[1,i]/dat[2,i],p2=dat[1,j]/dat[2,j])[[2]])
    }
  }
  
  par(mfrow=c(3,2))
  for(i in c(1:12)){
    dat <- my.all.sites[[i]][[3]][,-1]
    barplot(dat[1,]/dat[,2])
  }
  
  
  
  
  # Creat count matrix of how many
  motifs              <- unique(as.character(myPAS$rep.motifID..length.end..))
  myPASp              <- lapply(motifs,function(x)return(myPAS[which(as.character(myPAS$rep.motifID..length.end..)==x)]))
  myCounts            <- matrix(0,nrow=length(subGR),ncol=length(myPASp))
  colnames(myCounts)  <- motifs
  rownames(myCounts)  <- as.character(subGR$ID)
  
  for(i in c(1:ncol(myCounts))){
    mytar                      <- queryHits(findOverlaps(query=subGR,subject=myPASp[[i]],ignore.strand=FALSE))
    no.counts                  <- as.data.frame(table(mytar))
    myCounts[as.numeric(as.character(no.counts$mytar)),i]<-no.counts$Freq
    print(i)
  }
  
  # C. Statts on count
  count.tot     <- apply(myCounts,1,sum)
  
  sum(count.tot==0&as.character(subGR$isoform)=="0")
  sum(count.tot>0&as.character(subGR$isoform)=="0")
  sum(count.tot==0&as.character(subGR$isoform)=="0"&width(myGR)>10)/sum(as.character(subGR$isoform)=="0"&width(myGR)>10)
  %(30% w/0 site)

#Frequency of absence of a motif
barplot(c(sum(count.tot==0&as.character(subGR$isoform)=="0")/sum(as.character(subGR$isoform)=="0"),sum(count.tot==0&as.character(subGR$isoform)!="0")/sum(as.character(subGR$isoform)!="0")))
barplot(c(sum(count.tot==0&as.character(subGR$isoform)=="0"&width(myGR)>10)/sum(as.character(subGR$isoform)=="0"&width(myGR)>10),sum(count.tot==0&as.character(subGR$isoform)!="0")/sum(as.character(subGR$isoform)!="0")))



#For those which do not intersect with any motif, look for distnotnex
sel1        <- count.tot>0|as.character(myGR$isoform)=="0"#46'048
subGR       <- resize(myGR,width=6,fix="end")
no.pas      <- subGR[!sel1,]
mydist.pas  <-distanceToNearest(x=no.pas,myPAS,select="arbitrary",ignore.strand=FALSE)
tomodif

30948       29873        77

mydist[[3]]


L2  <- myGR[sel,]
export(L2,con="L2.gtf",format="gtf")




# 5. Characterisation of pA elements along 3' UTR
# --for mi1,mi2,mi3,mi4 only when no.mi>=4
# --

IsWithin25  <- do.call(IsWithinDist(mydist=25),what=cbind)
IsWithin50  <- do.call(IsWithinDist(mydist=50),what=cbind)
IsWithin100 <- do.call(IsWithinDist(mydist=100),what=cbind)
IsWithin150 <- do.call(IsWithinDist(mydist=150),what=cbind)
IsWithin200 <- do.call(IsWithinDist(mydist=200),what=cbind)

GetFreqAllMotif <- function(allMotifs=apply(IsWithin,1,function(x)return(sum(x)>0))){
  
  no.conservative            <- sum(allMotifs[anno_ngf$conservative&anno_ngf$newL>10])
  freq.conservative          <- no.conservative/sum(anno_ngf$conservative&anno_ngf$newL>10)
  no.new                     <- sum(allMotifs[!anno_ngf$conservative&anno_ngf$newL>10])
  freq.new                   <- no.new/sum(!anno_ngf$conservative&anno_ngf$newL>10)
  no.new.cat1                <- sum(allMotifs[!anno_ngf$conservative&abs(anno_ngf$dL)<150&abs(anno_ngf$dL)>50])
  freq.new.cat1              <- no.new.cat1/sum(!anno_ngf$conservative&abs(anno_ngf$dL)<150&abs(anno_ngf$dL)>50)
  no.new.cat2                <- sum(allMotifs[!anno_ngf$conservative&abs(anno_ngf$dL)>150])
  freq.new.cat2              <- no.new.cat2/sum(!anno_ngf$conservative&abs(anno_ngf$dL)>150)
  
  out <- cbind(c(no.conservative,no.new,no.new.cat1,no.new.cat2),c(freq.conservative,freq.new,freq.new.cat1,freq.new.cat2))
  colnames(out)<-c("no.","fraction with at least one motif")
  rownames(out)<-c("conservative","new","new.less.150","new.more.150")
  return(out)
}

freq1 <- GetFreqAllMotif(allMotifs=apply(IsWithin25,1,function(x)return(sum(x)>0)))
freq2 <- GetFreqAllMotif(allMotifs=apply(IsWithin50,1,function(x)return(sum(x)>0)))
freq3 <- GetFreqAllMotif(allMotifs=apply(IsWithin100,1,function(x)return(sum(x)>0)))
freq4 <- GetFreqAllMotif(allMotifs=apply(IsWithin150,1,function(x)return(sum(x)>0)))

freqall<- do.call(lapply(c(1:length(motif)),function(x)return(GetFreqAllMotif(allMotifs=IsWithin50[,x]))),what=cbind)
colnames(freqall)[grep(pattern="fraction",colnames(freqall))]<-paste("freq. ",unlist(lapply(motif,toString)),sep="")
colnames(freqall)[grep(pattern="no",colnames(freqall))]<-paste("no.UTR with ",unlist(lapply(motif,toString)),sep="")

pdf(paste(outdir,"freq_motif_new_after_filtering.pdf"))
par(mfrow=c(2,2),mar=c(2,2,2,2))
barplot(freq1[,2],main="within 25 nt")
mtext(side=2,line=2,text="fraction 3'UTR isoforms with at least one pA element")
barplot(freq2[,2],main="within 50 nt")
mtext(side=2,line=2,text="fraction 3'UTR isoforms with at least one pA element")
barplot(freq3[,2],main="within 100 nt")
mtext(side=2,line=2,text="fraction 3'UTR isoforms with at least one pA element")
barplot(freq4[,2],main="within 150 nt")
mtext(side=2,line=2,text="fraction 3'UTR isoforms with at least one pA element")

par(mfrow=c(2,5),mar=c(2,2,2,2))
for(i in c(1:length(motif))){
  barplot(freqall[,grep(pattern="freq",colnames(freqall))[i]])
  mtext(side=3,line=0,text=unlist(lapply(motif,toString))[i])
}

par(mfrow=c(2,2),mar=c(2,2,2,2))
for(i in c(1:4)){
  barplot(freqall[i,grep(pattern="freq",colnames(freqall))])
  mtext(side=3,line=0,text=rownames(freqall)[i])
}

par(mfrow=c(1,1),mar=c(2,2,2,2))
barplot(freqall[c(3,4),grep(pattern="freq",colnames(freqall))],beside=T,cex.names=0.5)

dev.off()


MyTestFrac <- function(n1,n2,p1,p2){
  p  <- (p1*n1+p2*n2)/(n1+n2)
  SE <- sqrt(p*(1-p)*(1/n1+1/n2))
  Z  <- (p1-p2)/SE
  PVal<- 2*pnorm(-abs(Z))#2sided
  return(list(Z,PVal))
}


no.new.cat1                <- sum(allMotifs[!anno_ngf$conservative&abs(anno_ngf$dL)<150&abs(anno_ngf$dL)>50])
freq.new.cat1              <- no.new.cat1/sum(!anno_ngf$conservative&abs(anno_ngf$dL)<150&abs(anno_ngf$dL)>50)
no.new.cat2                <- sum(allMotifs[!anno_ngf$conservative&abs(anno_ngf$dL)>150])
freq.new.cat2              <- no.new.cat2/sum(!anno_ngf$conservative&abs(anno_ngf$dL)>150)

my.no <- c(sum(!anno_ngf$conservative&abs(anno_ngf$dL)<150&abs(anno_ngf$dL)>50),sum(!anno_ngf$conservative&abs(anno_ngf$dL)>150))
my.p  <-freqall[c(3,4),grep(pattern="freq",colnames(freqall))]

z.val     <- p.val <-vector(length=length(motif))
for(i in c(1:length(motif))){
  z.val[i] <- MyTestFrac(n1=my.no[1],n2=my.no[2],p1=my.p[1,i],p2=my.p[2,i])[[1]]
  p.val[i] <- MyTestFrac(n1=my.no[1],n2=my.no[2],p1=my.p[1,i],p2=my.p[2,i])[[2]]
}
names(z.val)<-names(p.val)<-unlist(lapply(motif,toString))




#
# C. Motifs analysis now that the phasing is correct
#
IsWithinDist <- function(mydist,subRes=newRes){
  mytest <- list()
  for(i in c(1:length(subRes))){
    mytest[[i]]<- unlist(lapply(subRes[[i]]$dist_from_end,function(x)return(sum(x<=mydist&x>=(-mydist/2))>0)))
  }
  return(mytest)
}
IsWithin                <- lapply(c(25,50,100),function(X)return(do.call(IsWithinDist(X),what=cbind)))
iso_cons_merged         <- anno_ngf$iso_cons
iso_cons_merged         <- sapply(iso_cons_merged,function(X)return(gsub(X,pattern="[s,ss]",repl="")))
iso_cons_merged[which(!iso_cons_merged%in%c("-1","0","1","2","3","4"))]<-NA
anno_ngf$iso_cons_merged <- iso_cons_merged

iso_cons_focus          <- anno_ngf$iso_cons
for(i in c("-1s","0s","1s","2s","3s","4s")){
  iso_cons_focus[setdiff(grep(iso_cons_focus,pattern=i),grep(iso_cons_focus,pattern=paste("-",i,sep="")))] <- i
}
iso_cons_focus[setdiff(c(1:length(iso_cons_focus)),unlist(lapply(X=c("-1","0","1","2","3","4"),
                                                                 FUN=function(X)return(grep(iso_cons_focus,pattern=X)))))]<-NA
anno_ngf$iso_cons_focus  <- iso_cons_focus

Cols                     <- rainbow(length(levels(as.factor(iso_cons_merged))))
names(Cols)              <- levels(as.factor(iso_cons_merged))
myCols                   <- list()
myCols[[1]]              <- Cols
myCols[[2]]              <-  unlist(lapply(X=levels(factor(as.character(iso_cons_focus),
                                                           levels=c("-1s","-1", "0s","0", "1s","1","2s","2","3s","3","4s","4"))),
                                           FUN=function(x)return(Cols[which(names(Cols)==gsub(x,pattern="s",repl=""))])
)
)

### A. No difference between I/Is; focus until I4
torm              <- which(is.na(anno_ngf$iso_cons_merged))
FACT              <- factor(as.character(anno_ngf$iso_cons_merged[-torm]),levels=c("-1","0","1","2","3","4"))
no.per.cat        <- as.data.frame(table(FACT))[,2]
with.motif        <- lapply(IsWithin, function(X)return(apply(X,2,function(z)return(tapply(X=z[-torm],INDEX=FACT,FUN=function(x)return(sum(x)))))))
freq.motif        <- lapply(with.motif,function(X){
  for(i in c(1:ncol(X))){
    X[,i]  <- X[,i]/no.per.cat
  }
  colnames(X)<- lapply(motif,toString)
  return(X)}
)
freq.motif.var    <- lapply(with.motif,function(X){
  for(i in c(1:ncol(X))){
    X[,i]  <- X[,i]/no.per.cat
    X[,i]  <- X[,i]*(1-X[,i])/no.per.cat
    
  }
  colnames(X)<- lapply(motif,toString)
  return(X)}
)

Ntest            <- length(myGR)-length(torm)
freq.avg         <- lapply(freq.motif,function(X)return(apply(X,2,mean)))
freq.sd          <- lapply(freq.motif.var,function(X)return(apply(X,2,function(z)return(sum(z)/length(z)^2))))

ChiTest <- function(MOTIF=1,REGION=1){
  par(mar=c(3,3.5,2,1))
  dat            <- with.motif[[REGION]][,MOTIF]
  data           <- cbind(no.per.cat-dat,dat)
  dimnames(data) <- list("isoforms"=names(dat), "outcome"=c("w/0 motif","with motif"))
  result         <- chisq.test(data)$p.value
  Num            <- freq.motif[[REGION]][,MOTIF] - freq.avg[[REGION]][MOTIF]
  Denom          <- sqrt(freq.sd[[REGION]][MOTIF]+freq.motif.var[[REGION]][,MOTIF])
  
  
  Z.test        <- Num/Denom
  
  barplot(t(data), beside = TRUE, col = c("gray", "green"), las=1,cex.axis=0.5,cex=0.5)
  mtext(side=2,line=2.5,text="no.3'UTR isoforms",cex=0.5)
  mtext(side=3,line=0,text=paste("Chist-test P-value=",format(result,digits=2, scientific=TRUE)),cex=0.5)
  mtext(side=3,line=1,text=toString(motif[MOTIF]),cex=0.5)
  par(mar=c(0,0,0,0))
  pie(dat,col=Cols,labels=NA)
  par(mar=c(3,3,2,1))
  barplot(freq.motif[[REGION]][,MOTIF],col=Cols,las=1,cex=0.5,cex.axis=0.5)
  abline(h=freq.avg[[REGION]][MOTIF],col="black",lwd=1,lty=2)
  mtext(side=2,line=3,text="fraction with motif",cex=0.5)
  barplot(Z.test,col=Cols,las=1,cex=0.5,cex.axis=0.5)
  abline(h=c(-2.55,2.55),col="red",lty=2)
  mtext(side=2,line=2,text="z.value of enrichment",cex=0.5)
  return(result)
}

pdf(paste(outdir,"8-enrichment_motif_25nt_window.pdf",sep=""))
par(mfrow=c(6,4),mar=c(3,3.5,2,1))
for(i in c(1:12)){
  ChiTest(MOTIF=i,REGION=1)
}
dev.off()

pdf(paste(outdir,"8-enrichment_motif_50nt_window.pdf",sep=""))
par(mfrow=c(6,4),mar=c(3,3.5,2,1))
for(i in c(1:12)){
  ChiTest(MOTIF=i,REGION=2)
}
dev.off()

pdf(paste(outdir,"8-enrichment_motif_100nt_window.pdf",sep=""))
par(mfrow=c(6,4),mar=c(3,3.5,2,1))
for(i in c(1:12)){
  ChiTest(MOTIF=i,REGION=3)
}
dev.off()


myHeatmapMotifs <- function(mat=t(freq.motif[[1]]),centerCol=TRUE,centerRow=FALSE){
  
  if(centerCol){
    meanCol <- apply(mat,2,mean)
    sdCol   <- apply(mat,2,sd)
    newmat  <- mat
    for(i in c(1:ncol(mat))){newmat[,i]<- (mat[,i]-meanCol[i])/sdCol[i]}
    mat <- newmat
  }
  if(centerRow){
    meanRow <- apply(mat,1,mean)
    sdRow   <- apply(mat,1,sd)
    newmat  <- mat
    for(i in c(1:nrow(mat))){newmat[i,]<- (mat[i,]-meanRow[i])/sdRow[i]}
    mat <- newmat
  }
  
  require(gplots)
  mycols              <- greenred(255)
  
  mycols           <- colorRampPalette( c("cyan", "black", "magenta"), space="rgb")(31)
  
  b                   <- c(seq((-max(abs(mat))),0,length=16),seq(0,max(abs(mat)),length=16))
  heatmap.2(as.matrix(mat),keysize=1,col=mycols,breaks=b, scale="none",
            Rowv=TRUE,Colv=FALSE,dendr="row",key=TRUE,symkey=FALSE, density.info="none", trace="none",cexCol=1.0, cexRow=1.0, font.lab=1,las=2)
}

pdf(paste(outdir,"9-heatmap_motif_enrich_in_cat.pdf",sep=""))
myHeatmapMotifs(mat=t(freq.motif[[1]]),centerCol=FALSE,centerRow=TRUE)
myHeatmapMotifs(mat=t(freq.motif[[2]]),centerCol=FALSE,centerRow=TRUE)
myHeatmapMotifs(mat=t(freq.motif[[3]]),centerCol=FALSE,centerRow=TRUE)
dev.off()


### B. When considering the super short
torm              <- which(is.na(anno_ngf$iso_cons_focus))
FACT              <- factor(as.character(anno_ngf$iso_cons_focus[-torm]),levels=c("-1s","-1", "0s","0", "1s","1","2s","2","3s","3","4s","4"))
no.per.cat        <- as.data.frame(table(FACT))[,2]
with.motif        <- lapply(IsWithin, function(X)return(apply(X,2,function(z)return(tapply(X=z[-torm],INDEX=FACT,FUN=function(x)return(sum(x)))))))
freq.motif        <- lapply(with.motif,function(X){
  for(i in c(1:ncol(X))){
    X[,i]  <- X[,i]/no.per.cat
  }
  colnames(X)<- lapply(motif,toString)
  return(X)}
)
freq.motif.var    <- lapply(with.motif,function(X){
  for(i in c(1:ncol(X))){
    X[,i]  <- X[,i]/no.per.cat
    X[,i]  <- X[,i]*(1-X[,i])/no.per.cat
    
  }
  colnames(X)<- lapply(motif,toString)
  return(X)}
)

Ntest            <- length(myGR)-length(torm)
freq.avg         <- lapply(freq.motif,function(X)return(apply(X,2,mean)))
freq.sd          <- lapply(freq.motif.var,function(X)return(apply(X,2,function(z)return(sum(z)/length(z)^2))))

ChiTest <- function(MOTIF=1,REGION=1){
  par(mar=c(3,3.5,2,1))
  dat            <- with.motif[[REGION]][,MOTIF]
  data           <- cbind(no.per.cat-dat,dat)
  dimnames(data) <- list("isoforms"=names(dat), "outcome"=c("w/0 motif","with motif"))
  result         <- chisq.test(data)$p.value
  Num            <- freq.motif[[REGION]][,MOTIF] - freq.avg[[REGION]][MOTIF]
  Denom          <- sqrt(freq.sd[[REGION]][MOTIF]+freq.motif.var[[REGION]][,MOTIF])
  Z.test        <- Num/Denom
  
  barplot(t(data), beside = TRUE, col = c("gray", "green"), las=1,cex.axis=0.5,cex=0.5)
  mtext(side=2,line=2.5,text="no.3'UTR isoforms",cex=0.5)
  mtext(side=3,line=0,text=paste("Chist-test P-value=",format(result,digits=2, scientific=TRUE)),cex=0.5)
  mtext(side=3,line=1,text=toString(motif[MOTIF]),cex=0.5)
  par(mar=c(0,0,0,0))
  pie(dat,col=myCols[[2]],labels=NA)
  par(mar=c(3,3,2,1))
  barplot(freq.motif[[REGION]][,MOTIF],col=myCols[[2]],las=1,cex=0.5,cex.axis=0.5)
  abline(h=freq.avg[[REGION]][MOTIF],col="black",lwd=1,lty=2)
  mtext(side=2,line=3,text="fraction with motif",cex=0.5)
  barplot(Z.test,col=myCols[[2]],las=1,cex=0.5,cex.axis=0.5)
  abline(h=c(-2.55,2.55),col="red",lty=2)
  mtext(side=2,line=2,text="z.value of enrichment",cex=0.5)
  return(result)
}

pdf(paste(outdir,"8-enrichment_motif_25nt_window_ws.pdf",sep=""))
par(mfrow=c(6,4),mar=c(3,3.5,2,1))
for(i in c(1:12)){
  ChiTest(MOTIF=i,REGION=1)
}
dev.off()

pdf(paste(outdir,"8-enrichment_motif_50nt_window_ws.pdf",sep=""))
par(mfrow=c(6,4),mar=c(3,3.5,2,1))
for(i in c(1:12)){
  ChiTest(MOTIF=i,REGION=2)
}
dev.off()

pdf(paste(outdir,"8-enrichment_motif_100nt_window_ws.pdf",sep=""))
par(mfrow=c(6,4),mar=c(3,3.5,2,1))
for(i in c(1:12)){
  ChiTest(MOTIF=i,REGION=3)
}
dev.off()


myHeatmapMotifs <- function(mat=t(freq.motif[[1]]),centerCol=TRUE,centerRow=FALSE){
  
  if(centerCol){
    meanCol <- apply(mat,2,mean)
    sdCol   <- apply(mat,2,sd)
    newmat  <- mat
    for(i in c(1:ncol(mat))){newmat[,i]<- (mat[,i]-meanCol[i])/sdCol[i]}
    mat <- newmat
  }
  if(centerRow){
    meanRow <- apply(mat,1,mean)
    sdRow   <- apply(mat,1,sd)
    newmat  <- mat
    for(i in c(1:nrow(mat))){newmat[i,]<- (mat[i,]-meanRow[i])/sdRow[i]}
    mat <- newmat
  }
  
  require(gplots)
  mycols              <- greenred(255)
  
  mycols           <- colorRampPalette( c("cyan", "black", "magenta"), space="rgb")(31)
  
  b                   <- c(seq((-max(abs(mat))),0,length=16),seq(0,max(abs(mat)),length=16))
  heatmap.2(as.matrix(mat),keysize=1,col=mycols,breaks=b, scale="none",
            Rowv=TRUE,Colv=FALSE,dendr="row",key=TRUE,symkey=FALSE, density.info="none", trace="none",cexCol=1.0, cexRow=1.0, font.lab=1,las=2)
}

pdf(paste(outdir,"9-heatmap_motif_enrich_in_cat_ws.pdf",sep=""))
myHeatmapMotifs(mat=t(freq.motif[[1]]),centerCol=FALSE,centerRow=TRUE)
myHeatmapMotifs(mat=t(freq.motif[[2]]),centerCol=FALSE,centerRow=TRUE)
myHeatmapMotifs(mat=t(freq.motif[[3]]),centerCol=FALSE,centerRow=TRUE)
dev.off()




### C. Tset differences between short and long version


MyTestFrac <- function(n1,n2,p1,p2){
  p  <- (p1*n1+p2*n2)/(n1+n2)
  SE <- sqrt(p*(1-p)*(1/n1+1/n2))
  Z  <- (p1-p2)/SE
  PVal<- 2*pnorm(-abs(Z))#2sided
  return(list(Z,PVal))
}




Ntest   <- length(myGR)#n2
ptot25   <- apply(IsWithin25,2,sum)/Ntest
ptot50   <- apply(IsWithin50,2,sum)/Ntest
ptot100  <- apply(IsWithin100,2,sum)/Ntest

MyTestFrac <- function(n1,n2,p1,p2){
  p  <- (p1*n1+p2*n2)/(n1+n2)
  SE <- sqrt(p*(1-p)*(1/n1+1/n2))
  Z  <- (p1-p2)/SE
  PVal<- 2*pnorm(-abs(Z))#2sided
  return(list(Z,PVal))
}

z.val25      <- p.val25 <-matrix(ncol=length(motif),nrow=nrow(freq.motif25))
z.val50      <- p.val50 <-matrix(ncol=length(motif),nrow=nrow(freq.motif25))
z.val100     <- p.val100 <-matrix(ncol=length(motif),nrow=nrow(freq.motif25))

for(j in c(1:length(motif))){
  for(i in c(1:nrow(freq.motif25))){
    z.val25[i,j]  <- MyTestFrac(n1=no.per.cat[i],n2=Ntest,p1=freq.motif25[i,j],p2=ptot25[j])[[1]]
    p.val25[i,j]  <- -log10(MyTestFrac(n1=no.per.cat[i],n2=Ntest,p1=freq.motif25[i,j],p2=ptot25[j])[[2]])
    z.val50[i,j]  <- MyTestFrac(n1=no.per.cat[i],n2=Ntest,p1=freq.motif50[i,j],p2=ptot50[j])[[1]]
    p.val50[i,j]  <- -log10(MyTestFrac(n1=no.per.cat[i],n2=Ntest,p1=freq.motif50[i,j],p2=ptot50[j])[[2]])
    z.val100[i,j] <- MyTestFrac(n1=no.per.cat[i],n2=Ntest,p1=freq.motif100[i,j],p2=ptot100[j])[[1]]
    p.val100[i,j] <- -log10(MyTestFrac(n1=no.per.cat[i],n2=Ntest,p1=freq.motif100[i,j],p2=ptot100[j])[[2]])
  }
}

colnames(z.val25)<-colnames(z.val50)<-colnames(z.val100)<-colnames(p.val25)<-colnames(p.val50)<-colnames(p.val100)<- motif
rownames(z.val25)<-rownames(z.val50)<-rownames(z.val100)<-rownames(p.val25)<-rownames(p.val50)<-rownames(p.val100)<- rownames(freq.motif25)
ord  <- match(c("-1ss","-1s","-1",
                "0ss","0s","0",
                "1ss","1s","1",
                "2ss","2s","2",
                "3ss","3s","3",
                "4ss","4s","4",
                "5ss","5s","5",
                "6ss","6s","6",
                "7ss","7s","7",
                "8ss","8s","8"),rownames(z.val50))
ord <- ord[!is.na(ord)]
p.val25 <- p.val25[ord,]
p.val50 <- p.val50[ord,]
p.val100 <- p.val100[ord,]
no.per.cat <- no.per.cat[ord]
names(no.per.cat)<- rownames(p.val25)



