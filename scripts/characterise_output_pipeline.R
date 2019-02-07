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

# A. Load data
load("./annotation/rn5/GRanges_comprehensive_transcriptome_rat_24_nov_2015.RData")
myUTR        <- import.gff("./annotation/rn5/Lngf.gtf",format="gtf")
anno_ngf     <- read.table("./data/anno_ngf.tab",header=T,sep="\t")

# B. Characterise pipeline in terms of how many new annotated
no.txID.rn5     <- length(unique(mcols(g3utr)$X.transcript_id.))#22'845 as those for which 3' UTR length was zero were removed
no.isoforms     <- nrow(anno_ngf)
is.I0           <- tapply(anno_ngf$isoform,INDEX=anno_ngf$txID,function(Z)return(sum(Z==0)==0))
no.new.isoform  <- sum(anno_ngf$isoform!=0)#26'468
no.new.longer   <- sum(anno_ngf$newL-anno_ngf$initL>=30)
no.new.shorter  <- sum(-anno_ngf$newL+anno_ngf$initL>=30)
focus.iso       <- c(no.new.isoform,no.new.longer,no.new.shorter)

no.txID.modif   <- sum(tapply(abs(anno_ngf$newL-anno_ngf$initL),INDEX=anno_ngf$txID,function(Z)return(sum(Z>30)>0)))#10219 txID modified
no.txID.extended<- sum(unlist(lapply((anno_ngf$maxL-anno_ngf$initL)[match(unique(anno_ngf$txID),anno_ngf$txID)],function(Z)return(Z>30))))#7506 txID extended
no.txID.shortened<- sum(unlist(lapply((anno_ngf$initL-anno_ngf$minL)[match(unique(anno_ngf$txID),anno_ngf$txID)],function(Z)return(Z>30))))#4721 txID shortened
no.txID.both    <- sum(
  unlist(lapply((anno_ngf$initL-anno_ngf$minL)[match(unique(anno_ngf$txID),anno_ngf$txID)],function(Z)return(Z>30)))&
  unlist(lapply((anno_ngf$maxL-anno_ngf$initL)[match(unique(anno_ngf$txID),anno_ngf$txID)],function(Z)return(Z>30))))#2'008 txID extended
focus.txID       <- c(no.txID.modif,no.txID.extended,no.txID.shortened,no.txID.both)


# A. Figure 1C-D
# The analysis of the frequancy of motifs, overlap with Xeno and so on is performed on this one; however then we will focus our analysos on isoform separated by at least 500nt (Lngf_sub)
myno.new       <- list( no.novel = nrow(anno_ngf)-length(grep(anno_ngf$uniqueID,pattern="\\.0")),
                        no.extended= sum(anno_ngf$newL-anno_ngf$initL>0),
                        no.shortened= sum(anno_ngf$newL-anno_ngf$initL<0))
myno.new.tX   <- list( no.novel     = length(unique(as.character(anno_ngf$txID)[setdiff(c(1:nrow(anno_ngf)),grep(anno_ngf$uniqueID,pattern="\\.0"))])),
                       no.extended  = length(unique(as.character(anno_ngf$txID)[which(anno_ngf$newL-anno_ngf$initL>0)])),
                       no.shortened = length(unique(as.character(anno_ngf$txID)[which(anno_ngf$newL-anno_ngf$initL<0)])),
                       no.both      = length(intersect(unique(as.character(anno_ngf$txID)[which(anno_ngf$newL-anno_ngf$initL<0)]),
                                                       unique(as.character(anno_ngf$txID)[which(anno_ngf$newL-anno_ngf$initL>0)]))))

#Supplementary Figure 2b
par(mfrow=c(1,2))
mp<-barplot(unlist(myno.new),ylim=c(0,30000),las=1,col="white",cex.axis=0.7,cex.lab=0.7,cex.names=0.7,ylab="")
mtext(side=3,line=0,text=unlist(myno.new),at=mp,cex=0.5)
mtext(side=2,line=2,text="#newly annotated 3' UTR isoforms",cex=0.7)
mp<-barplot(unlist(myno.new.tX)[-1],ylim=c(0,8000),las=1,col="white",cex.axis=0.7,cex.lab=0.7,cex.names=0.7,ylab="")
mtext(side=3,line=0,text=unlist(myno.new.tX)[-1],at=mp,cex=0.5)
mtext(side=2,line=2,text="#newly Ensembl txID with novel 3'UTR",cex=0.7)



# C. Remove redundant transcripts and select reliably expressed genes
myUTR    <- import.gff("./annotation/rn5/Lngf_sub.gtf",format="gtf")
id       <- paste(start(myUTR),end(myUTR),strand(myUTR),sep=".")#2215 duplicated
myUTR    <- myUTR[-which(duplicated(id)),]
anno_ngf <- anno_ngf[match(myUTR$ID,anno_ngf$uniqueID),]

#Figure S3c
par(mfrow=c(2,2))
mydat                <- anno_ngf[which(anno_ngf$is.conservative),match(c("NGF.axon.1.raw","NGF.axon.2.raw","NGF.cb.1.raw","NGF.cb.2.raw"),colnames(anno_ngf))]
mydatall             <- anno_ngf[,match(c("NGF.axon.1.raw","NGF.axon.2.raw","NGF.cb.1.raw","NGF.cb.2.raw"),colnames(anno_ngf))]
for(i in c(1:ncol(mydat))){
  mydat[,i] <- log2(mydat[,i]+1)
}
lims      <- apply(mydat,2,function(Z)2^SelectExpressed(dat=Z,frac.bg=0.95,frac.fg=0.1))
tempsel   <- apply(mydatall,2,function(Z)return(Z>=mean(lims[1,c(1,2)])))
soft.sel  <- cbind(tempsel[,1]&tempsel[,2],tempsel[,3]&tempsel[,4])
tempsel   <- do.call(what=cbind,lapply(c(1:4),function(Z)return(mydatall[,Z]>=lims[2,Z])))
hard.sel  <- cbind(tempsel[,1]|tempsel[,2],tempsel[,3]|tempsel[,4])
final.sel <- cbind(soft.sel[,1]|hard.sel[,1],soft.sel[,2]|hard.sel[,2])

NGF.axon.is.expressed.iso                   <- final.sel[,1]
NGF.cb.is.expressed.iso                     <- final.sel[,2]
anno_ngf$NGF.axon.is.expressed.iso          <- final.sel[,1]
anno_ngf$NGF.cb.is.expressed.iso            <- final.sel[,2]

#FG limit adapted to all samples
#BG limit from axonal samples (better fit)

# D. Characterise the expression in each compartment in terms of the numer of expressed transcripts
axons                                              <- unique(anno_ngf$txID[anno_ngf$NGF.axon.is.expressed.iso])
cb                                                 <- unique(anno_ngf$txID[anno_ngf$NGF.cb.is.expressed.iso])
neurons                                            <- union(axons,cb)
cb.only                                            <- setdiff(cb,axons)
axon.only                                          <- setdiff(axons,cb)
both                                               <- intersect(axons,cb)
detect.txID.ngf                                    <- rep("no",nrow(anno_ngf))
detect.txID.ngf[which(anno_ngf$txID%in%both)]      <- "both"
detect.txID.ngf[which(anno_ngf$txID%in%cb.only)]   <- "cb.only"
detect.txID.ngf[which(anno_ngf$txID%in%axon.only)] <- "axon.only"
anno_ngf$detect.txID.ngf                           <- detect.txID.ngf
#uniqueID.expression
axons                                             <- anno_ngf$uniqueID[anno_ngf$NGF.axon.is.expressed.iso]
cb                                                <- anno_ngf$uniqueID[anno_ngf$NGF.cb.is.expressed.iso]
neur.neur                                         <- as.character(anno_ngf$uniqueID)[anno_ngf$uniqueID%in%intersect(cb,axons)&anno_ngf$txID%in%both]
neur.cb                                           <- as.character(anno_ngf$uniqueID)[anno_ngf$uniqueID%in%setdiff(cb,axons)&anno_ngf$txID%in%both]
neur.ax                                           <- as.character(anno_ngf$uniqueID)[anno_ngf$uniqueID%in%setdiff(axons,cb)&anno_ngf$txID%in%both]
ax.ax                                             <- as.character(anno_ngf$uniqueID)[anno_ngf$uniqueID%in%setdiff(axons,cb)&anno_ngf$txID%in%axon.only]
cb.cb                                             <- as.character(anno_ngf$uniqueID)[anno_ngf$uniqueID%in%setdiff(cb,axons)&anno_ngf$txID%in%cb.only]

detect.iso.ngf                                    <- rep("no",nrow(anno_ngf))
detect.iso.ngf[which(anno_ngf$txID%in%neur.neur)] <- "neur.neur"
detect.iso.ngf[which(anno_ngf$txID%in%neur.cb)]   <- "neur.cb"
detect.iso.ngf[which(anno_ngf$txID%in%neur.ax)]   <- "neur.ax"
detect.iso.ngf[which(anno_ngf$txID%in%ax.ax)]     <- "ax.ax"
detect.iso.ngf[which(anno_ngf$txID%in%cb.cb)]     <- "cb.cb"
anno_ngf$detect.iso.ngf                           <- detect.iso.ngf


no.txID <- c(cb=length(unique(anno_ngf$txID[anno_ngf$NGF.cb.is.expressed.iso])),
             axons=length(unique(anno_ngf$txID[anno_ngf$NGF.axon.is.expressed.iso]) ))
no.iso  <- c(cb=length(unique(anno_ngf$uniqueID[anno_ngf$NGF.cb.is.expressed.iso])),
             axons=length(unique(anno_ngf$uniqueID[anno_ngf$NGF.axon.is.expressed.iso]) ))

#Figure S3d
mycols   <- c(rgb(23/255,71/255,120/255),rgb(119/255,192/255,68/255))
par(mfrow=c(1,2))
mp<- barplot(no.iso,las=1,col=mycols)
mtext(side=3,line=0,text=no.iso,cex=0.7,at=mp)
mp<- barplot(no.txID,las=1,col=mycols)
mtext(side=3,line=0,text=no.txID,cex=0.7,at=mp)

# E. Characterise the expression in each compartment in terms of the enrichment in GO terms
t2g                     <- read.csv("./annotation/rn5/t2g_biomaRt.csv")
txID2GO                 <- tapply(t2g$ensembl_transcript_id,INDEX=t2g$go_id,FUN=function(x)return(x))
noNodes                 <- 300

myInterestingGenes     <- list(cb=unique(anno_ngf$txID[anno_ngf$NGF.cb.is.expressed.iso]),
                               axons=unique(anno_ngf$txID[anno_ngf$NGF.axon.is.expressed.iso]),
                               neurons=union(unique(anno_ngf$txID[anno_ngf$NGF.cb.is.expressed.iso]),unique(anno_ngf$txID[anno_ngf$NGF.axon.is.expressed.iso])),
                               cb.only=setdiff(unique(anno_ngf$txID[anno_ngf$NGF.cb.is.expressed.iso]),unique(anno_ngf$txID[anno_ngf$NGF.axon.is.expressed.iso]))
                               )
mysampleGO           <- lapply(myInterestingGenes,function(Z){
  geneNames              <- unique(t2g$ensembl_transcript_id)
  geneList               <- factor(as.integer(geneNames %in% Z))
  names(geneList)        <- geneNames
  return(new("topGOdata",description = "Simple session", ontology = "BP",allGenes = geneList, geneSel = Z,nodeSize = 10,annot = annFUN.GO2genes,GO2gene=txID2GO))
                        })
myenrich            <- lapply(mysampleGO,function(Z){
  resultFisher.weight01   <- runTest(Z, algorithm = "weight01", statistic = "fisher")
  allRes                 <- GenTable(Z, weight0Fisher=resultFisher.weight01,orderBy = "weight0Fisher", ranksOf = "weight0Fisher", topNodes = 100)
  return(list(allRes))})


outdir <- "./GOenrichment/"

write.csv(myenrich[[1]],paste(outdir,"enrichment_cb.csv",sep=""))
write.csv(myenrich[[2]],paste(outdir,"enrichment_axons.csv",sep=""))
write.csv(myenrich[[3]],paste(outdir,"enrichment_neurons.csv",sep=""))
write.csv(myenrich[[4]],paste(outdir,"enrichment_cb.only.csv",sep=""))

#Filter GO to remove redundant terms
enr1<- read.csv(paste(outdir,"enrichment_neurons_f.csv",sep=""))
enr2<- read.csv(paste(outdir,"enrichment_axons_f.csv",sep=""))
enr3<- read.csv(paste(outdir,"enrichment_cb.only_f.csv",sep=""))

#Figure S3d
par(mfrow=c(1,3))
vec<- -log10(enr1$weight0Fisher)
names(vec)<- as.character(enr1$Term)
barplot(rev(vec),horiz=T,cex.names=0.6,las=1)
mtext(side=3,line=0,text="neurones")
mtext(side=1,line=2,text="-log10(P-Value)")

vec<- -log10(enr3$weight0Fisher)
names(vec)<- as.character(enr3$Term)
barplot(rev(vec),horiz=T,cex.names=0.6,las=1)
mtext(side=3,line=0,text="cell body only")
mtext(side=1,line=2,text="-log10(P-Value)")

vec<- -log10(enr2$weight0Fisher)
names(vec)<- as.character(enr2$Term)
barplot(rev(vec),horiz=T,cex.names=0.6,las=1)
mtext(side=3,line=0,text="axons and cell body")
mtext(side=1,line=2,text="-log10(P-Value)")

# F. Characterise the expression in each compartment in terms of the 3' UTR length
Axons.txID <- unique(as.character(anno_ngf$txID)[anno_ngf$NGF.axon.is.expressed.iso])
CB.txID    <- unique(as.character(anno_ngf$txID)[anno_ngf$NGF.cb.is.expressed.iso])
Axons.iso  <- unique(as.character(anno_ngf$uniqueID)[anno_ngf$NGF.axon.is.expressed.iso])
CB.iso     <- unique(as.character(anno_ngf$uniqueID)[anno_ngf$NGF.cb.is.expressed.iso])

L1 <- list(
  cb.only=setdiff(CB.txID,Axons.txID),
  cb.axons=intersect(CB.txID,Axons.txID),
  axon.only=setdiff(Axons.txID,CB.txID)
)

L2 <- list(
  cb.only=setdiff(CB.iso,Axons.iso),
  cb.axons=intersect(CB.iso,Axons.iso),
  axon.only=setdiff(Axons.iso,CB.iso)
)


maxL        <- tapply(anno_ngf$newL[anno_ngf$NGF.axon.is.expressed.iso|anno_ngf$NGF.cb.is.expressed.iso],INDEX=factor(as.character(anno_ngf$txID[anno_ngf$NGF.axon.is.expressed.iso|anno_ngf$NGF.cb.is.expressed.iso])),FUN=max)
initL       <- tapply(anno_ngf$initL[anno_ngf$NGF.axon.is.expressed.iso|anno_ngf$NGF.cb.is.expressed.iso],INDEX=factor(as.character(anno_ngf$txID[anno_ngf$NGF.axon.is.expressed.iso|anno_ngf$NGF.cb.is.expressed.iso])),FUN=function(X)return(X[1]))
no.iso      <- tapply(anno_ngf$initL[anno_ngf$NGF.axon.is.expressed.iso|anno_ngf$NGF.cb.is.expressed.iso],INDEX=factor(as.character(anno_ngf$txID[anno_ngf$NGF.axon.is.expressed.iso|anno_ngf$NGF.cb.is.expressed.iso])),FUN=function(X)return(length(X)))

anno_ngf$maxL  <- maxL[match(anno_ngf$txID,names(maxL))]
anno_ngf$initL <- initL[match(anno_ngf$txID,names(initL))]
anno_ngf$no.iso<- no.iso[match(anno_ngf$txID,names(no.iso))]
with.multiple <- unlist(lapply(L1,function(Z)return(sum(no.iso[match(Z,names(no.iso))]>1)/length(Z))))[-3]
distr.no.iso  <- rbind(
              as.vector(table(no.iso[match(L1[[1]],names(no.iso))]))/length(L1[[1]]),
              (as.vector(table(no.iso[match(L1[[2]],names(no.iso))]))[c(1:8)])/length(L1[[2]]))


utrL          <- list(
  Ensembl=initL[match(unique(c(Axons.txID,CB.txID)),names(initL))],
  cb.only=maxL[match(L1[[1]],names(maxL))],
  axons=maxL[match(L1[[2]],names(maxL))]
)

#Figures S3e
par(mfrow=c(2,3),mar=c(3,4,3,3))
mycols <- c(rgb(254/255,218/255,0),rgb(150/255,150/255,150/255),"black")
mp<- barplot(unlist(lapply(L1,length)),las=1,col=mycols,frame=F)
mtext(side=3,line=0,text=unlist(lapply(L1,length)),at=mp,cex=0.6)
mtext(side=2,line=3,text="no.txID")

mycols <- c(rgb(254/255,218/255,0),rgb(150/255,150/255,150/255),"black")
mp<- barplot(unlist(lapply(L2,length)),las=1,col=mycols,frame=F)
mtext(side=3,line=0,text=unlist(lapply(L2,length)),at=mp,cex=0.5)
mtext(side=2,line=3,text="no.3' UTR isoforms")

#Figures 2c,d
barplot(with.multiple*100,col=mycols,las=1,frame=F)
mtext(side=2,line=3,text="fraction of txID with multiple iso")
barplot(distr.no.iso*100,beside=T,col=mycols[c(1,2)],las=1,,frame=F)
mtext(side=2,line=3,text="fraction of txID with multiple iso")
boxplot(utrL,outline=F,col=c("white",mycols[c(1,2)]),las=1,frame=F)
mtext(side=3,line=0,text=unlist(lapply(utrL,length)),at=c(1,2,3),cex=0.5)


#Statistical test to to test whether the difference in 3'UTR length is significant
wilcox.test(utrL[[2]],utrL[[3]], paired=F)$p.value
#Statistical test to test whether the fraction of txID with mult iso is different
Is.axons    <- unique(union(CB.txID,Axons.txID))%in%L1[[2]]
Is.multiple <- no.iso[match(unique(union(CB.txID,Axons.txID)),names(no.iso))]>1

fisher.test(Is.axons,Is.multiple)$p.value
temp<-as.matrix(table(Is.axons,Is.multiple))
temp<-temp/length(Is.axons)


#Test of significance between 3'UTR length between the CB and axonal compartment  using GLM
myMaxL <- tapply(anno_ngf$newL[anno_ngf$NGF.axon.is.expressed.iso|anno_ngf$NGF.cb.is.expressed.iso],INDEX=factor(as.character(anno_ngf$txID[anno_ngf$NGF.axon.is.expressed.iso|anno_ngf$NGF.cb.is.expressed.iso])),FUN=max)
myComp <- as.character(unlist(tapply(anno_ngf$detect.txID.ngf[anno_ngf$NGF.axon.is.expressed.iso|anno_ngf$NGF.cb.is.expressed.iso],INDEX=as.character(anno_ngf$txID)[anno_ngf$NGF.axon.is.expressed.iso|anno_ngf$NGF.cb.is.expressed.iso],FUN=function(x)return(as.character(unique(x))))))
torm   <- myComp=="axon.only"
myMaxL <- myMaxL[!torm]
myComp <- myComp[!torm]

datL         <-  data.frame(length=myMaxL,compartment=myComp,lengthL=log10(myMaxL),cb=myComp=="cb.only",neur=myComp=="both")
datL         <- datL[-which(datL$compartment=="no"),]


myComp        <- c(rep("cb.only",length(L1[[1]])),rep("axons",length(L1[[2]])))
utrL <- list(
  Ensembl=initL[match(unique(c(Axons.txID,CB.txID)),names(initL))],
  cb.only=maxL[match(L1[[1]],names(maxL))],
  axons=maxL[match(L1[[2]],names(maxL))]
)
myL1          <- c(utrL[[2]],utrL[[3]])
myL2          <- log10(myL1)
datL          <-  data.frame(length=myL1,compartment=myComp,lengthL=myL2,cb=myComp=="cb.only",neur=myComp=="axons")


with(datL, tapply(lengthL, compartment, function(x) {
  sprintf("M (SD) = %1.2f (%1.2f)", mean(x), sd(x))
}))

with(datL, tapply(lengthL, compartment, function(x) {
  sprintf("Median (SD) = %1.2f (%1.2f)", median(x), sd(x))
}))

with(datL, tapply(length, compartment, function(x) {
  sprintf("Median (SD) = %1.2f (%1.2f)", median(x), sd(x))
}))

require(MASS)
summary(m1 <- glm.nb(length ~1+compartment, data = datL))
m1          <- glm.nb(length ~compartment, data = datL)
test.fitted <-  predict(m1, type = "response")
confint(m1)




