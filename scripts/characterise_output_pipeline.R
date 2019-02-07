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



### NESTED FUNCTIONS

GetOI <- function(mygoID="GO:0006412",sampleGO){
  go.genes    <- genesInTerm(sampleGO, mygoID)[[1]]#To extract all genes related to this term
  sig.genes   <- sigGenes(sampleGO)
  goi         <- intersect(sig.genes,go.genes)
  return(goi)
}

CompareBP <- function(enr1=test1,enr2=test2,PLOT=TRUE,no=10,lab1="remodelled",lab2="transport"){
  
  myterms            <- unique(c(enr1$Term,enr2$Term))
  temp1              <- -log10(enr1[match(myterms,enr1$Term),]$P.DE)
  temp2              <- -log10(enr2[match(myterms,enr2$Term),]$P.DE)
  temp1[is.na(temp1)]<- 0
  names(temp1)      <-  myterms
  temp2[is.na(temp2)]<- 0
  names(temp2)       <- myterms
  out                <- data.frame(term=myterms,val1=temp1,val2=temp2)
  
  
  if(PLOT){
    
    selterms <- unique(c(names(sort(temp1,decreasing=T))[c(1:no)],names(sort(temp2,decreasing=T))[c(1:no)]))
    dat      <- out[match(selterms,out$term),]
    dat      <- dat[sort(dat$val1-dat$val2,decreasing=T,index.return=T)$ix,]
    L        <- nrow(dat)
    val1     <- as.vector(-dat$val1)
    val2     <- as.vector(dat$val2)
    
    par(mar=c(3,10,2,3),cex=0.7)
    plot(c(L+3,0),xlim=c(min(-dat$val1)-10,max(dat$val2)),type = "n",frame=F,yaxt="n",ylab="")
    mp=barplot(height = val1,add = TRUE,axes = FALSE,horiz=T,col="lightsteelblue3")
    barplot(height = val2,add = TRUE,axes = FALSE,horiz=T,col="midnightblue")
    text(x=(min(-dat$val1)-5),y=mp,lab=as.character(dat$term),cex=0.6)
    mtext(side=1,line=2,text="-log10[P-Value]",cex=0.6)
    legend("top",pch=15,col=c("lightsteelblue3","midnightblue"),bty="n",ncol=2,leg=c("axonal remodelling","facilitated transport"),cex=0.5)
  }
  return(out)
}


CompareBP_improved <- function(enr1=enrichLong[[2]][[3]],enr2=enrichShort[[2]][[3]],PLOT=TRUE,no=10,lab1="long",lab2="short",coi="weight0Fisher"){
  
  myterms            <- unique(c(enr1$Term[enr1$Significant>=5&as.numeric(enr1[,match(coi,colnames(enr1))])<=0.05],enr2$Term[enr1$Significant>=5&as.numeric(enr2[,match(coi,colnames(enr2))])<=0.05]))
  temp1              <- -log10(as.numeric(enr1[match(myterms,enr1$Term),match(coi,colnames(enr1))]))
  temp2              <- -log10(as.numeric(enr2[match(myterms,enr2$Term),match(coi,colnames(enr2))]))
  temp1[is.na(temp1)]<- 0
  names(temp1)      <-  myterms
  temp2[is.na(temp2)]<- 0
  names(temp2)       <- myterms
  out                <- data.frame(term=myterms,val1=temp1,val2=temp2)
  
  if(PLOT){
    
    selterms <- unique(c(names(sort(temp1,decreasing=T))[c(1:no)],names(sort(temp2,decreasing=T))[c(1:no)]))
    dat      <- out[match(selterms,out$term),]
    dat      <- dat[sort(dat$val1-dat$val2,decreasing=T,index.return=T)$ix,]
    L        <- nrow(dat)
    val1     <- as.vector(-dat$val1)
    val2     <- as.vector(dat$val2)
    
    par(mar=c(3,10,2,3),cex=0.7)
    plot(c(L+3,0),xlim=c(min(-dat$val1)-10,max(dat$val2)),type = "n",frame=F,yaxt="n",ylab="")
    mp=barplot(height = val1,add = TRUE,axes = FALSE,horiz=T,col="lightsteelblue3")
    barplot(height = val2,add = TRUE,axes = FALSE,horiz=T,col="midnightblue")
    text(x=(min(-dat$val1)-5),y=mp,lab=as.character(dat$term),cex=0.6)
    mtext(side=1,line=2,text="-log10[P-Value]",cex=0.6)
    legend("top",pch=15,col=c("lightsteelblue3","midnightblue"),bty="n",ncol=2,leg=c(lab1,lab2),cex=0.5)
  }
  colnames(out)      <- c("terms",lab1,lab2)
  
  return(out)
}

#Enrichment with GOANA
GetEnrich <- function(selection=selRUD[[5]],value=dRUD){
  mydat           <- subanno[selection,]
  if(!is.na(value)){
    print("with value")
    mydat <- data.frame(mydat,t2g[match(mydat$geneSymbol,t2g$external_gene_name),],myval=value[selection])
    go    <- goana(unique(mydat$entrezgene), trend=mydat$myval, universe = t2g$entrezgene, species = "Rn", prior.prob = NULL, covariate=NULL,plot=FALSE,FDR=0.01)
    BP    <- topGO(go, ontology="BP",number=Inf)
    return(BP[BP$P.DE<=0.01,])
  }
  else{
    mydat <- data.frame(mydat,t2g[match(mydat$geneSymbol,t2g$external_gene_name),])
    go    <- goana(unique(mydat$entrezgene), trend=NULL, universe = t2g$entrezgene, species = "Rn", prior.prob = NULL, covariate=NULL,plot=FALSE,FDR=0.01)
    BP    <- topGO(go, ontology="BP",number=Inf)
    return(BP[BP$P.DE<=0.01,])
  }
}


CreateSampleGO <- function(mysel){
  sampleGO               <- list()
  
  geneNames              <- myBG$txID
  myInterestingGenes     <- unique(as.character(subanno$txID)[mysel])
  myInterestingGenes     <- myInterestingGenes[!is.na(myInterestingGenes)]
  geneList               <- factor(as.integer(geneNames %in% myInterestingGenes))
  names(geneList)        <- geneNames
  sampleGO[[1]]         <- new("topGOdata",description = "Simple session", ontology = "BP",allGenes = geneList, geneSel = myInterestingGenes,nodeSize = 10,annot = annFUN.GO2genes,GO2gene=txID2GO)
  
  
  geneNames              <- unique(t2g$ensembl_transcript_id)
  myInterestingGenes     <- unique(as.character(subanno$txID)[mysel])
  myInterestingGenes     <- myInterestingGenes[!is.na(myInterestingGenes)]
  geneList               <- factor(as.integer(geneNames %in% myInterestingGenes))
  names(geneList)        <- geneNames
  sampleGO[[2]]     <- new("topGOdata",description = "Simple session", ontology = "BP",allGenes = geneList, geneSel = myInterestingGenes,nodeSize = 10,annot = annFUN.GO2genes,GO2gene=txID2GO)
  return(sampleGO)
}


CreateSampleGOList <- function(myInterestingGenes){
  sampleGO               <- list()
  geneNames              <- myBG$txID
  geneList               <- factor(as.integer(geneNames %in% myInterestingGenes))
  names(geneList)        <- geneNames
  sampleGO[[1]]         <- new("topGOdata",description = "Simple session", ontology = "BP",allGenes = geneList, geneSel = myInterestingGenes,nodeSize = 10,annot = annFUN.GO2genes,GO2gene=txID2GO)
  
  geneNames              <- unique(t2g$ensembl_transcript_id)
  geneList               <- factor(as.integer(geneNames %in% myInterestingGenes))
  names(geneList)        <- geneNames
  sampleGO[[2]]     <- new("topGOdata",description = "Simple session", ontology = "BP",allGenes = geneList, geneSel = myInterestingGenes,nodeSize = 10,annot = annFUN.GO2genes,GO2gene=txID2GO)
  return(sampleGO)
}

#Enrichment with topGO
getEnrich     <- function(mysampleGO=sampleGOdata1){
  resultFisher            <- runTest(mysampleGO, algorithm = "classic", statistic = "fisher")
  resultFisher.weight01   <- runTest(mysampleGO, algorithm = "weight01", statistic = "fisher")
  
  allRes1.1                 <- GenTable(mysampleGO, classicFisher = resultFisher,weight0Fisher=resultFisher.weight01,orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = noNodes)
  allRes1.2                 <- GenTable(mysampleGO, classicFisher = resultFisher, weight0Fisher=resultFisher.weight01,orderBy = "weight0Fisher", ranksOf = "weight0Fisher", topNodes = noNodes)
  return(list(allRes1.1,allRes1.2))
}

plotEnrich1 <- function(dat=elim){
  dat$pval.unique <- -log10(dat$pval.unique)
  dat$pval.pluri  <- -log10(dat$pval.pluri)
  dat[is.na(dat)] <- 0
  dat             <- dat[sort(dat$pval.unique-dat$pval.pluri,decreasing=T,index.return=T)$ix,]
  L               <- nrow(dat)
  val1            <- as.vector(-dat$pval.unique)
  val2            <- as.vector(dat$pval.pluri)
  
  par(mar=c(3,10,2,3),cex=0.7)
  plot(c(L+3,0),xlim=c(min(-dat$pval.unique)-10,max(dat$pval.pluri)),type = "n",frame=F,yaxt="n",ylab="")
  mp=barplot(height = val1,add = TRUE,axes = FALSE,horiz=T,col="lightsteelblue3")
  barplot(height = val2,add = TRUE,axes = FALSE,horiz=T,col="midnightblue")
  text(x=(min(-dat$pval.unique)-5),y=mp,lab=as.character(dat$Term),cex=0.6)
  mtext(side=1,line=2,text="-log10[P-Value]",cex=0.6)
  legend("top",pch=15,col=c("lightsteelblue3","midnightblue"),bty="n",ncol=2,leg=c("#iso=1","#iso>1"),cex=0.5)
}

MyEnrichPlot <- function(dat=wF,mytitle="BP"){
  
  par(mar=c(3,10,2,3),cex=0.7)
  dat                       <- dat[sort(dat[,1]-dat[,2],decreasing=T,index.return=T)$ix,]
  L                         <- nrow(dat)
  
  val1 <- as.vector(-dat[,1])
  val2 <- as.vector(dat[,2])
  
  plot(c(L+10,0),xlim=c(min(-dat[,1])-10,max(dat)),type = "n",frame=F,yaxt="n",ylab="")
  mp=barplot(height = val1,add = TRUE,axes = FALSE,horiz=T,col="lightsteelblue3")
  barplot(height = val2,add = TRUE,axes = FALSE,horiz=T,col="midnightblue")
  text(x=(min(-dat[,1])-5),y=mp,lab=rownames(dat),cex=0.6)
  mtext(side=1,line=2,text="-log10[P-Value]",cex=0.6)
  mtext(side=3,line=0,text=mytitle,cex=0.6)
  
}

PlotScatterRUD <- function(selS=selShort,selL=selLong){
  #par(mfrow=c(1,1))
  mycols             <- rep("grey",nrow(mRUD))
  mycols[selS] <-rgb(154/255,162/255,197/255)
  mycols[selL] <-rgb(27/255,35/255,83/255)
  plot(mRUD[!(selS|selL),c(2,1)],pch=19,col=mycols[!(selS|selL)],cex=0.3,frame=F,las=1,xlab="",ylab="",cex.axis=0.8)
  points(mRUD[selS|selL,c(2,1)],pch=19,col=mycols[(selS|selL)],cex=0.3)
  abline(v=0,lty=2,col="grey")
  abline(h=0,lty=2,col="grey")
  mtext(side=1,line=2,text="log2 proximal-to-distal poly(A) site ratio",cex=0.8)
  mtext(side=2,line=3,text="log2 proximal-to-distal poly(A) site ratio",cex=0.8)
  mtext(side=2,line=2,text="axonal compartment",cex=0.8)
  mtext(side=1,line=3,text="cell body compartment",cex=0.8)
  subGS   <- as.character(subanno$geneSymbol)
  text(x=-8,y=15,col="lightsteelblue3",labels=paste(length(unique(subGS[selS]))," proximal shifts in axons",sep=""),cex=0.7)
  text(x=8,y=-15,col="midnightblue",labels=paste(length(unique(subGS[selL]))," distal shifts in axons",sep=""),cex=0.7)
  text(x=-8,y=14,col="black",labels=paste("n=",length(unique(subGS))," tandem 3' UTR",sep=""),cex=0.7)
  text(x=-8,y=13,col="black",labels=paste("r=",round(cor(mRUD[,1],mRUD[,2],method="spearman"),digit=2),"(spearman)",sep=""),cex=0.7)
  #  abline(a=0,b=1,lty=2,col="red")
  #  abline(v=0,lty=2,col="red")
  #  abline(h=0,lty=2,col="red")
}

# An isoforms is considered to be reliably expressed if in all the samples the probability of it belonging to the non-expressed class is below 0.05 (less than 5% chance to below to the background in both replicates) -- soft threshold
#OR
# An isoforms is considered to be expressed if in at least one of the replicate the probability of it belonging to the expressed class is above 0.1 (more than 10% chance to belong to the foreground in at least one replicate) -- hard threshold

SelectExpressed <- function(dat=log2(htseq[,1]+1),frac.bg=0.6,frac.fg=0.1){
  bimdens       <- densityMclust(data=dat,G=2)
  x   <- seq(from=0, to=max(dat),length=100)
  if(is.na(bimdens$parameters$variance$sigmasq[2])){
    Lim.fg           <- qnorm(frac.fg,mean=bimdens$parameters$mean[2],sd=sqrt(bimdens$parameters$variance$sigmasq[1]))
    Lim.bg           <- qnorm(frac.bg,mean=bimdens$parameters$mean[1],sd=sqrt(bimdens$parameters$variance$sigmasq[1]))
    hx2              <- dnorm(x,mean=bimdens$parameters$mean[2],sd=sqrt(bimdens$parameters$variance$sigmasq[1]))
  }
  if(!is.na(bimdens$parameters$variance$sigmasq[2])){
    Lim.fg           <- qnorm(frac.fg,mean=bimdens$parameters$mean[2],sd=sqrt(bimdens$parameters$variance$sigmasq[2]))
    Lim.bg           <- qnorm(frac.bg,mean=bimdens$parameters$mean[1],sd=sqrt(bimdens$parameters$variance$sigmasq[1]))
    hx2              <- dnorm(x,mean=bimdens$parameters$mean[2],sd=sqrt(bimdens$parameters$variance$sigmasq[2]))
  }
  
  hx1 <- dnorm(x,mean=bimdens$parameters$mean[1],sd=sqrt(bimdens$parameters$variance$sigmasq[1]))
  hist(dat,breaks=50,col=rgb(0,0,0,alpha=0.2),freq=FALSE,xlab="",ylab="",main="",ylim=c(0,0.4),las=1,xlim=c(0,20))
  lines(x,hx2 , lwd=2, col="blue")
  lines(x,hx1 , lwd=2, col="red")
  abline(v=Lim.fg,col="blue",lty=2)
  abline(v=Lim.bg,col="red",lty=2)
  return(c(Lim.bg,Lim.fg))
}

###


# A. Load data
load("./annotation/rn5/GRanges_comprehensive_transcriptome_rat_24_nov_2015.RData")
myUTR        <- import.gff("./annotation/rn5/Lngf.gtf",format="gtf")#To obtain the number isoforms in figure 1
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

#Supplementary Figure 2b
par(mfrow=c(1,2))
mp<-barplot(focus.iso,las=1,frame=F,col="white",ylim=c(0,30000),ylab="# novel 3' UTR isoforms")
mtext(side=3,line=0,text=focus.iso,at=mp,cex=0.5)
mp<-barplot(focus.txID,las=1,frame=F,col="white",ylim=c(0,12000),ylab="# Ensembl txID")
mtext(side=3,line=0,text=focus.txID,at=mp,cex=0.5)

maxdL.pos        <- (anno_ngf$maxL-anno_ngf$initL)[match(unique(anno_ngf$txID),anno_ngf$txID)]
maxdL.pos        <- maxdL.pos[maxdL.pos>0]
maxdL.neg        <- abs(anno_ngf$minL-anno_ngf$initL)[match(unique(anno_ngf$txID),anno_ngf$txID)]
maxdL.neg        <- maxdL.neg[maxdL.neg>0]


# C. Remove redundant transcripts and select reliably expressed genes
id       <- paste(start(myUTR),end(myUTR),strand(myUTR),sep=".")#2215 duplicated
myUTR    <- myUTR[-which(duplicated(id)),]
anno_ngf <- anno_ngf[-which(duplicated(id)),]

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
Coefficients:
  Estimate Std. Error z value Pr(>|z|)
(Intercept)          7.7101     0.0127 607.078  < 2e-16 ***
  compartmentcb.only  -0.1383     0.0184  -7.513  5.8e-14 ***

m1          <- glm.nb(length ~compartment, data = datL)
test.fitted <-  predict(m1, type = "response")
confint(m1)



# G. Alternative 3' UTR usage between the two compartments
no.genes.expressed.neurons              <- length(unique(anno_ngf$geneSymbol[NGF.axon.is.expressed.iso|NGF.cb.is.expressed.iso]))
no.genes.expressed.axons                <- length(unique(anno_ngf$geneSymbol[NGF.axon.is.expressed.iso]))

#txID.expression
axons                                   <- unique(anno_ngf$txID[anno_ngf$NGF.axon.is.expressed.iso])
cb                                      <- unique(anno_ngf$txID[anno_ngf$NGF.cb.is.expressed.iso])
neurons                                 <- union(axons,cb)
cb.only                                 <- setdiff(cb,axons)
axon.only                               <- setdiff(axons,cb)
both                                    <- intersect(axons,cb)

#uniqueID.expression
axons                                  <- anno_ngf$uniqueID[anno_ngf$NGF.axon.is.expressed.iso]
cb                                     <- anno_ngf$uniqueID[anno_ngf$NGF.cb.is.expressed.iso]
neur.neur                              <- as.character(anno_ngf$uniqueID)[anno_ngf$uniqueID%in%intersect(cb,axons)&anno_ngf$txID%in%both]
neur.cb                                <- as.character(anno_ngf$uniqueID)[anno_ngf$uniqueID%in%setdiff(cb,axons)&anno_ngf$txID%in%both]
neur.ax                                <- as.character(anno_ngf$uniqueID)[anno_ngf$uniqueID%in%setdiff(axons,cb)&anno_ngf$txID%in%both]
ax.ax                                  <- as.character(anno_ngf$uniqueID)[anno_ngf$uniqueID%in%setdiff(axons,cb)&anno_ngf$txID%in%axon.only]
cb.cb                                  <- as.character(anno_ngf$uniqueID)[anno_ngf$uniqueID%in%setdiff(cb,axons)&anno_ngf$txID%in%cb.only]

temp       <-cbind(anno_ngf$uniqueID%in%neur.neur,anno_ngf$uniqueID%in%neur.cb,anno_ngf$uniqueID%in%neur.ax,anno_ngf$uniqueID%in%cb.cb,anno_ngf$uniqueID%in%ax.ax)
iso.class <- rep("no",length(anno_ngf$uniqueID))
iso.class[anno_ngf$uniqueID%in%neur.neur] <- "neur.neur"
iso.class[anno_ngf$uniqueID%in%neur.cb]   <- "neur.cb"
iso.class[anno_ngf$uniqueID%in%neur.ax]   <- "neur.ax"
iso.class[anno_ngf$uniqueID%in%cb.cb]     <- "cb.cb"
iso.class[anno_ngf$uniqueID%in%ax.ax]     <- "ax.ax"
anno_ngf$iso.class <- iso.class


#Create Backgroun with what is expressed in neurons
mytx <- unique(as.character(anno_ngf$txID[anno_ngf$NGF.axon.is.expressed.iso|anno_ngf$NGF.cb.is.expressed.iso]))
mygs <- unique(as.character(anno_ngf$geneSymbol[anno_ngf$NGF.axon.is.expressed.iso|anno_ngf$NGF.cb.is.expressed.iso]))#12'529
myBG <- data.frame(txID=mytx,GS=c(mygs,rep(NA,length(mytx)-length(mygs))))

#Remove those transcript ID which have only one isoforms
anno_ngf     <- anno_ngf[anno_ngf$NGF.axon.is.expressed.iso|NGF.cb.is.expressed.iso,]
temp         <- as.data.frame(table(as.character(anno_ngf$txID)))
no.iso       <- temp[match(anno_ngf$txID,as.character(temp[,1])),2]
anno_ngf     <- anno_ngf[anno_ngf$no.iso>1,]
names(myUTR) <- myUTR$ID
ix           <- match(anno_ngf$uniqueID,names(myUTR))
myUTR        <- myUTR[ix[!is.na(ix)],]
anno_ngf     <- anno_ngf[!is.na(ix),]


#Check that the closest neighbour is not within 300 nt
#iso_ordered
tempL                     <- anno_ngf$newL
names(tempL)              <- anno_ngf$uniqueID
tempG                     <- as.factor(as.character(anno_ngf$txID))
myord                     <- tapply(tempL,INDEX=tempG,function(x)return(cbind(names(x)[sort(as.numeric(as.character(x)),index.return=T,decreasing=F)$ix],c(1:length(x)))))
test                      <- do.call(what=rbind,args=myord)
iso_ordered               <- as.numeric(test[match(anno_ngf$uniqueID,test[,1]),2])
#nextIsoform
tempid                    <- paste(anno_ngf$txID,iso_ordered,sep=".")
next.iso                  <- paste(anno_ngf$txID,(as.numeric(as.character(iso_ordered))+1),sep=".")
nextID                    <- anno_ngf$uniqueID[match(next.iso,tempid)]
#distTonext
distonext                 <- rep(NA,nrow(anno_ngf))
sel                       <- !is.na(nextID)
d1                        <- as.numeric(as.character(anno_ngf$newL))[sel]
d2                        <- as.numeric(as.character(anno_ngf$newL))[match(nextID[sel],anno_ngf$uniqueID)]
distonext[sel]            <- d2-d1
#Fproxi
Fproxi                    <- which(distonext<300)#should be none

# Get ordering of the isoforms
tempL                     <- anno_ngf$newL
names(tempL)              <- anno_ngf$uniqueID
tempG                     <- as.factor(as.character(anno_ngf$txID))
myord                     <- tapply(tempL,INDEX=tempG,function(x)return(cbind(names(x)[sort(as.numeric(as.character(x)),index.return=T,decreasing=F)$ix],c(1:length(x)))))
test                      <- do.call(what=rbind,args=myord)
test                      <- test[match(anno_ngf$uniqueID,test[,1]),]
temp                      <- as.data.frame(table(as.character(anno_ngf$txID)))
no.iso                    <- temp[match(anno_ngf$txID,as.character(temp[,1])),2]
anno_ngf$no.iso           <- no.iso
anno_ngf$iso_ordered      <- test[,2]
tempL                     <- as.character(anno_ngf$iso_ordered)
names(tempL)              <- as.character(anno_ngf$uniqueID)
imID1                     <- tapply(tempL,INDEX=as.factor(as.character(anno_ngf$txID)),function(x)return(names(x)[x=="1"]))
imIDp                     <- as.character(imID1[match(anno_ngf$txID,names(imID1))])
anno_ngf$imIDp            <- imIDp


# Get the selection on which to focus for the analysis
# At least one of the pair (distal or proximal) must be detected in axons
sela1                <- (anno_ngf$NGF.axon.is.expressed.iso|anno_ngf$NGF.axon.is.expressed.iso[match(anno_ngf$imIDp,anno_ngf$uniqueID)])
# At least one of the pair (distal or proximal) must be detected in cb
sela2                <- (anno_ngf$NGF.cb.is.expressed.iso|anno_ngf$NGF.cb.is.expressed.iso[match(anno_ngf$imIDp,anno_ngf$uniqueID)])
# Remove all the isoform I0
sela3                 <- anno_ngf$uniqueID!=anno_ngf$imIDp
#Proximal is expressed in at least one of the 2 samples
sela4                 <- (anno_ngf$NGF.axon.is.expressed.iso|anno_ngf$NGF.cb.is.expressed.iso)[match(anno_ngf$imIDp,anno_ngf$uniqueID)]
#Distal is expressed in at least one of the 2 samples
sela5                 <- (anno_ngf$NGF.axon.is.expressed.iso|anno_ngf$NGF.cb.is.expressed.iso)

sela                  <- sela1&sela2&sela3&sela4&sela5
subanno               <- anno_ngf[sela,]

#
# C. IDENTIFICATION OF DIFFERENTIALLY EXPRESSED ISOFORMS
#

# C.1    Proximal to distal site usage AND fisher test
ix1                  <- c(1:nrow(anno_ngf))
ix2                  <- match(as.character(anno_ngf$imIDp),as.character(anno_ngf$uniqueID))
mydat                <- anno_ngf[,match(c("NGF.axon.1.raw","NGF.axon.2.raw","NGF.cb.1.raw","NGF.cb.2.raw"),colnames(anno_ngf))]
for(i in c(1:ncol(mydat))){mydat[,i]<-log2(1+mydat[,i])}
rownames(mydat)      <- anno_ngf$uniqueID
myProximal           <- apply(mydat,2,function(x)return(x[ix2]))
myDistal             <- apply(mydat,2,function(x)return(x[ix1]))
rownames(myProximal) <-rownames(myDistal)<-anno_ngf$uniqueID[ix1]
RUD                  <- do.call(lapply(c(1:4),function(x)return(myProximal[,x]-myDistal[,x])),what=cbind)
rownames(RUD)        <- anno_ngf$uniqueID[ix1]
colnames(RUD)        <- c("axon.1","axon.2","cb.1","cb.2")
RUD                  <- RUD[sela,]#8'853
sdRUD                 <- t(apply(RUD,1,function(x)return(tapply(x,INDEX=factor(c("axon","axon","cb","cb")),FUN=sd))))
mRUD                  <- t(apply(RUD,1,function(x)return(tapply(x,INDEX=factor(c("axon","axon","cb","cb")),FUN=mean))))
dRUD                  <- mRUD[,2]-mRUD[,1]

#Compare replicates
par(mfrow=c(2,2))
smoothScatter(RUD[,1],RUD[,2],frame=FALSE,xlab="ax1",ylab="ax2")
mtext(side=3,line=0,text=cor(RUD[,1],RUD[,2],method="spearman"))
smoothScatter(RUD[,3],RUD[,4],frame=FALSE,xlab="cb.1",ylab="cb.2")
mtext(side=3,line=0,text=cor(RUD[,3],RUD[,4],method="spearman"))

sumRUD               <- cbind(
  apply(anno_ngf[,match(c("NGF.axon.1.raw","NGF.axon.2.raw","NT3.axon.1.raw","NT3.axon.2.raw"),colnames(anno_ngf))],1,sum),
  apply(anno_ngf[,match(c("NGF.cb.1.raw","NGF.cb.2.raw","NT3.cb.1.raw","NT3.cb.2.raw"),colnames(anno_ngf))],1,sum)
)

#Fisher test on sum of RAW count data
sumRUD               <- cbind(
  apply(anno_ngf[,match(c("NGF.axon.1.raw","NGF.axon.2.raw"),colnames(anno_ngf))],1,sum),
  apply(anno_ngf[,match(c("NGF.cb.1.raw","NGF.cb.2.raw"),colnames(anno_ngf))],1,sum)
)
rownames(sumRUD)<- anno_ngf$uniqueID
colnames(sumRUD)<- c("axon","cb")
test.apa <- function(ix.proximal=ix1[1],ix.distal=ix2[1]){
  return(fisher.test(round(cbind(sumRUD[ix.proximal,],sumRUD[ix.distal,])))$p.value)
}
fisherRUD   <- unlist(lapply(c(1:nrow(sumRUD)),function(x)return(test.apa(ix.proximal=ix1[x],ix.distal=ix2[x]))))


my.proximal <- sumRUD[ix2,]
my.distal   <- sumRUD[ix1,]
rel.proximal.usage <- cbind(my.proximal[,1]/(my.distal[,1]+my.proximal[,1]),
                            my.proximal[,2]/(my.distal[,2]+my.proximal[,2]))

colnames(rel.proximal.usage) <- c("axon","cb")
rownames(rel.proximal.usage)  <- rownames(sumRUD)

sumRUD                 <- sumRUD[sela,]
fisherRUD              <- fisherRUD[sela]
rel.proximal.usage     <- rel.proximal.usage[sela,]
diff.rel.proximal.usage<- rel.proximal.usage[,2]-rel.proximal.usage[,1]
# Correct for multitest (remove non-important ones)
padjRUD              <-  p.adjust(fisherRUD,method="fdr")
mRUD                 <-  as.data.frame(mRUD)

# Selection
id.proxi             <-  anno_ngf$imIDp[match(names(dRUD),anno_ngf$uniqueID)]
seld1                <-  anno_ngf$NGF.axon.is.expressed.iso[match(id.proxi,anno_ngf$uniqueID)]
seld2                <-  anno_ngf$NGF.axon.is.expressed.iso[match(names(dRUD),anno_ngf$uniqueID)]
# Significant P-value of difference
selc                 <- padjRUD<0.01


#Relative Proximal to distal site usage
mydat                <- anno_ngf[,match(c("NGF.axon.1.raw","NGF.axon.2.raw","NGF.cb.1.raw","NGF.cb.2.raw"),colnames(anno_ngf))]
ix.distal            <- c(1:nrow(anno_ngf))
ix.proximal          <- match(as.character(anno_ngf$imIDp),as.character(anno_ngf$uniqueID))
psi                  <- mydat
for(i in c(1:ncol(psi))){
  psi[,i]                                                              <- mydat[ix.proximal,i]/(mydat[ix.distal,i]+mydat[ix.proximal,i])
  psi[is.na(psi[,i]),i]<-0
}
PUD                 <- psi
rownames(PUD)       <- anno_ngf$uniqueID
colnames(PUD)       <- c("axon.1","axon.2","cb.1","cb.2")
mPUD                <- t(apply(PUD,1,function(x)return(tapply(x,INDEX=factor(c("axon","axon","cb","cb")),FUN=mean))))
PUD                 <- PUD[sela,]#8'853
mPUD                <- mPUD[sela,]
diffPUD             <-  mPUD[,"axon"]-mPUD[,"cb"]
mPUD <- data.frame(mPUD)


Lngf       <- import.gff("~/Desktop/DataAnalsyisRiccio/Dec2016/utrid/APA_stringent/Lngf.gtf",format="gtf")
Lngf.sub.1 <- Lngf[setdiff(which(Lngf$is.pas=="FALSE"),grep(Lngf$ID,pattern="\\.0"))]
Lngf.sub.2 <- Lngf[intersect(which(Lngf$is.pas=="FALSE"),grep(Lngf$ID,pattern="\\.0")),]


anno_ngf$NGF.axon.is.expressed.txID <- anno_ngf$txID%in%unique(as.character(anno_ngf$txID))[anno_ngf$NGF.axon.is.expressed.iso]
anno_ngf$NGF.cb.is.expressed.txID   <- anno_ngf$txID%in%unique(as.character(anno_ngf$txID))[anno_ngf$NGF.cb.is.expressed.iso]

temp1 <- anno_ngf$utrL[anno_ngf$NGF.axon.is.expressed.iso]
temp2 <- anno_ngf$utrL[anno_ngf$NGF.cb.is.expressed.iso|!anno_ngf$NGF.axon.is.expressed.iso]
boxplot(list(maxL2,maxL1),outline=F,col=c("yellow","grey"))

#When alternative, what is the longest used; distal shift versus proximal shift

temp1 <- anno_ngf[anno_ngf$NGF.axon.is.expressed.iso,]
temp2 <- anno_ngf[anno_ngf$NGF.cb.is.expressed.iso&!anno_ngf$NGF.axon.is.expressed.iso,]
maxL1 <- tapply(temp1$newL,INDEX=factor(as.character(temp1$txID)),FUN=max)#4'192
maxL2 <- tapply(temp2$newL,INDEX=factor(as.character(temp2$txID)),FUN=max)#6'478
boxplot(list(maxL2,maxL1),outline=F,col=c("yellow","grey"))


temp1 <- anno_ngf[anno_ngf$NGF.axon.is.expressed.txID,]
temp2 <- anno_ngf[anno_ngf$NGF.cb.is.expressed.txID&!anno_ngf$NGF.axon.is.expressed.txID,]
maxL1 <- tapply(temp1$newL,INDEX=factor(as.character(temp1$txID)),FUN=max)#4'192
maxL2 <- tapply(temp2$newL,INDEX=factor(as.character(temp2$txID)),FUN=max)#6'478
boxplot(list(maxL2,maxL1),outline=F,col=c("yellow","grey"))

# C.3 Selection of genes of interest
opp.1 <- mRUD$axon>0&mRUD$cb<0&mPUD$axon>0.5&mPUD$cb<0.5
opp.2 <- mRUD$axon<0&mRUD$cb>0&mPUD$axon<0.5&mPUD$cb>0.5
Lim1  <- 0.2
Lim2  <- 0.8

selRUD1             <- list()
IX=1
selRUD1[[IX]]       <- padjRUD<0.01&dRUD<(-1)&seld1#980; proximal shift so proximal is expressed in axons
IX=IX+1
selRUD1[[IX]]       <- padjRUD<0.01&dRUD>(1)&seld2#1065; distal shift so distal expressed in axons

IX=IX+1
selRUD1[[IX]]      <- padjRUD<0.01&dRUD<(-1)&seld1&diffPUD>(0.15)#689-proximal shift so proximal is expressed in axons
IX=IX+1
selRUD1[[IX]]      <- padjRUD<0.01&dRUD>(1)&seld2&diffPUD<(-0.15)#737-distal shift so distal expressed in axons

IX=IX+1
selRUD1[[IX]]      <- padjRUD<0.01&dRUD<(-1)&seld1&diff.rel.proximal.usage<(-0.15)&rel.proximal.usage[,2]<=Lim1#178, less that 20% usage of the proximal in CB; therefore can be considered as absent
IX=IX+1
selRUD1[[IX]]      <- padjRUD<0.01&dRUD>(1)&seld2&diff.rel.proximal.usage>(0.15)&rel.proximal.usage[,2]>=Lim2#86,less than 20% usage of the distal in CB; therefore can be considered as absent

View(cbind(subanno[selRUD1[[5]],c("uniqueID","geneSymbol","imIDp","newL")],mPUD[selRUD1[[5]],]))
View(cbind(subanno[selRUD1[[6]],c("uniqueID","geneSymbol","imIDp","newL")],mPUD[selRUD1[[6]],]))

write.csv(cbind(subanno[selRUD1[[5]],c("uniqueID","geneSymbol","imIDp","newL")],mPUD[selRUD1[[5]],]),"~/Desktop/DataAnalsyisRiccio/Dec2016/axonal_remodelling/proximal_shifts.csv")
write.csv(cbind(subanno[selRUD1[[6]],c("uniqueID","geneSymbol","imIDp","newL")],mPUD[selRUD1[[6]],]),"~/Desktop/DataAnalsyisRiccio/Dec2016/axonal_remodelling/distal_shifts.csv")


my.no <- unlist(lapply(selRUD1,function(Z)return(length(unique(subanno$txID[Z])))))

pdf("~/Desktop/DataAnalysisRiccio/Dec2016/axonal_remodelling/alternative/no.uniqs_80.pdf")
mp<- barplot(my.no[c(5,6,9,10)],las=1)
mtext(side=3,line=0,text=my.no[c(5,6,9,10)],at=mp)
dev.off()



mylist <- lapply(selRUD1,function(Z)return(unique(subanno$txID[Z])))

setdiff(cand$txID,mylist[[5]])
intersect(cand$txID,mylist[[5]])

View(cbind(as.character(subanno$GS),diff.rel.proximal.usage,rel.proximal.usage)[subanno$txID%in%setdiff(cand$txID,mylist[[9]]),])
cand <- read.csv("~/Desktop/DataAnalysisRiccio/Riccio_exp1_analysis/candidates_rem.csv")



ix <- which(subanno$txID=="ENSRNOT00000016995")
ix <- which(subanno$txID=="ENSRNOT00000065177")
ix <- which(subanno$GS=="Zmynd11")
rel.proximal.usage[ix,]
(padjRUD<0.01&dRUD<(-1)&seld1&diffPUD>(0.15))[ix]
rel.proximal.usage[ix,2]

ix <- which(anno_ngf$GS=="Zmynd11")
my.proximal[ix,]
my.distal[ix,]
(cbind(my.proximal[,1]/(my.distal[,1]+my.proximal[,1]),
my.proximal[,2]/(my.distal[,2]+my.proximal[,2])))[ix,]


View(anno_ngf[,match(c("uniqueID","txID","GS","NGF.cb.1.raw","NGF.cb.2.raw","NGF.axon.1.raw","NGF.axon.2.raw"),colnames(anno_ngf))])

View(anno_ngf[which(anno_ngf$GS=="Zmynd11"),match(c("uniqueID","txID","GS","NGF.cb.1.raw","NGF.cb.2.raw","NGF.axon.1.raw","NGF.axon.2.raw"),colnames(anno_ngf))])


View(subanno[selRUD1[[9]],])
View(subanno[selRUD1[[10]],])

rel.proximal.usage
diff.rel.proximal.usage

table(anno_ngf$iso.class[match(unique(subanno$imIDp[selRUD1[[1]]]),anno_ngf$uniqueID)])
table(subanno$iso.class[selRUD1[[2]]])

View(anno_ngf[which(anno_ngf$uniqueID%in%subanno$imIDp[padjRUD<0.01&dRUD>(1)&seld2&diffPUD<(-0.15)&rel.proximal.usage[,2]<0.2]),match(c("uniqueID","txID","GS","NGF.cb.1.raw","NGF.cb.2.raw","NGF.axon.1.raw","NGF.axon.2.raw"),colnames(anno_ngf))])

View(anno_ngf[which(anno_ngf$uniqueID%in%subanno$imIDp[selRUD1[[1]]]&anno_ngf$iso.class=="neur.ax"),match(c("uniqueID","txID","GS","NGF.cb.1.raw","NGF.cb.2.raw","NGF.axon.1.raw","NGF.axon.2.raw"),colnames(anno_ngf))])
View(subanno[which(selRUD1[[2]]&subanno$iso.class=="neur.ax"),match(c("uniqueID","txID","GS","NGF.cb.1.raw","NGF.cb.2.raw","NGF.axon.1.raw","NGF.axon.2.raw"),colnames(subanno))])

s1p <- anno_ngf[which(anno_ngf$uniqueID%in%subanno$imIDp[selRUD1[[5]]]),match(c("uniqueID","txID","GS","NGF.cb.1.raw","NGF.cb.2.raw","NGF.axon.1.raw","NGF.axon.2.raw"),colnames(anno_ngf))]
s2p <- anno_ngf[which(anno_ngf$uniqueID%in%subanno$imIDp[selRUD1[[7]]]),match(c("uniqueID","txID","GS","NGF.cb.1.raw","NGF.cb.2.raw","NGF.axon.1.raw","NGF.axon.2.raw"),colnames(anno_ngf))]

s1d  <- subanno[which(selRUD1[[6]]),match(c("uniqueID","txID","GS","NGF.cb.1.raw","NGF.cb.2.raw","NGF.axon.1.raw","NGF.axon.2.raw"),colnames(subanno))]
s2d  <- subanno[which(selRUD1[[8]]),match(c("uniqueID","txID","GS","NGF.cb.1.raw","NGF.cb.2.raw","NGF.axon.1.raw","NGF.axon.2.raw"),colnames(subanno))]

write.csv(s2p,paste(outdir,"candidates_axonal_remodelling.csv",sep=""))
write.csv(s2d,paste(outdir,"candidates_faciliated_transport.csv",sep=""))



sampleGO.RUD1      <- lapply(selRUD1[c(3,4)],CreateSampleGO)
enrich.RUD1        <- lapply(sampleGO.RUD1,FUN=function(X)return(lapply(X,getEnrich)))
F1.w0     <- CompareBP_improved(enr1=enrich.RUD1[[1]][[2]][[2]],enr2=enrich.RUD1[[2]][[2]][[2]],PLOT=TRUE,no=20,lab2="long",lab1="short",coi="weight0Fisher")
F1.fisher <- CompareBP_improved(enr1=enrich.RUD1[[1]][[2]][[1]],enr2=enrich.RUD1[[2]][[2]][[1]],PLOT=TRUE,no=20,lab2="long",lab1="short",coi="classicFisher")
write.csv(file=paste(outdir,"F1.w0.r.csv",sep=""),F1.w0)
write.csv(file=paste(outdir,"F1.fisher.r.csv",sep=""),F1.fisher)



pdf(paste(outdir,"characterisation_APA_stringent.pdf",sep=""))
layout(matrix(c(1,1,1,1), 2, 2, byrow = FALSE))
par(mar=c(3,3,3,3))
final_enrich <- read.csv(paste(outdir,"F1.w0.r_f.csv",sep=""))
PlotScatterRUD(selS=selRUD1[[3]],selL=selRUD1[[4]])
#par(mfrow=c(2,2),mar=c(3,10,3,3))
vec1 <- final_enrich[,2]
vec2 <- final_enrich[,3]
names(vec1)<-names(vec2)<-as.character(final_enrich[,1])
vec1  <- sort(vec1,decreasing=F)
vec2  <- sort(vec2,decreasing=F)
vec1  <- vec1[vec1>=-log10(0.05)]
vec2  <- vec2[vec2>=-log10(0.05)]
barplot(vec1,horiz=TRUE,col=rgb(154/255,162/255,197/255),las=1)
mtext(side=3,line=0,text=foi[IX])
mtext(side=1,line=2,text="-log10(P-Value)")
barplot(vec2,horiz=TRUE,col=rgb(27/255,35/255,83/255),las=1)
mtext(side=1,line=2,text="-log10(P-Value)")

dat <- c(my.nos,length(unique(subanno$geneSymbol[selRUD1[[3]]])),length(unique(subanno$geneSymbol[selRUD1[[4]]])))
names(dat)[c(4,5)] <- c("proximal","distal")
mp  <- barplot(dat,las=1,ylab="no.GS")
mtext(side=3,line=0,text=dat,at=mp)
dev.off()

#Study how many express only the distal compared to all in each pairs: additional analysis for Antonella
ix.distal            <- c(1:nrow(anno_ngf))
ix.proximal          <- match(as.character(anno_ngf$imIDp),as.character(anno_ngf$uniqueID))
myIDs                <- cbind(anno_ngf$imIDp,as.character(anno_ngf$uniqueID))[sela,]#6'949 rows
cor.txID             <- anno_ngf$txID[sela]#4'191
cor.GS               <- anno_ngf$GS[sela]#4'102

#Analysis in axons
expr.axons           <- anno_ngf$uniqueID[anno_ngf$NGF.axon.is.expressed.iso]
expr.id.axons        <- apply(myIDs,2,function(Z)return(Z%in%expr.axons))
torm1                <- apply(expr.id.axons,1,sum)==0
colnames(expr.id.axons)<- c("proximal","distal")

only.proxi.axons     <- expr.id.axons[,1]&!expr.id.axons[,2]
only.prox.GS         <- tapply(only.proxi.axons,INDEX=factor(cor.GS),FUN=function(Z)sum(!Z)==0)#691
only.prox.GS         <- names(only.prox.GS)[only.prox.GS]
only.prox.txID       <- tapply(only.proxi.axons,INDEX=factor(cor.txID),FUN=function(Z)sum(!Z)==0)#708
only.prox.txID       <- names(only.prox.txID)[only.prox.txID]

only.dist.axons      <- !expr.id.axons[,1]&expr.id.axons[,2]
only.dist.GS         <- unique(cor.GS[only.dist.axons])#1'872
only.dist.txID       <- unique(cor.txID[only.dist.axons])#1'900

both.axons           <- expr.id.axons[,1]&expr.id.axons[,2]
both.GS              <- unique(cor.GS[both.axons])#1'545
both.txID            <- unique(cor.txID[both.axons])#1'583

#Please note that some GS are dual due to multiple transcript ID, and therefore the number of txID is more reliable
axons.compare.GS     <- list(both=both.GS,distal.only=only.dist.GS,short.only=only.prox.GS)
axons.compare.txID   <- list(both=both.txID,distal.only=only.dist.txID,short.only=only.prox.txID)


# Analysis in CB
expr.cb             <- anno_ngf$uniqueID[anno_ngf$NGF.cb.is.expressed.iso]
expr.id.cb          <- apply(myIDs,2,function(Z)return(Z%in%expr.cb))
torm1               <- apply(expr.id.cb,1,sum)==0
colnames(expr.id.cb)<- c("proximal","distal")

only.proxi.cb        <- expr.id.cb[,1]&!expr.id.cb[,2]
only.prox.GS         <- tapply(only.proxi.cb,INDEX=factor(cor.GS),FUN=function(Z)sum(!Z)==0)#3
only.prox.GS         <- names(only.prox.GS)[only.prox.GS]
only.prox.txID       <- tapply(only.proxi.cb,INDEX=factor(cor.txID),FUN=function(Z)sum(!Z)==0)#3
only.prox.txID       <- names(only.prox.txID)[only.prox.txID]

only.dist.cb      <- !expr.id.cb[,1]&expr.id.cb[,2]
only.dist.GS         <- unique(cor.GS[only.dist.cb])#5
only.dist.txID       <- unique(cor.txID[only.dist.cb])#5

both.cb           <- expr.id.cb[,1]&expr.id.cb[,2]
both.GS              <- unique(cor.GS[both.cb])#4095
both.txID            <- unique(cor.txID[both.cb])#4183

#Please note that some GS are dual due to multiple transcript ID, and therefore the number of txID is more reliable
cb.compare.GS     <- list(both=both.GS,distal.only=only.dist.GS,short.only=only.prox.GS)
cb.compare.txID   <- list(both=both.txID,distal.only=only.dist.txID,short.only=only.prox.txID)

mydat                <- anno_ngf[,match(c("NGF.axon.1.raw","NGF.axon.2.raw","NGF.cb.1.raw","NGF.cb.2.raw"),colnames(anno_ngf))]
for(i in c(1:ncol(mydat))){mydat[,i]<-log2(1+mydat[,i])}
rownames(mydat)      <- anno_ngf$uniqueID
myMean <- t(apply(mydat,1,function(Z)return(tapply(Z,INDEX=factor(c("axon","axon","cb","cb")),FUN=mean))))


myshifts.txID.1            <- list(proxi=unique(subanno$txID[selRUD1[[3]]]),distal=unique(subanno$txID[selRUD1[[4]]]))
myshifts.txID.2            <- list(proxi=unique(subanno$txID[selRUD1[[7]]]),distal=unique(subanno$txID[selRUD1[[8]]]))

myfrac.1 <- c(distal.shift=100*length(intersect(axons.compare.txID[[2]],myshifts.txID.1[[2]]))/length(myshifts.txID.1[[2]]),
              proximal.shift=100*length(intersect(axons.compare.txID[[3]],myshifts.txID.1[[1]]))/length(myshifts.txID.1[[1]]))
myfrac.2 <- c(distal.shift=100*length(intersect(axons.compare.txID[[2]],myshifts.txID.2[[2]]))/length(myshifts.txID.2[[2]]),
              proximal.shift=100*length(intersect(axons.compare.txID[[3]],myshifts.txID.2[[1]]))/length(myshifts.txID.2[[1]]))

pdf("~/Desktop/temp2.pdf")
par(mfrow=c(2,2))
barplot(myfrac.1,col=c(rgb(30/255,36/255,83/255),rgb(154/255,162/255,197/255)),frame=F,las=1,ylim=c(0,100))
barplot(myfrac.,col=c(rgb(30/255,36/255,83/255),rgb(154/255,162/255,197/255)),frame=F,las=1,ylim=c(0,100))


pdf("~/Desktop/DataAnalysisRiccio/Dec2016/axonal_remodelling/alternative/only.short.analysis.pdf")

#Eith compare with proximal shift in axons = selRUD[[3]] and distal shift in axons =selRUD[[4]]; however maybe some are simply due to the diffusion
myshifts.txID            <- list(proxi=unique(subanno$txID[selRUD1[[3]]]),distal=unique(subanno$txID[selRUD1[[4]]]))
myshifts.GS              <- list(proxi=unique(subanno$GS[selRUD1[[3]]]),distal=unique(subanno$GS[selRUD1[[4]]]))
proxi.focus              <- list( proxi.and.both=intersect(myshifts.txID[[1]],axons.compare.txID[[1]]),
                                  proxi.and.distal=intersect(myshifts.txID[[1]],axons.compare.txID[[2]]),
                                  proxi.and.proxi=intersect(myshifts.txID[[1]],axons.compare.txID[[3]]))

distal.focus             <- list( dist.and.both=intersect(myshifts.txID[[2]],axons.compare.txID[[1]]),
                                  dist.and.distal=intersect(myshifts.txID[[2]],axons.compare.txID[[2]]),
                                  dist.and.proxi=intersect(myshifts.txID[[2]],axons.compare.txID[[3]]))

lapply(proxi.focus,length)#54 of which proximal shift in axons AND only express short AND distal more expressed in CB
lapply(distal.focus,length)#75 of which distal shift in axons AND only express long AND proximal more in CB



#OR compare with those which have dramatic shifts
selt1      <- selRUD1[[7]]#185-proximal shift so proximal is expressed in axons
selt2      <- selRUD1[[8]]#134-distal shift so distal expressed in axons
myshifts.txID            <- list(proxi=unique(subanno$txID[selt1]),distal=unique(subanno$txID[selt2]))
myshifts.GS              <- list(proxi=unique(subanno$GS[selt1]),distal=unique(subanno$GS[selt2]))

proxi.focus              <- list( proxi.and.both=intersect(myshifts.txID[[1]],axons.compare.txID[[1]]),
                                  proxi.and.distal=intersect(myshifts.txID[[1]],axons.compare.txID[[2]]),
                                  proxi.and.proxi=intersect(myshifts.txID[[1]],axons.compare.txID[[3]]))

distal.focus             <- list( dist.and.both=intersect(myshifts.txID[[2]],axons.compare.txID[[1]]),
                                  dist.and.distal=intersect(myshifts.txID[[2]],axons.compare.txID[[2]]),
                                  dist.and.proxi=intersect(myshifts.txID[[2]],axons.compare.txID[[3]]))

lapply(proxi.focus,length)#54 of which proximal shift in axons AND only express short AND distal more expressed in CB
lapply(distal.focus,length)#75 of which distal shift in axons AND only express long AND proximal more in CB

par(mfrow=c(2,2),mar=c(3,4,3,3))
#layout(matrix(c(1,4,4,2,4,4,3,4,4),4,3,byrow = TRUE))
mp=barplot(c(distalshift=length(myshifts.txID[[2]]),proximashift=length(myshifts.txID[[1]])),col=c(rgb(224/255,40/255,38/255),rgb(194/255,164/255,152/255)),las=1,ylab="# Ensembl txID",ylim=c(0,140))
mtext(side=3,line=0,at=mp,text=c(length(myshifts.txID[[2]]),length(myshifts.txID[[1]])))
mp=barplot(c(distal.and.distalshift=lapply(distal.focus,length)[[2]],prox.only.and.proxi.shift=lapply(proxi.focus,length)[[3]]),col=c(rgb(30/255,36/255,83/255),rgb(154/255,162/255,197/255)),las=1,ylab="# Ensembl txID",ylim=c(0,140))
mtext(side=3,line=0,at=mp,text=c(lapply(distal.focus,length)[[2]],prox.only.and.proxi.shift=lapply(proxi.focus,length)[[3]]))

#Check RUD in CB for those which express either only long or only short with significant distal shift
comp.RUD         <- list(all.cb=mRUD[,2],all.ax=mRUD[,1],
                         proxi.and.proxi.cb=mRUD[which(subanno$txID%in%proxi.focus[[3]]),2],proxi.and.proxi.ax=mRUD[which(subanno$txID%in%proxi.focus[[3]]),1],
                         dist.and.dist.cb=mRUD[which(subanno$txID%in%distal.focus[[2]]),2],dist.and.dist.ax=mRUD[which(subanno$txID%in%distal.focus[[2]]),1]
)

my.no.txID       <- list(unique(subanno$txID),unique(proxi.focus[[3]]),unique(distal.focus[[2]]))
comp.utrL        <- lapply(my.no.txID,function(Z)return(anno_ngf$maxL[match(Z,anno_ngf$txID)]))

boxplot(comp.RUD,outline=F,las=1,frame=FALSE,col=c(rgb(22/255,70/255,120/255),rgb(112/255,191/255,68/255)),ylab="log2 proximal-to-distal poly(A) site ratio")
mtext(side=3,at=c(1.5,3.5,5.5),text=unlist(lapply(my.no.txID,length)))
legend("topright",col=c(rgb(22/255,70/255,120/255),rgb(112/255,191/255,68/255)),pch=15,bty="n",leg=c("cb","axons"))
boxplot(comp.utrL,outline=F,las=1,frame=FALSE,col=c("grey",rgb(154/255,162/255,197/255),rgb(30/255,36/255,83/255)),ylab="3' UTR length")

selS <- selRUD1[[3]]
selL <- selRUD1[[4]]
mycols             <- rep("grey",nrow(mRUD))
mycols[selS]       <- rgb(154/255,162/255,197/255)
mycols[selL]       <- rgb(30/255,36/255,83/255)
mycols[selS&subanno$txID%in%proxi.focus[[3]]&selt1] <- "red"
mycols[selL&subanno$txID%in%distal.focus[[2]]&selt2]<- "red"


plot(mRUD[!(selS|selL),c(2,1)],pch=19,col=mycols[!(selS|selL)],cex=0.3,frame=F,las=1,xlab="",ylab="",cex.axis=0.8)
points(mRUD[selS|selL,c(2,1)],pch=19,col=mycols[(selS|selL)],cex=0.3)
abline(v=0,lty=2,col="grey")
mtext(side=1,line=2,text="log2 proximal-to-distal poly(A) site ratio",cex=0.8)
abline(h=0,lty=2,col="grey")
mtext(side=2,line=2,text="axonal compartment",cex=0.8)
mtext(side=2,line=3,text="log2 proximal-to-distal poly(A) site ratio",cex=0.8)
mtext(side=1,line=3,text="cell body compartment",cex=0.8)
subGS   <- as.character(subanno$txID)
text(x=-8,y=15,col="lightsteelblue3",labels=paste(length(unique(subGS[selS]))," proximal shifts in axons",sep=""),cex=0.7)
text(x=8,y=-15,col="midnightblue",labels=paste(length(unique(subGS[selL]))," distal shifts in axons",sep=""),cex=0.7)
text(x=-8,y=14,col="black",labels=paste("n=",length(unique(subGS))," tandem 3' UTR",sep=""),cex=0.7)
text(x=-8,y=13,col="black",labels=paste("r=",round(cor(mRUD[,1],mRUD[,2],method="spearman"),digit=2),"(spearman)",sep=""),cex=0.7)




dev.off()











#Perform Gene Enrichment Analysis
sampleGO.txOI      <- lapply(tx.oi,CreateSampleGOList)
enrich.txOI        <- lapply(sampleGO.txOI,FUN=function(X)return(lapply(X,getEnrich)))
write.csv(enrich.txOI[[1]][[2]][[2]],file=paste(outdir,"only.short.enrich.csv",sep=""))
write.csv(enrich.txOI[[3]][[2]][[2]],file=paste(outdir,"only.long.enrich.csv",sep=""))

enr1 <- read.csv(paste(outdir,"only.short.enrich_f.csv",sep=""))
enr2 <- read.csv(paste(outdir,"only.long.enrich_f.csv",sep=""))



pdf(paste(outdir,"only.long.short.pdf",sep=""))
par(mfrow=c(2,2))
layout(matrix(c(1,1,2,3), 2, 2, byrow = FALSE))
mp<-barplot(unlist(lapply(tx.oi,length))[c(1,3,5,2,4,6)],col=c("lightsteelblue3","midnightblue","grey"),las=1)
mtext(side=3,at=mp,text=unlist(lapply(tx.oi,length))[c(1,3,5,2,4,6)],cex=0.4)

mp<-barplot(unlist(lapply(gs.oi,length))[c(1,3,5,2,4,6)],col=c("lightsteelblue3","midnightblue","grey"),las=1)
mtext(side=3,at=mp,text=unlist(lapply(gs.oi,length))[c(1,3,5,2,4,6)],cex=0.4)


dat<- -log10(enr1$weight0Fisher)
names(dat)<- as.character(enr1$Term)
dat       <- sort(dat)
barplot(dat,col="lightsteelblue3",horiz=TRUE,las=1)
mtext(side=1,line=3,text="-log10(P-Value)")
mtext(side=3,line=0,text="Only short in axons")

dat       <- -log10(enr2$weight0Fisher)
names(dat)<- as.character(enr2$Term)
dat       <- sort(dat)
barplot(dat,col="midnightblue",horiz=TRUE,las=1)
mtext(side=1,line=3,text="-log10(P-Value)")
mtext(side=3,line=0,text="Only long in axons")
dev.off()








