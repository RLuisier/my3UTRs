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




#Working on MacBook -- March
load("~/Desktop/DataAnalysisRiccio/GRanges_comprehensive_transcriptome_v78_rn5_23022014.RData")
load("~/Desktop/DataAnalysisRiccio/Dec2016/utrid/anno_ngf_stringent/anno_ngf_March_12.RData")#Contain "anno_ngf","ngfGRS"
myUTR        <- import.gff("~/Desktop/DataAnalysisRiccio/Dec2016/utrid/APA_stringent/Lngf_sub.gtf",format="gtf")
names(myUTR) <- as.character(myUTR$ID)
ix1          <- match(names(myUTR),anno_ngf$uniqueID)
myUTR        <- myUTR[!is.na(ix1),]
ix1          <- match(names(myUTR),anno_ngf$uniqueID)
anno_ngf    <- anno_ngf[ix1,]

#Working on iMac/new iMac
load("~/Desktop/DataAnalsyisRiccio/GRanges_comprehensive_transcriptome_v78_rn5_23022014.RData")
load("~/Desktop/DataAnalsyisRiccio/Dec2016/utrid/anno_ngf_stringent/anno_ngf_March_12.RData")#Contain "anno_ngf","ngfGRS"
myUTR        <- import.gff("~/Desktop/DataAnalsyisRiccio/Dec2016/utrid/APA_stringent/Lngf.gtf",format="gtf")#To obtain the number isoforms in figure 1
myUTR        <- import.gff("~/Desktop/DataAnalsyisRiccio/Dec2016/utrid/APA_stringent/Lngf_sub.gtf",format="gtf")
names(myUTR) <- as.character(myUTR$ID)
ix1          <- match(names(myUTR),anno_ngf$uniqueID)
myUTR        <- myUTR[!is.na(ix1),]
ix1          <- match(names(myUTR),anno_ngf$uniqueID)
anno_ngf    <- anno_ngf[ix1,]
outdir       <-  "~/Desktop/DataAnalsyisRiccio/Dec2016/axonal_remodelling/alternative/"



myUTR        <- import.gff("~/Desktop/DataAnalsyisRiccio/Dec2016/utrid/APA_stringent/Lngf.gtf",format="gtf")#To obtain the number isoforms in figure 1

#Characterise pipeline in terms of how many new annotated: longer, shorter,...
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

pdf("~/Desktop/DataAnalysisRiccio/Dec2016/axonal_remodelling/alternative/output.pipeline.pdf")
par(mfrow=c(1,2))
mp<-barplot(focus.iso,las=1,frame=F,col="white",ylim=c(0,30000),ylab="# novel 3' UTR isoforms")
mtext(side=3,line=0,text=focus.iso,at=mp,cex=0.5)
mp<-barplot(focus.txID,las=1,frame=F,col="white",ylim=c(0,12000),ylab="# Ensembl txID")
mtext(side=3,line=0,text=focus.txID,at=mp,cex=0.5)
dev.off()

maxdL.pos        <- (anno_ngf$maxL-anno_ngf$initL)[match(unique(anno_ngf$txID),anno_ngf$txID)]
maxdL.pos        <- maxdL.pos[maxdL.pos>0]
maxdL.neg        <- abs(anno_ngf$minL-anno_ngf$initL)[match(unique(anno_ngf$txID),anno_ngf$txID)]
maxdL.neg        <- maxdL.neg[maxdL.neg>0]


t2g                     <- read.csv("~/Desktop/DataAnalysisRiccio/t2g_biomaRt.csv")
t2g                     <- read.csv("~/Desktop/DataAnalsyisRiccio/t2g_biomaRt.csv")
t2g                     <- read.csv("/home/rluisier/data/Riccio/Exp_1/t2g_biomaRt.csv")
txID2GO                 <- tapply(t2g$ensembl_transcript_id,INDEX=t2g$go_id,FUN=function(x)return(x))
noNodes                 <- 300

#
# START NESTED FUNCTIONS
#
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

#
# END NESTED FUNCTIONS
#


#
# A. FILTERING
#

# A.1 Remove redundant transcript
id       <- paste(start(myUTR),end(myUTR),strand(myUTR),sep=".")#674 duplicated
myUTR    <- myUTR[-which(duplicated(id)),]
anno_ngf <- anno_ngf[-which(duplicated(id)),]

pdf(paste(outdir,"select_expressed_genes.pdf",sep=""))
par(mfrow=c(2,2))
mydat                <- anno_ngf[which(anno_ngf$is.conservative),match(c("NGF.axon.1.raw","NGF.axon.2.raw","NGF.cb.1.raw","NGF.cb.2.raw"),colnames(anno_ngf))]
mydatall             <- anno_ngf[,match(c("NGF.axon.1.raw","NGF.axon.2.raw","NGF.cb.1.raw","NGF.cb.2.raw"),colnames(anno_ngf))]
for(i in c(1:ncol(mydat))){
  mydat[,i] <- log2(mydat[,i]+1)
}

lims      <- apply(mydat,2,function(Z)2^SelectExpressed(dat=Z,frac.bg=0.95,frac.fg=0.1))
dev.off()

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

# A. Characterise the expression in each compartment

#txID.expression -- this cannot be used for counting given that there are duplicates in
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


pdf(paste(outdir,"no.genes.expressed_jul_2016.pdf",sep=""))

no.txID <- c(cb=length(unique(anno_ngf$txID[anno_ngf$NGF.cb.is.expressed.iso])),
             axons=length(unique(anno_ngf$txID[anno_ngf$NGF.axon.is.expressed.iso]) ))
no.iso  <- c(cb=length(unique(anno_ngf$uniqueID[anno_ngf$NGF.cb.is.expressed.iso])),
             axons=length(unique(anno_ngf$uniqueID[anno_ngf$NGF.axon.is.expressed.iso]) ))

mycols   <- c(rgb(23/255,71/255,120/255),rgb(119/255,192/255,68/255))
par(mfrow=c(1,2))
mp<- barplot(no.iso,las=1,col=mycols)
mtext(side=3,line=0,text=no.iso,cex=0.7,at=mp)
mp<- barplot(no.txID,las=1,col=mycols)
mtext(side=3,line=0,text=no.txID,cex=0.7,at=mp)
dev.off()

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

write.csv(myenrich[[1]],paste(outdir,"enrichment_cb.csv",sep=""))
write.csv(myenrich[[2]],paste(outdir,"enrichment_axons.csv",sep=""))
write.csv(myenrich[[3]],paste(outdir,"enrichment_neurons.csv",sep=""))
write.csv(myenrich[[4]],paste(outdir,"enrichment_cb.only.csv",sep=""))

enr1<- read.csv(paste(outdir,"enrichment_neurons_f.csv",sep=""))
enr2<- read.csv(paste(outdir,"enrichment_axons_f.csv",sep=""))
enr3<- read.csv(paste(outdir,"enrichment_cb.only_f.csv",sep=""))

pdf(paste(outdir,"enrichmt_compartmentstringent_jul_2016.pdf",sep=""))
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
dev.off()

# B. Length of the reliably expressed in each compartment
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


#In this version I consider all of them
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


pdf(paste(outdir,"forfig7_all_loose_july2016.pdf",sep=""))
par(mfrow=c(2,3),mar=c(3,4,3,3))
mycols <- c(rgb(254/255,218/255,0),rgb(150/255,150/255,150/255),"black")
mp<- barplot(unlist(lapply(L1,length)),las=1,col=mycols,frame=F)
mtext(side=3,line=0,text=unlist(lapply(L1,length)),at=mp,cex=0.6)
mtext(side=2,line=3,text="no.txID")

mycols <- c(rgb(254/255,218/255,0),rgb(150/255,150/255,150/255),"black")
mp<- barplot(unlist(lapply(L2,length)),las=1,col=mycols,frame=F)
mtext(side=3,line=0,text=unlist(lapply(L2,length)),at=mp,cex=0.5)
mtext(side=2,line=3,text="no.3' UTR isoforms")

barplot(with.multiple*100,col=mycols,las=1,frame=F)
mtext(side=2,line=3,text="fraction of txID with multiple iso")
barplot(distr.no.iso*100,beside=T,col=mycols[c(1,2)],las=1,,frame=F)
mtext(side=2,line=3,text="fraction of txID with multiple iso")
boxplot(utrL,outline=F,col=c("white",mycols[c(1,2)]),las=1,frame=F)
mtext(side=3,line=0,text=unlist(lapply(utrL,length)),at=c(1,2,3),cex=0.5)

dev.off()

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


no.genes.expressed.neurons              <- length(unique(anno_ngf$geneSymbol[NGF.axon.is.expressed.iso|NGF.cb.is.expressed.iso]))
no.genes.expressed.axons                <- length(unique(anno_ngf$geneSymbol[NGF.axon.is.expressed.iso]))


#txID.expression
axons                                   <- unique(anno_ngf$txID[anno_ngf$NGF.axon.is.expressed.iso])
cb                                      <- unique(anno_ngf$txID[anno_ngf$NGF.cb.is.expressed.iso])
neurons                                 <- union(axons,cb)
cb.only                                 <- setdiff(cb,axons)
axon.only                               <- setdiff(axons,cb)
both                                    <- intersect(axons,cb)

View(anno_ngf[which(anno_ngf$txID%in%axon.only),match(c("NGF.cb.1.raw","NGF.cb.2.raw","NGF.axon.1.raw","NGF.axon.2.raw"),colnames(anno_ngf))])
View(anno_ngf[setdiff(which(anno_ngf$uniqueID%in%setdiff(Axons.iso,CB.iso)),which(anno_ngf$txID%in%axon.only)),match(c("NGF.cb.1.raw","NGF.cb.2.raw","NGF.axon.1.raw","NGF.axon.2.raw"),colnames(anno_ngf))])

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

View(anno_ngf[which(anno_ngf$uniqueID%in%ax.ax),match(c("NGF.cb.1.raw","NGF.cb.2.raw","NGF.axon.1.raw","NGF.axon.2.raw"),colnames(anno_ngf))])
View(anno_ngf[which(anno_ngf$uniqueID%in%neur.ax),match(c("NGF.cb.1.raw","NGF.cb.2.raw","NGF.axon.1.raw","NGF.axon.2.raw"),colnames(anno_ngf))])


#Create Backgroun with what is expressed in neurons
mytx <- unique(as.character(anno_ngf$txID[anno_ngf$NGF.axon.is.expressed.iso|anno_ngf$NGF.cb.is.expressed.iso]))
mygs <- unique(as.character(anno_ngf$geneSymbol[anno_ngf$NGF.axon.is.expressed.iso|anno_ngf$NGF.cb.is.expressed.iso]))#12'529
myBG <- data.frame(txID=mytx,GS=c(mygs,rep(NA,length(mytx)-length(mygs))))

# A.2 Remove those transcript ID which have only one isoforms
anno_ngf     <- anno_ngf[anno_ngf$NGF.axon.is.expressed.iso|NGF.cb.is.expressed.iso,]
temp         <- as.data.frame(table(as.character(anno_ngf$txID)))
no.iso       <- temp[match(anno_ngf$txID,as.character(temp[,1])),2]
anno_ngf     <- anno_ngf[anno_ngf$no.iso>1,]#21377 rows
names(myUTR) <- myUTR$ID
myUTR        <- myUTR[match(anno_ngf$uniqueID,names(myUTR)),]#8'884 with only 1 isoforms

# A.3 Check that the closest neighbour is not within 300 nt
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
Fproxi                    <- which(distonext<300)#none; ok what expected
Fproxi

#
# B. ADD INFORMATION
#


# B.1 Remake ordering and imID
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


# B.2 Get the selection on which to focus for the analysis
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

length(unique(subanno$geneSymbol))

no.genes.expressed.neurons              <- length(unique(anno_ngf$geneSymbol[anno_ngf$NGF.axon.is.expressed.iso|anno_ngf$NGF.cb.is.expressed.iso]))
no.genes.expressed.axons                <- length(unique(anno_ngf$geneSymbol[anno_ngf$NGF.axon.is.expressed.iso]))
no.genes.expressed.axons.multiple.iso   <- length(unique(subanno$geneSymbol))
my.nos                                  <- c(no.genes.expressed.neurons,no.genes.expressed.axons,no.genes.expressed.axons.multiple.iso)
names(my.nos)                           <- c("neurons","axons","axons&no.iso>1")

no.genes.expressed.neurons.t              <- length(unique(anno_ngf$txID[anno_ngf$NGF.axon.is.expressed.iso|anno_ngf$NGF.cb.is.expressed.iso]))
no.genes.expressed.axons.t                <- length(unique(anno_ngf$txID[anno_ngf$NGF.axon.is.expressed.iso]))
no.genes.expressed.axons.multiple.iso.t   <- length(unique(subanno$txID))
my.nos.t                                  <- c(no.genes.expressed.neurons.t,no.genes.expressed.axons.t,no.genes.expressed.axons.multiple.iso.t)
names(my.nos.t)                           <- c("neurons","axons","axons&no.iso>1")


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

pdf(paste(outdir,"comparison_RUD_replicates.pdf",sep=""))
par(mfrow=c(2,2))
smoothScatter(RUD[,1],RUD[,2],frame=FALSE,xlab="ax1",ylab="ax2")
mtext(side=3,line=0,text=cor(RUD[,1],RUD[,2],method="spearman"))
smoothScatter(RUD[,3],RUD[,4],frame=FALSE,xlab="cb.1",ylab="cb.2")
mtext(side=3,line=0,text=cor(RUD[,3],RUD[,4],method="spearman"))
dev.off()



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

#my.proximal <- log2(1+sumRUD[ix2,])
#my.distal   <- log2(1+sumRUD[ix1,])
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


# C.2    Relative Proximal to distal site usage
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

unlist(lapply(selRUD1,function(Z)return(length(unique(subanno$geneSymbol[Z])))))
unlist(lapply(selRUD1,function(Z)return(sum(Z))))

#IX=IX+1
#selRUD1[[IX]]      <- padjRUD<0.01&dRUD<(-1)&seld1&diffPUD>(0.15)&opp.1#131-proximal shift so proximal is expressed in axons
#IX=IX+1
#selRUD1[[IX]]      <- padjRUD<0.01&dRUD>(1)&seld2&diffPUD<(-0.15)&opp.2#83-distal shift so distal expressed in axons

#IX=IX+1
#selRUD1[[IX]]      <- padjRUD<0.01&dRUD<(-1)&seld1&1&diff.rel.proximal.usage<(-0.15)&rel.proximal.usage[,2]<=Lim1&opp.1#78, less that 20% usage of the proximal in CB; therefore can be considered as absent
#IX=IX+1
#selRUD1[[IX]]      <- padjRUD<0.01&dRUD>(1)&seld2&diff.rel.proximal.usage>(0.15)&rel.proximal.usage[,2]>=Lim2&opp.2#34,less than 20% usage of the distal in CB; therefore can be considered as absent

my.no <- unlist(lapply(selRUD1,function(Z)return(length(unique(subanno$txID[Z])))))

pdf("~/Desktop/DataAnalysisRiccio/Dec2016/axonal_remodelling/alternative/no.uniqs_80.pdf")
mp<- barplot(my.no[c(5,6,9,10)],las=1)
mtext(side=3,line=0,text=my.no[c(5,6,9,10)],at=mp)
dev.off()

View(cbind(as.character(subanno$GS),diff.rel.proximal.usage,rel.proximal.usage)[selRUD1[[5]],])
View(cbind(as.character(subanno$GS),diff.rel.proximal.usage,rel.proximal.usage)[selRUD1[[6]],])
View(cbind(as.character(subanno$GS),diff.rel.proximal.usage,rel.proximal.usage))


mylist <- lapply(selRUD1,function(Z)return(unique(subanno$txID[Z])))

setdiff(cand$txID,mylist[[5]])
intersect(cand$txID,mylist[[5]])

View(cbind(as.character(subanno$GS),diff.rel.proximal.usage,rel.proximal.usage)[subanno$txID%in%setdiff(cand$txID,mylist[[9]]),])
cand <- read.csv("~/Desktop/DataAnalysisRiccio/Riccio_exp1_analysis/candidates_rem.csv")

View(cbind(as.character(subanno$GS),diff.rel.proximal.usage,rel.proximal.usage)[subanno$txID%in%setdiff(cand$txID,mylist[[9]]),])


View(subanno[selRUD1[[9]],])
View(subanno[selRUD1[[10]],])


ix <- which(subanno$GS=="Sub1")
ix <- which(subanno$GS=="Vps36")


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

#Check the level of expression of the proximal in CB for those which are not expressed in axons; they are way less expressed
ix.proxi               <- match(myIDs[,1],rownames(myMean))
proi                   <- lapply(distal.focus,function(Z)return(anno_ngf$imIDp[match(Z,anno_ngf$txID)]))
ix.proxi.distal        <- lapply(proi,function(Z)return(match(Z,rownames(myMean))))
cb.e                   <- list(all=myMean[ix.proxi,2],proxi.expressed=myMean[ix.proxi.distal[[1]],2],proxi.not.expressed=myMean[ix.proxi.distal[[2]],2])
boxplot(cb.e)
multidensity(cb.e)
abline(v=8)
abline(v=10)
plot(pch=19,myMean[,c(2,1)],cex=0.3,col=rgb(0,0,0,0.1))
abline(v=8)
abline(v=10)

par(mfrow=c(2,2),mar=c(3,4,3,3))
#layout(matrix(c(1,4,4,2,4,4,3,4,4),4,3,byrow = TRUE))
mp=barplot(c(distal.and.distalshift=lapply(distal.focus,length)[[2]],prox.only.and.proxi.shift=lapply(proxi.focus,length)[[3]]),col=c(rgb(30/255,36/255,83/255),rgb(154/255,162/255,197/255)),las=1,ylab="# Ensembl txID")
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
mycols[selS&subanno$txID%in%proxi.focus[[3]]] <- "red"
mycols[selL&subanno$txID%in%distal.focus[[2]]]<- "red"


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



### Compare with with Mass-Spec data
foi   <- list.files("~/Desktop/DataAnalsyisRiccio/MSdata")
foi   <- foi[grep(foi,pattern=".csv")]
data  <- lapply(paste("~/Desktop/DataAnalsyisRiccio/MSdata/",foi,sep=""),function(x)return(read.csv(x)))
descr <- c(as.character(data[[1]][,2]),as.character(data[[3]][,2]))
ExtractGeneName <- function(test){
  temp<- unlist(strsplit(test,split=" "))[grep(pattern="GN=",unlist(strsplit(test,split=" ")))]
  return(gsub(temp,pattern="GN=",repl=""))
}
GS    <- unlist(lapply(descr,ExtractGeneName))
for(i in c(1:length(data))){
  data[[i]]$GS <- GS[match(as.character(data[[i]][,2]),descr)]
  print(i)
}

#Some statistics on the number of protein detected in each compartments
ax.prot                         <- unique(as.character(data[[1]]$GS))
cb.prot                         <- unique(as.character(data[[3]]$GS))
neur                            <- intersect(ax.prot,cb.prot)
ax                              <- setdiff(ax.prot,cb.prot)
cb                              <- setdiff(cb.prot,ax.prot)
pdf(paste(outdir,"no.prot.detected.pdf",sep=""))
par(mfrow=c(2,2))
no.prot                         <- c(length(cb.prot),length(ax.prot))
mycols                          <- c(rgb(22/255,70/255,120/255),rgb(117/255,191/255,68/255))
mp<- barplot(no.prot,col=mycols,las=1)
mtext(side=3,line=0,text=no.prot,cex=0.6,at=mp)
mtext(side=2,line=3,text="# protein detected",cex=0.6)

no.prot                         <- c(length(neur),length(cb),length(ax))
mycols                          <- c(rgb(253/255,218/255,0/255),rgb(150/255,150/255,150/255),"black")
mp                              <- barplot(no.prot,col=mycols,las=1)
mtext(side=3,line=0,text=no.prot,cex=0.6,at=mp)
mtext(side=2,line=3,text="# protein detected",cex=0.5)

dev.off()

whereProt                        <- rep("neur",length(unique(GS)))
names(whereProt)                 <- unique(GS)
whereProt[names(whereProt)%in%ax]<- "axons"
whereProt[names(whereProt)%in%cb]<- "cb"

detect.prot                      <- whereProt[match(anno_ngf$geneSymbol,names(whereProt))]
detect.prot[is.na(detect.prot)]  <- "not detected"
anno_ngf$detect.prot             <- detect.prot


myGS <- unique(anno_ngf$geneSymbol[anno_ngf$])


pie(table(anno_ngf$detect.prot[match(subanno$txID[selRUD1[[1]]],anno_ngf$txID)]))
pie(table(anno_ngf$detect.prot[match(subanno$txID[selRUD1[[2]]],anno_ngf$txID)]))

pie(table(anno_ngf$detect.prot[match(subanno$txID[selRUD1[[3]]],anno_ngf$txID)]))
pie(table(anno_ngf$detect.prot[match(subanno$txID[selRUD1[[4]]],anno_ngf$txID)]))

pie(table(anno_ngf$detect.prot[match(subanno$txID[selRUD1[[5]]],anno_ngf$txID)]))
pie(table(anno_ngf$detect.prot[match(subanno$txID[selRUD1[[6]]],anno_ngf$txID)]))

table(myGS%in%subanno$geneSymbol[selRUD1[[3]]],myGS%in%neur)

strong <- read.csv("~/Desktop/DataAnalsyisRiccio/Dec2016/axonal_remodelling/alternative/candidates.csv",header=F,colClasses = "character")

vec<-data.frame(table(anno_ngf$detect.prot[match(strong[,1],anno_ngf$txID)]))

cbind(anno_ngf$detect.prot,as.character(anno_ngf$txID))[match(subanno$txID[selRUD1[[5]]],anno_ngf$txID),]

mySel.IDs   <- do.call(what=cbind,args=lapply(selRUD1,function(X)return(as.character(subanno$uniqueID[X]))))
for(i in c(1:ncol(mySel.IDs))){
  mySel.IDs[duplicated(mySel.IDs[,i]),i]<- NA
}
mySel.GS   <- do.call(what=cbind,args=lapply(selRUD1,function(X)return(as.character(subanno$geneSymbol[X]))))
for(i in c(1:ncol(mySel.GS))){
  mySel.GS[duplicated(mySel.GS[,i]),i]<- NA
}
shiftType                                            <- rep("null",nrow(anno_ngf))
shiftType[which(anno_ngf$uniqueID%in%mySel.IDs[,3])] <- "proximal.shift"
shiftType[which(anno_ngf$uniqueID%in%mySel.IDs[,4])] <- "distal.shift"
shiftTypeGS                                          <- rep("null",nrow(anno_ngf))
shiftTypeGS[which(anno_ngf$GS%in%mySel.GS[[7]])]     <- "proximal.shift"
shiftTypeGS[which(anno_ngf$GS%in%mySel.GS[[8]])]     <- "distal.shift"
shiftTypeGS[which(anno_ngf$GS%in%intersect(mySel.GS[[7]],mySel.GS[[8]]))]<-"both"
anno_ngf$shiftType                                   <- shiftType
anno_ngf$shiftTypeGS                                 <- shiftTypeGS


ix                               <- match(unique(anno_ngf$geneSymbol),anno_ngf$geneSymbol)
table(anno_ngf$detect.prot[ix],anno_ngf$shiftTypeGS[ix])






#Prepare for the output
mySelOverall.GS       <- do.call(what=cbind,args=lapply(selRUD1,function(X)return(unique(as.character(subanno$geneSymbol[X])))))
mySelOverall.txID     <- do.call(what=cbind,args=lapply(selRUD1,function(X)return(unique(as.character(subanno$txID[X])))))
mySelOverall.uniqueID <- do.call(what=cbind,args=lapply(selRUD1,function(X)return(unique(as.character(subanno$uniqueID[X])))))

for(i in c(1:ncol(mySelOverall))){
  mySelOverall.GS[duplicated(mySelOverall.GS[,i]),i]             <-NA
  mySelOverall.txID[duplicated(mySelOverall.txID[,i]),i]         <-NA
  mySelOverall.uniqueID[duplicated(mySelOverall.uniqueID[,i]),i] <-NA
}

colnames(mySelOverall.GS) <-colnames(mySelOverall.txID) <-colnames(mySelOverall.uniqueID) <- paste(rep(c("selShort","selLong"),3),c("sel0","sel0","sel1","sel2","sel3","sel3"),sep=".")
write.table(mySelOverall.GS,file=paste(outdir,"mySelOverall.GS.txt",sep=""),col.names = TRUE,row.names = FALSE,quote=F,sep="\t")
write.table(mySelOverall.txID,file=paste(outdir,"mySelOverall.txID.txt",sep=""),col.names = TRUE,row.names = FALSE,quote=F,sep="\t")
write.table(mySelOverall.uniqueID,file=paste(outdir,"mySelOverall.uniqueID.txt",sep=""),col.names = TRUE,row.names = FALSE,quote=F,sep="\t")

pdf("~/Desktop/DataAnalysisRiccio/Dec2016/axonal_remodelling/final_analysis/plot_final_enrich_p.pdf")
foi <- list.files("~/Desktop/DataAnalysisRiccio/Dec2016/axonal_remodelling/final_analysis/")
foi <- foi[intersect(grep(foi,pattern="_f"),grep(foi,pattern="F0"))][c(2:3)]
final_enrich <- lapply(foi,function(Z)return(read.csv(paste("~/Desktop/DataAnalysisRiccio/Dec2016/axonal_remodelling/final_analysis/",Z,sep=""))))
PlotFinalEnrich(final_enrich)

foi <- list.files("~/Desktop/DataAnalysisRiccio/Dec2016/axonal_remodelling/final_analysis/")
foi <- foi[intersect(grep(foi,pattern="_f"),grep(foi,pattern="F1"))][c(2:3)]
final_enrich <- lapply(foi,function(Z)return(read.csv(paste("~/Desktop/DataAnalysisRiccio/Dec2016/axonal_remodelling/final_analysis/",Z,sep=""))))
PlotFinalEnrich(final_enrich)

foi <- list.files("~/Desktop/DataAnalysisRiccio/Dec2016/axonal_remodelling/final_analysis/")
foi <- foi[intersect(grep(foi,pattern="_f"),grep(foi,pattern="F2"))][c(2:3)]
final_enrich <- lapply(foi,function(Z)return(read.csv(paste("~/Desktop/DataAnalysisRiccio/Dec2016/axonal_remodelling/final_analysis/",Z,sep=""))))
PlotFinalEnrich(final_enrich)

foi <- list.files("~/Desktop/DataAnalysisRiccio/Dec2016/axonal_remodelling/final_analysis/")
foi <- foi[intersect(grep(foi,pattern="_f"),grep(foi,pattern="F3"))][c(2:3)]
final_enrich <- lapply(foi,function(Z)return(read.csv(paste("~/Desktop/DataAnalysisRiccio/Dec2016/axonal_remodelling/final_analysis/",Z,sep=""))))
PlotFinalEnrich(final_enrich)

dev.off()


--> I should select some examples to illustrate such as Rab, translation, ubiquitination

#Long
GO:0016197 activation of JUN kinase activity
GetOI(mygoID="GO:0016197",sampleGO.RUD1[[7]][[2]])

GO:0048812 neuron projection morphogenesis
GetOI(mygoID="GO:0048812",sampleGOLong[[2]])
GO:0007265 Ras protein signal transduction
GetOI(mygoID="GO:0007265",sampleGOLong[[2]])
GO:0015031 protein transport
GetOI(mygoID="GO:0015031",sampleGOLong[[2]])


#Short
GO:0016197 endosomal transport
GetOI(mygoID="GO:0016197",sampleGOShort[[2]])
GO:0060675 ureteric bud morphogenesis
GetOI(mygoID="GO:0060675",sampleGOShort[[2]])
GO:0006417 regulation of translation
GetOI(mygoID="GO:0006417",sampleGOShort[[2]])
GO:0032482 Rab protein signal transduction
GetOI(mygoID="GO:0032482",sampleGOShort[[2]])
GO:0051301 cell division
GetOI(mygoID="GO:0051301",sampleGOShort[[2]])
GO:0006886 intracellular protein transport
GetOI(mygoID="GO:0006886",sampleGOShort[[2]])
GO:0006915 apoptotic process
GetOI(mygoID="GO:0006915",sampleGOShort[[2]])


###### To plot, get the values of interest
ix1                  <- c(1:nrow(anno_ngf))
ix2                  <- match(as.character(anno_ngf$imIDp),as.character(anno_ngf$uniqueID))
mydat                <- anno_ngf[,match(c("NGF.axon.1.raw","NGF.axon.2.raw","NGF.cb.1.raw","NGF.cb.2.raw"),colnames(anno_ngf))]
myProximal           <- apply(mydat,2,function(x)return(x[ix2]))
myDistal             <- apply(mydat,2,function(x)return(x[ix1]))
myProxi.avg          <- t(apply(myProximal,1,FUN=function(X)return(tapply(X,INDEX=factor(c("axon","axon","cb","cb"),levels=c("axon","cb")),FUN=mean))))
myDistal.avg        <- t(apply(myDistal,1,FUN=function(X)return(tapply(X,INDEX=factor(c("axon","axon","cb","cb"),levels=c("axon","cb")),FUN=mean))))
myProxi.sd          <- t(apply(myProximal,1,FUN=function(X)return(tapply(X,INDEX=factor(c("axon","axon","cb","cb"),levels=c("axon","cb")),FUN=sd))))
myDistal.sd        <- t(apply(myDistal,1,FUN=function(X)return(tapply(X,INDEX=factor(c("axon","axon","cb","cb"),levels=c("axon","cb")),FUN=sd))))

RDU                <- do.call(what=cbind,lapply(c(1:ncol(myProximal)),FUN=function(X)return(myDistal[,X]/(myDistal[,X]+myProximal[,X]))))
RDU.avg            <- t(apply(RDU,1,FUN=function(X)return(tapply(X,INDEX=factor(c("axon","axon","cb","cb"),levels=c("axon","cb")),FUN=mean))))
RDU.sd            <- t(apply(RDU,1,FUN=function(X)return(tapply(X,INDEX=factor(c("axon","axon","cb","cb"),levels=c("axon","cb")),FUN=sd))))

myDat <- data.frame(GS=anno_ngf$geneSymbol,
                    txID=anno_ngf$uniqueID,
                    proxiID=anno_ngf$uniqueID[ix2],
                    distalID=anno_ngf$uniqueID[ix1],

                    mean.cb.proxi=myProxi.avg[,2],
                    sd.cb.proxi=myProxi.sd[,2],
                    mean.cb.distal=myDistal.avg[,2],
                    sd.cb.distal=myDistal.sd[,2],
                    mean.rdu.cb=RDU.avg[,2],
                    sd.rdu.ax=RDU.sd[,2],

                    mean.ax.proxi=myProxi.avg[,1],
                    sd.ax.proxi=myProxi.sd[,1],
                    mean.ax.distal=myDistal.avg[,1],
                    sd.ax.distal=myDistal.sd[,1],
                    mean.rdu.ax=RDU.avg[,1],
                    sd.rdu.ax=RDU.sd[,1]
)


ix1                  <- c(1:nrow(anno_ngf))
ix2                  <- match(as.character(anno_ngf$imIDp),as.character(anno_ngf$uniqueID))
mydat                <- anno_ngf[,match(c("NGF.axon.1.raw","NGF.axon.2.raw","NGF.cb.1.raw","NGF.cb.2.raw"),colnames(anno_ngf))]
myProximal           <- apply(mydat,2,function(x)return(log2(1+x[ix2])))
myDistal             <- apply(mydat,2,function(x)return(log2(1+x[ix1])))
myProxi.avg          <- t(apply(myProximal,1,FUN=function(X)return(tapply(X,INDEX=factor(c("axon","axon","cb","cb"),levels=c("axon","cb")),FUN=mean))))
myDistal.avg        <- t(apply(myDistal,1,FUN=function(X)return(tapply(X,INDEX=factor(c("axon","axon","cb","cb"),levels=c("axon","cb")),FUN=mean))))
myProxi.sd          <- t(apply(myProximal,1,FUN=function(X)return(tapply(X,INDEX=factor(c("axon","axon","cb","cb"),levels=c("axon","cb")),FUN=sd))))
myDistal.sd        <- t(apply(myDistal,1,FUN=function(X)return(tapply(X,INDEX=factor(c("axon","axon","cb","cb"),levels=c("axon","cb")),FUN=sd))))

RDU                <- do.call(what=cbind,lapply(c(1:ncol(myProximal)),FUN=function(X)return(myDistal[,X]/(myDistal[,X]+myProximal[,X]))))
RDU.avg            <- t(apply(RDU,1,FUN=function(X)return(tapply(X,INDEX=factor(c("axon","axon","cb","cb"),levels=c("axon","cb")),FUN=mean))))
RDU.sd            <- t(apply(RDU,1,FUN=function(X)return(tapply(X,INDEX=factor(c("axon","axon","cb","cb"),levels=c("axon","cb")),FUN=sd))))

myDat2 <- data.frame(GS=anno_ngf$geneSymbol,
                     txID=anno_ngf$uniqueID,
                     proxiID=anno_ngf$uniqueID[ix2],
                     distalID=anno_ngf$uniqueID[ix1],

                     mean.cb.proxi=myProxi.avg[,2],
                     sd.cb.proxi=myProxi.sd[,2],
                     mean.cb.distal=myDistal.avg[,2],
                     sd.cb.distal=myDistal.sd[,2],
                     mean.rdu.cb=RDU.avg[,2],
                     sd.rdu.ax=RDU.sd[,2],

                     mean.ax.proxi=myProxi.avg[,1],
                     sd.ax.proxi=myProxi.sd[,1],
                     mean.ax.distal=myDistal.avg[,1],
                     sd.ax.distal=myDistal.sd[,1],
                     mean.rdu.ax=RDU.avg[,1],
                     sd.rdu.ax=RDU.sd[,1]
)


oi <- which(myDat$GS=="Msantd3")
View(myDat[oi,])
View(myDat2[oi,])

###### To plot, get the values of interest






### Analysis when considering only 2 isoforms per txID

load("/Users/luisie01/Desktop/DataAnalsyisRiccio/Dec2016/utrid/anno_ngf_stringent/anno_ngf_March_12.RData")#Contain "anno_ngf","ngfGRS"
myUTR        <- import.gff("/Users/luisie01/Desktop/DataAnalsyisRiccio/Dec2016/utrid/APA_stringent/Lngf_sub.gtf",format="gtf",asRangedData = FALSE)
names(myUTR) <- as.character(myUTR$ID)
ix1          <- match(names(myUTR),anno_ngf$uniqueID)
myUTR        <- myUTR[!is.na(ix1),]
ix1          <- match(names(myUTR),anno_ngf$uniqueID)
anno_ngf    <- anno_ngf[ix1,]
outdir       <-  "/Users/luisie01/Desktop/DataAnalsyisRiccio/Dec2016/axonal_remodelling/"
load("/Users/luisie01/Desktop/DataAnalsyisRiccio/GRanges_comprehensive_transcriptome_v78_rn5_23022014.RData")



#
# A. FILTERING
#

# A.1 Remove redundant transcript
id       <- paste(start(myUTR),end(myUTR),strand(myUTR),sep=".")#674 duplicated
myUTR    <- myUTR[-which(duplicated(id)),]
anno_ngf <- anno_ngf[-which(duplicated(id)),]

# A.2 Remove those transcript ID which have only one isoforms
temp    <- as.data.frame(table(as.character(anno_ngf$txID)))
no.iso  <- temp[match(anno_ngf$txID,as.character(temp[,1])),2]
myUTR   <- myUTR[no.iso>1,]#8'884 with only 1 isoforms
anno_ngf<- anno_ngf[no.iso>1,]

# A.3 Check that the closest neighbour is not within 300 nt
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
Fproxi                    <- which(distonext<300)#none; ok what expected
Fproxi


# A.4 Remove those potential spurious ones
torm     <- which(mcols(myUTR)$is.pas=="FALSE"&mcols(myUTR)$isoform!=0)
myUTR    <- myUTR[-torm,]
anno_ngf <- anno_ngf[-torm,]
temp    <- as.data.frame(table(as.character(anno_ngf$txID)))
no.iso  <- temp[match(anno_ngf$txID,as.character(temp[,1])),2]
myUTR   <- myUTR[no.iso>1,]
anno_ngf<- anno_ngf[no.iso>1,]
no.iso  <- no.iso[no.iso>1]

# A.5 When several isoforms select the most expressed in CB
tempL                     <- anno_ngf$newL
names(tempL)              <- anno_ngf$uniqueID
tempG                     <- as.factor(as.character(anno_ngf$txID))
myord                     <- tapply(tempL,INDEX=tempG,function(x)return(cbind(names(x)[sort(as.numeric(as.character(x)),index.return=T,decreasing=F)$ix],c(1:length(x)))))
test                      <- do.call(what=rbind,args=myord)
iso_ordered               <- as.numeric(test[match(anno_ngf$uniqueID,test[,1]),2])
gp1                       <- anno_ngf[no.iso==2|iso_ordered==1,]
gp2                       <- anno_ngf[no.iso>2&iso_ordered>1,]
tempL                     <- gp2$NGF.cb.log2.sum
names(tempL)              <- gp2$uniqueID
mymax                     <- tapply(tempL,INDEX=factor(as.character(gp2$txID)),function(x)return(names(x)[which(x==max(x))[1]]))
sel                       <- match(mymax,gp2$uniqueID)
gp2                       <- gp2[sel,]

anno_ngf                  <- rbind(gp1,gp2)
myUTR                     <- myUTR[match(anno_ngf$uniqueID,myUTR$ID),]
temp                      <- as.data.frame(table(as.character(anno_ngf$txID)))
no.iso                    <- temp[match(anno_ngf$txID,as.character(temp[,1])),2]

#
# B. ADD INFORMATION
#


# B.1 Remake ordering and imID
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

# B.2 Add key column in anno_ngf:expressed in at least one sample
#NGF.cb.is.expressed.iso                    <- anno_ngf$detect.iso.ngf%in%c("cb.cb","neur.cb","neur.neur")
#NGF.axon.is.expressed.iso                  <- anno_ngf$detect.iso.ngf%in%c("neur.neur","neur.ax","ax.ax")
NGF.axon.is.expressed.iso                   <- anno_ngf$NGF.axon.1.isexpressed|anno_ngf$NGF.axon.2.isexpressed
NGF.cb.is.expressed.iso                     <- anno_ngf$NGF.cb.1.isexpressed|anno_ngf$NGF.cb.2.isexpressed
anno_ngf                                    <- data.frame(anno_ngf,NGF.cb.is.expressed.iso,NGF.axon.is.expressed.iso)

# B.3 Get the order and the ID of the longest detected one
temp1                 <- tapply(anno_ngf$iso_ordered[anno_ngf$NGF.axon.is.expressed.iso],factor(as.character(anno_ngf$txID[anno_ngf$NGF.axon.is.expressed.iso])),FUN=function(x)return(max(as.numeric(x))))
long.ax                <- as.integer(temp1[match(anno_ngf$txID,names(temp1))])
long.ax[is.na(long.ax)]<- 1
temp1                <- tapply(anno_ngf$iso_ordered[anno_ngf$NGF.cb.is.expressed.iso],factor(as.character(anno_ngf$txID[anno_ngf$NGF.cb.is.expressed.iso])),FUN=max)
long.cb              <- as.integer(temp1[match(anno_ngf$txID,names(temp1))])
long.cb[is.na(long.cb)]<- 1

anno_ngf$ord.longest <- apply(cbind(long.ax,long.cb) ,1,max)
temp                 <- anno_ngf[,match(c("ord.longest","iso_ordered","uniqueID"),colnames(anno_ngf))]
myid<- tapply(c(1:nrow(temp)),INDEX=anno_ngf$txID,FUN=function(IX){
  return(as.character((temp$uniqueID[IX])[match(unique(temp$ord.longest[IX]),as.integer(temp$iso_ordered[IX]))]))
})
myid               <- myid[match(anno_ngf$txID,names(myid))]
anno_ngf$longestID <- myid

# B.4 Get distance to longest expressed
anno_ngf$distToLongest <- anno_ngf$newL[match(anno_ngf$longestID,anno_ngf$uniqueID)]-anno_ngf$newL

