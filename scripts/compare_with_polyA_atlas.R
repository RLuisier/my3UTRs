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


load("~/Desktop/DataAnalsyisRiccio/GRanges_comprehensive_transcriptome_v78_rn5_23022014.RData")
outdir <-  "~/Desktop/DataAnalsyisRiccio/Rebuttal_May_2017/"
dir.create(outdir)
load("~/Desktop/DataAnalsyisRiccio/Dec2016/utrid/anno_ngf_stringent/anno_ngf_March_12.RData")#anno_ngf[41204]

#load("/home/rluisier/data/Riccio/Exp_1/Dec2016/Coverage/utrCov_stringent/myCov_March_12_2016_500.RData")

myUTR        <- import.gff("~/Desktop/DataAnalsyisRiccio/Dec2016/utrid/APA_stringent/Lngf.gtf",format="gtf")
names(myUTR) <- as.character(myUTR$ID)
ix1          <- match(names(myUTR),anno_ngf$uniqueID) 
myUTR        <- myUTR[!is.na(ix1),]
ix1          <- match(names(myUTR),anno_ngf$uniqueID) 
anno_ngf     <- anno_ngf[ix1,]
myUTR.new    <- myUTR[-grep(myUTR$ID,pattern="\\.0"),]


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

SelectSurroundingEnd <- function(myGR=myUTR.new,before=10,after=10){
  NEG         <- as.character(strand(myGR))=="-"
  mystart     <- end(myGR)-before
  myend       <- end(myGR)+after
  mystart[NEG]<- start(myGR)[NEG]-after
  myend[NEG]  <- start(myGR)[NEG]+before
  end(myGR)    <- myend
  start(myGR)  <- mystart
  return(myGR)
}

GetFractionOverlap <- function(utrGR=myUTR.1,target=polya_db){
  return(length(unique(subjectHits(findOverlaps(query=target,subject=utrGR,ignore.strand=F))))/length(utrGR))
}

GetFractionRecovered <- function(utrGR=SelectSurroundingEnd(myGR=myUTR,before=10,after=10),target=polya_db){
  known <- utrGR[grep(utrGR$ID,pattern="\\.0"),]
  novel <- utrGR[-grep(utrGR$ID,pattern="\\.0"),] 
  absent<- target[-unique(subjectHits(findOverlaps(query=known,subject=target,ignore.strand=F))),]
  return(length(unique(subjectHits(findOverlaps(query=novel,subject=absent,ignore.strand=F))))/length(absent))
}

Extract_utr <- function(xeno=rn6){
  myEnd.pos       <- as.numeric(end(xeno))
  myStart.neg     <- as.numeric(start(xeno))
  mymaxEnd.pos    <- tapply(myEnd.pos,xeno$transcript_id,max)+10
  myStart.pos     <- mymaxEnd.pos-21
  myminStart.neg  <- tapply(myStart.neg,xeno$transcript_id,min)-10
  myEnd.neg       <- myminStart.neg+21
  POS             <- tapply(as.character(strand(xeno)),xeno$transcript_id,function(x)return(sum(x=="+")!=0))
  Chrom           <- as.character(seqnames(xeno))[match(names(POS),mcols(xeno)$transcript_id)]
  Strand           <- rep("+",length(POS))
  Strand[!POS]     <- "-"
  myStart         <- myStart.pos
  myStart[!POS]   <- myminStart.neg[!POS]
  myEnd           <- mymaxEnd.pos
  myEnd[!POS]     <- myEnd.neg[!POS]
  subxeno         <- GRanges(seqnames =Rle(Chrom),ranges =IRanges(start=myStart,end=myEnd),strand = Rle(Strand),transcript_id=names(POS)
  )
  return(subxeno)
}

MergeGRangesObject <- function(gr1=polya_db,gr2=refseq.rn5){
  Chrom   <- c(seqnames(gr1),seqnames(gr2))
  Strand  <- Rle(c(as.character(strand(gr1)),as.character(strand(gr2))))
  myRanges<- c(ranges(gr1),ranges(gr2))
  return(GRanges(seqnames =Chrom,ranges =myRanges,strand = Strand))
}

mydat                <- anno_ngf[which(anno_ngf$is.conservative),match(c("NGF.axon.1.raw","NGF.axon.2.raw","NGF.cb.1.raw","NGF.cb.2.raw"),colnames(anno_ngf))]
mydatall             <- anno_ngf[,match(c("NGF.axon.1.raw","NGF.axon.2.raw","NGF.cb.1.raw","NGF.cb.2.raw"),colnames(anno_ngf))]
for(i in c(1:ncol(mydat))){
  mydat[,i] <- log2(mydat[,i]+1)
}
require(mclust)
lims                <- apply(mydat,2,function(Z)2^SelectExpressed(dat=Z,frac.bg=0.95,frac.fg=0.1))
tempsel             <- apply(mydatall,2,function(Z)return(Z>=mean(lims[1,c(1,2)])))
soft.sel            <- cbind(tempsel[,1]&tempsel[,2],tempsel[,3]&tempsel[,4])
tempsel             <- do.call(what=cbind,lapply(c(1:4),function(Z)return(mydatall[,Z]>=lims[2,Z])))
hard.sel            <- cbind(tempsel[,1]|tempsel[,2],tempsel[,3]|tempsel[,4])
final.sel           <- cbind(soft.sel[,1]|hard.sel[,1],soft.sel[,2]|hard.sel[,2])

NGF.axon.is.expressed.iso                   <- final.sel[,1]
NGF.cb.is.expressed.iso                     <- final.sel[,2]
anno_ngf$NGF.axon.is.expressed.iso          <- final.sel[,1]
anno_ngf$NGF.cb.is.expressed.iso            <- final.sel[,2]
txoi                                        <- unique(anno_ngf$txID[anno_ngf$NGF.axon.is.expressed.iso|anno_ngf$NGF.cb.is.expressed.iso])



AnalysisFractionOverlap <- function(siteoi)return(
  c(100*GetFractionOverlap(utrGR=SelectSurroundingEnd(myGR=myUTR.new,before=10,after=10),target=siteoi),
    100*GetFractionOverlap(utrGR=SelectSurroundingEnd(myGR=myUTR.new,before=50,after=50),target=siteoi),
    100*GetFractionOverlap(utrGR=SelectSurroundingEnd(myGR=myUTR.new,before=100,after=100),target=siteoi),
    100*GetFractionOverlap(utrGR=SelectSurroundingEnd(myGR=myUTR.new,before=150,after=150),target=siteoi),
    100*GetFractionOverlap(utrGR=SelectSurroundingEnd(myGR=myUTR.new,before=200,after=200),target=siteoi)))


AnalysisFractionRecovered <- function(siteoi)return(
  c(100*GetFractionRecovered(utrGR=SelectSurroundingEnd(myGR=myUTR,before=10,after=10),target=siteoi),
    100*GetFractionRecovered(utrGR=SelectSurroundingEnd(myGR=myUTR,before=50,after=50),target=siteoi),
    100*GetFractionRecovered(utrGR=SelectSurroundingEnd(myGR=myUTR,before=100,after=100),target=siteoi)
  )
)



PlotFractionValidated <- function(siteoi=polya_db,name="polyAdb"){
  frac<- c(100*GetFractionOverlap(utrGR=SelectSurroundingEnd(myGR=myUTR.new,before=10,after=10),target=siteoi),
           100*GetFractionOverlap(utrGR=SelectSurroundingEnd(myGR=myUTR.new,before=50,after=50),target=siteoi),
           100*GetFractionOverlap(utrGR=SelectSurroundingEnd(myGR=myUTR.new,before=100,after=100),target=siteoi),
           100*GetFractionOverlap(utrGR=SelectSurroundingEnd(myGR=myUTR.new,before=150,after=150),target=siteoi),
           100*GetFractionOverlap(utrGR=SelectSurroundingEnd(myGR=myUTR.new,before=200,after=200),target=siteoi))
  mp<- barplot(frac,las=1,frame=FALSE,main="",ylab="",ylim=c(0,round(max(frac)+5)),cex.axis=0.7,cex=0.7,col=c("white","#E0E1E0","#989998","#5A5B5A","black"))
  text(x=mp,cex=0.7,y=ceiling(frac)+1,label=round(frac,digits=1))
  mtext(side=2,line=2,cex=0.7,text=paste("% validated by ",name))
}

#Import files of interest
#1. PolyA_db
polya_db                <- import.gff("~/Desktop/DataAnalsyisRiccio/Dec2016/rn.gtf")
AnalysisFractionRecovered(polya_db)#What is plotted on Figure 1D -- fraction of polyA absent from Rn5 and recovered


#2. APADB
test                    <-read.table("~/igv/Rn5/APADB/hg19.apadb_v2_final_rn5.bed",sep="\t")[,c(1,2,3,6)]
colnames(test)          <-c("seqnames","start","end","strand")
apadb.h19              <- makeGRangesFromDataFrame(test,ignore.strand=FALSE,seqnames.field="seqnames",start.field="start",end.field="end",strand.field="strand")
test                    <-read.table("~/igv/Rn5/APADB/mm10.apadb_v2_final_rn5.bed",sep="\t")[,c(1,2,3,6)]
colnames(test)          <-c("seqnames","start","end","strand")
apadb.mm10              <- makeGRangesFromDataFrame(test,ignore.strand=FALSE,seqnames.field="seqnames",start.field="start",end.field="end",strand.field="strand")
apadb <- c(apadb.h19,apadb.mm10)


#3. RefSeq
refseq.rn5              <- import.gff("~/igv/Riccio/RefSeq_rn5_UCSC_27062016.gtf",format="gtf")
refseq.rn6              <- import.gff("~/igv/Rn6/RefSeq/refSeq_rn6.Rn5.sorted.gtf",format="gtf")
refseq                  <- MergeGRangesObject(Extract_utr(refseq.rn5),Extract_utr(refseq.rn6))

#4.Ensembl Rn6
rn6                     <- import.gff("~/igv/Rn6/Ensembl/Rattus_norvegicus.Rnor_6.0.87.copy.sorted.Rn5.sorted.gtf")



#5. RefSeq xeno.musmusc
xeno.musmusc            <- import.gff("~/igv/xenoRefSeqrn5MusMusculus.gtf",format="gtf")
#6. RefSeq xeno.h19
xeno.h19                <- import.gff("~/igv/xenoRefSeqrn5HomoSapiens.gtf",format="gtf")
#7. RefSeq xeno.all
xeno.all                <- import.gff("~/igv/xenoRefSeqrn5.sorted.gtf",format="gtf")#all xeno

genome.annotation       <- MergeGRangesObject(MergeGRangesObject(Extract_utr(xeno.all),Extract_utr(rn6)),refseq)


#8. PolyA_unibasel -- mm10
test                    <-read.table("~/igv/polyA/polyA/clusters_withTissueInfo_mm10torn5.sorted.bed",sep="\t")[,c(1,2,3,6)]
colnames(test)          <-c("seqnames","start","end","strand")
polya.unibas.musmusc    <- makeGRangesFromDataFrame(test,ignore.strand=FALSE,seqnames.field="seqnames",start.field="start",end.field="end",strand.field="strand")
test                    <-read.table("~/igv/polyA/polyA/clusters_withTissueInfo_h19torn5.sorted.bed",sep="\t")[,c(1,2,3,6)]
colnames(test)          <-c("seqnames","start","end","strand")
polya.unibas.h19    <- makeGRangesFromDataFrame(test,ignore.strand=FALSE,seqnames.field="seqnames",start.field="start",end.field="end",strand.field="strand")
polya.sites             <- resize(c(polya.unibas.h19,polya.unibas.musmusc),width=5,fix="end")


#9. Derti sites from rat (testis and brain)
brain.sites  <- import.bed("~/Desktop/DataAnalsyisRiccio/Derti/GSM747486_rat_brain.sites.clustered.rn5.sorted.bed")
brain.algn   <- import.bed("~/Desktop/DataAnalsyisRiccio/Derti/GSM747486_rat_brain.alignments.sum.rn5.sorted.bed")
testis.sites <- import.bed("~/Desktop/DataAnalsyisRiccio/Derti/GSM747487_rat_testis.sites.clustered.rn5.sorted.bed")
testis.algn  <- import.bed("~/Desktop/DataAnalsyisRiccio/Derti/GSM747487_rat_testis.alignments.sum.rn5.sorted.bed")
derti.sites <- c(testis.sites,brain.sites)

polyA.sites <- MergeGRangesObject(polya.sites,derti.sites)

polyA.atlas <- MergeGRangesObject(polyA.sites,genome.annotation)

pdf("~/Desktop/DataAnalsyisRiccio/Rebuttal_May_2017/validation_polyadenylation_atlas.pdf")
par(mfrow=c(3,3))
PlotFractionValidated(siteoi=polya_db,name="polyAdb")
PlotFractionValidated(siteoi=apadb,name="APADB")
PlotFractionValidated(siteoi=MergeGRangesObject(polya_db,apadb),name="APADB-polyA_dB")
PlotFractionValidated(siteoi=refseq,name="RefSeq (Rn5-Rn6)")
PlotFractionValidated(siteoi=Extract_utr(rn6),name="Ensembl Rn6")
PlotFractionValidated(siteoi=Extract_utr(xeno.musmusc),name="Xeno RefSeq mm10")
PlotFractionValidated(siteoi=Extract_utr(xeno.h19),name="Xeno RefSeq h19")
PlotFractionValidated(siteoi=Extract_utr(xeno.all),name="Xeno RefSeq all")
PlotFractionValidated(siteoi=genome.annotation,name="(xeno)RefSeq -- Ensembl")
PlotFractionValidated(siteoi=polya.sites,name="polyA-unibasel")
PlotFractionValidated(siteoi=derti.sites,name="Derti data")
PlotFractionValidated(siteoi=polyA.sites,name="Derti +unibasel")
PlotFractionValidated(siteoi=polyA.atlas,name="polyA atlas")
par(mfrow=c(1,3))
PlotFractionValidated(siteoi=polyA.atlas,name="polyA atlas")
plot(density(log10(1+data.frame(distanceToNearest(x=SelectSurroundingEnd(myGR=myUTR.new,before=1,after=0),subject=polyA.atlas,ignore.strand=FALSE))[,3])),frame=FALSE)
plot(ecdf(log10(1+data.frame(distanceToNearest(x=SelectSurroundingEnd(myGR=myUTR.new,before=1,after=0),subject=polyA.atlas,ignore.strand=FALSE))[,3])),frame=FALSE)
dev.off()


100*GetFractionOverlap(utrGR=SelectSurroundingEnd(myGR=myUTR.new,before=150,after=150),target=polyA.atlas)
100*GetFractionOverlap(utrGR=SelectSurroundingEnd(myGR=myUTR.new,before=100,after=50),target=polyA.atlas)


myfun<-ecdf(log10(1+data.frame(distanceToNearest(x=SelectSurroundingEnd(myGR=myUTR.new,before=1,after=0),subject=polyA.atlas,ignore.strand=FALSE))[,3]))

AnalysisFractionOverlap(MergeGRangesObject(c(testis.sites,brain.sites),c(polya.unibas.h19,polya.unibas.musmusc)))#30.51609 51.68505 65.57352 73.76455 79.45821
#Derti+polyAunibasel+polyA_db (does not change anything)
AnalysisFractionOverlap(MergeGRangesObject(polya_db,MergeGRangesObject(c(testis.sites,brain.sites),c(polya.unibas.h19,polya.unibas.musmusc))))#30.66344 51.85129 65.78510 73.93834 79.67357
#Derti+xenoRefSeq (does not change anything)
AnalysisFractionOverlap(MergeGRangesObject(Extract_utr(xeno.all),MergeGRangesObject(polya_db,MergeGRangesObject(c(testis.sites,brain.sites),c(polya.unibas.h19,polya.unibas.musmusc)))))# 33.21369 55.13072 69.02675 76.99864 82.34094 --plot this one as add in Supplementary material the individual comparisons
#Derti+xenoRefSeqH19 (does not change anything)
AnalysisFractionOverlap(MergeGRangesObject(Extract_utr(xeno.h19),MergeGRangesObject(polya_db,MergeGRangesObject(c(testis.sites,brain.sites),c(polya.unibas.h19,polya.unibas.musmusc)))))#31.40396 52.99229 66.97522 75.07934 80.70878
#Derti+xenoRefSeqmm9
AnalysisFractionOverlap(MergeGRangesObject(Extract_utr(xeno.musmusc),MergeGRangesObject(polya_db,MergeGRangesObject(c(testis.sites,brain.sites),c(polya.unibas.h19,polya.unibas.musmusc)))))#31.57020 52.79961 66.74853 74.89043 80.49343

AnalysisFractionOverlap(MergeGRangesObject(Extract_utr(refseq.rn5),MergeGRangesObject(Extract_utr(xeno.all),MergeGRangesObject(polya_db,MergeGRangesObject(c(testis.sites,brain.sites),c(polya.unibas.h19,polya.unibas.musmusc))))))# 33.34215 55.34230 69.28366 77.24044 82.59030

#All sources (Ensembl Rn6, RefSeq Rn5, PolyA_db, Derti, poly_unibasel)
AnalysisFractionOverlap(MergeGRangesObject(Extract_utr(rn6),MergeGRangesObject(Extract_utr(refseq.rn5),MergeGRangesObject(Extract_utr(xeno.all),MergeGRangesObject(polya_db,MergeGRangesObject(c(testis.sites,brain.sites),c(polya.unibas.h19,polya.unibas.musmusc)))))))#33.81442 55.83724 69.79749 77.75805 83.02101

#All sources (Ensembl Rn6, RefSeq Rn5-Rn6, PolyA_db, Derti, poly_unibasel)
AnalysisFractionOverlap(MergeGRangesObject(Extract_utr(refseq.rn6),MergeGRangesObject(Extract_utr(rn6),MergeGRangesObject(Extract_utr(refseq.rn5),MergeGRangesObject(Extract_utr(xeno.all),MergeGRangesObject(polya_db,MergeGRangesObject(c(testis.sites,brain.sites),c(polya.unibas.h19,polya.unibas.musmusc))))))))#33.81820 55.84102 69.80505 77.76560 83.02478

xeno.all.p                <- import.gff("~/Downloads/XenoRefGene_ucs.gtf",format="gtf")#all xeno
AnalysisFractionOverlap(MergeGRangesObject(Extract_utr(refseq.rn6),MergeGRangesObject(Extract_utr(rn6),MergeGRangesObject(Extract_utr(refseq.rn5),MergeGRangesObject(Extract_utr(xeno.all.p),MergeGRangesObject(polya_db,MergeGRangesObject(c(testis.sites,brain.sites),c(polya.unibas.h19,polya.unibas.musmusc))))))))#33.84086 55.85613 69.82016 77.78072 83.05123

export.gff(xeno.all.p,"~/igv/Rn5/Xeno_refseq_ucsc.gtf",format="gtf")

#All together (XenoRefSeq, PolyA_unibasel,polya_db )
AnalysisFractionOverlap(MergeGRangesObject(gr1=polya_db,gr2=MergeGRangesObject(gr1=c(testis.sites,brain.sites),gr2=MergeGRangesObject(gr2=Extract_utr(xeno=xeno.all),gr1=c(polya.unibas.h19,polya.unibas.musmusc)))))#33.21369 55.13072 69.02675 76.99864 82.34094

#All together (XenoRefSeq, PolyA_unibasel,polya_db,apadb)
AnalysisFractionOverlap(MergeGRangesObject(apadb,MergeGRangesObject(gr1=polya_db,gr2=MergeGRangesObject(gr1=c(testis.sites,brain.sites),gr2=MergeGRangesObject(gr2=Extract_utr(xeno=xeno.all),gr1=c(polya.unibas.h19,polya.unibas.musmusc))))))#34.08267 56.18861 69.79371 77.67871 82.86610

all.gr <- MergeGRangesObject(apadb,MergeGRangesObject(gr1=polya_db,gr2=MergeGRangesObject(gr1=c(testis.sites,brain.sites),gr2=MergeGRangesObject(gr2=Extract_utr(xeno=xeno.all),gr1=c(polya.unibas.h19,polya.unibas.musmusc)))))

plot(density(log10(1+data.frame(distanceToNearest(x=SelectSurroundingEnd(myGR=myUTR.new,before=1,after=0),subject=all.gr,ignore.strand=FALSE))[,3])),frame=FALSE)
plot(ecdf(log10(1+data.frame(distanceToNearest(x=SelectSurroundingEnd(myGR=myUTR.new,before=1,after=0),subject=all.gr,ignore.strand=FALSE))[,3])),frame=FALSE)
barplot(AnalysisFractionOverlap(all.gr))

GetFractionOverlap(utrGR=SelectSurroundingEnd(myGR=myUTR.new,before=2,after=2),target=all.gr)

#Comparison with known polyA sites
#Coordinates of known transcript termini (Fig. 3B) were compiled from RefSeq Rn5 and Rn6, and Ensembl (Rn6) transcript models (downloaded from UCSC) (Fujita et al. 2011), taking the 3'-most site when transcript isoforms overlapped on the same genomic strand. 
#
#Known 3' UTR coordinates (Fig. 3D,E) were compiled from RefSeq (Rn5 and Rn6) Ensembl (Rn6) transcripts by collapsing overlapping isoforms (on the same strand) into a single 3' UTR model that represented the union of genomic coordinates.

#Previously reported polyA sites (Fig. 3G) were compiled from PolyA-DB2 (Lee et al. 2007), dbEST (Boguski et al. 1993), GenBank mRNAs (Wheeler et al. 2003), RefSeq genes (Wheeler et al. 2003), and Ensembl genes (Flicek et al. 2011). Coordinates were obtained from UCSC (Fujita et al. 2011) when available, or by aligning against the respective genomes using BLAT (Kent 2002) (using default settings and taking the best alignment).

# A. Validation
# A.1 RefSeq (Rn5, Rn6)
# A.2 PolyA_db
# A.3 Ensembl Rn6
# A.4 Xeno mus musculus
# A.5 Xeno homo sapiens
# A.6 Xeno 

# B. Fraction recovered annotation
# B.1 RefSeq
# B.2 PolyA_db
# B.3 Rn6

# Assembling the Polyadenylation Site Atlas. 
# First 3' UTR annotations from RefSeq (Rn5 and Rn6) and Ensembl (Rn6) are combined, where overlapping regions are merged. Given that rat genome is poorly annotated compared to other species such as human and mouse, and because 3' end are highly conserved across species, we futher combined existing rat 3' UTR annotation with RefSeq annotation from other species (xenoRefseq). 

# Then, to generate a comprehensive atlas of rat PAS in 3' UTR, polyadenylation annotations from PolyA_DB 2 (37), and APADB (39) were used together with the polyadyenlation atlas from unibasel (both human and mouse). Finally PAS obtained from rat 3'-Seq data (Derti et al.) were also used to expand the repertoire of PAS.
# 
# Genome annotation files from Rn6 assembly or other species have been converted to Rn5 using Crossmap.  
# 




myu1           <- resize(myUTR.1,width=2,fix="end")
myu2           <- resize(myUTR.2,width=2,fix="end")
myu3           <- resize(myUTR.3,width=2,fix="end")

#Identification of problematic once
EnsembltoNew <- data.frame(distanceToNearest(x=myu2,subject=myu1,ignore.strand=FALSE))
myu2[EnsembltoNew[EnsembltoNew[,3]==0,1],]



#Compute fraction recovered in function of distance
GetFracRecov <- function(ct=myu1,denov=myu2,target=subrefseq.p,w=1){
  ct                 <- resize(x=ct,fix="start",width=w,use.names=T)
  ct                 <- resize(x=ct,fix="end",width=w*2,use.names=T)
  denov              <- resize(x=denov,fix="start",width=w,use.names=T)
  denov              <- resize(x=denov,fix="end",width=w*2,use.names=T)
  
  in_Ens            <- unique(subjectHits(findOverlaps(query=ct,subject=target,ignore.strand=F)))
  in_New            <- unique(subjectHits(findOverlaps(query=denov,subject=target,ignore.strand=F)))
  ix.all            <- c(1:length(target))
  
  #ntot                <- length(target)
  n_not_in_Ens        <- length(setdiff(ix.all,in_Ens))
  #n_not_in_any        <- length(setdiff(ix.all,unique(c(in_Ens,in_New))))
  n_in_denov_not_ens  <- length(setdiff(in_New,in_Ens))
  frac.improved       <- n_in_denov_not_ens/n_not_in_Ens
  
  return(round(frac.improved,digit=2))
}


NEAT1_S chr1:228133821-228137077 -
  NEAT1_L1 chr1:228132548-228133609 -
  NEAT1_L2 chr1:228129563-228130675 -
  
  start        <- c(228133821,228132548,228129563)
end          <- c(228137077,228133609,228130675)
chr          <- c("chr1","chr1","chr1")
strand       <- c("-","-","-")
name         <- c("NEAT1_S","NEAT1_L1","NEAT1_L2")
out          <- cbind(chr, start,end,name,score=rep(1,length(start)),strand)
write.table(x=out,file="~/Desktop/NEAT1.bed", append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", row.names = FALSE,col.names = FALSE)
bedtools getfasta -fi  /Users/raphaelle/igv/genomes/rn5.fa -bed ~/Desktop/NEAT1.bed -s -name -fo ~/Desktop/NEAT1_Rn.fasta



require("rtracklayer")
malat1 <- import.gff("~/Desktop/MALAT1_Rn.gtf")
start        <- as.numeric(start(malat1))
end          <- as.numeric(end(malat1))
chr          <- as.character(seqnames(malat1))
strand       <- as.character(strand(malat1))
name         <- mcols(malat1)$exon_id
out          <- cbind(chr, start,end,name,score=rep(1,length(start)),strand)
write.table(x=out,file="~/Desktop/MALAT1_Rn.bed", append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", row.names = FALSE,col.names = FALSE)
bedtools getfasta -fi  /Users/raphaelle/igv/genomes/rn5.fa -bed ~/Desktop/MALAT1_Rn.bed -s -name -fo ~/Desktop/MALAT1_Rn.fasta



NEAT1 chr1:228133805-228137098 - reads in axonal compartments
Malat1_mm9 chr1:228084241-228091179 no reads in axonal compartments
MALAT1_h19 chr1:228084316-228090094 reads in axonal compartments
Pex3 chr1:9313349-9354376
Odz4 chr1:168082724-168094331
Btaf1 chr1:263055333-263130789
Exosc3 chr5:65375285-65380475 -
  Egr3 chr15:55482586-55487962 +
  Atpaf2 chr10:46523882-46539221 -
  
  Tug1_whole_locus chr14:84563853-84570887 -  
  Tug1_whole_utr chr14:84563866-84566455 -
  Xist chrX:75129831-75147353 -
  
  Xist_exon_1  
start        <- c(84563853,84563866,75129831)
end          <- c(84570887,84566455,75147353)
chr          <- c("chr14","chr14","chrX")
strand       <- c("-","-","-")
name         <- c("Tug1_whole_locus","Tug1_whole_utr","Xist")
out          <- cbind(chr, start,end,name,score=rep(1,length(start)),strand)[-7,]
write.table(x=out,file="~/Desktop/DataAnalsyisRiccio/Dec2016/axonal_remodelling/xist_Tug1.bed", append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", row.names = FALSE,col.names = FALSE)

bedtools getfasta -fi  /Users/raphaelle/igv/genomes/rn5.fa -bed ~/Desktop/DataAnalsyisRiccio/Dec2016/axonal_remodelling/xist_Tug1.bed -s -name -fo ~/Desktop/DataAnalsyisRiccio/Dec2016/axonal_remodelling/xist_Tug1.fasta

chrX:75138032-75138182
chrX:75129846-75137192

start        <- c(75138032,75129846)
end          <- c(75138182,75137192)
chr          <- c("chrX","chrX")
strand       <- c("-","-")
name         <- c("Xist_E1","Xist_E2")
out          <- cbind(chr, start,end,name,score=rep(1,length(start)),strand)[-7,]
write.table(x=out,file="~/Desktop/DataAnalsyisRiccio/Dec2016/axonal_remodelling/xist.bed", append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", row.names = FALSE,col.names = FALSE)
bedtools getfasta -fi  /Users/raphaelle/igv/genomes/rn5.fa -bed ~/Desktop/DataAnalsyisRiccio/Dec2016/axonal_remodelling/xist.bed -s -name -fo ~/Desktop/DataAnalsyisRiccio/Dec2016/axonal_remodelling/xist.fa


start        <- c(228133805,228084241,228084316,9313349,168082724,263055333,272540500,65375285,55482586,46523882)
end          <- c(228137098,228091179,228090094,9354376,168094331,263130789,272705421,65380475,55487962,46539221)
chr          <- c("chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr5","chr15","chr10")
strand       <- c("-","-","-","-","+","+","+","-","+","-")
name         <- c("NEAT1_mm9","Malat1_mm9","MALAT1_h19","Pex3","Odz4","Btaf1","Btrc","Exosc3","Egr3","Atpaf2")
out          <- cbind(chr, start,end,name,score=rep(1,length(start)),strand)[-7,]

write.table(x=out,file="~/Desktop/DataAnalsyisRiccio/Dec2016/axonal_remodelling/target_fraction_whole_locus.bed", append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", row.names = FALSE,col.names = FALSE)


start        <- c(228133805,228084241,228084316,9313349,168093701,263055333,272540500,65375285,55485494,46523882)
end          <- c(228135577,228085376,228085376,9314069,168094331,263127957,272705421,65375647,55487962,46524504)
chr          <- c("chr1","chr1","chr1","chr1","chr1","chr1","chr1","chr5","chr15","chr10")
strand       <- c("-","-","-","-","+","+","+","-","+","-")
name         <- c("NEAT1_mm9","Malat1_mm9","MALAT1_h19","Pex3","Odz4","Btaf1","Btrc","Exosc3","Egr3","Atpaf2")
out          <- cbind(chr, start,end,name,score=rep(1,length(start)),strand)[-7,]

write.table(x=out,file="~/Desktop/DataAnalsyisRiccio/Dec2016/axonal_remodelling/target_fraction_3utr.bed", append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", row.names = FALSE,col.names = FALSE)


Impa1 chr2 113451413-113451534 +
  Atp5f1 chr2 227828728-227829082 -
  
  start        <- c(113451413,227828728)
end          <- c(113451534,227829082)
chr          <- c("chr2","chr2")
strand       <- c("+","-")
name         <- c("Impa1","Atp5f1")
out          <- cbind(chr, start,end,name,score=rep(1,length(start)),strand)[-7,]
write.table(x=out,file="~/Desktop/DataAnalsyisRiccio/Dec2016/axonal_remodelling/impa1_atp5f1.bed", append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", row.names = FALSE,col.names = FALSE)


#In BASH
bedtools getfasta -fi  /Users/raphaelle/igv/genomes/rn5.fa -bed ~/Desktop/DataAnalsyisRiccio/Dec2016/axonal_remodelling/target_fraction_whole_locus.bed -s -name -fo ~/Desktop/DataAnalsyisRiccio/Dec2016/axonal_remodelling/target_fraction_whole_locus.fasta
bedtools getfasta -fi  /Users/raphaelle/igv/genomes/rn5.fa -bed ~/Desktop/DataAnalsyisRiccio/Dec2016/axonal_remodelling/target_fraction_3utr.bed -s -name -fo ~/Desktop/DataAnalsyisRiccio/Dec2016/axonal_remodelling/target_fraction_3utr.fasta

bedtools getfasta -fi  /Users/raphaelle/igv/genomes/rn5.fa -bed ~/Desktop/DataAnalsyisRiccio/Dec2016/axonal_remodelling/impa1_atp5f1.bed -s -name -fo ~/Desktop/DataAnalsyisRiccio/Dec2016/axonal_remodelling/impa1_atp5f1.fasta


#First compute distance to closest one in order to label whether or not in RefSeq
EnsembltoRefseq <- data.frame(distanceToNearest(x=subrefseq.p,subject=myu1,ignore.strand=FALSE))
newtoRefSeq     <- data.frame(distanceToNearest(x=subrefseq.p,subject=myu2,ignore.strand=FALSE))
newtoRefSeq     <- newtoRefSeq[match(EnsembltoRefseq[,1],newtoRefSeq[,1]),]

ix1             <- newtoRefSeq[,3]==0|EnsembltoRefseq[,3]==0

pdf("~/Desktop/DataAnalysisRiccio/Dec2016/utrid/anno_ngf_stringent/comparison_with_refseq.pdf")
par(mfrow=c(2,3),mar=c(4,4,3,3))
multidensity(list(Ensembl=log10(EnsembltoRefseq[!ix1,3]),New=log10(newtoRefSeq[!ix1,3])),frame=F,main=paste("#RefSeq not in Ensembl=",sum(!ix1),sep=""),xlab="distance t to RefSeq")
mp<-multiecdf(list(Ensemb=log10(EnsembltoRefseq[!ix1,3]),New=log10(newtoRefSeq[!ix1,3])),frame=F,xlab="distance t to RefSeq")
abline(v=c(1,2,3,4),col="grey",lty=2)
x        <- seq(0,4,by=0.05)
improved <- mp[[2]](x)-mp[[1]](x)
plot(x,improved,type="p",pch=19,main="increase fraction of New distant t to Refseq",frame=F)
plot(log10(EnsembltoRefseq[!ix1,3]),log10(newtoRefSeq[!ix1,3]),frame=F,cex=0.4,pch=19,las=1,main=paste("# dNew<dEnsembl=",sum(log10(EnsembltoRefseq[!ix1,3])-log10(newtoRefSeq[!ix1,3])>0),sep=""))
mtext(side=3,line=0,text=paste(round(100*sum(log10(EnsembltoRefseq[!ix1,3])-log10(newtoRefSeq[!ix1,3])>0)/sum(!ix1),digits=2),"% improved",sep=""))
tempval <- numeric()
for(i in c(1:10^4)){
  ix1             <- newtoRefSeq[,3]>=i
  tempval         <- c(tempval,round(100*sum(log10(EnsembltoRefseq[ix1,3])-log10(newtoRefSeq[ix1,3])>0)/sum(ix1)))
}
plot(log10(c(1:10^4)),tempval,ylim=c(0,60),las=1,type="l",frame=F,ylab="% improved",xlab="distance Ensembl to RefSeq")
dev.off()


#Then compute distance to closest Xeno one in order to label whether or not in Xeno
EnsembltoRefseq <- data.frame(distanceToNearest(x=subxeno.exp,subject=myu1,ignore.strand=FALSE))
newtoRefSeq     <- data.frame(distanceToNearest(x=subxeno.exp,subject=myu2,ignore.strand=FALSE))
newtoRefSeq     <- newtoRefSeq[match(EnsembltoRefseq[,1],newtoRefSeq[,1]),]

ix1             <- newtoRefSeq[,3]==0

pdf("~/Desktop/DataAnalysisRiccio/Dec2016/utrid/anno_ngf_stringent/comparison_with_xeno.pdf")
par(mfrow=c(2,3),mar=c(4,4,3,3))
multidensity(list(Ensembl=log10(EnsembltoRefseq[!ix1,3]),New=log10(newtoRefSeq[!ix1,3])),frame=F,main=paste("#XenoRefseq not in Ensembl=",sum(!ix1),sep=""),xlab="distance t to XenoRefseq")
mp<-multiecdf(list(Ensemb=log10(EnsembltoRefseq[!ix1,3]),New=log10(newtoRefSeq[!ix1,3])),frame=F,xlab="distance t to XenoRefseq")
abline(v=c(1,2,3,4),col="grey",lty=2)
x        <- seq(0,4,by=0.05)
improved <- mp[[2]](x)-mp[[1]](x)
plot(x,improved,type="p",pch=19,main="increase fraction of New distant t to XenoRefseq",frame=F)
plot(log10(EnsembltoRefseq[!ix1,3]),log10(newtoRefSeq[!ix1,3]),frame=F,cex=0.4,pch=19,las=1,main=paste("# dNew<dEnsembl=",sum(log10(EnsembltoRefseq[!ix1,3])-log10(newtoRefSeq[!ix1,3])>0),sep=""))
mtext(side=3,line=0,text=paste(round(100*sum(log10(EnsembltoRefseq[!ix1,3])-log10(newtoRefSeq[!ix1,3])>0)/sum(!ix1),digits=2),"% improved",sep=""))
tempval <- numeric()
for(i in c(1:10^4)){
  ix1             <- newtoRefSeq[,3]>=i
  tempval         <- c(tempval,round(100*sum(log10(EnsembltoRefseq[ix1,3])-log10(newtoRefSeq[ix1,3])>0)/sum(ix1)))
}
plot(log10(c(1:10^4)),tempval,ylim=c(0,60),las=1,type="l",frame=F,ylab="% improved",xlab="distance Ensembl to XenoRefseq")
dev.off()

#Then compute distance to closest polyA.db one in order to label whether or not in Xeno
padb <- import.gff("~/igv/Riccio/NewExp1/paDB.sorted.gtf",format="gtf")
EnsembltoRefseq <- data.frame(distanceToNearest(x=padb,subject=myu1,ignore.strand=FALSE))
newtoRefSeq     <- data.frame(distanceToNearest(x=padb,subject=myu2,ignore.strand=FALSE))
newtoRefSeq     <- newtoRefSeq[match(EnsembltoRefseq[,1],newtoRefSeq[,1]),]
ix1             <- newtoRefSeq[,3]==0

pdf("~/Desktop/DataAnalysisRiccio/Dec2016/utrid/anno_ngf_stringent/comparison_with_polyAdb.pdf")
par(mfrow=c(2,3),mar=c(4,4,3,3))
multidensity(list(Ensembl=log10(EnsembltoRefseq[!ix1,3]),New=log10(newtoRefSeq[!ix1,3])),frame=F,main=paste("#PolyA.db not in Ensembl=",sum(!ix1),sep=""),xlab="distance t to polyA.db")
mp<-multiecdf(list(Ensemb=log10(EnsembltoRefseq[!ix1,3]),New=log10(newtoRefSeq[!ix1,3])),frame=F,xlab="distance t to PolyA.db")
abline(v=c(1,2,3,4),col="grey",lty=2)
x        <- seq(0,4,by=0.05)
improved <- mp[[2]](x)-mp[[1]](x)
plot(x,improved,type="p",pch=19,main="increase fraction of New distant t to PolyA.db",frame=F)
plot(log10(EnsembltoRefseq[!ix1,3]),log10(newtoRefSeq[!ix1,3]),frame=F,cex=0.4,pch=19,las=1,main=paste("# dNew<dEnsembl=",sum(log10(EnsembltoRefseq[!ix1,3])-log10(newtoRefSeq[!ix1,3])>0),sep=""))
mtext(side=3,line=0,text=paste(round(100*sum(log10(EnsembltoRefseq[!ix1,3])-log10(newtoRefSeq[!ix1,3])>0)/sum(!ix1),digits=2),"% improved",sep=""))
tempval <- numeric()
for(i in c(1:10^4)){
  ix1             <- newtoRefSeq[,3]>=i
  tempval         <- c(tempval,round(100*sum(log10(EnsembltoRefseq[ix1,3])-log10(newtoRefSeq[ix1,3])>0)/sum(ix1)))
}
plot(log10(c(1:10^4)),tempval,ylim=c(0,60),las=1,type="l",frame=F,ylab="% improved",xlab="distance Ensembl to PolyA.db")
dev.off()
























reoi<- subrefseq.p[which(EnsembltoRefseq[,3]>20&newtoRefSeq[,3]<1),]
myDistToRefSeq     <- list(EnsembltoRefseq=log10(1+data.frame(distanceToNearest(x=myu1,subject=subrefseq,ignore.strand=FALSE))$distance),
                           newtoRefSeq=log10(1+data.frame(distanceToNearest(x=myu2,subject=subrefseq,ignore.strand=FALSE))$distance))

myDistToXeno        <- list(EnsembltoRefseq=log10(1+data.frame(distanceToNearest(x=myu1,subject=subxeno,ignore.strand=FALSE))$distance),
                            newtoRefSeq=log10(1+data.frame(distanceToNearest(x=myu2,subject=subxeno,ignore.strand=FALSE))$distance))



layout(matrix(c(1,1,2,3,2,3), 2, 3, byrow = FALSE))
plot(density(log10(newtoRefSeq[,3])),col="red",frame=F,las=1)
lines(density(log10(EnsembltoRefseq[,3])),col="grey")
abline(v=log10(100),col="black",lty=2)

plot(ecdf(EnsembltoRefseq[EnsembltoRefseq[,3]>100,3]),xlim=c(0,1000),frame=F,lwd=1,las=1,xlab="distance to nearest Ensembl 3' UTR",ylab="% RefSeq >= x",main=paste("n=",sum(EnsembltoRefseq[,3]>100), " RefSeq txID with d>100 to Ensembl"))
temp            <- ecdf(EnsembltoRefseq[EnsembltoRefseq[,3]>100,3])
abline(v=c(1:1000)[which(round(temp(c(1:1000)),digit=1)==0.5)[1]],col="red",lty=2)
abline(h=0.5,col="red",lty=2)
text(x=c(1:1000)[which(round(temp(c(1:1000)),digit=1)==0.5)[1]],y=0.6,lab=c(1:1000)[which(round(temp(c(1:1000)),digit=1)==0.5)[1]],col="red")


plot(ecdf(newtoRefSeq[EnsembltoRefseq[,3]>100,3]),xlim=c(0,1000),frame=F,lwd=1,las=1,xlab="distance to nearest novel 3' UTR",ylab="% RefSeq >= x",main=paste("n=",sum(EnsembltoRefseq[,3]>100), " RefSeq txID with d>100 to Ensembl"))
temp            <- ecdf(newtoRefSeq[EnsembltoRefseq[,3]>100,3])
abline(v=100,col="red",lty=2)
abline(h=temp(100),col="red",lty=2)
text(x=105,y=temp(100)+0.1,lab=round(temp(100),digits=2),col="red")
text(x=10,y=temp(10)+0.1,lab=round(temp(10),digits=2),col="green")
abline(v=10,col="green",lty=2)
abline(h=temp(10),col="green",lty=2)


layout(matrix(c(1,1,2,3,2,3), 2, 3, byrow = FALSE))
plot(density(log10(newtoRefSeq[,3])),col="red",frame=F,las=1)
lines(density(log10(EnsembltoRefseq[,3])),col="grey")
abline(v=log10(1000),col="black",lty=2)

plot(ecdf(EnsembltoRefseq[EnsembltoRefseq[,3]>1000,3]),xlim=c(0,10000),frame=F,lwd=1,las=1,xlab="distance to nearest Ensembl 3' UTR",ylab="% RefSeq >= x",main=paste("n=",sum(EnsembltoRefseq[,3]>1000), " RefSeq txID with d>100 to Ensembl"))
temp            <- ecdf(EnsembltoRefseq[EnsembltoRefseq[,3]>1000,3])
abline(v=c(1:10000)[temp(c(1:10000))==0.5],col="red",lty=2)
abline(h=0.5,col="red",lty=2)
text(x=c(1:10000)[temp(c(1:10000))==0.5],y=0.6,lab=c(1:10000)[temp(c(1:10000))==0.5],col="red")

plot(ecdf(newtoRefSeq[EnsembltoRefseq[,3]>1000,3]),xlim=c(0,1000),frame=F,lwd=1,las=1,xlab="distance to nearest novel 3' UTR",ylab="% RefSeq >= x",main=paste("n=",sum(EnsembltoRefseq[,3]>1000), " RefSeq txID with d>1000 to Ensembl"))
temp            <- ecdf(newtoRefSeq[EnsembltoRefseq[,3]>1000,3])
abline(v=100,col="red",lty=2)
abline(h=temp(100),col="red",lty=2)
text(x=105,y=temp(100)+0.1,lab=round(temp(100),digits=2),col="red")
text(x=10,y=temp(10)+0.1,lab=round(temp(10),digits=2),col="green")
abline(v=10,col="green",lty=2)
abline(h=temp(10),col="green",lty=2)



mycol <- rep("grey",nrow(newtoRefSeq))
mycol[EnsembltoRefseq[,3]>20]<-"red"
plot(log10(1+EnsembltoRefseq[,3]),log10(1+newtoRefSeq[,3]),col=mycol,pch=19,cex=0.6)
abline(v=log10(50))
abline(h=log10(50))

hist(log10(EnsembltoRefseq[EnsembltoRefseq[,3]>20,3]+1),col=rgb(0,0,0),breaks=75)
hist(log10(newtoRefSeq[EnsembltoRefseq[,3]>20,3]+1),col=rgb(1,0,0,0.2),add=T,breaks=100)

multidensity(list(refseqTOensemb=log10(EnsembltoRefseq[,3]),refseTOnew=log10(newtoRefSeq[,3])))

n1 <- sum(EnsembltoRefseq[,3]>20&newtoRefSeq[,3]<10)
n2 <- sum(EnsembltoRefseq[,3]>20)
n3 <- sum(EnsembltoRefseq[,3]>20&newtoRefSeq[,3]>10)


reoi<- subrefseq[which(EnsembltoRefseq[,3]>10&newtoRefSeq[,3]<10),]
temp<-ecdf(newtoRefSeq[EnsembltoRefseq[,3]>10,3]),xlim=c(0,500)
n1/(n1+n2)
16% of RefSeq absent of Ensembl AND expressed in these samples are recovered by the new pipeline-->I should check how many are indeed expressed


n1 <- sum(EnsembltoRefseq[,3]>10&newtoRefSeq[,3]<10)
n2 <- sum(EnsembltoRefseq[,3]>10)
n3 <- sum(EnsembltoRefseq[,3]>10&newtoRefSeq[,3]>10)

par(mfrow=c(2,1))
multidensity(myDistToRefSeq,frame=F,xlab="distance to nearest [log10(bp)]",ylab="density")
multidensity(myDistToXeno,frame=F,xlab="distance to nearest [log10(bp)]",ylab="density")


sum(myDistToRefSeq[[1]]>log10(30))
8107
sum(myDistToRefSeq[[1]]<=log10(30))
14679
abline(v=log10(30))




