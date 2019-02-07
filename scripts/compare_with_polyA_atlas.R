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
names(myUTR) <- as.character(myUTR$ID)
ix1          <- match(names(myUTR),anno_ngf$uniqueID) 
myUTR        <- myUTR[!is.na(ix1),]
ix1          <- match(names(myUTR),anno_ngf$uniqueID) 
anno_ngf     <- anno_ngf[ix1,]
myUTR.new    <- myUTR[-grep(myUTR$ID,pattern="\\.0"),]

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


#1. PolyA_db
polya_db                <- import.gff("./annotation/rn5/rn.gtf")
AnalysisFractionRecovered(polya_db)#What is plotted on Figure 1D -- fraction of polyA absent from Rn5 and recovered

#2. APADB
test                    <-read.table("./annotation/rn5/hg19.apadb_v2_final_rn5.bed",sep="\t")[,c(1,2,3,6)]
colnames(test)          <-c("seqnames","start","end","strand")
apadb.h19              <- makeGRangesFromDataFrame(test,ignore.strand=FALSE,seqnames.field="seqnames",start.field="start",end.field="end",strand.field="strand")
test                    <-read.table("./annotation/rn5/mm10.apadb_v2_final_rn5.bed",sep="\t")[,c(1,2,3,6)]
colnames(test)          <-c("seqnames","start","end","strand")
apadb.mm10              <- makeGRangesFromDataFrame(test,ignore.strand=FALSE,seqnames.field="seqnames",start.field="start",end.field="end",strand.field="strand")
apadb <- c(apadb.h19,apadb.mm10)


#3. RefSeq
refseq.rn5              <- import.gff("./annotation/rn5/RefSeq_rn5_UCSC_27062016.gtf",format="gtf")
refseq.rn6              <- import.gff("./annotation/rn5/refSeq_rn6.Rn5.sorted.gtf",format="gtf")
refseq                  <- MergeGRangesObject(Extract_utr(refseq.rn5),Extract_utr(refseq.rn6))

#4.Ensembl Rn6 -- file too large >100M for Github
rn6                     <- import.gff("./annotation/rn5/Rattus_norvegicus.Rnor_6.0.87.copy.sorted.Rn5.sorted.gtf")



#5. RefSeq xeno.musmusc
xeno.musmusc            <- import.gff("./annotation/rn5/xenoRefSeqrn5MusMusculus.gtf",format="gtf")
#6. RefSeq xeno.h19
xeno.h19                <- import.gff("./annotation/rn5/xenoRefSeqrn5HomoSapiens.gtf",format="gtf")
#7. RefSeq xeno.all -- file too large >100M for Github
xeno.all                <- import.gff("/annotation/rn5/xenoRefSeqrn5.sorted.gtf",format="gtf")#all xeno

genome.annotation       <- MergeGRangesObject(MergeGRangesObject(Extract_utr(xeno.all),Extract_utr(rn6)),refseq)


#8. PolyA_unibasel -- mm10
test                    <-read.table("./annotation/rn5/clusters_withTissueInfo_mm10torn5.sorted.bed",sep="\t")[,c(1,2,3,6)]
colnames(test)          <-c("seqnames","start","end","strand")
polya.unibas.musmusc    <- makeGRangesFromDataFrame(test,ignore.strand=FALSE,seqnames.field="seqnames",start.field="start",end.field="end",strand.field="strand")
test                    <-read.table("./annotation/rn5/clusters_withTissueInfo_h19torn5.sorted.bed",sep="\t")[,c(1,2,3,6)]
colnames(test)          <-c("seqnames","start","end","strand")
polya.unibas.h19    <- makeGRangesFromDataFrame(test,ignore.strand=FALSE,seqnames.field="seqnames",start.field="start",end.field="end",strand.field="strand")
polya.sites             <- resize(c(polya.unibas.h19,polya.unibas.musmusc),width=5,fix="end")


#9. Derti sites from rat (testis and brain)
brain.sites  <- import.bed("./annotation/rn5/GSM747486_rat_brain.sites.clustered.rn5.sorted.bed")
brain.algn   <- import.bed("./annotation/rn5/GSM747486_rat_brain.alignments.sum.rn5.sorted.bed")


testis.sites <- import.bed("./annotation/rn5/GSM747487_rat_testis.sites.clustered.rn5.sorted.bed")
testis.algn  <- import.bed("~/Desktop/DataAnalsyisRiccio/Derti/GSM747487_rat_testis.alignments.sum.rn5.sorted.bed")
derti.sites <- c(testis.sites,brain.sites)

polyA.sites <- MergeGRangesObject(polya.sites,derti.sites)

polyA.atlas <- MergeGRangesObject(polyA.sites,genome.annotation)

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
