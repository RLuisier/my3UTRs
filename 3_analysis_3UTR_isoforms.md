This file describes the downstream analysis of the 3' UTR expression in cell body and axonal compartment.

1. Create GTF file which focus on the last 500 nt fragment of 3' UTR

```R
library(grDevices)
library(Rsamtools)
library(IRanges)
library(GenomicRanges)
library(rtracklayer)

FocusEnd <- function(gr=L3,before=300,after=0){
    POS                <- as.character(strand(gr))=="+"
    grn                <- gr
    start(grn)[POS]  <- end(gr[POS])-before
    end(grn)[POS]    <- end(gr[POS])+after
    end(grn)[!POS]   <- start(gr[!POS])+before
    start(grn)[!POS] <- start(gr[!POS])-after
    start(grn)[start(grn)<0]<-1
    return(grn)
}

L2                 <- import.gff("./utrid/APA/L2.gtf",format="gtf")
L2_500             <- FocusEnd(gr=L2,before=500,after=0)
export.gff(L2_500,"./utrid/APA/L2_500.gtf",format="gtf")

```


2. Extract coverage of the last 500 nt fragment

```
OUTDIR=./Coverage/utrCov/500/
ANNOTATION=./utrid/APA/L2_500.gtf
SF=/SAN/luscombelab/Riccio/Exp_1/FASTQ/P2P3P4/scaling_factor.txt
FRAC=0.05
awk '{print $1}' ./annotation/rn5/chrom_oi.txt | while read CHR
do
    cat ./SAMPLES.txt  | while read SAMPLE
    do
        BAM=${BAMDIR}/SAMPLE.bam
        ./scripts/extractCovFirstSegmentPerChrImproved.R ${BAM} ${OUTDIR} ${CHR} ${SAMPLE} ${ANNOTATION} ${SF} ${FRAC}
    done
done

```

3. Collect coverage

```R
covdir               <- "./Coverage/utrCov/500/"
foi                  <- list.files(covdir)
names                <- gsub(foi,pattern=".txt",repl="")
types                <- do.call(lapply(names,function(x)return(length(unlist(strsplit(x,split="[\\_,\\.,-]"))))),what=c)
names_norms          <- names[types=="5"]
names_norms          <- do.call(lapply(names[types=="5"],function(x)return(unlist(strsplit(x,split="[\\_,\\.,-]")))),what=rbind)
names_norms          <- data.frame(foi[types=="5"],names_norms)
colnames(names_norms )<-c("id","treatment","origin","replicate","chr","type")
names_raw           <- names[intersect(which(types=="4"),grep(names,pattern=".chr"))]
names_raw           <- do.call(lapply(names[intersect(which(types=="4"),grep(names,pattern=".chr"))],function(x)return(unlist(strsplit(x,split="[\\_,\\.,-]")))),what=rbind)
names_raw           <- data.frame(foi[intersect(which(types=="4"),grep(names,pattern=".chr"))],names_raw)
colnames(names_raw )<-c("id","treatment","origin","replicate","chr")

samples                 <- factor(paste(names_raw$treatment,names_raw$origin,names_raw$replicate,sep="."))
covraw_500              <- lapply(c(1:length(unique(samples))),function(X)return(do.call(lapply(which(samples==unique(samples)[X]), function(Z)return(read.table(paste(covdir,names_raw$id[Z],sep="")))),what=rbind)))

covraw_500.sum          <- lapply(covraw_500,function(Z)return(tapply(Z[,2],INDEX=factor(as.character(Z[,1])),FUN=sum)))
covraw_500              <- covraw_500.sum

rown500                 <-  unique(do.call(lapply(covraw_500,function(x)return(names(x))),what=c))
myCovRaw500             <- matrix(0,nrow=length(rown500),ncol=length(covraw_500))
for(i in c(1:length(covraw_500))){
  myCovRaw500[,i]       <- covraw_500[[i]][match(rown500,names(covraw_500[[i]]))]
}
colnames(myCovRaw500)   <- unlist(lapply(unique(samples),function(x)return(paste(x,".raw",sep=""))))
rownames(myCovRaw500)   <- rown500


covnorm_500             <- lapply(c(1:length(unique(samples))),function(X)return(do.call(lapply(which(samples==unique(samples)[X]), function(Z)return(read.table(paste(covdir,names_norms$id[Z],sep="")))),what=rbind)))

covnorm_500.sum          <- lapply(covnorm_500,function(Z)return(tapply(Z[,2],INDEX=factor(as.character(Z[,1])),FUN=sum)))
covnorm_500              <- covraw_500.sum

rown500                 <-  unique(do.call(lapply(covnorm_500,function(x)return(names(x))),what=c))
myCovNorm500            <- matrix(0,nrow=length(rown500),ncol=length(covnorm_500))
for(i in c(1:length(covnorm_500))){
  myCovNorm500[,i]      <- covraw_500[[i]][match(rown500,names(covraw_500[[i]]))]
}
colnames(myCovNorm500)  <- unlist(lapply(unique(samples),function(x)return(paste(x,".norm",sep=""))))
rownames(myCovNorm500)  <- rown500

myCov500                <- myCovRaw500
myCov500                <- cbind(myCov500,myCovNorm500[match(rownames(myCov500),rownames(myCovNorm500)),])
myCov500.RpM            <- myCov500

save(list="myCov500",file= "./data/myCov500.RData")
```


4. Characterise 3'--5' end bias in coverage as described in [characterise_tx_read_density.md](./scripts/characterise_tx_read_density.md)



#********************************************************
# C. Create annotation file (RUD, newL, initL, phasing of the 3' UTR isoforms, distance to prior)
#
#********************************************************
code in   --> dependencies/import_coverage_stringent.R
stored in --> /home/rluisier/data/Riccio/Exp_1/Dec2016/Coverage/utrCov_stringent/coverage_detection.RData
stored in --> /home/rluisier/data/Riccio/Exp_1/Dec2016/Coverage/utrCov_stringent/myCov_March_12_2016_500.RData

#Sep2017
code in --> dependencies/import_coverage_stringent_Sep2017.R
stored in --> /home/rluisier/data/Riccio/Exp_1/Dec2016/utrid/APA_stringent/Ltot_Sep2017.gtf
             /home/rluisier/data/Riccio/Exp_1/Dec2016/utrid/APA_stringent/Ltotp_Sep2017.gtf
             /home/rluisier/data/Riccio/Exp_1/Dec2016/utrid/APA_stringent/coverage_detection_Sep2017.RData
code in --> dependencies/create_anno_tot_Sep2017.R
stored in --> /home/rluisier/data/Riccio/Exp_1/Dec2016/utrid/anno_tot_stringent/anno_tot_Sep2017.RData
             /home/rluisier/data/Riccio/Exp_1/Dec2016/utrid/APA_stringent/Ltot_sub_Sep2017.gtf
             ~/Desktop/DataAnalsyisRiccio/Dec2016/utrid/APA_stringent/anno_tot_stringent/anno_tot_corr_Sep2017.RData


#Nov2017
code in --> AxonalTransport/import_coverage_Nov2017.R
stored in --> ~/Desktop/DataAnalsyisRiccio/Nov2017/Ltot_Nov2017.gtf
          --> ~/Desktop/DataAnalsyisRiccio/Nov2017/Ltot_sub_Nov2017.gtf
          --> ~/Desktop/DataAnalsyisRiccio/Nov2017/anno_tot_Nov2017.RData


code in   --> dependencies/create_anno_ngf.R
stored in --> Coverage/utrCov/coverageUTR_AUG_25.RData
stored in --> /farm/home/luisie01/Riccio/Exp_1/WithNewReads/utrid/anno_ngf/anno_ngf_Aug_27.RData
stored in --> /home/rluisier/data/Riccio/Exp_1/Nov2015/Coverage/utrCov/coverageUTR_DEC_21.RData
stored in --> /home/rluisier/data/Riccio/Exp_1/Dec2016/Coverage/utrCov/myCov_Jan_26_2016_500.RData which contain myOut
stored in --> /home/rluisier/data/Riccio/Exp_1/Dec2016/Coverage/utrCov/coverage_detection.RData which contain tot.detect, myselG, myCov500, Lngf, Lnt3 and Ltot
stored in --> /home/rluisier/data/Riccio/Exp_1/Dec2016/utrid/anno_ngf/anno_ngf_Jan_27.RData
stored in -->  /home/rluisier/data/Riccio/Exp_1/Dec2016/utrid/Lngf_Jan_27.gtf
stored in --> /home/rluisier/data/Riccio/Exp_1/Dec2016/utrid/anno_ngf_stringent/anno_ngf_March_12.RData which contains "anno_ngf","ngfGRS"


code in   --> dependencies/create_anno_tot.R
stored in --> /home/rluisier/data/Riccio/Exp_1/Dec2016/utrid/anno_tot_stringent/anno_tot_March_12.RData which contains "anno_tot","totGRS"


#********************************************************
# D. Characterise output of pipeline (no.iso, length)
#
#********************************************************
code in   --> dependencies/characteriseOutput.R
stored in --> /home/rluisier/data/Riccio/Exp_1/Dec2016/utrid/anno_ngf_stringent
remark    --> should be continued


#********************************************************
# E. Comparison with PA.db  (with distribution of distance to next PAS)
#
#********************************************************
code in   --> dependencies/characterisePAdb.R #done
stored in -->  /home/rluisier/data/Riccio/Exp_1/Dec2016/utrid/padb_stringent

#********************************************************
# F. Content in term of PAS (enrichment towards high iso.no; distribution of distance between PAS)
#
#********************************************************
code in   --> dependencies/characterisePAS.R
stored in --> /home/rluisier/data/Riccio/Exp_1/Dec2016/utrid/pas_stringent
remark    --> TO BE CONTINUED; especially in regards of the positioning of the PAS motifs along the 3' UTR


#********************************************************
# H. Conservation score analysis
#
#********************************************************
code in    --> dependencies/PhastConsAnalysis.R
stored in  --> /home/rluisier/data/Riccio/Exp_1/Dec2016/utrid/PhastCons
remark    --> TO BE CONTINUED



