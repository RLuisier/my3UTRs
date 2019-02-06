#!/bin/R-3.0.2/bin/Rscript

library(grDevices)
library(Rsamtools)
library(IRanges)
library(GenomicRanges)
library(rtracklayer)
library(gplots)
require(EBS)
require(Segmentor3IsBack)



# A. Load arguments
args       <- commandArgs(TRUE)
outdir     <- args[1]
tardir     <- args[2]
chr        <- args[3]
region     <- args[4]
method     <- args[5]
annotation <- args[6]
limD       <- as.numeric(as.character(args[7]))
winD       <- as.numeric(as.character(args[8]))
noSeg      <- as.numeric(as.character(args[9]))
isLog      <- as.logical(as.character(args[10]))
minCov     <- as.numeric(as.character(args[11]))

print(args)

print("load annotation")
load(annotation)
print("load dependencies")
source("./scripts/nested_FUN_identify_pAs.R")



#B.Prepare annotation file
print("Prepare annotation file")
files        <- list.files(tardir)
files        <- files[grep(files,pattern=".RData")]
foi          <- data.frame(files,do.call(what=rbind,args=lapply(files,function(x)return(unlist(strsplit(x,split="[-.\\_]+"))))))[,c(1:5)]
rm("files")
colnames(foi)<- c("file","condition","origin","replicate","chr")

# D. Import Gross Segments
print("Import gross segment")
print(region)
myregion                        <- import.gff(region,format="gtf")
print("region imported")
myUTR                           <- myregion[which(seqnames(myregion)==chr),]
names(myUTR)                    <- as.character(myUTR$txID)
mcols(myUTR)$X.transcript_id.   <- myUTR$txID
out                             <- paste(outdir,"grUTR.",chr,".RData",sep="")


#C. Load data
print("Load data")
sel         <- which(foi$chr==chr)
foi         <- foi[sel,]#should result in 4 samples
myres       <-  list()
for(i in c(1:length(sel))){
    load(paste(tardir,foi$file[i],sep=""))
    myres[[i]] <- myoutput
}
names(myres)<- paste(foi$condition,foi$origin,foi$replicate,sep=".")



#E. Identify pA sites
print("Identify pA sites")
ids      <- do.call(args=mapply(A=myres,B=names(myres),function(A,B)return(cbind(A$txID,rep(B,length(A$txID))))),what=rbind)
tokeep   <- ids[,1]%in%as.character(mcols(myUTR)$X.transcript_id.)[width(myUTR)>0]
ids      <- ids[tokeep,]
groups   <- as.factor(ids[,1])
myiso    <- list()

for(s in c(1:length(levels(groups)))){
    print(paste("txID=",s,sep=""))
    myiso[[s]] <- mergePAsTissue(tempTX=levels(groups)[s],method=method,windowsize=winD,NumSeg=noSeg,doLog=isLog,limDist=limD,MinCov=minCov)
}

mytest <- mergePAsTissue(tempTX="ENSRNOT00000016995",method=method,windowsize=winD,NumSeg=noSeg,doLog=isLog,limDist=limD,MinCov=minCov)


# F. Create GRange object
myGRobject <- do.call(args=myiso,what=c)

#export.gff(object=myGRobject,con="tempGFFchr6.gtf",format="gtf")


# G. Save file
save(list="myGRobject",file=out)
print("Finished!")
