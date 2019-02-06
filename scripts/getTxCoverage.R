#!/bin/R-3.0.2/bin/Rscript

library(grDevices)
library(Rsamtools)
library(IRanges)
library(GenomicRanges)
library(rtracklayer)

# A. Load arguments
args        <- commandArgs(TRUE)
outdir      <- args[1]
chr         <- args[2]
annotation  <- args[3]
pathdir     <- args[4]
print(args)


# A. Load data
load(annotation)
tX    <- e75f[intersect(which(as.character(seqnames(e75f))==chr),which(as.character(mcols(e75f)$type)%in%"exon")),]
ID    <- as.character(mcols(tX)$X.transcript_id.)
tX$ID <- ID


# B. Import BAM files
files   <- system(paste("ls -d ",pathdir,"*.bam",sep=""),intern=T)
gr      <- list()
for(i in c(1:length(files))){
    bamfile          <- files[i]
    chrlength       <- scanBamHeader(bamfile)[[1]]$targets[[chr]]
    chrgr           <- GRanges(seqnames=chr,ranges=IRanges(start=1, end=chrlength))
    gr[[i]]         <- readGAlignmentPairsFromBam(file=bamfile, param=ScanBamParam(which=chrgr))
    print(i)
}

# C. Get Coverage from BAM file
cov          <- list()
myviews_exon <- list()

for(i in c(1:length(gr))){
    cov[[i]]                    <-  coverage(gr[[i]],method="sort")[[chr]]
    myviews_exon[[i]]           <-  lapply(Views(subject=cov[[i]], as(tX,"RangesList")[[chr]]),function(x)return(as.data.frame(x)[,1]))
    print(i)
}


# D. Create a single vector per tX isoform

uniqueID    <- unique(tX$ID)
allVecRes   <- list()

for(SAMPLE in c(1:length(myviews_exon))){

    myVec    <- list()
    IX       <- 1
    POS      <- as.character(strand(tX))[match(uniqueID,as.character(tX$ID))]=="+"
    id_pos   <- uniqueID[POS]
    for(i in c(1:length(id_pos))){
        ix.ID           <- which(as.character(tX$ID)%in%id_pos[i])
        myVec[[IX]]     <- do.call(myviews_exon[[SAMPLE]][ix.ID],what=c)
        IX              <- IX+1
        print(i)
    }

    IX       <- length(myVec)+1
    id_neg   <- uniqueID[!POS]
    for(i in c(1:length(id_neg))){

        ix.ID       <- which(as.character(tX$ID)%in%id_neg[i])
        myVec[[IX]] <- rev(do.call(myviews_exon[[SAMPLE]][rev(ix.ID)],what=c))
        IX          <- IX+1
     print(i)
    }
    names(myVec)        <- c(id_pos,id_neg)
    allVecRes[[SAMPLE]] <- myVec
 }


# D.1.2 Create BINS
myL  <- list()
BINS <- list()
for(i in c(1:length(allVecRes))){
    myL[[i]]  <- unlist(lapply(allVecRes[[i]],length))
    BINS[[i]] <- cut(x=myL[[i]],breaks=10^3*c(0.5,1,2,3,4,7,20),include.lowest=T)
    print(i)
}

# D.1.3 Sum up read coverage from all the tX in the same range of length
maxL    <- (10^3*c(0.5,1,2,3,4,7,20))[-1]
mySum <- list()
for(i in c(1:length(allVecRes))){
    #maxL       <- unlist(tapply(myL[[i]],INDEX=BINS[[i]],FUN=max))
    mySum[[i]] <- list()
    LEV        <- levels(BINS[[i]])
    print(i)
    for(j in c(1:length(maxL))){

        myList <- lapply(X=allVecRes[[i]][which(BINS[[i]]==LEV[j])],FUN=function(x)return(c(rep(0,maxL[j]-length(x)),x)))
        mySum[[i]][[j]] <- apply(
                                    do.call(
                                        args=myList,                                          what=rbind
                                            )
                                  ,2,sum
                                  )
        print(j)
    }
}



# D.1.4 Create density vector of length 100 for each range
GetDensityCov <- function(vec){
                    tempvec <- tapply(vec,INDEX=cut(c(1:length(vec)),100),FUN=sum)
                    return(tempvec/sum(tempvec))
            }

DensityCov   <- lapply(X=mySum,FUN=function(X)return(do.call(args=lapply(X=X,FUN=GetDensityCov),what=rbind)))

# E. Save
out <- paste(outdir,"txCoverage_",chr,".RData",sep="")
save(list=c("tX","allVecRes","DensityCov","mySum"),file=out)

print("I am done yeah!")
