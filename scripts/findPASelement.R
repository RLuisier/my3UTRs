library(grDevices)
library(Rsamtools)
library(IRanges)
library(GenomicRanges)
library(rtracklayer)
library(Biostrings)

# A. Load arguments
args       <- commandArgs(TRUE)
outdir     <- args[1]
IX         <- args[2]
myfasta    <- args[3]
region     <- args[4]
print(args)


# B. Load Fasta regions
target     <- readDNAStringSet(myfasta)
myGR       <- import.gff(region,format="gtf")
names(myGR)<- as.character(myGR$txID)
sel        <- match(names(target),names(myGR))
myGR       <- myGR[sel,]

# C. Load motif of interest
motif       <- RNAStringSet(as.character(read.delim("./annotation/PAS_paDB_v2.txt",sep="\t")[IX,2]))
print("loaded 1")
motif       <- DNAStringSet(motif)
print("loaded motif")
motifID     <- as.character(motif)

# D. Find motif
print("start finding motif")
my.res1        <- lapply(target,function(x)return(matchPattern(pattern=motif[[1]],subject = x, fixed = TRUE)))
print("finished with myRes1")

# E. Select the longest isoform per tX
uniqueID    <- names(target)
myL         <- nchar(target)
names(myL)  <- names(target)
txID        <- unlist(lapply(names(myL),function(x)unlist(strsplit(x,split="\\."))[1]))
mygroups    <- as.factor(as.character(txID))
mymax       <- tapply(myL,INDEX=mygroups,FUN=function(x)return(names(x)[x==max(x)[1]]))
oi          <- as.vector(unlist(mymax))

# F. Extract all positions along the transcriptome
selLong       <-  names(myL)%in%oi
selLong       <-  selLong[match(names(my.res1),names(myL))]
myGR          <-  myGR[match(names(my.res1),names(myGR)),]

subGR        <-  myGR[selLong,]
myStart      <-  lapply(my.res1[selLong],function(x)return(start(x)))
length.start <-  unlist(lapply(myStart,length))

subGR        <- subGR[length.start>0,]
myStart      <- myStart[length.start>0]
length.start <- length.start[length.start>0]


POS           <- as.character(strand(subGR))=="+"

start.gr       <- start(subGR)[POS]
myStartp       <- myStart[POS]
length.startp  <- length.start[POS]
chr.gr.p       <- as.character(seqnames(subGR))[POS]

end.gr         <- end(subGR)[!POS]
myStartn       <- myStart[!POS]
length.startn  <- length.start[!POS]
chr.gr.n       <- as.character(seqnames(subGR))[!POS]


start.p   <- unlist(lapply(c(1:length(length.startp)),function(x)return(rep(start.gr[x],length.startp[x]))))
start.p   <- start.p + unlist(myStartp)
end.p     <- start.p + nchar(motif) - 1
chr.p     <- unlist(lapply(c(1:length(length.startp)),function(x)return(rep(chr.gr.p[x],length.startp[x]))))
strand.p  <- rep("+",length(chr.p))

end.n     <- unlist(lapply(c(1:length(length.startn)),function(x)return(rep(end.gr[x],length.startn[x]))))
end.n     <- end.n -  unlist(myStartn)
start.n   <- end.n - nchar(motif) + 1
chr.n     <- unlist(lapply(c(1:length(length.startn)),function(x)return(rep(chr.gr.n[x],length.startn[x]))))
strand.n  <- rep("-",length(chr.n))

chr       <- c(chr.p,chr.n)
strand    <- c(strand.p,strand.n)
start     <- c(start.p,start.n)
end       <- c(end.p,end.n)
ID        <- data.frame(rep(motifID,length(end)))

myPAS      <- GRanges(
                   seqnames = Rle(chr),
                   ranges   = IRanges(start=start,end=end),
                   strand   = Rle(strand),
                   ID
              )

export.gff(object=myPAS,con=paste(outdir,motifID,".gtf",sep=""),format="gtf")
save(list="myPAS",file=paste(outdir,motifID,".RData",sep=""))











