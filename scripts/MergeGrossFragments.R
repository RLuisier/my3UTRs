library(grDevices)
library(Rsamtools)
library(IRanges)
library(GenomicRanges)
library(rtracklayer)

setwd("./utrid/GrossSegments/Extended")

# A. Import new annotation
samples <- sort(as.character(read.delim("./SAMPLES.txt",header=F)[,1]))
foi1    <- unlist(lapply(samples,function(x)return(paste("./utrid/GrossSegments/GFF/",x,".6.gff",sep=""))))
foi2    <- unlist(lapply(samples,function(x)return(paste("./utrid/GrossSegments/Extended/",x,".extended.gff",sep=""))))
fois    <- c(foi1,foi2)
myS     <- do.call(lapply(fois,function(x)return(import.gff(con=x, version = c("1"), genome = "rn5"))),what=c)


# B. Import known 3' UTR and annotation of transcript
load("./annotation/rn5/GRanges_comprehensive_transcriptome_rat_24_nov_2015.RData")


# C. Extract the longest 3'UTR isoforms from new annotation
ids          <- as.character(mcols(g3utr)$X.transcript_id)
group        <- as.character(mcols(g3utr)$group)
txID         <- ids[match(as.character(mcols(myS)$group),group)]
myS$txID     <- txID
dL           <- width(myS)-width(g3utr)[match(myS$txID,ids)]

newS         <- GRanges(seqnames = c(seqnames(myS),seqnames(g3utr)),
                       ranges   = c(ranges(myS),ranges(g3utr)),
                       strand   = c(strand(myS),strand(g3utr)),
                       txID     = c(myS$txID,as.character(mcols(g3utr)$X.transcript_id)),
                       dL       = c(dL,rep(0,length(g3utr)))
           )

no.iso       <- as.vector(table(newS$txID))[match(newS$txID,names(table(newS$txID)))]
newS$no.iso  <- no.iso


sub1        <- newS[newS$no.iso==1,]
sub2        <- newS[newS$no.iso>1,]
myTX        <- unique(sub2$txID)
temp        <- cbind(sub2$dL,sub2$txID,c(1:length(sub2)))
theLongest  <- do.call(lapply(myTX, function(x){
                            ix  <- which(temp[,2]==x)
                            sub <- temp[ix,]
                            return(sub[sort(as.numeric(sub[,1]),decreasing=T,index.return=T)$ix[1],3])
                         }
                         ),what=c)

myLongest              <- c(sub1,sub2[as.numeric(theLongest),])
POS                    <- as.character(strand(myLongest))=="+"
end(myLongest)[POS]    <- end(myLongest)[POS] + 100
start(myLongest)[!POS] <- start(myLongest)[!POS] - 100
start(myLongest)[start(myLongest)<1]<-1

export.gff(object=myLongest,con="./annotation/rn5/GrossFragments.gtf",format="gtf")
