print("start merging neighbors")
## A.1 Read regions of interest (per chromosome file)
print("Start read region of interest")
if(ensembl=="TRUE"){
    altchr <- gsub(chr,pattern="chr",repl="")
}
if(ensembl=="FALSE"){
    altchr <- chr
}
roi <- paste(altchr,".txt",sep="")
reg  <- read.table(roi, header=T,colClasses=c("numeric","numeric"))
nR   <- nrow(reg)
groi <- GRanges(seqnames =
            Rle(replicate(nR,chr)),
            ranges =
            IRanges(start=reg[,1],end=reg[,2]),
            strand =
            Rle(replicate(nR,strand))# TO BE MODIFIED!!!
            )
groi  <- reduce(groi)
nR    <- length(groi)
print("Finished read region of interest")

# A.2 Create GRanges from RepeatMask
print("Start RepeatMasker")
RpM    <- read.table(paste(RPM,chr,".RP",sep=""),header=T,comment.char = "&")
gRp <-
    GRanges(seqnames =
            Rle(replicate(nrow(RpM),chr)),
            ranges =
            IRanges(start=RpM$genoStart,end=RpM$genoEnd),
            strand =
            Rle(RpM$strand)
            )
# As there are some overlapping regions,
gRp <- reduce(gRp)
# A.3 Create GRanges from Gaps between the identified covered regions
print("Create GRanges from Gaps between the identified covered regions")
grga <-
    GRanges(seqnames =
            Rle(replicate(nR-1,chr)),
            ranges =
            IRanges(start=end(groi)[c(1:(nR-1))],end=start(groi)[c(2:nR)]),
            strand =
            Rle(replicate(nR-1,strand))
            )

# B. Merge regions separated by less than 150 nt of non-repetitive sequences
#
#
#   Conditions to merge two windows:
#                1/ w1---|-RRRRRRRRRRRRRRRRRR-|----w2 i.e. gap is only repeat sequence
#                2/ w1---|----RRR--RRR---RRR--|----w2 i.e. the number of non-repeat nt <150nt
#                3/ w1---|--------<150nt------|----w2 i.e. gap composed of non-repeat less than 150nt
#
#     BUT this may lead to merging two regions which are really distant given that 50% of genome is
#     covered by RepeatMask!

# B.1 Regions which are less than 150 nt
#
#
print("Identify gaps which are less than 150 nt")
TestA                                 <- replicate(length(grga),FALSE)
TestA[which(width(ranges(grga))<150)] <- TRUE

# B.2 Identify gaps that overlap with repetitive elemts (or where there is no more than 150nt in between)
#
#
print("Identify gaps that overlap with repetitive elemts (or where there is no more than 150nt in between)")
overlaps <- findOverlaps(query=grga,subject=gRp,ignore.strand=FALSE)
idxGa    <- queryHits(overlaps)
idxRp    <- subjectHits(overlaps)
#Percentage of gaps that overlap with Repeat Regions:  0.82 = length(unique(idxGa))/length(grga)
Lev    <- unique(idxGa)
Test1  <- replicate(length(Lev),FALSE)

for(j in 1:length(Lev)){
    i <- Lev[j]
    #Range of the gap of interest
    gap   <- ranges(grga[idxGa[i]])
    #Ranges of the overlapping repeated regions
    overs <- ranges(gRp[idxRp[which(idxGa==idxGa[i])]])
    #Force the extreme of the ranges of the overlapping Repeated Regions to be in range of the gap
    start(overs)[start(overs)<start(gap)]<-start(gap)
    end(overs)[end(overs)>end(gap)]      <-end(gap)
    #Test the two conditions
    Test1[j] <- (width(gap)-sum(width(overs)) <150)
}

TestB             <- replicate(length(grga),FALSE)
TestB[Lev[Test1]] <-TRUE

# B.3 Check that the gaps are spanned by at least one paired-end read IF paired-end!!!
#
#
if(pairs=="TRUE"){
    print("Check that the gaps are spanned by at least one paired-end read")
    chrgr       <- groi
    print("start reading BAM file with corresponding covered regions (take ages...)")

    gbam        <- readGAlignmentsFromBam(
        file=bamfile,
        param=ScanBamParam(which=chrgr,what=c("flag", "mrnm", "mpos","cigar")),
        use.names=TRUE
    )
    print("finished reading BAM file with corresponding covered regions")

    #Check first that names are not written with additional stuff
    names(gbam) <- sub("_.*$", "", names(gbam))
    print("finished changing names")
    gbam2      <- makeGAlignmentPairs(gbam)
    print("finished making GAlignment BAM file with corresponding covered regions")
    gbam2          <- gbam2[isProperPair(gbam2),]
    gbamp         <- granges(gbam2)
    print("finished making GRanges BAM file with corresponding covered regions")


    gOver         <- findOverlaps(query=groi,subject=gbamp,ignore.strand=FALSE)
    grouping      <- factor(queryHits(gOver))
    Lev           <- as.numeric(levels(grouping))
    print("start looking for gaps covered by pair...")
    TestC         <- replicate(length(grga),FALSE)
    idx1          <- queryHits(gOver)
    idx2          <- subjectHits(gOver)

    for(i in 1:(length(Lev)-1)){
        TestC[Lev[i]]<- sum(idx2[which(idx1==Lev[i]+1)]%in%idx2[which(idx1==Lev[i])])>0
    }
}


if(pairs=="FALSE"){
    TestC <- rep(TRUE,length(TestA))
}

TestF             <- which((TestA|TestB)&TestC)
print(paste("N0 of regions merged according to size smaller than 150nt only:", sum(TestA)))
print(paste("N0 of regions merged according to repetitive regions:", sum(TestB)))
print(paste("Gaps covered by at least one paired-end reads:", sum(TestC)))
print(paste("N0 of regions merged according to both test:", sum((TestA|TestB)&TestC)))
print(paste("N0 of gaps:", length(grga)))


# B.4 Merge these gapps
print("mergeGaps")
if(length(TestF)>0){
    tempS  <- start(groi)[-(TestF+1)]
    tempE  <- end(groi)[-(TestF)]
    nR     <- length(tempS)

    print(length(tempS))
    print(length(tempE))
    print(chr)

    groiv1 <- GRanges(seqnames =
            Rle(replicate(nR,chr)),
            ranges =
            IRanges(start=tempS,end=tempE),
            strand =
            Rle(replicate(nR,strand))
        )
}

if(length(TestF)==0){
    groiv1<- groi
}

library(rtracklayer)
setwd(tempdir)
export.gff(object=groi,con="groi.gff",format="gff")
export.gff(object=groiv1,con="groiv1.gff",format="gff")
print("terminated merging Neighbours")
