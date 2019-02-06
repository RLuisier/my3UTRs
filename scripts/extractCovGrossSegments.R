#!/bin/R-3.0.2/bin/Rscript

library(grDevices)
library(Rsamtools)
library(IRanges)
library(GenomicRanges)
library(rtracklayer)

# A. Load arguments
args       <- commandArgs(TRUE)
bamfile    <- args[1]
outdir     <- args[2]
chr        <- args[3]
mysample   <- args[4]
grosseg    <- args[5]
withcor    <- args[6]
pairs      <- args[7]
frac       <- as.numeric(as.character(args[8]))
protocol   <- args[9]
print(args)

# B. Load annotation file and extract Exons

# B.1 Import annotation of the gross segments
mySegmets           <- import.gff(grosseg,format="gtf")
myUTR               <- mySegmets[which(seqnames(mySegmets)==chr),]
out                 <- paste(mysample,".",chr,".RData",sep="")
names(myUTR)        <- as.character(myUTR$txID)

print("loaded annotation")

# C. Compute coverage

if(pairs=="TRUE"){
    # B.2 Load BAM file
    chrlength       <- scanBamHeader(bamfile)[[1]]$targets[[chr]]
    chrgr           <- GRanges(seqnames=chr,ranges=IRanges(start=1, end=chrlength))
    gr              <- readGAlignmentPairsFromBam(file=bamfile, param=ScanBamParam(which=chrgr))
    print("BAM loaded")
    dir.create(outdir)
    setwd(outdir)
    # C.1 Select only those reads from "+" strand
    m1              <- first(gr)
    m2              <- last(gr)
    if(protocol=="reverse"){
        sel1            <- which(strand(m1)=="-")
        sel2            <- which(strand(m2)=="+")
    }
    if(protocol!="reverse"){
        sel1            <- which(strand(m1)=="+")
        sel2            <- which(strand(m2)=="-")
    }
    if(sum(sel1!=sel2)>0){print("hey then there is an issue!")}
    m1              <- m1[sel1]
    m2              <- m2[sel2]
    mtot            <- c(m1,m2)
    covpos          <- coverage(mtot,method="sort")[[chr]]

    # C.2 Focus on those ranges that overlap with the one of interest
    utrPos   <- myUTR[which(as.character(strand(myUTR))=="+"),]
    over     <- findOverlaps(query=utrPos,subject=mtot,ignore.strand=T)
    tX       <- names(utrPos)[unique(queryHits(over))]
    subUTR   <- utrPos[which(names(utrPos)%in%tX)]

    # C.3 Create one vector of coverage per tX
    myviews1       <- Views(subject=covpos, as(subUTR,"RangesList")[[chr]])
    myCovPos       <- lapply(myviews1,function(x)return(as.data.frame(x)[,1]))
    names(myCovPos)<- names(subUTR)

    # C.4 Select only those reads from "-" strand
    m1              <- first(gr)
    m2              <- last(gr)
    if(protocol=="reverse"){
        sel1 <- which(strand(m1)=="+")
        sel2 <- which(strand(m2)=="-")
    }
    if(protocol!="reverse"){
        sel1 <- which(strand(m1)=="-")
        sel2 <- which(strand(m2)=="+")
    }
    if(sum(sel1!=sel2)>0){print("hey then there is an isse!")}
    m1              <- m1[sel1]
    m2              <- m2[sel2]
    mtot            <- c(m1,m2)
    covneg          <- coverage(mtot,method="sort")[[chr]]

    # C.5 Focus on those ranges that overlap with the one of interest
    utrNeg   <- myUTR[which(as.character(strand(myUTR))=="-"),]
    over     <- findOverlaps(query=utrNeg,subject=mtot,ignore.strand=T)
    tX       <- names(utrNeg)[unique(queryHits(over))]
    subUTR   <- utrNeg[which(names(utrNeg)%in%tX)]

    # C.6 Create one vector of coverage per tX
    myviews2       <- Views(subject=covneg, as(subUTR,"RangesList")[[chr]])
    myCovNeg       <- lapply(myviews2,function(x)return(as.data.frame(x)[,1]))
    names(myCovNeg)<- names(subUTR)
    #l1 <- unlist(lapply(myCovNeg,function(x)return(length(x))))
    #l2 <- width(subUTR)
    #sum(l1-l2) -->ok

    if(withcor=="FALSE"){
        # D. Merge POS and NEG and create output
        myCov          <- append(myCovNeg,myCovPos)
        idx            <- match(names(myCov),names(myUTR))
        myoutput       <- list(
                                chr=rep(chr,length(idx)),
                                start=start(myUTR)[idx],
                                end=end(myUTR)[idx],
                                strand=as.character(strand(myUTR))[idx],
                                txID=names(myCov),
                                cov=myCov
                )

        print("finished with creation of outut")

        # E. Save file
        dir.create(outdir)
        setwd(outdir)
        save(list="myoutput",file=out)
        print("Finished!")
    }

    if(withcor=="TRUE"){
        # C.1 Select only those reads from "-" strand
        m1              <- first(gr)
        m2              <- last(gr)
        if(protocol=="reverse"){
            sel1 <- which(strand(m1)=="+")
            sel2 <- which(strand(m2)=="-")
        }

        if(protocol!="reverse"){
            sel1 <- which(strand(m1)=="-")
            sel2 <- which(strand(m2)=="+")
        }
        if(sum(sel1!=sel2)>0){print("hey then there is an issue!")}
        m1              <- m1[sel1]
        m2              <- m2[sel2]
        mtot            <- c(m1,m2)
        covposf         <- coverage(mtot,method="sort")[[chr]]

        # C.2 Focus on those ranges that overlap with the one of interest
        overf     <- findOverlaps(query=utrPos,subject=mtot,ignore.strand=T)
        tXf       <- names(utrPos)[unique(queryHits(overf))]
        subUTRf   <- utrPos[which(names(utrPos)%in%tXf)]

        # C.3 Create one vector of coverage per tX
        myviews1f       <- Views(subject=covposf, as(subUTRf,"RangesList")[[chr]])
        myCovPosf       <- lapply(myviews1f,function(x)return(as.data.frame(x)[,1]))
        names(myCovPosf)<- names(subUTRf)
        #l1 <- unlist(lapply(myCovPosf,function(x)return(length(x))))
        #l2 <- width(subUTRf)
        #sum(l1-l2) -->ok


        # Positive strand correction
        tocorr      <- match(names(myCovPosf),names(myCovPos))
        myCovPosf   <- myCovPosf[!is.na(tocorr)]
        myCovPost   <- myCovPos[tocorr[!is.na(tocorr)]]

        CorCovPos <- function(IX){
           temp        <- round(myCovPost[[IX]]-(frac/(1-frac^2))*(myCovPosf[[IX]]-frac*myCovPost[[IX]]))
           temp[temp<0]<-0
           return(temp)
        }
        myCovPos[tocorr[!is.na(tocorr)]] <- lapply(c(1:length(myCovPost)),CorCovPos)


        # C.4 Select only those reads from "+" strand
        m1              <- first(gr)
        m2              <- last(gr)
        if(protocol=="reverse"){
            sel1 <- which(strand(m1)=="-")
            sel2 <- which(strand(m2)=="+")
        }
        if(protocol!="reverse"){
            sel1 <- which(strand(m1)=="+")
            sel2 <- which(strand(m2)=="-")
        }
        if(sum(sel1!=sel2)>0){print("hey then there is an isse!")}
        m1              <- m1[sel1]
        m2              <- m2[sel2]
        mtot            <- c(m1,m2)
        covnegf         <- coverage(mtot,method="sort")[[chr]]

        # C.5 Focus on those ranges that overlap with the one of interest
        overf     <- findOverlaps(query=utrNeg,subject=mtot,ignore.strand=T)
        tXf       <- names(utrNeg)[unique(queryHits(overf))]
        subUTRf   <- utrNeg[which(names(utrNeg)%in%tXf)]

        # C.6 Create one vector of coverage per tX
        myviews2f       <- Views(subject=covnegf, as(subUTRf,"RangesList")[[chr]])
        myCovNegf       <- lapply(myviews2f,function(x)return(as.data.frame(x)[,1]))
        names(myCovNegf)<- names(subUTRf)


        # Negative strand correction
        tocorr      <- match(names(myCovNegf),names(myCovNeg))
        myCovNegf   <- myCovNegf[!is.na(tocorr)]
        myCovNegt   <- myCovNeg[tocorr[!is.na(tocorr)]]

        CorCovNeg <- function(IX){
           temp        <- round(myCovNegt[[IX]]-(frac/(1-frac^2))*(myCovNegf[[IX]]-frac*myCovNegt[[IX]]))
           temp[temp<0]<-0
           return(temp)
        }
        myCovNeg[tocorr[!is.na(tocorr)]] <- lapply(c(1:length(myCovNegt)),CorCovNeg)


        # D. Merge POS and NEG and create output
        myCov          <- append(myCovNeg,myCovPos)
        idx            <- match(names(myCov),names(myUTR))
        myoutput       <- list(
                                chr=rep(chr,length(idx)),
                                start=start(myUTR)[idx],
                                end=end(myUTR)[idx],
                                strand=as.character(strand(myUTR))[idx],
                                txID=names(myCov),
                                cov=myCov
                )

        print("finished with creation of outut")

        # E. Save file
        save(list="myoutput",file=out)
        print("Finished!")

    }


}

if(pairs=="FALSE"){
    # B.2 Load BAM file
    chrlength       <- scanBamHeader(bamfile)[[1]]$targets[[chr]]
    chrgr           <- GRanges(seqnames=chr,ranges=IRanges(start=1, end=chrlength))
    gr              <- readGAlignmentFromBam(file=bamfile, param=ScanBamParam(which=chrgr))
    print("BAM loaded")
    dir.create(outdir)
    setwd(outdir)

    # C.1 Coverage positive strand -- I should check this to
    if(protocol=="reverse"){
        mtot            <- gr[which(strand(gr)=="-")]
        covpos          <- coverage(mtot,method="sort")[[chr]]
    }

    if(protocol!="reverse"){
        mtot            <- gr[which(strand(gr)=="+")]
        covpos          <- coverage(mtot,method="sort")[[chr]]
    }


    # C.2 Focus on those ranges that overlap with the one of interest
    utrPos   <- myUTR[which(as.character(strand(myUTR))=="+"),]
    over     <- findOverlaps(query=utrPos,subject=mtot,ignore.strand=T)
    tX       <- names(utrPos)[unique(queryHits(over))]
    subUTR   <- utrPos[which(names(utrPos)%in%tX)]

    # C.3 Create one vector of coverage per tX
    myviews1       <- Views(subject=covpos, as(subUTR,"RangesList")[[chr]])
    myCovPos       <- lapply(myviews1,function(x)return(as.data.frame(x)[,1]))
    names(myCovPos)<- names(subUTR)


    # C.4 Select only those reads from "-" strand
    if(protocol=="reverse"){
        mtot            <- gr[which(strand(gr)=="+")]
        covneg          <- coverage(mtot,method="sort")[[chr]]
    }

    if(protocol!="reverse"){
        mtot            <- gr[which(strand(gr)=="-")]
        covneg          <- coverage(mtot,method="sort")[[chr]]
    }

    # C.5 Focus on those ranges that overlap with the one of interest
    utrNeg   <- myUTR[which(as.character(strand(myUTR))=="-"),]
    over     <- findOverlaps(query=utrNeg,subject=mtot,ignore.strand=T)
    tX       <- names(utrNeg)[unique(queryHits(over))]
    subUTR   <- utrNeg[which(names(utrNeg)%in%tX)]

    # C.6 Create one vector of coverage per tX
    myviews2       <- Views(subject=covneg, as(subUTR,"RangesList")[[chr]])
    myCovNeg       <- lapply(myviews2,function(x)return(as.data.frame(x)[,1]))
    names(myCovNeg)<- names(subUTR)


    if(withcor=="FALSE"){
        # D. Merge POS and NEG and create output
        myCov          <- append(myCovNeg,myCovPos)
        idx            <- match(names(myCov),names(myUTR))
        myoutput       <- list(
                        chr=rep(chr,length(idx)),
                        start=start(myUTR)[idx],
                        end=end(myUTR)[idx],
                        strand=as.character(strand(myUTR))[idx],
                        txID=names(myCov),
                        cov=myCov
                )

        print("finished with creation of outut")

        # E. Save file
        dir.create(outdir)
        setwd(outdir)
        save(list="myoutput",file=out)
        print("Finished!")
    }

    if(withcor=="TRUE"){
        # C.1 Select only those reads from "-" strand
        if(protocol=="reverse"){
            mtot            <- gr[which(strand(gr)=="+")]
            covposf         <- coverage(mtot,method="sort")[[chr]]
        }

        if(protocol!="reverse"){
            mtot            <- gr[which(strand(gr)=="-")]
            covposf          <- coverage(mtot,method="sort")[[chr]]
        }

        # C.2 Focus on those ranges that overlap with the one of interest
        overf     <- findOverlaps(query=utrPos,subject=mtot,ignore.strand=T)
        tXf       <- names(utrPos)[unique(queryHits(overf))]
        subUTRf   <- utrPos[which(names(utrPos)%in%tXf)]

        # C.3 Create one vector of coverage per tX
        myviews1f       <- Views(subject=covposf, as(subUTRf,"RangesList")[[chr]])
        myCovPosf       <- lapply(myviews1f,function(x)return(as.data.frame(x)[,1]))
        names(myCovPosf)<- names(subUTRf)

        # Positive strand correction
        tocorr      <- match(names(myCovPosf),names(myCovPos))
        myCovPosf   <- myCovPosf[!is.na(tocorr)]
        myCovPost   <- myCovPos[tocorr[!is.na(tocorr)]]

        CorCovPos <- function(IX){
           temp        <- round(myCovPost[[IX]]-(frac/(1-frac^2))*(myCovPosf[[IX]]-frac*myCovPost[[IX]]))
           temp[temp<0]<-0
           return(temp)
        }
        myCovPos[tocorr[!is.na(tocorr)]] <- lapply(c(1:length(myCovPost)),CorCovPos)


        # C.4 Select only those reads from "-" strand
        if(protocol=="reverse"){
             mtot            <- gr[which(strand(gr)=="-")]
             covnegf          <- coverage(mtot,method="sort")[[chr]]
        }

        if(protocol!="reverse"){
            mtot            <- gr[which(strand(gr)=="+")]
            covnegf         <- coverage(mtot,method="sort")[[chr]]
        }


        # C.5 Focus on those ranges that overlap with the one of interest
        overf     <- findOverlaps(query=utrNeg,subject=mtot,ignore.strand=T)
        tXf       <- names(utrNeg)[unique(queryHits(overf))]
        subUTRf   <- utrNeg[which(names(utrNeg)%in%tXf)]

        # C.6 Create one vector of coverage per tX
        myviews2f       <- Views(subject=covnegf, as(subUTRf,"RangesList")[[chr]])
        myCovNegf       <- lapply(myviews2f,function(x)return(as.data.frame(x)[,1]))
        names(myCovNegf)<- names(subUTRf)

        # Negative strand correction
        tocorr      <- match(names(myCovNegf),names(myCovNeg))
        myCovNegf   <- myCovNegf[!is.na(tocorr)]
        myCovNegt   <- myCovNeg[tocorr[!is.na(tocorr)]]

        CorCovNeg <- function(IX){
           temp        <- round(myCovNegt[[IX]]-(frac/(1-frac^2))*(myCovNegf[[IX]]-frac*myCovNegt[[IX]]))
           temp[temp<0]<-0
           return(temp)
        }
        myCovNeg[tocorr[!is.na(tocorr)]] <- lapply(c(1:length(myCovNegt)),CorCovNeg)


        # D. Merge POS and NEG and create output
        myCov          <- append(myCovNeg,myCovPos)
        idx            <- match(names(myCov),names(myUTR))
        myoutput       <- list(
                        chr=rep(chr,length(idx)),
                        start=start(myUTR)[idx],
                        end=end(myUTR)[idx],
                        strand=as.character(strand(myUTR))[idx],
                        txID=names(myCov),
                        cov=myCov
                )

        print("finished with creation of outut")

        # E. Save file
        save(list="myoutput",file=out)
        print("Finished!")

    }


}

