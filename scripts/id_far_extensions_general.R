#!/bin/R-3.0.2/bin/Rscript
library(grDevices)
library(Rsamtools)
library(IRanges)
library(GenomicRanges)
library(rtracklayer)
library(mclust)
library(geneplotter)


#Load arguments
args        <- commandArgs(TRUE)
infile      <- as.character(args[1])
wd          <- as.character(args[2])
name        <- as.character(args[3])
minFS       <- as.numeric(as.character(args[4]))
checkfile   <- as.character(args[5])
maxDist     <- as.numeric(as.character(args[6]))
annotation  <- as.character(args[7])


print("here are the arguments:")
print(args)
print("So now that everything is ok, let's start!")

setwd(wd)
tempdir <-paste("new_extensions_",name,sep="")
dir.create(tempdir)
print("Start identification")
setwd(tempdir)

# A. Create GRanges objects
#
#
# A.1 Continue from gross merged segments: groiv1
#	  groiv1	-- Granges of regions of interest which have been merged
groiv1 <- import.gff(infile)
groiv1 <- groiv1[width(groiv1)>minFS,]
check  <- import.gff(con=checkfile,version="1")
# A.2 Load Granges from various sources:
#   e75p	--	GRange object containing all transcripts from Ensembl  without any filtering
#			i.e. initial gff file download from Ensembl encoded in GRange object.
#   e75f	--	GRange object containing all transcript encoding a PROTEIN with exons from v75 Ensbl.
#   g3utr	-- 	GRange object containing 3'UTR for all protein coding transcripts from Ensembl v75
#			i.e. from end of stop codon until end of the transcript).
#			PLease note however that thes UTRs have introns included.
#   g5utr	--	GRange object which contains 5'UTR for all protein coding transcripts from Ensembl v75
#			i.e. from start of the transcript until start of start codon.
#			PLease note however that thes UTRs have introns included.
#   g3introns 	--      GRange object of introns contained in g3utr (3'UTRs).
#   g5introns	--	GRange object of introns contained in g5utr (5'UTRs).
#       ---takes some time!
print("Load required GRanges: UTRs and transcriptomes form Ensembl.")
load(annotation)

# Step 1: remove all segments that overlap with ANYTHING.
print("Step 1: remove all segments that overlap with ANYTHING.")
uOver <- findOverlaps(query=groiv1,subject=e75p,ignore.strand=FALSE)
gF1   <- groiv1[-unique(queryHits(uOver)),]
export.gff(object=gF1,con=paste("gF1_",name,".gff",sep=""),format="gff")



# Step 2: identify nearest feature using strand specificity
print("Step 2: find nearest neighbour and select only those fragments for which nearest neighbour is 3' UTR wihtin acceptable distance.")
# A. Identify feature which precedes the extension or which feature is follwed by the extension
# ``follow'' works like this: for each range of 'x', returns index of subject IF there is a range in subject that is BEFORE the range in question.
gOver             <- nearest(x=gF1,subject=e75p,select="all",ignore.strand=FALSE)
idxU              <- queryHits(gOver)
idxG              <- subjectHits(gOver)
mydist            <- as.data.frame(distance(x=gF1[idxU],y=e75p[idxG]))#Distance to all features to get an idea


IDX    <- which(seqnames(gF1)=="chr19"&start(gF1)>54410236&end(gF1)<54414699)
IDX.p  <- which(idxU==IDX)

# B. Select only features that are indeed 3'UTR

#B.1 Compute distance to any nearest feature
gOver             <- nearest(x=gF1,subject=e75p,select="all",ignore.strand=FALSE)
idxU              <- queryHits(gOver)
idxG              <- subjectHits(gOver)
tempdist          <- as.vector(distance(x=gF1[idxU],y=e75p[idxG]))
names(tempdist)   <- c(1:length(tempdist))
tokeep            <- as.numeric(tapply(tempdist,INDEX=factor(idxU,levels=as.character(c(1:length(gF1)))),
                                                FUN=function(Z)return(names(Z)[which(Z==min(Z))[1]])))
d_all             <- tempdist[tokeep]

#B.2 Compute distance to nearest 3' UTR
gOver             <- nearest(x=gF1,subject=g3utr,select="all",ignore.strand=FALSE)
idxU              <- queryHits(gOver)
idxG              <- subjectHits(gOver)
tempdist          <- as.vector(distance(x=gF1[idxU],y=g3utr[idxG]))
names(tempdist)   <- c(1:length(tempdist))
tokeep            <- as.numeric(tapply(tempdist,INDEX=factor(idxU,levels=as.character(c(1:length(gF1)))),
                                                FUN=function(Z)return(names(Z)[which(Z==min(Z))[1]])))
d_utr             <- tempdist[tokeep]
idxU              <- idxU[tokeep]
idxG              <- idxG[tokeep]

print(sum(diff(rev(idxU))>1))

#B.3 Select those fragment which closest feature is a 3' UTR
sel1 <- d_utr<=d_all
sel2 <- d_utr<=maxDist
self <- sel1&sel2
idxG <- idxG[self]


# C. Compute distance to closest feature in general
print("Filter n0.2 : Nearest feature is within acceptable distance and that are 3'UTR")
mydistsub         <- as.data.frame(distance(x=gF1[idxU],y=e75p[idxG]))#Distance to all features to get an idea
gF2  <- gF1[self,]
export.gff(object=gF2,con=paste("gF2_",name,".gff",sep=""),format="gff")


# D. Select the threshold of 'acceptable distance' assuming very far elements are non-annotated new gene while 'close' distance are more likely to be 3'UTR extensions.
data     <- d_all+1
data_log <- log10(data)
bimdens  <- densityMclust(data=data_log,G=2)
Lim      <- 10^qnorm(0.9,mean=bimdens$parameters$mean[1],sd=sqrt(bimdens$parameters$variance$sigmasq[1]))

pdf(paste(name,"_plot_fit_2_components_distance_to_3UTR.pdf",sep=""))
x  <- seq(from=0, to=ceiling(max(data_log)), length=100)
if(length(bimdens$parameters$variance$sigmasq)>1){
    hx1 <- dnorm(x,mean=bimdens$parameters$mean[1],sd=sqrt(bimdens$parameters$variance$sigmasq[1]))
    hx2 <- dnorm(x,mean=bimdens$parameters$mean[2],sd=sqrt(bimdens$parameters$variance$sigmasq[2]))
    hist(data_log,breaks=100,col=rgb(0,0,0,alpha=0.2),freq=FALSE,ylim=c(0,max(c(hx1,hx2))),xlab="",ylab="",main="",xaxt="n")
    lines(x,hx2 , lwd=2, col="blue")
    lines(x,hx1 , lwd=2, col="red")
    abline(v=qnorm(0.9,mean=bimdens$parameters$mean[1],sd=sqrt(bimdens$parameters$variance$sigmasq[1])),lty=2,col="black")
    mtext(side=1,line=2,text="distance to next 3UTR")
    mtext(side=2,line=2,text="density")
    #Plot nice log-log axes
    full=c("1","10","100","1000","10000","100000","1000000","10000000","100000000")
    y1   <- range(x)
    pow  <- seq(y1[1], y1[2]+1)
    ticksaty <- as.vector(sapply(pow, function(p) (1:10)*10^p))
    laby <- replicate(length(ticksaty),"")
    idy = match(as.numeric(full),ticksaty)
    idy = idy[!is.na(idy)]
    laby[idy]<-full[c(1:length(match(as.numeric(full),ticksaty)))]
    axis(side=1, at=ticksaty, labels=TRUE, tcl=-0.25, lwd=0, lwd.ticks=0.5,cex=0.8)
}
if(length(bimdens$parameters$variance$sigmasq)<1){
    hx1 <- dnorm(x,mean=bimdens$parameters$mean[1],sd=sqrt(bimdens$parameters$variance$sigmasq))
    hx2 <- dnorm(x,mean=bimdens$parameters$mean[2],sd=sqrt(bimdens$parameters$variance$sigmasq))
    hist(data_log,breaks=100,col=rgb(0,0,0,alpha=0.2),freq=FALSE,ylim=c(0,max(c(hx1,hx2))),xlab="",ylab="",main="",xaxt="n")
    lines(x,hx2 , lwd=2, col="blue")
    lines(x,hx1 , lwd=2, col="red")
    abline(v=qnorm(0.9,mean=bimdens$parameters$mean[1],sd=sqrt(bimdens$parameters$variance$sigmasq[1])),lty=2,col="black")
    mtext(side=1,line=2,text="coverage")
    mtext(side=2,line=2,text="density")
    #Plot nice log-log axes
    full=c("1","10","100","1000","10000","100000","1000000","10000000","100000000")
    y1   <- range(x)
    pow  <- seq(y1[1], y1[2]+1)
    ticksaty <- as.vector(sapply(pow, function(p) (1:10)*10^p))
    laby <- replicate(length(ticksaty),"")
    idy = match(as.numeric(full),ticksaty)
    idy = idy[!is.na(idy)]
    laby[idy]<-full[c(1:length(match(as.numeric(full),ticksaty)))]
    axis(side=1, at=ticksaty, labels=TRUE, tcl=-0.25, lwd=0, lwd.ticks=0.5,cex=0.8)
}
dev.off()




# Step 3: Prepare GRanges object composed of 3' UTRs regions only (no gene body):
#   --	one SEGMENT can have as nearest neighbour a single 3'UTR but several transcripts
#		--> due to transcript isoforms that either share stop codons or part of the UTR
#		--> initially we will elongate every overlapped 3'UTRs
#	 	--> however the extensions that overlap 2 DIFFERENT stop codons will be removed in subsequent step
#
#   --	one TRANSCRIPT ID can have serveal nearest FRAGMENTS
#		--> in case of very large UTRs
#		--> in this case the fragment that contribute to longest UTRs is retained
print("Prepare GRanges object composed of UTRs regions only (no gene body)")
gF3               <- g3utr[idxG,]
POS               <- as.character(strand(gF3))=="+"
end(gF3)[POS]     <- end(gF2)[POS]
start(gF3)[!POS]  <- start(gF2)[!POS]
print("Filter n0.3 : only 3'UTR")
export.gff(object=gF3,con=paste("gF3_",name,".gff",sep=""),format="gff")

#
#
###



## Step 4 : Overlaps zero 5' exon boundaries
#
#    STEP 1 -- check for any overlap with 5'UTR
#    STEP 2 -- check for any overlap with any non-coding transcript


FilterExonBoundary <- function(gr){
  print("start filter no overlap with 5'exon boundaries")
  #Identify which covers a 5'UTR
  myStartCodons <- e75f[e75f$type=="start_codon"]
  gOver 		<- findOverlaps(query=gr,subject=myStartCodons,ignore.strand=FALSE)
  mygene1       <- gr$gene_id.[queryHits(gOver)]
  mygene2       <- myStartCodons$gene_id.[subjectHits(gOver)]
  F2a  		    <- unique(queryHits(gOver)[mygene2!=mygene1])

  #Identify which overlaps with any non-coding transcript
  myNonCoding <- e75p[e75p$source!="protein_coding"]
  gOver 	  <- findOverlaps(query=gr,subject=myNonCoding,ignore.strand=FALSE)
  mygene1     <- gr$gene_id.[queryHits(gOver)]
  mygene2     <- myNonCoding$gene_id.[subjectHits(gOver)]
  F2b  		  <- unique(queryHits(gOver)[mygene2!=mygene1])

  F2          <- unique(c(F2a,F2b))

  if(length(F2)!=0){
    print(paste("N0. of UTRs not passing the filter:", length(F2), sep=""))
    print("terminated filter no overlap with 5'exon boundaries")
    return(gr[-F2,])
  }
  else{
    print(paste("N0. of UTRs not passing the filter:", length(F2), sep=""))
    print("terminated filter no overlap with 5'exon boundaries")
    return(gr)
  }
}
print("Filter n0.2 : Overlaps zero 5' exon boundaries")
gF4 <- FilterExonBoundary(gF3)
export.gff(object=gF4,con=paste("gF4_",name,".gff",sep=""),format="gff")

#
###



# Filter n0.4 : Low-mappable region accounts for < 20% full extension length
#
#

FilterLowMappable <- function(gr){
    RpM     <- import.gff("/farm/home/luisie01/RepeatMask/rn5/rmsk_rn5_ucsc.gff")
    # As there are some overlapping regions,
    gRp <- reduce(RpM)
    print("start filter low mappable region")
    #test1       <- subsetByOverlaps(subject=gr,query=gRp,ignore.strand=FALSE)#returns the subset of ‘query’ that has an overlap hit with a range in ‘subject’ using the specified ‘findOverlaps’ parameters.
    gOver       <- findOverlaps(query=gr,subject=gRp,ignore.strand=FALSE)
    idxU        <- queryHits(gOver)
    idxG        <- subjectHits(gOver)

   if(length(idxU)>0){
        grouping    <- factor(idxU)
        Lev         <- as.numeric(levels(grouping))
        F5          <- replicate(length(Lev),FALSE)
        for(i in 1:length(Lev)){
            F5[i]<- sum(width(ranges(gRp[idxG[which(idxU==Lev[i])],])))/width(ranges(gr[Lev[i],]))>0.2
        }
        if(sum(F5)!=0){
            print(paste("N0. of UTRs not passing the filter:", sum(F5), sep=""))
            print("terminated filter low mappable region")
            return(gr[-Lev[F5],])
        }
        else{
        print(paste("N0. of UTRs not passing the filter:0", sep=""))
        print("terminated filter low mappable region")
        return(gr)
        }
   }
    else{
        print(paste("N0. of UTRs not passing the filter:0", sep=""))
        print("terminated filter low mappable region")
        return(gr)
    }
}

print("Filter n0.5 : Low-mappable region accounts for < 20% full extension length")
#gF5  <- FilterLowMappable(gF4)
gF5  <- gF4
export.gff(object=gF5,con=paste("gF5_",name,".gff",sep=""),format="gff")


#
#
###

# Filter n0.6:  Overlaps a single STOP codon
#
#
FilterSingleCodon <- function(gr){
   print("start filter single codon")
  gtemp                     <- e75p[which(e75p$type=="stop_codon"),]
  gtemp                     <- reduce(gtemp) #MERGING overlapping STOP codons as several transcripts may share STOP codons
  gtemp                     <- gtemp[match(unique(gtemp$gene_id.),gtemp$gene_id.),]
  POS                       <- as.character(strand(gr))=="+"
  start(ranges(gr))[POS]    <- start(ranges(gr))[POS]-1
  start(ranges(gr))[!POS]   <- start(ranges(gr))[!POS]+1
  gOver                     <- findOverlaps(query=gr,subject=gtemp,ignore.strand = FALSE)
  idxU                      <- queryHits(gOver)
  idxG                      <- subjectHits(gOver)
  grouping                  <- factor(idxU)
  F5                        <- unlist(tapply(idxG,grouping,
                                            function(x)return(length(x)>1)
                                            )
                                     )

  if(sum(F5)!=0){
    print(paste("N0. of UTRs not passing the filter:", sum(F5), sep=""))
    print("terminated filter single codon")
    return(gr[-as.numeric(names(F5))[F5],])
  }
  else{
    print(paste("N0. of UTRs not passing the filter:", sum(F5), sep=""))
    print("terminated filter single codon")
    return(gr)
  }
}
gF6 <- FilterSingleCodon(gF5)
print("Filter n0.6:  Overlaps a single STOP codon")
export.gff(object=gF6,con=paste("gF6_",name,".gff",sep=""),format="gff")

#
#
###


# Filter n0.7: Eliminate potential retained introns
#
#              Check for overlap with introns not contained in 3'UTR.
#

FilterIntrons <- function(gr){
  print("start filter introns")
  #Identify which overlaps with transcripts of any kind (BUT one extension can overlap with several transcript from same gene...)
  gtx         <- e75p[which(e75p$type%in%c("transcript")),]
  gexons      <- e75p[which(e75p$type%in%c("start_codon","stop_codon","exon","UTR","CDS")),]
  gin1        <- setdiff(gtx,gexons)
  ginf        <- setdiff(gin1,g3utr)

  #Add gene_id. to introns
  got                            <- findOverlaps(query=ginf,subject=gtx,ignore.strand=FALSE)
  ginf$gene_id.                  <- gtx$gene_id.[subjectHits(got)[match(c(1:length(ginf)),queryHits(got))]]

  # Then simply check that the UTR does not cover a single feature
  gOver                         <- findOverlaps(query=gr,subject=ginf,ignore.strand=FALSE)
  mygene1                       <- gr$gene_id.[queryHits(gOver)]
  mygene2                       <- ginf$gene_id.[subjectHits(gOver)]
  F7   						    <- unique(queryHits(gOver)[mygene2!=mygene1])


  if(length(F7)>=1){
    print(paste("N0. of UTRs not passing the filter F7:", length(F7), sep=""))
    print("terminated filter introns")
    return(gr[-F7,])
  }
  else{
    print(paste("N0. of UTRs not passing the filter:", length(F7), sep=""))
    print("terminated filter introns")
    return(gr)
  }

}
print("Filter n0.7: Eliminate potential retained introns")
gF7 <- FilterIntrons(gF6)
export.gff(object=gF7,con=paste("gF7_",name,".gff",sep=""),format="gff")


#
#
###

# Filter n0.8: Check that the last exon of the extension is also expressed, assuming that the shortest 3'UTR is allways expressed -- I do not agree with this given the current protocol; if 3'UTR super long, then we may never reach the end. check gr is output from 9-steps-....rat.r
#
#
#

LastCheck <- function(gr){
    print("start last check and merge")
    gOver <- findOverlaps(query=gr,subject=check,ignore.strand=FALSE,select="all")
    idxU  <- queryHits(gOver)
    idxG  <- subjectHits(gOver)
    #Check that there is no duplicated idxG as this would mean that there are more than one extension
    if(sum(duplicated(idxG))!=0){
        print("oups, more than one extension per tX id!")

        oi <- idxG[which(duplicated(idxG))]

        torm <- duplicated(idxG)
        idxU <- idxU[-torm]
        idxG <- idxG[-torm]
    }
    outgr <- check[idxG]
   #Now extend temp
    if(length(idxU)>=1){
        print(paste("N0. of UTRs not passing the last check:", length(gr)-length(unique(idxU)), sep=""))
        POS                              <- as.character(strand(outgr))=="+"
        end(ranges(outgr[POS]))          <- end(ranges(gr[idxU[POS]]))
        start(ranges(outgr[!POS]))       <- start(ranges(gr[idxU[!POS]]))
        return(outgr)
    }
    else{
        print(paste("N0. of UTRs not passing the filter:", length(idxU), sep=""))
        print("terminated filter introns")
        return(outgr)
    }

}
print("Filter n0.8: start last check and merge")
#gF8      <- LastCheck(gF7)
#export.gff(object=gF8,con=paste("gF8_",name,".gff",sep=""),format="gff")
#gF8   <- gF7
#gOver <- findOverlaps(query=gF8,subject=check,ignore.strand=FALSE,select="all")
#idxU  <- queryHits(gOver)
#idxG  <- subjectHits(gOver)
#toadd <- c(1:length(check))[-unique(idxG)]
#gF9   <- sort(c(gF8,check[toadd]))

gF9    <- sort(gF7)
export.gff(object=gF9,con=paste("gF9_",name,".gff",sep=""),format="gff")


#
#
###


print("Filtering terminated!")


