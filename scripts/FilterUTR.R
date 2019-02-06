library(grDevices)
library(Rsamtools)
library(IRanges)
library(GenomicRanges)
library(rtracklayer)
library(gplots)
require(EBS)
require(Segmentor3IsBack)


#
# 0. Import key files
#
# 1/ Import init annotated 3'UTR
load("./annotation/rn5/GRanges_comprehensive_transcriptome_rat_24_nov_2015.RData")
myPAS <- import.gff("./annotation/rn5/myPAS.gtf",format="gtf",asRangedData = FALSE)

#
# 1. Create myGR_L1
#
foi<- list.files("./utrid/APA")
foi <- foi[grep(foi,pattern=".RData")]
load(paste("/home/rluisier/data/Riccio/Exp_1/Dec2016/utrid/APA/",foi[1],sep=""))


mynewGR  <-myGRobject
start   <- start(myGRobject)
end     <- end(myGRobject)
chr     <- as.character(seqnames(myGRobject))
strand  <- as.character(strand(myGRobject))
txID    <- myGRobject$txID
isoform <- myGRobject$isoform

for(i in c(2:length(foi))){
  load(paste("./utrid/APA/",foi[i],sep=""))
  start  <- c(start,start(myGRobject))
  end    <- c(end,end(myGRobject))
  chr    <- c(chr,as.character(seqnames(myGRobject)))
  strand <- c(strand,as.character(strand(myGRobject)))
  txID   <- c(txID,myGRobject$txID)
  isoform<- c(isoform,myGRobject$isoform)
}

start       <- c(start(g3utr),start)
end         <- c(end(g3utr),end)
chr         <- c(as.character(seqnames(g3utr)),chr)
strand      <- c(as.character(strand(g3utr)),strand)
txID        <- c(as.character(mcols(g3utr)$X.transcript_id.),txID)
isoform     <- c(rep(0,length(g3utr)),isoform)

mynewGR    <-  GRanges(seqnames = Rle(chr),
                       ranges   = IRanges(start=start,end=end),
                       strand   = Rle(strand),
                       data.frame(txID=txID,isoform=isoform)
)


txID            <- mynewGR$txID
myIso           <- mynewGR$isoform
uniqueID        <- paste(txID,myIso,sep=".")
names(mynewGR)  <- uniqueID
dups            <- duplicated(paste(mynewGR$txID,start(mynewGR),end(mynewGR),as.character(strand(mynewGR)),sep="."))
mynewGR         <- mynewGR[!dups,]
txID            <- mynewGR$txID

myGR            <- mynewGR
ids            <- paste(as.character(seqnames(myGR)),as.character(start(myGR)),as.character(end(myGR)),as.character(strand(myGR)),sep=".")
no.dups        <- as.data.frame(table(ids))
ids.dups       <- unique(no.dups$ids[no.dups$Freq>1])
i1             <- match(ids.dups,ids)
tocor          <- ids[i1[myGR$isoform[i1]==0]]
i2             <- which(ids%in%tocor&myGR$isoform!=0)
torm           <- unique(c(i2,i1[myGR$isoform[i1]!=0]))
myGR           <- myGR[-torm,]
myGR$ID        <- names(myGR)


ComputeWidth <- function(gr){

  Lnew                         <- width(gr)
  # Remove introns from length of final 3'UTR
  over                         <-  as.matrix(findOverlaps(query=g3introns,subject=gr,ignore.strand=FALSE))
  torm                         <-  tapply(width(pintersect(g3introns[over[,1]],gr[over[,2]])),INDEX=as.factor(over[,2]),FUN=sum)
  # Correct for introns
  Lnew[as.numeric(names(torm))]<-  Lnew[as.numeric(names(torm))]-torm

  return(Lnew)
}

RemoveCloseNeighbour <- function(gr){

  txID      <- factor(as.character(gr$txID))
  newL      <- ComputeWidth(gr)
  temp      <- as.data.frame(table(txID))
  no.iso    <- temp[match(as.character(txID),as.character(temp[,1])),2]

  #iso_ordered
  tempL                     <- newL
  uniqueID                  <- paste(gr$txID,gr$isoform,sep=".")
  names(tempL)              <- uniqueID
  gr$uniqueID               <- uniqueID
  tempG                     <- as.factor(as.character(gr$txID))
  myord                     <- tapply(tempL,INDEX=tempG,function(x)return(cbind(names(x)[sort(as.numeric(as.character(x)),index.return=T,decreasing=F)$ix],c(1:length(x)))))
  test                      <- do.call(what=rbind,args=myord)
  test                      <- test[match(gr$uniqueID,test[,1]),]
  gr$iso_ordered            <- test[,2]
  #no.iso
  temp                      <- as.data.frame(table(as.character(gr$txID)))
  temp                      <- temp[match(gr$txID,as.character(temp[,1])),2]
  gr$no.iso                 <- temp
  #priorIsoform
  tempid                     <- paste(gr$txID,gr$iso_ordered,sep=".")
  prior.iso                  <- paste(gr$txID,(as.numeric(gr$iso_ordered)-1),sep=".")
  priorID                    <- gr$uniqueID[match(prior.iso,tempid)]
  gr$priorID                 <- priorID
  #distToprior
  distoprior                            <- rep(NA,length(gr))
  distoprior[!is.na(gr$priorID)]        <- as.numeric(as.character(newL))[!is.na(gr$priorID)] - as.numeric(as.character(newL))[match(gr$priorID[!is.na(gr$priorID)],gr$uniqueID)]
  gr$distoprior                         <-distoprior

  #distTonext
  tempid                               <- paste(gr$txID,gr$iso_ordered,sep=".")
  next.iso                             <- paste(gr$txID,(as.numeric(gr$iso_ordered)+1),sep=".")
  nextID                               <- gr$uniqueID[match(next.iso,tempid)]
  gr$nextID                            <- nextID
  distonext                            <- rep(NA,length(gr))
  sel                                  <- !is.na(gr$nextID)
  d1                                   <- as.numeric(as.character(newL))[sel]
  d2                                   <- as.numeric(as.character(newL))[match(gr$nextID[sel],gr$uniqueID)]
  distonext[sel]                       <- d2-d1
  gr$distonext                         <- distonext

  # Recursively remove closest when whithin 50 nt
  oi             <- gr$distonext<=50
  oi[is.na(oi)]  <- FALSE
  toi            <- unique(gr$txID[oi])
  if(length(toi)>0){
    torm1          <- which(oi&as.character(gr$isoform)!="0")
    torm2          <- match(gr$nextID[oi&as.character(gr$isoform)=="0"],gr$ID)
    torm           <- c(torm1,torm2)
    out             <- gr[-torm,]
    return(out)
  }
  else{
    print("no remaining close neighbours")
    return(gr)
  }

}



test  <- RemoveCloseNeighbour(myGR)
test2 <- RemoveCloseNeighbour(test)
test3 <- RemoveCloseNeighbour(test2)


export.gff(object=test2,con="./utrid/APA/L1.gtf",format="gtf")
myGR  <- test2


#
# 2. Remove if amount of RpM is higher than 30% -- I should here not ignore the strand
#
RpM                     <- import.gff("./annotation/rn5/rmsk_rn5_ucsc.gff",asRangedData = FALSE)
myRegion                <- myGR
POS                     <- as.character(strand(myRegion))=="+"
start(myRegion)[POS]    <- end(myRegion)[POS]-499
end(myRegion)[!POS]     <- start(myRegion)[!POS]+499
gOver                   <- findOverlaps(query=myRegion,subject=RpM,ignore.strand=TRUE)
idxU                    <- queryHits(gOver)
idxG                    <- subjectHits(gOver)

myStart                 <- apply(cbind(start(myRegion[idxU]),start(RpM[idxG])),1,max)
myEnd                   <- apply(cbind(end(myRegion[idxU]),end(RpM[idxG])),1,min)
mywid                   <- myEnd-myStart+1
myfracRpM               <- rep(0,length(myRegion))
myfracRpM[unique(idxU)] <- tapply(mywid,INDEX=idxU,FUN=sum)/500

torm                    <- myfracRpM>=0.3&as.character(myGR$isoform)!="0"
myGR                    <- myGR[!torm,]



#
# 3. Filter when overlap with antisense transcripts
#
myTX                <- e75f[which(e75f$type=="transcript"),]
POS                 <- as.character(strand(myTX))=="+"
strand(myTX)[POS]   <- "-"
strand(myTX)[!POS]  <- "+"

# B. Check for no overlap with the new isoforms (however this should not remove those which already overlap)
gOver 		<- findOverlaps(query=myGR,subject=myTX,ignore.strand = FALSE)
F1   		<- c(1:length(myGR))%in%unique(queryHits(gOver))#Check whether the new isoforms overlap with antisense transcripts; I may need to check whether the antisense transcritp is expressed?
corr_init   <- g3utr[match(myGR$txID,g3utr$X.transcript_id.),]#Check whether the intial I0 already overlaps with antisense transcript
gOver 		<- findOverlaps(query=corr_init,subject=myTX,ignore.strand = FALSE)
#Check that they both do not overlap with the same
F2          <- !c(1:length(myGR))%in%unique(queryHits(gOver))
F3          <- !c(1:length(myGR))%in%grep(myGR$uniqueID,pattern="\\.0")
L2          <- myGR[!(F1&F2&F3),]
myGR        <- L2



#
# 3. Add flags:
# 3.1 PAS which do not contain PAS within ]-50; +50[; keep always Iso.0
#

#Focus on the regions -100;+50 around the PAS (Gruber et al, positioning)
POS                             <- as.character(strand(myGR))=="+"
subGR                           <- myGR
start(subGR)[POS]               <- end(myGR[POS])-100
end(subGR)[POS]                 <- end(myGR[POS])+50
end(subGR)[!POS]                <- start(myGR[!POS])+100
start(subGR)[!POS]              <- start(myGR[!POS])-50
start(subGR)[start(subGR)<0]    <- 1
temp                            <- as.data.frame(queryHits(findOverlaps(query=subGR,subject=myPAS,ignore.strand=FALSE)))
is.motif                        <- rep(FALSE,length(subGR))
is.motif[unique(temp[,1])]      <- TRUE
myGR$is.pas                     <- is.motif


#
# 3. Add flags:
# 3.2 Contains polydT at the end -- check for internal priming
#
dT7 <- import.gff("./annotation/rn5/dT7.gtf",format="gtf",asRangedData=FALSE)#Flase positive

#Focus on the regions -10;+50 around the PAS
POS                             <- as.character(strand(myGR))=="+"
subGR                           <- myGR
start(subGR)[POS]               <- end(myGR[POS])-10
end(subGR)[POS]                 <- end(myGR[POS])+50
end(subGR)[!POS]                <- start(myGR[!POS])+10
start(subGR)[!POS]              <- start(myGR[!POS])-50
start(subGR)[start(subGR)<0]    <- 1
temp                            <- as.data.frame(queryHits(findOverlaps(query=subGR,subject=dT7,ignore.strand=FALSE)))
is.motif                        <- rep(FALSE,length(subGR))
is.motif[unique(temp[,1])]      <- TRUE
myGR$internal.priming           <- is.motif



#
# 3. Add flags:
# 3.3 Potential false positive due to internal priming on the other strand
#
dU7 <- import.gff("./annotation/rn5/dU7.gtf",format="gtf",asRangedData=FALSE)#Flase positive

#Focus on the regions -200;-150 around the PAS
POS                             <- as.character(strand(myGR))=="+"
subGR                           <- myGR
start(subGR)[POS]               <- end(myGR[POS])-200
end(subGR)[POS]                 <- end(myGR[POS])-100
end(subGR)[!POS]                <- start(myGR[!POS])+200
start(subGR)[!POS]              <- start(myGR[!POS])+100
start(subGR)[start(subGR)<0]    <- 1
temp                            <- as.data.frame(queryHits(findOverlaps(query=subGR,subject=dU7,ignore.strand=FALSE)))
is.motif                        <- rep(FALSE,length(subGR))
is.motif[unique(temp[,1])]      <- TRUE
myGR$internal.priming.reverse   <- is.motif

export.gff(object=myGR,con="./utrid/APA/L2.gtf",format="gtf")


