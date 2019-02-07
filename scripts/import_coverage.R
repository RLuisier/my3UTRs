library(grDevices)
library(Rsamtools)
library(IRanges)
library(GenomicRanges)
library(rtracklayer)

rm(list=ls())

### Nested Functions

#Selection on those which are separated by at least 500 nt starting from the most distant (should be improved; really not perfect)
SelectBasedOnDistance <- function(famtxID=anno_ngf[anno_ngf$txID=="ENSRNOT00000019181",]){
  famtxID <- famtxID[sort(famtxID$newL,index.return=T,decreasing=FALSE)$ix,]
  while(sum(diff(famtxID$newL)<300)>0){
    if(nrow(famtxID)==2){
      return(famtxID[2,])
    }
    if(nrow(famtxID)==3){
      famtxID <- famtxID[-2,]
    }
    if(nrow(famtxID)>3){
      famtxID <- famtxID[-(which(abs(diff(famtxID$newL))<300)+1)[1],]
    }
  }
  return(famtxID)
}

#Remove the duplicates
plotMclust <- function(mydata=mydat[[1]],mybim=bimdens[[1]],myLim=log2(Lim[1]),myname=names(Lim)[1],mytitle="TPM"){
  x  <- seq(from=0, to=floor(max(mydata)),length=100)
  if(length(mybim$parameters$variance$sigmasq)>1){
    hx1 <- dnorm(x,mean=mybim$parameters$mean[1],sd=sqrt(mybim$parameters$variance$sigmasq[1]))
    hx2 <- dnorm(x,mean=mybim$parameters$mean[2],sd=sqrt(mybim$parameters$variance$sigmasq[2]))
    hist(mydata,breaks=100,col=rgb(0,0,0,alpha=0.2),freq=FALSE,ylim=c(0,max(c(hx1,hx2))),xlab="",ylab="",main="",xaxt="n")
    lines(x,hx2 , lwd=2, col="blue")
    lines(x,hx1 , lwd=2, col="red")
    abline(v=myLim,lty=2,col="red")
    mtext(side=1,line=2,text=mytitle,cex=0.5)
    mtext(side=2,line=2,text="density",cex=0.5)
    mtext(side=3,line=1,text=myname,cex=0.5)
    #Plot nice log-log axes
    full=c("1","2","4","8","16","32","64","128","256")
    y1   <- floor(range(x))
    pow  <- seq(y1[1], y1[2]+1)
    ticksaty <- as.vector(sapply(pow, function(p) (1:10)*2^p))
    laby <- replicate(length(ticksaty),"")
    idy = match(as.numeric(full),ticksaty)
    idy = idy[!is.na(idy)]
    laby[idy]<-full[c(1:length(match(as.numeric(full),ticksaty)))]
    axis(side=1, at=ticksaty, labels=TRUE, tcl=-0.25, lwd=0, lwd.ticks=0.5,cex=0.8)
  }
  if(length(mybim$parameters$variance$sigmasq)==1){
    hx1 <- dnorm(x,mean=mybim$parameters$mean[1],sd=sqrt(mybim$parameters$variance$sigmasq[1]))
    hx2 <- dnorm(x,mean=mybim$parameters$mean[2],sd=sqrt(mybim$parameters$variance$sigmasq[1]))
    hist(mydata,breaks=100,col=rgb(0,0,0,alpha=0.2),freq=FALSE,ylim=c(0,max(c(hx1,hx2))),xlab="",ylab="",main="",xaxt="n")
    lines(x,hx2 , lwd=2, col="blue")
    lines(x,hx1 , lwd=2, col="red")
    abline(v=myLim,lty=2,col="red")
    mtext(side=1,line=2,text=mytitle,cex=0.7)
    mtext(side=2,line=2,text="density",cex=0.6)
    mtext(side=3,line=1,text=myname,cex=0.5)
    #Plot nice log-log axes
    full=c("1","2","4","8","16","32","64","128","256")
    y1   <- floor(range(x))
    pow  <- seq(y1[1], y1[2]+1)
    ticksaty <- as.vector(sapply(pow, function(p) (1:10)*2^p))
    laby <- replicate(length(ticksaty),"")
    idy = match(as.numeric(full),ticksaty)
    idy = idy[!is.na(idy)]
    laby[idy]<-full[c(1:length(match(as.numeric(full),ticksaty)))]
    axis(side=1, at=ticksaty, labels=TRUE, tcl=-0.25, lwd=0, lwd.ticks=0.5,cex=0.8)
  }
}


RemoveIds <- function(gr){
  
  uniqs   <- paste(as.character(seqnames(gr)),as.character(start(gr)),as.character(end(gr)),as.character(strand(gr)),sep=".")
  freq    <- as.data.frame(table(uniqs))
  todeal  <- unique(as.character(freq[freq[,2]>1,1]))
  
  
  #On garde soit le premier ...
  sel1    <- uniqs%in%todeal&as.character(gr$isoform)==0
  sel1n   <- uniqs%in%unique(uniqs[sel1])&as.character(gr$isoform)!=0
  
  #... soit un autre si le premier n'est pas isoform 0.
  todeal1 <- setdiff(todeal,uniqs[sel1])
  sel2    <- !c(1:length(uniqs))%in%which((uniqs%in%todeal1)&!c(1:length(uniqs))%in%match(todeal1,uniqs))
  
  gr      <- gr[sel1|sel1n|sel2,]
  
  return(gr)
  
}

DeltaLp <- function(gr){
  ix       <- match(as.character(gr$txID),as.character(mcols(g3utr)$X.transcript_id.))
  Lnew     <- width(gr)
  goi      <- g3utr[ix]
  
  # Remove introns from length of init 3'UTR
  over                         <-  as.matrix(findOverlaps(query=g3introns,subject=goi,ignore.strand=FALSE))
  torm                         <-  tapply(width(pintersect(g3introns[over[,1]],goi[over[,2]])),INDEX=as.factor(over[,2]),FUN=sum)
  Lutr                         <-  width(goi)
  Lutr[as.numeric(names(torm))]<-  Lutr[as.numeric(names(torm))]-torm
  
  
  # Remove introns from length of final 3'UTR
  over                         <-  as.matrix(findOverlaps(query=g3introns,subject=gr,ignore.strand=FALSE))
  torm                         <-  tapply(width(pintersect(g3introns[over[,1]],gr[over[,2]])),INDEX=as.factor(over[,2]),FUN=sum)
  Lnew                         <-  width(gr)
  Lnew[as.numeric(names(torm))]<-  Lnew[as.numeric(names(torm))]-torm
  
  # Compute difference in length
  dL                           <-  Lnew-Lutr
  return(dL)
}


ComputeWidth <- function(gr){
  
  Lnew                         <- width(gr)
  # Remove introns from length of final 3'UTR
  over                         <-  as.matrix(findOverlaps(query=g3introns,subject=gr,ignore.strand=FALSE))
  torm                         <-  tapply(width(pintersect(g3introns[over[,1]],gr[over[,2]])),INDEX=as.factor(over[,2]),FUN=sum)
  #?keep when the end of the 3'utr is in the middle of the introns
  
  # Correct for introns
  Lnew[as.numeric(names(torm))]<-  Lnew[as.numeric(names(torm))]-torm
  
  return(Lnew)
}

ComputeWidthTx <- function(gr){
  
  Lnew     <- width(gr)
  # Remove introns from length
  over                         <-  as.matrix(findOverlaps(query=gintrons,subject=gr,ignore.strand=FALSE))
  torm                         <-  tapply(width(pintersect(gintrons[over[,1]],gr[over[,2]])),INDEX=as.factor(over[,2]),FUN=sum)
  Lnew                         <-  width(gr)
  Lnew[as.numeric(names(torm))]<-  Lnew[as.numeric(names(torm))]-torm
  
  return(Lnew)
}
#
#
###


# A. Load data
load("./annotation/rn5/GRanges_comprehensive_transcriptome_rat_24_nov_2015.RData")
ngfGRS         <- import.gff("./utrid/APA/L2.gtf",format="gtf")
covdir         <- "./Coverage/utrCov"

# B. Import coverage
foi                  <- list.files(paste(covdir,"/500",sep=""))
names                <- gsub(foi,pattern=".txt",repl="")
types                <- do.call(lapply(names,function(x)return(length(unlist(strsplit(x,split="[\\_,\\.,-]"))))),what=c)
names_raw           <- names[intersect(which(types=="4"),grep(names,pattern=".chr"))]
names_raw           <- do.call(lapply(names[intersect(which(types=="4"),grep(names,pattern=".chr"))],function(x)return(unlist(strsplit(x,split="[\\_,\\.,-]")))),what=rbind)
names_raw           <- data.frame(foi[intersect(which(types=="4"),grep(names,pattern=".chr"))],names_raw)
colnames(names_raw )<-c("id","treatment","origin","replicate","chr")
samples                 <- factor(paste(names_raw$treatment,names_raw$origin,names_raw$replicate,sep="."))
covraw_500              <- lapply(c(1:length(unique(samples))),function(X)return(do.call(lapply(which(samples==unique(samples)[X]), function(Z)return(read.table(paste(covdir,"/500/",names_raw$id[Z],sep="")))),what=rbind)))
rown500                 <-  unique(do.call(lapply(covraw_500,function(x)return(as.character(x[,1]))),what=c))
myCovRaw500             <- matrix(0,nrow=length(rown500),ncol=length(covraw_500))
for(i in c(1:length(covraw_500))){
  myCovRaw500[,i]       <- covraw_500[[i]][match(rown500,as.character(covraw_500[[i]][,1])),2]
}
colnames(myCovRaw500)   <- unlist(lapply(unique(samples),function(x)return(paste(x,".raw",sep=""))))
rownames(myCovRaw500)   <- rown500
myCov500                <- myCovRaw500

write.table(x=myCov500,file="./data/myCov500.tab",row.names = TRUE,col.names = TRUE,quote=FALSE,sep="\t")


# C. Select only those transcripts covered by at least 10 reads
dat1               <- myCov500[grep(rownames(myCov500),pattern="\\.0"),]
mylistCov500       <- lapply(c(1:4),function(x)return(log2(dat1[dat1[,x]>=1,x])))
names(mylistCov500)<-colnames(myCov500)
mydat              <- mylistCov500
require(mclust)
bimdens <- list()
for(i in c(1:length(mydat))){
  bimdens[[i]]<- densityMclust(data=mydat[[i]],G=2)
}
names(bimdens)<-names(mydat)
Lim         <- rep(10,4)#Please note that this limits will not be used at a later stage; only for this initial
Limp        <- log2(Lim)

par(mfrow=c(2,2),mar=c(3,3,3,3))
for(i in c(1:length(Lim))){
  no.expressed <- sum(mydat[[i]]>=Limp[i])
  plotMclust(mydata=mydat[[i]],mybim=bimdens[[i]],myLim=Limp[i],myname=names(bimdens)[i],mytitle="Coverage")
    mtext(side=3,line=0,text=paste("no.expressed=",no.expressed),cex=0.6)

}

mydata  <- log2(myCov500)
mydata[is.na(mydata)]<-0
myselG  <- matrix(FALSE,nrow=nrow(mydata),ncol=ncol(mydata))
for(i in c(1:ncol(myselG))){
  myselG[,i] <- mydata[,i]>=Limp[i]
}
colnames(myselG) <- colnames(mydata)
rownames(myselG) <- rownames(mydata)

info                 <- do.call(lapply(colnames(myCov500),function(x)return(unlist(strsplit(x,split="\\.")))),what=rbind)
info                 <- data.frame(colnames(myCov500),info)
colnames(info)       <-c("id","treatment","origin","replicate","type")
is.expressed.iso            <- t(apply(myselG,1,function(x)return(tapply(x,INDEX=factor(info$origin),FUN=function(z)return(sum(z)==2)))))

#Selection based on the expression and I0
selNGF <- union(which(ngfGRS$ID%in%rownames(is.expressed.iso)[apply(is.expressed.iso,1,function(x)return(sum(x)>0))]),
            grep("\\.0",ngfGRS$ID))
Lngf                        <- ngfGRS[selNGF,]
export.gff(Lngf,con="./annotation/rn5/Lngf.gtf",format="gtf")
Lngf      <- RemoveIds(Lngf)


# F. Create matrix containing annotation of the isoforms together with coverage -- anno_ngf
myOut           <- myCov500[match(Lngf$ID,rownames(myCov500)),]
ngfGRS          <- Lngf
#I should remove from here all the isoforms for which utrL<=10
tokeep   <- width(ngfGRS)>10
ngfGRS   <- ngfGRS[tokeep,]
myOut    <- myOut[tokeep,]

#Is conservative
is.conservative           <- ngfGRS$isoform==0
#initL
Linit        <- ComputeWidth(g3utr)
ids          <- ngfGRS$ID
txID         <- ngfGRS$txID#  22'786 unique tX IDs for 43'203 isoforms
iso          <- ngfGRS$isoform
initL        <- Linit[match(txID,g3utr$X.transcript_id.)]
newL         <- ComputeWidth(ngfGRS)
GS           <- myOut$GS
#txLength
alltx               <- e75f[which(mcols(e75f)$type=="transcript")]
ix                  <- match(txID,as.character(mcols(alltx)$X.transcript_id.))
alltx               <- alltx[ix,]
gtx                 <- e75f[which(e75f$type%in%c("transcript")),]
gid                 <- e75f[which(e75f$type%in%c("gene")),]
gexons              <- e75f[which(e75f$type%in%c("start_codon","stop_codon","exon","UTR","CDS")),]
gintrons            <- setdiff(gtx,gexons)
gintergenic         <- setdiff(gtx,gexons)
txLength            <- ComputeWidthTx(alltx)
#maxL
tempMax      <- tapply(newL,INDEX=factor(txID),max)
maxL         <- tempMax[match(txID,names(tempMax))]
#minL
tempMin      <- tapply(newL,INDEX=factor(txID),min)
minL         <- tempMin[match(txID,names(tempMin))]
#GeneSymbol
geneSymbol <- as.character(mcols(e75f)$X.gene_name.)[match(txID,as.character(mcols(e75f)$X.transcript_id.))]
#
anno_ngf     <- data.frame(uniqueID=ids,txID=txID,isoform=iso,txLength,is.conservative,initL,newL,maxL,minL,geneSymbol,no.iso,myOut)
#no.iso
temp                      <- as.data.frame(table(as.character(anno_ngf$txID)))
anno_ngf$no.iso           <- temp[match(anno_ngf$txID,as.character(temp[,1])),2]

tomodif <- anno_ngf$no.iso>1
sub1    <- anno_ngf[!tomodif,]
sub2    <- anno_ngf[tomodif,]
Lev     <- unique(sub2$txID)
modifed <- SelectBasedOnDistance(famtxID=sub2[sub2$txID==Lev[1],])
for(i in c(2:length(Lev))){
    modifed <- rbind(modifed,SelectBasedOnDistance(famtxID=sub2[sub2$txID==Lev[i],]))
    print(i)

}
selected                        <- c(as.character(sub1$uniqueID),as.character(modifed$uniqueID))
anno_ngf$selected.for.distance  <- anno_ngf$uniqueID%in%selected


write.table(x=anno_ngf,file="./data/anno_ngf.tab",row.names = TRUE,col.names = TRUE,quote=FALSE,sep="\t")
ngfGRS <- ngfGRS[which(ngfGRS$uniqueID%in%anno_ngf$uniqueID[anno_ngf$selected.for.distance]),]
export.gff(object=ngfGRS,con="./annotation/rn5/Lngf_sub.gtf",format="gtf")






