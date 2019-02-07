library(grDevices)
library(Rsamtools)
library(IRanges)
library(GenomicRanges)
library(rtracklayer)

rm(list=ls())

# A. Load data
load("./annotation/rn5/GRanges_comprehensive_transcriptome_rat_24_nov_2015.RData")
ngfGRS         <- import.gff("./utrid/APA/L2.gtf",format="gtf")
covdir         <- "./Coverage/utrCov"

# B. Import coverage
foi                  <- list.files(paste(covdir,"/500",sep=""))
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
covraw_500              <- lapply(c(1:length(unique(samples))),function(X)return(do.call(lapply(which(samples==unique(samples)[X]), function(Z)return(read.table(paste(covdir,"/500/",names_raw$id[Z],sep="")))),what=rbind)))
rown500                 <-  unique(do.call(lapply(covraw_500,function(x)return(as.character(x[,1]))),what=c))
myCovRaw500             <- matrix(0,nrow=length(rown500),ncol=length(covraw_500))
for(i in c(1:length(covraw_500))){
  myCovRaw500[,i]       <- covraw_500[[i]][match(rown500,as.character(covraw_500[[i]][,1])),2]
}
colnames(myCovRaw500)   <- unlist(lapply(unique(samples),function(x)return(paste(x,".raw",sep=""))))
rownames(myCovRaw500)   <- rown500

samples                 <- factor(paste(names_norms$treatment,names_norms$origin,names_norms$replicate,sep="."))
covnorm_500             <- lapply(c(1:length(unique(samples))),function(X)return(do.call(lapply(which(samples==unique(samples)[X]), function(Z)return(read.table(paste(covdir,"/500/",names_norms$id[Z],sep="")))),what=rbind)))
rown500                 <-  unique(do.call(lapply(covraw_500,function(x)return(as.character(x[,1]))),what=c))
myCovNorm500            <- matrix(0,nrow=length(rown500),ncol=length(covnorm_500))
for(i in c(1:length(covnorm_500))){
  myCovNorm500[,i]      <- covnorm_500[[i]][match(rown500,as.character(covnorm_500[[i]][,1])),2]
}
colnames(myCovNorm500)  <- unlist(lapply(unique(samples),function(x)return(paste(x,".norm",sep=""))))
rownames(myCovNorm500)  <- rown500

myCov500                <- myCovRaw500
myCov500                <- cbind(myCov500,myCovNorm500[match(rownames(myCov500),rownames(myCovNorm500)),])
id                      <- match(as.character(ngfGRS$uniqueID),rownames(myCov500))
ngfGRS                  <- ngfGRS[!is.na(id),]
id                      <- match(as.character(ngfGRS$uniqueID),rownames(myCov500))
myCov500                <- myCov500[id,]

write.csv(x=myCov500,file="./data/myCov500.csv")





# D. Identify expressed genes according to coverage: 500 nt. Please note that it seems more accurate to use the 22'000 tX rather than the entire pool to determine the threshold.
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


dat1  <- myCov500[grep(rownames(myCov500),pattern="\\.0"),c(1:8)]



mylistCov500  <- lapply(c(1:8),function(x)return(log2(dat1[dat1[,x]>=5.012531e-02,x])))

mylistCov500  <- lapply(c(1:8),function(x)return(log2(dat1[dat1[,x]>=1,x])))
mylistCov300  <- lapply(c(1:8),function(x)return(log2(dat2[dat2[,x]>=1,x])))

names(mylistCov500)<-colnames(myCov500)[c(1:8)]
names(mylistCov300)<-colnames(myCov300)[c(1:8)]


mydat <- mylistCov500
require(mclust)
bimdens <- list()
for(i in c(1:length(mydat))){
  bimdens[[i]]<- densityMclust(data=mydat[[i]],G=2)
}
names(bimdens)<-names(mydat)

#Select expression based on Mclust clustering
SelectExpressed1 <- function(dat=log2(htseq[,1]+1)){
    mod1   <- Mclust(dat,G=2)
    mod1dr <- MclustDR(mod1)
    #mod1$classification
    is.expressed<-mod1dr$class==2
    lim=min(dat[is.expressed])
    return(list(is.expressed,lim))
}

Lim <- unlist(lapply(c(1:length(mydat)),FUN=function(x)return(2^SelectExpressed1(dat=mydat[[x]])[[2]])))
Limp<- log2(Lim)

pdf(paste(outdir,"detect_expressed_genes_cov_500_all.pdf",sep=""))
par(mfrow=c(2,2),mar=c(3,3,3,3))
for(i in c(1:length(Lim))){
  no.expressed <- sum(mydat[[i]]>=Limp[i])
  plotMclust(mydata=mydat[[i]],mybim=bimdens[[i]],myLim=log2(Lim[i]),myname=names(bimdens)[i],mytitle="Coverage")
  mtext(side=3,line=0,text=paste("no.expressed=",no.expressed),cex=0.6)
}

Lim         <- rep(c(10,10,10,10),2)#Please note that this limits will not be used at a later stage; only for this initial
Limp        <- log2(Lim)

par(mfrow=c(2,2),mar=c(3,3,3,3))
for(i in c(1:length(Lim))){
  no.expressed <- sum(mydat[[i]]>=Limp[i])
  plotMclust(mydata=mydat[[i]],mybim=bimdens[[i]],myLim=Limp[i],myname=names(bimdens)[i],mytitle="Coverage")
    mtext(side=3,line=0,text=paste("no.expressed=",no.expressed),cex=0.6)

}
dev.off()



mydata  <- log2(myCov500[,which(info$type=="raw")])
mydata[is.na(mydata)]<-0
myselG  <- matrix(FALSE,nrow=nrow(mydata),ncol=ncol(mydata))
for(i in c(1:ncol(myselG))){
  myselG[,i] <- mydata[,i]>=Limp[i]
}
colnames(myselG) <- colnames(mydata)
rownames(myselG) <- rownames(mydata)

groups                      <- factor(paste(info$treatment[info$type=="raw"],info$origin[info$type=="raw"],sep="."))
is.expressed.iso            <- t(apply(myselG,1,function(x)return(tapply(x,INDEX=groups,FUN=function(z)return(sum(z)==2)))))

#Selection based on the expression AND on I0
selNGF <- union(which(ngfGRS$ID%in%rownames(is.expressed.iso)[apply(is.expressed.iso[,c(1,2)],1,function(x)return(sum(x)>0))]),
            grep("\\.0",ngfGRS$ID))




#I should indeed keep also the isoform I0 for downstream analysis
Lngf                        <- ngfGRS[selNGF,]
export.gff(Lngf,con="/home/rluisier/data/Riccio/Exp_1/Dec2016/utrid/APA_stringent/Lngf.gtf",format="gtf")#50'080


#Remove the duplicates
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

Lngf <- RemoveIds(Lngf)

export.gff(Lngf,con="/home/rluisier/data/Riccio/Exp_1/Dec2016/utrid/APA_stringent/Lngfp.gtf",format="gtf")#PLease note that this version does not contain those 3' end that do not overlap with a polyA site however without duplicates

#is.expressed.iso --> must be expressed in both samples
allGS                       <- as.character(mcols(e75p)$X.gene_name.)
allTX                       <- as.character(mcols(e75p)$X.transcript_id.)
myselG                      <- myselG[match(Ltot$uniqueID,rownames(myselG)),]
myCov500                    <- myCov500[match(Ltot$uniqueID,rownames(myCov500)),]
is.expressed.iso            <- is.expressed.iso[match(Ltot$uniqueID,rownames(is.expressed.iso)),]

groups                      <- factor(paste(info$treatment[info$type=="raw"],info$origin[info$type=="raw"],sep="."))
colnames(myselG)            <- paste(unlist(lapply(colnames(myselG),function(x)return(gsub(x,pattern=".raw",repl="")))),"isexpressed",sep=".")
colnames(is.expressed.iso)  <- unlist(lapply(colnames(is.expressed.iso),function(x)paste(x,"is.expressed.iso",sep=".")))

txID                        <- unlist(lapply(rownames(myselG),function(x)return(unlist(strsplit(x,split="\\."))[1])))
GS                          <- allGS[match(txID,allTX)]
temp                        <- apply(is.expressed.iso,2,function(X)return(tapply(X,INDEX=factor(txID),FUN=function(z)return(sum(z)>0))))
is.expressed.txID           <- temp[match(txID,rownames(temp)),]
temp                        <- apply(is.expressed.iso,2,function(X)return(tapply(X,INDEX=factor(GS),FUN=function(z)return(sum(z)>0))))
is.expressed.GS             <- temp[match(GS,rownames(temp)),]

colnames(is.expressed.txID) <- unlist(lapply(colnames(is.expressed.txID),function(x)return(gsub(x,pattern="iso",repl="txID"))))
colnames(is.expressed.GS)   <- unlist(lapply(colnames(is.expressed.GS),function(x)return(gsub(x,pattern="iso",repl="GS"))))

rownames(is.expressed.txID) <- rownames(myselG)
rownames(is.expressed.iso)  <- rownames(myselG)
rownames(is.expressed.GS)   <- rownames(myselG)

is.expressed.iso            <- data.frame(is.expressed.iso)
is.expressed.txID           <- data.frame(is.expressed.txID)
is.expressed.GS             <- data.frame(is.expressed.GS)




# F. Create single matrix
myOut   <- data.frame(tot.detect,myCov500,myselG)


load("/home/rluisier/data/Riccio/Exp_1/Dec2016/Coverage/utrCov_stringent/myCov_March_12_2016_500.RData")
myOut   <- myOut[myOut$detect.iso.ngf!="no"|c(1:nrow(myOut))%in%grep(rownames(myOut),pattern="\\.0"),]
myOut   <- myOut[rownames(myOut)%in%ngfGRS$uniqueID,]
ngfGRS  <- ngfGRS[match(rownames(myOut),ngfGRS$ID),]
outdir <- "/home/rluisier/data/Riccio/Exp_1/Dec2016/utrid/anno_ngf_stringent/"


ngfGRS         <- import.gff("/home/rluisier/data/Riccio/Exp_1/Dec2016/utrid/APA_stringent/Lngf.gtf",format="gtf",asRangedData=FALSE)#Takes the one with the replicates
load("/home/rluisier/data/Riccio/Exp_1/Dec2016/Coverage/utrCov_stringent/myCov_March_12_2016_500.RData")
myOut   <- myOut[myOut$detect.iso.ngf!="no"|c(1:nrow(myOut))%in%grep(rownames(myOut),pattern="\\.0"),]
myOut   <- myOut[rownames(myOut)%in%ngfGRS$uniqueID,]
ngfGRS  <- ngfGRS[match(rownames(myOut),ngfGRS$ID),]
outdir <- "/home/rluisier/data/Riccio/Exp_1/Dec2016/utrid/anno_ngf_stringent/"

#I should remove from here all the isoforms for which utrL<=10
torm   <- which(width(ngfGRS)<=10)
ngfGRS <- ngfGRS[-torm,]
myOut  <- myOut[-torm,]
#
# 1. Create basis of anno_ngf
#

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
ix                  <- match(txID,as.character(mcols(alltx)$X.transcript_id))
alltx               <- alltx[ix,]
gtx                 <- e75f[which(e75f$type%in%c("transcript")),]
gid                 <- e75f[which(e75f$type%in%c("gene")),]
gexons              <- e75f[which(e75f$type%in%c("start_codon","stop_codon","exon","UTR","CDS")),]
gintrons            <- setdiff(gtx,gexons)
gintergenic         <- setdiff(gtx,gexons)
txLength            <- ComputeWidthTx(alltx)
BINS                 <- cut(x=txLength,breaks=10^3*c(0.5,1,2,3,4,7,20),include.lowest=T)


pdf(paste(outdir,"txLength.pdf",sep=""))
par(mfrow=c(2,2))
hist(txLength,breaks=50,las=1,frame=F)
plot(ecdf(txLength),las=1,frame=F)
hist(log10(txLength),breaks=30,las=1,frame=F)
plot(ecdf(log10(txLength)),las=1,frame=F)
dev.off()

txLength_minus_initL<- txLength-initL
txLength_minus_newL <- txLength-newL
utrL_to_txL         <- initL/(txLength)

#maxL
tempMax      <- tapply(newL,INDEX=factor(txID),max)
maxL         <- tempMax[match(txID,names(tempMax))]
#minL
tempMin      <- tapply(newL,INDEX=factor(txID),min)
minL         <- tempMin[match(txID,names(tempMin))]

pdf(paste(outdir,"txLength_utrLength_ration.pdf",sep=""))
par(mfrow=c(2,2))
plot(density(log10(txLength)),ylab="",xlab="",frame=F,las=1,lwd=2,main=NA)
lines(density(log10(txLength_minus_initL[txLength_minus_initL>0])),col="grey",lty=2,lwd=2)
mtext(side=1,line=2,text="txLength [log10]")
mtext(side=2,line=2,text="density")

plot(density(log10(initL[match(unique(txID),txID)])),ylab="",xlab="",frame=F,las=1,lwd=2,main=NA,ylim=c(0,1))
lines(density(log10(maxL[match(unique(txID),txID)])),col="red",lty=1,lwd=2)
mtext(side=1,line=2,text="utr Length [log10]")
mtext(side=2,line=2,text="density")
legend("topright",ncol=1,bty="n",col=c("black","red"),leg=c("init","new"),pch=15)

boxplot(list(log10(initL[match(unique(txID),txID)]),log10(maxL[match(unique(txID),txID)])),outline=F,col=c("white","red"),frame=F)
mtext(side=2,line=2,text="utr Length [log10]")


hist(utrL_to_txL[utrL_to_txL>0],ylab="",xlab="",las=1,main=NA)
mtext(side=2,line=3,text="txL/utrL")

dev.off()

geneSymbol <- as.character(mcols(e75f)$X.gene_name.)[match(txID,as.character(mcols(e75f)$X.transcript_id.))]

#
#
anno_ngf     <- data.frame(uniqueID=ids,txID=txID,isoform=iso,txLength,txLength_minus_initL,txLength_minus_newL,utrL_to_txL,is.conservative,initL,newL,maxL,minL,geneSymbol)
#
#

#iso_ordered
tempL                     <- anno_ngf$newL
names(tempL)              <- anno_ngf$uniqueID
tempG                     <- as.factor(as.character(anno_ngf$txID))
myord                     <- tapply(tempL,INDEX=tempG,function(x)return(cbind(names(x)[sort(as.numeric(as.character(x)),index.return=T,decreasing=F)$ix],c(1:length(x)))))
test                      <- do.call(what=rbind,args=myord)
test                      <- test[match(anno_ngf$uniqueID,test[,1]),]

#no.iso
temp                      <- as.data.frame(table(as.character(anno_ngf$txID)))
no.iso                    <- temp[match(anno_ngf$txID,as.character(temp[,1])),2]

#
#
anno_ngf                  <- data.frame(anno_ngf,iso_ordered=test[,2],no.iso)
#
#

#myProximal
temp                      <- cbind(as.character(anno_ngf$uniqueID[anno_ngf$iso_ordered==1]),as.character(anno_ngf$txID[anno_ngf$iso_ordered==1]))
myProximal                <- temp[match(anno_ngf$txID,as.character(temp[,2])),1]
anno_ngf                  <- data.frame(anno_ngf,myProximal)

#distalToproximal
distalToproximal          <- anno_ngf$newL - anno_ngf$newL[match(as.character(anno_ngf$myProximal),as.character(anno_ngf$uniqueID))]
anno_ngf     <- data.frame(anno_ngf,distalToproximal)

#priorIsoform
tempid     <- paste(anno_ngf$txID,anno_ngf$iso_ordered,sep=".")
prior.iso  <- paste(anno_ngf$txID,(as.numeric(as.character(anno_ngf$iso_ordered))-1),sep=".")
priorID    <- anno_ngf$uniqueID[match(prior.iso,tempid)]
anno_ngf   <- data.frame(anno_ngf,priorID)

#distToprior
distoprior                           <- rep(NA,nrow(anno_ngf))
distoprior[!is.na(anno_ngf$priorID)] <- as.numeric(as.character(anno_ngf$newL))[!is.na(anno_ngf$priorID)] - as.numeric(as.character(anno_ngf$newL))[match(anno_ngf$priorID[!is.na(anno_ngf$priorID)],anno_ngf$uniqueID)]
anno_ngf                             <- data.frame(anno_ngf,distoprior=distoprior)


#nextIsoform
tempid     <- paste(anno_ngf$txID,anno_ngf$iso_ordered,sep=".")
next.iso   <- paste(anno_ngf$txID,(as.numeric(as.character(anno_ngf$iso_ordered))+1),sep=".")
nextID     <- anno_ngf$uniqueID[match(next.iso,tempid)]
anno_ngf   <- data.frame(anno_ngf,nextID)

#distTonext
distonext                            <- rep(NA,nrow(anno_ngf))
sel                                  <- !is.na(anno_ngf$nextID)
d1                                   <- as.numeric(as.character(anno_ngf$newL))[sel]
d2                                   <- as.numeric(as.character(anno_ngf$newL))[match(anno_ngf$nextID[sel],anno_ngf$uniqueID)]
distonext[sel]                       <- d2-d1
anno_ngf                             <- data.frame(anno_ngf,distonext)


# iso_cons
iso_cons                    <- rep(NA,nrow(anno_ngf))
names(iso_cons)             <- anno_ngf$uniqueID

#
# Check that the super-short are less than 122 nt
#
require(mclust)
BINS1          <- cut(anno_ngf$newL,breaks=quantile(anno_ngf$newL,prob=seq(from=0,to=1.0,by=0.1),na.rm=T),include.lowest=T)
BINS2          <- cut(log10(anno_ngf$newL),breaks=quantile(log10(anno_ngf$newL),prob=seq(from=0,to=1.0,by=0.1),na.rm=T),include.lowest=T)
BINS3          <- cut(log10(anno_ngf$newL),breaks=c(1,1.5,2,2.25,2.75,3.25,3.75,4.25),include.lowest=T)
mydist1        <- tapply(anno_ngf$distonext,INDEX=BINS1,function(x)return(log10(x[!is.na(x)])))
mydist2        <- tapply(anno_ngf$distonext,INDEX=BINS2,function(x)return(log10(x[!is.na(x)])))
mydist3        <- tapply(anno_ngf$distonext,INDEX=BINS3,function(x)return(log10(x[!is.na(x)])))
mydist4        <- tapply(anno_ngf$distonext,INDEX=factor(anno_ngf$iso_ordered,levels=as.character(c(1,2,3,4,5,6))),function(x)return(log10(x[!is.na(x)])))

Dist_to_Next        <- log10(anno_ngf$distonext[!is.na(anno_ngf$distonext)])


myutrsize <- tapply(anno_ngf$newL[anno_ngf$no.iso>2],INDEX=factor(anno_ngf$iso_ordered[anno_ngf$no.iso>2],levels=as.character(c(1,2,3,4,5,6))),function(x)return(log10(x[!is.na(x)&x>10])))
bimdens <- list()
for(i in c(1:length(myutrsize))){
  bimdens[[i]]<- densityMclust(data=myutrsize[[i]],G=1)
}
Lim       <- do.call(lapply(bimdens,function(x){
  return(10^qnorm(c(0.1,0.5,0.9),mean=x$parameters$mean[1],sd=sqrt(x$parameters$variance$sigmasq[1])))
}),what=cbind)
colnames(Lim)<-names(myutrsize)
bim2 <- densityMclust(data=myutrsize[[1]],G=3)
Lim2       <- cbind(10^qnorm(c(0.1,0.5,0.9),mean=bim2$parameters$mean[1],sd=sqrt(bim2$parameters$variance$sigmasq[1])),
                    10^qnorm(c(0.1,0.5,0.9),mean=bim2$parameters$mean[2],sd=sqrt(bim2$parameters$variance$sigmasq[1])),
                    10^qnorm(c(0.1,0.5,0.9),mean=bim2$parameters$mean[3],sd=sqrt(bim2$parameters$variance$sigmasq[1])))

colnames(Lim)<-names(myutrsize)


pdf(paste(outdir,"size_utr_per_iso_ordered_new.pdf",sep=""))
#A. All super-imposed
par(mfrow=c(1,1))
multidensity(myutrsize)

par(mfrow=c(2,2))
#A. Iso no.1
x   <- seq(from=0, to=max(myutrsize[[1]]), length=100)
hx1 <- dnorm(x,mean=bim2$parameters$mean[1],sd=sqrt(bim2$parameters$variance$sigmasq[1]))
hx2 <- dnorm(x,mean=bim2$parameters$mean[2],sd=sqrt(bim2$parameters$variance$sigmasq[1]))
hx3 <- dnorm(x,mean=bim2$parameters$mean[3],sd=sqrt(bim2$parameters$variance$sigmasq[1]))
hist(myutrsize[[1]],breaks=100,col=rgb(0,0,0,alpha=0.2),freq=FALSE,ylim=c(0,max(c(hx1,hx2))),xlab="",ylab="",main="",xaxt="n",xlim=c(1,5))
lines(x,hx3 , lwd=2, col="green")
lines(x,hx2 , lwd=2, col="blue")
lines(x,hx1 , lwd=2, col="red")
abline(v=log10(Lim2[,3]),lty=2,col="green",lwd=2)
abline(v=log10(Lim2[,1]),lty=2,col="red",lwd=2)
abline(v=log10(Lim2[,2]),lty=2,col="blue",lwd=2)
mtext(side=1,line=2,text="3' UTR length")
mtext(side=2,line=2,text="density")
#Plot nice log-log axes
full=c("1","10","100","1000","10000","100000","1000000","10000000","100000000")
y1   <- range(x)
pow  <- seq(y1[1], y1[2]+1)
ticksaty <- as.vector(sapply(pow, function(p) (1:10)*10^p))
laby <- replicate(length(ticksaty),"")
idy = match(as.numeric(full),ticksaty)
idy = idy[!is.na(idy)]
laby[idy]<-full[c(1:length(idy))]
axis(side=1, at=ticksaty, labels=TRUE, tcl=-0.25, lwd=0, lwd.ticks=0.5,cex=0.8)

#B. Remaining
for(i in c(2:length(myutrsize))){
  x   <- seq(from=0, to=max(myutrsize[[i]]), length=100)
  hx1 <- dnorm(x,mean=bimdens[[i]]$parameters$mean[1],sd=sqrt(bimdens[[i]]$parameters$variance$sigmasq[1]))
  hist(myutrsize[[i]],breaks=100,col=rgb(0,0,0,alpha=0.2),freq=FALSE,ylim=c(0,max(c(hx1,hx2))),xlab="",ylab="",main="",xaxt="n",xlim=c(1,5))
  lines(x,hx1 , lwd=2, col="red")
  abline(v=log10(Lim[,i]),lty=2,col="red",lwd=2)
  mtext(side=1,line=2,text="3' UTR length")
  mtext(side=2,line=2,text="density")
  mtext(side=3,line=0,text=paste("isoform no.",i,sep=""))
  #Plot nice log-log axes
  full=c("1","10","100","1000","10000","100000","1000000","10000000","100000000")
  y1   <- range(x)
  pow  <- seq(y1[1], y1[2]+1)
  ticksaty <- as.vector(sapply(pow, function(p) (1:10)*10^p))
  laby <- replicate(length(ticksaty),"")
  idy = match(as.numeric(full),ticksaty)
  idy = idy[!is.na(idy)]
  laby[idy]<-full[c(1:length(idy))]
  axis(side=1, at=ticksaty, labels=TRUE, tcl=-0.25, lwd=0, lwd.ticks=0.5,cex=0.8)
}
dev.off()


require(geneplotter)
pdf(paste(outdir,"distance_to_next_per_UTR_length_new.pdf",sep=""))
par(mfrow=c(2,2))
multidensity(mydist1)
multidensity(mydist2)
multidensity(mydist3)
multidensity(mydist4)
par(mfrow=c(2,1))
multidensity(myutrsize)
hist(Dist_to_Next,breaks=200,xlab="Distance to Next",freq=F)
lines(density(Dist_to_Next),col="red")
dev.off()
#


# super-short == newL<150 == -1
is.super.short                        <- anno_ngf$newL<150 #3'535
iso_cons[is.super.short]              <- -1

# then ordered BUT if distTonex<150, then Iis
bim2 <- densityMclust(data=log10(anno_ngf$distonext[!is.na(anno_ngf$distonext)]),G=2)

Lim2       <- cbind(10^qnorm(c(0.1,0.5,0.9),mean=bim2$parameters$mean[1],sd=sqrt(bim2$parameters$variance$sigmasq[1])),
                    10^qnorm(c(0.1,0.5,0.9),mean=bim2$parameters$mean[2],sd=sqrt(bim2$parameters$variance$sigmasq[2]))
)

is.short.version                          <- anno_ngf$distonext<150
is.short.version[is.na(is.short.version)] <- FALSE

# Find ordered of those which are neither short version or super short
tempL                                       <- anno_ngf$newL[!(is.short.version|is.super.short)]
names(tempL)                                <- anno_ngf$uniqueID[!(is.short.version|is.super.short)]
tempG                                       <- as.factor(as.character(anno_ngf$txID[!(is.short.version|is.super.short)]))
myord                                       <- tapply(tempL,INDEX=tempG,function(x)return(cbind(names(x)[sort(as.numeric(as.character(x)),index.return=T,decreasing=F)$ix],c(0:(length(x)-1)))))
test                                        <- do.call(what=rbind,args=myord)
iso_cons[match(test[,1],names(iso_cons))]   <- test[,2]
# Finally is.short.version becomes is
idsnext                    <- anno_ngf$nextID[is.short.version]
# check if next is not itself the next of a next
is.next.short              <- is.short.version[match(idsnext,names(iso_cons))]#1'332

iso_cons[which(is.short.version)[!is.next.short]] <- paste(iso_cons[match(idsnext[!is.next.short],names(iso_cons))],"s",sep="")
iso_cons[which(is.short.version)[is.next.short]]  <- paste(iso_cons[match(idsnext[is.next.short],names(iso_cons))],"s",sep="")
anno_ngf$iso_cons         <- iso_cons

# ioID
tempL                                       <- as.character(anno_ngf$iso_cons)
names(tempL)                                <- as.character(anno_ngf$uniqueID)
tempL[grep(tempL,pattern="s")]              <-"none"
ioID                                        <- tapply(tempL,INDEX=as.factor(as.character(anno_ngf$txID)),function(x)return(names(x)[grep(x,pattern="0")]))
ioID                                        <- as.character(ioID[match(anno_ngf$txID,names(ioID))])
anno_ngf                                    <- data.frame(anno_ngf,ioID)


# imID
tempL                                         <- as.character(anno_ngf$iso_cons)
names(tempL)                                  <- as.character(anno_ngf$uniqueID)
tempL[grep(tempL,pattern="s")]                <-"none"

imID1                                         <- tapply(tempL,INDEX=as.factor(as.character(anno_ngf$txID)),function(x)return(names(x)[grep(x,pattern="-1")]))
imID1                                         <- as.character(imID1[match(anno_ngf$txID,names(imID1))])

imID2                                         <- tapply(tempL,INDEX=as.factor(as.character(anno_ngf$txID)),function(x)return(names(x)[grep(x,pattern="0")[1]]))
imID2                                         <- as.character(imID2[match(anno_ngf$txID,names(imID2))])

imID                                          <- imID1
imID[is.na(match(imID1,anno_ngf$uniqueID))]   <- imID2[is.na(match(imID1,anno_ngf$uniqueID))]
anno_ngf                                    <- data.frame(anno_ngf,imID)


#iso_merged
iso_cons_merged            <- anno_ngf$iso_cons
iso_cons_merged            <- sapply(iso_cons_merged,function(X)return(gsub(X,pattern="[s,ss]",repl="")))
iso_cons_merged[which(!iso_cons_merged%in%c("-1","0","1","2","3","4"))]<-NA
anno_ngf$iso_cons_merged   <- iso_cons_merged
iso_cons_focus             <- anno_ngf$iso_cons
for(i in c("-1s","0s","1s","2s","3s","4s")){
  iso_cons_focus[setdiff(grep(iso_cons_focus,pattern=i),grep(iso_cons_focus,pattern=paste("-",i,sep="")))] <- i
}
iso_cons_focus[setdiff(c(1:length(iso_cons_focus)),unlist(lapply(X=c("-1","0","1","2","3","4"),                                                               FUN=function(X)return(grep(iso_cons_focus,pattern=X)))))]<-NA
anno_ngf$iso_cons_focus         <- iso_cons_focus




#
# 2. Merge with coverage -- think a bit before moving forward of the meaning of summing-up before or after log2
#
anno_ngf   <- data.frame(anno_ngf,myOut[match(anno_ngf$uniqueID,rownames(myOut)),])

#Check data are log-normal
Removepattern <- function(x,pat)return(unlist(lapply(x,function(x)return(gsub(x,pattern=pat,repl="")))))

mydat   <- anno_ngf[,grep(colnames(anno_ngf),pattern=".raw")]
groups  <- factor(Removepattern(x=Removepattern(x=colnames(mydat),pat="\\.1.raw"),pat="\\.2.raw"))
mymean  <- t(apply(mydat,1,function(x)return(tapply(x,INDEX=groups,FUN=mean))))
mysd    <- t(apply(mydat,1,function(x)return(tapply(x,INDEX=groups,FUN=sd))))
myvar   <- t(apply(mydat,1,function(x)return(tapply(x,INDEX=groups,FUN=function(Z)return(sd(Z)^2)))))

mymeanL  <- t(apply(mydat,1,function(x)return(tapply(x,INDEX=groups,FUN=function(Z)return(mean(log2(1+Z)))))))
mysdL    <- t(apply(mydat,1,function(x)return(tapply(x,INDEX=groups,FUN=function(Z)return(sd(log2(1+Z)))))))
mysvarL  <- t(apply(mydat,1,function(x)return(tapply(x,INDEX=groups,FUN=function(Z)return(sd(log2(1+Z))^2)))))



pdf(paste(outdir,"mean_sd.pdf",sep=""))
par(mfrow=c(2,2))
for(i in c(1:ncol(myvar))){
    smoothScatter(mymean[,i],myvar[,i],xlab="mean",ylab="variance",main=colnames(myvar)[i])
}
par(mfrow=c(2,2))
for(i in c(1:ncol(myvar))){
    smoothScatter(log10(mymean[,i]),log10(myvar[,i]),xlab="mean",ylab="variance",main=colnames(myvar)[i])
}
for(i in c(1:ncol(myvar))){
    smoothScatter(mymeanL[mymeanL[,i]>0,i],mysvarL[mymeanL[,i]>0,i],xlab="mean",ylab="variance",main=colnames(myvar)[i])
}

dev.off()

#Create data of interest: sum of the isoform across biological replicates
dat                         <- anno_ngf[,grep(colnames(anno_ngf),pattern=".raw")]
groups                      <- factor(Removepattern(x=Removepattern(x=colnames(dat),pat="\\.1.raw"),pat="\\.2.raw"))
anocol                       <- do.call(lapply(colnames(dat),function(x)return(unlist(strsplit(x,split="[_,\\.]")))),what=rbind)

myMean.iso                  <- t(apply(dat,1,function(x)return(tapply(x,INDEX=groups,FUN=function(Z)return(mean(log2(Z+1)))))))
colnames(myMean.iso)        <- unlist(lapply(colnames(myMean.iso),function(x)return(paste(x,"mean.log2",sep="."))))

mySD.iso                    <- t(apply(dat,1,function(x)return(tapply(x,INDEX=groups,FUN=function(Z)return(sd(log2(Z+1)))))))
colnames(mySD.iso)          <- unlist(lapply(colnames(mySD.iso),function(x)return(paste(x,"sd.log2",sep="."))))

mySum.iso                   <- t(apply(dat,1,function(x)return(tapply(x,INDEX=groups,FUN=function(Z)return(log2(1+sum(Z)))))))
colnames(mySum.iso)         <- unlist(lapply(colnames(mySum.iso),function(x)return(paste(x,"log2.sum",sep="."))))

mySum.iso.p                 <- t(apply(dat,1,function(x)return(tapply(x,INDEX=groups,FUN=function(Z)return(sum(log2(1+Z)))))))
colnames(mySum.iso.p)       <- unlist(lapply(colnames(mySum.iso.p),function(x)return(paste(x,"sum.log2",sep="."))))

mySum.iso.raw               <- t(apply(dat,1,function(x)return(tapply(x,INDEX=groups,FUN=sum))))
colnames(mySum.iso.raw)     <- unlist(lapply(colnames(mySum.iso.raw),function(x)return(paste(x,"sum.raw",sep="."))))

anno_ngf                    <- data.frame(anno_ngf,myMean.iso,mySD.iso,mySum.iso,mySum.iso.p,mySum.iso.raw)



#maxTxID (maximum of the average)
coloi                       <- grep(colnames(anno_ngf),pattern="mean.log2")
subdat                      <- anno_ngf[,coloi]

maxTxID                     <- apply(subdat,2,function(x)return(tapply(x,INDEX=as.factor(as.character(anno_ngf$txID)),FUN=function(x)return(max(x,na.rm=T)))))
colnames(maxTxID)           <- paste(colnames(subdat),".maxTx",sep="")
maxTxID                     <- maxTxID[match(anno_ngf$txID,rownames(maxTxID)),]
maxTxID[maxTxID==-Inf]      <- NA
anno_ngf                    <- data.frame(anno_ngf,maxTxID)


#minTxID (minimum of the average)
minTxID                     <- apply(subdat,2,function(x)return(tapply(x,INDEX=as.factor(as.character(anno_ngf$txID)),FUN=function(x)return(min(x,na.rm=T)))))
colnames(minTxID)           <- paste(colnames(subdat),".minTx",sep="")
minTxID                     <- minTxID[match(anno_ngf$txID,rownames(minTxID)),]
minTxID[minTxID==-Inf]      <- NA
anno_ngf                   <- data.frame(anno_ngf,minTxID)

#sumTxID (of the average)
sumTxID                     <- apply(subdat,2,function(x)return(tapply(x,INDEX=as.factor(as.character(anno_ngf$txID)),FUN=function(x)return(sum(x,na.rm=T)))))
colnames(sumTxID)           <- paste(colnames(subdat),".sumTx",sep="")
sumTxID                     <- sumTxID[match(anno_ngf$txID,rownames(sumTxID)),]
anno_ngf                   <- data.frame(anno_ngf,sumTxID)



#sumTxIDRawPerSample (of the iso)
coloi                      <- grep(colnames(anno_ngf),pattern=".raw")[-c(9:12)]
sumTxRaw                   <- apply(anno_ngf[,coloi],2,function(x)return(tapply(x,INDEX=as.factor(as.character(anno_ngf$txID)),FUN=function(x)return(sum(x,na.rm=T)))))
colnames(sumTxRaw)         <- paste(Removepattern(colnames(sumTxRaw),pat="\\.raw"),".sumTxRaw",sep="")
sumTxRaw                   <- sumTxRaw[match(anno_ngf$txID,rownames(sumTxRaw)),]
anno_ngf                   <- data.frame(anno_ngf,sumTxRaw)

#sumlog2TxIDPerSample (of the iso)
coloi                      <- grep(colnames(anno_ngf),pattern=".raw")[-c(9:12)]
sumTxRaw                   <- apply(anno_ngf[,coloi],2,function(x)return(tapply(x,INDEX=as.factor(as.character(anno_ngf$txID)),FUN=function(x)return(sum(log2(1+x),na.rm=T)))))

colnames(sumTxRaw)         <- paste(Removepattern(colnames(sumTxRaw),pat="\\.raw"),".sumlog2Tx",sep="")
sumTxRaw                   <- sumTxRaw[match(anno_ngf$txID,rownames(sumTxRaw)),]
anno_ngf                   <- data.frame(anno_ngf,sumTxRaw)

dat <- anno_ngf[match(unique(anno_ngf$txID),anno_ngf$txID),grep(colnames(anno_ngf),pattern="sumlog2Tx")]
datL<- apply(dat,2,function(x)return(x[x>0]))

pdf(paste(outdir,"analysis_tx_log2.pdf",sep=""))
par(mfrow=c(2,2))
samples <- c("NGF.axon","NGF.cb","NT3.axon","NT3.cb")
lapply(samples,function(x)return(smoothScatter(x=dat[,grep(colnames(dat),pattern=x)[1]],
                                               y=dat[,grep(colnames(dat),pattern=x)[2]],
                                               main=x
                                               )))
layout(matrix(c(1,2,1,3), 2, 2, byrow = TRUE))
boxplot(datL,outline=F)
multidensity(datL[grep(colnames(dat),pattern="axon")])
multidensity(datL[grep(colnames(dat),pattern="cb")])

par(mfrow=c(2,4))
lapply(datL,function(x)return(hist(x,breaks=150)))
dev.off()


#MaxID and MinID
tempL                      <- anno_ngf$newL
names(tempL)               <- anno_ngf$uniqueID
minID                      <- tapply(tempL,INDEX=as.factor(as.character(anno_ngf$txID)),function(x)return(names(x)[x==min(x)]))
minID                      <- minID[match(anno_ngf$txID,names(minID))]
maxID                      <- tapply(tempL,INDEX=as.factor(as.character(anno_ngf$txID)),function(x)return(names(x)[x==max(x)]))
maxID                      <- maxID[match(anno_ngf$txID,names(maxID))]
anno_ngf$MaxID             <- maxID
anno_ngf$minID             <- minID








#
# Compare length before and after pipeline
#

L3              <- ngfGRS
myL_init        <- ComputeWidth(g3utr)
names(myL_init) <- g3utr$X.transcript_id.
myL_init        <- myL_init[match(L3$txID,names(myL_init))]
myL_new         <- ComputeWidth(L3)
names(myL_new)  <- L3$ID
mydL            <- DeltaLp(L3)
names(mydL)     <- L3$ID
myDat           <- data.frame(txID=names(myL_init),uniqueID=names(myL_new),myL_init,myL_new,mydL,detect=anno_ngf$detect.txID.ngf,max_eCB=anno_ngf$NGF.cb.mean.log2.maxTx,max_eAx=anno_ngf$NGF.axon.mean.log2.maxTx)

#Before computing which one is the longest, you must remove the iso which are not expressed
sel1            <- union(which(anno_ngf$detect.iso.ngf!="no"),grep(anno_ngf$uniqueID,pattern="\\.0"))
myDat           <- myDat[sel1,]

#Focus on the maxL per category
temp            <- myDat$myL_new
names(temp)     <- as.character(myDat$uniqueID)
selmaxL         <- tapply(temp,INDEX=as.factor(myDat$txID),FUN=function(x)return(names(x)[which(x==max(x))[1]]))
ix              <- match(selmaxL,as.character(myDat$uniqueID))
maxL            <- myDat$myL_new[ix]
maxL            <- maxL[match(as.character(myDat$txID),as.character(myDat$txID)[ix])]
myDat$maxL      <- maxL
maxdL           <- myDat$mydL[ix]
maxdL           <- maxdL[match(as.character(myDat$txID),as.character(myDat$txID)[ix])]
myDat$maxdL     <- maxdL

#Take single txID
myDat           <- myDat[match(unique(as.character(myDat$txID)),as.character(myDat$txID)),]#17'418 txID of which 14'116 are expressed in neurones


#Focus on extended and expressed (8'617 txID)
elong           <- which(myDat$maxdL>=10)#8'619
focus_elong     <- list(rn5_init=myDat$myL_init[elong],denovo=myDat$maxL[elong])

#Focus on expressed
lappend <- function (lst, ...){
  lst <- c(lst, list(...))
  return(lst)
}
byComp          <- tapply(myDat$maxL,INDEX=myDat$detect,function(x)return(x))
byComp          <- lappend(
                    byComp,
                    neurones_rn5=myDat$myL_init[myDat$detect!="no"],
                    neurones_new=myDat$maxL[myDat$detect!="no"]
                    )

focus_expressd  <- byComp[match(c("neurones_rn5","neurones_new","cb","neur"),names(byComp))]
names(focus_expressd)[c(3,4)]<- c("cb-restricted","targetd to axons")


pdf(paste(outdir,"elongation_pipeline_expressed_elong.pdf",sep=""))
par(mfrow=c(1,2))
boxplot(focus_elong,outline=FALSE,las=1,frame=F)
mtext(side=3,line=0,text=unlist(lapply(focus_elong,FUN=length)),at=c(1:length(focus_elong)))
mtext(side=3,line=2,text="maximum length of transcript elongated")
boxplot(focus_expressd,outline=FALSE,las=1,frame=F)
mtext(side=3,line=0,text=unlist(lapply(focus_expressd,FUN=length)),at=c(1:length(focus_expressd)))
mtext(side=3,line=2,text="maximum length of transcript expressed")
dev.off()



#Test whether there is a correlation between level of expression and level of lenghtening
byComp          <- tapply(myDat$max_eCB[myDat$myL_init>=10],INDEX=myDat$detect[myDat$myL_init>=10],function(x)return(x))

pdf(paste(outdir,"correlation_expression_extension.pdf",sep=""))
par(mfrow=c(2,2))
boxplot(byComp[c(2,3)],outline=FALSE)
smoothScatter(myDat$maxL[myDat$myL_init>=10],myDat$max_eCB[myDat$myL_init>=10],xlab="3' UTR length",ylab="Maximum coverage CB")
smoothScatter(log10(myDat$maxL[myDat$myL_init>=10]),myDat$max_eCB[myDat$myL_init>=10],xlab="3' UTR length [log10]",ylab="Maximum coverage CB")
smoothScatter(log10(1+myDat$maxdL[myDat$myL_init>=10]),myDat$max_eCB[myDat$myL_init>=10],xlab="3max dL",ylab="Maximum coverage CB")

#The best check is with I0


byComp          <- tapply(myDat$myL_init[myDat$myL_init>=10],INDEX=myDat$detect[myDat$myL_init>=10],function(x)return(x))
byComp          <- lappend(
                    byComp,
                    neurones_rn5=myDat$myL_init[myDat$detect!="no"]
                    )

focus_expressd  <- byComp[match(c("neurones_rn5","neurones_new","cb","neur"),names(byComp))]
names(focus_expressd)[c(3,4)]<- c("cb-restricted","targetd to axons")


par(mfrow=c(1,2))
boxplot(focus_elong,outline=FALSE,las=1,frame=F)
mtext(side=3,line=0,text=unlist(lapply(focus_elong,FUN=length)),at=c(1:length(focus_elong)))
mtext(side=3,line=2,text="maximum length of transcript elongated")
boxplot(focus_expressd,outline=FALSE,las=1,frame=F)
mtext(side=3,line=0,text=unlist(lapply(focus_expressd,FUN=length)),at=c(1:length(focus_expressd)))
mtext(side=3,line=2,text="maximum length of transcript expressed")
dev.off()


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

SelectBasedOnDistance(famtxID=anno_ngf[anno_ngf$txID=="ENSRNOT00000014274",])
SelectBasedOnDistance(famtxID=anno_ngf[anno_ngf$txID=="ENSRNOT00000016995",])

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


save(list=c("anno_ngf","ngfGRS"),file="/home/rluisier/data/Riccio/Exp_1/Dec2016/utrid/anno_ngf_stringent/anno_ngf_March_12.RData")
#save(list=c("anno_ngf","ngfGRS"),file="/home/rluisier/data/Riccio/Exp_1/Dec2016/utrid/anno_ngf_stringent/anno_ngf_Feb_02.RData")

ngfGRS <- ngfGRS[which(ngfGRS$uniqueID%in%anno_ngf$uniqueID[anno_ngf$selected.for.distance]),]
export.gff(object=ngfGRS,con="/home/rluisier/data/Riccio/Exp_1/Dec2016/utrid/APA_stringent/Lngf_sub.gtf",format="gtf")






