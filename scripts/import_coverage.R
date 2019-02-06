library(grDevices)
library(Rsamtools)
library(IRanges)
library(GenomicRanges)
library(rtracklayer)

load("./annotation/GRanges_comprehensive_transcriptome_rat_24_nov_2015.RData")
ngfGRS         <- import.gff("./utrid/APA/L2.gtf",format="gtf",asRangedData=FALSE)
covdir         <- "./Coverage/utrCov"


# A. Import coverage

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


#B. Log2-transform
IX.norm.500 <- grep(colnames(myCov500),pattern="norm")
IX.raw.500  <- grep(colnames(myCov500),pattern="raw")
MIN1        <-min(as.vector(myCov500[,IX.norm.500])[as.vector(myCov500[,IX.norm.500]>0)])
myCov500l  <- cbind(log2(myCov500[,IX.raw.500]+1),
                    log2(myCov500[,IX.norm.500]+MIN1)
)
myCov500l <- myCov500l[,match(colnames(myCov500),colnames(myCov500l))]
myCov500  <- myCov500[,match(colnames(myCov500l),colnames(myCov500))]


#C. Create annotation
info                 <- do.call(lapply(colnames(myCov500l),function(x)return(unlist(strsplit(x,split="\\.")))),what=rbind)
info                 <- data.frame(colnames(myCov500l),info)
colnames(info)       <-c("id","treatment","origin","replicate","type")
save(list=c("myCov500","myCov500l","info"),file="/home/rluisier/data/Riccio/Exp_1/Dec2016/Coverage/utrCov_stringent/coverageUTR_March_12.RData")



# D. Identify reliably expressed genes. Please note that it seems more accurate to use the 22'000 tX rather than the entire pool to determine the threshold.

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
mylistCov300  <- lapply(c(1:8),function(x)return(log2(dat2[dat2[,x]>=5.012531e-02,x])))

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
selNT3 <- union(which(ngfGRS$ID%in%rownames(is.expressed.iso)[apply(is.expressed.iso[,c(3,4)],1,function(x)return(sum(x)>0))]),
            grep("\\.0",ngfGRS$ID))
seltot <- union(which(ngfGRS$ID%in%rownames(is.expressed.iso)[apply(is.expressed.iso,1,function(x)return(sum(x)>0))]),
            grep("\\.0",ngfGRS$ID))

#Selection based on I0


#I should indeed keep also the isoform I0 for downstream analysis
Lngf                        <- ngfGRS[selNGF,]
Lnt3                        <- ngfGRS[selNT3,]
Ltot                        <- ngfGRS[seltot,]
export.gff(Lngf,con="/home/rluisier/data/Riccio/Exp_1/Dec2016/utrid/APA_stringent/Lngf.gtf",format="gtf")#50'080
export.gff(Lnt3,con="/home/rluisier/data/Riccio/Exp_1/Dec2016/utrid/APA_stringent/Lnt3.gtf",format="gtf")#49'853
export.gff(Ltot,con="/home/rluisier/data/Riccio/Exp_1/Dec2016/utrid/APA_stringent/Ltot.gtf",format="gtf")#50'975


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
Lnt3 <- RemoveIds(Lnt3)
Ltot <- RemoveIds(Ltot)

export.gff(Lngf,con="/home/rluisier/data/Riccio/Exp_1/Dec2016/utrid/APA_stringent/Lngfp.gtf",format="gtf")#PLease note that this version does not contain those 3' end that do not overlap with a polyA site however without duplicates
export.gff(Lnt3,con="/home/rluisier/data/Riccio/Exp_1/Dec2016/utrid/APA_stringent/Lnt3p.gtf",format="gtf")
export.gff(Ltot,con="/home/rluisier/data/Riccio/Exp_1/Dec2016/utrid/APA_stringent/Ltotp.gtf",format="gtf")


#is.expressed.txid and GS
#is.expressed.iso --> must be expressed in both samples
allGS                       <- as.character(mcols(e75p)$X.gene_name.)
allTX                       <- as.character(mcols(e75p)$X.transcript_id.)
myselG                      <- myselG[match(Ltot$uniqueID,rownames(myselG)),]
myCov500                    <- myCov500[match(Ltot$uniqueID,rownames(myCov500)),]
myCov300                    <- myCov300[match(Ltot$uniqueID,rownames(myCov300)),]
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



#Detection per compartment per condition: NGF
sel1                        <- is.expressed.iso$NGF.cb.is.expressed.iso#all isoforms which are detected in CB
sel2                        <- is.expressed.iso$NGF.axon.is.expressed.iso#all isoforms detected in axons

txID.cb                     <- unique(txID[sel1])#all txID detected in CB --> 14'698
txID.ax                     <- unique(txID[sel2])#all txID detected in axons --> 5'722
GS.cb                       <- unique(GS[sel1])#all GS detected in CB --> 13'084
GS.ax                       <- unique(GS[sel2])#all GS detected in axons --> 5'109

selAll                      <- txID%in%union(txID.cb,txID.ax) #txID detected in axons and/or cell body
selCB                       <- txID%in%setdiff(txID.cb,txID.ax) #txID detected in CB only
selNeurons                  <- txID%in%intersect(txID.cb,txID.ax) #txID detected in both axons AND cell body
selAxonsonly                <- txID%in%setdiff(txID.ax,txID.cb) #txID detected in both axons and cell body

detect.txID                 <-  rep("no",length(txID))
detect.txID[selAxonsonly]   <-  "ax"
detect.txID[selNeurons]     <-  "neur"
detect.txID[selCB]          <-  "cb"

selAx.Iso.axons             <- selNeurons&sel2
selAx.Iso.CB                <- selNeurons&sel1&!sel2
selAx.Iso.axons.only        <- selNeurons&sel2&!sel1

detect.iso                                      <- rep("no",length(txID))
detect.iso[detect.txID=="cb"&sel1]              <- "cb.cb"
detect.iso[detect.txID=="ax"&sel2]              <- "ax.ax"
detect.iso[detect.txID=="neur"&sel1&sel2]       <- "neur.neur"
detect.iso[detect.txID=="neur"&sel1&!sel2]      <- "neur.cb"
detect.iso[detect.txID=="neur"&sel2&!sel1]      <- "neur.ax"

tot.detect.iso                                  <- data.frame(detect.iso)
colnames(tot.detect.iso)                        <- "detect.iso.ngf"
tot.detect.txID                                 <- data.frame(detect.txID)
colnames(tot.detect.txID)                       <- "detect.txID.ngf"


selAll                      <- GS%in%union(GS.cb,GS.ax) #txID detected in axons and/or cell body
selCB                       <- GS%in%setdiff(GS.cb,GS.ax) #txID detected in CB only
selNeurons                  <- GS%in%intersect(GS.cb,GS.ax) #txID detected in both axons AND cell body
selAxonsonly                <- GS%in%setdiff(GS.ax,GS.cb) #txID detected in both axons and cell body

detect.GS                 <-  rep("no",length(GS))
detect.GS[selAxonsonly]   <-  "ax"
detect.GS[selNeurons]     <-  "neur"
detect.GS[selCB]          <-  "cb"
tot.detect.gs             <- data.frame(detect.GS)
colnames(tot.detect.gs)   <- "detect.GS.ngf"

#Detection per compartment per condition: NT3
sel1                        <- is.expressed.iso$NT3.cb.is.expressed.iso#all isoforms which are detected in CB
sel2                        <- is.expressed.iso$NT3.axon.is.expressed.iso#all isoforms detected in axons

txID.cb                     <- unique(txID[sel1])#all txID detected in CB --> 14'108
txID.ax                     <- unique(txID[sel2])#all txID detected in axons --> 5'129
GS.cb                       <- unique(GS[sel1])#all GS detected in CB --> 19'295
GS.ax                       <- unique(GS[sel2])#all GS detected in axons --> 4'754

selAll                      <- txID%in%union(txID.cb,txID.ax) #txID detected in axons and/or cell body
selCB                       <- txID%in%setdiff(txID.cb,txID.ax) #txID detected in CB only
selNeurons                  <- txID%in%intersect(txID.cb,txID.ax) #txID detected in both axons AND cell body
selAxonsonly                <- txID%in%setdiff(txID.ax,txID.cb) #txID detected in both axons and cell body

detect.txID                 <-  rep("no",length(txID))
detect.txID[selAxonsonly]   <-  "ax"
detect.txID[selNeurons]     <-  "neur"
detect.txID[selCB]          <-  "cb"

selAx.Iso.axons             <- selNeurons&sel2
selAx.Iso.CB                <- selNeurons&sel1&!sel2
selAx.Iso.axons.only        <- selNeurons&sel2&!sel1

detect.iso                                      <- rep("no",length(txID))
detect.iso[detect.txID=="cb"&sel1]              <- "cb.cb"
detect.iso[detect.txID=="ax"&sel2]              <- "ax.ax"
detect.iso[detect.txID=="neur"&sel1&sel2]       <- "neur.neur"
detect.iso[detect.txID=="neur"&sel1&!sel2]      <- "neur.cb"
detect.iso[detect.txID=="neur"&sel2&!sel1]      <- "neur.ax"

tot.detect.iso                                  <- data.frame(tot.detect.iso,detect.iso)
colnames(tot.detect.iso)[2]                    <- "detect.iso.nt3"
tot.detect.txID                                 <- data.frame(tot.detect.txID ,detect.txID)
colnames(tot.detect.txID)[2]                    <- "detect.txID.nt3"


selAll                      <- GS%in%union(GS.cb,GS.ax) #txID detected in axons and/or cell body
selCB                       <- GS%in%setdiff(GS.cb,GS.ax) #txID detected in CB only
selNeurons                  <- GS%in%intersect(GS.cb,GS.ax) #txID detected in both axons AND cell body
selAxonsonly                <- GS%in%setdiff(GS.ax,GS.cb) #txID detected in both axons and cell body

detect.GS                 <-  rep("no",length(GS))
detect.GS[selAxonsonly]   <-  "ax"
detect.GS[selNeurons]     <-  "neur"
detect.GS[selCB]          <-  "cb"

tot.detect.gs             <- data.frame(tot.detect.gs,detect.GS)
colnames(tot.detect.gs)[2]   <- "detect.GS.nt3"

tot.detect                  <- data.frame(tot.detect.iso ,tot.detect.txID ,tot.detect.gs )
rownames(tot.detect)<- rownames(myselG)


txID                        <- unlist(lapply(rownames(myselG),function(x)return(unlist(strsplit(x,split="\\."))[1])))
GS                          <- allGS[match(txID,allTX)]

tot.detect <- data.frame(txID,GS,tot.detect)

save(list=c("tot.detect","myselG","myCov500","myCov300","Lngf","Lnt3","Ltot","txID","GS"),file= "/home/rluisier/data/Riccio/Exp_1/Dec2016/Coverage/utrCov_stringent/coverage_detection.RData")

unique(GS[tot.detect$detect.iso.ngf=="ax.ax"])

pdf(paste(outdir,"stat_expression_compartments_500_v_cor.pdf",sep=""))

layout(matrix(c(1,2,3,1,4,5), 2, 3, byrow = TRUE))

#Analysis no.expressed transcript per compartments I
temp1 <- as.data.frame(table(tot.detect$detect.txID.ngf[match(unique(tot.detect$txID),tot.detect$txID)]))[c(2,3,1),]
temp2 <- as.data.frame(table(tot.detect$detect.txID.nt3[match(unique(tot.detect$txID),tot.detect$txID)]))[c(2,3,1),]
myval <- c(  sum(temp1$Freq),sum(temp2$Freq),
    sum(temp1$Freq[c(1,2)]),sum(temp2$Freq[c(1,2)]),
    sum(temp1$Freq[c(3,2)]),sum(temp2$Freq[c(3,2)]))
names(myval)<- rep("",6)
mp          <- barplot(myval,las=1,col=c("grey","grey","blue","blue","green","green"))
mtext(at=mp,side=3,line=0,text=myval,cex=0.5)
mtext(at=mp,side=1,line=0,text=c("all.ngf","all.nt3","cb.ngf","cb.nt3","ax.ngf","ax.nt3"),cex=0.5)
mtext(side=2,line=2,text="no.tX")

#Analysis no.expressed transcript per compartments II
myval           <- tot.detect$detect.txID.ngf[match(unique(tot.detect$txID),tot.detect$txID)]
mp              <- barplot(table(myval)[c(2,3,1)],col=c("orange","green","black"),las=1,cex.names=0.5,ylim=c(0,8000))
mtext(at=mp,side=3,line=0,text=as.vector(table(myval))[c(2,3,1)],cex=0.5)
mtext(side=2,line=3,text="no.tx",cex=0.5)
mtext(side=3,line=1,text="NGF",cex=0.5)
myval           <- tot.detect$detect.txID.nt3[match(unique(tot.detect$txID),tot.detect$txID)]
mp              <- barplot(table(myval)[c(2,3,1)],col=c("orange","green","black"),las=1,cex.names=0.5,ylim=c(0,8000))
mtext(at=mp,side=3,line=0,text=as.vector(table(myval))[c(2,3,1)],cex=0.5)
mtext(side=2,line=3,text="no.tx",cex=0.5)
mtext(side=3,line=1,text="NT3",cex=0.5)

#Analysis no.expressed isoforoms per compartments
myval   <- tot.detect$detect.iso.ngf
mp<- barplot(table(myval)[c(2,4,5,3,1)],col=c("orange","orange","green","black","black"),las=1,cex.names=0.5,ylim=c(0,20000))
mtext(at=mp,side=3,line=0,text=as.vector(table(myval))[c(2,4,5,3,1)],cex=0.5)
mtext(side=2,line=3,text="no.isoforms",cex=0.5)
mtext(side=3,line=1,text="NGF",cex=0.5)
myval   <- tot.detect$detect.iso.nt3
mp<- barplot(table(myval)[c(2,4,5,3,1)],col=c("orange","orange","green","black","black"),las=1,cex.names=0.5,ylim=c(0,20000))
mtext(at=mp,side=3,line=0,text=as.vector(table(myval))[c(2,4,5,3,1)],cex=0.5)
mtext(side=2,line=3,text="no.isoforms",cex=0.5)
mtext(side=3,line=1,text="NGF",cex=0.5)






