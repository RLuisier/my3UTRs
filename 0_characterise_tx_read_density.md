1. In BASH get the coverage from all chromosomes

```bash
OUTDIR=./Coverage/txCoverage/
ANNOTATION=./annotation/GRanges_comprehensive_transcriptome_rat_24_nov_2015.RData
awk '{ print $1}' ./annotation/rn5/chrom_oi.txt  | while read myline
do
    ./scripts/getTxCoverage.R $OUTDIR $CHR $ANNOTATION $BAMDIR
    echo " "
    echo " "
done
````


2. Import all chromosome files
```R
library(grDevices)
library(Rsamtools)
library(IRanges)
library(GenomicRanges)
library(rtracklayer)

foi <- list.files("./Coverage/txCoverage")
foi <- foi[grep(foi,pattern=".RData")]

load(paste("./Coverage/txCoverage/",foi[1],sep=""))
myCoverage <- mySum
for(i in c(2:length(foi))){
    load(paste("./Coverage/txCoverage/",foi[i],sep=""))
    for(i in c(1:length(mySum))){
        for(j in c(1:length(mySum[[i]]))){
            myCoverage[[i]][[j]] <- myCoverage[[i]][[j]]+mySum[[i]][[j]]
            print(j)
        }
    }
}


names   <- unlist(lapply(unlist(lapply(system(paste("ls -d ",bamdir,"*.bam",sep=""),intern=T),function(x)return(unlist(strsplit(x,split="/"))[10]))),function(z)return(gsub(z,pattern=".sorted.bam",repl=""))))
names(myCoverage)<-names
```


3. Plot the coverage density focussing on the last 1000 nucletotides

```R
for(i in c(1:length(myCoverage))){
    for(j in c(1:length(myCoverage[[i]]))){
        sf                    <- 100/sum(myCoverage[[i]][[j]])
        print(sf)
        myCoverage[[i]][[j]] <- sf*myCoverage[[i]][[j]]
    }
}


myCumSum <- list()
for(i in c(1:length(myCoverage))){
    myCumSum[[i]] <- list()
    for(j in c(1:length(myCoverage[[i]]))){
        myCumSum[[i]][[j]]    <- 100-cumsum(myCoverage[[i]][[j]])
    }
}

outdir <- "./Coverage/txCoverage/"


pdf(paste(outdir,"perSamplepertxLength.pdf",sep=""))
for(i in c(1:length(myCoverage))){
    par(mfrow=c(3,4))
    for(j in c(1:length(myCoverage[[i]]))){
        plot(myCumSum[[i]][[j]],type="l",frame=F,mian="",xlab="",ylab="")
        abline(h=70,col="red",lty=2)
        IX<-which(myCumSum[[i]][[j]]>70)
        IX<- IX[length(IX)]
        abline(v=IX,col="red",lty=2)
        myd <- length(myCumSum[[i]][[j]])-IX
        mtext(side=3,line=0,text=paste("d = ",myd," nt",sep=""),cex=0.7)
        mtext(side=1,line=2,text="distance from 3'end",cex=0.7)
        mtext(side=2,line=2,text="cumulative read density",cex=0.7)

        plot(myCoverage[[i]][[j]],type="l",frame=F,mian="",xlab="",ylab="")
        abline(v=IX,col="red",lty=2)
        mtext(side=3,line=0,text=paste("d = ",myd," nt",sep=""),cex=0.7)
        mtext(side=1,line=2,text="distance from 3'end",cex=0.7)
        mtext(side=2,line=2,text="read coverage density",cex=0.7)


    }
}
dev.off()


avgCumSum  <- list()
sdCumSum   <- list()
avgCovDens <- list()
sdCovDens  <- list()

groups <- unlist(lapply(names,function(x)return(gsub(x,pattern="[_1,_2]",repl=""))))
Lev    <- unique(groups)
for(i in c(1:length(Lev))){
    A=myCumSum[[which(groups==Lev[i])[1]]]
    B=myCumSum[[which(groups==Lev[i])[2]]]
    avgCumSum[[i]] <- lapply(c(1:length(A)),function(X){ return( apply(cbind(A[[X]],B[[X]]),1,mean)) })
    sdCumSum[[i]]  <- lapply(c(1:length(A)),function(X){ return( apply(cbind(A[[X]],B[[X]]),1,sd)) })
    A=myCoverage[[which(groups==Lev[i])[1]]]
    B=myCoverage[[which(groups==Lev[i])[2]]]
    avgCovDens[[i]] <- lapply(c(1:length(A)),function(X){ return( apply(cbind(A[[X]],B[[X]]),1,mean)) })
    sdCovDens[[i]]  <- lapply(c(1:length(A)),function(X){ return( apply(cbind(A[[X]],B[[X]]),1,sd)) })
    print(i)
}

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
      if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
      stop("vectors must be same length")
      arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
    }



pdf(paste(outdir,"merged_samples_of_interest.pdf",sep=""))

par(mfrow=c(3,4))
for(i in c(1:length(avgCumSum[[1]]))){
        MIN <- min(c(min(avgCovDens[[1]][[i]]-sdCovDens[[1]][[i]]),min(avgCovDens[[2]][[i]]-sdCovDens[[2]][[i]])))
        MAX <- max(c(max(avgCovDens[[1]][[i]]+sdCovDens[[1]][[i]]),max(avgCovDens[[2]][[i]]+sdCovDens[[2]][[i]])))

        plot(avgCovDens[[1]][[i]],type="l",frame=F,mian="",xlab="",ylab="",ylim=c(MIN,MAX),col="green",lwd=2,las=1)
        lines(avgCovDens[[1]][[i]]+sdCovDens[[1]][[i]],col="green",lty=2,lwd=0.5)
        lines(avgCovDens[[1]][[i]]-sdCovDens[[1]][[i]],col="green",lty=2,lwd=0.5)

        lines(avgCovDens[[2]][[i]],col="orange",lwd=2)
        lines(avgCovDens[[2]][[i]]+sdCovDens[[2]][[i]],col="orange",lty=2,lwd=0.5)
        lines(avgCovDens[[2]][[i]]-sdCovDens[[2]][[i]],col="orange",lty=2,lwd=0.5)

        IX<- length(myCumSum[[2]][[i]])-500
        abline(v=IX,col="red",lty=2)
        IX<- length(myCumSum[[2]][[i]])-300
        abline(v=IX,col="grey",lty=2)

        mtext(side=1,line=2,text="distance from 3'end",cex=0.7)
        mtext(side=2,line=2,text="cumulative read density",cex=0.7)

        MIN <- min(c(min(avgCumSum[[1]][[i]]-sdCumSum[[1]][[i]]),min(avgCumSum[[2]][[i]]-sdCumSum[[2]][[i]])))
        MAX <- max(c(max(avgCumSum[[1]][[i]]+sdCumSum[[1]][[i]]),max(avgCumSum[[2]][[i]]+sdCumSum[[2]][[i]])))

        plot(avgCumSum[[1]][[i]],type="l",frame=F,mian="",xlab="",ylab="",ylim=c(MIN,MAX),col="green",lwd=2,las=1)
        lines(avgCumSum[[1]][[i]]+sdCovDens[[1]][[i]],col="green",lty=2,lwd=0.5)
        lines(avgCumSum[[1]][[i]]-sdCovDens[[1]][[i]],col="green",lty=2,lwd=0.5)

        lines(avgCumSum[[2]][[i]],col="orange",lwd=2)
        lines(avgCumSum[[2]][[i]]+sdCumSum[[2]][[i]],col="orange",lty=2,lwd=0.5)
        lines(avgCumSum[[2]][[i]]-sdCumSum[[2]][[i]],col="orange",lty=2,lwd=0.5)

        IX<- length(myCumSum[[2]][[i]])-500
        abline(v=IX,col="red",lty=2)
        myd  <- round(avgCumSum[[1]][[i]][IX],digit=0)
        mtext(side=3,line=0,text=paste("cumsum = ",myd," % cov (500)",sep=""),cex=0.7)
        abline(h=myd,col="red",lty=2)

        IX<- length(myCumSum[[2]][[i]])-300
        abline(v=IX,col="grey",lty=2)
        myd  <- round(avgCumSum[[1]][[i]][IX],digit=0)
        mtext(side=3,line=1,text=paste("cumsum = ",myd," % cov (300)",sep=""),cex=0.7)
        abline(h=myd,col="grey",lty=2)

        mtext(side=1,line=2,text="distance from 3'end",cex=0.7)
        mtext(side=2,line=2,text="read coverage density",cex=0.7)
}


par(mfrow=c(3,4))
for(i in c(1:length(avgCumSum[[1]]))){
        MIN <- min(c(min(avgCovDens[[3]][[i]]-sdCovDens[[3]][[i]]),min(avgCovDens[[4]][[i]]-sdCovDens[[4]][[i]])))
        MAX <- max(c(max(avgCovDens[[3]][[i]]+sdCovDens[[3]][[i]]),max(avgCovDens[[4]][[i]]+sdCovDens[[4]][[i]])))

        plot(avgCovDens[[3]][[i]],type="l",frame=F,mian="",xlab="",ylab="",ylim=c(MIN,MAX),col="green",lwd=2,las=1)
        lines(avgCovDens[[3]][[i]]+sdCovDens[[3]][[i]],col="green",lty=2,lwd=0.5)
        lines(avgCovDens[[3]][[i]]-sdCovDens[[3]][[i]],col="green",lty=2,lwd=0.5)

        lines(avgCovDens[[4]][[i]],col="orange",lwd=2)
        lines(avgCovDens[[4]][[i]]+sdCovDens[[4]][[i]],col="orange",lty=2,lwd=0.5)
        lines(avgCovDens[[4]][[i]]-sdCovDens[[4]][[i]],col="orange",lty=2,lwd=0.5)

        IX<- length(myCumSum[[4]][[i]])-500
        abline(v=IX,col="red",lty=2)
        IX<- length(myCumSum[[4]][[i]])-300
        abline(v=IX,col="grey",lty=2)

        mtext(side=1,line=2,text="distance from 3'end",cex=0.7)
        mtext(side=2,line=2,text="cumulative read density",cex=0.7)

        MIN <- min(c(min(avgCumSum[[3]][[i]]-sdCumSum[[3]][[i]]),min(avgCumSum[[4]][[i]]-sdCumSum[[4]][[i]])))
        MAX <- max(c(max(avgCumSum[[3]][[i]]+sdCumSum[[3]][[i]]),max(avgCumSum[[4]][[i]]+sdCumSum[[4]][[i]])))

        plot(avgCumSum[[3]][[i]],type="l",frame=F,mian="",xlab="",ylab="",ylim=c(MIN,MAX),col="green",lwd=2,las=1)
        lines(avgCumSum[[3]][[i]]+sdCovDens[[3]][[i]],col="green",lty=2,lwd=0.5)
        lines(avgCumSum[[3]][[i]]-sdCovDens[[3]][[i]],col="green",lty=2,lwd=0.5)

        lines(avgCumSum[[4]][[i]],col="orange",lwd=2)
        lines(avgCumSum[[4]][[i]]+sdCumSum[[4]][[i]],col="orange",lty=2,lwd=0.5)
        lines(avgCumSum[[4]][[i]]-sdCumSum[[4]][[i]],col="orange",lty=2,lwd=0.5)

        IX<- length(myCumSum[[4]][[i]])-500
        abline(v=IX,col="red",lty=2)
        myd  <- round(avgCumSum[[1]][[i]][IX],digit=0)
        mtext(side=3,line=0,text=paste("cumsum = ",myd," % cov (500)",sep=""),cex=0.7)
        abline(h=myd,col="red",lty=2)

        IX<- length(myCumSum[[4]][[i]])-300
        abline(v=IX,col="grey",lty=2)
        myd  <- round(avgCumSum[[1]][[i]][IX],digit=0)
        mtext(side=3,line=1,text=paste("cumsum = ",myd," % cov (300)",sep=""),cex=0.7)
        abline(h=myd,col="grey",lty=2)

        mtext(side=1,line=2,text="distance from 3'end",cex=0.7)
        mtext(side=2,line=2,text="read coverage density",cex=0.7)
}

dev.off()
```
