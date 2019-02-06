
###                                                                                                 ###
#  NESTED FUNCTIONS                                                                                   #
#                                                                                                     #

RemoveToSmall <- function(myvec,mylim){
     mydist                      <- diff(myvec)
     while(sum(mydist<mylim)>0){
            myvec <- myvec[-which(mydist<mylim)[1]]
            mydist <- diff(myvec)
      }
     return(myvec)
}

#limDisit=minimum distance between pA site; always take the longest when merging.
mergePAsTissue <- function(tempTX=levels(groups)[1],method="ebs",limDist=50,windowsize=winD,NumSeg=noSeg,doLog=isLog,MinCov=minCov){
    rm("new.isoforms","soi","mystart","myend","myPa","tempres","mystrd")
    soi                     <- ids[which(ids[,1]==tempTX),2]
    mystart                 <- 0
    myend                   <- 0
    myPa                    <- list()
    myPaf                   <- list()
    k=1#maximum 16 given that 4 different samples but raw and log2 count
    for(i in c(1:length(soi))){
        print(i)
        tempres                 <-  myres[[match(soi[i],names(myres))]]#select all coverage which correspond to txID
        myix                    <-  match(tempTX,tempres$txID)
        if(method=="SEG"){
                            # wd   = width of window for running median algo
                            # Llim = minimal distance to consider between 2 consecutive fragments
                            # Km   = max no. segments
                            # MIN  =  minimum coverage to consider; otherwise set to zero
                            myPa[[k]]<-  getpASeg(myobj=tempres,Ix=myix,wd=windowsize,Llim=limDist,Km=NumSeg,MIN=MinCov,LOG=doLog)
                            if(!is.na(myPa[[k]]$new.isoforms)){
                                tokeep                  <-  width(myPa[[k]]$new.isoforms)<10^5
                                myPa[[k]]$new.isoforms  <-  myPa[[k]]$new.isoforms[tokeep,]
                            }

                            if(length(myPa[[k]]$new.isoforms)==0|is.na(myPa[[k]]$new.isoforms)){
                                myPa[[k]]$new.isoforms    <- g3utr[match(tempTX,as.character(mcols(g3utr)$X.transcript_id)),]
                            }
                            k=k+1

                            #Je devrai filtrer les samples ici; enlever les voisins trop proches; filterer entre les k
                            if(length(myPaf)==0){myPaf <- myPa[[1]]$new.isoforms}
                            if(k>1){

                                     mystrd                  <-  tempres$strand[myix]
                                        L1                      <-  myPa[[k-1]]$new.isoforms
                                        L2                      <-  myPaf
                                        if(mystrd=="+"){
                                            start(L1) <-  end(L1)-1
                                            start(L2) <-  end(L2)-1
                                        }
                                        if(mystrd=="-"){
                                            end(L1)<- start(L1)+1
                                            end(L2)<- start(L2)+1
                                        }

                                         test     <-  as.data.frame(distanceToNearest(L1,L2,ignore.strand=FALSE))
                                         sel1     <-  test$distance<=windowsize
                                         sel2     <-  width(myPa[[k-1]]$new.isoforms)[test$queryHits]>width(myPaf)[test$subjectHits]
                                         tokeep1  <-  sel1&sel2|!sel1 #either keep the longest when 2 neighbours OR when not close
                                         tokeep2  <-  !(c(1:length(myPaf))%in%test$subjectHits[sel1&sel2])

                                         myPaf      <- GRanges( seqnames  = Rle(chr,sum(c(tokeep1,tokeep2))),
                                                                ranges    = IRanges(
                                                                    start    = c(start(myPaf[tokeep2]),start(myPa[[k-1]]$new.isoforms[tokeep1])),
                                                                    end      = c(end(myPaf[tokeep2]),end(myPa[[k-1]]$new.isoforms[tokeep1]))
                                                                    ),
                                                                    strand   = Rle(mystrd, sum(c(tokeep1,tokeep2)))
                                                               )
                            }
                            if(doLog){
                                myPa[[k]]<-  getpASeg(myobj=tempres,Ix=myix,wd=windowsize,Llim=limDist,Km=NumSeg,MIN=MinCov,LOG=FALSE)

                                if(!is.na(myPa[[k]]$new.isoforms)){
                                    tokeep                  <-  width(myPa[[k]]$new.isoforms)<10^5
                                    myPa[[k]]$new.isoforms  <-  myPa[[k]]$new.isoforms[tokeep,]
                                }


                                if(length(myPa[[k]]$new.isoforms)==0|is.na(myPa[[k]]$new.isoforms)){
                                    myPa[[k]]$new.isoforms    <- g3utr[match(tempTX,as.character(mcols(g3utr)$X.transcript_id)),]
                                }

                                k=k+1

                                #Je devrai filtrer les samples ici; enlever les voisins trop proches; filterer entre les k
                                if(k>1){
                                        mystrd                  <-  tempres$strand[myix]
                                        L1                      <-  myPa[[k-1]]$new.isoforms
                                        L2                      <-  myPaf
                                        if(mystrd=="+"){
                                            start(L1) <-  end(L1)-1
                                            start(L2) <-  end(L2)-1
                                        }
                                        if(mystrd=="-"){
                                            end(L1)<- start(L1)+1
                                            end(L2)<- start(L2)+1
                                        }

                                         test     <-  as.data.frame(distanceToNearest(L1,L2,ignore.strand=FALSE))
                                         sel1     <-  test$distance<=windowsize
                                         sel2     <-  width(myPa[[k-1]]$new.isoforms)[test$queryHits]>width(myPaf)[test$subjectHits]
                                         tokeep1  <-  sel1&sel2|!sel1 #either keep the longest when 2 neighbours OR when not close
                                         tokeep2  <-  !(c(1:length(myPaf))%in%test$subjectHits[sel1&sel2])


                                         myPaf      <- GRanges( seqnames  = Rle(chr,sum(c(tokeep1,tokeep2))),
                                                                ranges    = IRanges(
                                                                    start    = c(start(myPaf[tokeep2]),start(myPa[[k-1]]$new.isoforms[tokeep1])),
                                                                    end      = c(end(myPaf[tokeep2]),end(myPa[[k-1]]$new.isoforms[tokeep1]))
                                                                    ),
                                                                    strand   = Rle(mystrd, sum(c(tokeep1,tokeep2)))
                                                               )

                                }
                            }


        print(width(myPaf))
        }
        if(method=="EBS"){#Beaucoup trop lent cette routine!
                            # wd   = width of window for running median algo
                            # Llim = minimal distance to consider between 2 consecutive fragments
                            # Km   = max no. segments
                            # MIN  =  minimum coverage to consider; otherwise set to zero
                            myPa[[k]]<-  getpAEBS(myobj=tempres,Ix=myix,wd=windowsize,Llim=limDist,Km=NumSeg,MIN=MinCov,LOG=doLog)
                            if(!is.na(myPa[[k]]$new.isoforms)){
                                tokeep                  <-  width(myPa[[k]]$new.isoforms)<10^5
                                myPa[[k]]$new.isoforms  <-  myPa[[k]]$new.isoforms[tokeep,]
                            }

                            if(length(myPa[[k]]$new.isoforms)==0|is.na(myPa[[k]]$new.isoforms)){
                                myPa[[k]]$new.isoforms    <- g3utr[match(tempTX,as.character(mcols(g3utr)$X.transcript_id)),]
                            }

                            k=k+1

                            if(doLog){
                                myPa[[k]]<-  getpAEBS(myobj=tempres,Ix=myix,wd=windowsize,Llim=limDist,Km=NumSeg,MIN=MinCov,LOG=FALSE)

                                if(!is.na(myPa[[k]]$new.isoforms)){
                                    tokeep                  <-  width(myPa[[k]]$new.isoforms)<10^5
                                    myPa[[k]]$new.isoforms  <-  myPa[[k]]$new.isoforms[tokeep,]
                                }


                                if(length(myPa[[k]]$new.isoforms)==0|is.na(myPa[[k]]$new.isoforms)){
                                    myPa[[k]]$new.isoforms    <- g3utr[match(tempTX,as.character(mcols(g3utr)$X.transcript_id)),]
                                }

                                k=k+1
                            }
        }
        mystrd                  <-  tempres$strand[myix]
        myend                   <-  max(myend,tempres$end[myix])
        mystart                 <-  min(myend,tempres$start[myix])
    }

    print("computed myPAs")
    if(mystrd=="+"){
        print("I am +")



        myends                      <- sort(unique(end(myPaf)))
        new.isoforms                <- GRanges(seqnames = Rle(chr,length(myends)),
                                               ranges   = IRanges(start=rep(mystart,length(myends)),end=myends),
                                               strand   = Rle(mystrd, length(myends))
                                                )
        mcols(new.isoforms)$isoform <- c(1:length(myends))
        names(new.isoforms)         <- paste(tempTX,c(1:length(myends)),sep=".")
        mcols(new.isoforms)$txID    <- rep(tempTX,length(myends))
    }

    if(mystrd=="-"){
        print("I am -")
        mystart                     <- sort(unique(start(myPaf)))
        Ord                         <- rev(c(1:length(mystart)))
        new.isoforms                <- GRanges(seqnames = Rle(chr,length(mystart)),
                                               ranges   = IRanges(start=mystart[Ord],end=rep(myend,length(mystart))),
                                               strand   = Rle(mystrd, length(mystart))
                                                )
        mcols(new.isoforms)$isoform <- c(1:length(mystart))
        names(new.isoforms)         <- paste(tempTX,c(1:length(mystart)),sep=".")
        mcols(new.isoforms)$txID    <- rep(tempTX,length(mystart))
    }
    rm(method,limDist,windowsize,NumSeg,doLog,MinCov,myPA)
    return(new.isoforms)
}



getpASeg <- function(myobj=myres[[1]],Ix=3,wd=75,Llim=50,Km=10,MIN=7,LOG=TRUE){
    #Select data of interest
    myvalues <- myobj$cov[[Ix]]
    myID     <- myobj$txID[Ix]
    print(myID)
    myrange  <- myUTR[which(as.character(mcols(myUTR)$X.transcript_id.)==myID),]
    if(length(myrange)>1){
        myrange <- myrange[which(width(myrange)==length(myvalues))[1],]
    }
    strd     <- as.character(strand(myrange))
    vec      <- myvalues
    L        <- length(vec)

    #Check which of these extensions contain introns
    gOver    <- findOverlaps(query=myrange,subject=g3introns)
    idxU     <- queryHits(gOver)
    idxG     <- subjectHits(gOver)
    test     <- length(gOver)>0
    if(test){
        test <- width(g3introns[idxG,])<width(myrange)
    }

    if(!test&length(vec)>10){
        if(strd=="-"){vec <- vec[c(length(vec):1)]}
        #Compute running median
        if(LOG){
                mymed  <- as.integer(log2(1+as.vector(runmed(x=vec,k=min(wd,length(vec)-10),endrule="median"))))
                myLim  <- log2(MIN)
                mymed[mymed<=MIN] <- 0
        }
        if(!LOG){
                mymed  <- as.integer(as.vector(runmed(x=vec,k=min(wd,length(vec)-10),endrule="median")))
                myLim  <- MIN
                mymed[mymed<=MIN] <- 0
        }
        maxK   <- length(rle(mymed)$lengths)
        Km     <- min(maxK,Km)

        if(sum(mymed>myLim)>Llim&Km>1){
            #Segment
            Seg     <- Segmentor(mymed,model=3,Kmax=Km)
            Kchoose <- SelectModel(Seg, penalty="oracle")
            brks    <- as.vector(getBreaks(Seg)[Kchoose, 1:Kchoose])
            brks    <- unique(c(1,brks,length(mymed)))#the first and the last must always stay

#pdf("test.pdf")
#plot(mymed,type="l")
#lines(vec,col="green")
#abline(v=brks,col="red")
#abline(v=pAs,col="blue")
#dev.off()

            #brks    <- unique(c(1,brks))
            mydist  <- diff(brks)
            while(sum(mydist<Llim)>0){
                    torm   <- which(mydist<Llim)[1]
                    if(torm==1){torm<-2}
                    if(!is.na(torm)){
                        brks   <- brks[-torm]
                        mydist <- diff(brks)
                        }
            }
            brks   <- unique(c(brks,length(mymed)))
            #check that last segment is covered; if not remove it!
            no.brks <- length(brks)
            if(median(vec[c(brks[no.brks-1]:brks[no.brks])])<=myLim&(no.brks>1)){
                brks <- brks[-no.brks]
            }

            if(length(brks)<2){
                #print(paste("length brks=",length(brks),"easy!"))
                new.isoforms   <- NA
                pAs            <- 0

            }
            if(length(brks)>=2){
                #print(paste("length brks=",length(brks),"more difficult!"))
                mysegs  <- list()
                for(i in c(2:length(brks))){
                 mysegs[[i]] <- c(brks[(i-1)]:(brks[i]))
                }
                mysegs <- mysegs[-1]
                #Compute median of each segment
                mym    <-  unlist(lapply(mysegs,function(x){return(median(mymed[x]))}))
                mym[mym<=myLim]<-0
                #Identify pAs
                if(length(which(mym>myLim))>0){
                    sel <- union(which(diff(mym)<0),max(which(mym>myLim)))
                    pAs <- sapply(sel,function(x)return(max(mysegs[[x]])))
                    #Create real range
                    new.isoforms                <- myrange[rep(1,length(pAs)),]
                    mcols(new.isoforms)$isoform <- c(1:length(pAs))
                    names(new.isoforms)         <- paste(myID,c(1:length(pAs)),sep=".")
                    if(strd=="+"){end(new.isoforms)  <- start(new.isoforms)+pAs-1}
                    if(strd=="-"){start(new.isoforms)<- end(new.isoforms)-pAs-1}
                }
                if(length(which(mym>MIN))==0){
                        new.isoforms   <- NA
                        pAs            <- 0
                }
            }
        }
        if(Km==1){
            new.isoforms   <- myrange
            pAs            <- 0
        }
        if(sum(mymed>10)<=Llim){
            new.isoforms   <- NA
            pAs            <- 0
        }
            return(list(new.isoforms=new.isoforms,pAs=pAs))
    }

    #CASE with INTRONS
    if(test&length(vec)>10){
    ##                                                                          ##
    #Then should not account for the intron                                      #
    #
    if(strd=="+"){
        #Check that first intron is the 'first'
        mytemp <- end(g3introns[idxG,])
        ord    <- sort(mytemp,index.return=T,decreasing=F)$ix
        idxG   <- idxG[ord]
        #Compute relative start and end of intron relative to the end of the 3' UTR
        startR <- start(g3introns[idxG,]) - start(myrange)+1
        endR   <- end(g3introns[idxG,]) - start(myrange)+1
        wR     <- width(g3introns[idxG,])
        #Regions to remove
        torm <- do.call(args=lapply(c(1:length(idxG)),function(x)return(seq(from=startR[x],to=endR[x]))),what=c)
        vec  <- vec[-torm]
        #Identify position of introns along new shortened 3' UTR
        ri     <- rep(0,length(idxG))
        for(i in c(1:length(idxG))){
            if(i==1){
                ri[i]<- startR[i]
            }
            if(i>1){
                ri[i]<- startR[i] - sum(wR[c(1:i-1)])
            }
        }
     }

    if(strd=="-"){
        #Check that first intron is the 'last'
        mytemp <- end(g3introns[idxG,])
        ord    <- sort(mytemp,index.return=T,decreasing=T)$ix
        idxG   <- idxG[ord]
        #Compute relative start and end of intron relative to the end of the 3' UTR
        startR <- end(myrange)-end(g3introns[idxG,]) +1 #distance from start of the 3' UTR
        endR   <- end(myrange)-start(g3introns[idxG,])+1
        wR     <- width(g3introns[idxG,])
        #Inverse vector
        vec    <- vec[c(length(vec):1)]
        #Regions to remove
        torm <- do.call(args=lapply(c(1:length(idxG)),function(x)return(seq(from=startR[x],to=endR[x]))),what=c)
        vec  <- vec[-torm]
        #Identify position of introns along new shortened 3' UTR
        ri     <- rep(0,length(idxG))
        for(i in c(1:length(idxG))){
            if(i==1){
                ri[i]<- startR[i]
            }
            if(i>1){
                ri[i]<- startR[i] - sum(wR[c(1:i-1)])
            }
        }
     }
        #Compute running median
        if(LOG){
                mymed  <- as.integer(log2(1+as.vector(runmed(x=vec,k=min(wd,length(vec)-10),endrule="median"))))
                myLim  <- log2(MIN)
                mymed[mymed<=MIN] <- 0
        }
        if(!LOG){
                mymed  <- as.integer(as.vector(runmed(x=vec,k=min(wd,length(vec)-10),endrule="median")))
                myLim  <- MIN
                mymed[mymed<=MIN] <- 0

        }
        maxK   <- length(rle(mymed)$lengths)
        Km     <- min(maxK,Km)

        if(sum(mymed>myLim)>Llim&Km>1){
            #Segment
            Seg     <- Segmentor(mymed,model=3,Kmax=Km)
            Kchoose <- SelectModel(Seg, penalty="oracle")
            brks    <- as.vector(getBreaks(Seg)[Kchoose, 1:Kchoose])
            #brks    <- unique(c(1,brks))
            brks    <- unique(c(1,brks,length(mymed)))#the first and the last must always stay
            no.brks <- length(brks)
            if(median(vec[c(brks[no.brks-1]:brks[no.brks])])<=myLim&(no.brks>1)){
                brks <- brks[-no.brks]
            }

            mydist  <- diff(brks)
            while(sum(mydist<Llim)>0){
                    torm   <- which(mydist<Llim)[1]
                    if(torm==1){torm<-2}
                    if(!is.na(torm)){
                        brks   <- brks[-torm]
                        mydist <- diff(brks)
                        }
            }
            brks   <- unique(c(brks,length(mymed)))
            if(length(brks)<2){
                        new.isoforms   <- NA
                        pAs            <- 0
            }
            if(length(brks)>=2){
                mysegs  <- list()
                for(i in c(2:length(brks))){
                    mysegs[[i]] <- c(brks[(i-1)]:(brks[i]))
                }
                mysegs <- mysegs[-1]
                #Compute median of each segment
                mym    <-  unlist(lapply(mysegs,function(x){return(median(mymed[x]))}))
                mym[mym<=myLim]<-0

                #Identify pAs
                if(length(which(mym>myLim))>0){
                    sel <- union(which(diff(mym)<0),max(which(mym>myLim)))
                    pAs <- sapply(sel,function(x)return(max(mysegs[[x]])))#Distance from the start of the 3' UTR (therefore the end coordinate when negative)
                    #Correct pAs with considering the removed part due to introns
                    for(k in c(1:length(pAs))){
                          ixoi <- which(ri<pAs[k])
                          if(sum(ri<pAs[k])>0){
                            pAs[k] <- pAs[k]+sum(wR[ixoi])
                          }
                    }

                    #Create real range
                    #new.isoforms                <- myrange[rep(1,length(pAs)),]
                    #mcols(new.isoforms)$isoform <- c(1:length(pAs))
                    #names(new.isoforms)         <- paste(myID,c(1:length(pAs)),sep=".")
                    #if(strd=="+"){end(new.isoforms)  <- start(new.isoforms)+pAs-1}
                    #if(strd=="-"){start(new.isoforms)<- end(new.isoforms)-pAs-1}

                    if(strd=="+"){
                        #Create real range
                        new.isoforms                <- myrange[rep(1,length(pAs)),]
                        mcols(new.isoforms)$isoform <- c(1:length(pAs))
                        names(new.isoforms)         <- paste(myID,c(1:length(pAs)),sep=".")
                        end(new.isoforms)           <- start(new.isoforms)+pAs-1
                    }
                    if(strd=="-"){#There is a problem HERE!
                        #Create real range
                        new.isoforms                <- myrange[rep(1,length(pAs)),]#I here take the range of
                        mcols(new.isoforms)$isoform <- c(1:length(pAs))
                        names(new.isoforms)         <- paste(myID,c(1:length(pAs)),sep=".")
                        start(new.isoforms)         <- end(new.isoforms)-pAs+1
                    }
                }
                if(length(which(mym>MIN))==0){
                        new.isoforms   <- NA
                        pAs            <- 0
                }
            }
        }
        if(Km<2){
            new.isoforms   <- myrange
            pAs            <- 0
        }
        if(sum(mymed>10)<Llim){
            new.isoforms   <- NA
            pAs            <- 0
        }
        return(list(new.isoforms=new.isoforms,pAs=pAs))
    }
    if(length(vec)<=10){
            new.isoforms   <- myrange
            pAs            <- length(myvalues)
            return(list(new.isoforms=new.isoforms,pAs=pAs))
    }

}

getpAEBS <- function(myobj=myres[[1]],Ix=3,wd=75,Llim=50,Km=10,MIN=7,LOG=TRUE){
    #Select data of interest
    myvalues <- myobj$cov[[Ix]]
    myID     <- myobj$txID[Ix]
    print(myID)
    myrange  <- myUTR[which(as.character(mcols(myUTR)$X.transcript_id.)==myID),]
    if(length(myrange)>1){
        myrange <- myrange[which(width(myrange)==length(myvalues))[1],]
    }
    strd     <- as.character(strand(myrange))
    vec      <- myvalues
    L        <- length(vec)

    #Check which of these extensions contain introns
    gOver    <- findOverlaps(query=myrange,subject=g3introns)
    idxU     <- queryHits(gOver)
    idxG     <- subjectHits(gOver)
    test     <- length(gOver)>0
    if(test){
        test <- width(g3introns[idxG,])<width(myrange)
    }

    if(!test&length(vec)>10){
        if(strd=="-"){vec <- vec[c(length(vec):1)]}
        #Compute running median
        if(LOG){
                mymed  <- as.integer(log2(1+as.vector(runmed(x=vec,k=min(wd,length(vec)-10),endrule="median"))))
                myLim  <- log2(MIN)
                mymed[mymed<=MIN] <- 0
        }
        if(!LOG){
                mymed  <- as.integer(as.vector(runmed(x=vec,k=min(wd,length(vec)-10),endrule="median")))
                myLim  <- MIN
                mymed[mymed<=MIN] <- 0
        }
        maxK   <- length(rle(mymed)$lengths)
        Km     <- min(maxK,Km)

        if(sum(mymed>myLim)>Llim&Km>1){
            #Segment
            out   <- EBSegmentation(mymed,Kmax=Km)
            mybic <- EBSBIC(out)
            bic   <- (EBSBIC(out)$NbBIC)-1
            GetHighestDensity <- function(ki){
                yi<-EBSDistrib(out,k=ki,Kk=bic+1)
                return(which(yi==max(yi)))
            }
            mypos  <- unlist(lapply(c(1:bic),GetHighestDensity))
            brks   <- unique(c(1,mypos,length(mymed)))
#pdf("test.pdf")
#plot(mymed,type="l")
#lines(vec,col="green")
#abline(v=brks,col="red")
#abline(v=pAs,col="blue")
#dev.off()

            #brks    <- unique(c(1,brks))
            mydist  <- diff(brks)
            while(sum(mydist<Llim)>0){
                    torm   <- which(mydist<Llim)[1]
                    if(torm==1){torm<-2}
                    if(!is.na(torm)){
                        brks   <- brks[-torm]
                        mydist <- diff(brks)
                        }
            }
            brks   <- unique(c(brks,length(mymed)))
            #check that last segment is covered; if not remove it!
            no.brks <- length(brks)
            if(median(vec[c(brks[no.brks-1]:brks[no.brks])])<=myLim&(no.brks>1)){
                brks <- brks[-no.brks]
            }

            if(length(brks)<2){
                #print(paste("length brks=",length(brks),"easy!"))
                new.isoforms   <- NA
                pAs            <- 0

            }
            if(length(brks)>=2){
                #print(paste("length brks=",length(brks),"more difficult!"))
                mysegs  <- list()
                for(i in c(2:length(brks))){
                 mysegs[[i]] <- c(brks[(i-1)]:(brks[i]))
                }
                mysegs <- mysegs[-1]
                #Compute median of each segment
                mym    <-  unlist(lapply(mysegs,function(x){return(median(mymed[x]))}))
                mym[mym<=myLim]<-0
                #Identify pAs
                if(length(which(mym>myLim))>0){
                    sel <- union(which(diff(mym)<0),max(which(mym>myLim)))
                    pAs <- sapply(sel,function(x)return(max(mysegs[[x]])))
                    #Create real range
                    new.isoforms                <- myrange[rep(1,length(pAs)),]
                    mcols(new.isoforms)$isoform <- c(1:length(pAs))
                    names(new.isoforms)         <- paste(myID,c(1:length(pAs)),sep=".")
                    if(strd=="+"){end(new.isoforms)  <- start(new.isoforms)+pAs-1}
                    if(strd=="-"){start(new.isoforms)<- end(new.isoforms)-pAs-1}
                }
                if(length(which(mym>MIN))==0){
                        new.isoforms   <- NA
                        pAs            <- 0
                }
            }
        }
        if(Km==1){
            new.isoforms   <- myrange
            pAs            <- 0
        }
        if(sum(mymed>10)<=Llim){
            new.isoforms   <- NA
            pAs            <- 0
        }
            return(list(new.isoforms=new.isoforms,pAs=pAs))
    }

    #CASE with INTRONS
    if(test&length(vec)>10){
    ##                                                                          ##
    #Then should not account for the intron                                      #
    #
    if(strd=="+"){
        #Check that first intron is the 'first'
        mytemp <- end(g3introns[idxG,])
        ord    <- sort(mytemp,index.return=T,decreasing=F)$ix
        idxG   <- idxG[ord]
        #Compute relative start and end of intron relative to the end of the 3' UTR
        startR <- start(g3introns[idxG,]) - start(myrange)+1
        endR   <- end(g3introns[idxG,]) - start(myrange)+1
        wR     <- width(g3introns[idxG,])
        #Regions to remove
        torm <- do.call(args=lapply(c(1:length(idxG)),function(x)return(seq(from=startR[x],to=endR[x]))),what=c)
        vec  <- vec[-torm]
        #Identify position of introns along new shortened 3' UTR
        ri     <- rep(0,length(idxG))
        for(i in c(1:length(idxG))){
            if(i==1){
                ri[i]<- startR[i]
            }
            if(i>1){
                ri[i]<- startR[i] - sum(wR[c(1:i-1)])
            }
        }
     }

    if(strd=="-"){
        #Check that first intron is the 'last'
        mytemp <- end(g3introns[idxG,])
        ord    <- sort(mytemp,index.return=T,decreasing=T)$ix
        idxG   <- idxG[ord]
        #Compute relative start and end of intron relative to the end of the 3' UTR
        startR <- end(myrange)-end(g3introns[idxG,]) +1 #distance from start of the 3' UTR
        endR   <- end(myrange)-start(g3introns[idxG,])+1
        wR     <- width(g3introns[idxG,])
        #Inverse vector
        vec    <- vec[c(length(vec):1)]
        #Regions to remove
        torm <- do.call(args=lapply(c(1:length(idxG)),function(x)return(seq(from=startR[x],to=endR[x]))),what=c)
        vec  <- vec[-torm]
        #Identify position of introns along new shortened 3' UTR
        ri     <- rep(0,length(idxG))
        for(i in c(1:length(idxG))){
            if(i==1){
                ri[i]<- startR[i]
            }
            if(i>1){
                ri[i]<- startR[i] - sum(wR[c(1:i-1)])
            }
        }
     }
        #Compute running median
        if(LOG){
                mymed  <- as.integer(log2(1+as.vector(runmed(x=vec,k=min(wd,length(vec)-10),endrule="median"))))
                myLim  <- log2(MIN)
                mymed[mymed<=MIN] <- 0
        }
        if(!LOG){
                mymed  <- as.integer(as.vector(runmed(x=vec,k=min(wd,length(vec)-10),endrule="median")))
                myLim  <- MIN
                mymed[mymed<=MIN] <- 0

        }
        maxK   <- length(rle(mymed)$lengths)
        Km     <- min(maxK,Km)

        if(sum(mymed>myLim)>Llim&Km>1){
            #Segment
            out   <- EBSegmentation(mymed,Kmax=Km)
            mybic <- EBSBIC(out)
            bic   <- (EBSBIC(out)$NbBIC)-1
            GetHighestDensity <- function(ki){
                yi<-EBSDistrib(out,k=ki,Kk=bic+1)
                return(which(yi==max(yi)))
            }
            mypos  <- unlist(lapply(c(1:bic),GetHighestDensity))
            brks   <- c(1,mypos,length(mymed))#the first and the last must always stay

            no.brks <- length(brks)
            if(median(vec[c(brks[no.brks-1]:brks[no.brks])])<=myLim&(no.brks>1)){
                brks <- brks[-no.brks]
            }

            mydist  <- diff(brks)
            while(sum(mydist<Llim)>0){
                    torm   <- which(mydist<Llim)[1]
                    if(torm==1){torm<-2}
                    if(!is.na(torm)){
                        brks   <- brks[-torm]
                        mydist <- diff(brks)
                        }
            }
            brks   <- unique(c(brks,length(mymed)))
            if(length(brks)<2){
                        new.isoforms   <- NA
                        pAs            <- 0
            }
            if(length(brks)>=2){
                mysegs  <- list()
                for(i in c(2:length(brks))){
                    mysegs[[i]] <- c(brks[(i-1)]:(brks[i]))
                }
                mysegs <- mysegs[-1]
                #Compute median of each segment
                mym    <-  unlist(lapply(mysegs,function(x){return(median(mymed[x]))}))
                mym[mym<=myLim]<-0

                #Identify pAs
                if(length(which(mym>myLim))>0){
                    sel <- union(which(diff(mym)<0),max(which(mym>myLim)))
                    pAs <- sapply(sel,function(x)return(max(mysegs[[x]])))#Distance from the start of the 3' UTR (therefore the end coordinate when negative)
                    #Correct pAs with considering the removed part due to introns
                    for(k in c(1:length(pAs))){
                          ixoi <- which(ri<pAs[k])
                          if(sum(ri<pAs[k])>0){
                            pAs[k] <- pAs[k]+sum(wR[ixoi])
                          }
                    }

                    #Create real range
                    #new.isoforms                <- myrange[rep(1,length(pAs)),]
                    #mcols(new.isoforms)$isoform <- c(1:length(pAs))
                    #names(new.isoforms)         <- paste(myID,c(1:length(pAs)),sep=".")
                    #if(strd=="+"){end(new.isoforms)  <- start(new.isoforms)+pAs-1}
                    #if(strd=="-"){start(new.isoforms)<- end(new.isoforms)-pAs-1}

                    if(strd=="+"){
                        #Create real range
                        new.isoforms                <- myrange[rep(1,length(pAs)),]
                        mcols(new.isoforms)$isoform <- c(1:length(pAs))
                        names(new.isoforms)         <- paste(myID,c(1:length(pAs)),sep=".")
                        end(new.isoforms)           <- start(new.isoforms)+pAs-1
                    }
                    if(strd=="-"){#There is a problem HERE!
                        #Create real range
                        new.isoforms                <- myrange[rep(1,length(pAs)),]#I here take the range of
                        mcols(new.isoforms)$isoform <- c(1:length(pAs))
                        names(new.isoforms)         <- paste(myID,c(1:length(pAs)),sep=".")
                        start(new.isoforms)         <- end(new.isoforms)-pAs+1
                    }
                }
                if(length(which(mym>MIN))==0){
                        new.isoforms   <- NA
                        pAs            <- 0
                }
            }
        }
        if(Km<2){
            new.isoforms   <- myrange
            pAs            <- 0
        }
        if(sum(mymed>10)<Llim){
            new.isoforms   <- NA
            pAs            <- 0
        }
        return(list(new.isoforms=new.isoforms,pAs=pAs))
    }
    if(length(vec)<=10){
            new.isoforms   <- myrange
            pAs            <- length(myvalues)
            return(list(new.isoforms=new.isoforms,pAs=pAs))
    }

}


getpASegOld <- function(myobj=myres[[1]],Ix=3,wd=75,Llim=50,Km=10,MIN=7,LOG=TRUE){
    #Select data of interest
    
    myvalues <- myobj$cov[[Ix]]
    myvalneg <- myobj$covn[[Ix]]
    if(length(myvalneg)>0){
        print("I will correct for opposite strand reads")
        myvalues            <- myvalues-myvalneg
        myvalues[myvalues<0]<- 0
    }
    myID     <- myobj$txID[Ix]
    myrange  <- myUTR[which(as.character(mcols(myUTR)$X.transcript_id.)==myID),]
    if(length(myrange)>1){
        myrange <- myrange[which(width(myrange)==length(myvalues))[1],]
    }
    strd     <- as.character(strand(myrange))
    vec      <- myvalues
    L        <- length(vec)

    #Check which of these extensions contain introns
    gOver    <- findOverlaps(query=myrange,subject=g3introns)
    idxU     <- queryHits(gOver)
    idxG     <- subjectHits(gOver)
    test     <- length(gOver)>0
    if(test){
        test <- width(g3introns[idxG,])<width(myrange)
    }

    if(!test&length(vec)>10){
        if(strd=="-"){vec <- vec[c(length(vec):1)]}
        #Compute running median
        if(LOG){
                mymed  <- as.integer(log2(1+as.vector(runmed(x=vec,k=min(wd,length(vec)-10),endrule="median"))))
                myLim  <- log2(MIN)
        }
        if(!LOG){
                mymed  <- as.integer(as.vector(runmed(x=vec,k=min(wd,length(vec)-10),endrule="median")))
                myLim  <- MIN
        }
        maxK   <- length(rle(mymed)$lengths)
        Km     <- min(maxK,Km)

        if(sum(mymed>myLim)>Llim&Km>1){
            #Segment
            Seg     <- Segmentor(mymed,model=3,Kmax=Km)
            Kchoose <- SelectModel(Seg, penalty="oracle")
            brks    <- as.vector(getBreaks(Seg)[Kchoose, 1:Kchoose])
            brks    <- unique(c(1,brks,length(mymed)))#the first and the last must always stay
            #brks    <- unique(c(1,brks))
            mydist  <- diff(brks)
            while(sum(mydist<Llim)>0){
                    torm   <- which(mydist<Llim)[1]
                    if(torm==1){torm<-2}
                    brks   <- brks[-torm]
                    mydist <- diff(brks)
            }
            brks   <- unique(c(brks,length(mymed)))
            if(length(brks)==2){
                #print(paste("length brks=",length(brks),"easy!"))
                new.isoforms   <- myrange
                pAs            <- length(myvalues)
            }
            if(length(brks)>2){
                #print(paste("length brks=",length(brks),"more difficult!"))
                mysegs  <- list()
                for(i in c(2:length(brks))){
                 mysegs[[i]] <- c(brks[(i-1)]:(brks[i]))
                }
                mysegs <- mysegs[-1]
                #Compute median of each segment
                mym    <-  unlist(lapply(mysegs,function(x){return(mean(mymed[x]))}))
                #Identify pAs
                if(LOG){MIN<-log2(MIN)}
                if(length(which(mym>MIN))>0){
                    sel <- union(which(diff(mym)<0),max(which(mym>MIN)))
                    pAs <- sapply(sel,function(x)return(max(mysegs[[x]])))
                    #Create real range
                    new.isoforms                <- myrange[rep(1,length(pAs)),]
                    mcols(new.isoforms)$isoform <- c(1:length(pAs))
                    names(new.isoforms)         <- paste(myID,c(1:length(pAs)),sep=".")
                    if(strd=="+"){end(new.isoforms)  <- start(new.isoforms)+pAs-1}
                    if(strd=="-"){start(new.isoforms)<- end(new.isoforms)-pAs-1}
                }
                if(length(which(mym>MIN))==0){
                        new.isoforms   <- NA
                        pAs            <- 0
                }
            }
        }
        if(Km==1){
            new.isoforms   <- myrange
            pAs            <- 0
        }
        if(sum(mymed>10)<=Llim){
            new.isoforms   <- NA
            pAs            <- 0
        }
            return(list(new.isoforms=new.isoforms,pAs=pAs))
    }

    #CASE with INTRONS
    if(test&length(vec)>10){
    ##                                                                          ##
    #Then should not account for the intron                                      #
    #                                                                            #
        startR <- start(g3introns[idxG,]) - start(myrange)+1
        endR   <- end(g3introns[idxG,]) - start(myrange)+1
        wR     <- endR-startR+1
        ri     <- rep(0,length(endR))
        for(i in c(1:length(startR))){
            if(i>1){
                ri[i] <- startR[i]-sum(wR[c(1:(i-1))])
                vec   <- vec[-seq(from=ri[i],to=(endR[i]-sum(wR[c(1:(i-1))])))]
            }
            if(i==1){
                ri[i] <- startR[i]
                vec   <- vec[-seq(from=ri[i],to=endR[i])]
            }
        }

        if(strd=="-"){vec <- vec[c(length(vec):1)]}
        #Compute running median
        if(LOG){
                mymed  <- as.integer(log2(1+as.vector(runmed(x=vec,k=min(wd,length(vec)-10),endrule="median"))))
                myLim  <- log2(MIN)
        }
        if(!LOG){
                mymed  <- as.integer(as.vector(runmed(x=vec,k=min(wd,length(vec)-10),endrule="median")))
                myLim  <- MIN
        }
        maxK   <- length(rle(mymed)$lengths)
        Km     <- min(maxK,Km)
        if(LOG){MIN<-log2(MIN)}
        if(sum(mymed>myLim)>Llim&Km>1){
            #Segment
            Seg     <- Segmentor(mymed,model=3,Kmax=Km)
            Kchoose <- SelectModel(Seg, penalty="oracle")
            brks    <- as.vector(getBreaks(Seg)[Kchoose, 1:Kchoose])
            #brks    <- unique(c(1,brks))
            brks    <- unique(c(1,brks,length(mymed)))#the first and the last must always stay
            mydist  <- diff(brks)
            while(sum(mydist<Llim)>0){
                    torm   <- which(mydist<Llim)[1]
                    if(torm==1){torm<-2}
                    brks   <- brks[-torm]
                    mydist <- diff(brks)
            }
            brks   <- unique(c(brks,length(mymed)))
            if(length(brks)==2){
                new.isoforms   <- myrange
                pAs            <- length(myvalues)
            }
            if(length(brks)>2){
                mysegs  <- list()
                for(i in c(2:length(brks))){
                    mysegs[[i]] <- c(brks[(i-1)]:(brks[i]))
                }
                mysegs <- mysegs[-1]
                #Compute median of each segment
                mym    <-  unlist(lapply(mysegs,function(x){return(mean(mymed[x]))}))
                #Identify pAs

                if(length(which(mym>MIN))>0){
                    sel <- union(which(diff(mym)<0),max(which(mym>MIN)))
                    pAs <- sapply(sel,function(x)return(max(mysegs[[x]])))
                    #Create real range
                    new.isoforms                <- myrange[rep(1,length(pAs)),]
                    mcols(new.isoforms)$isoform <- c(1:length(pAs))
                    names(new.isoforms)         <- paste(myID,c(1:length(pAs)),sep=".")
                    if(strd=="+"){end(new.isoforms)  <- start(new.isoforms)+pAs-1}
                    if(strd=="-"){start(new.isoforms)<- end(new.isoforms)-pAs-1}

                    if(strd=="+"){
                        #Correct pAs with considering the removed part due to introns
                        for(i in c(1:length(ri))){
                            tocor <- pAs>=ri[i]
                            if(sum(tocor)>0){pAs[tocor]=pAs[tocor]+wR[i]}
                    }
                    #Create real range
                    new.isoforms                <- myrange[rep(1,length(pAs)),]
                    mcols(new.isoforms)$isoform <- c(1:length(pAs))
                    names(new.isoforms)         <- paste(myID,c(1:length(pAs)),sep=".")
                    end(new.isoforms)           <- start(new.isoforms)+pAs-1
                    }
                    if(strd=="-"){
                        #Correct pAs with considering the removed part due to introns
                        pAs <- (length(vec)-pAs+1)[c(length(pAs):1)]
                        for(i in c(1:length(ri))){
                            tocor <- pAs>=ri[i]
                            if(sum(tocor)>0){pAs[tocor]=pAs[tocor]+wR[i]}
                        }
                        #Create real range
                        new.isoforms                <- myrange[rep(1,length(pAs)),]
                        mcols(new.isoforms)$isoform <- c(1:length(pAs))
                        names(new.isoforms)         <- paste(myID,c(1:length(pAs)),sep=".")
                        start(new.isoforms)         <- start(new.isoforms)+pAs-1
                    }
                }
                if(length(which(mym>MIN))==0){
                        new.isoforms   <- NA
                        pAs            <- 0
                }
            }
        }
        if(Km<2){
            new.isoforms   <- myrange
            pAs            <- 0
        }
        if(sum(mymed>10)<Llim){
            new.isoforms   <- NA
            pAs            <- 0
        }
        return(list(new.isoforms=new.isoforms,pAs=pAs))
    }
    if(length(vec)<=10){
            new.isoforms   <- myrange
            pAs            <- length(myvalues)
            return(list(new.isoforms=new.isoforms,pAs=pAs))
    }

}

### GET PAsites with EBS -- much more time
getpAEBS_old <- function(myobj=myres[[1]],Ix=3,wd=75,Llim=50,Km=5){
    #Select data of interest
    myvalues <- myobj$cov[[Ix]]
    myID     <- myobj$txID[Ix]
    myrange  <- myUTR[which(as.character(mcols(myUTR)$X.transcript_id.)==myID),]
    if(length(myrange)>1){
        myrange <- myrange[which(width(myrange)==length(myvalues))[1],]
    }
    strd     <- as.character(strand(myrange))
    vec      <- myvalues
    L        <- length(vec)

    #Check which of these extensions contain introns
    gOver    <- findOverlaps(query=myrange,subject=g3introns)
    idxU     <- queryHits(gOver)
    idxG     <- subjectHits(gOver)
    test     <- length(gOver)>0

    if(!test&length(vec)>10){
        if(strd=="-"){vec <- vec[c(length(vec):1)]}
        #Compute running median
        mymed  <- as.vector(runmed(x=vec,k=min(wd,length(vec)-10),endrule="median"))
        if(sum(mymed>10)>Llim){
            #Segment
            out   <- EBSegmentation(mymed,Kmax=5)
            mybic <- EBSBIC(out)
            bic   <- (EBSBIC(out)$NbBIC)-1
            GetHighestDensity <- function(ki){
                yi<-EBSDistrib(out,k=ki,Kk=bic+1)
                return(which(yi==max(yi)))
            }
            mypos  <- unlist(lapply(c(1:bic),GetHighestDensity))
            brks   <- c(1,mypos,length(mymed))#the first and the last must always stay
            mydist  <- diff(brks)
            while(sum(mydist<Llim)>0){
                    torm   <- which(mydist<Llim)[1]
                    if(torm==1){torm<-2}
                    brks   <- brks[-torm]
                    mydist <- diff(brks)
            }
            brks   <- unique(c(brks,length(mymed)))
            if(length(brks)==2){
                #print(paste("length brks=",length(brks),"easy!"))
                new.isoforms   <- myrange
                pAs            <- length(myvalues)
            }
            if(length(brks)>2){
                #print(paste("length brks=",length(brks),"more difficult!"))
                mysegs  <- list()
                for(i in c(2:length(brks))){
                 mysegs[[i]] <- c(brks[(i-1)]:(brks[i]))
                }
                mysegs <- mysegs[-1]
                #Compute median of each segment
                mym    <-  unlist(lapply(mysegs,function(x){
                                                            tempval <- mymed[x]
                                                            tempval <- tempval[tempval!=0]
                                                            return(mean(tempval))
                                                            }
                                                            )
                                  )
                #Identify pAs
                sel <- union(which(diff(mym)<0),max(which(mym!=0)))
                pAs <- sapply(sel,function(x)return(max(mysegs[[x]])))
                #Create real range
                new.isoforms                <- myrange[rep(1,length(pAs)),]
                mcols(new.isoforms)$isoform <- c(1:length(pAs))
                names(new.isoforms)         <- paste(myID,c(1:length(pAs)),sep=".")

                if(strd=="+"){end(new.isoforms)  <- start(new.isoforms)+pAs-1}
                if(strd=="-"){start(new.isoforms)<- end(new.isoforms)-pAs-1}

            }
        }
        if(sum(mymed>10)<=Llim){
            new.isoforms   <- myrange
            pAs            <- length(myvalues)
        }
            return(list(new.isoforms=new.isoforms,pAs=pAs))
    }

    #CASE with INTRONS
    if(test&length(vec)>10){
    ##                                                                          ##
    #Then should not account for the intron                                      #
    #                                                                            #
        startR <- start(g3introns[idxG,]) - start(myrange)+1
        endR   <- end(g3introns[idxG,]) - start(myrange)+1
        wR     <- endR-startR+1
        ri     <- rep(0,length(endR))
        for(i in c(1:length(startR))){
            if(i>1){
                ri[i] <- startR[i]-sum(wR[c(1:(i-1))])
                vec   <- vec[-seq(from=ri[i],to=(endR[i]-sum(wR[c(1:(i-1))])))]
            }
            if(i==1){
                ri[i] <- startR[i]
                vec   <- vec[-seq(from=ri[i],to=endR[i])]
            }
        }
        #if(strd=="-"){
        #k=1
        #  for(i in c(length(startR):1)){
        #    if(i>1){
        #        ri[k] <- end(myrange)-end(g3introns[idxG[i],])-sum(wR[c(length(startR):(k-1))]) + 1
        #        k=k+1
        #    }
        #    if(i==1){
        #        ri[k] <- end(myrange)-end(g3introns[idxG[i],]) + 1
        #        k=k+1
        #    }
        #  }
        #}
    #                                                                            #
    ##                                                                          ##
        if(strd=="-"){vec <- vec[c(length(vec):1)]}
        #Compute running median
        mymed  <- as.vector(runmed(x=vec,k=min(wd,length(vec)-10),endrule="median"))
        if(sum(mymed>10)>Llim){
            #Segment
            out   <- EBSegmentation(mymed,Kmax=5)
            mybic <- EBSBIC(out)
            bic   <- (EBSBIC(out)$NbBIC)-1
            GetHighestDensity <- function(ki){
                yi<-EBSDistrib(out,k=ki,Kk=bic+1)
                return(which(yi==max(yi)))
            }
            mypos  <- unlist(lapply(c(1:bic),GetHighestDensity))
            brks   <- c(1,mypos,length(mymed))#the first and the last must always stay
            mydist  <- diff(brks)
            while(sum(mydist<Llim)>0){
                    torm   <- which(mydist<Llim)[1]
                    if(torm==1){torm<-2}
                    brks   <- brks[-torm]
                    mydist <- diff(brks)
            }
            brks   <- unique(c(brks,length(mymed)))
            if(length(brks)==2){
                new.isoforms   <- myrange
                pAs            <- length(myvalues)
            }
            if(length(brks)>2){
                mysegs  <- list()
                for(i in c(2:length(brks))){
                    mysegs[[i]] <- c(brks[(i-1)]:(brks[i]))
                }
                mysegs <- mysegs[-1]
                #Compute median of each segment
                mym    <-  unlist(lapply(mysegs,function(x){
                                                            tempval <- mymed[x]
                                                            tempval <- tempval[tempval!=0]
                                                            return(mean(tempval))
                                                            }
                                                            )
                                  )
                #Identify pAs
                sel <- union(which(diff(mym)<0),max(which(mym!=0)))
                pAs <- sapply(sel,function(x)return(max(mysegs[[x]])))

                if(strd=="+"){
                    #Correct pAs with considering the removed part due to introns
                    for(i in c(1:length(ri))){
                        tocor <- pAs>=ri[i]
                        if(sum(tocor)>0){pAs[tocor]=pAs[tocor]+wR[i]}
                    }
                    #Create real range
                    new.isoforms                <- myrange[rep(1,length(pAs)),]
                    mcols(new.isoforms)$isoform <- c(1:length(pAs))
                    names(new.isoforms)         <- paste(myID,c(1:length(pAs)),sep=".")
                    end(new.isoforms)           <- start(new.isoforms)+pAs-1
                }
                if(strd=="-"){
                    #Correct pAs with considering the removed part due to introns
                    pAs <- (length(vec)-pAs+1)[c(length(pAs):1)]
                    for(i in c(1:length(ri))){
                        tocor <- pAs>=ri[i]
                        if(sum(tocor)>0){pAs[tocor]=pAs[tocor]+wR[i]}
                    }
                    #Create real range
                    new.isoforms                <- myrange[rep(1,length(pAs)),]
                    mcols(new.isoforms)$isoform <- c(1:length(pAs))
                    names(new.isoforms)         <- paste(myID,c(1:length(pAs)),sep=".")
                    start(new.isoforms)         <- start(new.isoforms)+pAs-1
                }
            }
        }
        if(sum(mymed>10)<Llim){
            new.isoforms   <- myrange
            pAs            <- length(myvalues)
        }
            return(list(new.isoforms=new.isoforms,pAs=pAs))
    }
    if(length(vec)<=10){
            new.isoforms   <- myrange
            pAs            <- length(myvalues)
            return(list(new.isoforms=new.isoforms,pAs=pAs))
    }


}



getpAs <- function(myobj=myres[[1]],Ix=3,wd=75,Llim=50){

    #Select data of interest
    myvalues <- myobj$cov[[Ix]]
    myID     <- myobj$txID[Ix]
    myrange  <- myUTR[which(as.character(mcols(myUTR)$X.transcript_id.)==myID),]
    if(length(myrange)>1){
        myrange <- myrange[which(width(myrange)==length(myvalues))[1],]
    }
    strd     <- as.character(strand(myrange))
    vec      <- myvalues
    L        <- length(vec)

    #Check which of these extensions contain introns
    gOver    <- findOverlaps(query=myrange,subject=g3introns)
    idxU     <- queryHits(gOver)
    idxG     <- subjectHits(gOver)
    test     <- length(gOver)>0

    if(!test){
        if(strd=="-"){vec <- vec[c(length(vec):1)]}

        #Compute running median
        mymed  <- as.vector(runmed(x=vec,k=wd,endrule="median"))
        if(sum(mymed>10)>Llim){
            #Segment
            out   <- EBSegmentation(mymed,Kmax=5)
            mybic <- EBSBIC(out)
            bic   <- (EBSBIC(out)$NbBIC)-1
            GetHighestDensity <- function(ki){
                yi<-EBSDistrib(out,k=ki,Kk=bic+1)
                return(which(yi==max(yi)))
            }
            mypos  <- unlist(lapply(c(1:bic),GetHighestDensity))
            brks   <- c(1,mypos,length(mymed))#the first and the last must always stay
            mydist <- diff(brks)
            while(sum(mydist<Llim)>0){
                    brks   <- brks[-which(mydist<Llim)[1]]
                    mydist <- diff(brks)
            }
            brks   <- unique(c(brks,length(mymed)))
            if(length(brks)==2){
                new.isoforms   <- myrange
                pAs            <- length(myvalues)
            }
            if(length(brks>2)){
                mysegs  <- list()
                for(i in c(2:length(brks))){
                 mysegs[[i]] <- c(brks[(i-1)]:(brks[i]))
                }
                mysegs <- mysegs[-1]
                #Compute median of each segment
                mym    <-  unlist(lapply(mysegs,function(x)return(median(mymed[x]))))
                #Identify pAs
                sel <- union(which(diff(mym)<0),max(which(mym!=0)))
                pAs <- sapply(sel,function(x)return(max(mysegs[[x]])))
                #Create real range
                new.isoforms                <- myrange[rep(1,length(pAs)),]
                mcols(new.isoforms)$isoform <- c(1:length(pAs))
                names(new.isoforms)         <- paste(myID,c(1:length(pAs)),sep=".")

                if(strd=="+"){end(new.isoforms)  <- start(new.isoforms)+pAs-1}
                if(strd=="-"){start(new.isoforms)<- end(new.isoforms)-pAs-1}

            }
        }
        if(sum(mymed>10)<Llim){
            new.isoforms   <- myrange
            pAs            <- length(myvalues)
            mybic          <- NA
        }
    }

    #CASE with INTRONS
    if(test){
    ##                                                                          ##
    #Then should not account for the intron                                      #
    #                                                                            #
        startR <- start(g3introns[idxG,]) - start(myrange)+1
        endR   <- end(g3introns[idxG,]) - start(myrange)+1
        wR     <- endR-startR+1
        ri     <- rep(0,length(endR))
        for(i in c(1:length(startR))){
            if(i>1){
                ri[i] <- startR[i]-sum(wR[c(1:(i-1))])
                vec   <- vec[-seq(from=ri[i],to=(endR[i]-sum(wR[c(1:(i-1))])))]
            }
            if(i==1){
                ri[i] <- startR[i]
                vec   <- vec[-seq(from=ri[i],to=endR[i])]
            }
        }
        #if(strd=="-"){
        #k=1
        #  for(i in c(length(startR):1)){
        #    if(i>1){
        #        ri[k] <- end(myrange)-end(g3introns[idxG[i],])-sum(wR[c(length(startR):(k-1))]) + 1
        #        k=k+1
        #    }
        #    if(i==1){
        #        ri[k] <- end(myrange)-end(g3introns[idxG[i],]) + 1
        #        k=k+1
        #    }
        #  }
        #}
    #                                                                            #
    ##                                                                          ##
        if(strd=="-"){vec <- vec[c(length(vec):1)]}

        #Compute running median
        mymed  <- as.vector(runmed(x=vec,k=wd,endrule="median"))
        if(sum(mymed>10)>Llim){
            #Segment
            out     <- EBSegmentation(mymed,Kmax=5)
            mybic   <- EBSBIC(out)
            bic     <- (EBSBIC(out)$NbBIC)-1
            GetHighestDensity <- function(ki){
                yi<-EBSDistrib(out,k=ki,Kk=bic+1)
                return(which(yi==max(yi)))
            }
            mypos   <- unlist(lapply(c(1:bic),GetHighestDensity))
            brks   <- c(1,mypos,length(mymed))#the first and the last must always stay
            mydist <- diff(brks)
            while(sum(mydist<Llim)>0){
                    brks   <- brks[-which(mydist<Llim)[1]]
                    mydist <- diff(brks)
            }
            brks   <- unique(c(brks,length(mymed)))
            if(length(brks)==2){
                new.isoforms   <- myrange
                pAs            <- length(myvalues)
            }
            if(length(brks>2)){
                mysegs  <- list()
                for(i in c(2:length(brks))){
                    mysegs[[i]] <- c(brks[(i-1)]:(brks[i]))
                }
                mysegs <- mysegs[-1]
                #Compute median of each segment
                mym    <-  unlist(lapply(mysegs,function(x)return(median(mymed[x]))))
                #Identify pAs
                sel <- union(which(diff(mym)<0),max(which(mym!=0)))
                pAs <- sapply(sel,function(x)return(max(mysegs[[x]])))

                if(strd=="+"){
                    #Correct pAs with considering the removed part due to introns
                    for(i in c(1:length(ri))){
                        tocor <- pAs>=ri[i]
                        if(sum(tocor)>0){pAs[tocor]=pAs[tocor]+wR[i]}
                    }
                    #Create real range
                    new.isoforms                <- myrange[rep(1,length(pAs)),]
                    mcols(new.isoforms)$isoform <- c(1:length(pAs))
                    names(new.isoforms)         <- paste(myID,c(1:length(pAs)),sep=".")
                    end(new.isoforms)           <- start(new.isoforms)+pAs-1
                }
                if(strd=="-"){
                    #Correct pAs with considering the removed part due to introns
                    pAs <- (length(vec)-pAs+1)[c(length(pAs):1)]
                    for(i in c(1:length(ri))){
                        tocor <- pAs>=ri[i]
                        if(sum(tocor)>0){pAs[tocor]=pAs[tocor]+wR[i]}
                    }
                    #Create real range
                    new.isoforms                <- myrange[rep(1,length(pAs)),]
                    mcols(new.isoforms)$isoform <- c(1:length(pAs))
                    names(new.isoforms)         <- paste(myID,c(1:length(pAs)),sep=".")
                    start(new.isoforms)         <- start(new.isoforms)+pAs-1
                }
            }
        }
        if(sum(mymed>10)<Llim){
            new.isoforms   <- myrange
            pAs            <- length(myvalues)
            mybic          <- NA
        }
    }

    return(list(new.isoforms=new.isoforms,pAs=pAs,mybic=mybic))

}




getpAsOld <- function(myobj=mytX[[is]][[2]],Ix,myID,wd=151,Llim=50,Plot=TRUE,nsample=NA){

    #Select data of interest
    myvalues <- myobj[[1]][[Ix]]
    mylims   <- myobj[[2]][[Ix]]
    myID     <- names(myobj[[1]])[Ix]
    myrange  <- mergedGRS[mergedGRS$ID==myID]
    strd     <- as.character(strand(mergedGRS[mergedGRS$ID==myID]))

    #Check which of these extensions contain introns
    gOver    <- findOverlaps(query=myrange,subject=g3introns)
    idxU     <- queryHits(gOver)
    idxG     <- subjectHits(gOver)
    test     <- length(as.matrix(gOver))>2


    #Focus on 3'UTR
    Lbefore<-  mylims[nrow(mylims),1]-1
    vec    <- as.vector(myvalues)[c(mylims[nrow(mylims),1]:mylims[nrow(mylims),2])]
    L      <- length(vec)


    if(!test){


        #Compute running median
        mymed  <- runmed(x=vec,k=wd,endrule="median")
        #Segment
        Seg     <-Segmentor(as.vector(mymed),model=3,Kmax=10)
        Kchoose <-SelectModel(Seg, penalty="oracle")
        brks    <- as.vector(getBreaks(Seg)[Kchoose, 1:Kchoose])
        brks    <- c(1,brks[-which(diff(brks)<Llim)-1])
        mysegs  <- list()
        for(i in c(2:length(brks))){
            mysegs[[i]] <- c(brks[(i-1)]:brks[i])
        }
        mysegs <- mysegs[-1]
        #Compute median of each segment
        mym    <-  unlist(lapply(mysegs,function(x)return(median(mymed[x]))))
        #Identify pAs
        sel <- union(which(diff(mym)<0),max(which(mym!=0)))
        pAs <- sapply(sel,function(x)return(max(mysegs[[x]])))
        #Compute width of first stable peak
        myw <- length(mysegs[[max(which(mym!=0))]])
        pAs <- pAs+Lbefore

        #Plot
        if(Plot){
            header <- myID
            txCov  <- myvalues
            txLim  <- mylims
            require(grDevices)

            H    <- 0.025*(max(c(1,max(txCov)-min(txCov))))
            dh   <- 0.01*H
            lLm  <- -0.7*H
            plot(txCov,type="l",col="black",lwd=2,xlab="",ylab="",main="",frame=FALSE)
            mypalette             <- colorRampPalette(colors=c("white","black"),bias = 1, space = c("rgb"),interpolate = c("linear"))
            mycols                <- mypalette(n=nrow(txLim))
            mycols[1]             <- "red"
            mycols[length(mycols)]<- "green"
            for(k in c(1:nrow(txLim))){
                rect(txLim[k,1],lLm,txLim[k,2],lLm+H,border=NA,col=mycols[k])
            }
            lines(txCov,col="black",lwd=2)
            lines(runmed(x=txCov,k=wd,endrule="median"),col="red",lwd=2)
            mtext(side=1,line=2,text="position [nt]")
            mtext(side=2,line=2,text="coverage")
            mtext(side=3,line=0.5,text=header)
            abline(v=pAs,col="green",lwd=2,lty=2)
            if(!is.na(nsample)){mtext(side=3,line=1.5,text=nsample)}
        }


        return(list(myw,pAs))

    }


    if(test){
    ##                                                                          ##
    #Then should not account for the intron                                      #
    #                                                                            #
        startR <- start(g3introns[idxG[which(idxU==ix)],]) - start(totRegion[ix,])
        endR   <- end(g3introns[idxG[which(idxU==ix)],]) - start(totRegion[ix,])
        intoRM <- seq(from=startR,to=endR)
        Ls     <- length(intoRM)
        vec    <- vec[-intoRM]
        L      <- length(vec)
    #                                                                            #
    ##                                                                          ##
        tempv <- c(rep(median(vec[c(1:floor(wd/2))]),floor(wd/2)),vec,rep(median(vec[c((L-floor(wd/2)+1):L)]),floor(wd/2)))
        Lp    <- length(tempv)
        idx1  <-  c(1:(Lp-2*floor(wd/2)))
        idx2  <-  c((1+2*floor(wd/2)):Lp)
        if(strd=="+"){
            temp <- (tempv[idx2]-tempv[idx1])/wd
            Ll   <- length(temp)
            #Re-add last values to introns for plotting purpose
            lastV <-  temp[(startR-1)]
            temp  <- c(temp[c(1:(startR-1))],rep(lastV,Ls),temp[c(startR:Ll)])
        }
        else{
            tempv <- rev(tempv)
            temp <- (tempv[idx2]-tempv[idx1])/wd
            temp <- rev(temp)
            Ll   <- length(temp)
            lastV <-  temp[startR]
            temp <- c(temp[c(1:(startR-1))],rep(lastV,Ls),temp[c(startR:Ll)])
         }
    }

    return(temp)
}



