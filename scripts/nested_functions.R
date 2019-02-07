### NESTED FUNCTIONS
StudyFrequencyPAS <- function(utrGR=focusGR,utrAnno=focusAnno,target=myGRpaf){
  # a) Find overlap
  gOver               <- findOverlaps(query=target,subject=utrGR,ignore.strand=F)
  idxU                <- queryHits(gOver)
  idxG                <- subjectHits(gOver)
  # b) Compute frequency
  Lev                    <- c("-1","0","1","2","3","4")
  subAnno                <- utrAnno[unique(idxG),]
  all                    <- unlist(lapply(Lev,function(x)return(sum(as.character(utrAnno$iso_cons_merged)[utrAnno$detect.iso.ngf!="no"]==x,na.rm=T))))
  all                    <- c(sum(utrAnno$detect.iso.ngf[utrAnno$is.conservative]!="no"),all)
  names(all)             <- c("conservative",Lev)
  
  indb                   <- unlist(lapply(Lev,function(x)return(sum(as.character(subAnno$iso_cons_merged)[subAnno$detect.iso.ngf!="no"]==x,na.rm=T))))
  indb                   <- c(sum(subAnno$detect.txID.ngf[subAnno$is.conservative]!="no"),indb)
  names(indb)            <- c("conservative",Lev)
  
  # c) Create output
  out                     <- data.frame(all=all,indb=indb,freq=indb/all)
  return(out)
}


StudyFrequencyPAS_v2 <- function(utrGR=focusGR,utrAnno=focusAnno,target=myGRpaf,version="all"){
  
  if(version=="all"){
    print("no modif required")
  }
  if(version=="conservative"){
    utrGR   <- utrGR[utrAnno$is.conservative,]
    utrAnno <- utrAnno[utrAnno$is.conservative,]
  }
  if(version=="new"){
    utrGR   <- utrGR[!utrAnno$is.conservative,]
    utrAnno <- utrAnno[!utrAnno$is.conservative,]
  }
  
  # a) Find overlap
  gOver               <- findOverlaps(query=target,subject=utrGR,ignore.strand=F)
  idxU                <- queryHits(gOver)
  idxG                <- subjectHits(gOver)
  subAnno             <- utrAnno[unique(idxG),]
  
  #b) Compute frequency per bin of length
  #BINS2               <- cut(log10(utrAnno$newL),breaks=quantile(log10(utrAnno$newL),prob=seq(from=0,to=1.0,by=0.1),na.rm=T),include.lowest=T)
  BINS1                <- cut(log10(utrAnno$newL),breaks=c(1,1.5,2,2.25,2.75,3.25,3.75,4.25),include.lowest=T)
  BINS2               <- cut(log10(subAnno$newL),breaks=c(1,1.5,2,2.25,2.75,3.25,3.75,4.25),include.lowest=T)
  Lev                  <- levels(BINS)
  
  all                 <- as.vector(table(BINS1))
  all                 <- c(sum(utrAnno$is.conservative),all)
  names(all)          <- c("conservative",Lev)
  
  indb                <- as.vector(table(BINS2))
  indb                <- c(sum(subAnno$is.conservative),indb)
  names(indb)         <- c("conservative",Lev)
  
  # c) Create output
  out                 <- data.frame(all=all,indb=indb,freq=indb/all)
  return(out)
}



StudyNumberPAS <- function(utrGR=focusGR,utrAnno=focusAnno,target=myPA[[1]]){
  
  Lev                    <- c("-1","0","1","2","3","4")
  sel1                   <- utrAnno$detect.iso.ngf!="no"
  mysel                  <- lapply(Lev,function(x){
    temp<-as.character(utrAnno$iso_cons_merged)==x
    temp[is.na(temp)]<-FALSE
    return(temp)})
  myListUTR              <- unlist(lapply(c(1:length(mysel)),function(x)return(utrGR[sel1&mysel[[x]],])))
  
  # Compute number of sites per category
  GetNumberSites <- function(subUTR=myListUTR[[1]]){
    myHits      <- subjectHits(findOverlaps(query=target,subject=subUTR,ignore.strand=F))
    no.missing  <- sum(!c(1:length(subUTR))%in%unique(myHits))
    out         <- as.data.frame(table(table(myHits)))#get number of sites per isoform
    temp        <- sum(as.numeric(as.character(out$Var1))*as.numeric(as.character(out$Freq)))/(no.missing+sum(as.numeric(as.character(out$Freq))))
    Out         <- c(rep(0,no.missing),do.call(lapply(c(1:nrow(out)),function(x)return(rep(out$Var1[x],out$Freq[x]))),what=c))
    
    return(list(temp,Out))
  }
  
  no.sites    <-  lapply(myListUTR,GetNumberSites)
  
  
  my.all.sites <- lapply(c(1:length(no.sites)),function(x)return(no.sites[[x]][[2]]))
  my.all.sites[[length(my.all.sites)+1]] <- GetNumberSites(utrGR[sel1&utrAnno$is.conservative,])[[2]]
  my.all.sites <- my.all.sites[c(7,c(1:6))]
  names(my.all.sites)<- c("conservative",Lev)
  
  my.no.with.sites <- rbind(unlist(lapply(my.all.sites,function(x)return(sum(x!=0)))),
                            unlist(lapply(my.all.sites,function(x)return(length(x))))
  )
  
  rownames(my.no.with.sites)<- c("with sites", "tot")
  
  out          <- unlist(lapply(c(1:length(no.sites)),function(x)return(no.sites[[x]][[1]])))
  out          <- c(GetNumberSites(utrGR[sel1&utrAnno$is.conservative,])[[1]],out)
  names(out)   <-c("conservative",Lev)
  
  return(list(out,my.all.sites,my.no.with.sites))
}



StudyNumberPASPerLength <- function(utrGR=focusGR,utrAnno=focusAnno,target=myPA[[1]],Lev=BINS1){
  
  sel1                   <- utrAnno$detect.iso.ngf!="no"
  mysel                  <- lapply(Lev,function(x){
    temp             <-as.character(utrAnno$iso_cons_merged)==x
    temp[is.na(temp)]<-FALSE
    return(temp)
  }
  )
  
  myListUTR              <- unlist(lapply(c(1:length(mysel)),function(x)return(utrGR[sel1&mysel[[x]],])))
  
  # Compute number of sites per category
  GetNumberSites <- function(subUTR=myListUTR[[1]]){
    myHits      <- subjectHits(findOverlaps(query=target,subject=subUTR,ignore.strand=F))
    no.missing  <- sum(!c(1:length(subUTR))%in%unique(myHits))
    out         <- as.data.frame(table(table(myHits)))#get number of sites per isoform
    
    temp       <- sum(as.numeric(as.character(out$Var1))*as.numeric(as.character(out$Freq)))/(no.missing+sum(as.numeric(as.character(out$Freq))))
    
    Out <- c(rep(0,no.missing),do.call(lapply(c(1:nrow(out)),function(x)return(rep(out$Var1[x],out$Freq[x]))),what=c))
    
    return(list(temp,Out))
  }
  
  no.sites    <-  lapply(myListUTR,GetNumberSites)
  
  
  my.all.sites <- lapply(c(1:length(no.sites)),function(x)return(no.sites[[x]][[2]]))
  my.all.sites[[length(my.all.sites)+1]] <- GetNumberSites(utrGR[sel1&utrAnno$is.conservative,])[[2]]
  my.all.sites <- my.all.sites[c(7,c(1:6))]
  names(my.all.sites)<- c("conservative",Lev)
  
  my.no.with.sites <- rbind(unlist(lapply(my.all.sites,function(x)return(sum(x!=0)))),
                            unlist(lapply(my.all.sites,function(x)return(length(x))))
  )
  
  rownames(my.no.with.sites)<- c("with sites", "tot")
  
  out          <- unlist(lapply(c(1:length(no.sites)),function(x)return(no.sites[[x]][[1]])))
  out          <- c(GetNumberSites(utrGR[sel1&utrAnno$is.conservative,])[[1]],out)
  names(out)   <-c("conservative",Lev)
  
  return(list(out,my.all.sites,my.no.with.sites))
}

StudyNumberPAS_v2 <- function(utrGR=focusGR,utrAnno=focusAnno,target=myPA[[1]],version="all"){
  
  if(version=="all"){
    print("no modif required")
  }
  if(version=="conservative"){
    utrGR   <- utrGR[utrAnno$is.conservative,]
    utrAnno <- utrAnno[utrAnno$is.conservative,]
  }
  if(version=="new"){
    utrGR   <- utrGR[!utrAnno$is.conservative,]
    utrAnno <- utrAnno[!utrAnno$is.conservative,]
  }
  
  BINS1                <- cut(log10(utrAnno$newL),breaks=c(1,1.5,2,2.25,2.75,3.25,3.75,4.25),include.lowest=T)
  print("debug1")
  Lev                  <- levels(BINS1)
  print("debug2")
  myListUTR            <- tapply(c(1:length(utrGR)),INDEX=BINS1,FUN=function(x)return(utrGR[x]))
  print("debug3")
  print(length(myListUTR))
  
  # Compute number of sites per category
  GetNumberSites <- function(subUTR=myListUTR[[1]]){
    if(length(subUTR)>10){
      myHits      <- subjectHits(findOverlaps(query=target,subject=subUTR,ignore.strand=F))
      print("debug4")
      no.missing  <- sum(!c(1:length(subUTR))%in%unique(myHits))
      print("debug5")
      if(length(myHits)==0){
        out         <- rep(0,length(subUTR))
        Out         <- rep(0,no.missing)
      }
      
      
      
      if(length(myHits>0)){
        out         <- as.data.frame(table(table(myHits)))#get number of sites per isoform
        print("debug6")
        temp       <- mean(c(as.data.frame(table(myHits))$Freq,rep(0,no.missing)))#Average number of PAS per isofor
        print("debug7")
        Out        <- c(rep(0,no.missing),do.call(lapply(c(1:nrow(out)),function(x)return(rep(out$Var1[x],out$Freq[x]))),what=c))#Total number of motifs
      }
    }
    else{
      temp=NA
      Out=NA
    }
    return(list(temp,Out))
  }
  
  no.sites                                <-  lapply(myListUTR,GetNumberSites)
  
  
  my.all.sites                            <- lapply(c(1:length(no.sites)),function(x)return(no.sites[[x]][[2]]))
  my.all.sites[[length(my.all.sites)+1]]  <- GetNumberSites(utrGR[utrAnno$is.conservative,])[[2]]
  my.all.sites                            <- my.all.sites[c((length(Lev)+1),c(1:length(Lev)))]
  names(my.all.sites)                     <- c("conservative",Lev)
  
  my.no.with.sites                       <- rbind(unlist(lapply(my.all.sites,function(x)return(sum(x!=0)))),
                                                  unlist(lapply(my.all.sites,function(x)return(length(x))))
  )
  
  rownames(my.no.with.sites)<- c("with sites", "tot")
  
  out          <- unlist(lapply(c(1:length(no.sites)),function(x)return(no.sites[[x]][[1]])))
  out          <- c(GetNumberSites(utrGR[utrAnno$is.conservative,])[[1]],out)#Average number of PAS per isoform
  names(out)   <-c("conservative",Lev)
  
  return(list(out,my.all.sites,my.no.with.sites))
}


#This one needs corrections as performed with the motif analysis
GetZScore <- function(){
  
  # A. FG
  #1/ Sum over all uniqueID
  mysumid <- tapply(myfc[,3],INDEX=as.factor(as.character(myfc[,5])),FUN=sum)
  mysumid <- data.frame(ID=rownames(mysumid),iso=anno_ngf$iso_cons[match(rownames(mysumid),as.character(anno_ngf$uniqueID))],val=as.vector(mysumid))
  #2/ Group by iso_cons
  if(withsum){fgcg      <- tapply(mysumid$val,INDEX=as.factor(mysumid$iso),FUN=sum)}
  if(!withsum){fgcg     <- tapply(mysumid$val,INDEX=as.factor(mysumid$iso),FUN=mean)}
  scorec                <- fgcg[match(c("-1","0","1","2","3"),names(fgcg))]
  
  # B. BG
  #shuffle labels 500 times and compute global score
  myscoreg  <- list()
  myscorec  <- list()
  
  # Create matrix of randomised index/labels: when I shuffle the iso_cons, I need to shuffle the uniqueID together otherwise the same txID will be find in several places; donc randomisation en pair.
  #1/ Reshuffle the labels
  myIdx <- matrix(NA,nrow=nrow(mysumid),ncol=500)
  for(k in c(1:500)){
    myIdx[,k]<- sample(as.character(mysumid$iso))
  }
  #2/ Sum over all the labels
  for(j in c(1:500)){
    if(withsum){myscoreg[[j]]   <- tapply(mysumid$val,INDEX=as.factor(myIdx[,j]),FUN=sum)}
    if(!withsum){myscoreg[[j]]  <- tapply(mysumid$val,INDEX=as.factor(myIdx[,j]),FUN=mean)}
    myscoreg[[j]]               <- myscoreg[[j]][match(c("-1","0","1","2","3"),names(myscoreg[[j]]))]
  }
  
  
  myscorec  <-  do.call(myscoreg,what=rbind)
  
  mymbg      <- apply(myscorec,2,mean)
  msdbg      <- apply(myscorec,2,sd)
  zval       <- (scorec-mymbg)/msdbg
  
  pval.more  <- vector(length=length(zval))
  pval.less  <- vector(length=length(zval))
  for(i in c(1:ncol(myscorec))){
    pval.more[i]<- sum(myscorec[,i]>scorec[i])/length(myscorec[,i])
    pval.less[i]<- sum(myscorec[,i]<scorec[i])/length(myscorec[,i])
  }
  
  log.p.more <- -log10(pval.more+0.001)
  log.p.less <- -log10(pval.less+0.001)
  log.p.tot  <- log.p.more
  log.p.tot[abs(log.p.more)<abs(log.p.less)]<- -log.p.less[abs(log.p.more)<abs(log.p.less)]
  
}

GetZScore_v2 <- function(){
  # A. FG
  #1/ Sum over all uniqueID
  mysumid <- tapply(myfc[,3],INDEX=as.factor(as.character(myfc[,5])),FUN=sum)
  mysumid <- data.frame(ID=rownames(mysumid),iso=anno_ngf$iso_cons[match(rownames(mysumid),as.character(anno_ngf$uniqueID))],val=as.vector(mysumid))
  #2/ Group by iso_cons
  if(withsum){fgcg      <- tapply(mysumid$val,INDEX=as.factor(mysumid$iso),FUN=sum)}
  if(!withsum){fgcg     <- tapply(mysumid$val,INDEX=as.factor(mysumid$iso),FUN=mean)}
  scorec                <- fgcg[match(c("-1","0","1","2","3"),names(fgcg))]
  
  # B. BG
  #shuffle labels 500 times and compute global score
  myscoreg  <- list()
  myscorec  <- list()
  
  # Create matrix of randomised index/labels: when I shuffle the iso_cons, I need to shuffle the uniqueID together otherwise the same txID will be find in several places; donc randomisation en pair.
  #1/ Reshuffle the labels
  myIdx <- matrix(NA,nrow=nrow(mysumid),ncol=500)
  for(k in c(1:500)){
    myIdx[,k]<- sample(as.character(mysumid$iso))
  }
  #2/ Sum over all the labels
  for(j in c(1:500)){
    if(withsum){myscoreg[[j]]   <- tapply(mysumid$val,INDEX=as.factor(myIdx[,j]),FUN=sum)}
    if(!withsum){myscoreg[[j]]  <- tapply(mysumid$val,INDEX=as.factor(myIdx[,j]),FUN=mean)}
    myscoreg[[j]]               <- myscoreg[[j]][match(c("-1","0","1","2","3"),names(myscoreg[[j]]))]
  }
  
  
  myscorec  <-  do.call(myscoreg,what=rbind)
  
  mymbg      <- apply(myscorec,2,mean)
  msdbg      <- apply(myscorec,2,sd)
  zval       <- (scorec-mymbg)/msdbg
  
  pval.more  <- vector(length=length(zval))
  pval.less  <- vector(length=length(zval))
  for(i in c(1:ncol(myscorec))){
    pval.more[i]<- sum(myscorec[,i]>scorec[i])/length(myscorec[,i])
    pval.less[i]<- sum(myscorec[,i]<scorec[i])/length(myscorec[,i])
  }
  
  log.p.more <- -log10(pval.more+0.001)
  log.p.less <- -log10(pval.less+0.001)
  log.p.tot  <- log.p.more
  log.p.tot[abs(log.p.more)<abs(log.p.less)]<- -log.p.less[abs(log.p.more)<abs(log.p.less)]
  
}

MyTestFrac <- function(n1,n2,p1,p2){
  p  <- (p1*n1+p2*n2)/(n1+n2)
  SE <- sqrt(p*(1-p)*(1/n1+1/n2))
  Z  <- (p1-p2)/SE
  PVal<- 2*pnorm(-abs(Z))#2sided
  return(list(Z,PVal))
}


error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}



#THIS MUST BE MODIFIED
ComputeEnrichmentPAS <- function(ID=unique(as.character(mydat[,5])),withsum=TRUE,Nrand=500,PLOT=TRUE,BG=anno_tot$uniqueID[anno_tot$isoform=="0"]){
  
  ID <- ID[!is.na(ID)]
  ix <- match(ID,rownames(mysumpost))
  # A. FG
  #1/ Group by iso_cons
  if(withsum){fgcg       <- apply(mysumpost[ix,],2,FUN=sum)}
  if(!withsum){fgcg      <- apply(mysumpost[ix,],2,FUN=mean)}
  
  # B. BG
  #Pick labels 500 times and compute global score
  if(withsum){
    myscorec       <- do.call(what=rbind,
                              args=lapply(c(1:Nrand),FUN=function(Z)return(apply(mysumpost[match(sample(x=BG,size=length(ID)),rownames(mysumpost)),],2,FUN=sum))))
  }
  
  
  if(!withsum){
    myscorec       <- do.call(what=rbind,
                              args=lapply(c(1:Nrand),FUN=function(Z)return(apply(mysumpost[match(sample(x=BG,size=length(ID)),rownames(mysumpost)),],2,FUN=mean))))
  }
  
  scorec     <-  fgcg
  mymbg      <- apply(myscorec,2,mean)
  msdbg      <- apply(myscorec,2,sd)
  zval       <- (fgcg-mymbg)/msdbg
  
  pval.more  <- vector(length=length(zval))
  pval.less  <- vector(length=length(zval))
  for(i in c(1:ncol(myscorec))){
    pval.more[i]<- sum(myscorec[,i]>scorec[i])/length(myscorec[,i])
    pval.less[i]<- sum(myscorec[,i]<scorec[i])/length(myscorec[,i])
  }
  
  log.p.more <- -log10(pval.more+0.001)
  log.p.less <- -log10(pval.less+0.001)
  log.p.tot  <- log.p.more
  log.p.tot[abs(log.p.more)<abs(log.p.less)]<- -log.p.less[abs(log.p.more)<abs(log.p.less)]
  
}

GetOI <- function(mygoID="GO:0006412",sampleGO){
  go.genes    <- genesInTerm(sampleGO, mygoID)[[1]]#To extract all genes related to this term
  sig.genes   <- sigGenes(sampleGO)
  goi         <- intersect(sig.genes,go.genes)
  return(goi)
}

CompareBP <- function(enr1=test1,enr2=test2,PLOT=TRUE,no=10,lab1="remodelled",lab2="transport"){
  
  myterms            <- unique(c(enr1$Term,enr2$Term))
  temp1              <- -log10(enr1[match(myterms,enr1$Term),]$P.DE)
  temp2              <- -log10(enr2[match(myterms,enr2$Term),]$P.DE)
  temp1[is.na(temp1)]<- 0
  names(temp1)      <-  myterms
  temp2[is.na(temp2)]<- 0
  names(temp2)       <- myterms
  out                <- data.frame(term=myterms,val1=temp1,val2=temp2)
  
  
  if(PLOT){
    
    selterms <- unique(c(names(sort(temp1,decreasing=T))[c(1:no)],names(sort(temp2,decreasing=T))[c(1:no)]))
    dat      <- out[match(selterms,out$term),]
    dat      <- dat[sort(dat$val1-dat$val2,decreasing=T,index.return=T)$ix,]
    L        <- nrow(dat)
    val1     <- as.vector(-dat$val1)
    val2     <- as.vector(dat$val2)
    
    par(mar=c(3,10,2,3),cex=0.7)
    plot(c(L+3,0),xlim=c(min(-dat$val1)-10,max(dat$val2)),type = "n",frame=F,yaxt="n",ylab="")
    mp=barplot(height = val1,add = TRUE,axes = FALSE,horiz=T,col="lightsteelblue3")
    barplot(height = val2,add = TRUE,axes = FALSE,horiz=T,col="midnightblue")
    text(x=(min(-dat$val1)-5),y=mp,lab=as.character(dat$term),cex=0.6)
    mtext(side=1,line=2,text="-log10[P-Value]",cex=0.6)
    legend("top",pch=15,col=c("lightsteelblue3","midnightblue"),bty="n",ncol=2,leg=c("axonal remodelling","facilitated transport"),cex=0.5)
  }
  return(out)
}


CompareBP_improved <- function(enr1=enrichLong[[2]][[3]],enr2=enrichShort[[2]][[3]],PLOT=TRUE,no=10,lab1="long",lab2="short",coi="weight0Fisher"){
  
  myterms            <- unique(c(enr1$Term[enr1$Significant>=5&as.numeric(enr1[,match(coi,colnames(enr1))])<=0.05],enr2$Term[enr1$Significant>=5&as.numeric(enr2[,match(coi,colnames(enr2))])<=0.05]))
  temp1              <- -log10(as.numeric(enr1[match(myterms,enr1$Term),match(coi,colnames(enr1))]))
  temp2              <- -log10(as.numeric(enr2[match(myterms,enr2$Term),match(coi,colnames(enr2))]))
  temp1[is.na(temp1)]<- 0
  names(temp1)      <-  myterms
  temp2[is.na(temp2)]<- 0
  names(temp2)       <- myterms
  out                <- data.frame(term=myterms,val1=temp1,val2=temp2)
  
  if(PLOT){
    
    selterms <- unique(c(names(sort(temp1,decreasing=T))[c(1:no)],names(sort(temp2,decreasing=T))[c(1:no)]))
    dat      <- out[match(selterms,out$term),]
    dat      <- dat[sort(dat$val1-dat$val2,decreasing=T,index.return=T)$ix,]
    L        <- nrow(dat)
    val1     <- as.vector(-dat$val1)
    val2     <- as.vector(dat$val2)
    
    par(mar=c(3,10,2,3),cex=0.7)
    plot(c(L+3,0),xlim=c(min(-dat$val1)-10,max(dat$val2)),type = "n",frame=F,yaxt="n",ylab="")
    mp=barplot(height = val1,add = TRUE,axes = FALSE,horiz=T,col="lightsteelblue3")
    barplot(height = val2,add = TRUE,axes = FALSE,horiz=T,col="midnightblue")
    text(x=(min(-dat$val1)-5),y=mp,lab=as.character(dat$term),cex=0.6)
    mtext(side=1,line=2,text="-log10[P-Value]",cex=0.6)
    legend("top",pch=15,col=c("lightsteelblue3","midnightblue"),bty="n",ncol=2,leg=c(lab1,lab2),cex=0.5)
  }
  colnames(out)      <- c("terms",lab1,lab2)
  
  return(out)
}

#Enrichment with GOANA
GetEnrich <- function(selection=selRUD[[5]],value=dRUD){
  mydat           <- subanno[selection,]
  if(!is.na(value)){
    print("with value")
    mydat <- data.frame(mydat,t2g[match(mydat$geneSymbol,t2g$external_gene_name),],myval=value[selection])
    go    <- goana(unique(mydat$entrezgene), trend=mydat$myval, universe = t2g$entrezgene, species = "Rn", prior.prob = NULL, covariate=NULL,plot=FALSE,FDR=0.01)
    BP    <- topGO(go, ontology="BP",number=Inf)
    return(BP[BP$P.DE<=0.01,])
  }
  else{
    mydat <- data.frame(mydat,t2g[match(mydat$geneSymbol,t2g$external_gene_name),])
    go    <- goana(unique(mydat$entrezgene), trend=NULL, universe = t2g$entrezgene, species = "Rn", prior.prob = NULL, covariate=NULL,plot=FALSE,FDR=0.01)
    BP    <- topGO(go, ontology="BP",number=Inf)
    return(BP[BP$P.DE<=0.01,])
  }
}


CreateSampleGO <- function(mysel){
  sampleGO               <- list()
  
  geneNames              <- myBG$txID
  myInterestingGenes     <- unique(as.character(subanno$txID)[mysel])
  myInterestingGenes     <- myInterestingGenes[!is.na(myInterestingGenes)]
  geneList               <- factor(as.integer(geneNames %in% myInterestingGenes))
  names(geneList)        <- geneNames
  sampleGO[[1]]         <- new("topGOdata",description = "Simple session", ontology = "BP",allGenes = geneList, geneSel = myInterestingGenes,nodeSize = 10,annot = annFUN.GO2genes,GO2gene=txID2GO)
  
  
  geneNames              <- unique(t2g$ensembl_transcript_id)
  myInterestingGenes     <- unique(as.character(subanno$txID)[mysel])
  myInterestingGenes     <- myInterestingGenes[!is.na(myInterestingGenes)]
  geneList               <- factor(as.integer(geneNames %in% myInterestingGenes))
  names(geneList)        <- geneNames
  sampleGO[[2]]     <- new("topGOdata",description = "Simple session", ontology = "BP",allGenes = geneList, geneSel = myInterestingGenes,nodeSize = 10,annot = annFUN.GO2genes,GO2gene=txID2GO)
  return(sampleGO)
}


CreateSampleGOList <- function(myInterestingGenes){
  sampleGO               <- list()
  geneNames              <- myBG$txID
  geneList               <- factor(as.integer(geneNames %in% myInterestingGenes))
  names(geneList)        <- geneNames
  sampleGO[[1]]         <- new("topGOdata",description = "Simple session", ontology = "BP",allGenes = geneList, geneSel = myInterestingGenes,nodeSize = 10,annot = annFUN.GO2genes,GO2gene=txID2GO)
  
  geneNames              <- unique(t2g$ensembl_transcript_id)
  geneList               <- factor(as.integer(geneNames %in% myInterestingGenes))
  names(geneList)        <- geneNames
  sampleGO[[2]]     <- new("topGOdata",description = "Simple session", ontology = "BP",allGenes = geneList, geneSel = myInterestingGenes,nodeSize = 10,annot = annFUN.GO2genes,GO2gene=txID2GO)
  return(sampleGO)
}

#Enrichment with topGO
getEnrich     <- function(mysampleGO=sampleGOdata1){
  resultFisher            <- runTest(mysampleGO, algorithm = "classic", statistic = "fisher")
  resultFisher.weight01   <- runTest(mysampleGO, algorithm = "weight01", statistic = "fisher")
  
  allRes1.1                 <- GenTable(mysampleGO, classicFisher = resultFisher,weight0Fisher=resultFisher.weight01,orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = noNodes)
  allRes1.2                 <- GenTable(mysampleGO, classicFisher = resultFisher, weight0Fisher=resultFisher.weight01,orderBy = "weight0Fisher", ranksOf = "weight0Fisher", topNodes = noNodes)
  return(list(allRes1.1,allRes1.2))
}

plotEnrich1 <- function(dat=elim){
  dat$pval.unique <- -log10(dat$pval.unique)
  dat$pval.pluri  <- -log10(dat$pval.pluri)
  dat[is.na(dat)] <- 0
  dat             <- dat[sort(dat$pval.unique-dat$pval.pluri,decreasing=T,index.return=T)$ix,]
  L               <- nrow(dat)
  val1            <- as.vector(-dat$pval.unique)
  val2            <- as.vector(dat$pval.pluri)
  
  par(mar=c(3,10,2,3),cex=0.7)
  plot(c(L+3,0),xlim=c(min(-dat$pval.unique)-10,max(dat$pval.pluri)),type = "n",frame=F,yaxt="n",ylab="")
  mp=barplot(height = val1,add = TRUE,axes = FALSE,horiz=T,col="lightsteelblue3")
  barplot(height = val2,add = TRUE,axes = FALSE,horiz=T,col="midnightblue")
  text(x=(min(-dat$pval.unique)-5),y=mp,lab=as.character(dat$Term),cex=0.6)
  mtext(side=1,line=2,text="-log10[P-Value]",cex=0.6)
  legend("top",pch=15,col=c("lightsteelblue3","midnightblue"),bty="n",ncol=2,leg=c("#iso=1","#iso>1"),cex=0.5)
}

MyEnrichPlot <- function(dat=wF,mytitle="BP"){
  
  par(mar=c(3,10,2,3),cex=0.7)
  dat                       <- dat[sort(dat[,1]-dat[,2],decreasing=T,index.return=T)$ix,]
  L                         <- nrow(dat)
  
  val1 <- as.vector(-dat[,1])
  val2 <- as.vector(dat[,2])
  
  plot(c(L+10,0),xlim=c(min(-dat[,1])-10,max(dat)),type = "n",frame=F,yaxt="n",ylab="")
  mp=barplot(height = val1,add = TRUE,axes = FALSE,horiz=T,col="lightsteelblue3")
  barplot(height = val2,add = TRUE,axes = FALSE,horiz=T,col="midnightblue")
  text(x=(min(-dat[,1])-5),y=mp,lab=rownames(dat),cex=0.6)
  mtext(side=1,line=2,text="-log10[P-Value]",cex=0.6)
  mtext(side=3,line=0,text=mytitle,cex=0.6)
  
}

PlotScatterRUD <- function(selS=selShort,selL=selLong){
  #par(mfrow=c(1,1))
  mycols             <- rep("grey",nrow(mRUD))
  mycols[selS] <-rgb(154/255,162/255,197/255)
  mycols[selL] <-rgb(27/255,35/255,83/255)
  plot(mRUD[!(selS|selL),c(2,1)],pch=19,col=mycols[!(selS|selL)],cex=0.3,frame=F,las=1,xlab="",ylab="",cex.axis=0.8)
  points(mRUD[selS|selL,c(2,1)],pch=19,col=mycols[(selS|selL)],cex=0.3)
  abline(v=0,lty=2,col="grey")
  abline(h=0,lty=2,col="grey")
  mtext(side=1,line=2,text="log2 proximal-to-distal poly(A) site ratio",cex=0.8)
  mtext(side=2,line=3,text="log2 proximal-to-distal poly(A) site ratio",cex=0.8)
  mtext(side=2,line=2,text="axonal compartment",cex=0.8)
  mtext(side=1,line=3,text="cell body compartment",cex=0.8)
  subGS   <- as.character(subanno$geneSymbol)
  text(x=-8,y=15,col="lightsteelblue3",labels=paste(length(unique(subGS[selS]))," proximal shifts in axons",sep=""),cex=0.7)
  text(x=8,y=-15,col="midnightblue",labels=paste(length(unique(subGS[selL]))," distal shifts in axons",sep=""),cex=0.7)
  text(x=-8,y=14,col="black",labels=paste("n=",length(unique(subGS))," tandem 3' UTR",sep=""),cex=0.7)
  text(x=-8,y=13,col="black",labels=paste("r=",round(cor(mRUD[,1],mRUD[,2],method="spearman"),digit=2),"(spearman)",sep=""),cex=0.7)
  #  abline(a=0,b=1,lty=2,col="red")
  #  abline(v=0,lty=2,col="red")
  #  abline(h=0,lty=2,col="red")
}

# An isoforms is considered to be reliably expressed if in all the samples the probability of it belonging to the non-expressed class is below 0.05 (less than 5% chance to below to the background in both replicates) -- soft threshold
#OR
# An isoforms is considered to be expressed if in at least one of the replicate the probability of it belonging to the expressed class is above 0.1 (more than 10% chance to belong to the foreground in at least one replicate) -- hard threshold

SelectExpressed <- function(dat=log2(htseq[,1]+1),frac.bg=0.6,frac.fg=0.1){
  bimdens       <- densityMclust(data=dat,G=2)
  x   <- seq(from=0, to=max(dat),length=100)
  if(is.na(bimdens$parameters$variance$sigmasq[2])){
    Lim.fg           <- qnorm(frac.fg,mean=bimdens$parameters$mean[2],sd=sqrt(bimdens$parameters$variance$sigmasq[1]))
    Lim.bg           <- qnorm(frac.bg,mean=bimdens$parameters$mean[1],sd=sqrt(bimdens$parameters$variance$sigmasq[1]))
    hx2              <- dnorm(x,mean=bimdens$parameters$mean[2],sd=sqrt(bimdens$parameters$variance$sigmasq[1]))
  }
  if(!is.na(bimdens$parameters$variance$sigmasq[2])){
    Lim.fg           <- qnorm(frac.fg,mean=bimdens$parameters$mean[2],sd=sqrt(bimdens$parameters$variance$sigmasq[2]))
    Lim.bg           <- qnorm(frac.bg,mean=bimdens$parameters$mean[1],sd=sqrt(bimdens$parameters$variance$sigmasq[1]))
    hx2              <- dnorm(x,mean=bimdens$parameters$mean[2],sd=sqrt(bimdens$parameters$variance$sigmasq[2]))
  }
  
  hx1 <- dnorm(x,mean=bimdens$parameters$mean[1],sd=sqrt(bimdens$parameters$variance$sigmasq[1]))
  hist(dat,breaks=50,col=rgb(0,0,0,alpha=0.2),freq=FALSE,xlab="",ylab="",main="",ylim=c(0,0.4),las=1,xlim=c(0,20))
  lines(x,hx2 , lwd=2, col="blue")
  lines(x,hx1 , lwd=2, col="red")
  abline(v=Lim.fg,col="blue",lty=2)
  abline(v=Lim.bg,col="red",lty=2)
  return(c(Lim.bg,Lim.fg))
}

SelectSurroundingEnd <- function(myGR=myUTR.new,before=10,after=10){
  NEG         <- as.character(strand(myGR))=="-"
  mystart     <- end(myGR)-before
  myend       <- end(myGR)+after
  mystart[NEG]<- start(myGR)[NEG]-after
  myend[NEG]  <- start(myGR)[NEG]+before
  end(myGR)    <- myend
  start(myGR)  <- mystart
  return(myGR)
}

GetFractionOverlap <- function(utrGR=myUTR.1,target=polya_db){
  return(length(unique(subjectHits(findOverlaps(query=target,subject=utrGR,ignore.strand=F))))/length(utrGR))
}

GetFractionRecovered <- function(utrGR=SelectSurroundingEnd(myGR=myUTR,before=10,after=10),target=polya_db){
  known <- utrGR[grep(utrGR$ID,pattern="\\.0"),]
  novel <- utrGR[-grep(utrGR$ID,pattern="\\.0"),] 
  absent<- target[-unique(subjectHits(findOverlaps(query=known,subject=target,ignore.strand=F))),]
  return(length(unique(subjectHits(findOverlaps(query=novel,subject=absent,ignore.strand=F))))/length(absent))
}

Extract_utr <- function(xeno=rn6){
  myEnd.pos       <- as.numeric(end(xeno))
  myStart.neg     <- as.numeric(start(xeno))
  mymaxEnd.pos    <- tapply(myEnd.pos,xeno$transcript_id,max)+10
  myStart.pos     <- mymaxEnd.pos-21
  myminStart.neg  <- tapply(myStart.neg,xeno$transcript_id,min)-10
  myEnd.neg       <- myminStart.neg+21
  POS             <- tapply(as.character(strand(xeno)),xeno$transcript_id,function(x)return(sum(x=="+")!=0))
  Chrom           <- as.character(seqnames(xeno))[match(names(POS),mcols(xeno)$transcript_id)]
  Strand           <- rep("+",length(POS))
  Strand[!POS]     <- "-"
  myStart         <- myStart.pos
  myStart[!POS]   <- myminStart.neg[!POS]
  myEnd           <- mymaxEnd.pos
  myEnd[!POS]     <- myEnd.neg[!POS]
  subxeno         <- GRanges(seqnames =Rle(Chrom),ranges =IRanges(start=myStart,end=myEnd),strand = Rle(Strand),transcript_id=names(POS)
  )
  return(subxeno)
}

MergeGRangesObject <- function(gr1=polya_db,gr2=refseq.rn5){
  Chrom   <- c(seqnames(gr1),seqnames(gr2))
  Strand  <- Rle(c(as.character(strand(gr1)),as.character(strand(gr2))))
  myRanges<- c(ranges(gr1),ranges(gr2))
  return(GRanges(seqnames =Chrom,ranges =myRanges,strand = Strand))
}

AnalysisFractionOverlap <- function(siteoi)return(
  c(100*GetFractionOverlap(utrGR=SelectSurroundingEnd(myGR=myUTR.new,before=10,after=10),target=siteoi),
    100*GetFractionOverlap(utrGR=SelectSurroundingEnd(myGR=myUTR.new,before=50,after=50),target=siteoi),
    100*GetFractionOverlap(utrGR=SelectSurroundingEnd(myGR=myUTR.new,before=100,after=100),target=siteoi),
    100*GetFractionOverlap(utrGR=SelectSurroundingEnd(myGR=myUTR.new,before=150,after=150),target=siteoi),
    100*GetFractionOverlap(utrGR=SelectSurroundingEnd(myGR=myUTR.new,before=200,after=200),target=siteoi)))


AnalysisFractionRecovered <- function(siteoi)return(
  c(100*GetFractionRecovered(utrGR=SelectSurroundingEnd(myGR=myUTR,before=10,after=10),target=siteoi),
    100*GetFractionRecovered(utrGR=SelectSurroundingEnd(myGR=myUTR,before=50,after=50),target=siteoi),
    100*GetFractionRecovered(utrGR=SelectSurroundingEnd(myGR=myUTR,before=100,after=100),target=siteoi)
  )
)



PlotFractionValidated <- function(siteoi=polya_db,name="polyAdb"){
  frac<- c(100*GetFractionOverlap(utrGR=SelectSurroundingEnd(myGR=myUTR.new,before=10,after=10),target=siteoi),
           100*GetFractionOverlap(utrGR=SelectSurroundingEnd(myGR=myUTR.new,before=50,after=50),target=siteoi),
           100*GetFractionOverlap(utrGR=SelectSurroundingEnd(myGR=myUTR.new,before=100,after=100),target=siteoi),
           100*GetFractionOverlap(utrGR=SelectSurroundingEnd(myGR=myUTR.new,before=150,after=150),target=siteoi),
           100*GetFractionOverlap(utrGR=SelectSurroundingEnd(myGR=myUTR.new,before=200,after=200),target=siteoi))
  mp<- barplot(frac,las=1,frame=FALSE,main="",ylab="",ylim=c(0,round(max(frac)+5)),cex.axis=0.7,cex=0.7,col=c("white","#E0E1E0","#989998","#5A5B5A","black"))
  text(x=mp,cex=0.7,y=ceiling(frac)+1,label=round(frac,digits=1))
  mtext(side=2,line=2,cex=0.7,text=paste("% validated by ",name))
}


lappend <- function (lst, ...){
  lst <- c(lst, list(...))
  return(lst)
}
###