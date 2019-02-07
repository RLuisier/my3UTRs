### NESTED FUNCTIONS

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

###