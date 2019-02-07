source("http://bioconductor.org/biocLite.R")
require("GO.db")
require("limma")
require("topGO")
require("biomaRt")
require("org.Rn.eg.db")
library(grDevices)
library(Rsamtools)
library(IRanges)
library(GenomicRanges)
library(rtracklayer)
library(geneplotter)
require("multtest")
require("mclust")
rm(list = ls())

source("./scripts/nested_functions.R")

# A. Load data
load("./annotation/rn5/GRanges_comprehensive_transcriptome_rat_24_nov_2015.RData")
myUTR    <- import.gff("./annotation/rn5/Lngf_sub.gtf",format="gtf")
anno_ngf     <- read.table("./data/anno_ngf.tab",header=T,sep="\t")

# B. Remove redundant transcripts and select reliably expressed genes
id       <- paste(start(myUTR),end(myUTR),strand(myUTR),sep=".")#2215 duplicated
myUTR    <- myUTR[-which(duplicated(id)),]
anno_ngf <- anno_ngf[match(myUTR$ID,anno_ngf$uniqueID),]

par(mfrow=c(2,2))
mydat                <- anno_ngf[which(anno_ngf$is.conservative),match(c("NGF.axon.1.raw","NGF.axon.2.raw","NGF.cb.1.raw","NGF.cb.2.raw"),colnames(anno_ngf))]
mydatall             <- anno_ngf[,match(c("NGF.axon.1.raw","NGF.axon.2.raw","NGF.cb.1.raw","NGF.cb.2.raw"),colnames(anno_ngf))]
for(i in c(1:ncol(mydat))){
  mydat[,i] <- log2(mydat[,i]+1)
}
lims      <- apply(mydat,2,function(Z)2^SelectExpressed(dat=Z,frac.bg=0.95,frac.fg=0.1))
tempsel   <- apply(mydatall,2,function(Z)return(Z>=mean(lims[1,c(1,2)])))
soft.sel  <- cbind(tempsel[,1]&tempsel[,2],tempsel[,3]&tempsel[,4])
tempsel   <- do.call(what=cbind,lapply(c(1:4),function(Z)return(mydatall[,Z]>=lims[2,Z])))
hard.sel  <- cbind(tempsel[,1]|tempsel[,2],tempsel[,3]|tempsel[,4])
final.sel <- cbind(soft.sel[,1]|hard.sel[,1],soft.sel[,2]|hard.sel[,2])

NGF.axon.is.expressed.iso                   <- final.sel[,1]
NGF.cb.is.expressed.iso                     <- final.sel[,2]
anno_ngf$NGF.axon.is.expressed.iso          <- final.sel[,1]
anno_ngf$NGF.cb.is.expressed.iso            <- final.sel[,2]

# C. Sub-selection of those 3' UTR isoforms for down-stream analysis of differential expression
no.genes.expressed.neurons              <- length(unique(anno_ngf$geneSymbol[NGF.axon.is.expressed.iso|NGF.cb.is.expressed.iso]))
no.genes.expressed.axons                <- length(unique(anno_ngf$geneSymbol[NGF.axon.is.expressed.iso]))

#txID.expression
axons                                   <- unique(anno_ngf$txID[anno_ngf$NGF.axon.is.expressed.iso])
cb                                      <- unique(anno_ngf$txID[anno_ngf$NGF.cb.is.expressed.iso])
neurons                                 <- union(axons,cb)
cb.only                                 <- setdiff(cb,axons)
axon.only                               <- setdiff(axons,cb)
both                                    <- intersect(axons,cb)

#uniqueID.expression
axons                                  <- anno_ngf$uniqueID[anno_ngf$NGF.axon.is.expressed.iso]
cb                                     <- anno_ngf$uniqueID[anno_ngf$NGF.cb.is.expressed.iso]
neur.neur                              <- as.character(anno_ngf$uniqueID)[anno_ngf$uniqueID%in%intersect(cb,axons)&anno_ngf$txID%in%both]
neur.cb                                <- as.character(anno_ngf$uniqueID)[anno_ngf$uniqueID%in%setdiff(cb,axons)&anno_ngf$txID%in%both]
neur.ax                                <- as.character(anno_ngf$uniqueID)[anno_ngf$uniqueID%in%setdiff(axons,cb)&anno_ngf$txID%in%both]
ax.ax                                  <- as.character(anno_ngf$uniqueID)[anno_ngf$uniqueID%in%setdiff(axons,cb)&anno_ngf$txID%in%axon.only]
cb.cb                                  <- as.character(anno_ngf$uniqueID)[anno_ngf$uniqueID%in%setdiff(cb,axons)&anno_ngf$txID%in%cb.only]

temp       <-cbind(anno_ngf$uniqueID%in%neur.neur,anno_ngf$uniqueID%in%neur.cb,anno_ngf$uniqueID%in%neur.ax,anno_ngf$uniqueID%in%cb.cb,anno_ngf$uniqueID%in%ax.ax)
iso.class <- rep("no",length(anno_ngf$uniqueID))
iso.class[anno_ngf$uniqueID%in%neur.neur] <- "neur.neur"
iso.class[anno_ngf$uniqueID%in%neur.cb]   <- "neur.cb"
iso.class[anno_ngf$uniqueID%in%neur.ax]   <- "neur.ax"
iso.class[anno_ngf$uniqueID%in%cb.cb]     <- "cb.cb"
iso.class[anno_ngf$uniqueID%in%ax.ax]     <- "ax.ax"
anno_ngf$iso.class <- iso.class


#Create Backgroun with what is expressed in neurons
mytx <- unique(as.character(anno_ngf$txID[anno_ngf$NGF.axon.is.expressed.iso|anno_ngf$NGF.cb.is.expressed.iso]))
mygs <- unique(as.character(anno_ngf$geneSymbol[anno_ngf$NGF.axon.is.expressed.iso|anno_ngf$NGF.cb.is.expressed.iso]))#12'529
myBG <- data.frame(txID=mytx,GS=c(mygs,rep(NA,length(mytx)-length(mygs))))

#Remove those transcript ID which have only one isoforms
anno_ngf     <- anno_ngf[anno_ngf$NGF.axon.is.expressed.iso|NGF.cb.is.expressed.iso,]
temp         <- as.data.frame(table(as.character(anno_ngf$txID)))
no.iso       <- temp[match(anno_ngf$txID,as.character(temp[,1])),2]
anno_ngf     <- anno_ngf[anno_ngf$no.iso>1,]
names(myUTR) <- myUTR$ID
ix           <- match(anno_ngf$uniqueID,names(myUTR))
myUTR        <- myUTR[ix[!is.na(ix)],]
anno_ngf     <- anno_ngf[!is.na(ix),]


#Check that the closest neighbour is not within 300 nt
#iso_ordered
tempL                     <- anno_ngf$newL
names(tempL)              <- anno_ngf$uniqueID
tempG                     <- as.factor(as.character(anno_ngf$txID))
myord                     <- tapply(tempL,INDEX=tempG,function(x)return(cbind(names(x)[sort(as.numeric(as.character(x)),index.return=T,decreasing=F)$ix],c(1:length(x)))))
test                      <- do.call(what=rbind,args=myord)
iso_ordered               <- as.numeric(test[match(anno_ngf$uniqueID,test[,1]),2])
#nextIsoform
tempid                    <- paste(anno_ngf$txID,iso_ordered,sep=".")
next.iso                  <- paste(anno_ngf$txID,(as.numeric(as.character(iso_ordered))+1),sep=".")
nextID                    <- anno_ngf$uniqueID[match(next.iso,tempid)]
#distTonext
distonext                 <- rep(NA,nrow(anno_ngf))
sel                       <- !is.na(nextID)
d1                        <- as.numeric(as.character(anno_ngf$newL))[sel]
d2                        <- as.numeric(as.character(anno_ngf$newL))[match(nextID[sel],anno_ngf$uniqueID)]
distonext[sel]            <- d2-d1
#Fproxi
Fproxi                    <- which(distonext<300)#should be none

# Get ordering of the isoforms
tempL                     <- anno_ngf$newL
names(tempL)              <- anno_ngf$uniqueID
tempG                     <- as.factor(as.character(anno_ngf$txID))
myord                     <- tapply(tempL,INDEX=tempG,function(x)return(cbind(names(x)[sort(as.numeric(as.character(x)),index.return=T,decreasing=F)$ix],c(1:length(x)))))
test                      <- do.call(what=rbind,args=myord)
test                      <- test[match(anno_ngf$uniqueID,test[,1]),]
temp                      <- as.data.frame(table(as.character(anno_ngf$txID)))
no.iso                    <- temp[match(anno_ngf$txID,as.character(temp[,1])),2]
anno_ngf$no.iso           <- no.iso
anno_ngf$iso_ordered      <- test[,2]
tempL                     <- as.character(anno_ngf$iso_ordered)
names(tempL)              <- as.character(anno_ngf$uniqueID)
imID1                     <- tapply(tempL,INDEX=as.factor(as.character(anno_ngf$txID)),function(x)return(names(x)[x=="1"]))
imIDp                     <- as.character(imID1[match(anno_ngf$txID,names(imID1))])
anno_ngf$imIDp            <- imIDp

# Get the selection on which to focus for the analysis
# At least one of the pair (distal or proximal) must be detected in axons
sela1                <- (anno_ngf$NGF.axon.is.expressed.iso|anno_ngf$NGF.axon.is.expressed.iso[match(anno_ngf$imIDp,anno_ngf$uniqueID)])
# At least one of the pair (distal or proximal) must be detected in cb
sela2                <- (anno_ngf$NGF.cb.is.expressed.iso|anno_ngf$NGF.cb.is.expressed.iso[match(anno_ngf$imIDp,anno_ngf$uniqueID)])
# Remove all the isoform I0
sela3                 <- anno_ngf$uniqueID!=anno_ngf$imIDp
#Proximal is expressed in at least one of the 2 samples
sela4                 <- (anno_ngf$NGF.axon.is.expressed.iso|anno_ngf$NGF.cb.is.expressed.iso)[match(anno_ngf$imIDp,anno_ngf$uniqueID)]
#Distal is expressed in at least one of the 2 samples
sela5                 <- (anno_ngf$NGF.axon.is.expressed.iso|anno_ngf$NGF.cb.is.expressed.iso)

sela                  <- sela1&sela2&sela3&sela4&sela5
subanno               <- anno_ngf[sela,]


# D. IDENTIFICATION OF DIFFERENTIALLY EXPRESSED ISOFORMS
#Proximal to distal site usage AND fisher test
ix1                  <- c(1:nrow(anno_ngf))
ix2                  <- match(as.character(anno_ngf$imIDp),as.character(anno_ngf$uniqueID))
mydat                <- anno_ngf[,match(c("NGF.axon.1.raw","NGF.axon.2.raw","NGF.cb.1.raw","NGF.cb.2.raw"),colnames(anno_ngf))]
for(i in c(1:ncol(mydat))){mydat[,i]<-log2(1+mydat[,i])}
rownames(mydat)      <- anno_ngf$uniqueID
myProximal           <- apply(mydat,2,function(x)return(x[ix2]))
myDistal             <- apply(mydat,2,function(x)return(x[ix1]))
rownames(myProximal) <-rownames(myDistal)<-anno_ngf$uniqueID[ix1]
RUD                  <- do.call(lapply(c(1:4),function(x)return(myProximal[,x]-myDistal[,x])),what=cbind)
rownames(RUD)        <- anno_ngf$uniqueID[ix1]
colnames(RUD)        <- c("axon.1","axon.2","cb.1","cb.2")
RUD                  <- RUD[sela,]#8'853
sdRUD                 <- t(apply(RUD,1,function(x)return(tapply(x,INDEX=factor(c("axon","axon","cb","cb")),FUN=sd))))
mRUD                  <- t(apply(RUD,1,function(x)return(tapply(x,INDEX=factor(c("axon","axon","cb","cb")),FUN=mean))))
dRUD                  <- mRUD[,2]-mRUD[,1]

#Compare replicates
par(mfrow=c(2,2))
smoothScatter(RUD[,1],RUD[,2],frame=FALSE,xlab="ax1",ylab="ax2")
mtext(side=3,line=0,text=cor(RUD[,1],RUD[,2],method="spearman"))
smoothScatter(RUD[,3],RUD[,4],frame=FALSE,xlab="cb.1",ylab="cb.2")
mtext(side=3,line=0,text=cor(RUD[,3],RUD[,4],method="spearman"))


#Fisher test on sum of RAW count data
sumRUD               <- cbind(
  apply(anno_ngf[,match(c("NGF.axon.1.raw","NGF.axon.2.raw"),colnames(anno_ngf))],1,sum),
  apply(anno_ngf[,match(c("NGF.cb.1.raw","NGF.cb.2.raw"),colnames(anno_ngf))],1,sum)
)
rownames(sumRUD)<- anno_ngf$uniqueID
colnames(sumRUD)<- c("axon","cb")
test.apa <- function(ix.proximal=ix1[1],ix.distal=ix2[1]){
  return(fisher.test(round(cbind(sumRUD[ix.proximal,],sumRUD[ix.distal,])))$p.value)
}
fisherRUD   <- unlist(lapply(c(1:nrow(sumRUD)),function(x)return(test.apa(ix.proximal=ix1[x],ix.distal=ix2[x]))))


my.proximal        <- sumRUD[ix2,]
my.distal          <- sumRUD[ix1,]
rel.proximal.usage <- cbind(my.proximal[,1]/(my.distal[,1]+my.proximal[,1]),
                            my.proximal[,2]/(my.distal[,2]+my.proximal[,2]))

colnames(rel.proximal.usage) <- c("axon","cb")
rownames(rel.proximal.usage)  <- rownames(sumRUD)

sumRUD                 <- sumRUD[sela,]
fisherRUD              <- fisherRUD[sela]
rel.proximal.usage     <- rel.proximal.usage[sela,]
diff.rel.proximal.usage<- rel.proximal.usage[,2]-rel.proximal.usage[,1]
# Correct for multitest (remove non-important ones)
padjRUD              <-  p.adjust(fisherRUD,method="fdr")
mRUD                 <-  as.data.frame(mRUD)

# Selection
id.proxi             <-  anno_ngf$imIDp[match(names(dRUD),anno_ngf$uniqueID)]
seld1                <-  anno_ngf$NGF.axon.is.expressed.iso[match(id.proxi,anno_ngf$uniqueID)]
seld2                <-  anno_ngf$NGF.axon.is.expressed.iso[match(names(dRUD),anno_ngf$uniqueID)]
# Significant P-value of difference
selc                 <- padjRUD<0.01

#Relative Proximal to distal site usage
mydat                <- anno_ngf[,match(c("NGF.axon.1.raw","NGF.axon.2.raw","NGF.cb.1.raw","NGF.cb.2.raw"),colnames(anno_ngf))]
ix.distal            <- c(1:nrow(anno_ngf))
ix.proximal          <- match(as.character(anno_ngf$imIDp),as.character(anno_ngf$uniqueID))
psi                  <- mydat
for(i in c(1:ncol(psi))){
  psi[,i]                                                              <- mydat[ix.proximal,i]/(mydat[ix.distal,i]+mydat[ix.proximal,i])
  psi[is.na(psi[,i]),i]<-0
}
PUD                 <- psi
rownames(PUD)       <- anno_ngf$uniqueID
colnames(PUD)       <- c("axon.1","axon.2","cb.1","cb.2")
mPUD                <- t(apply(PUD,1,function(x)return(tapply(x,INDEX=factor(c("axon","axon","cb","cb")),FUN=mean))))
PUD                 <- PUD[sela,]#8'853
mPUD                <- mPUD[sela,]
diffPUD             <-  mPUD[,"axon"]-mPUD[,"cb"]
mPUD <- data.frame(mPUD)

#Selection of genes of interest
opp.1 <- mRUD$axon>0&mRUD$cb<0&mPUD$axon>0.5&mPUD$cb<0.5
opp.2 <- mRUD$axon<0&mRUD$cb>0&mPUD$axon<0.5&mPUD$cb>0.5
Lim1  <- 0.2
Lim2  <- 0.8

selRUD1             <- list()
IX=1
selRUD1[[IX]]       <- padjRUD<0.01&dRUD<(-1)&seld1#980; proximal shift so proximal is expressed in axons
IX=IX+1
selRUD1[[IX]]       <- padjRUD<0.01&dRUD>(1)&seld2#1065; distal shift so distal expressed in axons

IX=IX+1
selRUD1[[IX]]      <- padjRUD<0.01&dRUD<(-1)&seld1&diffPUD>(0.15)#689-proximal shift so proximal is expressed in axons
IX=IX+1
selRUD1[[IX]]      <- padjRUD<0.01&dRUD>(1)&seld2&diffPUD<(-0.15)#737-distal shift so distal expressed in axons

IX=IX+1
selRUD1[[IX]]      <- padjRUD<0.01&dRUD<(-1)&seld1&diff.rel.proximal.usage<(-0.15)&rel.proximal.usage[,2]<=Lim1#178, less that 20% usage of the proximal in CB; therefore can be considered as absent
IX=IX+1
selRUD1[[IX]]      <- padjRUD<0.01&dRUD>(1)&seld2&diff.rel.proximal.usage>(0.15)&rel.proximal.usage[,2]>=Lim2#86,less than 20% usage of the distal in CB; therefore can be considered as absent


write.csv(cbind(subanno[selRUD1[[5]],c("uniqueID","geneSymbol","imIDp","newL")],mPUD[selRUD1[[5]],]),"./differential_expression/proximal_shifts.csv")
write.csv(cbind(subanno[selRUD1[[6]],c("uniqueID","geneSymbol","imIDp","newL")],mPUD[selRUD1[[6]],]),"./differential_expression/distal_shifts.csv")


#Inset Figure 2e
my.no <- unlist(lapply(selRUD1,function(Z)return(length(unique(subanno$txID[Z])))))
mp<- barplot(my.no[c(5,6,9,10)],las=1)
mtext(side=3,line=0,text=my.no[c(5,6,9,10)],at=mp)
mylist <- lapply(selRUD1,function(Z)return(unique(subanno$txID[Z])))


s1p <- anno_ngf[which(anno_ngf$uniqueID%in%subanno$imIDp[selRUD1[[3]]]),match(c("uniqueID","txID","geneSymbol","NGF.cb.1.raw","NGF.cb.2.raw","NGF.axon.1.raw","NGF.axon.2.raw"),colnames(anno_ngf))]
s2p <- anno_ngf[which(anno_ngf$uniqueID%in%subanno$imIDp[selRUD1[[5]]]),match(c("uniqueID","txID","geneSymbol","NGF.cb.1.raw","NGF.cb.2.raw","NGF.axon.1.raw","NGF.axon.2.raw"),colnames(anno_ngf))]

s1d  <- subanno[which(selRUD1[[4]]),match(c("uniqueID","txID","geneSymbol","NGF.cb.1.raw","NGF.cb.2.raw","NGF.axon.1.raw","NGF.axon.2.raw"),colnames(subanno))]
s2d  <- subanno[which(selRUD1[[6]]),match(c("uniqueID","txID","geneSymbol","NGF.cb.1.raw","NGF.cb.2.raw","NGF.axon.1.raw","NGF.axon.2.raw"),colnames(subanno))]

write.csv(s2p,"./differential_expression/candidates_axonal_remodelling.csv")
write.csv(s2d,"./differential_expression/candidates_faciliated_transport.csv")

sampleGO.RUD1      <- lapply(selRUD1[c(3,4)],CreateSampleGO)
enrich.RUD1        <- lapply(sampleGO.RUD1,FUN=function(X)return(lapply(X,getEnrich)))
F1.w0     <- CompareBP_improved(enr1=enrich.RUD1[[1]][[2]][[2]],enr2=enrich.RUD1[[2]][[2]][[2]],PLOT=TRUE,no=20,lab2="long",lab1="short",coi="weight0Fisher")
F1.fisher <- CompareBP_improved(enr1=enrich.RUD1[[1]][[2]][[1]],enr2=enrich.RUD1[[2]][[2]][[1]],PLOT=TRUE,no=20,lab2="long",lab1="short",coi="classicFisher")
write.csv(file="./GOenrichment/F1.w0.r.csv",F1.w0)
write.csv(file="./GOenrichment/F1.fisher.r.csv",F1.fisher)


layout(matrix(c(1,1,1,1), 2, 2, byrow = FALSE))
par(mar=c(3,3,3,3))
final_enrich <- read.csv("./GOenrichment/F1.w0.r_f.csv")
#Figure 2e 
PlotScatterRUD(selS=selRUD1[[3]],selL=selRUD1[[4]])
vec1 <- final_enrich[,2]
vec2 <- final_enrich[,3]
names(vec1)<-names(vec2)<-as.character(final_enrich[,1])
vec1  <- sort(vec1,decreasing=F)
vec2  <- sort(vec2,decreasing=F)
vec1  <- vec1[vec1>=-log10(0.05)]
vec2  <- vec2[vec2>=-log10(0.05)]
#Figure 2f
barplot(vec1,horiz=TRUE,col=rgb(154/255,162/255,197/255),las=1)
mtext(side=3,line=0,text=foi[IX])
mtext(side=1,line=2,text="-log10(P-Value)")
barplot(vec2,horiz=TRUE,col=rgb(27/255,35/255,83/255),las=1)
mtext(side=1,line=2,text="-log10(P-Value)")

dat <- c(my.no,length(unique(subanno$geneSymbol[selRUD1[[3]]])),length(unique(subanno$geneSymbol[selRUD1[[4]]])))
names(dat)[c(4,5)] <- c("proximal","distal")
mp  <- barplot(dat,las=1,ylab="no.GS")
mtext(side=3,line=0,text=dat,at=mp)

