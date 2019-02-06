#!/bin/R-3.0.2/bin/Rscript
library(grDevices)
library(Rsamtools)
library(IRanges)
library(GenomicRanges)
library(rtracklayer)

# A. Load arguments
args       <- commandArgs(TRUE)
bamfile    <- args[1]
outdir     <- args[2]
chr        <- args[3]
mysample   <- args[4]
annotation <- args[5]
SF         <- args[6]
frac       <- as.numeric(as.character(args[7]))
print(args)


out                      <- paste(mysample,chr,"txt",sep=".")
out.norm                 <- paste(mysample,chr,"norm","txt",sep=".")
# B. Load annotation file and extract Exons

# B.1 Import 3' UTR new annotation with focus on latest 750 nt; contains object 'bothGRSfocus' that is a GRange object
myUTR               <- import.gff(annotation,format="gtf",asRangedData=FALSE)
myUTR               <- myUTR[as.character(seqnames(myUTR))==chr,]
print("loaded annotation")

# B.2 Load BAM file
chrlength       <- scanBamHeader(bamfile)[[1]]$targets[[chr]]
chrgr           <- GRanges(seqnames=chr,ranges=IRanges(start=1, end=chrlength))
gr              <- readGAlignmentPairsFromBam(file=bamfile, param=ScanBamParam(which=chrgr))
print("BAM loaded")
setwd(outdir)

# C. Compute coverage

# C.1 Select only those reads from "+" strand
m1              <- first(gr)
m2              <- last(gr)
sel1            <- which(strand(m1)=="-")
sel2            <- which(strand(m2)=="+")
if(sum(sel1!=sel2)>0){print("hey then there is an issue!")}
m1              <- m1[sel1]
m2              <- m2[sel2]
mtot            <- c(m1,m2)
covpos          <- assays(summarizeOverlaps(features=myUTR, reads=mtot,ignore.strand=FALSE,inter.feature=FALSE))$counts

#I should maybe ignore the strand

# C.2 Select only those reads from "-" strand
m1              <- first(gr)
m2              <- last(gr)
sel1            <- which(strand(m1)=="+")
sel2            <- which(strand(m2)=="-")
if(sum(sel1!=sel2)>0){print("hey then there is an isse!")}
m1              <- m1[sel1]
m2              <- m2[sel2]
mtot            <- c(m1,m2)
covneg          <- assays(summarizeOverlaps(features=myUTR, reads=mtot,ignore.strand=FALSE,inter.feature=FALSE))$counts


# C.3 Select only counts from appropriate strand and substract those from opposite strand
POS                 <- which(as.character(strand(myUTR))=="+")
NEG                 <- which(as.character(strand(myUTR))=="-")
covfinal            <- covpos
covfinal[NEG]       <- covneg[NEG]
#Correction for the amount of contamination
covfinal[POS]       <- covfinal[POS]-(frac/(1-frac^2))*(covneg[POS]-frac*covfinal[POS])
covfinal[NEG]       <- (covfinal[NEG]-frac*covpos[NEG])/(1-frac^2)
covfinal[covfinal<0]<-0
rownames(covfinal)  <- as.character(myUTR$ID)



# C.4 Normalise with scale factor
if(!is.na(SF)){
    scalef   <- read.table(SF,header=F)
    scalef   <- as.numeric(as.character(scalef[match(mysample,scalef[,1]),3]))
    covfinalnorm <- covfinal*scalef
}

# D. Save file
write.table(x=covfinal,file=out,append=F,quote=F,row.names=T,col.names=F)
write.table(x=covfinalnorm,file=out.norm,append=F,quote=F,row.names=T,col.names=F)
print("Finished!")


