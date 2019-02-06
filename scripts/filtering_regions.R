print("Start Filtering Script")
# A. Create GRanges objects
#
#

# A.1 Continue from gross merged segments: groiv1
#     groi    -- Granges of regions of interest (regions which are covered by more than 10 reads)
#	  grga		-- Granges of the the gaps between regions of interest
#	  groiv1	-- Granges of regions of interest which have been merged
#     TestA    	-- Gaps that are less than 150nt long
#	  TestB		-- Gaps that are composed of RepeatRegions + less than 150nt of non-repeat regions
#	  TestC		-- Gaps that are spanned by at least one pair on reads
#     TestF    	-- (TestA|TestB)&TestC
#setwd("/farm/home/luisie01/2011_Keane/tissues/raw/done/fastqfiles/brain/all")
#load("tempFile.RData")


# A.2 Load Granges from various sources:
#   e75p	--	GRange object containing all transcripts from Ensembl  without any filtering
#			i.e. initial gff file download from Ensembl encoded in GRange object.
#   e75f	--	GRange object containing all transcript encoding a PROTEIN with exons from v75 Ensbl.
#   g3utr	-- 	GRange object containing 3'UTR for all protein coding transcripts from Ensembl v75
#			i.e. from end of stop codon until end of the transcript).
#			PLease note however that thes UTRs have introns included.
#   g5utr	--	GRange object which contains 5'UTR for all protein coding transcripts from Ensembl v75
#			i.e. from start of the transcript until start of start codon.
#			PLease note however that thes UTRs have introns included.
#   g3introns 	--      GRange object of introns contained in g3utr (3'UTRs).
#   g5introns	--	GRange object of introns contained in g5utr (5'UTRs).
#       ---takes some time!
print("Load required GRanges: UTRs and transcriptomes form Ensembl.")
load(annotation)

####                THE FITLERING
####
####

# Step INIT == identification of which of the identified regions ( segments) overlaps with known UTRs from ensembl
# Create GR of extensions only
#
print("Step INIT == identification of which of the identified regions ( segments) overlaps with known UTRs from ensembl")
uOver <- findOverlaps(query=groiv1,subject=g3utr,ignore.strand=FALSE)#Use strand information
gover <- groiv1[unique(queryHits(uOver)),]
export.gff(object=gover,con="gover.gff",format="gff")


## Filter n0.1 : Overlaps a single gene (a single 3' UTR in this case)
#
#
FilterSingleGene <- function(gr){
  print("start filter for single gene overlap")
  gOver    <- findOverlaps(query=gr,subject=g3utr,ignore.strand=FALSE)
  idxU     <- queryHits(gOver)
  idxG     <- subjectHits(gOver)
  grouping <-factor(idxU)
  F1    <- unlist(tapply(idxG,grouping,
                         function(x)return(length(levels(factor(g3utr[x,]$gene_id.)))==1)
  )
  )
  gr <- gr[as.numeric(names(F1))[F1],]
  print(paste("N0. of UTRs passing the filter:", sum(F1), sep=""))
  print("Terminated filter for single gene overlap")
  return(gr)
}
print("Filter n0.1 : Overlaps a single gene")
gF1<-FilterSingleGene(gover)
export.gff(object=gF1,con="gF1_init.gff",format="gff")

# Prepare GRanges object composed of UTRs regions only (no gene body):
#   --	one SEGMENT can overlap with several TRANSCRIPT ID (but from same gene)
#		--> due to transcript isoforms that either share stop codons or part of the UTR
#		--> initially we will elongate every overlapped 3'UTRs
#	 	--> however the extensions that overlap 2 stop codons will be removed in subsequent step
#
#   --	one TRANSCRIPT ID can be overlapped by several FRAGMENTS
#		--> in case of very large UTRs
#		--> in this case the fragment that contribute to longest UTRs is retained
print("Prepare GRanges object composed of UTRs regions only (no gene body)")
gOver <- findOverlaps(query=gF1,subject=g3utr,ignore.strand=FALSE)
idxU  <- queryHits(gOver)
idxG  <- subjectHits(gOver)


# --STEP 1: select 3'UTR that are overlapped by fragments
gInit                           <- g3utr[idxG,]
# --STEP 2: elongate 3'UTR when necessary
#            (1) if only one fragment  --> update end of the UTR
#            (2) if several framgments --> update with UTR that is most extreme (depends on strand)
Pos                             <- as.character(strand(gInit))=="+"
end(ranges(gInit))[Pos]         <- end(ranges(gF1))[idxU[Pos]]
start(ranges(gInit))[!Pos]      <- start(ranges(gF1))[idxU[!Pos]]
gF1                             <- gInit
dups				            <- as.character(gInit$X.transcript_id.[duplicated(gInit$X.transcript_id.)])
Pos                             <- as.character(strand(gF1))=="+"
for(i in dups){
  #ix	<-	ranges(gInit[which(gInit$X.transcript_id.==i),])
  ix 	<- 	which(gF1$X.transcript_id.==i)
  if(Pos[ix][1]){
    end(ranges(gF1[ix[1],]))   <- max(end(ranges(gF1[ix,])))
  }
  else{
    start(ranges(gF1[ix[1],])) <- min(start(ranges(gF1[ix,])))
  }
}
gF1 <- gF1[!duplicated(gInit$X.transcript_id.),]
export.gff(object=gF1,con="gF1.gff",format="gff")

#
#
###



## Filter n0.2 : Overlaps zero 5' exon boundaries
#
#    STEP 1 -- check for any overlap with 5'UTR
#    STEP 2 -- check for any overlap with any non-coding transcript

FilterExonBoundary <- function(gr){
  print("start filter no overlap with 5'exon boundaries")
  #Identify which covers a 5'UTR
  myStartCodons <- e75f[e75f$type=="start_codon"]
  gOver 		<- findOverlaps(query=gr,subject=myStartCodons,ignore.strand=FALSE)
  mygene1       <- gr$gene_id.[queryHits(gOver)]
  mygene2       <- myStartCodons$gene_id.[subjectHits(gOver)]
  F2a  		    <- unique(queryHits(gOver)[mygene2!=mygene1])

  #Identify which overlaps with any non-coding transcript
  myNonCoding <- e75p[e75p$source!="protein_coding"]
  gOver 	  <- findOverlaps(query=gr,subject=myNonCoding,ignore.strand=FALSE)
  mygene1     <- gr$gene_id.[queryHits(gOver)]
  mygene2     <- myNonCoding$gene_id.[subjectHits(gOver)]
  F2b  		  <- unique(queryHits(gOver)[mygene2!=mygene1])
  
  F2          <- unique(c(F2a,F2b))
  
  if(length(F2)!=0){
    print(paste("N0. of UTRs not passing the filter:", length(F2), sep=""))
    print("terminated filter no overlap with 5'exon boundaries")
    return(gr[-F2,])
  }
  else{
    print(paste("N0. of UTRs not passing the filter:", length(F2), sep=""))
    print("terminated filter no overlap with 5'exon boundaries")
    return(gr)
  }
}
print("Filter n0.2 : Overlaps zero 5' exon boundaries")
gF2 <- FilterExonBoundary(gF1)
export.gff(object=gF2,con="gF2.gff",format="gff")



#
###

# Filter n0.3 : At least 500 nt of space before next annotated feature (txStart, txStop, or any codon) -- I should consider this filter carefully; I do not do this for this experiment given the strong bias i 3' end
#
#
FilterSpaceBeforeNext <- function(gr){
  print("start filter 500nt space before next annotated feature")
  # Create new GRanges object which start from UTR and terminates 500nt after end of the UTR
  gtemp                       <- gr
  POS                         <- as.character(strand(gr))=="+"
  end(ranges(gtemp))[POS]     <- end(ranges(gr))[POS]+500
  start(ranges(gtemp))[POS]   <- end(ranges(gr))[POS]+1
  start(ranges(gtemp))[!POS]  <- start(ranges(gr))[!POS]-500
  end(ranges(gtemp))[!POS]    <- start(ranges(gr))[!POS]-1
  
  # Then simply check that the UTR does not cover a single feature
  gOver  	  <- findOverlaps(query=gtemp,subject=e75p,ignore.strand=FALSE)
  mygene1     <- gr$gene_id.[queryHits(gOver)]
  mygene2     <- myNonCoding$gene_id.[subjectHits(gOver)]
  F3   	      <- unique(queryHits(gOver)[mygene2!=mygene1])
  
  if(length(F3)!=0){
    print(paste("N0. of UTRs not passing the filter:", length(F3), sep=""))
    print("terminated filter 500nt space before next annotated feature")
    return(gr[-F3,])
  }
  else{
    print(paste("N0. of UTRs not passing the filter:", length(F3), sep=""))
    print("terminated filter 500nt space before next annotated feature")
    return(gr)
  }
}
#print("Filter n0.3 : At least 500 nt of space before next annotated feature (txStart, txStop, or any codon)")
#gF3 <- FilterSpaceBeforeNext(gF2)
#export.gff(object=gF3,con="gF3.gff",format="gff")

#For the moment I prefer to ignore this filter given that the strandedness is already relatively stringent
gF3 <- gF2


#
#
###

# Filter n0.4 : Low-mappable region accounts for < 20% full extension length
#
#

#Low mappable region are tested by TestB&TestC

FilterLowMappable <- function(gr){
  print("start filter low mappable region")
  gaps        <- grga[TestF,]
  gOver       <- findOverlaps(query=gF3,subject=gaps,ignore.strand=FALSE)
  idxU        <- queryHits(gOver)
  idxG        <- subjectHits(gOver)
  if(length(idxU)>0){
    grouping    <- factor(idxU)
    Lev         <- as.numeric(levels(grouping))
    F4          <- replicate(length(Lev),FALSE)
    for(i in 1:length(Lev)){
      F4[i]<- sum(width(ranges(gaps[idxG[which(idxU==Lev[i])],])))/width(ranges(gr[Lev[i],]))>0.2
    }
    if(sum(F4)!=0){
      print(paste("N0. of UTRs not passing the filter:", sum(F4), sep=""))
      print("terminated filter low mappable region")
      return(gr[-Lev[F4],])
    }
    else{
      print(paste("N0. of UTRs not passing the filter:0", sep=""))
      print("terminated filter low mappable region")
      return(gr)
    }
  }
  else{
    print(paste("N0. of UTRs not passing the filter:0", sep=""))
    print("terminated filter low mappable region")
    return(gr)
  }
}


print("Filter n0.4 : Low-mappable region accounts for < 20% full extension length")
gF4  <- FilterLowMappable(gF3)
export.gff(object=gF4,con="gF4.gff",format="gff")


#
#
###

# Filter n0.5:  Overlaps a single STOP codon; here I should select the longest stop codon WITHIN THE GENE FAMILY
#
#
FilterSingleCodon <- function(gr){
  print("start filter single codon")
  gtemp                     <- e75p[which(e75p$type=="stop_codon"),]
  gtemp                     <- reduce(gtemp) #MERGING overlapping STOP codons as several transcripts may share STOP codons
  gtemp                     <- gtemp[match(unique(gtemp$gene_id.),gtemp$gene_id.),]
  POS                       <- as.character(strand(gr))=="+"
  start(ranges(gr))[POS]    <- start(ranges(gr))[POS]-1
  start(ranges(gr))[!POS]   <- start(ranges(gr))[!POS]+1
  gOver                     <- findOverlaps(query=gr,subject=gtemp,ignore.strand = FALSE)
  idxU                      <- queryHits(gOver)
  idxG                      <- subjectHits(gOver)
  grouping                  <- factor(idxU)
  F5                        <- unlist(tapply(idxG,grouping,
                                            function(x)return(length(x)>1)
                                            )
                                     )
  
  if(sum(F5)!=0){
    print(paste("N0. of UTRs not passing the filter:", sum(F5), sep=""))
    print("terminated filter single codon")
    return(gr[-as.numeric(names(F5))[F5],])
  }
  else{
    print(paste("N0. of UTRs not passing the filter:", sum(F5), sep=""))
    print("terminated filter single codon")
    return(gr)
  }
}
gF5 <- FilterSingleCodon(gF4)
print("Filter n0.5:  Overlaps a single STOP codon")
export.gff(object=gF5,con="gF5.gff",format="gff")
#
#
###

# Filter n0.6: Expressed above a minimum threshold:above 1 FPKM for mouse and 1.5 FPKM for human
#
#

# I have a big issue with this one given the normalisation method...
# I guess this step is indeed performed at initial stage. Just bizarre to write in the method as at stage 5.
# In order to test whether this is TRUE regarding FPKM, I would need to only consider a read when mapped in pairs

#
#
###
#print("Filter n0.6: Expressed above a minimum threshold:above 1 FPKM for mouse and 1.5 FPKM for human")
gF6=gF5
#export.gff(object=gF6,con="gF6.gff",format="gff")

# Filter n0.7: Eliminate potential retained introns
#
#              Check for overlap with introns not contained in 3'UTR.
#

FilterIntrons <- function(gr){
  print("start filter introns")
  #Identify which overlaps with transcripts of any kind (BUT one extension can overlap with several transcript from same gene...)
  gtx         <- e75p[which(e75p$type%in%c("transcript")),]
  gexons      <- e75p[which(e75p$type%in%c("start_codon","stop_codon","exon","UTR","CDS")),]
  gin1        <- setdiff(gtx,gexons)
  ginf        <- setdiff(gin1,g3utr)

  #Add gene_id. to introns
  got                            <- findOverlaps(query=ginf,subject=gtx,ignore.strand=FALSE)
  ginf$gene_id.                  <- gtx$gene_id.[subjectHits(got)[match(c(1:length(ginf)),queryHits(got))]]

  # Then simply check that the UTR does not cover a single feature
  gOver                         <- findOverlaps(query=gr,subject=ginf,ignore.strand=FALSE)
  mygene1                       <- gr$gene_id.[queryHits(gOver)]
  mygene2                       <- ginf$gene_id.[subjectHits(gOver)]
  F7   						    <- unique(queryHits(gOver)[mygene2!=mygene1])


  if(length(F7)>=1){
    print(paste("N0. of UTRs not passing the filter F7:", length(F7), sep=""))
    print("terminated filter introns")
    return(gr[-F7,])
  }
  else{
    print(paste("N0. of UTRs not passing the filter:", length(F7), sep=""))
    print("terminated filter introns")
    return(gr)
  }
  
}
print("Filter n0.7: Eliminate potential retained introns")
gF7 <- FilterIntrons(gF6)
export.gff(object=gF7,con="gF7.gff",format="gff")


#
#
###

# Filter n0.8: Extends current model by >=500nt
#
#
FilterExtension <- function(gr,minExt=200){
  print("start filter for significant extension")
  gOver       <- findOverlaps(query=gr,subject=g3utr,ignore.strand=FALSE)
  idxU        <- queryHits(gOver)
  idxG        <- subjectHits(gOver)
  #Compare length of the original against length of the length of the new
  F8          <- unique(idxU[width(ranges(gr[idxU,]))-width(ranges(g3utr[idxG,]))<minExt])
  if(length(F8)>=1){
    print(paste("N0. of UTRs not passing the filter:", length(F8), sep=""))
    print("terminated filter for significant extension")
    return(gr[-F8,])
  }
  else{
    print(paste("N0. of UTRs not passing the filter:", sum(F8), sep=""))
    print("terminated filter for significant extension")
    return(gr)
  }
}
print("Filter n0.8: Extends current model by >=500nt")
gF8 <- FilterExtension(gF7)
export.gff(object=gF8,con="gF8.gff",format="gff")
#
#
###

print("Filtering terminated!")

