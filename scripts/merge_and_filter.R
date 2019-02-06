#!/bin/R-3.0.2/bin/Rscript
library(grDevices)
library(Rsamtools)
library(IRanges)
library(GenomicRanges)


annotation="./annotation/rn5/GRanges_comprehensive_transcriptome_rat_24_nov_2015.RData"
RPM="./annotation/rn5/RpM/"
ensembl="TRUE"

#Load arguments
args            <- commandArgs(TRUE)
wd              <- as.character(args[1])
chr             <- as.character(args[2])
bamfile         <- as.character(args[3])
out             <- as.character(args[4])
strand          <- as.character(args[5])
pairs           <- as.character(args[6])
annotation      <- as.character(args[7])
RPM             <- as.character(args[8])
ensembl         <- as.character(args[9])

print(args)

setwd(wd)
tempdir <-paste("merge_and_filter_",chr,sep="")
dir.create(tempdir)

#Merge neighbouring regions
source("./scripts/mergeNeighbors.R")

#Filter regions
source("./scripts/filtering_regions.R")

#Save files
export.gff(object=gF7,con=out,format="gff")

#Terminated
print("Merging and Filtering finished!")
