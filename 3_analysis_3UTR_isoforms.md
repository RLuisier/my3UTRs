This file describes the downstream analysis of the 3' UTR expression in cell body and axonal compartment.

1. Create GTF file which focus on the last 500 nt fragment of 3' UTR

```R
library(grDevices)
library(Rsamtools)
library(IRanges)
library(GenomicRanges)
library(rtracklayer)

FocusEnd <- function(gr=L3,before=300,after=0){
    POS                <- as.character(strand(gr))=="+"
    grn                <- gr
    start(grn)[POS]  <- end(gr[POS])-before
    end(grn)[POS]    <- end(gr[POS])+after
    end(grn)[!POS]   <- start(gr[!POS])+before
    start(grn)[!POS] <- start(gr[!POS])-after
    start(grn)[start(grn)<0]<-1
    return(grn)
}

L2                 <- import.gff("./utrid/APA/L2.gtf",format="gtf")
L2_500             <- FocusEnd(gr=L2,before=500,after=0)
export.gff(L2_500,"./utrid/APA/L2_500.gtf",format="gtf")

```


2. Extract coverage of the last 500 nt fragment

```
OUTDIR=./Coverage/utrCov/500/
ANNOTATION=./utrid/APA/L2_500.gtf
SF=/SAN/luscombelab/Riccio/Exp_1/FASTQ/P2P3P4/scaling_factor.txt
FRAC=0.05
awk '{print $1}' ./annotation/rn5/chrom_oi.txt | while read CHR
do
    cat ./SAMPLES.txt  | while read SAMPLE
    do
        BAM=${BAMDIR}/SAMPLE.bam
        ./scripts/extractCovFirstSegmentPerChrImproved.R ${BAM} ${OUTDIR} ${CHR} ${SAMPLE} ${ANNOTATION} ${SF} ${FRAC}
    done
done

```

3. Collect coverage
Described in file [import_coverage.R](./scripts/import_coverage.R), the outputs are the following:

*[myCov500](./data/myCov500.tab) the raw count table. 

*[anno_ngf](./data/anno_ngf.tab)

*[Lngf.gtf](./annotation/rn5/Lngf.gtf) the gtf file of the 3' UTR isoforms covered by at least 10 reads.

*[Lngf_sub.gtf](./annotation/rn5/Lngf_sub.gtf) the same as Lngf.gtf however where 3'UTR which are within 500nt of each other are filtered out, keeping the longest 3' UTR isoform. Furthermore isoforms of length smaller than 10 nt were remove.


4. Characterise the output of the pipeline


