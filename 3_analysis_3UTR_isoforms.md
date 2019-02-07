This file describes the downstream analysis of the 3' UTR expression in cell body and axonal compartment.

### 1. Create GTF file which focus on the last 500 nt fragment of 3' UTR

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


### 2. Extract coverage of the last 500 nt fragment

```bash
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

### 3. Collect coverage

Described in file [import_coverage.R](./scripts/import_coverage.R), the outputs are the following:
  * [myCov500.tab](./data/myCov500.tab) the raw count table.
  * [anno_ngf.tab](./data/anno_ngf.tab) a table containing both the raw counts together with additional information about the 3' UTR isoforms.
  * [Lngf.gtf](./annotation/rn5/Lngf.gtf) the gtf file of the 3' UTR isoforms covered by at least 10 reads.
  * [Lngf_sub.gtf](./annotation/rn5/Lngf_sub.gtf) the same as Lngf.gtf however where 3'UTR which are within 500nt of each other are filtered out, keeping the longest 3' UTR isoform. Furthermore isoforms of length smaller than 10 nt were remove.


### 4. Validate the pipeline by comparing with compiled polyA atlas

The reliability of the novel 3' end annotations was assessed by checking them against a compiled polyA Atlas. Assembling the Polyadenylation Site Atlas: first 3' UTR annotations from RefSeq (Rn5 and Rn6) and Ensembl (Rn6) are combined, where overlapping regions are merged. Given that rat genome is poorly annotated compared to other species such as human and mouse, and because 3' end are highly conserved across species, we futher combined existing rat 3' UTR annotation with RefSeq annotation from other species (xenoRefseq). Then, to generate a comprehensive atlas of rat PAS in 3' UTR, polyadenylation annotations from PolyA_DB 2 (37), and APADB (39) were used together with the polyadyenlation atlas from unibasel (both human and mouse). Finally PAS obtained from rat 3'-Seq data (Derti et al.) were also used to expand the repertoire of PAS. Genome annotation files from Rn6 assembly or other species have been converted to Rn5 using Crossmap. Downstream analysis is described in [compare_with_polyA_atlas.R](./scripts/compare_with_polyA_atlas.R).

### 5. Characterise the output of the pipeline

Characterise the output of the pipeline in terms of number of newly identified 3' UTR isoforms, 3' UTR expression in each compartment, number of isoforms expressed. Described in file [characterise_output_pipeline.R](./scripts/characterise_output_pipeline.R).

### 6. Differential 3’UTR isoforms expression analysis

Perform differential 3' UTR expression analysis between cell body and axonal compartment leading to the identification of candidate 3' UTR for axonal remodelling. Described in file [differential_expression_between_compartments.R](./scripts/differential_expression_between_compartments.R).

