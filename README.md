# my3UTR
This repository contains all scripts related to the manuscript [3′UTR cleavage of transcripts localized in axons of sympathetic neurons](https://www.biorxiv.org/content/10.1101/170100v2) which includes the custom pipeline that re-annotates and quantifies 3' end alternative transcripts using RNA-seq data together with down-stream analysis. It is organised as follows:

### A. Bioinformatic pipeline to identify alternative 3' ends from 3'end-seq data
This comprises both the identification of the longest 3' UTR expressed per Ensembl transcript ID as described in file [1_identify_gross_fragments.md](./1_identify_gross_fragments.md), as well as the identification of alternative 3' end as described in file [2_identify_APA.md](./2_identify_APA.md). The output ([L2.gtf](./utrid/APA/L2.gtf)) is a gtf file of alternative 3' UTR isoforms as identified by the pipeline.

### B. Downstream analysis of 3' UTR isoform expression in cell body and distal axons
The down-stream analysis includes 1) the quantification of the expression of each alternative 3' UTR isoform in each of the samples, 2) differential 3’UTR isoforms expression analysis; 3) Gene Ontology enrichment analysis.

## Samples description
mRNA was obtained in duplicates from either cell body or distal axons of rat sympathetic neurons which distal compartment has been exposed to NGF. In this model system, distal axons are separated from the cell bodies by a 1 mm wide Teflon divider, allowing the isolation of mRNA from distinct cellular compartments. Prior to sequencing, mRNA was subjected to two rounds of linear amplification, which similarly to Poly(A)-Seq
(Shepard et al., 2011) led to the accumulation of reads at the 3’ end of transcripts


## Dependencies
TopHat2, BEDTools suite,

## Project inventory

After trimming and QC with fastqc, the RNA-seq data were aligned to the rat genome Rn5 with TopHat using the following parameters:
`tophat2 --mate-inner-dist 1 --max-multihits 1 -N 2 -p 8 -G $geneModel --library-type fr-firststrand -o ${out}${data} $genome ${paths}${file1} ${paths}${file2}` where `$geneModel` is the path to the Ensembl Rn5 (v78) gtf annotation file, `$genome` is the path to the indexed rat genome. Since several run of sequencing have been performed, we then merged all BAM files from same samples together using `samtools merge`

### [0_characterise_tx_read_density.md](./characterise_tx_read_density.md)
This file describes the steps to characterise the 3'--5' end bias in coverage, demonstrating that the data exhibits a similar coverage profile as those from 3'end seq data. 


### [1_identify_gross_fragments.md](./1_identify_gross_fragments.md)
This files describes the steps to extract the longest expressed 3' UTR in each samples, starting from the Ensembl 3' en annotation v78 (Rn4) using a combinations of Bash, Python and R scripts.

Briefly continuously transcribed regions were identified using a sliding window across the genome requiring a minimum coverage of 7 reads in more than 80 positions per window of 100 bp. Neighbouring regions separated by lowmappable regions were merged as described in (Miura et al., 2013). Expressed fragments were associated with matching strand overlapping 3’UTR using Ensembl Rn5 version 78 (v78) (Flicek et al., 2013). Isolated expressed fragments that did not overlap with any feature were associated with the closest 3’UTR if (1) the 3’UTR was <10kb and (2) there were no intervening annotations. We filtered assigned expressed fragments to exclude potential intragenic transcription, overlapping transcripts, and retained introns as described in (Miura et al., 2013). If the expressed fragment continued beyond the end of the annotated 3’UTR, we took the fragment as the new 3’ end.

### [2_identify_APA.md](./2_identify_APA.md)
This file describes the steps to identify expressed alternative 3' UTR isoforms.

Briefly a marked change in the level of coverage in the 3’ to 5’ end direction is expected to occur at the boundaries of alternative polyadenylation sites due to the accumulation of the reads in 3’ termini of the transcripts. To identify alternative 3’UTR isoforms we used Segmentor3IsBack R package (Cleynen et al., 2014) that modelled nucleotide level read coverage along 3’UTRs with a negative binomial distribution and identified segments of spatially coherent coverage along individual 3’UTRs. Prior to coverage segmentation, nucleotidelevel read coverage was smoothed using a running median of width 150 nt. As previously reported, each 3’ end can have multiple cleavage positions in a small window (Tian et al., 2005); we thus clustered together change-points located within 25 nt from one another and selected the most promoter-distal change-point. We applied the algorithm to the raw coverage and log2-scaled coverage of both cell body and axonderived samples. We then merged all 4 annotations (cell body and axon samples, linear and log scale) and clustered 3’ end located within 50 nt distance, selecting the most promoter-distal annotation. We searched the -150 nt to +100 nt region surrounding the 3’ end termini of Ensembl annotated and newly annotated 3’UTR isoforms for 12 canonical and non-canonical PAS motifs listed in PolyA_db (Lee et al., 2007) (AATACA, ATTAAA, TATAAA, AATATA, AATAGA, AGTAAA, AATGAA, ACTAAA, CATAAA, GATAAA, AAGAAA, and AATAAA) using the matchPattern function from the Biostrings R package (Pages et al, 2008.). We tested for the statistical enrichment of the PAS motifs in 3’UTR isoforms using the Fisher’s exact test.

### [3_analysis_3UTR_isoforms.md](./3_analysis_3UTR_isoforms.md)
This files describes the steps to analyse the 3' UTR isoforms usage in cell body and distal axons. This includes 3’UTR isoform quantification and identification of transcripts localized to axons; differential 3’UTR isoforms expression analysis; Gene Ontology enrichment analysis.
