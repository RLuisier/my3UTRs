#Identification of Polyadenylation signals along the longest expressed 3' UTR
library(grDevices)
library(Rsamtools)
library(IRanges)
library(GenomicRanges)
library(rtracklayer)

# 1. Create BED regions
utrs         <- import.gff("./annotation/GrossFragments.gtf",format="gtf")
start        <- start(utrs)
end          <- end(utrs)
chr          <- as.character(seqnames(utrs))
name         <- utrs$txID
strand       <- as.character(strand(utrs))
out          <- cbind(chr, start,end,name,score=rep(1,length(utrs)),strand)
write.table(x=out,file="/home/rluisier/data/Riccio/Exp_1/Nov2015/utrid/GrossSegments/GrossFragments.bed", append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", row.names = FALSE,col.names = FALSE)


# B.2 Extract FASTA sequences
bedtools getfasta -fi /farm/home/luisie01/Conservative_reference_annotation_rat/rn5/genome/rn5.fa -bed /farm/home/luisie01/Riccio/Exp_1/WithNewReads/utrid/GrossSegments/GrossFragments.bed -s -name -fo /farm/home/luisie01/Riccio/Exp_1/WithNewReads/utrid/GrossSegments/GrossFragments.fasta

#November
bedtools getfasta -fi /home/rluisier/data/Conservative_reference_annotation_rat/rn5/genome/rn5.fa -bed /home/rluisier/data/Riccio/Exp_1/Nov2015/utrid/GrossSegments/GrossFragments.bed -s -name -fo /home/rluisier/data/Riccio/Exp_1/Nov2015/utrid/GrossSegments/GrossFragments.fasta


# B.2 Find PAS element in 3' UTR
cd /farm/home/luisie01/Riccio/Exp_1/WithNewReads/utrid/GrossSegments
echo -e  "1\n 2\n 3\n 4\n 5\n 6\n 7\n 8\n 9\n 10\n 11\n 12"| while read IX
do
    OUTDIR=/farm/home/luisie01/Riccio/Exp_1/WithNewReads/utrid/GrossSegments/
    SCRIPT=/farm/home/luisie01/Scripts/3-UTR-pipeline/rat_study/findPaElement_improved.R
    FASTA=/farm/home/luisie01/Riccio/Exp_1/WithNewReads/utrid/GrossSegments/GrossFragments.fasta
    REGION=/farm/home/luisie01/Riccio/Exp_1/WithNewReads/utrid/GrossSegments/GrossFragments.gtf
    printf "echo \"  "$SCRIPT" "$OUTDIR" "$IX" "FASTA" "REGION" \"| msub -l nodes=1:ppn=1,walltime=5:00:00 -d $PWD" | sh
    echo " "
    echo " "
done

#November 2015
cd /home/rluisier/data/Riccio/Exp_1/Nov2015/utrid/GrossSegments
echo -e  "1\n 2\n 3\n 4\n 5\n 6\n 7\n 8\n 9\n 10\n 11\n 12"| while read IX
do
    OUTDIR=/home/rluisier/data/Riccio/Exp_1/Nov2015/utrid/GrossSegments/
    SCRIPT=/home/rluisier/Scripts/3-UTR-pipeline/rat_study_replicate/dependencies/findPaElement_improved.R
    FASTA=/home/rluisier/data/Riccio/Exp_1/Nov2015/utrid/GrossSegments/GrossFragments.fasta
    REGION=/home/rluisier/data/Riccio/Exp_1/Nov2015/utrid/GrossSegments/GrossFragments.gtf
    printf "echo \"  "$SCRIPT" "$OUTDIR" "$IX" "$FASTA" "$REGION" \"|  qsub -l h_vmem=10G,tmem=10G,h_rt=5:00:00 -R y -S /bin/bash -wd $PWD -o $OUTDIR/logp.$IX -j y " | sh
    echo " "
    echo " "
done



# B.3 Merge all 12 together
/SAN/luscombelab/general/rluisier/bin/R-3.0.2/bin/R
library(grDevices)
library(Rsamtools)
library(IRanges)
library(GenomicRanges)
library(rtracklayer)
library(Gviz)
foi <- list.files("/farm/home/luisie01/Riccio/Exp_1/WithNewReads/utrid/GrossSegments")
foi <- foi[grep(foi,pattern="A.gtf")]
myPAS <- do.call(lapply(foi,function(x)return(import.gff(con=paste("/farm/home/luisie01/Riccio/Exp_1/WithNewReads/utrid/GrossSegments/",x,sep=""),format="gtf"))),what=c)
export.gff(object=myPAS,con="/farm/home/luisie01/Riccio/Exp_1/WithNewReads/utrid/GrossSegments/myPAS.gtf",format="gtf")

#November version
/SAN/luscombelab/general/rluisier/bin/R-3.0.2/bin/R
library(grDevices)
library(Rsamtools)
library(IRanges)
library(GenomicRanges)
library(rtracklayer)
foi <- list.files("/home/rluisier/data/Riccio/Exp_1/Nov2015/utrid/GrossSegments")
foi <- foi[grep(foi,pattern="A.gtf")]
myPAS <- do.call(lapply(foi,function(x)return(import.gff(con=paste("/home/rluisier/data/Riccio/Exp_1/Nov2015/utrid/GrossSegments/",x,sep=""),format="gtf"))),what=c)

export.gff(object=myPAS,con="/home/rluisier/data/Riccio/Exp_1/Nov2015/utrid/GrossSegments/myPAS.gtf",format="gtf")