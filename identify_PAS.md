##Identification of Polyadenylation signals along the longest expressed 3' UTR

1. Create BED regions of gross fragments

``` R
library(grDevices)
library(Rsamtools)
library(IRanges)
library(GenomicRanges)
library(rtracklayer)

utrs         <- import.gff("./annotation/rn5/GrossFragments.gtf",format="gtf")
start        <- start(utrs)
end          <- end(utrs)
chr          <- as.character(seqnames(utrs))
name         <- utrs$txID
strand       <- as.character(strand(utrs))
out          <- cbind(chr, start,end,name,score=rep(1,length(utrs)),strand)
write.table(x=out,file="./annotation/rn5/GrossFragments.bed", append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", row.names = FALSE,col.names = FALSE)
```

2. Get FASTA sequences using bedtools
```
bedtools getfasta -fi ~/genome/rn5.fa -bed ./annotation/rn5/GrossFragments.bed -s -name -fo ./annotation/rn5/GrossFragments.fasta
```

3. Find PAS element in 3' UTR
```
echo -e  "1\n 2\n 3\n 4\n 5\n 6\n 7\n 8\n 9\n 10\n 11\n 12"| while read IX
do
    OUTDIR=./utrid/GrossSegments/
    FASTA=./annotation/rn5/GrossFragments.fasta
    REGION=./annotation/rn5/GrossFragments.gtf
    ./scripts/findPASelement.R $OUTDIR $IX $FASTA $REGION
    echo " "
    echo " "
done
```

4. Merge all 12 together

``` R
library(grDevices)
library(Rsamtools)
library(IRanges)
library(GenomicRanges)
library(rtracklayer)
foi <- list.files("./utrid/GrossSegments/")
foi <- foi[grep(foi,pattern="A.gtf")]
myPAS <- do.call(lapply(foi,function(x)return(import.gff(con=paste("./utrid/GrossSegments/",x,sep=""),format="gtf"))),what=c)
export.gff(object=myPAS,con="./annotation/rn5/myPAS.gtf",format="gtf")
```
Output stored in [myPAS.gtf](./annotation/rn5/myPAS.gtf)
