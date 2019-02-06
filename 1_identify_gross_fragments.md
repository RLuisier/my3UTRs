1. Get stranded coverage at single nucleotide level [PE_stranded_Coverage.sh](../scripts/PE_stranded_Coverage.sh)
```
OUTDIR=./Coverage
./scripts/PE_stranded_Coverage.sh  -p $OUTDIR -b $BAM -f 2 -n $OUTNAME
```
This script generates 2 files for each sample, for the positive (`$OUTNAME.pos.coverage`) and the negative strand (`$OUTNAME.neg.coverage`).



2. Remove positions which are not covered by at least 7 reads
```
TARDIR=./Coverage/
T=7
for FILES in `ls $TARDIR`; do
    awk '{if($3 > ${T}) print($0)}' ${TARDIR}${FILES}  > ${TARDIR}sub.${FILES}
done
```



3. Split coverage profile per chromosome
```
TARDIR=./Coverage/
cat ./SAMPLES.txt | while read SAMPLE
do
    mkdir ${TARDIR}${SAMPLE}
    mkdir ${TARDIR}${SAMPLE}/pos
    mkdir ${TARDIR}${SAMPLE}/neg

    cd ${TARDIR}${SAMPLE}/pos
    awk '{close(f);f=$1}{print $2 > f".bed"}' ${TARDIR}sub.${SAMPLE}.pos.coverage
    cd ${TARDIR}${SAMPLE}/neg
    awk '{close(f);f=$1}{print $2 > f".bed"}' ${TARDIR}sub.${SAMPLE}.neg.coverage
done
```


4. Identify segments of the genome larger or equal to 100nt covered by at leat 7 reads
```
TARDIR=./Coverage/
cat ./SAMPLES.txt| while read SAMPLE
do
    ls ${TARDIR}${SAMPLE}/pos/*.bed | while read FILE
    do
      ./scripts/findRegions.py -i ${TARDIR}${SAMPLE}/pos/${FILE} -s 100 -t 80 -c 0 -o ${TARDIR}${SAMPLE}/pos/${FILE%.bed}.txt
    done

    ls ${TARDIR}${SAMPLE}/neg/*.bed | while read FILE
    do
      ./scripts/findRegions.py -i ${TARDIR}${SAMPLE}/neg/${FILE} -s 100 -t 80 -c 0 -o ${TARDIR}${SAMPLE}/neg/${FILE%.bed}.txt
    done
done
```

5. Merge expressed regions gapped by low-mappable using script mergeNeighbors.R and filter

```
INDIR=./Coverage
OUTDIR=./utrid/GrossSegments
ANNOTATION=./annotation/rn5/GRanges_comprehensive_transcriptome_rat_24_nov_2015.RData
RPM=./annotation/rn5/RpM/

cat ./SAMPLES.txt| while read SAMPLE
do
    mkdir ${OUTDIR}/${SAMPLE}
    mkdir ${OUTDIR}/${SAMPLE}/neg
    mkdir ${OUTDIR}/${SAMPLE}/pos


    awk '{print $1}' ./annotation/rn5/chrom_oi.txt | while read CHROM
    do
        WD=${OUTDIR}/${SAMPLE}/pos
        ./scripts/merge_and_filter.R ${WD} ${CHROM} ${BAMDIR}/${SAMPLE}.bam ${OUTDIR}/${SAMPLE}/${CHROM}.pos.gff + TRUE ${ANNOTATION} ${RPM} TRUE
        WD=${OUTDIR}/${SAMPLE}/neg
        ./scripts/merge_and_filter.R ${WD} ${CHROM} ${BAMDIR}/${SAMPLE}.bam ${OUTDIR}/${SAMPLE}/${CHROM}.neg.gff - TRUE ${ANNOTATION} ${RPM} FALSE
    done
done
```


6. Merge all gff files
```
cd ./utrid/GrossSegments
cat ./SAMPLES.txt| while read SAMPLE
do
    cd ./utrid/GrossSegments/GrossSegments/${SAMPLE}/neg
    TARDIR=./utrid/GrossSegments/GFF
    cat merge_and*/groi.gff > ${TARDIR}/${SAMPLE}.1.gff
    cat merge_and*/groiv1.gff > ${TARDIR}/${SAMPLE}.2.gff
    cat merge_and*/gover.gff > ${TARDIR}/${SAMPLE}.3.gff
    cat merge_and*/gF1_init.gff > ${TARDIR}/${SAMPLE}.4.gff
    cat merge_and*/gF1.gff > ${TARDIR}/${SAMPLE}.5.gff
    cat merge_and*/gF2.gff > ${TARDIR}/${SAMPLE}.6.gff
    cat merge_and*/gF4.gff > ${TARDIR}/${SAMPLE}.7.gff
    cat merge_and*/gF5.gff > ${TARDIR}/${SAMPLE}.8.gff
    cat merge_and*/gF7.gff > ${TARDIR}/${SAMPLE}.9.gff
    cat merge_and*/gF8.gff > ${TARDIR}/${SAMPLE}.10.gff


    cd ./utrid/GrossSegments/GrossSegments/${SAMPLE}/pos
    cat merge_and*/groi.gff >> ${TARDIR}/${SAMPLE}.1.gff
    cat merge_and*/groiv1.gff >> ${TARDIR}/${SAMPLE}.2.gff
    cat merge_and*/gover.gff >> ${TARDIR}/${SAMPLE}.3.gff
    cat merge_and*/gF1_init.gff >> ${TARDIR}/${SAMPLE}.4.gff
    cat merge_and*/gF1.gff >> ${TARDIR}/${SAMPLE}.5.gff
    cat merge_and*/gF2.gff >> ${TARDIR}/${SAMPLE}.6.gff
    cat merge_and*/gF4.gff >> ${TARDIR}/${SAMPLE}.7.gff
    cat merge_and*/gF5.gff >> ${TARDIR}/${SAMPLE}.8.gff
    cat merge_and*/gF7.gff >> ${TARDIR}/${SAMPLE}.9.gff
    cat merge_and*/gF8.gff >> ${TARDIR}/${SAMPLE}.10.gff
done
```

7. Add extension

```
ANNOTATION=./annotation/rn5/GRanges_comprehensive_transcriptome_rat_24_nov_2015.RData
TARDIR=./utrid/GrossSegments/GFF
cat  ./SAMPLES.txt  | while read SAMPLE
do
    OUTDIR=./utrid/GrossSegments/Extended
    LIM=300
    MAXDIST=16342
    ./scripts/id_far_extensions_general.R ${TARDIR}/${SAMPLE}.2.gff ${OUTDIR} ${SAMPLE} ${LIM} ${TARDIR}/${SAMPLE}.9.gff ${MAXDIST} ${ANNOTATION}
    cp ./utrid/GrossSegments/Extended/new_extensions_${SAMPLE}/gF9_${SAMPLE}.gff ./utrid/GrossSegments/Extended/${SAMPLE}.extended.gff
done
```
Finally merge annotations in R:

``` R
library(grDevices)
library(Rsamtools)
library(IRanges)
library(GenomicRanges)
library(rtracklayer)

setwd("./utrid/GrossSegments/Extended")

# A. Import new annotation
samples <- sort(as.character(read.delim("./SAMPLES.txt",header=F)[,1]))
foi1    <- unlist(lapply(samples,function(x)return(paste("./utrid/GrossSegments/GFF/",x,".6.gff",sep=""))))
foi2    <- unlist(lapply(samples,function(x)return(paste("./utrid/GrossSegments/Extended/",x,".extended.gff",sep=""))))
fois    <- c(foi1,foi2)
myS     <- do.call(lapply(fois,function(x)return(import.gff(con=x, version = c("1"), genome = "rn5"))),what=c)


# B. Import known 3' UTR and annotation of transcript
load("./annotation/rn5/GRanges_comprehensive_transcriptome_rat_24_nov_2015.RData")


# C. Extract the longest 3'UTR isoforms from new annotation
ids          <- as.character(mcols(g3utr)$X.transcript_id)
group        <- as.character(mcols(g3utr)$group)
txID         <- ids[match(as.character(mcols(myS)$group),group)]
myS$txID     <- txID
dL           <- width(myS)-width(g3utr)[match(myS$txID,ids)]

newS         <- GRanges(seqnames = c(seqnames(myS),seqnames(g3utr)),
                       ranges   = c(ranges(myS),ranges(g3utr)),
                       strand   = c(strand(myS),strand(g3utr)),
                       txID     = c(myS$txID,as.character(mcols(g3utr)$X.transcript_id)),
                       dL       = c(dL,rep(0,length(g3utr)))
           )

no.iso       <- as.vector(table(newS$txID))[match(newS$txID,names(table(newS$txID)))]
newS$no.iso  <- no.iso


sub1        <- newS[newS$no.iso==1,]
sub2        <- newS[newS$no.iso>1,]
myTX        <- unique(sub2$txID)
temp        <- cbind(sub2$dL,sub2$txID,c(1:length(sub2)))
theLongest  <- do.call(lapply(myTX, function(x){
                            ix  <- which(temp[,2]==x)
                            sub <- temp[ix,]
                            return(sub[sort(as.numeric(sub[,1]),decreasing=T,index.return=T)$ix[1],3])
                         }
                         ),what=c)

myLongest              <- c(sub1,sub2[as.numeric(theLongest),])
POS                    <- as.character(strand(myLongest))=="+"
end(myLongest)[POS]    <- end(myLongest)[POS] + 100
start(myLongest)[!POS] <- start(myLongest)[!POS] - 100
start(myLongest)[start(myLongest)<1]<-1

export.gff(object=myLongest,con="./annotation/rn5/GrossFragments.gtf",format="gtf")
```
