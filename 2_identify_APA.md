
### Step no.2: identification of alternative 3' end along the longest expressed 3' UTR

1. Extract coverage along longest 3'UTR fragment using [extractCovGrossSegments.R](./scripts/extractCovGrossSegments.R)

```bash
OUTDIR=./Coverage/utrCov/
GROSSEG=./annotation/rn5/GrossFragments.gtf
WITHCOR=TRUE
PAIRS=TRUE
FRAC=0.05
PROTOCOL=reverse
awk '{ print $1 }' ./annotation/rn5/chrom_oi.txt | while read CHR
do
    cat ./SAMPLES.txt  | while read SAMPLE
    do
        BAM=${BAMDIR}/{SAMPLE}.bam
        SCRIPT=
        ./scripts/extractCovGrossSegments.R ${BAM} ${OUTDIR} ${CHR} ${SAMPLE} ${GROSSEG} ${WITHCOR} ${PAIRS} ${FRAC} ${PROTOCOL}
    done
done
```

2. Identify APA from coverage profile

```bash
COVDIR=./Coverage/utrCov/
OUTDIR=./utrid/APA/
REGION=./annotation/rn5/GrossFragments.gtf
METHOD=SEG
ANNOTATION=./annotation/GRanges_comprehensive_transcriptome_rat_24_nov_2015.RData
limD=50
winD=151
noSeg=10
isLog=TRUE
minCov=8

awk '{ print $1 }' ./annotation/rn5/chrom_oi.txt | while read CHR
do
    ./scripts/identify_pA_sites.R $OUTDIR $COVDIR $CHR $REGION $METHOD $ANNOTATION $limD $winD $noSeg $isLog $minCov
done
```
Outputs saved in RData objects that can be found in ./utrid/APA.

3. Filter 3' UTR
 
  3.1 Identify PAS along longest 3' UTR expressed fragment as described in [identify_PAS.md](./identify_PAS.md). The output is stored in [myPAS.gtf](./annotation/myPAS.gtf).
  
  3.2 Filter 3' UTR isoforms as described in [FilterUTR.R](./scripts/FilterUTR.R) to remove close neighbours and those overlapping with repeated regions. The output is stored in [L2.gtf](./utrid/APA/L2.gtf).
