#!/bin/bash

# Usage: ./PE_stranded_Coverage.sh -p [OUTDIR FOR TEMP FILES] -b [PATH TO BAM FILE] -f [PROTOCOL USED: 1 IF ++/-- AND 2 OTHERWISE] -n [PATH TO OUTPUT FILE]

while getopts p:b:f:n: option
do
        case "${option}"
        in
                p) OUTDIR=${OPTARG};;
                b) BAM=${OPTARG};;
                f) PROTOCOL=${OPTARG};;
                n) OUT=${OPTARG};;
        esac
done

mkdir $OUTDIR

# 1. Separate mates
if [ -f $OUTDIR/m1.bam ]; then
    echo No need to create BAM files
else
    echo Create BAM files
    samtools view -b -f 0x0040 -h -o $OUTDIR/m1.bam $BAM
    samtools view -b -f 0x0080 -h -o $OUTDIR/m2.bam $BAM
fi


echo "finished separation between left and right reads"
# 2. Calculate coverage separately for each strand

###
# M1
###

# Watson-strand:
genomeCoverageBed -ibam $OUTDIR/m1.bam -dz -split -strand + > $OUTDIR/m1.pos.coverage
# Crick-strand:
genomeCoverageBed -ibam $OUTDIR/m1.bam -dz -split -strand - > $OUTDIR/m1.neg.coverage

echo "finished reading coverage m1"

###
# M2
###

# Watson-strand:
genomeCoverageBed -ibam $OUTDIR/m2.bam -dz -split -strand + > $OUTDIR/m2.pos.coverage
# Crick-strand:
genomeCoverageBed -ibam $OUTDIR/m2.bam -dz -split -strand - > $OUTDIR/m2.neg.coverage

echo "finished reading coverage m2"


# 5. Modify all numbers that are bigger than 9999
ls -d $OUTDIR/*.coverage | while read Line
do
    echo $Line
    temp=$Line.temp
    /bin/sed 's/[0-9]*.[0-9]*e+0[1-9]$/99999/' $Line > $temp && mv $temp $Line
    echo " "
    echo " "
done

# 6. Merge m1.pos with m2.neg and m1.neg with m2.pos
#    PROTOCOL=1 --> ++,--; PROTOCOL=2 --> +-,-+
if [ $PROTOCOL != 2 ]; then
    cat $OUTDIR/m1.pos.coverage $OUTDIR/m2.neg.coverage > $OUT.pos.coverage
    cat $OUTDIR/m1.neg.coverage $OUTDIR/m2.pos.coverage > $OUT.neg.coverage
else
    cat $OUTDIR/m1.neg.coverage $OUTDIR/m2.pos.coverage > $OUT.pos.coverage
    cat $OUTDIR/m1.pos.coverage $OUTDIR/m2.neg.coverage > $OUT.neg.coverage
fi

# 7. Clean the directory

#/bin/rm -rf m1.pos.coverage m2.pos.coverage m1.neg.coverage m2.neg.coverage m1.bam m2.bam

echo "I am done with creating per nucleotide coverage file"

