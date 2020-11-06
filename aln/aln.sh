#!/bin/bash
set -x
REF=../ref

STAR=../sw/STAR-2.7.3a/bin/Linux_x86_64/STAR
SKEWER=../sw/skewer/skewer
CWD=$(pwd)


for FQZ1 in *R1*.fastq.gz ; do

FQZ2=$(echo $FQZ1 | sed 's/_1/_2/')
FQ1=$(echo $FQZ1 | sed 's/.gz$/-trimmed-pair1.fastq/')
FQ2=$(echo $FQZ1 | sed 's/.gz$/-trimmed-pair2.fastq/')
BASE=$(echo $FQZ1 | sed 's/.fq.gz//')
BAM=$BASE.bam

# quality trim with skewer
$SKEWER -t $(nproc) -q 20 \
 -x AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -y AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
$FQZ1 $FQZ2

# clip off umi
#fastx_trimmer -t 10 -i $FQ1 -o tmp && mv tmp $FQ1
#fastx_trimmer -f 13 -i $FQ2 -o tmp && mv tmp $FQ2

# map with STAR single end mode
$STAR --runThreadN 30 --quantMode GeneCounts --genomeLoad LoadAndKeep  \
  --outSAMtype None --genomeDir $REF --readFilesIn=$FQ1 \
  --outFileNamePrefix $BASE.

rm $FQ1 $FQ2

done

$STAR --genomeLoad Remove --genomeDir $REF

for TAB in *ReadsPerGene.out.tab ; do
  tail -n +5 $TAB | cut -f1,4 | sed "s/^/${TAB}\t/"
done > 3col.tsv


