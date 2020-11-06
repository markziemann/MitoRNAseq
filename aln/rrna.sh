#!/bin/bash
set -x
REF=../ref

for FQZ1 in *R1*.fastq.gz ; do
  BASE=$(echo $FQZ1 | sed 's/.fq.gz//')
  FQZ2=$(echo $FQZ1 | sed 's/_1/_2/')
  BASE=$(echo $FQZ1 | sed 's/.fq.gz//')
  gunzip -dc $FQZ1 | head -4000000 \
  | bwa mem -t 30 ../ref/rrna.fa - \
  | samtools view -@16 -bS - > $BASE.rrna.bam
done

for BAM in *bam ; do samtools flagstat $BAM > $BAM.flagstat.txt ; done
