#/bin/bash

# purpose is to quantify exonic (true) reads and intronic/flanking (noise) reads.
# we will then use this to filter true hits from noise.

# exonic reads summarised to genes

GTF=../ref2/Rattus_norvegicus.mRatBN7.2.106.gtf
SAF=../ref2/Rattus_norvegicus.mRatBN7.2.106.genes.saf

for BAM in *bam ; do

  OUT_EXON=$BAM.exon.tsv
  OUT_GENE=$BAM.gene.tsv
  OUT=$BAM.fc.tsv

  featureCounts -a $GTF \
    -F GTF -t exon -Q 20 -s 0 -T 16 -O \
    -o $OUT_EXON $BAM

  featureCounts -a $SAF \
    -F SAF -f -Q 20 -s 0 -T 16 -O \
    -o $OUT_GENE $BAM

  sed 1d $OUT_EXON | cut -f1,7 > tmp ; mv tmp $OUT_EXON
  sed 1d $OUT_GENE | cut -f1,7 > tmp ; mv tmp $OUT_GENE

  join -1 1 -2 1 $OUT_EXON $OUT_GENE | sed 's/ /\t/g' > $OUT
  sed -i "s/^/${BAM}\t/" $OUT_EXON
  sed -i "s/^/${BAM}\t/" $OUT_GENE

done


cat *exon.tsv | grep -v "Geneid" > EXON.tsv
cat *gene.tsv | grep -v "Geneid" > GENE.tsv

