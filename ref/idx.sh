#!/bin/bash
#need to have the ensembl GTF and FA file in current dir
STAR=/mnt/mziemann/adam_trewin/sw/STAR-2.7.3a/bin/Linux_x86_64/STAR

GTF=$(find . | grep gtf$)
FA=$(find . | egrep '(.fa)')
CWD=$(pwd)
GNAMES=$(echo $GTF | sed 's#.gtf#.gnames.txt#')

$STAR --runMode genomeGenerate \
--sjdbGTFfile $GTF \
--genomeDir $CWD  \
--genomeFastaFiles $CWD/$FA \
--runThreadN $(nproc)

#grep -w gene $GTF | cut -d '"' -f2,6 \
#| tr '"' '\t' | sort -k 1b,1 > $GNAMES

grep -w gene $GTF  | cut -d '"' -f2,6,10 \
| tr '"' '\t' > $GNAMES
