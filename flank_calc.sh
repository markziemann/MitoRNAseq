#!/bin/bash

grep -w gene Rattus_norvegicus.mRatBN7.2.106.gtf \
| cut -f1,4,5,9 \
| cut -d '"' -f-2,6 \
| sed 's/gene_id "//' \
| sed 's/"/_/' \
| bedtools slop -i - -g chrNameLength.txt -b 3000 > Rattus_norvegicus.mRatBN7.2.106.genes.bed

cut -d '_' -f1 Rattus_norvegicus.mRatBN7.2.106.genes.bed \
| awk '{OFS="\t"} {print $4,$1,$2,$3,$3-$2}' > Rattus_norvegicus.mRatBN7.2.106.genes.saf
