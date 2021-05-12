#!/bin/bash

zcat Rattus_norvegicus.Rnor_6.0.cdna+ncrna.fa.gz \
| grep '>' | sed 's/:/\n/2' | sed 's/:/ /' \
| sed 's/gene_symbol:/\n/' | sed 's/gene_biotype:/\n/' \
| sed 's/gene:/\n/' \
| cut -d ' ' -f1  |  paste - - - - - \
| sed 's/>//' \
| awk '{OFS="\t"} {print $1,$2,$3,$5,$4}' > Rattus_norvegicus.Rnor_6.0.cdna+ncrna.gene_names.tsv
