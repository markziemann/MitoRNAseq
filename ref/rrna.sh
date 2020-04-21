#!/bin/bash
for t in $(grep rRNA Rattus_norvegicus.Rnor_6.0.ncrna.fa | cut -d  ' ' -f1 | tr -d '>') ; do
  samtools faidx Rattus_norvegicus.Rnor_6.0.ncrna.fa $t
done > rrna.fa

