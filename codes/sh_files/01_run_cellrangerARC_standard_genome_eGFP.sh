#!/bin/bash

# multiome specific reference
GRCh38="path for refdata-cellranger-arc-GRCh38-2020-A-2.0.0"
#GRCz11=/lustre/home/mc900/RP_T121869/data1/references/danRer11/GRCz11
GRCz11="path for danRer11_UCSC_eGFP_crh_only"
# for bioseqfs02:
CRARC="path for cellranger-arc" 

while read sample; do 
  $CRARC count \
     --id ${sample} \
     --reference=${GRCz11} \
     --libraries=${sample}_libraries.csv
done < samples.txt

