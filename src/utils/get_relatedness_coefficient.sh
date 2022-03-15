#!/bin/bash

sed 's/ /\t/g' meansP0_AncientDNA_normalized | join -j1 -t$'\t' - READ_results | awk 'BEGIN{FS="\t"; OFS=","}{print $1, 2*(1-$2), $6}' | column -s, -t
