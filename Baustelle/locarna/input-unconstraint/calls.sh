#!/usr/bin/env bash

# get data
for f in ../../cofold/input/locARNA_*.fa; do grep -P "^[>ACGUT]|#\d" $f > $(basename $f); done


for f in locARNA_*.fa; do 
  awk '{if (substr($1,1,1)==">"){print $1}else{print substr($1,1,110)" "$2}}' $f > 5p.$f
  awk '{if (substr($1,1,1)==">"){print $1}else{print substr($1,118)" "$2}}' $f > 3p.$f
done


