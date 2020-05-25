#!/bin/bash

VCF_FILE=$1
MIDPOINT=$2
WIDTH=$3

# TODO : put in a usage message here
[ $# -eq 0 ] && { echo "Usage: $0 <vcf> <midpoint> <width>"; exit 1;}

# Running the plink inference for the r2 matrix 
plink --vcf $VCF_FILE --double-id --chr 22  --from-bp $(($MIDPOINT - $WIDTH/2)) --to-bp $(($MIDPOINT + $WIDTH/2)) --r2 triangle --out test

# Cleaning up the files because we don't acutally want to keep them around 
rm test.*
rm plink.log
