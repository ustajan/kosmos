#!/usr/bin/env bash

# Usage: ./analysis.sh threshold input_file
# Example: ./analysis.sh 0.6 test.dat


# --------- 1. Check input arguments ----------
if [[ $# -ne 2 ]]; then
  echo "Usage: $0 <threshold> <input_file>"
  echo "Example: $0 0.6 test.dat"
  exit 1
fi
# --------- 1. Input arguments ----------
THRESHOLD="$1"
INPUT_FILE="$2"

awk -v threshold="$THRESHOLD" 'BEGIN{e=0;n=0;}{if($2>threshold) {e+=$2;n++;}}END{print "counts=" n " period=" 1e+9/(31591000*n/NR)"ns" " Average E=" e/n "MeV";}' $INPUT_FILE