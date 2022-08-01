#!/bin/bash
set -e 

testf='test.fmt6'
oud=test_oud_dir
mkdir -p $oud

python3 -m AIMH \
    -i $testf -o $oud/fmt6.tsv -O $oud/AIMH.tsv \
    -s ',' -S $'\t' -l bn





