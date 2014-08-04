#!/bin/bash
cd /spacebase/data/alexji/karlsson-model
envname=atomiccoolinghalo
sfrname=fixTS500
for sfrname in fixTS500 10xTS500
do
for p in "0.00" "0.10" "0.30" "0.50" "0.70" "0.90" "1.00"
do
    python plot_histlist.py --save $envname $sfrname "p${p}"
done
done