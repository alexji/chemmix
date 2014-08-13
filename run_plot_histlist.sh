#!/bin/bash
cd /spacebase/data/alexji/karlsson-model
envname=atomiccoolinghalo
#envname=atomiccoolinghalo_lowDt
#envname=minihalo

for sfrname in fixTS500 10xTS500 
#for sfrname in burstfixTS500 burst10xTS500 
#for sfrname in fixTS150 #fixTS200
do
#python plot_histlist.py --save --plottwo $envname $sfrname
#python plot_histlist.py --save --plotthree $envname $sfrname
#python plot_histlist.py --save --plotfour $envname $sfrname

#for p in "0.00" "0.10" "0.30" "0.50" "0.70" "0.90" "1.00"
#do
#    python plot_histlist.py --save --plotone $envname $sfrname "p${p}"
#done
done