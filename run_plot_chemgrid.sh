#!/bin/bash

cd /spacebase/data/alexji/karlsson-model
envname=atomiccoolinghalo
sfrname=burstachsimTS500
#sfrname=fidTS500
#sfrname=flatfixTS500
for p in "0.50" #"0.00" "0.10" "0.30" "0.50" "0.70" "0.90" "1.00"
do
    python plot_chemgrid.py --save $envname $sfrname "p${p}"
    #python plot_cfrac.py --save $envname $sfrname "p${p}"
done

