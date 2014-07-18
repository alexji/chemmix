#!/bin/bash

cd ~/karlsson-model
envname=atomiccoolinghalo
for ts in 300 400 500  #600 700 800 900 1000
do
    sfrname="fidTS${ts}"
    python plot_fMkkp.py $envname $sfrname
done

