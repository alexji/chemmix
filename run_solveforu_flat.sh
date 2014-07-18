#!/bin/csh

cd ~/karlsson-model
set envname = atomiccoolinghalo
foreach ts (500 600 700 800 900 1000)
python solveforu.py $envname flatfixTS${ts}
python solveforu.py $envname flat10xTS${ts}
end
