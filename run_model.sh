#!/bin/bash
#python run_model.py -j 16 --mmixgrid atomiccoolinghalo dummysfrname
#python run_model.py -j 16 --sfr --ckk --fMkk atomiccoolinghalo fixTS500
#python run_model.py -j 16 --sfr --ckk --fMkk atomiccoolinghalo 10xTS500
#python run_model.py -j 16 --sfr --ckk --fMkk atomiccoolinghalo fixTS600
#python run_model.py -j 16 --sfr --ckk --fMkk atomiccoolinghalo 10xTS600

#python run_model.py -j 16 --mmixgrid minihalo dummysfrname
#python run_model.py -j 16 --sfr --ckk --fMkk minihalo fixTS150 
#python run_model.py -j 16 --sfr --ckk --fMkk minihalo fixTS200

python run_model.py -j 8 --mmixgrid minihalo_lowDt dummysfrname
python run_model.py -j 8 --sfr --ckk --fMkk minihalo_lowDt fixTS150 
python run_model.py -j 8 --sfr --ckk --fMkk minihalo_lowDt fixTS200

#python run_model.py -j 16 --mmixgrid atomiccoolinghalo_lowDt dummysfrname
#python run_model.py -j 16 --sfr --ckk --fMkk atomiccoolinghalo_lowDt 10xTS500
#python run_model.py -j 16 --sfr --ckk --fMkk atomiccoolinghalo_lowDt fixTS500
#python run_model.py -j 16 --sfr --ckk --fMkk atomiccoolinghalo_lowDt 10xTS600
#python run_model.py -j 16 --sfr --ckk --fMkk atomiccoolinghalo_lowDt fixTS600

#python run_model.py -j 16 --mmixgrid atomiccoolinghalo_lowmass dummysfrname
#python run_model.py -j 16 --sfr --ckk --fMkk atomiccoolinghalo_lowmass 10xTS500
#python run_model.py -j 16 --sfr --ckk --fMkk atomiccoolinghalo_lowmass fixTS500
#python run_model.py -j 16 --sfr --ckk --fMkk atomiccoolinghalo_lowmass 10xTS600
#python run_model.py -j 16 --sfr --ckk --fMkk atomiccoolinghalo_lowmass fixTS600

#python run_model.py -j 16 --mmixgrid minihalo dummysfrname
#python run_model.py -j 16 --mmixgrid atomiccoolinghalo dummysfrname

#python run_model.py --sfr minihalo fixTS150
#python run_model.py --sfr minihalo fixTS200
#python run_model.py --sfr atomiccoolinghalo fixTS500
#python run_model.py --sfr atomiccoolinghalo 10xTS500
#python run_model.py --sfr atomiccoolinghalo fixTS600
#python run_model.py --sfr atomiccoolinghalo 10xTS600
#python run_model.py --sfr atomiccoolinghalo fixTS700
#python run_model.py --sfr atomiccoolinghalo 10xTS700

#python run_model.py -j 16 --ckk --fMkk minihalo fixTS150
#python run_model.py -j 16 --ckk --fMkk minihalo fixTS200
#python run_model.py -j 16 --ckk --fMkk atomiccoolinghalo fixTS500
#python run_model.py -j 16 --ckk --fMkk atomiccoolinghalo 10xTS500
#python run_model.py -j 16 --ckk --fMkk atomiccoolinghalo fixTS600
#python run_model.py -j 16 --ckk --fMkk atomiccoolinghalo 10xTS600

