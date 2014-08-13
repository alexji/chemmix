#!/bin/bash
numprocs=8
#python run_model.py -j $numprocs --mmixgrid atomiccoolinghalo dummysfrname
#python run_model.py -j $numprocs --sfr --ckk --fMkk atomiccoolinghalo fixTS500
#python run_model.py -j $numprocs --sfr --ckk --fMkk atomiccoolinghalo 10xTS500
#python run_model.py -j $numprocs --sfr --ckk --fMkk atomiccoolinghalo fixTS600
#python run_model.py -j $numprocs --sfr --ckk --fMkk atomiccoolinghalo 10xTS600
#python run_model.py -j $numprocs --sfr --ckk --fMkk atomiccoolinghalo burstfixTS500
#python run_model.py -j $numprocs --sfr --ckk --fMkk atomiccoolinghalo burst10xTS500
#python run_model.py -j $numprocs --sfr --ckk --fMkk atomiccoolinghalo burstachsimTS500
#python run_model.py -j $numprocs --sfr --ckk --fMkk atomiccoolinghalo burstachsloTS500

#python run_model.py -j $numprocs --mmixgrid minihalo dummysfrname
#python run_model.py -j $numprocs --sfr --ckk --fMkk minihalo fixTS150 
#python run_model.py -j $numprocs --sfr --ckk --fMkk minihalo fixTS200
#python run_model.py -j $numprocs --sfr --ckk --fMkk minihalo mhsimTS150 

python run_model.py -j $numprocs --mmixgrid minihalo_lowE dummysfrname
#python run_model.py -j $numprocs --sfr --ckk --fMkk minihalo_lowE mhsimTS150 

#python run_model.py -j $numprocs --mmixgrid minihalo_lowDt dummysfrname
#python run_model.py -j $numprocs --sfr --ckk --fMkk minihalo_lowDt fixTS150 
#python run_model.py -j $numprocs --sfr --ckk --fMkk minihalo_lowDt fixTS200

#python run_model.py -j $numprocs --mmixgrid atomiccoolinghalo_lowDt dummysfrname
#python run_model.py -j $numprocs --sfr --ckk --fMkk atomiccoolinghalo_lowDt 10xTS500
#python run_model.py -j $numprocs --sfr --ckk --fMkk atomiccoolinghalo_lowDt fixTS500
#python run_model.py -j $numprocs --sfr --ckk --fMkk atomiccoolinghalo_lowDt 10xTS600
#python run_model.py -j $numprocs --sfr --ckk --fMkk atomiccoolinghalo_lowDt fixTS600

#python run_model.py -j $numprocs --mmixgrid atomiccoolinghalo_lowmass dummysfrname
#python run_model.py -j $numprocs --sfr --ckk --fMkk atomiccoolinghalo_lowmass 10xTS500
#python run_model.py -j $numprocs --sfr --ckk --fMkk atomiccoolinghalo_lowmass fixTS500
#python run_model.py -j $numprocs --sfr --ckk --fMkk atomiccoolinghalo_lowmass 10xTS600
#python run_model.py -j $numprocs --sfr --ckk --fMkk atomiccoolinghalo_lowmass fixTS600




#python run_model.py --sfr minihalo fixTS150
#python run_model.py --sfr minihalo fixTS200
#python run_model.py --sfr atomiccoolinghalo fixTS500
#python run_model.py --sfr atomiccoolinghalo 10xTS500
#python run_model.py --sfr atomiccoolinghalo fixTS600
#python run_model.py --sfr atomiccoolinghalo 10xTS600
#python run_model.py --sfr atomiccoolinghalo fixTS700
#python run_model.py --sfr atomiccoolinghalo 10xTS700

#python run_model.py -j $numprocs --ckk --fMkk minihalo fixTS150
#python run_model.py -j $numprocs --ckk --fMkk minihalo fixTS200
#python run_model.py -j $numprocs --ckk --fMkk atomiccoolinghalo fixTS500
#python run_model.py -j $numprocs --ckk --fMkk atomiccoolinghalo 10xTS500
#python run_model.py -j $numprocs --ckk --fMkk atomiccoolinghalo fixTS600
#python run_model.py -j $numprocs --ckk --fMkk atomiccoolinghalo 10xTS600

