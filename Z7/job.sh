#!/bin/bash
#PBS -N i2D_C
#PBS -q cpu
#PBS -l select=1:ncpus=1
#PBS -l walltime=06:00:00
#PBS -o i2D.scr.log
#PBS -e i2D.err.log
#PBS -M mberos@pmfst.hr
#PBS -m e

cd ${PBS_O_WORKDIR}

gcc i2D.c -o i2D -lm

./i2D

