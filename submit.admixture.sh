#!/bin/bash

#PBS -q bensasson_q
#PBS -N paint
#PBS -l walltime=48:00:00,vmem=10gb,nodes=1:ppn=1
#PBS -m abe
#PBS -M jhamlin@uga.edu

cd ~/jenna/scripts

admixture.methods.py -i Random.list.txt

