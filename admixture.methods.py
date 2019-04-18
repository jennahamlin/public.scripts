#!/usr/bin/env python

import argparse
import os
import sys
import re
import subprocess
import random
import time

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--inputFile', help = "Name of inputFile which is a txt file with each strain listed on one line")
args = parser.parse_args()

inputFile = args.inputFile
inputFile = list(open(inputFile))

#static location of whole genome files to make a list of strains to use
faDir = "/scratch/jhamlin/workDir/Admixture.manuscript/Backbone2/Confirmed.Whole.Genomes"

#need to define a variable for -M 
#M = 0.00145
#relaxedM = 0.00367

##### FUNCTION TO RANDOMLY SELECT FROM AN INPUT TEXT FILE #####
def makeRandom(inputFile):
	rand_item=random.sample(inputFile, 1)
	rand_item[:]=[item.rstrip('\n') for item in rand_item]
	return rand_item

##### FUNCTION TO GET DATA #####
def getData():
	cmd="""cp /scratch/jhamlin/workDir/Admixture.manuscript/Backbone2/Confirmed.Chr.Data/%s*.fa /scratch/jhamlin/workDir/Admixture.manuscript/Backbone2/Random.Test """ %(strain)
	subprocess.call(cmd, shell=True)
#	cmd=""" mv /scratch/jhamlin/workDir/Admixture.manuscript/Backbone2/Confirmed.Chr.Data/%s.fa /scratch/jhamlin/workDir/Admixture.manuscript/Backbone2/Confirmed.Whole.Genomes""" %(strain)
#	subprocess.call(cmd, shell=True)


##### FUNCTION TO EXTRACT LIST OF STRAIN NAMES FROM DIRECTORY #####
def strainList(faDir):  #from Rosemary
        strainlist = []
        for filename in os.listdir(faDir):
                try:
                        match = re.match('(.*).fa', filename)
                        strain = match.group(1)
                        strainlist.append(strain)
                except AttributeError:
                        pass
        return strainlist

##### FUNCTION TO COMPILE SPLIT CHROMOSOME FILES AND PAINT CHROMSOMES#####
def CompilePaint (strain):
	#strainlist = strainList(faDir)
	#strainlist= (' '.join(strainlist))
	command = """#!/bin/bash

#PBS -q bensasson_q
#PBS -N compile.paint
#PBS -l walltime=10:00:00,vmem=8gb,nodes=1:ppn=1
#PBS -m abe
#PBS -M jhamlin@uga.edu

module load Perl
module load R
 
cd /scratch/jhamlin/workDir/Admixture.manuscript/Backbone2/Random.Test

for chr in chrI chrII chrIII chrIV chrV chrVI chrVII chrVIII chrIX chrX chrXI chrXII chrXIII chrXIV chrXV chrXVI chrM
do
	cat *_$chr.fa >$chr.mfa
	perl /scratch/jhamlin/workDir/scripts/fastatools/fastaLC2n.pl -i $chr.mfa -o $chr.temp
	perl /scratch/jhamlin/workDir/scripts/fastatools/alcat.pl -i $chr.temp -o $chr.mfa
done

rm chr*.temp

mkdir /scratch/jhamlin/workDir/Admixture.manuscript/Backbone2/Random.Test/%s

mv chrI.mfa /scratch/jhamlin/workDir/Admixture.manuscript/Backbone2/Random.Test/%s/chr01.mfa
mv chrII.mfa /scratch/jhamlin/workDir/Admixture.manuscript/Backbone2/Random.Test/%s/chr02.mfa
mv chrIII.mfa /scratch/jhamlin/workDir/Admixture.manuscript/Backbone2/Random.Test/%s/chr03.mfa
mv chrIV.mfa /scratch/jhamlin/workDir/Admixture.manuscript/Backbone2/Random.Test/%s/chr04.mfa
mv chrV.mfa /scratch/jhamlin/workDir/Admixture.manuscript/Backbone2/Random.Test/%s/chr05.mfa
mv chrVI.mfa /scratch/jhamlin/workDir/Admixture.manuscript/Backbone2/Random.Test/%s/chr06.mfa
mv chrVII.mfa /scratch/jhamlin/workDir/Admixture.manuscript/Backbone2/Random.Test/%s/chr07.mfa
mv chrVIII.mfa /scratch/jhamlin/workDir/Admixture.manuscript/Backbone2/Random.Test/%s/chr08.mfa
mv chrIX.mfa /scratch/jhamlin/workDir/Admixture.manuscript/Backbone2/Random.Test/%s/chr09.mfa
mv chrX.mfa /scratch/jhamlin/workDir/Admixture.manuscript/Backbone2/Random.Test/%s/chr10.mfa
mv chrXI.mfa /scratch/jhamlin/workDir/Admixture.manuscript/Backbone2/Random.Test/%s/chr11.mfa
mv chrXII.mfa /scratch/jhamlin/workDir/Admixture.manuscript/Backbone2/Random.Test/%s/chr12.mfa
mv chrXIII.mfa /scratch/jhamlin/workDir/Admixture.manuscript/Backbone2/Random.Test/%s/chr13.mfa
mv chrXIV.mfa /scratch/jhamlin/workDir/Admixture.manuscript/Backbone2/Random.Test/%s/chr14.mfa
mv chrXV.mfa /scratch/jhamlin/workDir/Admixture.manuscript/Backbone2/Random.Test/%s/chr15.mfa
mv chrXVI.mfa /scratch/jhamlin/workDir/Admixture.manuscript/Backbone2/Random.Test/%s/chr16.mfa

cd /scratch/jhamlin/workDir/Admixture.manuscript/Backbone2/Random.Test/%s

AR=(%s)

for i in "${AR[@]}"
do 
perl /scratch/jhamlin/workDir/scripts/fastatools/faChrompaint.pl -I chr -r ${i} -W 50000 -M 0.00145 -c /scratch/jhamlin/workDir/Admixture.manuscript/Backbone2/clade.file.txt -C /scratch/jhamlin/workDir/Admixture.manuscript/Backbone2/color.file.txt -e %s
done 

cat >> /scratch/jhamlin/workDir/Admixture.manuscript/Backbone2/Random.Test/Random.All.Count.tsv << EOF
%s
EOF

cat %s.nearest.tsv | awk '{ print $5 }' | sort | uniq -c | head -n -1 >> /scratch/jhamlin/workDir/Admixture.manuscript/Backbone2/Random.Test/Random.All.Count.tsv 

""" % (strain, strain, strain, strain, strain, strain, strain, strain, strain, strain, strain, strain, strain, strain, strain, strain, strain, strain, strain, strain, strain, strain)	#only strain needed here to replace %s. If more output files placeholderes than more calls to outputfile here  #Can change the AR section to strainlist to do all files in the location of whole genomes folder 
	outF = open("Compile.Paint.sh", "w")
	for line in command:
		outF.write(line)
	outF.close()

##### FUNCTION TO COUNT NUMBER OF WINDOWS AND KEEP STRAIN IF DIVERGED OR REPRESENTIVE OF A POPULATION #####


##### FUNCTION TO RUN RAXML AND BUILD PHYLOGENY #####
def BuildPhylogeny (strain): 
	command = """#!/bin/bash

#PBS -q bensasson_q
#PBS -N RAxML
#PBS -l walltime=70:00:00,vmem=15gb,nodes=1:ppn=16
#PBS -m abe
#PBS -M jhamlin@uga.edu

module load Perl
module load RAxML/8.2.4-foss-2016b-pthreads-avx 

cd /scratch/jhamlin/workDir/Admixture.manuscript/Backbone2/Random.Test/%s

grep '>' chr01.mfa > strain.confirmation.tsv

perl /scratch/jhamlin/workDir/scripts/fastatools/alcat.pl -i 'chr01.mfa chr02.mfa chr03.mfa chr04.mfa chr05.mfa chr06.mfa chr07.mfa chr08.mfa chr09.mfa chr10.mfa chr11.mfa chr12.mfa chr13.mfa chr14.mfa chr15.mfa chr16.mfa'  -o %s.mfa

raxmlHPC-PTHREADS-AVX -T 16 -f a -x 12345 -p 12345 -m GTRGAMMA -N 100 -s %s.mfa -n %s.mfa

""" % (strain, strain, strain, strain) 
	outF = open("Build.Phylogeny.sh", "w")
	for line in command:
		outF.write(line)
	outF.close()

##### FUNCTION TO REMOVE DATA #####
def removeData():
        cmd="""rm /scratch/jhamlin/workDir/Admixture.manuscript/Backbone2/Random.Test/%s*.fa """ %(strain)
        subprocess.call(cmd, shell=True)

	#cmd="""rm /scratch/jhamlin/workDir/Admixture.manuscript/Backbone2/Confirmed.Whole.Genomes/%s*.fa """ %(strain)
        #subprocess.call(cmd, shell=True)

oldstrains = []

while len(oldstrains)!=len(inputFile):
	strain = makeRandom(inputFile)[0]
	if strain not in oldstrains:
		getData()
		CompilePaint(strain)
		BuildPhylogeny(strain)
		qsub = os.popen("""qsub Compile.Paint.sh""").read()
		qsub1= qsub.rstrip()[0:6]
		os.popen("""$(qsub -W depend=afterok:%s Build.Phylogeny.sh)"""%(qsub1))
		time.sleep(180) #holds the python script for X seconds allowing the first script to run and then performs the remove data funtion and then continuing on with the loop. Not sure if this is the best but is a solution
#		removeData()
		oldstrains.append(strain)
	else:
		continue
