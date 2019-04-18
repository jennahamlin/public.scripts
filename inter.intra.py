#in python/3.4.3
import os
import pandas as pd 

#check working directory
cwd = os.getcwd()

#my tsv files are combined per population (i.e. cat Y12.pDiffs.tsv YJM1388.pDiffs.tsv YJM1389.pDiffs.tsv > Sake.pDiffs.tsv) 
#and then run the below for each population adjusting the intra to be population specific as an updated dictionary

# read in each individual TSV file
tsv_file = 'wine.combined.tsv'

# use pd.read_csv to load the tsv_file
datafile = pd.read_csv(tsv_file,sep="\t")


#read in the dictionary which defines each strain as being inter or intra this has to be adjusted for each population
d = {'HUN9.1s1': 'inter', 'PYR4b.1.1': 'inter', 'ZP560': 'inter', 'ZP561': 'inter', 'ZP562': 'inter', 'ZP565': 'inter','ZP568': 'inter', 'ZP633': 'inter', 'ZP636': 'inter',

'UWOPS03.461.4': 'inter', 
'SDO1s': 'inter', 'SDO2s1': 'inter', 'SDO3s1': 'inter', 'SDO4s1': 'inter', 'SDO6s1': 'inter', 'SDO7s1': 'inter', 'SM12s1': 'inter', 'SM17s1': 'inter', 'SM66s1': 'inter', 'SM69s1': 'inter', 
 
'YPS128m':'inter', 'YPS139s':'inter', 'YPS396': 'inter', 'YPS400':'inter', 'YPS600': 'inter', 'YPS602': 'inter', 'YPS604':'inter', 'YPS606m': 'inter', 'YPS608':'inter',  'YPS610':'inter', 

'Y12':'inter', 
'DBVPG6044':'inter', 

'AN3e.1.1':'intra','DBVPG1106':'intra', 'DBVPG1373':'intra', 'DBVPG1788':'intra', 'DBVPG6765':'intra',  'L1374':'intra', 'L1528':'intra', 'YJM975':'intra', 'YJM978':'intra', 'YJM981':'intra', 'ZP577':'intra', 'ZP578':'intra', 'ZP579':'intra',


'YJM1463': 'inter' }


d = {       'YPS128':'inter', 
             'YJM1273':'inter',
             'YPS606':'inter',
             'UWOPS03.461.4':'inter',
             'UWOPS05.217.3':'inter',
             'UWOPS05.227.2':'inter',
             'Y12':'inter',
             'K11':'inter',
             'YJM1388':'inter',
             'DBVPG6044':'inter',
             'YJM1439':'inter',
             'YJM1248':'inter',
             'L1374':'inter',
             'BC187':'inter',
             'YJM1332':'inter',
             'EN14S01':'inter',
             'EM14S01.3B':'inter',
             'GE14S01_7B':'inter',
             'ZP560':'inter',
             'ZP565':'inter',
             'ZP636':'inter',
             'SDO2s1':'inter',
             'SDO3s1':'inter',
             'SM66s1':'inter',
             'YJM1463':'inter',
             'ZP653':'inter',
             'ZP654':'inter',
             'ZP657':'inter',
             'ZP778':'intra',
             'ZP779':'intra',
             'ZP785':'intra',
             'ZP784':'intra',
             'ZP786':'intra',
             'ZP823':'intra',
             'ZP793':'intra',
             'UFMG.CM.Y264':'inter',
             'UFMG.CM.Y266':'inter',
             'UFMG.CM.Y455':'inter',
             'UFMG.CM.Y639':'inter',
             'UFMG.CM.Y641':'inter',
             'UFMG.CM.Y642':'inter',
             'UFMG.CM.Y263':'inter',
             'UFMG.CM.Y269':'inter',
             'S8BM.32.4D.a':'inter',
             'CBS1594':'inter',
             'CBS3000':'inter',
             'CBS1419':'inter',
             'CBS4456':'inter',
             'CLIB219.2b':'inter',
             'CBS7959':'inter',
             'CBS2888.1b':'inter',
             'CBS1463':'inter',
             'Y10.1B':'inter',
             'YPS1000':'inter',
             'ZP780':'intra',
             'UWOPS83.787.3':'inter',
             'YJM1400':'inter',
             'YJM1479':'inter',
             'YJM1401':'inter'
     }





#on data file look at the colun strain and use the dictionary to map either inter or intra
datafile['intra/inter'] = datafile['strain'].map(d)

#output the data file 
datafile.to_csv('wine.combined.out.csv')


###after outputing the new csv file, I then select columns 2 - 11 using awk and specifing the delimiter as a comma 
awk  'BEGIN{FS=OFS=","}{print $2,$3,$4,$5,$6,$7,$8,$9}' wine.combined.out.csv > wine.combined.out.final.csv 

##remove duplicate rows 
awk '!seen[$0]++' combined.test.out.pDiffs2.csv > combined.test.out.pDiffs3.csv


