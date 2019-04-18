#!/bin/bash

#PBS -q bensasson_q
#PBS -N assemble
#PBS -l walltime=1:00:00,vmem=64gb,nodes=1:ppn=1
#PBS -m abe
#PBS -M jhamlin@uga.edu

cd /scratch/jhamlin/workDir

#Insert in each strains individual folder name 
for file in *_1.fastq
do
        strain_name=$(echo $file | cut -d'_' -f 1)
        ar+=($strain_name)
done

for strain in "${ar[@]}"
do
 #I have both paired end reads for each strain parsed into individual folders
 #so change into each folder unzip each read (if necessary) and then make the two scripts (X_fastqc.raw.trim.sh and X_assembly.sh)
mkdir ${strain} 
mv ${strain}_1.fastq ${strain}_2.fastq /scratch/jhamlin/workDir/${strain} 
cd ${strain}
####fastqc on raw fastq, clean and trim and fastqc post trimming #### uses fastqc and trimmomatic
 echo "#!/bin/bash">${strain}_fastqc.raw.trim.sh
 echo " ">>${strain}_fastqc.raw.trim.sh
 echo "#PBS -q bensasson_q" >>${strain}_fastqc.raw.trim.sh
 echo "#PBS -l walltime=5:00:00,vmem=2gb:nodes=1:ppn=1" >>${strain}_fastqc.raw.trim.sh
 echo "#PBS -M jhamlin@uga.edu" >>${strain}_fastqc.raw.trim.sh
 echo "#PBS -m abe" >>${strain}_fastqc.raw.trim.sh
 echo "#PBS -j oe" >>${strain}_fastqc.raw.trim.sh
 echo " ">> ${strain}_fastqc.raw.trim.sh
 echo "cd /scratch/jhamlin/workDir/${strain}">>${strain}_fastqc.raw.trim.sh
 echo " " >>${strain}_fastqc.raw.trim.sh
 echo "module load Java/1.8.0_144">>${strain}_fastqc.raw.trim.sh
 echo "module load FastQC/0.11.8-Java-1.8.0_144 ">>${strain}_fastqc.raw.trim.sh
 echo " " >>${strain}_fastqc.raw.trim.sh
 echo "mkdir fastqc.raw">>${strain}_fastqc.raw.trim.sh
 echo "mkdir fastqc.trimmed">>${strain}_fastqc.raw.trim.sh
 echo " " >>${strain}_fastqc.raw.trim.sh
 echo "fastqc -o fastqc.raw/ ${strain}_1.fastq">>${strain}_fastqc.raw.trim.sh
 echo "fastqc -o fastqc.raw/ ${strain}_2.fastq">>${strain}_fastqc.raw.trim.sh
 echo " " >>${strain}_fastqc.raw.trim.sh
 echo "module load Trimmomatic/0.33-Java-1.8.0_144">>${strain}_fastqc.raw.trim.sh
 echo "java -jar \$EBROOTTRIMMOMATIC/trimmomatic-0.33.jar PE -baseout ${strain}.fq /scratch/jhamlin/workDir/${strain}/${strain}_1.fastq /scratch/jhamlin/workDir/${strain}/${strain}_2.fastq SLIDINGWINDOW:4:20 MINLEN:36">>${strain}_fastqc.raw.trim.sh
 echo " " >>${strain}_fastqc.raw.trim.sh
 echo "fastqc -o fastqc.trimmed/ ${strain}_1P.fq">>${strain}_fastqc.raw.trim.sh
 echo "fastqc -o fastqc.trimmed/ ${strain}_2P.fq">>${strain}_fastqc.raw.trim.sh

####Mapping#### uses bwa mem this assumes you already have the reference genome sorted and indexed 
#and samtools for  mpileup to generate a vcf file individually
 echo "#!/bin/bash">${strain}_assembly.sh
 echo " ">>${strain}_assembly.sh
 echo "#PBS -q bensasson_q" >>${strain}_assembly.sh
 echo "#PBS -l walltime=10:00:00,vmem=8gb:nodes=1:ppn=1" >>${strain}_assembly.sh
 echo "#PBS -M jhamlin@uga.edu" >>${strain}_assembly.sh
 echo "#PBS -m abe" >>${strain}_assembly.sh
 echo "#PBS -j oe" >>${strain}_assembly.sh
 echo " ">> ${strain}_assembly.sh
 echo "cd /scratch/jhamlin/workDir/${strain}">>${strain}_assembly.sh
 echo "OD=/scratch/jhamlin/workDir/${strain}/vcf">>${strain}_assembly.sh
 echo " " >>${strain}_assembly.sh
 echo "mkdir vcf">>${strain}_assembly.sh
 echo " " >>${strain}_assembly.sh
 echo "module load BWA">>${strain}_assembly.sh
 echo "module load SAMtools/1.3.1-foss-2016b">>${strain}_assembly.sh
 echo "module load BCFtools/1.3.1-foss-2016b">>${strain}_assembly.sh
 echo " " >>${strain}_assembly.sh
 echo "bwa mem /scratch/jhamlin/workDir/sacCer3/sacCer3.fa ${strain}_1P.fq  ${strain}_2P.fq | samtools sort -o ${strain}.sorted -">>${strain}_assembly.sh
 echo " " >>${strain}_assembly.sh
 echo "samtools mpileup -d 100000 -uf /lustre1/jhamlin/workDir/sacCer3/sacCer3.fa ${strain}.sorted -I | bcftools call -c -f GQ > \$OD/${strain}.raw.vcf">>${strain}_assembly.sh
 
 chmod 755 ${strain}_fastqc.raw.trim.sh
 chmod 755 ${strain}_assembly.sh

#make a variable called qcheck which is the submission id generated after qsub ${strain}_fastqc.raw.trim.sh 
 qcheck=$(qsub ${strain}_fastqc.raw.trim.sh)
#strip the .pbs.scm from the job id and store in qcheck 2
 qcheck2=${qcheck%%.pbs.scm}
#qsub the next script (assembly.sh) after the first script (fastqc.raw.trim.sh) is completed
 qsub -W depend=afterok:$qcheck2 ${strain}\_assembly.sh

  #cd back up and move into the next strains folder and repeat
 cd ..

done
