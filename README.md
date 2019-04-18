The above scripts are a select subset of scripts which I have developed over the past ~1.5 years. 

**admixture.methods.py** this is a combined python pipeline which can be run on UGA Sapelo2 cluster and can be called like so:
`admixture.methods.py -i Random.list.txt` The -i is just a standard text document with one strain per row. 

**app.R** this is the shiny script which renders the chromopainting plot for a selected strain shown here: 
https://jhamlin.shinyapps.io/chromopaintApp/

**assemble.sh** this is a shell script to perform genome mapping for as many individuals as are in the directory. Essentially,
it makes subfolders for each pair of fastq output per individual and then runs fastqc, trimmomatic, fastqc on trimmed reads,
bwa-mem, and samtools ultimately generating a vcf file for each individual. 

**build.phylogeny.R** this is an R script to plot and color a phylogenetic tree produced by RAxML

**inter.intra.py** this is a python script to add a column specifying if the comparison is inter-population or intra-population
via a dictionary

**submit.admixture.sh** this is just a shell script to submit and run the admixture.methods.py script. Long run time 
is requested here as the admixture.method.py script runs until it has processed all of the files in the Random.list.txt
