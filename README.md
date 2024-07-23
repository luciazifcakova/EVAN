# EVAN
EVAN tool for imaging GC heterogeneity in genomes. Written by Marta Vohnoutova and Lucia Zifcakova.

![EvanLogo.gif](https://github.com/martavohnoutova/Evan_Octopus/blob/main/Evan%20logo2.png)

Python based tool can create Guanine Cytosine percentage distribution plots across soft-masked genomes in chromosomal resolutions. Prior runing EVAN, please, soft-mask your genomes using RepeatMasker (see more info in https://github.com/luciazifcakova/scripts_for_paper_Regional_distribution_GC_content). The output files are figures created separately for each chromosome, showing genomic coordinates on X axis and percentage of GC% on Y axis as a function of masked(repetitive)/unmasked(coding) proportions of genomic region. Additional output files show the plot of average GC% of each chromosome as a function of masked/unmasked proportions. Output csv files contains data about rep% (% of repetitive sequences), gc% (GC% of masked/repetitive sequences) and GC% (% of unmasked(coding) sequences). The tool consists of 4 configurable Jupyter programs and can operate on many samples together in cycle.

You can run Jypyter notebook scripts as described in individuall script files (e.g.: 1_GC_repeats_profile_PLOS_TRAINING_universal_v.0.2.ipynb etc.)

Otehrwise, there are 3 options how to run EVAN tool.

1. Download GC_with_args_for_HPC.py and run it with python3 like this:

"python3 /path/to/your/script/GC_with_args_for_HPC.py /path/to/your/fasta.fa species_name window_size"

E.G.: "python3 /Users/GC_with_args_for_HPC.py /Users/fasta.fa Nautilus 3000"

2. Run with Docker. Make sure you have Docker installed and running on you local machine.
2.a Either build your own Docker container from Dockerfile from this repository
or
2.b pull from Docker hub

"docker pull zifcakova/evan_hpc_ver2.2"

run Docker container with external data by mounting volume to container:

"docker run -p 4000:8080 -v /Users:/data evan_hpc_ver2.2 /data/fasta.fa Nautilus 3000"

3. Run with Singularity. Have Singularity installed on your machine.
Pull the image:

"singularity pull evan_hpc_ver2.2.sif docker://zifcakova/evan_hpc_ver2.2:latest"

Execute in the writtable location, as the script writes output into input directory. Have your fasta files in /writable/location/outside/of/container. 

"singularity exec --bind /writable/location/outside/of/container:/mnt \
/path/to/your/singularity/image/evan_hpc_ver2.2.sif \
python /app/GC_percentage_masked_unmasked_in_genome.py /mnt/fasta.fa Nautilus 3000"

Aternatevely, for high throughput, you can use slurm array to feed input into singularity container in this way (example is for 3 input fasta files):

#!/bin/bash

#SBATCH --job-name=singGC

#SBATCH --partition=your_partition_name

#SBATCH --time=1-0

#SBATCH --mem-per-cpu=10G

#SBATCH --cpus-per-task=20

#SBATCH --mail-user=your_email_registered_at_the_cluster

#SBATCH --mail-type=BEGIN,FAIL,END

#SBATCH --array=0-2 #3 input fasta files

#SBATCH --output=./singGC_%A-%a.out

files=(/writable/location/outside/of/container/*_soft_masked.fa)

file=${files[${SLURM_ARRAY_TASK_ID}]}

file1=$(basename $file)

out=$(basename -s .fa $file) 

ml bioinfo-ugrp-modules #load necessary module for your cluster

ml DebianMed

ml singularity

singularity exec --bind /writable/location/outside/of/container:/mnt \
/path/to/your/singularity/image/evan_hpc_ver2.2.sif \
python /app/GC_percentage_masked_unmasked_in_genome.py /mnt/${file1} ${out} 3000 #3000 is window size, you can change to any size you want (over 1000bp makes sense)
