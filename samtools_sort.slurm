#!/usr/bin/env bash
#SBATCH --partition=short    ### Partition
#SBATCH --job-name=sam_view.job  ### Job Name
#SBATCH --output=sam_view.job.out    ### File in which to store job output
#SBATCH --time=0-23:00:00   ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1           ### Number of nodes needed for the job
#SBATCH --ntasks-per-node=28 ### Number of tasks to be launcged per Node

module purge
module load samtools/1.5

cd /projects/bgmp/mdriskil/deduper

files=`ls -1 ./*.sam`
for file in ${files}; do
	 trimmedfilename=`echo ${file} | cut -d "/" -f2 | cut -d "." -f1`;
samtools sort -o sorted_${trimmedfilename}.sam -O sam ${file};

done
exit
