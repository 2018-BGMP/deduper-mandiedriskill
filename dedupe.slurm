#!/usr/bin/env bash
#SBATCH --partition=short    ### Partition
#SBATCH --job-name=dedupejob  ### Job Name
#SBATCH --output=dedupejob.out    ### File in which to store job output
#SBATCH --time=0-00:30:00   ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1           ### Number of nodes needed for the job
#SBATCH --ntasks-per-node=14 ### Number of tasks to be launcged per Node

cd /projects/bgmp/mdriskil/deduper

files=`ls -1 ./sorted_*`
for file in ${files}; do
         trimmedfilename=`echo ${file} | cut -d "/" -f2 | cut -d "_" -f2 | cut -d "." -f1`;
./driskill_deduper.py -f ${file} -o ${trimmedfilename} -u STL96.txt;
done
exit
