#!/bin/bash
#SBATCH --job-name=pc                   # Job name
#SBATCH --partition=sixhour             # Partition Name (Required)
#SBATCH --mail-type=ALL                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=gwwilson@ku.edu     # Where to send mail	
#SBATCH --ntasks=1                      # Run 1 task on one node
#SBATCH --cpus-per-task=12              # Number of threads to use
#SBATCH --mem=2gb                       # Job memory request
#SBATCH --time=0-06:00:00               # Time limit days-hrs:min:sec
#SBATCH --output=pc_%A_%a.log           # Standard output
#SBATCH --error=pc_%A_%a.err           # Standard error
#SBATCH --array=1-12                    # Job array

echo 'SLURM_ARRY_TASK_ID: '$SLURM_ARRAY_TASK_ID

LIST="MinBias2018_Rsim_GWW_ALL-J-"$SLURM_ARRAY_TASK_ID
NUMFILES=0
NTHREADS=12
VERSION=6

pwd
hostname
date

echo "PATH: "
echo $PATH
 
echo "LD_LIBRARY_PATH"
echo $LD_LIBRARY_PATH

echo "Run conversion code "

python2 --version

#Need to source ROOT

. ./setenv.sh

echo $ROOTSYS
 
MYWDIR=$PWD
cd ${MYWDIR}
echo "Now in directory "
pwd

echo "Start execution"

./runmyanalysis.sh ${LIST} ${NUMFILES} ${NTHREADS} ${VERSION}

date

exit
