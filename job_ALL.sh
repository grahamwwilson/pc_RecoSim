#!/bin/bash
#SBATCH --job-name=conv2                # Job name
#SBATCH --partition=sixhour             # Partition Name (Required)
#SBATCH --mail-type=FAIL,END            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=gwwilson@ku.edu     # Where to send mail	
#SBATCH --ntasks=1                      # Run 1 task on one node
#SBATCH --cpus-per-task=24              # Number of threads to use
#SBATCH --mem=8gb                       # Job memory request
#SBATCH --time=0-06:00:00               # Time limit days-hrs:min:sec
#SBATCH --output=conv2_%j.log           # Standard output and error log

LIST="MinBias2018_Rsim_GWW_ALL"
NUMFILES=0
NTHREADS=24
VERSION=1

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
