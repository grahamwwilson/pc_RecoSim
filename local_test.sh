#!/bin/sh
# runmyanalysis.sh

# This should work on my laptop

date

CODEDIR=$PWD

. ./setenv.sh

INPUT=${1-MinBias2018_Rsim_GWW_18}

echo ${INPUT}

INPUTLIST=${CODEDIR}/Lists/${INPUT}.list

#Execute this from execution directory, so that we can have several output files in parallel
EXEDIR=$PWD/pc_ExecutionDirectory/${INPUT}

rm -r ${EXEDIR}
# Check if it exists. If not make it.
mkdir ${EXEDIR}

cd ${EXEDIR}
pwd

#input is: num files, numthreads, yourdata.list, execution_path
# num files ==0 => read number from input list  
# Note any change to the thread count needs to also be in the job description file ..
python2 ${CODEDIR}/runmacrolocal.py 0 1 ${INPUTLIST} ${CODEDIR}

#Need to find some way of having several of these in parallel ...
#cp Outfile.root ${CODEDIR}/PC_${INPUT}.root

date

exit
