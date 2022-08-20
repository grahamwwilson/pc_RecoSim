#!/bin/sh
# runmyanalysis.sh

date

CODEDIR=$PWD

INPUT=${1-local_test3}
NUMFILES=${2-0}
NTHREADS=${3-1}
VERSION=${4-1}
MACRO=${5-runmacro.py}

echo "runmyanalysis.sh Input Args"
echo ${INPUT}
echo ${NUMFILES}
echo ${NTHREADS}
echo ${VERSION}
echo ${MACRO}

INPUTLIST=${CODEDIR}/Lists/${INPUT}.list

#Execute this from execution directory, so that we can have several output files in parallel
EXEDIR=$PWD/pc_ExecutionDirectory/${INPUT}-${NUMFILES}-${NTHREADS}-V${VERSION}

rm -r ${EXEDIR}
# Check if it exists. If not make it.
mkdir ${EXEDIR}

cd ${EXEDIR}
pwd

# input is: numfiles, numthreads,  yourdata.list, execution-path
# (numfiles==0 => read all files specfied in yourdata.list)
# Note any change to the thread count needs to also be in the job description file ..
python2 ${CODEDIR}/${MACRO} ${NUMFILES} ${NTHREADS} ${INPUTLIST} ${CODEDIR}

date

exit
