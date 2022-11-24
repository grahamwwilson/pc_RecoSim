#!/bin/bash
#
# Collate jobs and add output histograms
#
# Need to source ROOT for hadd step
. ./setenv.sh

mkdir Collating/MC-V5
cd Collating/MC-V5

execdir=../../pc_ExecutionDirectory
jobfile=../../Lists/jobs.list
masterlist="MinBias2018_Rsim_GWW_ALL"
masterlistfile=${masterlist}.list

command="hadd Outfile.root"

for job in $(cat ${jobfile})
do
# Make symbolic link in collation subdirectory
   echo 'job: '$job
   mydir=MinBias2018_Rsim_GWW_ALL-J-${job}-0-12-V5
   ofile=" Outfile${job}.root"
   ln -s ${execdir}/${mydir}/Outfile.root ${ofile}
   command+="${ofile}"
done

echo ${command}
${command}

exit
