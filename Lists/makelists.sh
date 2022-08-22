#!/bin/bash
#
# For convenience let's make 12 separate list files each 
# containing the next 83 root file names in the overall list 
# file (996 total). So all files are included and no file is repeated.
#

jobfile=jobs.list
masterlist="MinBias2018_Rsim_GWW_ALL"
masterlistfile=${masterlist}.list

first=1
last=83

for job in $(cat ${jobfile})
do
   echo 'job: '$job
   echo 'first: '$first
   echo 'last:  '$last
   outfile=${masterlist}-J-${job}.list
   sed -n "${first},${last} p" ${masterlistfile} > ${outfile}
   let first+=83
   let last+=83
done

exit
