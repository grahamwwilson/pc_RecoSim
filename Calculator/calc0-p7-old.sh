#!/bin/sh
# Here we try to reproduce our previous results with the current implementation.
# We set ppfrac = 1.0 (ignore CS) and use our old version for tratio namely pp/Tsai
# This is for |eta| < 0.70
echo ' '
echo '-----------------------------'
echo 'Results for |eta| < 0.70 (V2.5)'
echo '-----------------------------'
echo ' '

# Material budget tool numbers are with 'Realistic25ns13TeVEarly2018Collision' vertex smearing
# and are as reported on slide 24 of updated slides of July 14th, 2022.

python3 radlen.py -r "(BP   )" -n 30524  -d 18719557 -m 0.245801 -p 1.0 -t 0.86557 
python3 radlen.py -r "(BPIX1)" -n 189492 -d 18647028 -m 1.5203   -p 1.0 -t 0.89493
python3 radlen.py -r "(BPIX2)" -n 195393 -d 18389727 -m 1.5561   -p 1.0 -t 0.89493 
python3 radlen.py -r "(BPIX3)" -n 195789 -d 18125168 -m 1.5630   -p 1.0 -t 0.89493
python3 radlen.py -r "(BPIX4)" -n 201955 -d 17881366 -m 1.6268   -p 1.0 -t 0.89493
python3 radlen.py -r "(Other)" -n 473562 -d 17641702 -m 3.967    -p 1.0 -t 0.89493 

exit
