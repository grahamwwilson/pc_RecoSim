#!/bin/sh

# This is for |eta| < 0.70
echo ' '
echo '-----------------------------'
echo 'Results for |eta| < 0.70 (V3)'
echo '-----------------------------'
echo ' '

# Material budget tool numbers are with 'Realistic25ns13TeVEarly2018Collision' vertex smearing
# and are as reported on slide 24 of updated slides of July 14th, 2022.

python3 radlen.py -r "(BP   )" -n 30524  -d 18719557 -m 0.245801 -p 0.95888 -t 0.90246 
python3 radlen.py -r "(BPIX1)" -n 189492 -d 18647028 -m 1.5203   -p 0.97810 -t 0.91483
python3 radlen.py -r "(BPIX2)" -n 195393 -d 18389727 -m 1.5561   -p 0.97810 -t 0.91483 
python3 radlen.py -r "(BPIX3)" -n 195789 -d 18125168 -m 1.5630   -p 0.97810 -t 0.91483
python3 radlen.py -r "(BPIX4)" -n 201955 -d 17881366 -m 1.6268   -p 0.97810 -t 0.91483
python3 radlen.py -r "(Other)" -n 473562 -d 17641702 -m 3.967    -p 0.97810 -t 0.91483

echo 'Try Carbon'
python3 radlen.py -r "(BPIX1)" -n 189492 -d 18647028 -m 1.5203   -p 0.969614 -t 0.900854
echo 'Try Silicon' 
python3 radlen.py -r "(BPIX1)" -n 189492 -d 18647028 -m 1.5203   -p 0.984953 -t 0.929085

exit
