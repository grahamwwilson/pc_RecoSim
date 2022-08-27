#!/bin/sh

# This is for |eta| < 0.10
echo ' '
echo '------------------------'
echo 'Results for |eta| < 0.10'
echo '------------------------'
echo ' '

# Material budget tool numbers are with 'Realistic25ns13TeVEarly2018Collision' vertex smearing
# and are as reported on slide 23 of updated slides of July 14th, 2022.

python3 radlen.py -r "(BP   )"    -n 4115  -d 2679647 -m 0.227200 -c 4
python3 radlen.py -r "(BPIX1)" -n 24636 -d 2670018 -m 1.4008
python3 radlen.py -r "(BPIX2)" -n 25384 -d 2636163 -m 1.3997
python3 radlen.py -r "(BPIX3)" -n 24272 -d 2601274 -m 1.3631
python3 radlen.py -r "(BPIX4)" -n 24441 -d 2570358 -m 1.3647
python3 radlen.py -r "(Other)" -n 58079 -d 2541106 -m 3.227

exit
