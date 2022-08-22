#!/bin/sh

# This is for |eta| < 1.25
echo ' '
echo '------------------------'
echo 'Results for |eta| < 1.25'
echo ' (note do not yet have MBT results for this acceptance region) '
echo '------------------------'
echo ' '

# Material budget tool numbers are with 'Realistic25ns13TeVEarly2018Collision' vertex smearing
# and are as reported on slide 24 of updated slides of July 14th, 2022 for |eta| < 0.7!

python3 radlen.py -r "BP"    -n 64682  -d 33047557 -m 0.2458
python3 radlen.py -r "BPIX1" -n 408888 -d 32902522 -m 1.5203
python3 radlen.py -r "BPIX2" -n 434761 -d 32388557 -m 1.5561
python3 radlen.py -r "BPIX3" -n 448641 -d 31827269 -m 1.5630
python3 radlen.py -r "BPIX4" -n 438583 -d 30804492 -m 1.6268
python3 radlen.py -r "Other" -n 664990 -d 27115728 -m 3.967

exit
