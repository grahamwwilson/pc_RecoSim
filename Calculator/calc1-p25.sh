#!/bin/sh

# This is for |eta| < 1.25
echo ' '
echo '------------------------'
echo 'Results for |eta| < 1.25'
echo '------------------------'
echo ' '

# Material budget tool numbers are for VersionG.

python3 radlen.py -r "(BP   )" -n 64682  -d 33047557 -m 0.29085 -c 4
python3 radlen.py -r "(BPIX1)" -n 408888 -d 32902522 -m 1.9016
python3 radlen.py -r "(BPIX2)" -n 434761 -d 32388557 -m 2.0138
python3 radlen.py -r "(BPIX3)" -n 448641 -d 31827269 -m 2.1175
python3 radlen.py -r "(BPIX4)" -n 438583 -d 30804492 -m 2.4658
python3 radlen.py -r "(Other)" -n 664990 -d 27115728 -m 5.1844

exit
