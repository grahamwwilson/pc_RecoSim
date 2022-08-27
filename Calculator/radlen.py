# radlen.py
#
# Simple radiation length calculator
#
import myPythonCheck
import radlenArgs
import math

myPythonCheck.Check()                         # Enforce use of python3

corrtype, region, nobs, ngam, mbt = radlenArgs.getArguments(None)
radlenArgs.showArgs(corrtype, region, nobs, ngam, mbt)

# Calculate conversion probability

# First define cross-section ratio (sigmaPP/sigmaTsai)
factor = 0.8895                     # default based on Freon (lambda/X0 = 2.958)
if corrtype == 4:
    factor = 0.8596                 # based on Beryllium
elif corrtype == -1:
    factor = 1.0                    # No correction

print('Factor set to ',factor)

p = nobs/ngam
dp = math.sqrt(p*(1.0-p)/ngam)

print('Observed conversion probability in %',100.0*p,'+-',100.0*dp,' ( ',100.0*dp/p,')')

# Now translate this to radiation lengths using the naive formula, namely
# x/X0 = - (9/7)*log(1-pconv)

radl = -(9.0/7.0)*(1.0/factor)*math.log(1.0 - p)
dradl = (9.0/7.0)*(1.0/factor)*dp/(1.0 - p)
print('Corresponding to ',100.0*radl,'+-',100.0*dradl,' radiation lengths in %',' ( ',100.0*dradl/radl,')')
print('(RecoSIM/MBT ratio) for ',region,100.0*radl/mbt,'+-',100.0*dradl/mbt)
print(' ')
