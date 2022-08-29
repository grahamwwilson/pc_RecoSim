# radlen.py
#
# Radiation length calculator
#
# ppfrac = xsPair/xsTotal
# tratio = xsTotal/xsTsai
#
import myPythonCheck
import radlenArgs
import math

myPythonCheck.Check()                         # Enforce use of python3

ppfrac, tratio, region, nobs, ngam, mbt = radlenArgs.getArguments(None)
radlenArgs.showArgs(ppfrac, tratio, region, nobs, ngam, mbt)

p = nobs/ngam
dp = math.sqrt(p*(1.0-p)/ngam)

print('Observed conversion probability in %',100.0*p,'+-',100.0*dp,' ( ',100.0*dp/p,')')

# Now translate this to radiation lengths using increasingly more sophisticated formulae
# 1. x/X0 = - (9/7)*log(1-pconv)
# 2. x/X0 = - (9/7)*(1/tratio)*log(1-pconv)      (was initially implemented ignoring CS, ie. ppfrac=1).
# 3. x/X0 = - (9/7)*(1/tratio)*log(1 - (pconv/ppfrac))

radl = -(9.0/7.0)*(1.0/tratio)*math.log(1.0 - (p/ppfrac))
dradl = (9.0/7.0)*(1.0/(ppfrac*tratio))*dp/(1.0 - (p/ppfrac))
print('Corresponding to ',100.0*radl,'+-',100.0*dradl,' radiation lengths in %',' ( ',100.0*dradl/radl,')')
print('(RecoSIM/MBT ratio) for ',region,100.0*radl/mbt,'+-',100.0*dradl/mbt)
print(' ')
