# radlenArgs.py
from argparse import ArgumentParser

# See https://stackoverflow.com/questions/26785952/python-argparse-as-a-function

def getArgs(argv=None):
# Set command line configurable parameters. Do python3 program.py -h to see this in action.
    parser = ArgumentParser(description="Calculate radiation length from observations")
    parser.add_argument("-n", "--nobs", type=int, default=2, help="Number of observed photon conversions")    
    parser.add_argument("-d", "--ngam", type=float, default=2.0e6, help="Number of incident photons")
         
    args=parser.parse_args(argv)
    print('(pvalueArgs.getArgs     ) Found argument list: ',args)
    
    return args
    
def showArgs(nobs,ngam):
# Check these are what we want
    print('(pvalueArgs.ShowArgs    ) Program has set')
    print('nobs:   ',nobs)
    print('ngam:   ',ngam)
    return
        
def getArguments(argv=None):
# Do 2 things at once.
# i)   set defaults and parse them using getArgs above
# ii)  set values for our program

    args = getArgs(argv)

    print('(radlenArgs.getArguments) Assigning arguments to program variables')
    
    nobs = args.nobs
    ngam = args.ngam
   
    return nobs,ngam    
