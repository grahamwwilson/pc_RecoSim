# radlenArgs.py
from argparse import ArgumentParser

# See https://stackoverflow.com/questions/26785952/python-argparse-as-a-function

def getArgs(argv=None):
# Set command line configurable parameters. Do python3 program.py -h to see this in action.
    parser = ArgumentParser(description="Calculate radiation length from observations")
    
    parser.add_argument("-r", "--region", type=str, default="BP", help="Radial region")      
    parser.add_argument("-n", "--nobs", type=int, default=2, help="Number of observed photon conversions")    
    parser.add_argument("-d", "--ngam", type=int, default=100000, help="Number of incident photons")
    parser.add_argument("-m", "--mbt", type=float, default=1.0, help="MBT estimate of radiation lengths (in per cent)")    
  
    args=parser.parse_args(argv)
#    print('(pvalueArgs.getArgs     ) Found argument list: ',args)
    
    return args
    
def showArgs(region,nobs,ngam,mbt):
# Check these are what we want
#    print('(pvalueArgs.ShowArgs):')
    
    print('region: ',region)    
    print('nobs:   ',nobs)
    print('ngam:   ',ngam)
    print('mbt:    ',mbt)    
    

    return
        
def getArguments(argv=None):
# Do 2 things at once.
# i)   set defaults and parse them using getArgs above
# ii)  set values for our program

    args = getArgs(argv)

#    print('(radlenArgs.getArguments) Assigning arguments to program variables')
   
    region = args.region    
    nobs = args.nobs
    ngam = args.ngam
    mbt = args.mbt

    return region,nobs,ngam,mbt    
