// Adapt the GEANT4 code from G4BetheHeitlerModel.cc for the total cross-section
// Calculates the microscopic cross section in GEANT4 internal units.
// A parametrized formula from L. Urban is used to estimate
// the total cross section in ComputeCrossSectionPerAtom.
// It gives a good description of the data from 1.5 MeV to 100 GeV.
// below 1.5 MeV: sigma=sigma(1.5MeV)*(GammaEnergy-2electronmass)
//                                   *(GammaEnergy-2electronmass)

#include <iostream>
#include <math.h>
#include <string>
#include <iomanip>
using namespace std;

// Lrad function used in the radiation length definition of 
// Y-S. Tsai equation 3.66 on p827 of 
// Review of Modern Physics, Vol 46 (4) 1974 pp815-851.
double Lrad(double dZ){
    int Z;
    double val;
    const double ONETHIRD = 1.0/3.0;

    Z = int(dZ+0.5);

    if(Z>=5){
        val = log(184.15/(pow(dZ,ONETHIRD)));
    }
    else if(Z==1){
        val = 5.31;
    }
    else if(Z==2){
        val = 4.79;
    }
    else if(Z==3){
        val = 4.74;
    }
    else if(Z==4){
        val = 4.71;
    }
    return val;
}

// Lradprime function used in the radiation length definition of 
// Y-S. Tsai equation 3.66 on p827 of 
// Review of Modern Physics, Vol 46 (4) 1974 pp815-851.
double Lradprime(double dZ){
    int Z;
    double val;
    const double TWOTHIRDS = 2.0/3.0;

    Z = int(dZ+0.5);
    if(Z>=5){
        val = log(1194.0/(pow(dZ,TWOTHIRDS)));
    }
    else if(Z==1){
        val = 6.144;
    }
    else if(Z==2){
        val = 5.621;
    }
    else if(Z==3){
        val = 5.805;
    }
    else if(Z==4){
        val = 5.924;
    }
    return val;
}

double RadiationLengthConstant(){

//  Calculate the coefficient (approximately 716.408) in Tsai equation 3.66 
//  using modern values of fundamental constants.
//  Units here are mol cm^-2

    static const double alp=1.0/137.035999084;   // alpha
    static const double re=2.8179403262e-13;     // electron classical radius (cm)
    static const double NAvogadro=6.02214075e23; // Avogadro's number (mol^-1)
   
    static const double rlConstant = 1.0/(4.0*alp*re*re*NAvogadro);
   
    return rlConstant;
}   

double ComputeApproxCrossSectionPerAtom(double Z){
//
//  Define sigmaINF = 7/9 (A / (X0 * NA )) in barns/atom  
//  where the radiation length, X0, in g cm^-2 is taken from Tsai eqn 3.66.
//  This is Tsai equation 3.73, that neglects the -(2/21)*(Z**2 + Z) term of equation 3.61 
//  and is an approximation for the total pair production cross-section at infinite energy.  

//  Coulomb correction function
    static const double alp=1.0/137.035999084;   // alpha
    static const double re=2.8179403262e-13;     // electron classical radius (cm)
    double azsq = alp*alp*Z*Z;
    double fc = azsq*(0.20206 -0.0369*azsq +0.0083*azsq*azsq -0.0020*azsq*azsq*azsq + (1.0/(1.0+azsq)) );
   
    double sigmaINF = (7.0/9.0)*(alp*re*re)*4.0*( Z*Z*(Lrad(Z) - fc) + Z*Lradprime(Z) );
    sigmaINF = sigmaINF/1.0e-24;
    return sigmaINF;    //units of barns
}

double ComputeBetterCrossSectionPerAtom(double Z){
//
//  Define sigmaINF using equation 3.61 of Tsai. 
//  This is also an approximation for the total pair production cross-section at infinite energy 
//  but includes the correction term with coefficient -2/21 in eqn 3.61 
//  that is neglected in the radiation length formalism. 
//  While this is interesting to compute, and a more exact approximation in the 
//  complete screening high energy limit, the Tsai radiation length is defined excluding this 
//  correction term, so it is not directly relevant to extracting what is accepted as the 
//  definition of radiation length.

//  Coulomb correction function
    static const double alp=1.0/137.035999084;   // alpha
    static const double re=2.8179403262e-13;     // electron classical radius (cm)
    double azsq = alp*alp*Z*Z;
    double fc = azsq*(0.20206 -0.0369*azsq +0.0083*azsq*azsq -0.0020*azsq*azsq*azsq + (1.0/(1.0+azsq)) );
   
//   double l1p = 20.863 - (4.0/3.0)*log(Z) - 4.0*fc;
//   double l2p = 28.352 - (8.0/3.0)*log(Z);
   
    const double COEFF = 1.0/42.0;
    double sigmaINF = (7.0/9.0)*(alp*re*re)*4.0*( Z*Z*( (Lrad(Z) - fc) - COEFF) + Z*( Lradprime(Z) - COEFF) );
    sigmaINF = sigmaINF/1.0e-24;
    return sigmaINF;    //units of barns
}
 
double ComputeCrossSectionPerAtom(double gammaEnergy, double Z){ 
//
//  Returns cross-section per atom in barns. Basically from G4BetheHeitlerModel.cc in Geant4.
// 
//  gammaEnergy in MeV
//  Z: atomic number (number of protons for elemental materials
//
    double xSection = 0.0 ;
    // short versions
    static const double kMC2  = 0.51099906;  // MeV
    // zero cross section below the kinematical limit: Eg<2mc^2
    if (Z < 0.9 || gammaEnergy <= 2.0*kMC2) { return xSection; }
    //
    static const double gammaEnergyLimit = 1.5;  // MeV
    // set coefficients a, b c in units of microbarn
    static const double a0 =  8.7842e+2;   // microbarn
    static const double a1 = -1.9625e+3;   // microbarn 
    static const double a2 =  1.2949e+3;   // microbarn
    static const double a3 = -2.0028e+2;   // microbarn 
    static const double a4 =  1.2575e+1;   // microbarn 
    static const double a5 = -2.8333e-1;   // microbarn
  
    static const double b0 = -1.0342e+1;   // microbarn
    static const double b1 =  1.7692e+1;   // microbarn
    static const double b2 = -8.2381   ;   // microbarn
    static const double b3 =  1.3063   ;   // microbarn
    static const double b4 = -9.0815e-2;   // microbarn
    static const double b5 =  2.3586e-3;   // microbarn
  
    static const double c0 = -4.5263e+2;   // microbarn
    static const double c1 =  1.1161e+3;   // microbarn 
    static const double c2 = -8.6749e+2;   // microbarn
    static const double c3 =  2.1773e+2;   // microbarn 
    static const double c4 = -2.0467e+1;   // microbarn
    static const double c5 =  6.5372e-1;   // microbarn
    // check low energy limit of the approximation (1.5 MeV)
    double gammaEnergyOrg = gammaEnergy;
    if (gammaEnergy < gammaEnergyLimit) { gammaEnergy = gammaEnergyLimit; }
    // compute gamma energy variables
    const double x  = log(gammaEnergy/kMC2);
    const double x2 =  x*x; 
    const double x3 = x2*x;
    const double x4 = x3*x;
    const double x5 = x4*x;
    //
    const double F1 = a0 + a1*x + a2*x2 + a3*x3 + a4*x4 + a5*x5;
    const double F2 = b0 + b1*x + b2*x2 + b3*x3 + b4*x4 + b5*x5;
    const double F3 = c0 + c1*x + c2*x2 + c3*x3 + c4*x4 + c5*x5;     
    // compute the approximated cross section 
    xSection = (Z + 1.)*(F1*Z + F2*Z*Z + F3);
    // check if we are below the limit of the approximation and apply correction
    if (gammaEnergyOrg < gammaEnergyLimit) {
        const double dum = (gammaEnergyOrg-2.*kMC2)/(gammaEnergyLimit-2.*kMC2);
        xSection *= dum*dum;
    }
    // make sure that the cross section is never negative
    xSection = 1.0e-6*max(xSection, 0.); 
    return xSection; // units of barns
}
