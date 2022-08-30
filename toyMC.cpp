#include "PhotonCrossSections.h"
#include <random>
typedef std::mt19937 RandomNumberGenerator;
using namespace std;

int main(){

    unsigned long int seed;
/*
    cout << "Give seed value " << endl;
    cin >> seed;
    cout << "Base seed set to " << seed << endl;
*/
    seed = 13579135L;    
    RandomNumberGenerator ge(seed);
    RandomNumberGenerator gu(seed+1);

    const double GeVtoMeV = 1000.0;
    const double barnstocmsq = 1.0e-24;
    const int NPHOTONS = 10000000;
    const int NSLICES = 100;
    const double ECUT = 0.4;   
    const double Z = 6.0;             // Choose Carbon for now.
    const double rho = 2.210;         // g/cc
    const double A = 12.0107;         // Atomic weight g/mol
    const double NA = 6.02214075e23;  // Avogadro const (/mol)
    double xsTsai = ComputeApproxCrossSectionPerAtom(Z);  // barns/atom    
    double radlen = (7.0/9.0)*(A/(NA*xsTsai))*1.0e24;
    double radlength = radlen/rho;
    
    cout << "radlen    = " << radlen    << " g/cm^2 " << endl;
    cout << "radlength = " << radlength << " cm     " << endl;
    
// Set target thickness
    double targett = 0.015*radlength;    
    
    exponential_distribution<double> expo(3.0);
    uniform_real_distribution<double> uniform;
    
    double Esum = 0.0;
    double ppsum = 0.0;
    double tratiosum = 0.0;
    double tppsum = 0.0;
    
    double E, EinGeV;

    int nsurviving = 0;
    int nconversions = 0;
    int ncomptons = 0;
        
    for (int i=0; i<NPHOTONS; ++i) {
// Choose random energy according to particular distribution. Set minimum at 0.4 GeV.
        EinGeV = 0.0;
        while (EinGeV < ECUT){
            EinGeV = expo(ge);
        }
        Esum += EinGeV;
        E = EinGeV*GeVtoMeV;
        double xsPair    = ComputeCrossSectionPerAtom(E,Z);
        double xsCompton = ComputeComptonCrossSectionPerAtom(E,Z);
        double xsTotal   = xsPair + xsCompton;
        double ppFraction = xsPair/xsTotal;
        ppsum += ppFraction;
        tratiosum += (xsTotal/xsTsai);
        tppsum += (xsPair/xsTsai);        
        
// Divide target into many slices each with same thickness
        int iresult = 0;     // survival
        for (int j=0; j<NSLICES; ++j){
             double Lumi = (NA/A)*rho*targett/double(NSLICES);   //cm^-2
             double pint = xsTotal*Lumi*barnstocmsq;  // interaction probability per layer
//             if(j==0)cout << "pint = " << pint << endl;
             double r = uniform(gu);
// Interaction in this layer
             if(r <= pint){
                 if(r <= ppFraction*pint){
                    iresult = 1;   // conversion in this layer
                 }
                 else{
                    iresult = 2;   // Compton scattering
                 }
                 cout << "Interaction in layer " << j << " with result " << iresult << " E = " << EinGeV << endl;                 
             }
             if(iresult != 0)break; // Photon already interacted break out of loop
        }
        if(iresult == 0)nsurviving++;
        if(iresult == 1)nconversions++;
        if(iresult == 2)ncomptons++;
    }
    cout << " " << endl;
    cout << "<Energy> = "      << Esum/double(NPHOTONS) << endl;
    cout << "<PP fraction> = " << ppsum/double(NPHOTONS) << endl;
    cout << "<tratio>      = " << tratiosum/double(NPHOTONS) << endl;
    cout << "<t*PP>        = " << tppsum/double(NPHOTONS) << endl;    
    
    cout << "  " << endl;
    cout << "Nsurviving   " << nsurviving << endl;
    cout << "Nconversions " << nconversions << endl;
    cout << "N Comptons   " << ncomptons << endl;
    
    double pconv = double(nconversions)/double(NPHOTONS);
    double dpconv=sqrt(pconv*(1.0-pconv)/double(NPHOTONS));
    double ppfrac = ppsum/double(NPHOTONS);
    double tratio = tratiosum/double(NPHOTONS);
    double radlmeas = -(9.0/7.0)*(1.0/tratio)*log(1.0 - (pconv/ppfrac));
    double dradlmeas = (9.0/7.0)*(1.0/(ppfrac*tratio))*dpconv/(1.0 - (pconv/ppfrac));    

    cout << "Estimate of radiation length (%) = " << 100.0*radlmeas << " +- " << 100.0*dradlmeas << endl;

    return 0;
}
