#include "PhotonCrossSections.h"
#include <random>
#include <fstream>
typedef std::mt19937 RandomNumberGenerator;
using namespace std;

int main(){

    ofstream fout;
    fout.open("Conversions.dat");

    unsigned long int seed;
    seed = 124794357L;        
/*
    cout << "Give seed value " << endl;
    cin >> seed;
*/
    cout << "Base seed set to " << seed << endl;

    RandomNumberGenerator ge(seed);
    RandomNumberGenerator gu(seed+1);

    const double GeVtoMeV = 1000.0;
    const double barnstocmsq = 1.0e-24;
    const int NPHOTONS = 1000000;
    const int NSLICES = 1000;
    const double SEVENOVERNINE = 7.0/9.0;
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
    double targett = 2.00*radlength;    
    
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
    
    fout << NPHOTONS << endl;
        
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
        double tratio = xsTotal/xsTsai;
        ppsum += ppFraction;
        tratiosum += tratio;
        tppsum += (xsPair/xsTsai);        
        
// Divide target into many slices each with same thickness
        int iresult = 0;     // survival
        for (int j=0; j<NSLICES; ++j){
             double Lumi = (NA/A)*rho*targett/double(NSLICES);   //cm^-2
             double pint1 = xsTotal*Lumi*barnstocmsq;  // interaction probability per layer = sigma*L.
             double dx = targett/double(NSLICES);      // slice thickness [cm]
             double probInt = SEVENOVERNINE * tratio * (rho*dx/radlen); // written more explicitly as (7/9) tratio (x/X0) 
//             cout << "pint1 = " << pint1 << " pint2 = " << pint2 << endl;
//             if(j==0)cout << "pint = " << pint << endl;
             double r = uniform(gu);
// Interaction in this layer
             if(r <= probInt){
                 if(r <= ppFraction*probInt){
                    iresult = 1;   // conversion in this layer
                 }
                 else{
                    iresult = 2;   // Compton scattering
                 }
//                 cout << "Interaction in layer " << j << " with result " << iresult << " E = " << EinGeV << endl;                 
             }
             if(iresult != 0)break; // Photon already interacted break out of loop
        }
        if(iresult == 0)nsurviving++;
        if(iresult == 1)nconversions++;
        if(iresult == 2)ncomptons++;
        fout << iresult << " " << EinGeV << " " << ppFraction << " " << tratio << endl;
    }
    cout << " " << endl;
    cout << "<Energy> = "      << Esum/double(NPHOTONS) << endl;
    cout << "<PP fraction> = " << ppsum/double(NPHOTONS) << endl;
    cout << "<tratio>      = " << tratiosum/double(NPHOTONS) << endl;
    cout << "<t*PP>        = " << tppsum/double(NPHOTONS) << endl;    
    
    cout << "  " << endl;
    cout << "N survivors   " << nsurviving << endl;
    cout << "N conversions " << nconversions << endl;
    cout << "N Comptons    " << ncomptons << endl;
    
    double pconv = double(nconversions)/double(NPHOTONS);
    double dpconv=sqrt(pconv*(1.0-pconv)/double(NPHOTONS));
    double ppfrac = ppsum/double(NPHOTONS);
    double tratio = tratiosum/double(NPHOTONS);
    double tpp = tppsum/double(NPHOTONS);
    
// This is version 1 in radlen.py    
    double fradl1meas = -(9.0/7.0)*log(1.0 - pconv);
    double dfradl1meas = (9.0/7.0)*dpconv/(1.0 - pconv);    
    cout << "Naive estimate of fractional radiation length (%)   = " << 100.0*fradl1meas << " +- " << 100.0*dfradl1meas << endl;    
    
// This is version 2 in radlen.py    
    double fradl2meas = -(9.0/7.0)*(1.0/tpp)*log(1.0 - pconv);
    double dfradl2meas = (9.0/7.0)*(1.0/tpp)*dpconv/(1.0 - pconv);    
    cout << "Refined estimate of fractional radiation length (%) = " << 100.0*fradl2meas << " +- " << 100.0*dfradl2meas << endl;      
        
// This is version 3 in radlen.py
    double fradlmeas = -(9.0/7.0)*(1.0/tratio)*log(1.0 - (pconv/ppfrac));
    double dfradlmeas = (9.0/7.0)*(1.0/(ppfrac*tratio))*dpconv/(1.0 - (pconv/ppfrac));    
    cout << "Best estimate of fractional radiation length (%)    = " << 100.0*fradlmeas << " +- " << 100.0*dfradlmeas << endl;
    
// Let's try simply survival probability as the estimator
    double psurv = double(nsurviving)/double(NPHOTONS);
    double dpsurv=sqrt(psurv*(1.0-psurv)/double(NPHOTONS)); 
    double fradl4meas = -(9.0/7.0)*(1.0/tratio)*log(psurv);
    double dfradl4meas = (9.0/7.0)*(1.0/tratio)*dpsurv/(psurv);    
    cout << "Survival estimate of fractional radiation length (%)= " << 100.0*fradl4meas << " +- " << 100.0*dfradl4meas << endl;                    
    
    fout.close();

    return 0;
}
