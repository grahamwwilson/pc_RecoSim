#include "PairCrossSections.h"

int main(){

    double rlConstant = RadiationLengthConstant();
    cout << " " << endl;
    cout << "GEANT4 total pair-production cross-section estimates vs (E,Z)" << endl;
    cout << "These use the 18-parameter parametrization by L. Urban encoded in G4BetheHeitler.cpp" << endl;
    cout << "based on least-squares fits to the data of Hubbell, Gimm, Overbo (1980)." << endl;
    cout << "The range of validity is 1.5MeV < E < 100GeV and 1 <= Z <= 100" << endl;
    cout << "This may well be what is used in the CMS simulations (I need to check) " << endl;
    cout << " " << endl;
    cout << "The accuracy is estimated to be at worst 5% with a mean value of about 2.2% " << endl;
    cout << " " << endl;
    cout << "Radiation length constant calculated as " << setw(12) << rlConstant << endl;
    
// According to 
// https://lss.fnal.gov/archive/2021/slides/fermilab-slides-21-046-cms.pdf    
// CMS uses the FTFP_BERT_EMM physics list. So the electro-magnetic part is EMM.
// Also see CMS DP-2020/050
// Looks like gamma conversion model is BetheHeitlerLPM

    const double GeVtoMeV = 1000.0;
    const int NEL = 13;
    const int NENERGIES = 10;    
    string names[NEL] = { "Hydrogen", "Helium", "Beryllium", "Carbon", "Freon", "F/Ne", "Neon", "Aluminum", "Silicon", "Iron", "Copper", "Lead", "Uranium" };
    double elements[NEL] = { 1, 2, 4, 6, 9, 9.5, 10, 13, 14, 26, 29, 82, 92 };        
    double energies[NENERGIES] = {0.01, 0.03, 0.1, 0.3, 0.74, 1.0, 3.0, 10.0, 30.0, 100.0};
    
    for (int i=0; i<NEL; ++i) {
        double Z = elements[i];
        double sigmaINFA = ComputeApproxCrossSectionPerAtom(Z);
        double sigmaINFB = ComputeBetterCrossSectionPerAtom(Z);       
        cout << " " << endl;
        cout << " " << endl;
        cout << names[i] << " Z = " << elements[i] << " sigmaINFA, sigmaINFB = " 
             << setw(9) << sigmaINFA << " " << setw(9) << sigmaINFB << " barns/atom " << endl;
        cout << " " << endl;
        for (int j=0; j<NENERGIES; ++j) {
            double E = energies[j]*GeVtoMeV;
            double xs = ComputeCrossSectionPerAtom(E,Z);
            cout << "E = " << setw(4) << energies[j] << " GeV   Cross-section: " << setw(9) << xs << " barns/atom " 
                 << " RatioA: " << setw(9) << xs/sigmaINFA << " RatioB " << setw(9) << xs/sigmaINFB << endl;           
        } 
    }
    return 0;
}
