#include "TF2.h"
#include "TFile.h"
#include "Math/WrappedMultiTF1.h"
#include "Math/AdaptiveIntegratorMultiDim.h"
#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLorentzVector.h"
#include "TRandom.h"
#include "TGraph.h"
#include "Math/Integrator.h"
#include "Math/IntegratorMultiDim.h"
#include "Math/AllIntegrationTypes.h"
#include "Math/Functor.h"
#include "Math/GaussIntegrator.h"
#include "finalstate.h"
#include "nucleus.h"

static const double mass_electric            = 0.00051099907;  //(GeV/c^2)
static const double mass_chargedpion         = 0.93956995;     //(GeV/c^2)
static const double hbarc                    = 0.197327053;    //(GeV·fm)
static const double pi                       = TMath::Pi();
static const double alpha                    = 1/137.0359895;  //dimensionless, electromagnetic coupling const
//double Gamma                    = 106.6;          //dimensionless, Lorentz factor
double Gamma = 11.727079;// 11 GeV
//double Gamma = 2694.629; //5.02TeV
//double Gamma = 1471; // 2.76TeV
//double Gamma = 33.495; //62.4GeV
//double Gamma                    = 2665.2452;//5TeV
static const double mass_proton              = 0.938;          //(GeV/c^2)
static const double Z_proton                 = 1;
double              energy_proton            = mass_proton*Gamma;
double              bmin                     = 12.0;           //(fm)
double              bmax                     = 200;            //(fm)
const double 	    Emin                     = 1.E-3;          //(GeV)
const double	    Emax                     = 200;            //(GeV)
//const double 	    Emin                     = 1.E-5;          //(GeV)
//const double	    Emax                     = 200;            //(GeV)
const double        Eminl                    = TMath::Log(Emin);
const double        Emaxl                    = TMath::Log(Emax);
//const double        Ebins                    = 100;
const double        Ebins                    = 2000;
int                 NbinsR                   = 30;
int                 NbinsPhi                 = 20;
double              _channelMass             = 0.7685;
//static const double crosspp                  = 4.2;//fm^2
static const double crosspp = 3.8;
//static const double crosspp                  = 6.18;//2.76TeV
double              _wideWmax                = _channelMass+6*width_rho0;
double              _wideWmin                = 2.0*mass_chargedpion;
//double              NQ                       = 0.9;
//double              GT                       = 1;
double              NQ                       = 1.8;
double              GT                       = 0.5;
//double              NQ                       = 1.0/0.105;
//double              NQ                       = 1/0.3;
//double              NQ                       = 2.4;
//double              GT                       = 0.5;
double              Y_current                = 0;
double              rho_r_normalization_parameter;
double              breitWigner_normalization_parameter;
//double              prob_nohadron;
double              prob_photon;
double              Ep       = mass_proton*Gamma;
//double              thicknessnorm = 0.1604;//Pb
double              thicknessnorm = 0.16934601;//Au
//double              thicknessnorm = 0.169;
//double              thicknessnorm = 0.159;
double              Wignernorm = 0.0887;
double              rhonorm  = 0.087;



