#ifndef EICvalueconst_H
#define EICvalueconst_H

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
static const double mass_chargedpion         = 0.13956995;     //(GeV/c^2)
static const double hbarc                    = 0.197327053;    //(GeV·fm)
static const double pi                       = TMath::Pi();
static const double alpha                    = 1/137.0359895;  //dimensionless, electromagnetic coupling const
//double Gamma                    = 106.6;          //dimensionless, Lorentz factor
//double Gamma_electron = 6849.3151; //3.5 GeV
//double Gamma_proton = 21.321962; //20 GeV
//double Gamma_electron = 9784.7358; // 5 GeV
//double Gamma_electron = 6849.3151;//3.5 GeV
//double Gamma_electron = 19569.472; // 10GeV
double Gamma_electron = 35225.; // 18GeV
double Gamma_electron_target = 2.*Gamma_electron*Gamma_electron-1.;
//double Gamma_electron = 106;
//double Gamma_electron = 200;//test
//double Gamma_proton = 26.652452; //25 GeV
//double Gamma_proton = 106.60981; //100GeV
double Gamma_proton = 43.71; //41GeV
//double Gamma_proton = 17.057;//16 GeV
//double Gamma_proton = 11.034;//10.35 GeV
//double Gamma_proton = 27.718;//26 GeV
//double Gamma_proton = 11.727;//11 GeV
//double Gamma_proton = 21.322;//20 GeV
//double Gamma = 2694.629; //5.02TeV
//double Gamma = 1471; // 2.76TeV
//double Gamma = 33.495; //62.4GeV
//const double	    t_perp_x		     = 0.85;	       // t_perp x-component fraction of t_perp
const double        deltar_np                = 0.14817259;
static const double mass_proton              = 0.938;          //(GeV/c^2)
static const double Z_proton                 = 1;
double              energy_proton            = mass_proton*Gamma_proton;
double              bmin                     = 12.0;           //(fm)
double              bmax                     = 200;            //(fm)
const double 	    Emin                     = 1.E-3;          //(GeV)
const double	    Emax                     = 110;            //(GeV)
//const double 	    Emin                     = 1.E-5;          //(GeV)
//const double	    Emax                     = 200;            //(GeV)
const double        Eminl                    = TMath::Log(Emin);
const double        Emaxl                    = TMath::Log(Emax);
const double        Ebins                    = 100;
int                 NbinsR                   = 30;
int                 NbinsPhi                 = 20;
double              _channelMass             = 0.7685;
static const double crosspp                  = 4.2;//fm^2
double              _wideWmax                = _channelMass+6*width_rho0;
double              _wideWmin                = 2.0*mass_chargedpion;
//double              NQ                       = 0.9;
//double              GT                       = 1;
//double              NQ                       = 1.8;
//double              GT                       = 0.5;
//double              NQ                       = 1.0/0.105;
double              NQ                       = 1/0.3;
double              GT                       = 0.5;
double              Y_current                = 0;
double              rho_r_normalization_parameter;
double              breitWigner_normalization_parameter;
double              prob_nohadron;
double              prob_photon;
double              Ep       = mass_proton*Gamma_proton;
//double              thicknessnorm = 0.16039006;//Pb
double              thicknessnorm = 0.16934601;//Au
//double              thicknessnorm = 0.15348210;//nuetron skin
//double              thicknessnorm = 0.159;
double              Wignernorm;
double              rhonorm  = 0.087;
double              c1_rho0 = 2.09;
double              c2_rho0 = 0.73*1.e-2;  // GeV^-2
double              c1_Jpsi = 2.36;
double              c2_Jpsi = 0.29*1.e-2;  // GeV^-2
//double              c1_Jpsi = 2.45;
//double              c2_Jpsi = 0.00084;
double              c1_phi = 2.15;
double              c2_phi = 0.74*1.e-2;  // GeV^-2
					  
#endif
