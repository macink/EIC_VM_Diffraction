#ifndef FormFactor_resolution_add_wedge_1D2_H
#define FormFactor_resolution_add_wedge_1D2_H

#include "EICvalueconst.h"

using namespace std;




/********************************************** 
 
 * Form Factor with detector resolution and wedge cut 1D: |F(t)|^2
 
 * Parameters:
    Vo: Nuclear density (hard sphere approximation) [fm^-3]
    A: Atomic mass number
    R: Radius of nucleus [fm]
    a0: Range of Yukawa potential [fm]
 * Variables:
    q: Transverse momentum >> q = sqrt(qx^2+qy^2) [GeV]
    t: Four-momentum transfer squared >> t = tx+ty = q^2 [GeV^2]

**********************************************/
class FormFactor_resolution_add_wedge_1D2
{
public:
    FormFactor_resolution_add_wedge_1D2(double t_min, double t_max, 
        double bins, double phi_min, double phi_max)
    {
        ROOT::Math::IntegratorOneDimOptions::SetDefaultAbsTolerance(1e-4);
        ROOT::Math::IntegratorOneDimOptions::SetDefaultRelTolerance(1e-2);
        ROOT::Math::IntegratorOneDimOptions::SetDefaultIntegrator("Gauss");

        
        // Initialize |F(t)|^2 wedge cut TF1
        cutFFt2 = new TF1("Cut Form Factor: |F(t)|^{2}", [this, phi_min, phi_max] (double *var, double *par)
        {
            double t = var[0];
            TF1 fft2("", FF_cut_wRes, phi_min, phi_max);
            fft2.SetParameters(t);
            double integral = fft2.Integral(phi_min, phi_max, 1e-12);
            return integral;
        }, t_min, t_max, 0);

        // Initialize hisogram 
        hist1D = new TH1D("", "", bins, t_min, t_max);
        hist1D->Sumw2();
        double step3 = t_max/bins;
        for(int i=0;i<bins;i++)
        {
            double t = hist1D->GetXaxis()->GetBinCenter(i+1);
            //double t = (i+1)*step3;
            //double q = sqrt(t);
            //double val = wedge_FF(q, sigma, phi_min, phi_max);
            //hist1D->SetBinContent(i+1, val);
            hist1D->SetBinContent(i+1, cutFFt2->Eval(t)); 
        }

    }

    // Integrand for wedge cut
    static double FF_cut_wRes(double *x, double *par)
    {
        double t = par[0];
        double phi = x[0];  // integration variable
        double q = sqrt(t), qx_prime = q*sin(phi), qy = q*cos(phi); 
        // call form factor with resolution
	    return sqrt(qx_prime*qx_prime + qy*qy);
    }
   

    // Functions
    TF1 *getWedgeRes_fun_1D() const {return cutFFt2;}
    TH1D *getWedgeRes_hist_1D() const {return hist1D;}

private:
    double A, Vo, R, a0;
    TLorentzVector vmMC;
    double tx, ty; 
    TF1 *cutFFt2;
    TH1D *hist1D;
};



#endif






