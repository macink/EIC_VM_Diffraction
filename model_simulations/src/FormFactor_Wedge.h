#ifndef FormFactor_Wedge_H
#define FormFactor_Wedge_H

#include "EICvalueconst.h"

using namespace std;

/********************************************** 
 
 * Form Factor Wedge Cuts
 
 * Parameters:
    Vo: Nuclear density (hard sphere approximation) [fm^-3]
    A: Atomic mass number
    R: Radius of nucleus [fm]
    a0: Range of Yukawa potential [fm]
 * Variables:
    t: Four-momentum transfer squared >> t = tx+ty = q^2 [GeV^2]
    phi: angle of wedge to cut out from y-axis
        qx = q sin(phi)
        qy = q cos(phi)
        **Do in terms of q then convert back to t

**********************************************/
class FormFactor_Wedge
{
public:
    FormFactor_Wedge(double A_init, double Vo_init, double R_init, double a0_init, double t_min, double t_max, double phi_min, 
        double phi_max, double bins): A(A_init), Vo(Vo_init), R(R_init), a0(a0_init)
    {
        // Set integration preferences
        ROOT::Math::IntegratorOneDimOptions::SetDefaultAbsTolerance(1e-4);
        ROOT::Math::IntegratorOneDimOptions::SetDefaultRelTolerance(1e-2);
        ROOT::Math::IntegratorOneDimOptions::SetDefaultIntegrator("Gauss");

        // Initialize |F(t)|^2 wedge cut TF1
        wedge_cut_formFactor_TF1 = new TF1("", [this, phi_min, phi_max] (double *var, double *par)
        {
            double t = var[0];  // t is variable
            TF1 wedge_formFactor_integral("", calc_formFactor_wedge, phi_min, phi_max, 5);
            wedge_formFactor_integral.SetParameters(par[0], par[1], par[2], par[3], t);
            double integral = wedge_formFactor_integral.Integral(phi_min, phi_max, 1e-12);
            return integral;
        }, t_min, t_max, 4);
        wedge_cut_formFactor_TF1->SetParameters(A,Vo,R,a0);

        // Initialize histogram for TF1
        wedge_cut_formFactor_TH1D = new TH1D("", "", bins, t_min, t_max);
        wedge_cut_formFactor_TH1D->Sumw2();
        double step = t_max/bins;
        for(int i=0;i<bins;i++) 
        {
            double t = (i+1)*step;
            double q = sqrt(t);
            double val = wedge_integral(q, phi_min, phi_max, Vo, A, R, a0);
            wedge_cut_formFactor_TH1D->SetBinContent(i+1, val); 
        }
}


    // Calculations for form factor
    static double calc_formFactor_wedge(double *var, double *par)
    {
        double A = par[0], Vo = par[1], R = par[2], a0 = par[3], t = par[4];
        double phi = var[0];
        double q = sqrt(t), qx = q*sin(phi), qy = q*cos(phi);
	    const double arg1 = q*R / hbarc;  
	    const double arg2 = hbarc / q;
	    const double arg3  = (sin(arg1)-arg1*cos(arg1)) * 4*pi*Vo * arg2*arg2*arg2 / A; 
	    double result = 1/(2*pi)*arg3 / (1. + (a0*a0 * q*q/(hbarc*hbarc))); 
        return result*result; // [-]
    }
    // Integrand to calculate |F(q)|^2 for wedge cut
    static double wedge_integrand(double *var, double *par)
    {
        double A = par[0], Vo = par[1], R = par[2], a0 = par[3], t = par[4];
        //double qq = par[0];
        double phi = var[0];  // integration variable
        // put q in terms of theta
        double q = sqrt(t), qx = q*sin(phi), qy = q*cos(phi); 
        double vars[] = {phi};
        double pars[] = {A, Vo, R, a0, t};
	    return calc_formFactor_wedge(vars, pars);  // Integrand
    }
    // Integration over phi
    double wedge_integral(double q, double phi_min, double phi_max, double Vo, double A, double R, double a0)
    {
        TF1 *f = new TF1("f", wedge_integrand, phi_min, phi_max, 5);  
        f->SetParameters(q, Vo, A, R, a0);
        return f->Integral(phi_min, phi_max, 1e-12);
    }
   

    // Functions
    TF1 *getWedgeCut1D_fun() const {return wedge_cut_formFactor_TF1;}
    TH1D *getWedgeCut1D_hist() const {return wedge_cut_formFactor_TH1D;}

private:
    double A, Vo, R, a0;
    TF1 *wedge_cut_formFactor_TF1;
    TH1D *wedge_cut_formFactor_TH1D;
};

#endif
