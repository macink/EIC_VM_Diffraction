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
        ROOT::Math::IntegratorOneDimOptions::SetDefaultAbsTolerance(1e-4);
        ROOT::Math::IntegratorOneDimOptions::SetDefaultRelTolerance(1e-2);
        ROOT::Math::IntegratorOneDimOptions::SetDefaultIntegrator("Gauss");

        // Initialize |F(t)|^2 wedge cut TF1
        cutFFt2 = new TF1("Cut Form Factor: |F(t)|^{2}", [this, phi_min, phi_max] (double *var, double *par)
        {
            double t = var[0];  // t is variable
            TF1 fft2("fft2", calc_FormFactor_tWedge, phi_min, phi_max, 5);
            fft2.SetParameters(par[0], par[1], par[2], par[3], t);
            double integral = fft2.Integral(phi_min, phi_max, 1e-12);
            return integral;
        }, t_min, t_max, 4);
        cutFFt2->SetParameters(A,Vo,R,a0);

        // Initialize histogram for TF1
        hist1D = new TH1D("FF", "Form Factor", bins, t_min, t_max);
        hist1D->Sumw2();
        double step = t_max/bins;
        for(int i=0;i<bins;i++) 
        {
            double t = (i+1)*step;
            double q = sqrt(t);
            double val = wedge_FF2_inte(q, phi_min, phi_max);
            hist1D->SetBinContent(i+1, val); 
        }
}


    // Calculations for form factor
    static double calc_FormFactor_tWedge(double *var, double *par)
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
    static double calc_FormFactor_qWedge(double q)
    {
        double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7;
	    const double arg1 = q*R / hbarc;  
	    const double arg2 = hbarc / q;
	    const double arg3  = (sin(arg1)-arg1*cos(arg1)) * 4*pi*Vo * arg2*arg2*arg2 / A; 
	    double result = 1/(2*pi)*arg3 / (1. + (a0*a0 * q*q/(hbarc*hbarc))); 
        return result*result; // [-]
    }
    // Integrand to calculate |F(q)|^2 for wedge cut
    static double wedge_calcFF2(double *var, double *par)
    {
        double qq = par[0];
        double phi = var[0];  // integration variable
        // put q in terms of theta
        double qx = qq*sin(phi), qy = qq*cos(phi), q = sqrt(qx*qx+qy*qy); 
	    return calc_FormFactor_qWedge(q);  // Integrand
    }
    // Integration over phi
    double wedge_FF2_inte(double q, double phi_min, double phi_max)
    {
        TF1 *f = new TF1("f",wedge_calcFF2,phi_min,phi_max,1);  
        f->SetParameters(q);
        return f->Integral(phi_min,phi_max,1e-12);
    }

    // Functions
    TF1 *getCutFormFactor_t2_1D() const {return cutFFt2;}
    TH1D *getCutHist_1D() const {return hist1D;}

private:
    double A, Vo, R, a0;
    TF1 *cutFFt2;
    TH1D *hist1D;
};

#endif