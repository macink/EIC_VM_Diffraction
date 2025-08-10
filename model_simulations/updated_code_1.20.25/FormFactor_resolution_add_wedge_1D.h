#ifndef FormFactor_resolution_add_wedge_1D_H
#define FormFactor_resolution_add_wedge_1D_H

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
class FormFactor_resolution_add_wedge_1D
{
public:
    FormFactor_resolution_add_wedge_1D(double A_init, double Vo_init, double R_init, double a0_init, double t_min, double t_max, double bins, 
        double phi_min, double phi_max, double sigma, double r_min, double r_max): A(A_init), Vo(Vo_init), R(R_init), a0(a0_init) 
    {
        ROOT::Math::IntegratorOneDimOptions::SetDefaultAbsTolerance(1e-4);
        ROOT::Math::IntegratorOneDimOptions::SetDefaultRelTolerance(1e-2);
        ROOT::Math::IntegratorOneDimOptions::SetDefaultIntegrator("Gauss");

        // Initialize |F(t)|^2 wedge cut TF1
        cutFFt2 = new TF1("Cut Form Factor: |F(t)|^{2}", [this, phi_min, phi_max, sigma] (double *var, double *par)
        {
            double t = var[0];  // t is variable
            TF1 fft2("", FF_cut_wRes, phi_min, phi_max, 2);
            fft2.SetParameters(t, sigma);
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


    // Calculations to add resolution to form factor
    static double calcFF(double q)
    {
        double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7;
	    const double arg1 = q*R / hbarc;  
	    const double arg2 = hbarc / q;
	    const double arg3  = (sin(arg1)-arg1*cos(arg1)) * 4*pi*Vo * arg2*arg2*arg2 / A; 
	    double result = 1/(2*pi)*arg3 / (1. + (a0*a0 * q*q/(hbarc*hbarc))); 
        return result*result; // [-]
    }
    // Integrand to calcuate smearing
    static double guassian(double *x, double *par)
    {
        double qy = par[0], qx_prime = par[1], sigma = par[2];
        //double qy = sqrt(ty), qx_prime = sqrt(tx_prime);
        double qx = x[0];  // integration variable
        double q = sqrt(qx*qx+qy*qy);
	    return calcFF(q)*exp(-(qx-qx_prime)*(qx-qx_prime)/(2*sigma*sigma));
    }
    // Integral over qx -> new form factor with smearing
    static double FF_wRes(double qy, double qx_prime, double sigma) 
    {
        double qqx_min = 0, qqx_max = 20;
        TF1 *f = new TF1("f", guassian, qqx_min, qqx_max, 3);  
        f->SetParameters(qy, qx_prime, sigma);
        return f->Integral(qqx_min, qqx_max, 1e-12); // F(qx',qy) [GeV]
    }
    // Integrand for wedge cut
    static double FF_cut_wRes(double *x, double *par)
    {
        double t = par[0], sigma = par[1];
        double phi = x[0];  // integration variable
        double q = sqrt(t), qx_prime = q*sin(phi), qy = q*cos(phi);
        // call form factor with resolution
	    return FF_wRes(qy, qx_prime, sigma); // [GeV^2]
    }
    // Integral to integrate over theta
    static double wedge_FF(double t, double sigma, double phi_min, double phi_max)
    {
        double q = sqrt(t);
        TF1 *f = new TF1("", FF_cut_wRes, phi_min, phi_max, 2);  
        f->SetParameters(q, sigma);
        return f->Integral(phi_min, phi_max, 1e-12); //[GeV]
    }
   

    // Functions
    TF1 *getWedgeRes_fun_1D() const {return cutFFt2;}
    TH1D *getWedgeRes_hist_1D() const {return hist1D;}

private:
    double A, Vo, R, a0;
    TF1 *cutFFt2;
    TH1D *hist1D;
};

#endif






