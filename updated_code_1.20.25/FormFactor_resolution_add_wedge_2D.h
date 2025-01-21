#ifndef FormFactor_resolution_add_wedge_2D_H
#define FormFactor_resolution_add_wedge_2D_H

#include "EICvalueconst.h"

using namespace std;




/********************************************** 
 
 * Form Factor with detector resolution and wedge cut 2D: |F(t_x,t_y)|^2
 
 * Parameters:
    Vo: Nuclear density (hard sphere approximation) [fm^-3]
    A: Atomic mass number
    R: Radius of nucleus [fm]
    a0: Range of Yukawa potential [fm]
 * Variables:
    q: Transverse momentum >> q = sqrt(qx^2+qy^2) [GeV]
    t: Four-momentum transfer squared >> t = tx+ty = q^2 [GeV^2]

**********************************************/
class FormFactor_resolution_add_wedge_2D
{
public:
    FormFactor_resolution_add_wedge_2D(double A_init, double Vo_init, double R_init, double a0_init, double ty_min, double ty_max, 
        double tx_prime_min, double tx_prime_max, double bins, double phi_min, double phi_max, double sigma, double x_min, 
        double x_max, double y_min, double y_max): A(A_init), Vo(Vo_init), R(R_init), a0(a0_init) 
    {
        ROOT::Math::IntegratorOneDimOptions::SetDefaultAbsTolerance(1e-4);
        ROOT::Math::IntegratorOneDimOptions::SetDefaultRelTolerance(1e-2);
        ROOT::Math::IntegratorOneDimOptions::SetDefaultIntegrator("Gauss");

        // Initialize |F(tx,ty)|^2 wedge cut TF2
        cutFFt2d = new TF2("Cut Form Factor: |F(t_{x},t_{y})|^{2}", [this, phi_min, phi_max, sigma] (double *var, double *par)
        {
            double tx_prime = var[0], ty = var[1]; 
            TF1 fft2("fft2", FF_cut_wRes_2D, phi_min, phi_max, 3);
            fft2.SetParameters(tx_prime, ty, sigma);
            double integral = fft2.Integral(phi_min, phi_max, 1e-12);
            return integral;
        }, tx_prime_min, tx_prime_max, ty_min, ty_max, 0);

        // Initialize histogram for TF2
        hist2D = new TH2D("FF", "Form Factor", bins, tx_prime_min, tx_prime_max, bins, ty_min, ty_max);
        hist2D->Sumw2();
        double step2 = tx_prime_max/bins;
        for(int i=0;i<bins;i++) 
        {
            //double tx_prime = (i+1)*step2;
            for(int j=0;j<bins;j++) 
            {
                //double ty = (j+1)*step2;
                //double val = wedge_FF_2D(tx_prime, ty, sigma, phi_min, phi_max);
                //hist2D->SetBinContent(i+1, j+1, val); 
                double tx_prime = hist2D->GetXaxis()->GetBinCenter(i+1);
                double ty = hist2D->GetYaxis()->GetBinCenter(j+1);
                hist2D->SetBinContent(i+1, j+1, cutFFt2d->Eval(tx_prime, ty));
            }
        }
    }


    // Calculations to add resolution to form factor
    static double calcFF_2D(double qx, double qy)
    {
        double q = sqrt(qx*qx+qy*qy);
	    if(q==0){return 0;}
        else
        {
            double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7;
	        const double arg1 = q*R / hbarc;  
	        const double arg2 = hbarc /q;
	        const double sph  = (sin(arg1)-arg1*cos(arg1)) * 4*pi*Vo * arg2*arg2*arg2 / A; 
	        double result = sph * 1/(1. + (a0*a0 * q*q/(hbarc*hbarc)));  // [-]
            return result*result;
        }
    }
    // Integrand to calcuate smearing
    static double guassian_2D(double *x, double *par)
    {
        double sigma = par[0], qx_prime = par[1], qy = par[2];
        //double qy = sqrt(ty), qx_prime = sqrt(tx_prime);
        double qx = x[0];  // integration variable
	    return calcFF_2D(qx,qy)*exp(-(qx-qx_prime)*(qx-qx_prime)/(2*sigma*sigma));
    }
    // Integral over qx -> new form factor with smearing
    static double FF_wRes_2D(double qy, double qx_prime, double sigma) 
    {
        double qqx_min = 0, qqx_max = 5;
        TF1* f = new TF1("f", guassian_2D, qqx_min, qqx_max, 3);  
        f->SetParameters(sigma, qx_prime, qy);
        return f->Integral(qqx_min, qqx_max, 1e-12); // F(tx',ty) [GeV]
    }
    // Integrand for wedge cut
    static double FF_cut_wRes_2D(double *x, double *par)
    {
        double tx_prime = par[0], ty = par[1], sigma = par[2];
        double phi = x[0];  // integration variable
        double q = sqrt(tx_prime+ty), qx_prime = q*sin(phi), qy = q*cos(phi);
        // call form factor with resolution
	    return FF_wRes_2D(qy,qx_prime,sigma); // [GeV^2]
    }
    // Integral to integrate over theta
    static double wedge_FF_2D(double tx_prime, double ty, double sigma, double phi_min, double phi_max)
    {
        double qx_prime = sqrt(tx_prime), qy = sqrt(qy);
        TF1 *f = new TF1("f",FF_cut_wRes_2D,phi_min,phi_max,3);  
        f->SetParameters(qx_prime, qy, sigma);
        return f->Integral(phi_min,phi_max,1e-12); //[GeV]
    }


    // Functions
    TF2 *getWedgeRes_fun_2D() const {return cutFFt2d;}
    TH2D *getWedgeRes_hist_2D() const {return hist2D;}


private:
    double A, Vo, R, a0;
    TF2 *cutFFt2d;
    TH2D *hist2D;
};

#endif

