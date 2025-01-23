#ifndef FormFactor_transform_resolution_add_wedge_1D_H
#define FormFactor_transform_resolution_add_wedge_1D_H

#include "EICvalueconst.h"

using namespace std;




/********************************************** 
 
 * Form Factor with detector resolution and wedge cut 1D transform to position space (Woods-Saxon density distribution): |F(t)|^2 -> G(r)
 
 * Parameters:
    Vo: Nuclear density (hard sphere approximation) [fm^-3]
    A: Atomic mass number
    R: Radius of nucleus [fm]
    a0: Range of Yukawa potential [fm]
 * Variables:
    q: Transverse momentum >> q = sqrt(qx^2+qy^2) [GeV]
    t: Four-momentum transfer squared >> t = tx+ty = q^2 [GeV^2]
    r: radial distance from center of nucleus >> r = sqrt(x^2+y^2) [fm]

**********************************************/
class FormFactor_transform_resolution_add_wedge_1D
{
public:
/*
    FormFactor_transform_resolution_add_wedge_1D(double A_init, double Vo_init, double R_init, double a0_init, double bins, 
        double phi_min, double phi_max, double sigma, double r_min, double r_max): A(A_init), Vo(Vo_init), R(R_init), a0(a0_init) 
    {
        ROOT::Math::IntegratorOneDimOptions::SetDefaultAbsTolerance(1e-4);
        ROOT::Math::IntegratorOneDimOptions::SetDefaultRelTolerance(1e-2);
        ROOT::Math::IntegratorOneDimOptions::SetDefaultIntegrator("Gauss");

        // Initialize TF1
        double qq_min = 0, qq_max = 1;
        FF = new TF1("", [qq_min, qq_max, phi_min, phi_max, sigma] (double *var, double *par)
        {
            double r = var[0]; 
            TF1 fft2("fft2", G_integrand_1D, qq_min, qq_max, 4);
            fft2.SetParameters(r, phi_min, phi_max, sigma);
            double integral = fft2.Integral(qq_min, qq_max, 1e-12);
            return integral;
        }, r_min, r_max, 0);

        // Initialize histogram 
        hist1D = new TH1D("hist", "", bins, r_min, r_max);
        hist1D->Sumw2();
        double step = r_max/bins;
        for(int i=0;i<bins;i++)
        {
            double r = (i+1)*step;
            double val = G_integral_dq_1D(r, phi_min, phi_max, sigma);
            hist1D->SetBinContent(i+1, val);
        }
    }
*/

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
    static double guassian(double *x, double *par)
    {
        double qy = par[0], qx_prime = par[1], sigma = par[2];
        double qx = x[0];  // integration variable
        double q = sqrt(qx*qx+qy*qy);
	    return calcFF(q)*exp(-(qx-qx_prime)*(qx-qx_prime)/(2*sigma*sigma));
    }
    // Integral over qx, new 2D form factor with smearing
    static double FF_wRes(double qy, double qx_prime, double sigma)
    {
        double qx_min = 0, qx_max = 5;
        TF1 *f = new TF1("f", guassian, qx_min, qx_max, 3);  
        f->SetParameters(qy, qx_prime, sigma);
        return f->Integral(qx_min, qx_max, 1e-12); // F(qx',qy) [GeV]
    }
    // Integrand for wedge cut
    static double testwedge_calcFF2_2D(double *x, double *par)
    {
        double q = par[0], sigma = par[1];
        double phi = x[0];  // integration variable
        double qx_prime = q*sin(phi);
        double qy = q*cos(phi);
	    return FF_wRes(qy, qx_prime, sigma); // [GeV]
    }
    // Integral to integrate over phi
    static double testwedge_FF2_inte(double q, double phi_min, double phi_max, double sigma)
    {
        TF1 *f = new TF1("f", testwedge_calcFF2_2D, phi_min, phi_max, 2);  
        f->SetParameters(q, sigma);
        return f->Integral(phi_min, phi_max, 1e-12); //[GeV]
    }
    // Integrand: |F(q)|^2 ->G (r) 
    static double G_integrand_1D(double *x, double *par) 
    {
        double r = par[0], phi_min = par[1], phi_max = par[2], sigma = par[3];
        double q = x[0];  // integration variable
        return sqrt(2*pi/(r))*sqrt(2*hbarc/(pi*q*r))*sin(q*r/hbarc)*testwedge_FF2_inte(q,phi_min,phi_max,sigma)*sqrt(q)*1/hbarc*1/hbarc*1/sqrt(hbarc); // [GeV^-1 fm^-3]
    }
    // Integral for G(r)
    static double G_integral_dq_1D(double r,double phi_min, double phi_max, double sigma)
    {
        double qq_min = 0., qq_max = 5.;
        TF1 *f = new TF1("f", G_integrand_1D, qq_min, qq_max, 4);  
        f->SetParameters(r, phi_min, phi_max, sigma);
        return f->Integral(qq_min, qq_max, 1e-12); //[1/fm^3]
    }


    TF1 *transform_FF_to_G_fun(double phi_min, double phi_max, double sigma, double r_min, double r_max)
    {
        double qq_min = 0, qq_max = 1;
        TF1 *FF = new TF1("", [qq_min, qq_max, phi_min, phi_max, sigma] (double *var, double *par)
        {
            double r = var[0]; 
            TF1 fft2("fft2", G_integrand_1D, qq_min, qq_max, 4);
            fft2.SetParameters(r, phi_min, phi_max, sigma);
            double integral = fft2.Integral(qq_min, qq_max, 1e-12);
            return integral;
        }, r_min, r_max, 0);
        return FF;
    }

    TH1D *transform_FF_to_G_hist(double phi_min, double phi_max, double sigma, double r_min, double r_max, double bins)
    {
        TH1D *hist1D = new TH1D("hist", "", bins, r_min, r_max);
        hist1D->Sumw2();
        double step = r_max/bins;
        for(int i=0;i<bins;i++)
        {
            double r = (i+1)*step;
            double val = G_integral_dq_1D(r, phi_min, phi_max, sigma);
            hist1D->SetBinContent(i+1, val);
        }
        return hist1D;
    }



/*
    // Functions
    TF1 *getWedgeResTransform_fun_1D() const {return FF;}
    TH1D *getWedgeResTransform_hist_1D() const {return hist1D;}


private:

    double A, Vo, R, a0;
    TF1 *FF;
    TH1D *hist1D;

    */
};

#endif






