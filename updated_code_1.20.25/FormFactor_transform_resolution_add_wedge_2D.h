#ifndef FormFactor_transform_resolution_add_wedge_2D_H
#define FormFactor_transform_resolution_add_wedge_2D_H

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
class FormFactor_transform_resolution_add_wedge_2D
{
public:
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
    double qx_min = 0;
    double qx_max = 20;
    TF1 *f = new TF1("",guassian,qx_min,qx_max,3);  
    f->SetParameters(qy,qx_prime,sigma);
    return f->Integral(qx_min,qx_max,1e-12); // F(qx',qy) [GeV]
}

// Integrand for wedge cut
static double testwedge_calcFF2(double *x, double *par)
{
    double q_prime = par[0], sigma = par[1];
    // integration variable
    double phi = x[0];
    double qx_prime = q_prime*sin(phi);
    double qy = q_prime*cos(phi);
    // call form factor with resolution
	return FF_wRes(qy,qx_prime,sigma)*FF_wRes(qy,qx_prime,sigma); // [GeV^2]
}

// Integral to integrate over theta
static double testwedge_FF2_inte_2D(double q_prime, double phi_min, double phi_max, double sigma)
{
    TF1 *f = new TF1("",testwedge_calcFF2,phi_min,phi_max,2);  
    f->SetParameters(q_prime,sigma);
    return f->Integral(phi_min,phi_max,1e-12); //[GeV]
}




    // Integrand: |F(t)|^2 (wedge, no resolution)->G(x,y) 
static double G_integrand(double *x, double *par) 
{
    double xx = par[0];
    double y = par[1],  phi_min = par[2], phi_max = par[3], sigma = par[4];
    double r = sqrt(xx*xx + y*y);
    // integration variable
    double qx_prime = x[0], qy = x[1];
    double q_prime = sqrt(qx_prime*qx_prime+qy*qy);
    return sqrt(2*pi/(r))*sqrt(2*hbarc/(pi*q_prime*r))*sin(q_prime*r/hbarc)*testwedge_FF2_inte_2D(q_prime,phi_min,phi_max,sigma)*sqrt(q_prime)*1/hbarc*1/hbarc*1/sqrt(hbarc); // [GeV^-1 fm^-3]
}

// Integral for G(x,y)
static double G_integral_dq(double x, double y,double phi_min, double phi_max, double sigma)
{
    double qx_prime_min = 0., qx_prime_max = 1., qy_min = 0., qy_max = 1.;
    TF2 *f = new TF2("",G_integrand,qx_prime_min,qx_prime_max,qy_min,qy_max,5);  
    f->SetParameters(x,y,phi_min,phi_max,sigma);
    return f->Integral(qx_prime_min,qx_prime_max,qy_min,qy_max,1e-12); //[1/fm^3]
}


       
        TH2D *hist2(const char *name, const char *title, double bins, double xmin, double xmax, double ymin, double ymax, double phi_min, double phi_max, double sigma)
        {
            TH2D *hist2D = new TH2D(name, title, bins, xmin, xmax, bins, ymin, ymax);
            hist2D->Sumw2();
            double step = xmax/bins;
            for(int i=0;i<bins;i++)
            {
                double x = (i+1)*step;
                for(int j=0;j<bins;j++)
                {
                    double y = (j+1)*step;
                    double val = G_integral_dq(x, y, phi_min, phi_max, sigma);
                    hist2D->SetBinContent(i+1, j+1, val);
                }
                cout << "x: " << x << "\n";
            }
            return hist2D;
        }


    
};

#endif

