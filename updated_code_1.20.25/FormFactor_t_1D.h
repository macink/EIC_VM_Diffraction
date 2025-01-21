#ifndef FormFactor_t_1D_H
#define FormFactor_t_1D_H

#include "EICvalueconst.h"

using namespace std;

/********************************************** 
 
 * Form Factor 1D: F(t) [-]
 
 * Parameters:
    Vo: Nuclear density (hard sphere approximation) [fm^-3]
    A: Atomic mass number
    R: Radius of nucleus [fm]
    a0: Range of Yukawa potential [fm]
 * Variables:
    t: Four-momentum transfer squared >> t = tx+ty = q^2 [GeV^2]

**********************************************/
class FormFactor_t_1D
{
public:
    FormFactor_t_1D(double A_init, double Vo_init, double R_init, double a0_init,
        double t_min, double t_max, double bins, double r_min, double r_max): A(A_init), Vo(Vo_init), R(R_init), a0(a0_init)
    {
        ROOT::Math::IntegratorOneDimOptions::SetDefaultAbsTolerance(1e-4);
        ROOT::Math::IntegratorOneDimOptions::SetDefaultRelTolerance(1e-2);
        ROOT::Math::IntegratorOneDimOptions::SetDefaultIntegrator("Gauss");

        // Initialize |F(t)|^2 TF1
        FFt1 = new TF1("Form Factor: F(t)", calc_FormFactor_t1D, t_min, t_max, 4);
        FFt1->SetParameters(A, Vo, R, a0);

        // Initialize hist
        FF_hist = new TH1D("FF_hist","Form Factor", bins, t_min, t_max);
        FF_hist->Sumw2();
        double step2 = t_max/bins;
        for(int i=0;i<bins;i++)
        {
            double t = (i+1)*step2;
            double q = sqrt(t);
            double val = calcFF_1D(q);
            FF_hist->SetBinContent(i+1, val);
        }
        
        // Initialize TF1 for transform: |F(t)|^2 -> G(r)
        double qq_min = 0, qq_max = 5;
        transform_TF1 = new TF1("Fourier-Bessel Transformation: |F(t)|^{2} -> G(r)", [this, qq_min, qq_max] (double *var, double *par)
        {
            double r = var[0];
            TF1 trans_integral("trans_integral", trueG1_integrand_dq_1D, qq_min, qq_max, 1);
            trans_integral.SetParameters(r);
            return trans_integral.Integral(qq_min, qq_max, 1e-12);
        }, r_min, r_max, 0); 

       // Initialize hist for transform: |F(t)|^2 -> G(r)
        transform_hist = new TH1D("transform_hist","Transform Form Factor", bins, r_min, r_max);
        transform_hist->Sumw2();
        double step1 = r_max/bins;
        for(int i=0;i<bins;i++)
        {
            double r = (i+1)*step1;
            double val = trueG1_integral_dq_1D(r);
            transform_hist->SetBinContent(i+1, val);
        }
    }


    // Calculate Form Factor and G integrand
    static double calc_FormFactor_t1D(double *var, double *par) // F(t)
    {
        double A = par[0], Vo = par[1], R = par[2], a0 = par[3];
        double t = var[0];
        double q = sqrt(t);
        const double arg1 = q*R / hbarc;  
        const double arg2 = hbarc / q;
        const double arg3 = (sin(arg1)-arg1*cos(arg1)) * 4*pi*Vo * arg2*arg2*arg2 / A; 
        double result = arg3 / (1. + (a0*a0 * q*q / (hbarc * hbarc)));  
        return result*result; // [-]
    }
    static double calcFF_1D(double q)
    {
	    if(q==0){return 0;}
        else{
            double A = 197, Vo = 2.12, R = 6.38, a0 = .70;
	        const double arg1 = q*R / hbarc;  
	        const double arg2 = hbarc /q;
	        const double sph  = (sin(arg1)-arg1*cos(arg1)) * 4*pi*Vo * arg2*arg2*arg2 / A; 
            double result = sph * 1/(1. + (a0*a0 * q*q/(hbarc*hbarc)));
	        return result*result;  // [-]
        }
    }
    // Calculations for transformation |F(t)|^2 -> G(r)
    static double calcFF_1D_transform(double q)
    {
	    if(q==0){return 0;}
        else{
            double A = 197, Vo = 2.12, R = 6.38, a0 = .70;
	        const double arg1 = q*R / hbarc;  
	        const double arg2 = hbarc /q;
	        const double sph  = (sin(arg1)-arg1*cos(arg1)) * 4*pi*Vo * arg2*arg2*arg2 / A; 
            double result = sph * 1/(1. + (a0*a0 * q*q/(hbarc*hbarc)));
	        return result;  // [-]
        }
    }
    static double trueG1_integrand_dq_1D(double *var, double *par) 
    {
        double r = par[0];
        double q = var[0];
        return 2*pi*sqrt(2*pi/r)*sqrt(2*hbarc/(pi*q*r))*sin(q*r/hbarc)*calcFF_1D_transform(q)*sqrt(q)*q*1/hbarc*1/hbarc*1/sqrt(hbarc); //[1/fm^3GeV^2] 
    }   
    static double trueG1_integral_dq_1D(double r)
    {
        double qq_min = 0;
        double qq_max = 5;
        TF1 *f = new TF1("f", trueG1_integrand_dq_1D, qq_min, qq_max, 1);  
        f->SetParameters(r);
        return f->Integral(qq_min, qq_max, 1e-12); //[1/fm^3]
    }


    // Functions
    TF1 *getFormFactort_1D() const {return FFt1;}
    TF1 *getTransformedFF_TF1() const {return transform_TF1;}
    TH1D *getTransformedFF_hist() const {return transform_hist;}
    TH1D *getFF_hist() const {return FF_hist;}

private:
    double A, Vo, R, a0;
    TF1 *FFt1, *transform_TF1;
    TH1D *transform_hist, *FF_hist;
};


#endif
