#ifndef FF_WS_Transforms_H
#define FF_WS_Transforms_H

#include "EICvalueconst.h"

using namespace std;

/********************************************** 
 
 * Form Factor 1D: F(q) and F(t) [-]
 
 * Parameters:
    Vo: Nuclear density (hard sphere approximation) [fm^-3]
    A: Atomic mass number
    R: Radius of nucleus [fm]
    a0: Range of Yukawa potential [fm]
 * Variables:
    q: Transverse momentum >> q = sqrt(qx^2+qy^2) [GeV]
    t: Four-momentum transfer squared >> t = tx+ty = q^2 [GeV^2]

**********************************************/
class FF_WS_Transforms
{
public:
    FF_WS_Transforms(double A_init, double Vo_init, double R_init, double a0_init, double q_min, double q_max, double t_min, double t_max, 
        double phi_min, double phi_max, double bins, double r_min, double r_max): A(A_init), Vo(Vo_init), R(R_init), a0(a0_init)
    {
        ROOT::Math::IntegratorOneDimOptions::SetDefaultAbsTolerance(1e-4);
        ROOT::Math::IntegratorOneDimOptions::SetDefaultRelTolerance(1e-2);
        ROOT::Math::IntegratorOneDimOptions::SetDefaultIntegrator("Gauss");

       
        // Initialize TF1 for transform: |F(t)|^2 -> G(r)
        transform1 = new TF1("Fourier-Bessel Transformation: |F(t)|^{2} -> G(r)",[this,q_min,q_max](double *var, double *par){double r = var[0];
            TF1 trans_integral("trans_integral",trueG1_integrand_dq_1D,q_min,q_max,1);
            trans_integral.SetParameters(r);
            return trans_integral.Integral(q_min,q_max,1e-12);
        },r_min,r_max,0);  

        // Initialize hist for transform: |F(t)|^2 -> G(r)
        hankelTransformFF = new TH1D("hankelTransformFF","Transform Form Factor",bins,r_min,r_max);
        hankelTransformFF->Sumw2();
        double step = r_max/bins;
        for(int i=0;i<bins;i++)
        {
            double r = (i+1)*step;
            double val = trueG1_integral_dq_1D(r);
            hankelTransformFF->SetBinContent(i+1,val);
            //cout << "r: " << r << endl;
        }
        hankelTransformFF->SetTitle("Transformed WS: |F(t)|^{2} #rightarrow G(r)");
        hankelTransformFF->GetYaxis()->SetTitle("G(r) [fm^{-3}]");
	    hankelTransformFF->GetXaxis()->SetTitle("r [fm]");    
    }


    static double calcFF_1D(double q)
    {
	    if(q==0){return 0;}
        else{
            double A = 197, ro = 1.25, rho0 = 2.12, R = 6.38, a0 = .70;
	        const double arg1 = q * R / hbarc;  
	        const double arg2 = hbarc /q;
	        const double sph  = (sin(arg1) - arg1 * cos(arg1)) * 4*pi*rho0 * arg2 * arg2 *arg2/ double(A); 
	        return sph * 1/(1. + (a0 * a0 * q*q/(hbarc*hbarc)));  // [-]
        }
    }
    static double trueG1_integrand_dq_1D(double *x, double *par) 
    {
        double r = par[0];
        double q = x[0];
        //double q = sqrt(t);
        return 2*pi*sqrt(2*pi/r)*sqrt(2*hbarc/(pi*q*r))*sin(q*r/hbarc)*calcFF_1D(q)*sqrt(q)*q*1/hbarc*1/hbarc*1/sqrt(hbarc); //[1/fm^3GeV^2] 
    }
    static double trueG1_integral_dq_1D(double r)
    {
        double q_min = 0.000005;
        double q_max = 5;
        TF1 *f = new TF1("f",trueG1_integrand_dq_1D,q_min,q_max,1);  
        f->SetParameters(r);
        return f->Integral(q_min,q_max,1e-6); //[1/fm^3]
    }


// Method to transform a given TF1 object
TF1* transformTF1(TF1* inputTF1, double r_min, double r_max, double q_min, double q_max)
{
TF1 *transform1 = new TF1("Fourier-Bessel Transformation: |F(t)|^{2} -> G(r)",[this,q_min,q_max](double *var, double *par){double r = var[0];
            TF1 trans_integral("trans_integral",trueG1_integrand_dq_1D,q_min,q_max,1);
            trans_integral.SetParameters(r);
            return trans_integral.Integral(q_min,q_max,1e-12);
        },r_min,r_max,0);  

return transform1;
}

    

    // Functions
    TF1 *getFFtransform() const {return transform1;}
    TH1D *getTransformFFhist() const {return hankelTransformFF;}

private:
    double A, Vo, R, a0;
    TF1 *transform1;
    vector<double> x_vals_vec;
    vector<double> y_vals_vec;
    TH1D *hankelTransformFF;
};


#endif
