#ifndef FormFactor_q_2DwResCut_H
#define FormFactor_q_2DwResCut_H

#include "EICvalueconst.h"

using namespace std;

/********************************************** 
 
 * Form Factor 2D: F(qx,qy) [-]
 
 * Parameters:
    Vo: Nuclear density (hard sphere approximation) [fm^-3]
    A: Atomic mass number
    R: Radius of nucleus [fm]
    a0: Range of Yukawa potential [fm]
 * Variables:
    q: Transverse momentum >> q = sqrt(qx^2+qy^2) [GeV]

**********************************************/
class FormFactor_q_2DwResCut
{
public:
    FormFactor_q_2DwResCut(double A_init, double Vo_init, double R_init, double a0_init, double qy_min, double qy_max, 
        double qx_prime_min, double qx_prime_max, double bins, double sigma): A(A_init), Vo(Vo_init), R(R_init), a0(a0_init)
    {
        ROOT::Math::IntegratorOneDimOptions::SetDefaultAbsTolerance(1e-4);
        ROOT::Math::IntegratorOneDimOptions::SetDefaultRelTolerance(1e-2);
        ROOT::Math::IntegratorOneDimOptions::SetDefaultIntegrator("Gauss");

        // Initialize |F(tx,ty)|^2 with resolution TF2
        FFsmear_t = new TF2("", FF_wRes_2D, -qx_prime_max, qx_prime_max, -qy_max, qy_max, 1);
        double params[] = {sigma};
        FFsmear_t->SetParameters(params);

        // Initialize 2D hist
        FF = new TH2D("","", bins, qx_prime_min, qx_prime_max, bins, qy_min, qy_max);
        FF->Sumw2();
        double step = (qx_prime_max-qx_prime_min)/bins;
        for(int i=0;i<bins;i++)
        {
            double qx_prime = qx_prime_min+(i+0.5)*step;
            //double qx_prime = (i+1)*step;
            double tx_prime = qx_prime*qx_prime;
            for(int j=0;j<bins;j++)
            {
                double qy = qy_min+(j+0.5)*step;
                //double qy = (j+1)*step;
                double ty = qy*qy;
                double var[] = {tx_prime, ty};
                double par[] = {sigma};
                double val = FF_wRes_2D(var, par);
                FF->SetBinContent(i+1, j+1, val);
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
        double sigma = par[0], ty = par[2], tx_prime = par[1];
        double qy = sqrt(ty), qx_prime = sqrt(tx_prime);
        double qx = x[0];  // integration variable
	    return 1/sqrt(2*pi)*1/sigma*calcFF_2D(qx,qy)*exp(-(qx-qx_prime)*(qx-qx_prime)/(2*sigma*sigma));
    }
    // Integral over qx, new 2D form factor with smearing
    static double FF_wRes_2D(double *var, double *par) 
    {
        double nBins = 5000;
        double tx_prime = var[0], ty = var[1];
        double sigma = par[0];
        double qqx_min = -5, qqx_max = 5;
        double pars[] = {sigma, tx_prime, ty};
        double step = (qqx_max - qqx_min) / nBins;
        double integral = 0;
        for (int i = 0; i < nBins; i++) {
            double qx = qqx_min + (i + 0.5) * step;
            double x[] = { qx };
            integral += guassian_2D(x, pars) * step;
        }
        return integral; 
        //TF2* f = new TF2("", guassian_2D, qqx_min, qqx_max, qqx_min, qqx_max, 3);  
        //f->SetParameters(sigma, tx_prime, ty);
        //return f->Integral(qqx_min, qqx_max,qqx_min,qqx_max,1e-15); // F(tx',ty) [GeV] 
    }
    

    // Functions
    TF2 *getSmearedFormFactor() const {return FFsmear_t;}
    TH2D *getSmeared_hist() const {return FF;}


private:
    double A, Vo, R, a0;
    TF2 *FFsmear_t;
    TH2D *FF;
};
/*class FormFactor_q_2DwResCut
{
public:
    FormFactor_q_2DwResCut(double A_init, double Vo_init, double R_init, double a0_init, double qy_min, double qy_max, 
        double qx_prime_min, double qx_prime_max, double bins, double sigma): A(A_init), Vo(Vo_init), R(R_init), a0(a0_init)
    {
        ROOT::Math::IntegratorOneDimOptions::SetDefaultAbsTolerance(1e-4);
        ROOT::Math::IntegratorOneDimOptions::SetDefaultRelTolerance(1e-2);
        ROOT::Math::IntegratorOneDimOptions::SetDefaultIntegrator("Gauss");

        // Initialize |F(tx,ty)|^2 with resolution TF2
        FFsmear_t = new TF2("", FF_wRes_2D, qx_prime_min, qx_prime_max, qy_min, qy_max, 1);
        FFsmear_t->SetParameters(sigma);

        // Initialize 2D hist
        FF = new TH2D("","", bins, -qx_prime_max, qx_prime_max, bins, -qy_max, qy_max);
        FF->Sumw2();
        double step = qx_prime_max*2/bins;
        for(int i=0;i<bins;i++)
        {
            double qx_prime = -qx_prime_max+(i+0.5)*step;
            double tx_prime = qx_prime*qx_prime;
            for(int j=0;j<bins;j++)
            {
                double qy = -qy_max+(j+0.5)*step;
                double ty = qy*qy;
                double var[] = {tx_prime, ty};
                double par[] = {sigma};
                double val = FF_wRes_2D(var, par);
                FF->SetBinContent(i+1, j+1, val);
            }
        }
    }


    // Calculations to add resolution to form factor
    static double calcFF_2D(double qx, double qy)
    {
        double q = sqrt(qx*qx+qy*qy);
	    if(q==0){return 1e-10;}
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
        double sigma = par[0], ty = par[2], tx_prime = par[1];
        double qy = sqrt(ty), qx_prime = sqrt(tx_prime);
        double qx = x[0];  // integration variable
	    return calcFF_2D(qx,qy)*exp(-(qx-qx_prime)*(qx-qx_prime)/(2*sigma*sigma));
    }
    // Integral over qx, new 2D form factor with smearing
    static double FF_wRes_2D(double *var, double *par) 
    {
        double tx_prime = var[0], ty = var[1];
        double sigma = par[0];
        double qqx_min = 0, qqx_max = 5;
        TF1* f = new TF1("f", guassian_2D, qqx_min, qqx_max, 3);  
        f->SetParameters(sigma, tx_prime, ty);
        return f->Integral(qqx_min, qqx_max, 1e-9); // F(tx',ty) [GeV]
    }
    

    // Functions
    TF2 *getSmearedFormFactor() const {return FFsmear_t;}
    TH2D *getSmeared_hist() const {return FF;}


private:
    double A, Vo, R, a0;
    TF2 *FFsmear_t;
    TH2D *FF;
};*/

#endif
