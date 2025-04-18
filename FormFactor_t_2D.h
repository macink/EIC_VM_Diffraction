#ifndef FormFactor_t_2D_H
#define FormFactor_t_2D_H

#include "EICvalueconst.h"

using namespace std;

/********************************************** 
 
 * Form Factor 2D: F(tx,ty) [-]
 
 * Parameters:
    Vo: Nuclear density (hard sphere approximation) [fm^-3]
    A: Atomic mass number
    R: Radius of nucleus [fm]
    a0: Range of Yukawa potential [fm]
 * Variables:
    t: Four-momentum transfer squared >> t = tx+ty = q^2 [GeV^2]

**********************************************/
class FormFactor_t_2D
{
public:
    FormFactor_t_2D(double A_init, double Vo_init, double R_init, double a0_init, double tx_min, double tx_max,
        double ty_min, double ty_max, double bins, double x_min, double x_max, double y_min, double y_max): A(A_init), Vo(Vo_init), R(R_init), a0(a0_init)
    {
        ROOT::Math::IntegratorOneDimOptions::SetDefaultAbsTolerance(1e-4);
        ROOT::Math::IntegratorOneDimOptions::SetDefaultRelTolerance(1e-2);
        ROOT::Math::IntegratorOneDimOptions::SetDefaultIntegrator("Gauss");
        
        // Initialize |F(tx,ty)|^2 TF2
        FF2DtSquared = new TF2("", calc_FormFactor_2t2D, tx_min, tx_max, ty_min, ty_max, 4);
        FF2DtSquared->SetParameters(A, Vo, R, a0);

        // Initialize |F(tx,ty)^2| 2D histogram
        hist2D = new TH2D("", "|F(t_{x},t_{y})^{2}| 2D Histogram", bins, tx_min, tx_max, bins, ty_min, ty_max);
        hist2D->Sumw2();
        double step = tx_max/bins;
        for (int i=0;i<bins;i++) 
        {
            double tx = (i+1)*step;
            for (int j=0;j<bins;j++) 
            {
                double ty = (j+1)*step;
                double val = calc_FF_2t2D(tx, ty);
                hist2D->SetBinContent(i+1, j+1, val);
            }
        }

        // Initialize TF2 for transform: |F(tx,ty)|^2 -> G(x,y)
        transform_TF2 = new TF2("", [this, tx_min, tx_max, ty_min, ty_max] (double *var, double *par)
        {
            double x = var[0], y = var[1];
            TF2 trans_integral("", trueG1_integrand_dq_2D, tx_min, tx_max, ty_min, ty_max, 2);
            trans_integral.SetParameters(x, y);
            return trans_integral.Integral(tx_min, tx_max, ty_min, ty_max, 1e-12);
        }, x_min, x_max, y_min, y_max, 0);  

        // Initialize hist for transform: |F(tx,ty)|^2 -> G(x,y)
        transform_hist = new TH2D("", "Transform Form Factor", bins, x_min, x_max, bins, y_min, y_max);
        transform_hist->Sumw2();
        double step2 = x_max/bins;
        for(int i=0;i<bins;i++)
        {
            double x = (i+1)*step2;
            for(int j=0;j<bins;j++)
            {
                double y = (j+1)*step2;
                double val = trueG1_integral_dq_2D(x, y);
                transform_hist->SetBinContent(i+1, j+1, val);
            }
        }
    }


    // Calculate Form Factor and G integrand
    static double calc_FormFactor_2t2D(double *var, double *par) // |F(tx,ty)|^2
    {
        double A = par[0], Vo = par[1], R = par[2], a0 = par[3];
        double tx = var[0], ty = var[1];
        double t = tx+ty;
        double q = sqrt(t);
	    const double arg1 = q*R / hbarc;  
	    const double arg2 = hbarc / q;
	    const double arg3  = (sin(arg1)-arg1*cos(arg1)) * 4*pi*Vo * arg2*arg2*arg2 / A; 
	    double result = arg3 / (1. + (a0*a0 * q*q/(hbarc*hbarc)));  // [-]
        return result*result;
    }
    static double calc_FF_2t2D(double tx, double ty) // |F(qx,qy)|^2 [-]
    {
        double t = tx+ty, q = sqrt(t);
        double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7;
	    const double arg1 = q*R / hbarc;  
	    const double arg2 = hbarc / q;
	    const double arg3  = (sin(arg1)-arg1*cos(arg1)) * 4*pi*Vo * arg2*arg2*arg2 / A; 
        double result = arg3 / (1. + (a0*a0 * q*q/(hbarc*hbarc))); 
	    return result*result; 
    }
    // Calculations for transformation |F(tx,ty)|^2 -> G(x,y)
    static double calc_FF(double qx, double qy) // |F(qx,qy)|^2 [-]
    {
        double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7;
        double q = sqrt(qx*qx+qy*qy);
	    const double arg1 = q*R / hbarc;  
	    const double arg2 = hbarc / q;
	    const double arg3  = (sin(arg1)-arg1*cos(arg1)) * 4*pi*Vo * arg2*arg2*arg2 / A; 
        double result = arg3 / (1. + (a0*a0 * q*q/(hbarc*hbarc))); 
	    return result; 
    }
    static double trueG1_integrand_dq_2D(double *var, double *par) 
    {
        double x = par[0], y = par[1];
        double qx = var[0], qy = var[1];
        double r = sqrt(x*x+y*y), q = sqrt(qx*qx+qy*qy);
        return 2*pi*sqrt(2*pi/r)*sqrt(2*hbarc/(pi*q*r))*sin(q*r/hbarc)*calc_FF(qx,qy)*sqrt(q)*q*1/hbarc*1/hbarc*1/sqrt(hbarc); //[1/fm^3GeV^2] 
    }
    static double trueG1_integral_dq_2D(double x, double y)
    {
        double qq_min = 0, qq_max = 20;
        TF1 *f = new TF1("f", trueG1_integrand_dq_2D, qq_min, qq_max, 2);  
        f->SetParameters(x, y);
        return f->Integral(qq_min, qq_max, 1e-12); //[1/fm^3]
    }
    

    // Functions
    TF2 *getFormFactort2_2D() const {return FF2DtSquared;}
    TH2D *getFormFactort_hist() const {return hist2D;}
    TF2 *getTransformed_TF2() const {return transform_TF2;}
    TH2D *getTransformed_hist() const {return transform_hist;}


private:
    double A, Vo, R, a0;
    TF2 *FF2DtSquared, *transform_TF2;
    TH2D *hist2D, *transform_hist;
};


#endif