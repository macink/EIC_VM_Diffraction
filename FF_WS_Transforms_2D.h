#ifndef FF_WS_Transforms_2D_H
#define FF_WS_Transforms_2D_H

#include "EICvalueconst.h"

using namespace std;

/********************************************** 
 
 * 2D Fourier-Bessel Transformations of Form Factor: |F(tx,ty)|^2 and Woods-Saxon Distribution: G(x,y)
 
 * Parameters:
    Vo: Nuclear density (hard sphere approximation) [fm^-3]
    A: Atomic mass number
    R: Radius of nucleus [fm]
    a0: Range of Yukawa potential [fm]
 * Variables:
    t: Four-momentum transfer squared >> t = tx+ty = q^2 [GeV^2]
    r: position space >> r = sqrt(x^2+y^2)

**********************************************/
class FF_WS_Transforms_2D
{
public:

    // Calculations for form factor and G integrand
    static double calc_FF_2q2D(double tx, double ty) 
    {
        double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7;
        double q = sqrt(tx+ty);
	    const double arg1 = q*R / hbarc;  
	    const double arg2 = hbarc / q;
	    const double arg3  = (sin(arg1)-arg1*cos(arg1)) * 4*pi*Vo * arg2*arg2*arg2 / A; 
        double result = arg3 / (1. + (a0*a0 * q*q/(hbarc*hbarc))); 
	    return result*result; // [-]
    }
    static double trueG1_integrand_dq_2D(double *var, double *par) // |F(tx,ty)|^2 -> G(x,y)
    {
        double x = par[0], y = par[1];
        double tx = var[0], ty = var[1];
        double q = sqrt(tx+ty), r = sqrt(x*x+y*y);
        return 2*pi*sqrt(2*pi/r)*sqrt(2*hbarc/(pi*q*r))*sin(q*r/hbarc)*calc_FF_2q2D(tx,ty)*sqrt(q)*q*1/hbarc*1/hbarc*1/sqrt(hbarc); //[1/fm^3GeV^2] 
    }
    

    // Calculations for woods-saxon and F integrand
    static double calc_WS2D(double x, double y)
    {
        double Vo = 2.12, R = 6.38, a = 0.535; 
        double r = sqrt(x*x + y*y); 
        return Vo/(1. + exp((r-R)/a)); // [1/fm^3]
    }
    static double F2_integrand_dr(double *var, double *par)  // G(x,y) -> |F(tx,ty)|^2
    {
        double tx = par[0], ty = par[1];
        double x = var[0], y = var[1];
        double r = sqrt(x*x+y*y);
        double q = sqrt(tx+ty);
        return 2*pi*sqrt(2*pi/q)*sqrt(2*hbarc/(pi*q*r))*sin(q *r/hbarc)*calc_WS2D(x,y)*sqrt(hbarc)*sqrt(r)*r; //[1/fm]
    }
    

    // Functions to transform a given TF1 object
    TF2* transformTF2_FtoG(TF2* inputTF2, double x_min, double x_max, double y_min, double y_max)
    {
        ROOT::Math::IntegratorOneDimOptions::SetDefaultAbsTolerance(1e-4);
        ROOT::Math::IntegratorOneDimOptions::SetDefaultRelTolerance(1e-2);
        ROOT::Math::IntegratorOneDimOptions::SetDefaultIntegrator("Gauss");
    
        double qq_min = 0, qq_max = 5;
        TF2 *transform2 = new TF2("Fourier-Bessel Transformation: |F(t_{x},t_{y})|^{2} -> G(x,y)", [this, qq_min, qq_max] (double *var, double *par)
        {
            double x = var[0], y = var[1];
            TF1 trans_integral("trans_integral", trueG1_integrand_dq_2D, qq_min, qq_max, 2);
            trans_integral.SetParameters(x, y);
            return trans_integral.Integral(qq_min, qq_max, 1e-19);
        }, x_min, x_max, y_min, y_max, 0);  
        return transform2;
    }
    TF2* transformTF2_GtoF(TF2* inputTF2, double tx_min, double tx_max, double ty_min, double ty_max)
    {
        ROOT::Math::IntegratorOneDimOptions::SetDefaultAbsTolerance(1e-4);
        ROOT::Math::IntegratorOneDimOptions::SetDefaultRelTolerance(1e-2);
        ROOT::Math::IntegratorOneDimOptions::SetDefaultIntegrator("Gauss");
    
        double xx_min = 0, xx_max = 30, yy_min = 0, yy_max = 30;
        TF2 *transform2 = new TF2("Fourier-Bessel Transformation: G(x,y) -> |F(t_{x},t_{y})|^{2}", [this,xx_min, xx_max, yy_min, yy_max] (double *var, double *par)
        {
            double tx = var[0], ty = var[1];
            TF2 trans_integral("trans_integral", F2_integrand_dr, xx_min, xx_max, yy_min, yy_max, 2);
            trans_integral.SetParameters(tx, ty);
            return trans_integral.Integral(xx_min, xx_max, yy_min, yy_max, 1e-12)*trans_integral.Integral(xx_min, xx_max, yy_min, yy_max, 1e-12);
        }, tx_min, tx_max, ty_min, ty_max, 0);  
        return transform2;
    }

};


#endif
