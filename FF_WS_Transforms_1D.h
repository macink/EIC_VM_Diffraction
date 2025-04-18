#ifndef FF_WS_Transforms_1D_H
#define FF_WS_Transforms_1D_H

#include "EICvalueconst.h"

using namespace std;

/********************************************** 
 
 * 1D Fourier-Bessel Transformations for Form Factor |F(t)|^2 and Woods-Saxon Distribution: G(r)
    G(r)     ->  |F(t)|^2
    |F(t)|^2 ->  G(r)
 
 * Parameters:
    Vo: Nuclear density (hard sphere approximation) [fm^-3]
    A: Atomic mass number
    R: Radius of nucleus [fm]
    a0: Range of Yukawa potential [fm]
 * Variables:
    t: Four-momentum transfer squared >> t = tx+ty = q^2 [GeV^2]
    r: position space >> r = sqrt(x^2+y^2)

**********************************************/
class FF_WS_Transforms_1D
{
public:
    
    // Calculate form factor and G(r) integrand
    static double calcFF_1D(double q)
    {
	    if(q==0){return 0;}
        else{
            double A = 197, Vo = 2.12, R = 6.38, a0 = .70;
	        const double arg1 = q * R / hbarc;  
	        const double arg2 = hbarc /q;
	        const double sph  = (sin(arg1) - arg1 * cos(arg1)) * 4*pi*Vo * arg2 * arg2 *arg2/ double(A); 
	        return sph * 1/(1. + (a0 * a0 * q*q/(hbarc*hbarc)));  // [-]
        }
    }
    static double trueG1_integrand_dq_1D(double *x, double *par) // |F(t)|^2 -> G(r)
    {
        double r = par[0];
        double q = x[0];
        //double t = q*q;
        return 2*pi*sqrt(2*pi/r)*sqrt(2*hbarc/(pi*q*r))*sin(q*r/hbarc)*calcFF_1D(q)*sqrt(q)*q*1/hbarc*1/hbarc*1/sqrt(hbarc); //[1/fm^3GeV^2] 
    }
   
    
    // Calculate Woods-Saxon density distribution and F(t) integrand
    static double calc_WS1D(double r)
    {
        double Vo = 2.12, A = 197, a = depth_Au, R = 6.38;
        return Vo/(1. + exp((r-R)/a)); // [1/fm^3]
    }
    static double F_integrand_dr(double *var, double *par)  // G(r) -> |F(t)|^2
    {
        double t = par[0];
        double r = var[0];
        double q = sqrt(t);
        return 2*pi*sqrt(2*pi/q)*sqrt(2*hbarc/(pi*q*r))*sin(q*r/hbarc)*calc_WS1D(r)*sqrt(hbarc)*sqrt(r)*r; //[1/fm]
    }
    

};


#endif
