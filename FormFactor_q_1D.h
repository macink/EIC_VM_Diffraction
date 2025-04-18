#ifndef FormFactor_q_1D_H
#define FormFactor_q_1D_H

#include "EICvalueconst.h"

using namespace std;

/********************************************** 
 
 * Form Factor 1D: F(q) [-]
 
 * Parameters:
    Vo: Nuclear density (hard sphere approximation) [fm^-3]
    A: Atomic mass number
    R: Radius of nucleus [fm]
    a0: Range of Yukawa potential [fm]
 * Variables:
    q: Transverse momentum >> q = sqrt(qx^2+qy^2) [GeV]

**********************************************/
class FormFactor_q_1D
{
public:
    FormFactor_q_1D(double A_init, double Vo_init, double R_init, double a0_init, 
        double q_min, double q_max): A(A_init), Vo(Vo_init), R(R_init), a0(a0_init)
    {
        ROOT::Math::IntegratorOneDimOptions::SetDefaultAbsTolerance(1e-4);
        ROOT::Math::IntegratorOneDimOptions::SetDefaultRelTolerance(1e-2);
        ROOT::Math::IntegratorOneDimOptions::SetDefaultIntegrator("Gauss");

        // Initialize |F(q)|^2 TF1
        FFq2 = new TF1("Form Factor: |F(q)|^{2}", calc_FormFactor_2q1D, q_min, q_max, 4);
        FFq2->SetParameters(A, Vo, R, a0);
    }

    // Calculate Form Factor
    static double calc_FormFactor_2q1D(double *var, double *par) // |F(q)|^2 [-]
    {
        double A = par[0], Vo = par[1], R = par[2], a0 = par[3];
        double q = var[0];
	    const double arg1 = q*R / hbarc;  
	    const double arg2 = hbarc / q;
	    const double arg3  = (sin(arg1)-arg1*cos(arg1)) * 4*pi*Vo * arg2*arg2*arg2 / A; 
        double result = arg3 / (1. + (a0*a0 * q*q/(hbarc*hbarc)));
	    return result*result; 
    }
    

    // Functions
    TF1 *getFormFactorq_1D() const {return FFq2;}

private:
    double A, Vo, R, a0;
    TF1 *FFq2;
};


#endif
