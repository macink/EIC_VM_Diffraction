#ifndef FormFactor_q_2D_H
#define FormFactor_q_2D_H

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
class FormFactor_q_2D
{
public:
    FormFactor_q_2D(double A_init, double Vo_init, double R_init, double a0_init,double qy_min, double qy_max, double qx_min, 
        double qx_max, double bins): A(A_init), Vo(Vo_init), R(R_init), a0(a0_init)
    {
        // Initialize |F(qx,qy)|^2 TF2
        FF2DqSquared = new TF2("Form Factor", calc_FormFactor_2q2D, qx_min, qx_max, qy_min, qy_max, 4);
        FF2DqSquared->SetParameters(A, Vo, R, a0);

        // Initialize |F(qx,qy)^2| 2D histogram
        hist2D = new TH2D("", "|F(q_{x},q_{y})^{2}| 2D Histogram", bins, qx_min, qx_max, bins, qy_min, qy_max);
        hist2D->Sumw2();
        double pars[] = {A, Vo, R, a0};
        for (int i=1;i<=bins;i++) 
        {
            for (int j=1;j<=bins;j++) 
            {
                double qx = hist2D->GetXaxis()->GetBinCenter(i);
                double qy = hist2D->GetYaxis()->GetBinCenter(j);
                double vars[] = {qx, qy};
                double val = calc_FormFactor_2q2D(vars, pars);
                hist2D->SetBinContent(i, j, val);
            }
        }
    }


    // Calculate Form Factor
    static double calc_FormFactor_2q2D(double *var, double *par) // |F(qx,qy)|^2 [-]
    {
        double A = par[0], Vo = par[1], R = par[2], a0 = par[3];
        double qx = var[0], qy = var[1];
        double q = sqrt(qx*qx+qy*qy);
	    const double arg1 = q*R / hbarc;  
	    const double arg2 = hbarc / q;
	    const double arg3  = (sin(arg1)-arg1*cos(arg1)) * 4*pi*Vo * arg2*arg2*arg2 / A; 
        double result = arg3 / (1. + (a0*a0 * q*q/(hbarc*hbarc))); 
	    return result*result; 
    }
    

    // Functions
    TF2 *getFormFactorq2_2D() const {return FF2DqSquared;}
    TH2D *getFormFactorq_2D_hist() const {return hist2D;}

private:
    double A, Vo, R, a0;
    TF2 *FF2DqSquared;
    TH2D *hist2D;
};


#endif