#ifndef FormFactor_1D_H
#define FormFactor_1D_H

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
class FormFactor_1D
{
public:
    FormFactor_1D(double A_init, double Vo_init, double R_init, double a0_init, double q_min, double q_max, double t_min, double t_max, 
        double phi_min, double phi_max): A(A_init), Vo(Vo_init), R(R_init), a0(a0_init)
    {
        ROOT::Math::IntegratorOneDimOptions::SetDefaultAbsTolerance(1e-4);
        ROOT::Math::IntegratorOneDimOptions::SetDefaultRelTolerance(1e-2);
        ROOT::Math::IntegratorOneDimOptions::SetDefaultIntegrator("Gauss");

        // Initialize F(q) TF1
        FFq1 = new TF1("Form Factor: F(q)",calc_FormFactor_q1D,q_min,q_max,4);
        FFq1->SetParameters(A,Vo,R,a0);
        FFq1->SetParNames("A","Vo","R","a0");

        // Initialize |F(q)|^2 TF1
        FFq2 = new TF1("Form Factor: |F(q)|^{2}",calc_FormFactor_2q1D,q_min,q_max,4);
        FFq2->SetParameters(A,Vo,R,a0);
        FFq2->SetParNames("A","Vo","R","a0");

        // Initialize F(t) TF1
        FFt1 = new TF1("Form Factor: F(t)",calc_FormFactor_t1D,t_min,t_max,4);
        FFt1->SetParameters(A,Vo,R,a0);
        FFt1->SetParNames("A","Vo","R","a0");

        // Initialize |F(t)|^2 TF1
        FFt2 = new TF1("Form Factor: |F(t)|^{2}",calc_FormFactor_2t1D,t_min,t_max,4);
        FFt2->SetParameters(A,Vo,R,a0);
        FFt2->SetParNames("A","Vo","R","a0");

        // Initialize F(q) wedge cut TF1
        cutFFq = new TF1("Cut Form Factor: F(q)",[this,phi_min,phi_max](double *var, double *par){double q = var[0];  // q is variable
            TF1 ffq("ffq",calc_FormFactor_qWedge,phi_min,phi_max,5);
            ffq.SetParameters(par[0],par[1],par[2],par[3],q);
            return ffq.Integral(phi_min,phi_max,1e-12);
        },q_min,q_max,4);
        cutFFq->SetParameters(A,Vo,R,a0);
        cutFFq->SetParNames("A","Vo","R","a0");

        // Initialize |F(q)|^2 wedge cut TF1
        cutFFq2 = new TF1("Cut Form Factor: |F(q)|^{2}",[this,phi_min,phi_max](double *var, double *par){double q = var[0];  // q is variable
            TF1 ffq2("ffq",calc_FormFactor_qWedge,phi_min,phi_max,5);
            ffq2.SetParameters(par[0],par[1],par[2],par[3],q);
            double integral = ffq2.Integral(phi_min, phi_max,1e-12);
            return integral*integral;
        },q_min,q_max,4);
        cutFFq2->SetParameters(A,Vo,R,a0);
        cutFFq2->SetParNames("A","Vo","R","a0");

        // Initialize F(t) wedge cut TF1
        cutFFt = new TF1("Cut Form Factor: F(t)",[this,phi_min,phi_max](double *var, double *par){double t = var[0];  // t is variable
            TF1 fft("fft",calc_FormFactor_tWedge,phi_min,phi_max,5);
            fft.SetParameters(par[0],par[1],par[2],par[3],t);
            return fft.Integral(phi_min,phi_max,1e-12);
        },t_min,t_max,4);
        cutFFt->SetParameters(A,Vo,R,a0);
        cutFFt->SetParNames("A","Vo","R","a0");

        // Initialize |F(t)|^2 wedge cut TF1
        cutFFt2 = new TF1("Cut Form Factor: |F(t)|^{2}",[this,phi_min,phi_max](double *var, double *par){double t = var[0];  // t is variable
            TF1 fft2("fft2",calc_FormFactor_tWedge,phi_min,phi_max,5);
            fft2.SetParameters(par[0],par[1],par[2],par[3],t);
            double integral = fft2.Integral(phi_min,phi_max,1e-12);
            return integral*integral;
        },t_min,t_max,4);
        cutFFt2->SetParameters(A,Vo,R,a0);
        cutFFt2->SetParNames("A","Vo","R","a0");
    }

    // Calculate Form Factor
    static double calc_FormFactor_q1D(double *var, double *par) // F(q) [-]
    {
        double A = par[0], Vo = par[1], R = par[2], a0 = par[3];
        double q = var[0];
	    const double arg1 = q*R / hbarc;  
	    const double arg2 = hbarc / q;
	    const double arg3  = (sin(arg1)-arg1*cos(arg1)) * 4*pi*Vo*arg2*arg2*arg2 / A; 
        double result = arg3 / (1. + (a0 * a0 * q*q/(hbarc*hbarc)));
	    return result; 
    }
    static double calc_FormFactor_2q1D(double *var, double *par) // |F(q)|^2 [-]
    {
        double A = par[0], Vo = par[1], R = par[2], a0 = par[3];
        double q = var[0];
	    const double arg1 = q*R / hbarc;  
	    const double arg2 = hbarc / q;
	    const double arg3  = (sin(arg1)-arg1*cos(arg1)) * 4*pi*Vo*arg2*arg2*arg2 / A; 
        double result = arg3 / (1. + (a0 * a0 * q*q/(hbarc*hbarc)));
	    return result*result; 
    }
    static double calc_FormFactor_t1D(double *var, double *par) // F(t)
    {
        double A = par[0], Vo = par[1], R = par[2], a0 = par[3];
        double t = var[0];
        double q = sqrt(t);
        const double arg1 = q * R / hbarc;  
        const double arg2 = hbarc / q;
        const double arg3 = (sin(arg1) - arg1 * cos(arg1)) * 4 * pi * Vo * arg2 * arg2 * arg2 / A; 
        return arg3 / (1. + (a0 * a0 * q * q / (hbarc * hbarc)));  // [-]
    }
    static double calc_FormFactor_2t1D(double *var, double *par) // |F(t)|^2
    {
        double A = par[0], Vo = par[1], R = par[2], a0 = par[3];
        double t = var[0];
        double q = sqrt(t);
	    const double arg1 = q*R / hbarc;  
	    const double arg2 = hbarc / q;
	    const double arg3  = (sin(arg1)-arg1*cos(arg1)) * 4*pi*Vo*arg2*arg2*arg2 / A; 
	    return arg3 / (1. + (a0 * a0 * q*q/(hbarc*hbarc)))*arg3 / (1. + (a0 * a0 * q*q/(hbarc*hbarc)));  // [-]
    }
    

    // Calculations to integrate out wedge (phi from qy axis)
    static double calc_FormFactor_qWedge(double *var, double *par)
    {
        double A = par[0], Vo = par[1], R = par[2], a0 = par[3], q = par[4];
        double phi = var[0];
        double qx = q*sin(phi), qy = q*cos(phi), qq = sqrt(qx*qx+qy*qy);
	    const double arg1 = qq*R / hbarc;  
	    const double arg2 = hbarc / qq;
	    const double arg3  = (sin(arg1)-arg1*cos(arg1)) * 4*pi*Vo*arg2*arg2*arg2 / A; 
        double result = 1/(2*pi)*arg3 / (1. + (a0 * a0 * qq*qq/(hbarc*hbarc)));
        return result;  // [-]
    }
    static double calc_FormFactor2_qWedge(double *var, double *par)
    {
        double A = par[0], Vo = par[1], R = par[2], a0 = par[3], q = par[4];
        double phi = var[0];
        double qx = q*sin(phi), qy = q*cos(phi), qq = sqrt(qx*qx+qy*qy);
	    const double arg1 = qq*R / hbarc;  
	    const double arg2 = hbarc / qq;
	    const double arg3  = (sin(arg1)-arg1*cos(arg1)) * 4*pi*Vo*arg2*arg2*arg2 / A; 
        double result = 1/(2*pi)*arg3 / (1. + (a0 * a0 * qq*qq/(hbarc*hbarc)));
        return result*result;  // [-]
    }
    static double calc_FormFactor_tWedge(double *var, double *par)
    {
        double A = par[0], Vo = par[1], R = par[2], a0 = par[3], t = par[4];
        double phi = var[0];
        double q = sqrt(t), qx = q*sin(phi), qy = q*cos(phi), qq = sqrt(qx*qx+qy*qy);
	    const double arg1 = q*R / hbarc;  
	    const double arg2 = hbarc / q;
	    const double arg3  = (sin(arg1)-arg1*cos(arg1)) * 4*pi*Vo*arg2*arg2*arg2 / A; 
	    double result = 1/(2*pi)*arg3 / (1. + (a0 * a0 * q*q/(hbarc*hbarc))); 
        return result; // [-]
    }
    static double calc_FormFactor2_tWedge(double *var, double *par)
    {
        double A = par[0], Vo = par[1], R = par[2], a0 = par[3], t = par[4];
        double phi = var[0];
        double q = sqrt(t), qx = q*sin(phi), qy = q*cos(phi), qq = sqrt(qx*qx+qy*qy);
	    const double arg1 = q*R / hbarc;  
	    const double arg2 = hbarc / q;
	    const double arg3  = (sin(arg1)-arg1*cos(arg1)) * 4*pi*Vo*arg2*arg2*arg2 / A; 
	    double result = 1/(2*pi)*arg3 / (1. + (a0 * a0 * q*q/(hbarc*hbarc))); 
        return result*result; // [-]
    }


    // Functions
    TF1 *getFormFactorq2_1D() const {return FFq2;}
    TF1 *getFormFactort2_1D() const {return FFt2;}
    TF1 *getCutFormFactor_q() const {return cutFFq;}
    TF1 *getCutFormFactor_t() const {return cutFFt;}
    TF1 *getFormFactorq_1D() const {return FFq1;}
    TF1 *getFormFactort_1D() const {return FFt1;}
    TF1 *getCutFormFactor_q2() const {return cutFFq2;}
    TF1 *getCutFormFactor_t2() const {return cutFFt2;}

private:
    double A, Vo, R, a0;
    TF1 *FFq1, *FFt1, *cutFFq, *cutFFt, *FFq2, *FFt2, *cutFFq2, *cutFFt2, *FFsmearWedge_q, *FFsmearWedge_t;
};


#endif
