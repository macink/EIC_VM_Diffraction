#ifndef FormFactor_2D_H
#define FormFactor_2D_H

#include "EICvalueconst.h"

using namespace std;

/********************************************** 
 
 * Form Factor 2D: F(qx,qy) and F(tx,ty) [-]
 
 * Parameters:
    Vo: Nuclear density (hard sphere approximation) [fm^-3]
    A: Atomic mass number
    R: Radius of nucleus [fm]
    a0: Range of Yukawa potential [fm]
 * Variables:
    q: Transverse momentum >> q = sqrt(qx^2+qy^2) [GeV]
    t: Four-momentum transfer squared >> t = tx+ty = q^2 [GeV^2]

**********************************************/
class FormFactor_2D
{
public:
    FormFactor_2D(double A_init, double Vo_init, double R_init, double a0_init,double qy_min, double qy_max, double qx_min, double qx_max, double tx_min, 
        double tx_max, double ty_min, double ty_max, double bins, double qxBins, double qyBins, double txBins, double tyBins, double x_min,
        double x_max, double y_min, double y_max): A(A_init), Vo(Vo_init), R(R_init), a0(a0_init)
    {
        ROOT::Math::IntegratorOneDimOptions::SetDefaultAbsTolerance(1e-4);
        ROOT::Math::IntegratorOneDimOptions::SetDefaultRelTolerance(1e-2);
        ROOT::Math::IntegratorOneDimOptions::SetDefaultIntegrator("Gauss");

        // Initialize F(qx,qy) TF2
        FF2Dq = new TF2("Form Factor", calc_FormFactor_q2D, qx_min, qx_max, qy_min, qy_max, 4);
        FF2Dq->SetParameters(A,Vo,R,a0);
        FF2Dq->SetParNames("A","Vo","R","a0");

        // Initialize |F(qx,qy)|^2 TF2
        FF2DqSquared = new TF2("Form Factor", calc_FormFactor_2q2D, qx_min, qx_max, qy_min, qy_max, 4);
        FF2DqSquared->SetParameters(A,Vo,R,a0);
        FF2DqSquared->SetParNames("A","Vo","R","a0");

        // Initialize F(tx,ty) TF2
        FF2Dt = new TF2("Form Factor", calc_FormFactor_t2D, tx_min, tx_max, ty_min, ty_max, 4);
        FF2Dt->SetParameters(A,Vo,R,a0);
        FF2Dt->SetParNames("A","Vo","R","a0");

        // Initialize |F(tx,ty)|^2 TF2
        FF2DtSquared = new TF2("Form Factor", calc_FormFactor_2t2D, tx_min, tx_max, ty_min, ty_max, 4);
        FF2DtSquared->SetParameters(A,Vo,R,a0);
        FF2DtSquared->SetParNames("A","Vo","R","a0");

        // Initialize hist for transform: |F(t)|^2 -> G(r)
        hankelTransformFF = new TH2D("hankelTransformFF","Transform Form Factor",bins,x_min,x_max,bins,y_min,y_max);
        hankelTransformFF->Sumw2();
        double step = x_max/bins;
        for(int i=0;i<bins;i++)
        {
            double x = (i+1)*step;
            for(int j=0;j<bins;j++)
            {
                double y = (j+1)*step;
                double val = trueG1_integral_dq_1D(x,y);
                hankelTransformFF->SetBinContent(i+1,j+1,val);
            }
            //cout << "x: " << x << endl;
        }
        hankelTransformFF->SetTitle("Transformed WS: |F(tx,ty)|^{2} #rightarrow G(x,y)");
        hankelTransformFF->GetYaxis()->SetTitle("y [fm]");
	    hankelTransformFF->GetXaxis()->SetTitle("x [fm]");    
        hankelTransformFF->GetZaxis()->SetTitle("G(x,y) [fm^{-3}]");
    }


    // Calculate Form Factor
    static double calc_FormFactor_q2D(double *var, double *par) // F(qx,qy) [-]
    {
        double A = par[0], Vo = par[1], R = par[2], a0 = par[3];
        double qx = var[0], qy = var[1];
        double q = sqrt(qx*qx+qy*qy);
	    const double arg1 = q*R / hbarc;  
	    const double arg2 = hbarc / q;
	    const double arg3  = (sin(arg1)-arg1*cos(arg1)) * 4*pi*Vo*arg2*arg2*arg2 / A; 
        double result = arg3 / (1. + (a0 * a0 * q*q/(hbarc*hbarc))); 
	    return result;
    }
    static double calc_FormFactor_2q2D(double *var, double *par) // |F(qx,qy)|^2 [-]
    {
        double A = par[0], Vo = par[1], R = par[2], a0 = par[3];
        double qx = var[0], qy = var[1];
        double q = sqrt(qx*qx+qy*qy);
	    const double arg1 = q*R / hbarc;  
	    const double arg2 = hbarc / q;
	    const double arg3  = (sin(arg1)-arg1*cos(arg1)) * 4*pi*Vo*arg2*arg2*arg2 / A; 
        double result = arg3 / (1. + (a0 * a0 * q*q/(hbarc*hbarc))); 
	    return result*result; 
    }
    static double calc_FormFactor_t2D(double *var, double *par) // F(tx,ty)
    {
        double A = par[0], Vo = par[1], R = par[2], a0 = par[3];
        double tx = var[0], ty = var[1];
        double t = tx+ty;
        double q = sqrt(t);
	    const double arg1 = q*R / hbarc;  
	    const double arg2 = hbarc / q;
	    const double arg3  = (sin(arg1)-arg1*cos(arg1)) * 4*pi*Vo*arg2*arg2*arg2 / A; 
	    return arg3 / (1. + (a0 * a0 * q*q/(hbarc*hbarc)));  // [-]
    }
    static double calc_FormFactor_2t2D(double *var, double *par) // |F(tx,ty)|^2
    {
        double A = par[0], Vo = par[1], R = par[2], a0 = par[3];
        double tx = var[0], ty = var[1];
        double t = tx+ty;
        double q = sqrt(t);
	    const double arg1 = q*R / hbarc;  
	    const double arg2 = hbarc / q;
	    const double arg3  = (sin(arg1)-arg1*cos(arg1)) * 4*pi*Vo*arg2*arg2*arg2 / A; 
	    return arg3 / (1. + (a0 * a0 * q*q/(hbarc*hbarc)))*arg3 / (1. + (a0 * a0 * q*q/(hbarc*hbarc)));  // [-]
    }
    static double calc_FF_2q2D(double qx, double qy) // |F(qx,qy)|^2 [-]
    {
        double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7;
        double q = sqrt(qx*qx+qy*qy);
	    const double arg1 = q*R / hbarc;  
	    const double arg2 = hbarc / q;
	    const double arg3  = (sin(arg1)-arg1*cos(arg1)) * 4*pi*Vo*arg2*arg2*arg2 / A; 
        double result = arg3 / (1. + (a0 * a0 * q*q/(hbarc*hbarc))); 
	    return result*result; 
    }
    static double trueG1_integrand_dq_1D(double *var, double *par) 
    {
        double x = par[0], y = par[1];
        double qx = var[0], qy = var[1];
        double q = sqrt(qx*qx+qy*qy), r = sqrt(x*x+y*y);
        return 2*pi*sqrt(2*pi/r)*sqrt(2*hbarc/(pi*q*r))*sin(q*r/hbarc)*calc_FF_2q2D(qx,qy)*sqrt(q)*q*1/hbarc*1/hbarc*1/sqrt(hbarc); //[1/fm^3GeV^2] 
    }
    static double trueG1_integral_dq_1D(double x, double y)
    {
        double q_min = 0.000005;
        double q_max = 5;
        TF1 *f = new TF1("f",trueG1_integrand_dq_1D,q_min,q_max,2);  
        f->SetParameters(x,y);
        return f->Integral(q_min,q_max,1e-6); //[1/fm^3]
    }

    

    // Define hists functions
    TH2D *create_FormFactorq_2Dcontour(const char *name, const char *title, int qxBins, double qx_min, double qx_max, int qyBins, double qy_min, double qy_max)
    {
        TF2 *FF2 = getFormFactorq2_2D();
        // Create a 2D histogram to represent the function
        TH2D *hist2D = new TH2D(name, title, qxBins, qx_min, qx_max, qyBins, qy_min, qy_max);
        hist2D->Sumw2();
        for (int i=1;i<=qxBins;i++) 
        {
            for (int j=1;j<=qyBins;j++) 
            {
                double qx = hist2D->GetXaxis()->GetBinCenter(i);
                double qy = hist2D->GetYaxis()->GetBinCenter(j);
                hist2D->SetBinContent(i, j, FF2->Eval(qx, qy));
            }
            //cout << "bin: " << i << endl;
        }
        //hist2D->SetTitle("Form Factor 2D: F(qx,qy)");
        hist2D->GetYaxis()->SetTitle("qy [GeV]");
	    hist2D->GetXaxis()->SetTitle("qx [GeV]");
        hist2D->GetZaxis()->SetTitle("F(qx,qy)");
        return hist2D;
    }
    TH2D *create_FormFactort_2Dcontour(const char *name, const char *title, int txBins, double tx_min, double tx_max, int tyBins, double ty_min, double ty_max)
    {
        TF2 *FF2 = getFormFactort2_2D();
        // Create a 2D histogram to represent the function
        TH2D *hist2D = new TH2D(name, title, txBins, tx_min, tx_max, tyBins, ty_min, ty_max);
        hist2D->Sumw2();
        for (int i=1;i<=txBins;i++) 
        {
            for (int j=1;j<=tyBins;j++) 
            {
                double tx = hist2D->GetXaxis()->GetBinCenter(i);
                double ty = hist2D->GetYaxis()->GetBinCenter(j);
                hist2D->SetBinContent(i, j, FF2->Eval(tx, ty));
            }
            //cout << "bin: " << i << endl;
        }
        //hist2D->SetTitle("Form Factor 2D: F(tx,ty)");
        hist2D->GetYaxis()->SetTitle("ty [GeV^{2}]");
	    hist2D->GetXaxis()->SetTitle("tx [GeV^{2}]");
        hist2D->GetZaxis()->SetTitle("F(tx,ty)");
        return hist2D;
    }
    

    // Functions
    TF2 *getFormFactorq_2D() const {return FF2Dq;}
    TF2 *getFormFactorq2_2D() const {return FF2DqSquared;}
    TF2 *getFormFactort_2D() const {return FF2Dt;}
    TF2 *getFormFactort2_2D() const {return FF2DtSquared;}
    TH2D *getFFtransform_2D() const {return hankelTransformFF;}


private:
    double A, Vo, R, a0;
    TF2 *FF2Dq, *FF2DqSquared, *FF2Dt, *FF2DtSquared;
    TH2D *hankelTransformFF;
};


#endif
