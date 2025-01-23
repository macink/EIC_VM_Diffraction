#ifndef WoodsSaxon_1D_H
#define WoodsSaxon_1D_H

#include "EICvalueconst.h"

using namespace std;



/************************************************ 
 
 *  Woods-Saxon Dentsity Distribution 1D: G(r) [1/fm^-3]
 
 * Parameters:
    Vo: Nuclear density (hard sphere approximation) [fm^-3]
    R: Radius of nucleus [fm]
    a: Skin-depth [fm]
 * Variables:
    r: Radial distance from center of nucleus [fm]
        r = sqrt(x^2+y^2)
 
**************************************************/ 
class WoodsSaxon_1D
{
public:

    WoodsSaxon_1D(double Vo_init, double R_init, double a_init, double r_min, double r_max,
        double bins, double t_min, double t_max): Vo(Vo_init), R(R_init), a(a_init)
    {
        ROOT::Math::IntegratorOneDimOptions::SetDefaultAbsTolerance(1e-4);
        ROOT::Math::IntegratorOneDimOptions::SetDefaultRelTolerance(1e-2);
        ROOT::Math::IntegratorOneDimOptions::SetDefaultIntegrator("Gauss");
        
        // Initialize TF1
        WS1 = new TF1("Woods-Saxon 1D", calc_WoodsSaxon1D, r_min, r_max, 3);
        WS1->SetParameters(Vo, R, a);

        // Initialize TH1D
        WS1_hist = new TH1D("Woods-Saxon 1D", "Woods-Saxon", bins, r_min, r_max);
        WS1_hist->Sumw2();
        double step2 = r_max/bins;
        for(int i=0;i<bins;i++)
        {
            double r = (i+1)*step2;
            double vars[] = {r};
            double pars[] = {Vo,R,a};
            double val = calc_WoodsSaxon1D(vars, pars);
            WS1_hist->SetBinContent(i+1, val);
        }
    
        // Initialize TF1 for Transformation: G(r) -> |F(t)|^2
        double rr_min = 0, rr_max = 30;
        transform_TF1 = new TF1("Fourier-Bessel Transformation: G(r) -> |F(t)|^{2}", [this , rr_min, rr_max] (double *var, double *par)
        {
            double t = var[0];
            TF1 transform_integral("trans_integral", F_integrand_dr, rr_min, rr_max, 1);
            transform_integral.SetParameters(t);
            return transform_integral.Integral(rr_min, rr_max, 1e-12)*transform_integral.Integral(rr_min, rr_max, 1e-12);
        }, t_min, t_max, 0);  

        // Initialize hist for transform: G(r) -> |F(t)|^2
        transform_hist = new TH1D("hankelTransformWS", "Transform Woods-Saxon", bins, t_min, t_max);
        transform_hist->Sumw2();
        double step = t_max/bins;
        for(int i=0;i<bins;i++)
        {
            double t = (i+1)*step;
	        double val = F_integral_dr(t);
            double result = val*val;
            transform_hist->SetBinContent(i+1, result);
        }
        transform_hist->GetYaxis()->SetTitle("|F(t)|^{2}");
	    transform_hist->GetXaxis()->SetTitle("t [GeV^{2}]");

    }

    // Calculate Woods-Saxon density distribution
    static double calc_WoodsSaxon1D(double *var, double *par)
    {
        double Vo = par[0], R = par[1], a = par[2]; 
        double r = var[0]; 
        return Vo/(1. + exp((r-R)/a)); 
    }
    static double calc_WS1D(double r)
    {
        double Vo = 2.12, A = 197, a = depth_Au, R = 6.38;
        return Vo/(1. + exp((r-R)/a)); // [1/fm^3]
    }


    // Calculate Fourier-Bessel transformations
    static double F_integrand_dr(double *var, double *par)  // G(r) -> |F(t)|^2
    {
        double t = par[0];
        double r = var[0];
        double q = sqrt(t);
        return 2*pi*sqrt(2*pi/q)*sqrt(2*hbarc/(pi*q*r))*sin(q*r/hbarc)*calc_WS1D(r)*sqrt(hbarc)*sqrt(r)*r; //[1/fm]
    }
    static double F_integral_dr(double t)
    {
        double rr_min = 0., rr_max = 30;
        TF1 *f = new TF1("f", F_integrand_dr, rr_min, rr_max, 1);  
        f->SetParameters(t);
        return f->Integral(rr_min, rr_max, 1e-12); // [-]
    }


    TF1 *getWoodsSaxon1D() const {return WS1;}
    TF1 *getWStransform() const {return transform_TF1;}
    TH1D *getWStransformHist() const {return transform_hist;}
    TH1D *getWS_hist() const {return WS1_hist;}


private:
    double Vo, R, a;
    TF1 *WS1, *transform_TF1;
    TH1D *transform_hist, *WS1_hist;
};

#endif