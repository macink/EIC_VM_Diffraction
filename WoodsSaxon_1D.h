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
        // Set preferences for integration
        ROOT::Math::IntegratorOneDimOptions::SetDefaultAbsTolerance(1e-4);
        ROOT::Math::IntegratorOneDimOptions::SetDefaultRelTolerance(1e-2);
        ROOT::Math::IntegratorOneDimOptions::SetDefaultIntegrator("Gauss");
        
        // Initialize TF1
        WoodsSaxon_TF1 = new TF1("", calc_WoodsSaxon1D, r_min, r_max, 3);
        WoodsSaxon_TF1->SetParameters(Vo, R, a);

        // Initialize TH1D
        WoodsSaxon_TH1D = new TH1D("", "", bins, r_min, r_max);
        WoodsSaxon_TH1D->Sumw2();
        double r_step = r_max/bins;
        for(int i=0;i<bins;i++)
        {
            double r = (i+1)*r_step;
            double vars[] = {r};
            double pars[] = {Vo,R,a};
            double val = calc_WoodsSaxon1D(vars, pars);
            WoodsSaxon_TH1D->SetBinContent(i+1, val);
        }
    
        // Initialize TF1 for Transformation: G(r) -> |F(t)|^2
        /*double r_lower_bound = 0, r_upper_bound = 30;  //Limits of integration
        transform_WoodsSaxon_TF1 = new TF1("", [this , r_lower_bound, r_upper_bound] (double *var, double *par)
        {
            double t = var[0];
            TF1 transform_WoodsSaxon_integral("", F_integrand_dr, r_lower_bound, r_upper_bound, 4);
            transform_WoodsSaxon_integral.SetParameters(Vo, R, a, t);
            return transform_WoodsSaxon_integral.Integral(r_lower_bound, r_upper_bound, 1e-12);
        }, t_min, t_max, 3);  
        transform_WoodsSaxon_TF1->SetParameters(Vo, R, a);
        */

        // Initialize hist for transform: G(r) -> |F(t)|^2
        transform_WoodsSaxon_TH1D = new TH1D("", "", bins, t_min, t_max);
        transform_WoodsSaxon_TH1D->Sumw2();
        double t_step = t_max/bins;
        for(int i=0;i<bins;i++)
        {
            double t = (i+1)*t_step;
	        double val = F_integral_dr(t, Vo, R, a);
            double result = val*val;
            transform_WoodsSaxon_TH1D->SetBinContent(i+1, result);
        }
        transform_WoodsSaxon_TH1D->GetYaxis()->SetTitle("|F(t)|^{2}");
	    transform_WoodsSaxon_TH1D->GetXaxis()->SetTitle("t [GeV^{2}]");

    }

    // Calculate Woods-Saxon density distribution
    static double calc_WoodsSaxon1D(double *var, double *par)
    {
        double Vo = par[0], R = par[1], a = par[2]; 
        double r = var[0]; 
        return Vo/(1. + exp((r-R)/a)); 
    }
    // Calculate Fourier-Bessel transformation integrand
    static double F_integrand_dr(double *var, double *par) 
    {
        double Vo = par[0], R = par[1], a = par[2], t = par[3];
        double r = var[0];
        double q = sqrt(t);
        double vars[] = {r};
        double pars[] = {Vo, R, a};
        return 2*pi*sqrt(2*pi/q)*sqrt(2*hbarc/(pi*q*r))*sin(q*r/hbarc)*calc_WoodsSaxon1D(vars, pars)*sqrt(hbarc)*sqrt(r)*r; //[1/fm]
    }
     // G(r) -> |F(t)|^2
    static double F_integral_dr(double t, double Vo, double R, double a)
    {
        double r_lower_bound = 0., r_upper_bound = 100;
        TF1 *f = new TF1("", F_integrand_dr, r_lower_bound, r_upper_bound, 4); 
        double params[] = {t, Vo, R, a};
        f->SetParameters(params);
        return f->Integral(r_lower_bound, r_upper_bound, 1e-12); // [-]
    }


    TF1 *getWoodsSaxon1D_fun() const {return WoodsSaxon_TF1;}
    //TF1 *get_transformed_WoodsSaxon1D_fun() const {return transform_WoodsSaxon_TF1;}
    TH1D *get_transformed_WoodsSaxon1D_hist() const {return transform_WoodsSaxon_TH1D;}
    TH1D *getWoodsSaxon1D_hist() const {return WoodsSaxon_TH1D;}


private:
    double Vo, R, a;
    TF1 *WoodsSaxon_TF1, *transform_WoodsSaxon_TF1;
    TH1D *transform_WoodsSaxon_TH1D, *WoodsSaxon_TH1D;
};

#endif