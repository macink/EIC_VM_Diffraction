#ifndef WoodsSaxon_2D_H
#define WoodsSaxon_2D_H

#include "EICvalueconst.h"

using namespace std;



/************************************************ 
 
 *  Woods-Saxon Dentsity Distribution 2D: G(x,y) [1/fm^-3]
 
 * Parameters:
    Vo: Nuclear density (hard sphere approximation) [fm^-3]
    R: Radius of nucleus [fm]
    a: Skin-depth [fm]
 * Variables:
    r: Radial distance from center of nucleus [fm]
        r = sqrt(x^2+y^2)
 
**************************************************/ 
class WoodsSaxon_2D
{
public:

    WoodsSaxon_2D(double Vo, double R, double a, double x_min, double x_max, double y_min, double y_max, 
        double bins, double tx_min, double tx_max, double ty_min, double ty_max)
    {
        // Set default integration preferences
        ROOT::Math::IntegratorOneDimOptions::SetDefaultAbsTolerance(1.E-6);
        ROOT::Math::IntegratorOneDimOptions::SetDefaultRelTolerance(1.E-6);
        ROOT::Math::IntegratorOneDimOptions::SetDefaultIntegrator("Gauss");

        // Initialize TF2
        WoodsSaxon_TF2 = new TF2("", calc_WoodsSaxon2D, x_min, x_max, y_min, y_max, 3);
        WoodsSaxon_TF2->SetParameters(Vo, R, a);

        // Initialize Hist for G(x,y)
        WoodsSaxon_TH2D = new TH2D("", "", bins, x_min, x_max, bins, y_min, y_max);
        WoodsSaxon_TH2D->Sumw2();
        for(int i=1;i<=bins;i++) 
        {
            for(int j=1;j<=bins;j++) 
            {
            double x = WoodsSaxon_TH2D->GetXaxis()->GetBinCenter(i);
            double y = WoodsSaxon_TH2D->GetYaxis()->GetBinCenter(j);
            WoodsSaxon_TH2D->SetBinContent(i, j, WoodsSaxon_TF2->Eval(x, y)); 
            }
        }
        WoodsSaxon_TH2D->GetYaxis()->SetTitle("y [fm]");
        WoodsSaxon_TH2D->GetXaxis()->SetTitle("x [fm]");

        // Initialize TF2 for Transformation: G(x,y) -> |F(tx,ty)|^2
        /*transform_WoodsSaxon_TF2 = new TF2("", [this, x_min, x_max, y_min, y_max] (double *var, double *par)
        {
            double r_min = 0., r_max = 100;
            double tx = var[0], ty = var[1];
            TF1 transform_WoodsSaxon_integral("", F_integrand_dr, r_min, r_max, 5);
            transform_WoodsSaxon_integral.SetParameters(tx, ty, par[2], par[3], par[4]);
            return transform_WoodsSaxon_integral.Integral(r_min, r_max, 1e-12);
        }, tx_min, tx_max, ty_min, ty_max, 0);  
        */
       
        // Initialize hist for transform: G(x,y) -> |F(tx,ty)|^2
        transform_WoodsSaxon_TH2D = new TH2D("", "",bins ,tx_min, tx_max, bins, ty_min, ty_max);
        transform_WoodsSaxon_TH2D->Sumw2();
        double step = tx_max/bins;
        for(int i=0; i<bins; i++)
        {
            double tx = (i+1)*step;
            for(int j=0; j<bins; j++)
            {
                double ty = (j+1)*step;
	            double val = F_integral_dr(tx, ty, Vo, R, a);
                double result = val*val;
                transform_WoodsSaxon_TH2D->SetBinContent(i+1, j+1, result);
            }
        }
        transform_WoodsSaxon_TH2D->GetYaxis()->SetTitle("t_{y} [GeV^{2}]");
	    transform_WoodsSaxon_TH2D->GetXaxis()->SetTitle("t_{x} [GeV^{2}]");
    }


    // Calculations for Woods-Saxon density distribution
    static double calc_WoodsSaxon2D(double *var, double *par)
    {
        double Vo = par[0], R = par[1], a = par[2]; 
        double x = var[0], y = var[1];
        double r = sqrt(x*x + y*y); 
        return Vo/(1. + exp((r-R)/a)); // [1/fm^3]
    }
    // F integrand to transform Woods-Saxon to momentum space
    static double F_integrand_dr(double *var, double *par)
    {
        double tx = par[0], ty = par[1], Vo = par[2], R = par[3], a = par[4];
        double x = var[0], y = var[1];
        double r = sqrt(x*x+y*y);
        double q = sqrt(tx+ty);
        double pars[] = {Vo, R, a};
        double vars[] = {x, y};
        return 2*pi*sqrt(2*pi/q)*sqrt(2*hbarc/(pi*q*r))*sin(q *r/hbarc)*calc_WoodsSaxon2D(vars, pars)*sqrt(hbarc)*sqrt(r)*r; //[1/fm]
    }
    // G(x,y) -> |F(tx,ty)|^2
    static double F_integral_dr(double tx, double ty, double Vo, double R, double a)
    {
        double r_min = 0., r_max = 30;
        TF1 *f = new TF1("", F_integrand_dr, r_min, r_max, 5);  
        f->SetParameters(tx, ty, Vo, R, a);
        return f->Integral(r_min, r_max, 1e-12); // [-]
    }


    // Functions
    TF2 *getWoodsSaxon2D_fun() const {return WoodsSaxon_TF2;}
    TH2D *getWoodsSaxon2D_hist() const {return WoodsSaxon_TH2D;}
    //TF2 *get_transformed_WoodsSaxon2D_fun() const {return transform_WoodsSaxon_TF2;}
    TH2D *get_transformed_WoodsSaxon2D_hist() const {return transform_WoodsSaxon_TH2D;}

private:
    //double Vo, R, a;
    TF2 *WoodsSaxon_TF2, *transform_WoodsSaxon_TF2;
    TH2D *WoodsSaxon_TH2D, *transform_WoodsSaxon_TH2D;
};

#endif

