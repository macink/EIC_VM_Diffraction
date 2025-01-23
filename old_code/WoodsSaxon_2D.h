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

    WoodsSaxon_2D(double Vo_init, double R_init, double a_init, double x_min, double x_max, double y_min, double y_max, 
        double bins, double tx_min, double tx_max, double ty_min, double ty_max): Vo(Vo_init), R(R_init), a(a_init)
    {
        ROOT::Math::IntegratorOneDimOptions::SetDefaultAbsTolerance(1.E-6);
        ROOT::Math::IntegratorOneDimOptions::SetDefaultRelTolerance(1.E-6);
        ROOT::Math::IntegratorOneDimOptions::SetDefaultIntegrator("Gauss");

        // Initialize TF2
        WS2 = new TF2("", calc_WoodsSaxon2D, x_min, x_max, y_min, y_max, 3);
        WS2->SetParameters(Vo, R, a);

        // Initialize Hist for G(x,y)
        WS3 = new TH2D("", "Woods-Saxon Hist", bins, x_min, x_max, bins, y_min, y_max);
        WS3->Sumw2();
        for(int i=1;i<=bins;i++) 
        {
            for(int j=1;j<=bins;j++) 
            {
            double x = WS3->GetXaxis()->GetBinCenter(i);
            double y = WS3->GetYaxis()->GetBinCenter(j);
            WS3->SetBinContent(i, j, WS2->Eval(x, y)); 
            }
        }
        WS3->GetYaxis()->SetTitle("y [fm]");
        WS3->GetXaxis()->SetTitle("x [fm]");

        // Initialize TF2 for Transformation: G(x,y) -> |F(tx,ty)|^2
        transform_TF2 = new TF2("", [this, x_min, x_max, y_min, y_max] (double *var, double *par)
        {
            double rr_min = 0., rr_max = 100;
            double tx = var[0], ty = var[1];
            TF1 transform_integral("", F2_integrand_dr, rr_min, rr_max, 2);
            transform_integral.SetParameters(tx, ty);
            return transform_integral.Integral(rr_min, rr_max, 1e-12)*transform_integral.Integral(rr_min, rr_max, 1e-12);
        }, tx_min, tx_max, ty_min, ty_max,0);  

        // Initialize hist for transform: G(x,y) -> |F(tx,ty)|^2
        transform_hist = new TH2D("", "Transform Woods-Saxon ",bins ,tx_min, tx_max, bins, ty_min, ty_max);
        transform_hist->Sumw2();
        double step = tx_max/bins;
        for(int i=0; i<bins; i++)
        {
            double tx = (i+1)*step;
            for(int j=0; j<bins; j++)
            {
                double ty = (j+1)*step;
	            double val = F2_integral_dr(tx, ty);
                double result = val*val;
                transform_hist->SetBinContent(i+1, j+1, result);
            }
        }
        transform_hist->GetYaxis()->SetTitle("t_{y} [GeV^{2}]");
	    transform_hist->GetXaxis()->SetTitle("t_{x} [GeV^{2}]");
    }


    // Calculations for Woods-Saxon density distribution and F integrand
    static double calc_WoodsSaxon2D(double *var, double *par)
    {
        double Vo = par[0], R = par[1], a = par[2]; 
        double x = var[0], y = var[1];
        double r = sqrt(x*x + y*y); 
        return Vo/(1. + exp((r-R)/a)); // [1/fm^3]
    }
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
    static double F2_integral_dr(double tx, double ty)
    {
        double rr_min = 0., rr_max = 100;
        TF1 *f = new TF1("", F2_integrand_dr, rr_min, rr_max, 2);  
        f->SetParameters(tx,ty);
        return f->Integral(rr_min, rr_max, 1e-12); // [-]
    }


    // Functions
    TF2 *getWoodsSaxon2D() const {return WS2;}
    TH2D *getWoodsSaxonHist() const {return WS3;}
    TF2 *getWStransform_fun() const {return transform_TF2;}
    TH2D *getTransformed_hist() const {return transform_hist;}

private:
    double Vo, R, a;
    TF2 *WS2, *transform_TF2;
    TH2D *WS3, *transform_hist;
};

#endif

