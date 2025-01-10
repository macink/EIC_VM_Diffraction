#ifndef WoodsSaxon_1D_H
#define WoodsSaxon_1D_H

#include "EICvalueconst.h"

using namespace std;



/************************************************ 
 
 *  Woods-Saxon Dentsity Distribution: G(r) [1/fm^-3]
 
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
        WS1 = new TF1("Woods-Saxon 1D",calc_WoodsSaxon1D,r_min,r_max,3);
        WS1->SetParameters(Vo,R,a);
        WS1->SetParNames("Vo","R","a");

        // Initialize TF1 for Transformation: G(r) -> |F(t)|^2
        transform1 = new TF1("Fourier-Bessel Transformation: G(r) -> |F(t)|^{2}",[this,r_min,r_max](double *var, double *par){double t = var[0];
            TF1 trans_integral("trans_integral", F2_integrand_dr,r_min,r_max,1);
            trans_integral.SetParameters(t);
            return trans_integral.Integral(r_min,r_max,1e-12)*trans_integral.Integral(r_min,r_max,1e-12);
        },t_min,t_max,0);  

        // Initialize hist for transform: G(r) -> |F(t)|^2
        transform = new TH1D("hankelTransformWS","Transform Woods-Saxon",bins,t_min,t_max);
        transform->Sumw2();
        double step = t_max/bins;
        for(int i=0;i<bins;i++)
        {
            double t = (i+1)*step;
	        double val = F2_integral_dr(t);
            double result = val*val;
            transform->SetBinContent(i+1,result);
            //cout << "t:" << t << endl;
        }
	    transform->SetTitle("Transformed WS: G_{true} (r) #rightarrow |F(t)|^{2}");
        transform->GetYaxis()->SetTitle("|F(t)|^{2}");
	    transform->GetXaxis()->SetTitle("t [GeV^{2}]");

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
    static double F2_integrand_dr(double *var, double *par)  // G(r) -> |F(t)|^2
    {
        double t = par[0];
        double r = var[0];
        double q = sqrt(t);
        return 2*pi*sqrt(2*pi/q)*sqrt(2*hbarc/(pi*q*r))*sin(q *r/hbarc)*calc_WS1D(r)*sqrt(hbarc)*sqrt(r)*r; //[1/fm]
    }
    static double F2_integral_dr(double t)
    {
        double r_min = 0., r_max = 30;
        TF1 *f = new TF1("f",F2_integrand_dr,r_min,r_max,1);  
        f->SetParameters(t);
        return f->Integral(r_min,r_max,1e-12); // [-]
    }


    TF1 *getWoodsSaxon1D() const {return WS1;}
    TF1 *getWStransform() const {return transform1;}
    TH1D *getWStransformHist() const {return transform;}


private:
    double Vo, R, a;
    TF1 *WS1, *transform1;
    TH1D *transform;
};

#endif
