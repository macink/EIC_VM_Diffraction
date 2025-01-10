#ifndef WoodsSaxon_2D_H
#define WoodsSaxon_2D_H

#include "EICvalueconst.h"

using namespace std;



/************************************************ 
 
 *  Woods-Saxon Dentsity Distribution: G(x,y) [1/fm^-3]
 
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
        WS2 = new TF2("Woods-Saxon 2D",calc_WoodsSaxon2D,x_min,x_max,y_min,y_max,3);
        WS2->SetParameters(Vo,R,a);
        WS2->SetParNames("Vo","R","a");

        // Initialize Hist  
        WS3 = new TH2D("wsHist","Woods-Saxon Hist",bins,x_min,x_max,bins,y_min,y_max);
        WS3->Sumw2();
        for(int i=1;i<=bins;i++) 
        {
            for(int j=1;j<=bins;j++) 
            {
            double x = WS3->GetXaxis()->GetBinCenter(i);
            double y = WS3->GetYaxis()->GetBinCenter(j);
            WS3->SetBinContent(i,j,WS2->Eval(x, y)); 
            }
        }
        WS3->SetTitle("Woods-Saxon 2D: G(x,y)");
        WS3->GetYaxis()->SetTitle("y [fm]");
        WS3->GetXaxis()->SetTitle("x [fm]");
        WS3->GetZaxis()->SetTitle("G(x,y) [fm^{-3}]");

        // Initialize hist for transform: G(x,y) -> |F(tx,ty)|^2
        transform = new TH2D("hankelTransformWS","Transform Woods-Saxon",bins,tx_min,tx_max,bins,ty_min,ty_max);
        transform->Sumw2();
        double step = tx_max/bins;
        for(int i=0;i<bins;i++)
        {
            double tx = (i+1)*step, qx = sqrt(tx);
            for(int j=0;j<bins;j++)
            {
            double ty = (j+1)*step, qy = sqrt(ty);
	        double val = F2_integral_dr(qx,qy);
            double result = val*val;
            transform->SetBinContent(i+1,j+1,result);
            }
            //cout << "tx:" << tx << endl;
        }
	    transform->SetTitle("Transformed WS: G (x,y) #rightarrow |F(tx,ty)|^{2}");
        transform->GetYaxis()->SetTitle("ty [GeV^{2}]");
	    transform->GetXaxis()->SetTitle("tx [GeV^{2}]");
        transform->GetYaxis()->SetTitle("|F(tx,ty)|^{2}");
    }


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
        double qx = par[0], qy = par[1];
        double x = var[0], y = var[1];
        double r = sqrt(x*x+y*y);
        double q = sqrt(qx*qx+qy*qy);
        return 2*pi*sqrt(2*pi/q)*sqrt(2*hbarc/(pi*q*r))*sin(q *r/hbarc)*calc_WS2D(x,y)*sqrt(hbarc)*sqrt(r)*r; //[1/fm]
    }
    static double F2_integral_dr(double qx, double qy)
    {
        double r_min = -30., r_max = 30;
        TF1 *f = new TF1("f",F2_integrand_dr,r_min,r_max,2);  
        f->SetParameters(qx,qy);
        return f->Integral(r_min,r_max,1e-12); // [-]
    }



    TF2 *getWoodsSaxon2D() const {return WS2;}
    TH2D *getWoodsSaxonHist() const {return WS3;}
    TH2D *getWStransformHist() const {return transform;}

private:
    double Vo, R, a;
    TF2 *WS2;
    TH2D *WS3, *transform;
};

#endif

