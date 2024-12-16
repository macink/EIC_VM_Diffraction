#include <iostream>
#include "TF2.h"
#include "EICvalueconst.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"
#include <cmath>


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
class WoodsSaxon
{
public:
    // Calculations
    static double calc_WoodsSaxon1D(double *var, double *par)
    {
        double Vo = par[0], R = par[1], a = par[2]; 
        double r = var[0]; 
        return Vo/(1. + exp((r-R)/a)); 
    }
    static double calc_WoodsSaxon2D(double *var, double *par)
    {
        double Vo = par[0], R = par[1], a = par[2]; 
        double x = var[0], y = var[1];
        double r = sqrt(x*x + y*y); 
        return Vo/(1. + exp((r-R)/a)); // [1/fm^3]
    }

    // Functions
    TF1 *WoodsSaxon1D(double r_min, double r_max)
    {
        TF1 *WS1 = new TF1("Woods-Saxon 1D", calc_WoodsSaxon1D, r_min, r_max, 3);
        WS1->SetParameters(Vo, R, a);
        WS1->SetParNames("Vo", "R", "a");
        return WS1;
    }
    TF2 *WoodsSaxon2D(double x_min, double x_max, double y_min, double y_max)
    {
        TF2 *WS2 = new TF2("Woods-Saxon 2D", calc_WoodsSaxon2D, x_min, x_max, y_min, y_max, 3);
        WS2->SetParameters(Vo, R, a);
        WS2->SetParNames("Vo", "R", "a");
        return WS2;
    }

    // Histogram
    TH2D *create_WoodsSaxonHist2D(const char *name, const char *title, int xBins, double x_min, double x_max, int yBins, double y_min, double y_max)
    {
        TF2 *WS = WoodsSaxon2D(x_min,x_max,y_min,y_max);
        // Create a 2D histogram to represent the function
        TH2D *hist2D = new TH2D(name, title, xBins, x_min, x_max, yBins, y_min, y_max);
        hist2D->Sumw2();
        for (int i=1;i<=xBins;i++) 
        {
            for (int j=1;j<=yBins;j++) 
            {
                double x = hist2D->GetXaxis()->GetBinCenter(i);
                double y = hist2D->GetYaxis()->GetBinCenter(j);
                hist2D->SetBinContent(i, j, WS->Eval(x, y));
            }
            cout << "bin: " << i << endl;
        }
        hist2D->SetTitle("Woods-Saxon 2D: G(x,y)");
        hist2D->GetYaxis()->SetTitle("y [fm]");
	    hist2D->GetXaxis()->SetTitle("x [fm]");
        hist2D->GetZaxis()->SetTitle("G(x,y) [fm^{-3}]");
        return hist2D;
    }

private:
    double Vo = 2.12, R = 6.38, a = depth_Au;
};


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
    // Calculations
    static double calc_FormFactor_q1D(double *var, double *par) // F(q) [-]
    {
        double A = par[0], Vo = par[1], R = par[2], a0 = par[3];
        double q = var[0];
	    const double arg1 = q*R / hbarc;  
	    const double arg2 = hbarc / q;
	    const double arg3  = (sin(arg1)-arg1*cos(arg1)) * 4*pi*Vo*arg2*arg2*arg2 / A; 
        double result = arg3 / (1. + (a0 * a0 * q*q/(hbarc*hbarc)));
        cout << "q:  " << q << "  result: " << result << "\n";
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


    // Functions
    TF1 *get_FormFactor_q1D(double q_min, double q_max)
    {
        TF1 *FF1 = new TF1("Form Factor", calc_FormFactor_q1D, q_min, q_max, 4);
        FF1->SetParameters(197,2.12,6.38,0.70);
        FF1->SetParNames("A","Vo","R","a0");
        return FF1;
    }
    TF1 *get_FormFactor_2q1D(double q_min, double q_max)
    {
        TF1 *FF2 = new TF1("Form Factor", calc_FormFactor_2q1D, q_min, q_max, 4);
        FF2->SetParameters(197,2.12,6.38,0.70);
        FF2->SetParNames("A","Vo","R","a0");
        return FF2;
    }
    TF1 *get_FormFactor_t1D(double t_min, double t_max) 
    {
        TF1 *FF1 = new TF1("Form Factor", calc_FormFactor_t1D, t_min, t_max, 4);
        FF1->SetParameters(197,2.12,6.38,0.70);
        FF1->SetParNames("A", "Vo", "R", "a0");
        return FF1;
    }
    TF1 *get_FormFactor_2t1D(double t_min, double t_max)
    {
        TF1 *FF2 = new TF1("Form Factor", calc_FormFactor_2t1D, t_min, t_max, 4);
        FF2->SetParameters(197,2.12,6.38,0.70);
        FF2->SetParNames("A","Vo","R","a0");
        return FF2;
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


    // Wedge functions
    TF1 *qWedge_integral(double q_min, double q_max, double phi_max) 
    {
        double phi_min = 0;
        TF1 *cutFF = new TF1("Cut Form Factor", [phi_min, phi_max](double *var, double *par) {double q = var[0];  // q is variable
            TF1 f("f", calc_FormFactor_qWedge, phi_min, phi_max, 5);
            f.SetParameters(197, 2.12, 6.38, 0.70,q);
            return f.Integral(phi_min, phi_max, 1e-12);
        }, q_min, q_max, 0);
        return cutFF;
    }
    TF1 *q2Wedge_integral(double q_min, double q_max, double phi_max) 
    {
        double phi_min = 0;
        TF1 *cutFF = new TF1("Cut Form Factor", [phi_min, phi_max](double *var, double *par) {double q = var[0];  // q is variable
            TF1 f("f", calc_FormFactor_qWedge, phi_min, phi_max, 5);
            f.SetParameters(197, 2.12, 6.38, 0.70,q);
            return f.Integral(phi_min, phi_max, 1e-12);
        }, q_min, q_max, 0);

        TF1* cutFF2 = new TF1("Cut Form Factor Squared", [cutFF](double* var, double* par) {
            double q = var[0];
            return cutFF->Eval(q) * cutFF->Eval(q);
        }, q_min, q_max, 0);

        return cutFF2;
    }
    TF1 *tWedge_integral(double t_min, double t_max, double phi_max) 
    {
        double phi_min = 0;
        TF1 *cutFF = new TF1("Cut Form Factor", [phi_min, phi_max](double *var, double *par) {double t = var[0];  // t is variable
            TF1 f("f", calc_FormFactor_tWedge, phi_min, phi_max, 5);
            f.SetParameters(197, 2.12, 6.38, 0.70,t);
            return f.Integral(phi_min, phi_max, 1e-12);
        }, t_min, t_max, 0);
        return cutFF;
    }
    TF1 *t2Wedge_integral(double t_min, double t_max, double phi_max) 
    {
        double phi_min = 0;
        TF1 *cutFF = new TF1("Cut Form Factor", [phi_min, phi_max](double *var, double *par) {double t = var[0];  // t is variable
            TF1 f("f", calc_FormFactor_tWedge, phi_min, phi_max, 5);
            f.SetParameters(197, 2.12, 6.38, 0.70,t);
            return f.Integral(phi_min, phi_max, 1e-12);
        }, t_min, t_max, 0);

        TF1* cutFF2 = new TF1("Cut Form Factor Squared", [cutFF](double* var, double* par) {
            double t = var[0];
            return cutFF->Eval(t) * cutFF->Eval(t);
        }, t_min, t_max, 0);

        return cutFF2;
    }


    // Add resolution
    static double calc_FormFactor(double q)
    {
	    if(q==0){return 0;}
        else{
        double A = 197, rho0 = 2.12, R = 6.38, a0 = 0.7;
	    const double arg1 = q * R / hbarc;  
	    const double arg2 = hbarc /q;
	    const double sph  = (sin(arg1) - arg1 * cos(arg1)) * 4*pi*rho0 * arg2 * arg2 *arg2/ double(A); 
        const double result = sph * 1/(1. + (a0 * a0 * q*q/(hbarc*hbarc)));  // [-]
	    return result;
        }
    }
    static double combine_guassian_formFactor(double *var, double *par)
    {
        double sigma = 0.05;  // [GeV]
        double qy = par[0], qx_prime = par[1];
        double qx = var[0];  // integration variable
        double q = sqrt(qx*qx+qy*qy);
        double combined_funcs = calc_FormFactor(q)*exp(-(qx-qx_prime)*(qx-qx_prime)/(2*sigma*sigma));
	    return combined_funcs;
    }
    static double integrate_combined_funcs(double qy, double qx_prime)
    {
        double qx_min = 0, qx_max = 5;
        TF1 *new_ff = new TF1("new_ff",combine_guassian_formFactor,qx_min,qx_max,2);  
        new_ff->SetParameters(qy,qx_prime);
        return new_ff->Integral(qx_min,qx_max,1e-12); // F(qx',qy) [GeV]
    }
    static double FormFactor_wResolution_squared(double *var, double *par)
    {
        double qy = var[0], qx_prime = var[1];
        double smeared_formFactor = integrate_combined_funcs(qy, qx_prime)*integrate_combined_funcs(qy, qx_prime);
        return smeared_formFactor;
    }
    TF2 *define_new_formFactor()
    {
        TF2 *new_ff = new TF2("new_ff", FormFactor_wResolution_squared, 0, 5, 0, 5, 0);
        return new_ff;
    }
    
    TH2D *create_contour_q(const char *name, const char *title, int qy_bins, double qy_min, double qy_max, int qx_prime_bins, double qx_prime_min, double qx_prime_max)
    {
        TF2 *FF2 = new TF2("FF2", FormFactor_wResolution_squared, qx_prime_min, qx_prime_max, qy_min, qy_max, 0);
        // Create a 2D histogram to represent the function
        TH2D *hist2D = new TH2D(name, title, qy_bins, qy_min, qy_max, qx_prime_bins, qx_prime_min, qx_prime_max);
        hist2D->Sumw2();
        for (int i=1;i<=qy_bins;i++) 
        {
            for (int j=1;j<=qx_prime_bins;j++) 
            {
                double qy = hist2D->GetYaxis()->GetBinCenter(j);
                double qx_prime = hist2D->GetXaxis()->GetBinCenter(i);
                hist2D->SetBinContent(i, j, FF2->Eval(qx_prime, qy));
            }
            cout << "bin: " << i << endl;
        }
        hist2D->SetTitle("Form Factor 2D with Resolution: |F(qx',qy)|^{2}");
        hist2D->GetYaxis()->SetTitle("qy [GeV]");
	    hist2D->GetXaxis()->SetTitle("qx' [GeV]");
        hist2D->GetZaxis()->SetTitle("|F(qx',qy)|^{2}");
        return hist2D;
    }
    TH2D *create_contour_t(const char *name, const char *title, int ty_bins, double ty_min, double ty_max, int tx_prime_bins, double tx_prime_min, double tx_prime_max)
    {
        double qx_prime_min= sqrt(tx_prime_min), qx_prime_max = sqrt(tx_prime_max), qy_min = sqrt(ty_min), qy_max = sqrt(ty_max);
        TF2 *FF2 = new TF2("FF2", FormFactor_wResolution_squared, qx_prime_min, qx_prime_max, qy_min, qy_max, 0);
        // Create a 2D histogram to represent the function
        TH2D *hist2D = new TH2D(name, title, ty_bins, ty_min, ty_max, tx_prime_bins, tx_prime_min, tx_prime_max);
        hist2D->Sumw2();
        for (int i=1;i<=ty_bins;i++) 
        {
            for (int j=1;j<=tx_prime_bins;j++) 
            {
                double ty = hist2D->GetYaxis()->GetBinCenter(j);
                double tx_prime = hist2D->GetXaxis()->GetBinCenter(i);
                hist2D->SetBinContent(i, j, FF2->Eval(tx_prime, ty));
            }
            cout << "bin: " << i << endl;
        }
        hist2D->SetTitle("Form Factor 2D with Resolution: |F(tx,ty)|^{2}");
        hist2D->GetYaxis()->SetTitle("ty [GeV^{2}]");
	    hist2D->GetXaxis()->SetTitle("tx' [GeV^{2}]");
        hist2D->GetZaxis()->SetTitle("|F(tx,ty)|^{2}");
        return hist2D;
    }
    

private:
    double A = 197, Vo = 2.12, R = 6.38, a0 = 0.70;
};

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
    // Calculations
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

    // Functions
    TF2 *get_FormFactor_q2D(double qx_min, double qx_max, double qy_min, double qy_max)
    {
        TF2 *FF1 = new TF2("Form Factor", calc_FormFactor_q2D, qx_min, qx_max, qy_min, qy_max, 4);
        FF1->SetParameters(197,2.12,6.38,0.70);
        FF1->SetParNames("A","Vo","R","a0");
        return FF1;
    }
    TF2 *get_FormFactor_2q2D(double qx_min, double qx_max, double qy_min, double qy_max)
    {
        TF2 *FF2 = new TF2("Form Factor", calc_FormFactor_2q2D, qx_min, qx_max, qy_min, qy_max, 4);
        FF2->SetParameters(197,2.12,6.38,0.70);
        FF2->SetParNames("A","Vo","R","a0");
        return FF2;
    }
    TF2 *get_FormFactor_t2D(double tx_min, double tx_max, double ty_min, double ty_max)
    {
        TF2 *FF1 = new TF2("Form Factor", calc_FormFactor_t2D, tx_min, tx_max, ty_min, ty_max, 4);
        FF1->SetParameters(197,2.12,6.38,0.70);
        FF1->SetParNames("A","Vo","R","a0");
        return FF1;
    }
    TF2 *get_FormFactor_2t2D(double tx_min, double tx_max, double ty_min, double ty_max)
    {
        TF2 *FF2 = new TF2("Form Factor", calc_FormFactor_2t2D, tx_min, tx_max, ty_min, ty_max, 4);
        FF2->SetParameters(197,2.12,6.38,0.70);
        FF2->SetParNames("A","Vo","R","a0");
        return FF2;
    }


    // Histograms of |F(qx,qy)|^2 and |F(tx,ty)|^2
    TH2D *create_FormFactorq_2Dcontour(const char *name, const char *title, int qxBins, double qx_min, double qx_max, int qyBins, double qy_min, double qy_max)
    {
        TF2 *FF2 = get_FormFactor_2q2D(qx_min,qx_max,qy_min,qy_max);
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
            cout << "bin: " << i << endl;
        }
        hist2D->SetTitle("Form Factor 2D: F(qx,qy)");
        hist2D->GetYaxis()->SetTitle("qy [GeV]");
	    hist2D->GetXaxis()->SetTitle("qx [GeV]");
        hist2D->GetZaxis()->SetTitle("F(qx,qy)");
        return hist2D;
    }
    TH2D *create_FormFactort_2Dcontour(const char *name, const char *title, int txBins, double tx_min, double tx_max, int tyBins, double ty_min, double ty_max)
    {
        TF2 *FF2 = get_FormFactor_2t2D(tx_min,tx_max,ty_min,ty_max);
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
            cout << "bin: " << i << endl;
        }
        hist2D->SetTitle("Form Factor 2D: F(tx,ty)");
        hist2D->GetYaxis()->SetTitle("ty [GeV^{2}]");
	    hist2D->GetXaxis()->SetTitle("tx [GeV^{2}]");
        hist2D->GetZaxis()->SetTitle("F(tx,ty)");
        return hist2D;
    }


private:
    double A = 197, Vo = 2.12, R = 6.38, a0 = 0.70;
};


// test these functions
void test()
{
    WoodsSaxon ws;
    TF1 *WS1 = ws.WoodsSaxon1D(0,10);
    WS1->Draw();
}

void test2()
{
    WoodsSaxon ws;
    TF2 *WS2 = ws.WoodsSaxon2D(0,10,0,10);
    WS2->Draw();
}

void test3()
{
    WoodsSaxon ws;
    TH2D *WS_hist2D = ws.create_WoodsSaxonHist2D("2D hist", "Woods-Saxon Contour", 1000, 0, 10, 1000, 0, 10);
    WS_hist2D->Draw();
}

void test4()
{
    WoodsSaxon ws;
    TH2D *WS_hist2D = ws.create_WoodsSaxonHist2D("2D hist", "Woods-Saxon Contour: y-Projection", 1000, 0, 10, 1000, 0, 10);
    WS_hist2D->ProjectionY()->Draw();
}

void test5()
{
    FormFactor_1D ff;
    TF1 *test = ff.get_FormFactor_q1D(0,0.5);
    test->Draw();
}

void test6()
{
    FormFactor_1D ff;
    TF1 *test = ff.get_FormFactor_2q1D(0,0.5);
    test->Draw();
}

void test7()
{
    FormFactor_1D ff;
    TF1 *test = ff.get_FormFactor_2t1D(0,0.25);
    test->Draw();
    gPad->SetLogy(1);
}

void test8()
{
    FormFactor_1D ff;
    TF1 *test = ff.get_FormFactor_t1D(0,0.25);
    test->Draw();
}

void test9()
{
    FormFactor_2D ff;
    TF2 *test = ff.get_FormFactor_q2D(0,0.5,0,0.5);
    test->Draw();
}

void test10()
{
    FormFactor_2D ff;
    TF2 *test = ff.get_FormFactor_2q2D(0,0.5,0,0.5);
    test->Draw();
}

void test11()
{
    FormFactor_2D ff;
    TF2 *test = ff.get_FormFactor_2t2D(0,0.25,0,0.25);
    test->Draw();
}

void test12()
{
    FormFactor_2D ff;
    TF2 *test = ff.get_FormFactor_t2D(0,0.25,0,0.25);
    test->Draw();
}

void test13()
{
    FormFactor_2D ff;
    TH2D *FF_hist2D = ff.create_FormFactorq_2Dcontour("2D hist", "Form Factor Contour", 1000, 0, 0.5, 1000, 0, 0.5);
    FF_hist2D->Draw();
}

void test14()
{
    FormFactor_2D ff;
    TH2D *FF_hist2D = ff.create_FormFactort_2Dcontour("2D hist", "Form Factor Contour", 1000, 0, 0.25, 1000, 0, 0.25);
    FF_hist2D->Draw();
}

void test15()
{
    FormFactor_2D ff;
    TH2D *FF_hist2D = ff.create_FormFactort_2Dcontour("2D hist", "Form Factor Contour: y-Projection", 1000, 0, 0.25, 1000, 0, 0.25);
    FF_hist2D->ProjectionY()->Draw();
    gPad->SetLogy(1);
}


// test the wedges
void draw_qwedge()
{
    double q_min = 0, q_max = 0.5, phi_max = 2*pi;
    FormFactor_1D ff_wedge;

    TF1 *wedge = ff_wedge.qWedge_integral(q_min, q_max, phi_max);
    wedge->SetTitle("Form Factor Wedge Cut");
    wedge->GetYaxis()->SetTitle("F(q)");
    wedge->GetXaxis()->SetTitle("q [GeV]");
    wedge->Draw();
}

void compare_qwedge() 
{
    double q_min = 0, q_max = 0.5, phi_max = 2*pi;

    // Original form factor: F(q)
    FormFactor_1D ff;
    TF1 *FF1 = ff.get_FormFactor_q1D(q_min,q_max);
    FF1->Draw();

    // wedge integrated out
    FormFactor_1D ff_wedge;
    TF1 *FF2 = ff_wedge.qWedge_integral(q_min, q_max,phi_max);
    FF2->SetLineColor(kBlue);
    FF2->SetLineStyle(2);
    FF2->Draw("same");


    TLegend *legend = new TLegend(0.7, 0.8, 0.9, 0.9);
    legend->AddEntry(FF1, "no cut", "l");
    legend->AddEntry(FF2, "with wedge cut", "l");
    legend->Draw();
}

void draw_q2wedge()
{
    double q_min = 0, q_max = 0.5, phi_max = 2*pi;
    FormFactor_1D ff_wedge;
    TF1 *wedge = ff_wedge.q2Wedge_integral(q_min, q_max,phi_max);
    wedge->SetTitle("Form Factor Wedge Cut");
    wedge->GetYaxis()->SetTitle("|F(q)|^{2}");
    wedge->GetXaxis()->SetTitle("q [GeV^{q}]");
    wedge->Draw();
    gPad->SetLogy(1);
}

void compare_q2wedge() 
{
    double q_min = 0, q_max = 0.5, phi_max = 2*pi;

    // Original form factor: |F(q)|^2
    FormFactor_1D ff;
    TF1 *FF1 = ff.get_FormFactor_2q1D(q_min,q_max);
    FF1->Draw();

    // wedge integrated out
    FormFactor_1D ff_wedge;
    TF1 *FF2 = ff_wedge.q2Wedge_integral(q_min, q_max,phi_max);
    FF2->SetLineColor(kBlue);
    FF2->SetLineStyle(2);
    FF2->Draw("same");
    gPad->SetLogy(1);

    TLegend *legend = new TLegend(0.7, 0.8, 0.9, 0.9);
    legend->AddEntry(FF1, "no cut", "l");
    legend->AddEntry(FF2, "with wedge cut", "l");
    legend->Draw();
}


void draw_twedge()
{
    double t_min = 0, t_max = 0.25, phi_max = 2*pi;
    FormFactor_1D ff_wedge;

    TF1 *wedge = ff_wedge.tWedge_integral(t_min, t_max, phi_max);
    wedge->SetTitle("Form Factor Wedge Cut");
    wedge->GetYaxis()->SetTitle("|F(t)|^{2}");
    wedge->GetXaxis()->SetTitle("t [GeV^{2}]");
    wedge->Draw();
}

void compare_twedge() {
    double t_min = 0, t_max = 0.25, phi_max = 2*pi;

    // Original form factor: F(t)
    FormFactor_1D ff;
    TF1 *FF1 = ff.get_FormFactor_t1D(t_min,t_max);
    FF1->Draw();

    // wedge integrated out
    FormFactor_1D ff_wedge;
    TF1 *FF2 = ff_wedge.tWedge_integral(t_min,t_max,phi_max);
    FF2->SetLineColor(kBlue);
    FF2->SetLineStyle(2);
    FF2->Draw("same");

    TLegend *legend = new TLegend(0.7, 0.8, 0.9, 0.9);
    legend->AddEntry(FF1, "no cut", "l");
    legend->AddEntry(FF2, "with wedge cut", "l");
    legend->Draw();
}

void draw_t2wedge()
{
    double t_min = 0, t_max = 0.25, phi_max = 2*pi;
    FormFactor_1D ff_wedge;
    TF1 *wedge = ff_wedge.t2Wedge_integral(t_min, t_max,phi_max);
    wedge->SetTitle("Form Factor Wedge Cut");
    wedge->GetYaxis()->SetTitle("|F(t)|^{2}");
    wedge->GetXaxis()->SetTitle("t [GeV^{2}]");
    wedge->Draw();
    gPad->SetLogy(1);
}

void compare_t2wedge() 
{
    double t_min = 0, t_max = 0.25, phi_max = 2*pi;

    // Original form factor: |F(t)|^2
    FormFactor_1D ff;
    TF1 *FF1 = ff.get_FormFactor_2t1D(t_min,t_max);
    FF1->Draw();

    // wedge integrated out
    FormFactor_1D ff_wedge;
    TF1 *FF2 = ff_wedge.t2Wedge_integral(t_min, t_max,phi_max);
    FF2->SetLineColor(kBlue);
    FF2->SetLineStyle(2);
    FF2->Draw("same");
    gPad->SetLogy(1);

    TLegend *legend = new TLegend(0.7, 0.8, 0.9, 0.9);
    legend->AddEntry(FF1, "no cut", "l");
    legend->AddEntry(FF2, "with wedge cut", "l");
    legend->Draw();
}


void test16()
{
    FormFactor_1D test;
    TF2 *f2 = test.define_new_formFactor();
    f2->Draw();
}

void test17()
{
    FormFactor_1D test;
    TH2D *hist = test.create_contour_q("test", "test", 100, 0, 0.5, 100, 0 ,0.5);
    hist->Draw();
}

void test18()
{
    FormFactor_1D test;
    TH2D *hist = test.create_contour_q("test", "test", 100, 0, 0.5, 100, 0 ,0.5);
    TH1D *projY = hist->ProjectionY();
    projY->Draw();
}

void test19()
{
    FormFactor_1D test;
    TH2D *hist = test.create_contour_t("test", "test", 100, 0, 0.5, 100, 0 ,0.5);
    hist->Draw();
}

void test20()
{
    FormFactor_1D test;
    TH2D *hist = test.create_contour_t("test", "test", 1000, 0, 0.25, 1000, 0 ,0.25);
    TH1D *projY = hist->ProjectionY();
    projY->Draw();
}


