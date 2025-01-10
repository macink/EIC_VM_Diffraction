#ifndef FormFactor_saturation_data_H
#define FormFactor_saturation_data_H

#include "EICvalueconst.h"

using namespace std;


/********************************************** 
 
 * Form Factor saturation model: F(q) and F(t) [-]
 
 * Parameters:
    Vo: Nuclear density (hard sphere approximation) [fm^-3]
    A: Atomic mass number
    R: Radius of nucleus [fm]
    a0: Range of Yukawa potential [fm]
 * Variables:
    q: Transverse momentum >> q = sqrt(qx^2+qy^2) [GeV]
    t: Four-momentum transfer squared >> t = tx+ty = q^2 [GeV^2]

**********************************************/
class FormFactor_saturation_data
{
public:
    // Calculations to integrate out wedge (phi from qy axis)
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
    static double calc_FormFactor_noRes_wCut(const vector<double>& var, const vector<double>& par) 
    {
        double t = var[0];
        double phi_min = par[4], phi_max = par[5];

        TF1 cutFFt2("Cut Form Factor: |F(t)|^{2}",[phi_min,phi_max](double* var, double* par){double t = var[0]; 
            TF1 fft2("fft",calc_FormFactor2_tWedge,phi_min,phi_max,5);
            fft2.SetParameters(par[0],par[1],par[2],par[3],t);
            double integral = fft2.Integral(phi_min,phi_max,1e-12);
            return integral;
        },par[0],par[1],4);

        cutFFt2.SetParameters(par[0],par[1],par[2],par[3]);
        cutFFt2.SetParNames("A","Vo","R","a0");

        double result = cutFFt2.Eval(t);
        return result;
    }


    // Calculations to add resolution to form factor
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
    static double FormFactort_wResolution_squared(double *var, double *par)
    {
        double tx = var[0], ty = var[1];
        double qy = sqrt(ty), qx = sqrt(tx);
        double smeared_formFactor = integrate_combined_funcs(qy, qx)*integrate_combined_funcs(qy, qx);
        return smeared_formFactor;
    }
    static double calc_FormFactor_wResolution_noCut_data(const vector<double>& var, const vector<double>& par) 
    {
        double tx = var[0], ty = var[1];
        TF2 smearedFFt("Smeared Form Factor: |F(t)|^{2}",FormFactort_wResolution_squared,par[4],par[5],par[6],par[7],0);

        smearedFFt.SetParameters(par[0],par[1],par[2],par[3]);
        smearedFFt.SetParNames("A","Vo","R","a0");

        double result = smearedFFt.Eval(tx);
        return result;
    }
    static double calc_FormFactor_wResolution_WedgeCut_data(const vector<double>& var, const vector<double>& par) 
    {
        double tx = var[0], ty = var[1];
        double phi_min = par[4], phi_max = par[5];
        double t = tx+ty;
        // Initialize |F(t)|^2 with resolution TF2
        TF2 smearedFFq("Smeared Form Factor: |F(t)|^{2}",FormFactort_wResolution_squared,par[6],par[7],par[8],par[9], 0);

        // Initialize |F(t)|^2 wedge cut TF1    
        TF1 cutFFq2("Cut Form Factor: |F(t)|^{2}",[phi_min, phi_max](double* var, double* par){double t = var[0]; 
            TF1 ffq2("ffq",calc_FormFactor2_tWedge,phi_min,phi_max,5);
            ffq2.SetParameters(par[0],par[1],par[2],par[3],t);
            double integral = ffq2.Integral(phi_min,phi_max,1e-12);
            return integral;
        },par[0],par[1],4);

        cutFFq2.SetParameters(par[0],par[1],par[2],par[3]);
        cutFFq2.SetParNames("A","Vo","R","a0");

        // Calculate the form factor with resolution and wedge cut
        double result = cutFFq2.Eval(t);
        return result;
    }


    // Create plots
    TGraph *createGraph_data_noRes_noCut(const char *title, const char *x_axis, const char *y_axis, int bins, double* x_vals, double* y_vals) 
    {
        TCanvas *canvas = new TCanvas("canvas", "Graph from Data", 800, 600);
        TGraph *graph = new TGraph(bins, x_vals, y_vals);
        graph->SetTitle(title);
        graph->GetXaxis()->SetTitle("|t| [GeV^{2}]");
        graph->GetYaxis()->SetTitle("d#sigma/dt [nb/GeV^{2}]");
        graph->SetMarkerStyle(20);
        graph->Draw("AP");
        gPad->SetLogy(1);
        canvas->Update();
        return graph;
    }
    void createGraph_data_wRes_wCut(const vector<double>& var, const vector<double>& par) 
    {
        TGraph *graph = new TGraph();
        for (size_t i=0; i<var.size(); ++i) 
        {
            double result = calc_FormFactor_wResolution_WedgeCut_data({var[i]},par);
            graph->SetPoint(i,var[i],result);
        }
        TCanvas *c1 = new TCanvas("c1","Form Factor with Resolution and Wedge Cut",800,600);
        graph->SetTitle("Form Factor with Resolution and Wedge Cut; t; |F(t)|^{2}");
        graph->GetXaxis()->SetTitle("|t| [GeV^{2}]");
        graph->GetYaxis()->SetTitle("d#sigma/dt [nb/GeV^{2}]");
        gPad->SetLogy(1);
        graph->Draw("AL");
    }
    void createGraph_data_wRes_noCut(const vector<double>& var, const vector<double>& par) 
    {
        TGraph *graph = new TGraph();
        for (size_t i=0; i<var.size(); ++i) 
        {
            double result = calc_FormFactor_wResolution_noCut_data({var[i]},par);
            graph->SetPoint(i,var[i],result);
        }
        TCanvas *c1 = new TCanvas("c1","Form Factor with Resolution, No Cut",800,600);
        graph->SetTitle("Form Factor with Resolution, No Cut; t; |F(t)|^{2}");
        graph->GetXaxis()->SetTitle("|t| [GeV^{2}]");
        graph->GetYaxis()->SetTitle("d#sigma/dt [nb/GeV^{2}]");
        gPad->SetLogy(1);
        graph->Draw("AL");
    }
    void createGraph_data_noRes_wCut(const vector<double>& var, const vector<double>& par) 
    {
        TGraph *graph = new TGraph();
        for (size_t i=0; i<var.size(); ++i) 
        {
            double result = calc_FormFactor_noRes_wCut({var[i]},par);
            graph->SetPoint(i,var[i],result);
        }

        TCanvas *c1 = new TCanvas("c1","Form Factor with Cut, No Resolution",800,600);
        graph->SetTitle("Form Factor with Cut, No Resolution; t; |F(t)|^{2}");
        graph->GetXaxis()->SetTitle("|t| [GeV^{2}]");
        graph->GetYaxis()->SetTitle("d#sigma/dt [nb/GeV^{2}]");
        gPad->SetLogy(1);
        graph->Draw("AL");
    }
    

private:
    double A, Vo, R, a0;
    vector<double> x_vals_vec;
    vector<double> y_vals_vec;
};




#endif
