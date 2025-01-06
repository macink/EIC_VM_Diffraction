#include <iostream>
//#include "TComplex.h"
//#include "TVector2.h"
#include "TF2.h"
//#include "TVirtualFFT.h"
#include "EICvalueconst.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"
//#include <Riostream.h>
//#include "TLegend.h"
//#include "TLegendEntry.h"
//#include "Math/IFunction.h"
//#include <cmath>
//#include "TSystem.h"
//#include "TAxis.h"
//#include "TPaveLabel.h"
//#include "TFormula.h"
//#include "TRandom.h"
//#include "TGraph.h"



using namespace std;

/************************************************ 
 
 *  Woods-Saxon Dentsity Distribution: G(r)
 
 * Parameters:
    Vo: Nuclear density (hard sphere approximation)
    a: Skin-depth
    R: Radius of nucleus
 * Variables:
    r: Radial distance from center of nucleus
 
**************************************************/ 
double calc_WS1D(double *var, double *par)
{
    double Vo = par[0], a = par[1], R = par[2]; 
    double r = var[0]; 
    return Vo/(1. + exp((r-R)/a)); // [1/fm^3]
}

TF1 *Woods_Saxon1D(double r_min,double r_max) 
{
    TF1 *WS1 = new TF1("Woods-Saxon 1D", "calc_WS1D", r_min, r_max, 3);
    WS1->SetParameters(2.12,depth_Au,6.38);
    WS1->SetParNames("Vo","a","R");
    return WS1;
}

void draw_WS_1D()
{
    double r_min = 0, r_max = 10;
    TF1 *WS1 = Woods_Saxon1D(r_min,r_max);
    
    TCanvas *c1 = new TCanvas("c1", "Woods-Saxon", 800, 600);
    WS1->SetTitle("Woods-Saxon 1D");
    WS1->GetYaxis()->SetTitle("G(r) [fm^{-3}]");
	WS1->GetXaxis()->SetTitle("r [fm]");
    WS1->Draw();

    TFile *ws = new TFile("WS_1D.root","recreate");
    WS1->Write("WS1");
    ws->Close();
    c1->SaveAs("WS1D.png");
}

double calc_WS2D(double *var, double *par)
{
    double Vo = par[0], a = par[1], R = par[2]; 
    double x = var[0], y = var[1]; 
    double r = sqrt(x*x + y*y);
    return Vo/(1. + exp((r-R)/a)); // [1/fm^3]
}

TF2 *Woods_Saxon2D(double x_min,double x_max, double y_min, double y_max) 
{
    TF2 *WS2 = new TF2("Woods-Saxon 2D", "calc_WS2D", x_min, x_max, y_min, y_max, 3);
    WS2->SetParameters(2.12,depth_Au,6.38);
    WS2->SetParNames("Vo","a","R");
    return WS2;
}

void draw_WS_2D()
{
    double x_min = 0, x_max = 10, y_min = 0, y_max = 10;
    TF2 *WS2 = Woods_Saxon2D(x_min,x_max,y_min,y_max);

    TCanvas *c1 = new TCanvas("c1", "Woods-Saxon", 800, 600);
    WS2->SetTitle("Woods-Saxon 2D: G(x,y)");
    WS2->GetYaxis()->SetTitle("y [fm]");
	WS2->GetXaxis()->SetTitle("x [fm]");
    WS2->Draw();

    TFile *ws = new TFile("WS_2D.root","recreate");
    WS2->Write("WS2");
    ws->Close();
    c1->SaveAs("WS2D.png");

}

void draw_WS_contour()
{
    double x_min = 0, x_max = 10, y_min = 0, y_max = 10;
    TF2 *WS2 = Woods_Saxon2D(x_min,x_max,y_min,y_max);

    // Create a 2D histogram to represent the function
    int bins = 1000; 
    TH2D *hist2D = new TH2D("hist2D", "Woods-Saxon Potential 2D", bins, x_min, x_max, bins, y_min, y_max);
    hist2D->Sumw2();
    for (int i=1;i<=bins;i++) 
    {
        for (int j=1;j<=bins;j++) 
        {
            double x = hist2D->GetXaxis()->GetBinCenter(i);
            double y = hist2D->GetYaxis()->GetBinCenter(j);
            hist2D->SetBinContent(i, j, WS2->Eval(x, y));
            cout << "y: " << y << endl;
        }
    }
    TFile *WS = new TFile("WS_2D_contour.root","recreate");
    hist2D->SetTitle("Woods-Saxon 2D: G(x,y)");
    hist2D->GetYaxis()->SetTitle("y [fm]");
	hist2D->GetXaxis()->SetTitle("x [fm]");
    hist2D->Write();
    WS->Close();

    TCanvas *c1 = new TCanvas("Woods-Saxon", "Woods-Saxon 2D", 800, 600);
    hist2D->Draw();
    gStyle->SetOptStat(0); // remove stat box on histogram
    c1->SaveAs("WS2D_contour.png");
}

void draw_yproj_WS2D()
{
    TFile *f = new TFile("WS_2D_contour.root");
    TH2D *get = (TH2D *)f->Get("hist2D");
    TH1D *get1 = (TH1D *)get->ProjectionY();

    TCanvas *c1 = new TCanvas("Woods-Saxon", "Woods-Saxon 2D", 800, 600);
    get1->SetTitle("Woods-Saxon");
    get1->GetXaxis()->SetTitle("y [fm]");
    get1->GetYaxis()->SetTitle("G(x,y)"); 
    get1->Draw();
    gStyle->SetOptStat(0); // remove stat box on histogram
    c1->SaveAs("WS2D_yproj.png");
}


/********************************************** 
 
 * Form Factor: F(q)
 
 * Parameters:
    Vo: Nuclear density (hard sphere approximation)
    A: Atomic mass number
    R: Radius of nucleus
    a0: Range of Yukawa potential
 * Variables:
    q: Transverse momentum >> q = sqrt(qx^2+qy^2)
**********************************************/
double calc_FFq1D(double *var, double *par)
{
    double A = par[0], Vo = par[1], R = par[2], a0 = par[3];
    double q = var[0];
	const double arg1 = q*R / hbarc;  
	const double arg2 = hbarc / q;
	const double arg3  = (sin(arg1)-arg1*cos(arg1)) * 4*pi*Vo*arg2*arg2*arg2 / A; 
	return arg3 / (1. + (a0 * a0 * q*q/(hbarc*hbarc)));  // [-]
}

TF1 *form_factorq1D(double q_min, double q_max)
{
    TF1 *FF1 = new TF1("Form Factor", "calc_FFq1D", q_min, q_max, 4);
    FF1->SetParameters(197,2.12,6.38,0.70);
    FF1->SetParNames("A","Vo","R","a0");
    return FF1;
}

void draw_FFq_1D()
{
    double q_min = 0, q_max = 0.5;

    TCanvas *c1 = new TCanvas("c1", "Form Factor", 800, 600);
    TF1 *FF1 = form_factorq1D(q_min,q_max);
    FF1->SetTitle("Form Factor 1D");
    FF1->GetYaxis()->SetTitle("F(q)");
	FF1->GetXaxis()->SetTitle("q [GeV]");
    FF1->Draw();

    TFile *ff = new TFile("FFq_1D.root","recreate");
    FF1->Write("FF1");
    ff->Close();
}

double calc_FF2q1D(double *var, double *par)
{
    double A = par[0], Vo = par[1], R = par[2], a0 = par[3];
    double q = var[0];
	const double arg1 = q*R / hbarc;  
	const double arg2 = hbarc / q;
	const double arg3  = (sin(arg1)-arg1*cos(arg1)) * 4*pi*Vo*arg2*arg2*arg2 / A; 
	return arg3 / (1. + (a0 * a0 * q*q/(hbarc*hbarc)))*arg3 / (1. + (a0 * a0 * q*q/(hbarc*hbarc)));  // [-]
}

TF1 *form_factor2q1D(double q_min, double q_max)
{
    TF1 *FF2 = new TF1("Form Factor", "calc_FF2q1D", q_min, q_max, 4);
    FF2->SetParameters(197,2.12,6.38,0.70);
    FF2->SetParNames("A","Vo","R","a0");
    return FF2;
}

void draw_FF2q_1D()
{
    double q_min = 0, q_max = 0.5;

    TCanvas *c1 = new TCanvas("c1", "Form Factor", 800, 600);
    TFile *ff = new TFile("FF2q_1D.root","recreate");
    TF1 *FF2 = form_factor2q1D(q_min,q_max);
    FF2->SetTitle("Form Factor 1D");
    FF2->GetYaxis()->SetTitle("|F(q)|^{2}");
	FF2->GetXaxis()->SetTitle("q [GeV]");
    FF2->Draw();
    FF2->Write("FF2");
    ff->Close();
}

double calc_FFq2D(double *var, double *par)
{
    double A = par[0], Vo = par[1], R = par[2], a0 = par[3];
    double qx = var[0], qy = var[1];
    double q = sqrt(qx*qx+qy*qy);
	const double arg1 = q*R / hbarc;  
	const double arg2 = hbarc / q;
	const double arg3  = (sin(arg1)-arg1*cos(arg1)) * 4*pi*Vo*arg2*arg2*arg2 / A; 
	return arg3 / (1. + (a0 * a0 * q*q/(hbarc*hbarc)));  // [-]
}

TF2 *form_factorq2D(double qx_min, double qx_max, double qy_min, double qy_max)
{
    TF2 *FF1 = new TF2("Form Factor", "calc_FFq2D", qx_min, qx_max, qy_min, qy_max, 4);
    FF1->SetParameters(197,2.12,6.38,0.70);
    FF1->SetParNames("A","Vo","R","a0");
    return FF1;
}

void draw_FFq_2D()
{
    double qx_min = 0, qx_max = 0.5, qy_min = 0, qy_max = 0.5;

    TCanvas *c1 = new TCanvas("c1", "Form Factor", 800, 600);
    TF2 *FF1 = form_factorq2D(qx_min,qx_max,qy_min,qy_max);
    FF1->SetTitle("Form Factor 2D");
    FF1->GetYaxis()->SetTitle("q_{y} [GeV]");
	FF1->GetXaxis()->SetTitle("q_{x} [GeV]");
    FF1->Draw();

    TFile *ff = new TFile("FFq_2D.root","recreate");
    FF1->Write("FF1");
    ff->Close();
}

void draw_FFq_contour()
{
    double qx_min = 0, qx_max = 0.5, qy_min = 0, qy_max = 0.5;
    TF2 *FF2 = form_factorq2D(qx_min,qx_max,qy_min,qy_max);

    // Create a 2D histogram to represent the function
    int bins = 1000; 
    TH2D *hist2D = new TH2D("hist2D", "Form Factor 2D", bins, qx_min, qx_max, bins, qy_min, qy_max);
    hist2D->Sumw2();
    for (int i=1;i<=bins;i++) 
    {
        for (int j=1;j<=bins;j++) 
        {
            double qx = hist2D->GetXaxis()->GetBinCenter(i);
            double qy = hist2D->GetYaxis()->GetBinCenter(j);
            hist2D->SetBinContent(i, j, FF2->Eval(qx, qy));
            cout << "qx: " << qx << endl;
        }
    }
    TFile *FF = new TFile("FFq_2D_contour.root","recreate");
    hist2D->SetTitle("Form Factor 2D ");
    hist2D->GetYaxis()->SetTitle("q_{y} [GeV]");
	hist2D->GetXaxis()->SetTitle("q_{x} [GeV]");
    hist2D->Write();
    FF->Close();

    TCanvas *c1 = new TCanvas("Form Factor", "Form Factor 2D", 800, 600);
    hist2D->Draw();
    gStyle->SetOptStat(0);
}

void draw_yproj_FFq2D()
{
    TFile *f = new TFile("FFq_2D_contour.root");
    TH2D *get = (TH2D *)f->Get("hist2D");
    TH1D *get1 = (TH1D *)get->ProjectionY();
    get1->SetTitle("Form Factor");
    get1->GetXaxis()->SetTitle("q_{y} [GeV]");
    get1->GetYaxis()->SetTitle("F(q_{x},q_{y})"); 
    get1->Draw();
}

double calc_FF2q2D(double *var, double *par)
{
    double A = par[0], Vo = par[1], R = par[2], a0 = par[3];
    double qx = var[0], qy = var[1];
    double q = sqrt(qx*qx+qy*qy);
	const double arg1 = q*R / hbarc;  
	const double arg2 = hbarc / q;
	const double arg3  = (sin(arg1)-arg1*cos(arg1)) * 4*pi*Vo*arg2*arg2*arg2 / A; 
	return arg3 / (1. + (a0 * a0 * q*q/(hbarc*hbarc)))*arg3 / (1. + (a0 * a0 * q*q/(hbarc*hbarc)));  // [-]
}

TF2 *form_factor2q2D(double qx_min, double qx_max, double qy_min, double qy_max)
{
    TF2 *FF2 = new TF2("Form Factor", "calc_FF2q2D", qx_min, qx_max, qy_min, qy_max, 4);
    FF2->SetParameters(197,2.12,6.38,0.70);
    FF2->SetParNames("A","Vo","R","a0");
    return FF2;
}

void draw_FF2q_2D()
{
    double qx_min = 0, qx_max = 0.5, qy_min = 0, qy_max = 0.5;

    TCanvas *c1 = new TCanvas("c1", "Form Factor", 800, 600);
    TF2 *FF2 = form_factor2q2D(qx_min,qx_max,qy_min,qy_max);
    FF2->SetTitle("Form Factor 2D");
    FF2->GetYaxis()->SetTitle("q_{y} [GeV]");
	FF2->GetXaxis()->SetTitle("q_{x} [GeV]");
    FF2->Draw();

    TFile *ff = new TFile("FF2q_2D.root","recreate");
    FF2->Write("FF2");
    ff->Close();
}

void draw_FF2q_contour()
{
    double qx_min = 0, qx_max = 0.5, qy_min = 0, qy_max = 0.5;
    TF2 *FF2 = form_factor2q2D(qx_min,qx_max,qy_min,qy_max);

    // Create a 2D histogram to represent the function
    int bins = 1000; 
    TH2D *hist2D = new TH2D("hist2D", "Form Factor 2D", bins, qx_min, qx_max, bins, qy_min, qy_max);
    hist2D->Sumw2();
    for (int i=1;i<=bins;i++) 
    {
        for (int j=1;j<=bins;j++) 
        {
            double qx = hist2D->GetXaxis()->GetBinCenter(i);
            double qy = hist2D->GetYaxis()->GetBinCenter(j);
            hist2D->SetBinContent(i, j, FF2->Eval(qx, qy));
            cout << "qx: " << qx << endl;
        }
    }
    TFile *FF = new TFile("FF2q_2D_contour.root","recreate");
    hist2D->SetTitle("Form Factor 2D ");
    hist2D->GetYaxis()->SetTitle("q_{y} [GeV]");
	hist2D->GetXaxis()->SetTitle("q_{x} [GeV]");
    hist2D->Write();
    FF->Close();

    TCanvas *c1 = new TCanvas("Form Factor", "Form Factor 2D", 800, 600);
    hist2D->Draw();
    gStyle->SetOptStat(0);
}

void draw_yproj_FF2q2D()
{
    TFile *f = new TFile("FF2q_2D_contour.root");
    TH2D *get = (TH2D *)f->Get("hist2D");
    TH1D *get1 = (TH1D *)get->ProjectionY();
    get1->SetTitle("Form Factor");
    get1->GetXaxis()->SetTitle("q_{y} [GeV]");
    get1->GetYaxis()->SetTitle("|F(q_{x},q_{y})|^{2}"); 
    get1->Draw();
}


// Integrate out wedge (theta from qy axis)
double calc_FFq_wedge(double *var, double *par)
{
    double A = par[0], Vo = par[1], R = par[2], a0 = par[3], q = par[4];
    double theta = var[0];
    double qx = q*sin(theta), qy = q*cos(theta);
    double qq = sqrt(qx*qx+qy*qy);
	const double arg1 = qq*R / hbarc;  
	const double arg2 = hbarc / qq;
	const double arg3  = (sin(arg1)-arg1*cos(arg1)) * 4*pi*Vo*arg2*arg2*arg2 / A; 
    return 1/(2*pi)*arg3 / (1. + (a0 * a0 * qq*qq/(hbarc*hbarc)));  // [-]
}

TF1 *qwedge_inte(double q_min, double q_max) 
{
    double theta_min = 0, theta_max = 2*pi;
    TF1 *cutFF = new TF1("Cut Form Factor", [theta_min, theta_max](double *var, double *par) {double q = var[0];  // q is variable
        TF1 f("f", "calc_FFq_wedge", theta_min, theta_max, 5);
        f.SetParameters(197, 2.12, 6.38, 0.70,q);
        return f.Integral(theta_min, theta_max, 1e-12);
    }, q_min, q_max, 0);
    return cutFF;
}

void draw_qwedge()
{
    double q_min = 0, q_max = 0.5;
    TCanvas *c1 = new TCanvas("c1", "Form Factor", 800, 600);
    TF1 *FF2 = qwedge_inte(q_min, q_max);
    FF2->SetTitle("Form Factor Wedge Cut");
    FF2->GetYaxis()->SetTitle("|F(q)|^{2}");
    FF2->GetXaxis()->SetTitle("q [GeV]");
    FF2->Draw();

    TFile *ff = new TFile("wedge_q.root", "recreate");
    FF2->Write("FF2");
    ff->Close();
}

void compare_qwedge() 
{
    double q_min = 0, q_max = 0.5;
    TCanvas *c1 = new TCanvas("c1", "Form Factor", 800, 600);
    
    // Original form factor: F(q)
    //TF1 *FF1 = form_factorq1D(q_min,q_max);
    //TFile *f1 = TFile::Open("FFq_1D.root");
    TFile *f1 = new TFile("FFq_1D.root","READ");
    TF1 *FF1 = (TF1*)f1->Get("FF1");
    FF1->SetTitle("Form Factor 1D");
    FF1->GetYaxis()->SetTitle("F(q)");
	FF1->GetXaxis()->SetTitle("q [GeV]");
    FF1->SetLineColor(kRed);
    FF1->DrawCopy();

    // wedge integrated out: theta = 2pi
    TF1 *FF2 = qwedge_inte(q_min, q_max);
    FF2->SetLineColor(kBlue);
    FF2->SetLineStyle(2);
    FF2->Draw("same");


    TLegend *legend = new TLegend(0.7, 0.8, 0.9, 0.9);
    legend->AddEntry(FF1, "no cut", "l");
    legend->AddEntry(FF2, "with wedge cut", "l");
    legend->Draw();

    c1->SaveAs("compare_qwedges.png");
}

double calc_FF2q_wedge(double *var, double *par)
{
    double A = par[0], Vo = par[1], R = par[2], a0 = par[3], q = par[4];
    double theta = var[0];
    double qx = q*sin(theta), qy = q*cos(theta);
    double qq = sqrt(qx*qx+qy*qy);
	const double arg1 = qq*R / hbarc;  
	const double arg2 = hbarc / qq;
	const double arg3  = (sin(arg1)-arg1*cos(arg1)) * 4*pi*Vo*arg2*arg2*arg2 / A; 
    return 1/(2*pi)*arg3 / (1. + (a0 * a0 * qq*qq/(hbarc*hbarc)))*1/(2*pi)*arg3 / (1. + (a0 * a0 * qq*qq/(hbarc*hbarc)));  // [-]
}

TF1 *q2wedge_inte(double q_min, double q_max) 
{
    double theta_min = 0, theta_max = 2*pi;
    TF1 *cutFF = new TF1("Cut Form Factor", [theta_min, theta_max](double *var, double *par) {double q = var[0];  // q is variable
        TF1 f("f", "calc_FF2q_wedge", theta_min, theta_max, 5);
        f.SetParameters(197, 2.12, 6.38, 0.70,q);
        return f.Integral(theta_min, theta_max, 1e-12);
    }, q_min, q_max, 0);
    return cutFF;
}

void draw_q2wedge()
{
    double q_min = 0, q_max = 0.5;
    TCanvas *c1 = new TCanvas("c1", "Form Factor", 800, 600);
    TF1 *FF2 = q2wedge_inte(q_min, q_max);
    FF2->SetTitle("Form Factor Wedge Cut");
    FF2->GetYaxis()->SetTitle("|F(q)|^{2}");
    FF2->GetXaxis()->SetTitle("q [GeV^{q}]");
    FF2->Draw();

    TFile *ff = new TFile("wedge_q2.root", "recreate");
    FF2->Write("FF2");
    ff->Close();
}


/********************************************** 
 
 * Form Factor: F(t)
 
 * Parameters:
    Vo: Nuclear density (hard sphere approximation)
    A: Atomic mass number
    R: Radius of nucleus
    a0: Range of Yukawa potential
 * Variables:
    q: Transverse momentum >> q = sqrt(qx^2+qy^2)
    t: Four-momentum transfer squared >> t = tx+ty = q^2
**********************************************/
double calc_FFt1D(double *var, double *par)
{
    double A = par[0], Vo = par[1], R = par[2], a0 = par[3];
    double t = var[0];
    double q = sqrt(t);
	const double arg1 = q*R / hbarc;  
	const double arg2 = hbarc / q;
	const double arg3  = (sin(arg1)-arg1*cos(arg1)) * 4*pi*Vo*arg2*arg2*arg2 / A; 
	return arg3 / (1. + (a0 * a0 * q*q/(hbarc*hbarc)));  // [-]
}

TF1 *form_factort1D(double t_min, double t_max)
{
    TF1 *FF1 = new TF1("Form Factor", "calc_FFt1D", t_min, t_max, 4);
    FF1->SetParameters(197,2.12,6.38,0.70);
    FF1->SetParNames("A","Vo","R","a0");
    return FF1;
}

void draw_FFt_1D()
{
    double t_min = 0, t_max = 0.25;

    TCanvas *c1 = new TCanvas("c1", "Form Factor", 800, 600);
    TF1 *FF1 = form_factort1D(t_min,t_max);
    FF1->SetTitle("Form Factor 1D");
    FF1->GetYaxis()->SetTitle("F(t)");
	FF1->GetXaxis()->SetTitle("t [GeV^{2}]");
    FF1->Draw();

    TFile *ff = new TFile("FFt_1D.root","recreate");
    FF1->Write("FF1");
    ff->Close();
}

double calc_FF2t1D(double *var, double *par)
{
    double A = par[0], Vo = par[1], R = par[2], a0 = par[3];
    double t = var[0];
    double q = sqrt(t);
	const double arg1 = q*R / hbarc;  
	const double arg2 = hbarc / q;
	const double arg3  = (sin(arg1)-arg1*cos(arg1)) * 4*pi*Vo*arg2*arg2*arg2 / A; 
	return arg3 / (1. + (a0 * a0 * q*q/(hbarc*hbarc)))*arg3 / (1. + (a0 * a0 * q*q/(hbarc*hbarc)));  // [-]
}

TF1 *form_factor2t1D(double t_min, double t_max)
{
    TF1 *FF2 = new TF1("Form Factor", "calc_FF2t1D", t_min, t_max, 4);
    FF2->SetParameters(197,2.12,6.38,0.70);
    FF2->SetParNames("A","Vo","R","a0");
    return FF2;
}

void draw_FF2t_1D()
{
    double t_min = 0, t_max = 0.25;

    TCanvas *c1 = new TCanvas("c1", "Form Factor", 800, 600);
    TF1 *FF2 = form_factor2t1D(t_min,t_max);
    FF2->SetTitle("Form Factor 1D");
    FF2->GetYaxis()->SetTitle("|F(t)|^{2}");
	FF2->GetXaxis()->SetTitle("t [GeV^{2}]");
    FF2->Draw();
    gPad->SetLogy(1);

    TFile *ff = new TFile("FF2t_1D.root","recreate");
    FF2->Write("FF2");
    ff->Close();
}

double calc_FFt2D(double *var, double *par)
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

TF2 *form_factort2D(double tx_min, double tx_max, double ty_min, double ty_max)
{
    TF2 *FF1 = new TF2("Form Factor", "calc_FFt2D", tx_min, tx_max, ty_min, ty_max, 4);
    FF1->SetParameters(197,2.12,6.38,0.70);
    FF1->SetParNames("A","Vo","R","a0");
    return FF1;
}

void draw_FFt_2D()
{
    double tx_min = 0, tx_max = 0.25, ty_min = 0, ty_max = 0.25;
    
    TCanvas *c1 = new TCanvas("c1", "Form Factor", 800, 600);
    TF2 *FF1 = form_factort2D(tx_min,tx_max,ty_min,ty_max);
    FF1->SetTitle("Form Factor 2D");
    FF1->GetYaxis()->SetTitle("t_{y} [GeV^{2}]");
	FF1->GetXaxis()->SetTitle("t_{x} [GeV^{2}]");
    FF1->Draw();

    TFile *ff = new TFile("FFt_2D.root","recreate");
    FF1->Write("FF1");
    ff->Close();
}

void draw_FFt_contour()
{
    double tx_min = 0, tx_max = 0.25, ty_min = 0, ty_max = 0.25;

    TF2 *FF2 = form_factort2D(tx_min,tx_max,ty_min,ty_max);
    // Create a 2D histogram to represent the function
    int bins = 1000; 
    TH2D *hist2D = new TH2D("hist2D", "Form Factor 2D", bins, tx_min, tx_max, bins, ty_min, ty_max);
    hist2D->Sumw2();
    for (int i=1;i<=bins;i++) 
    {
        for (int j=1;j<=bins;j++) 
        {
            double tx = hist2D->GetXaxis()->GetBinCenter(i);
            double ty = hist2D->GetYaxis()->GetBinCenter(j);
            hist2D->SetBinContent(i, j, FF2->Eval(tx, ty));
            cout << "tx: " << tx << endl;
        }
    }
    TFile *FF = new TFile("FFt_2D_contour.root","recreate");
    hist2D->SetTitle("Form Factor 2D ");
    hist2D->GetYaxis()->SetTitle("t_{y} [GeV]");
	hist2D->GetXaxis()->SetTitle("t_{x} [GeV]");
    hist2D->Write();
    FF->Close();

    TCanvas *c1 = new TCanvas("Form Factor", "Form Factor 2D", 800, 600);
    hist2D->Draw();
    gStyle->SetOptStat(0);
}

void draw_yproj_FFt2D()
{
    TFile *f = new TFile("FFt_2D_contour.root");
    TH2D *get = (TH2D *)f->Get("hist2D");
    TH1D *get1 = (TH1D *)get->ProjectionY();
    get1->SetTitle("Form Factor");
    get1->GetXaxis()->SetTitle("t_{y} [GeV^{2}]");
    get1->GetYaxis()->SetTitle("F(t_{x},t_{y})"); 
    get1->Draw();
}

double calc_FF2t2D(double *var, double *par)
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

TF2 *form_factor2t2D(double tx_min, double tx_max, double ty_min, double ty_max)
{
    TF2 *FF2 = new TF2("Form Factor", "calc_FF2t2D", tx_min, tx_max, ty_min, ty_max, 4);
    FF2->SetParameters(197,2.12,6.38,0.70);
    FF2->SetParNames("A","Vo","R","a0");
    return FF2;
}

void draw_FF2t_2D()
{
    double tx_min = 0, tx_max = 0.25, ty_min = 0, ty_max = 0.25;
    
    TCanvas *c1 = new TCanvas("c1", "Form Factor", 800, 600);
    TF2 *FF2 = form_factor2t2D(tx_min,tx_max,ty_min,ty_max);
    FF2->SetTitle("Form Factor 2D");
    FF2->GetYaxis()->SetTitle("t_{y} [GeV^{2}]");
	FF2->GetXaxis()->SetTitle("t_{x} [GeV^{2}]");
    FF2->Draw();

    TFile *ff = new TFile("FF2t_2D.root","recreate");
    FF2->Write("FF2");
    ff->Close();
}

void draw_FF2t_contour()
{
    double tx_min = 0, tx_max = 0.25, ty_min = 0, ty_max = 0.25;

    TF2 *FF2 = form_factor2t2D(tx_min,tx_max,ty_min,ty_max);
    // Create a 2D histogram to represent the function
    int bins = 1000; 
    TH2D *hist2D = new TH2D("hist2D", "Form Factor 2D", bins, tx_min, tx_max, bins, ty_min, ty_max);
    hist2D->Sumw2();
    for (int i=1;i<=bins;i++) 
    {
        for (int j=1;j<=bins;j++) 
        {
            double tx = hist2D->GetXaxis()->GetBinCenter(i);
            double ty = hist2D->GetYaxis()->GetBinCenter(j);
            hist2D->SetBinContent(i, j, FF2->Eval(tx, ty));
            cout << "tx: " << tx << endl;
        }
    }
    TFile *FF = new TFile("FF2t_2D_contour.root","recreate");
    hist2D->SetTitle("Form Factor 2D ");
    hist2D->GetYaxis()->SetTitle("t_{y} [GeV]");
	hist2D->GetXaxis()->SetTitle("t_{x} [GeV]");
    hist2D->Write();
    FF->Close();

    TCanvas *c1 = new TCanvas("Form Factor", "Form Factor 2D", 800, 600);
    hist2D->Draw();
    gStyle->SetOptStat(0);
}

void draw_yproj_FF2t2D()
{
    TFile *f = new TFile("FF2t_2D_contour.root");
    TH2D *get = (TH2D *)f->Get("hist2D");
    TH1D *get1 = (TH1D *)get->ProjectionY();
    get1->SetTitle("Form Factor");
    get1->GetXaxis()->SetTitle("t_{y} [GeV^{2}]");
    get1->GetYaxis()->SetTitle("|F(t_{x},t_{y})|^{2}"); 
    get1->Draw();
    gPad->SetLogy(1);
}

// Integrate out wedge (theta from qy axis)
double calc_FFt_wedge(double *var, double *par)
{
    double A = par[0], Vo = par[1], R = par[2], a0 = par[3], t = par[4];
    double theta = var[0];
    double q = sqrt(t);
    double qx = q*sin(theta), qy = q*cos(theta);
    double qq = sqrt(qx*qx+qy*qy);
	const double arg1 = q*R / hbarc;  
	const double arg2 = hbarc / q;
	const double arg3  = (sin(arg1)-arg1*cos(arg1)) * 4*pi*Vo*arg2*arg2*arg2 / A; 
	return 1/(2*pi)*arg3 / (1. + (a0 * a0 * q*q/(hbarc*hbarc)));  // [-]
}

TF1 *twedge_inte(double t_min, double t_max) 
{
    double theta_min = 0, theta_max = 2*pi;
    TF1 *cutFF = new TF1("Cut Form Factor", [theta_min, theta_max](double *var, double *par) {double t = var[0];  // t is variable
        TF1 f("f", "calc_FFt_wedge", theta_min, theta_max, 5);
        f.SetParameters(197, 2.12, 6.38, 0.70,t);
        return f.Integral(theta_min, theta_max, 1e-12);
    }, t_min, t_max, 0);
    return cutFF;
}

void draw_twedge()
{
    double t_min = 0, t_max = 0.25;
    TCanvas *c1 = new TCanvas("c1", "Form Factor", 800, 600);
    TF1 *FF2 = twedge_inte(t_min, t_max);
    FF2->SetTitle("Form Factor Wedge Cut");
    FF2->GetYaxis()->SetTitle("|F(t)|^{2}");
    FF2->GetXaxis()->SetTitle("t [GeV^{2}]");
    FF2->Draw();

    TFile *ff = new TFile("wedge_t.root", "recreate");
    FF2->Write("FF2");
    ff->Close();
}

void compare_twedge() {
    double t_min = 0, t_max = 0.25;
    TCanvas *c1 = new TCanvas("c1", "Form Factor", 800, 600);
    
    // Original form factor: F(q)
    //TF1 *FF1 = form_factorq1D(q_min,q_max);
    //TFile *f1 = TFile::Open("FFq_1D.root");
    TFile *f1 = new TFile("FFt_1D.root","READ");
    TF1 *FF1 = (TF1*)f1->Get("FF1");
    FF1->SetTitle("Form Factor 1D");
    FF1->GetYaxis()->SetTitle("F(q)");
	FF1->GetXaxis()->SetTitle("q [GeV]");
    FF1->SetLineColor(kRed);
    FF1->DrawCopy();

    // wedge integrated out: theta = 0
    TF1 *FF2 = twedge_inte(t_min, t_max);
    FF2->SetLineColor(kBlue);
    FF2->SetLineStyle(2);
    FF2->Draw("same");

    TLegend *legend = new TLegend(0.7, 0.8, 0.9, 0.9);
    legend->AddEntry(FF1, "no cut", "l");
    legend->AddEntry(FF2, "with wedge cut", "l");
    legend->Draw();

    c1->SaveAs("compare_twedges.png");
}

double calc_FF2t_wedge(double *var, double *par)
{
    double A = par[0], Vo = par[1], R = par[2], a0 = par[3], t = par[4];
    double theta = var[0];
    double q = sqrt(t);
    double qx = q*sin(theta), qy = q*cos(theta);
    double qq = sqrt(qx*qx+qy*qy);
	const double arg1 = q*R / hbarc;  
	const double arg2 = hbarc / q;
	const double arg3  = (sin(arg1)-arg1*cos(arg1)) * 4*pi*Vo*arg2*arg2*arg2 / A; 
	return 1/(2*pi)*arg3 / (1. + (a0 * a0 * q*q/(hbarc*hbarc)))*1/(2*pi)*arg3 / (1. + (a0 * a0 * q*q/(hbarc*hbarc)));  // [-]
}

TF1 *t2wedge_inte(double t_min, double t_max)
{
    double theta_min = 0, theta_max = 2*pi;
    TF1 *cutFF = new TF1("Cut Form Factor", [theta_min, theta_max](double *var, double *par) {double t = var[0];  // t is variable
        TF1 f("f", "calc_FF2t_wedge", theta_min, theta_max, 5);
        f.SetParameters(197, 2.12, 6.38, 0.70,t);
        return f.Integral(theta_min, theta_max, 1e-12);
    }, t_min, t_max, 0);
    return cutFF;
}

void draw_t2wedge()
{
    double t_min = 0, t_max = 0.25;
    TCanvas *c1 = new TCanvas("c1", "Form Factor", 800, 600);
    TF1 *FF2 = t2wedge_inte(t_min, t_max);
    FF2->SetTitle("Form Factor Wedge Cut");
    FF2->GetYaxis()->SetTitle("|F(t)|^{2}");
    FF2->GetXaxis()->SetTitle("t [GeV^{2}]");
    FF2->Draw();
    gPad->SetLogy(1);

    TFile *ff = new TFile("wedge_t2.root", "recreate");
    FF2->Write("FF2");
    ff->Close();
}

/*
// Compare histograms |F(q)|^2 with different theta max angles integrated over
void compare_FF_qwedge()
{
    // Real Form Factor: |Ftrue(q)|^2
    TFile *f1 = new TFile("FF2q_1D.root");
    TH2D *get1_1 = (TH2D *)f1->Get("FF2");
    TH1D *get1 = (TH1D *)get1_1->ProjectionY();
    get1->SetTitle("Wedge Integration of Form Factor: |F(q)|^{2}");
    get1->GetYaxis()->SetTitle("|F(q)|^{2}");
    get1->GetXaxis()->SetTitle("q [GeV]");
    get1->Scale(197./get1->Integral(), "width");
    //get1->GetYaxis()->SetRangeUser(1e-5,1e5);
    get1->SetLineColor(kRed);

    // theta = 2pi
    TFile *f22 = new TFile("testFF_WS_Transformations_real_FF2_2D_wedge_2pi_1000ptq.root");
    TH1D *get22 = (TH1D *)f22->Get("FF");
    get22->Scale(197./get22->Integral(), "width");
    get22->SetLineColor(kBlue-2);
    get22->SetLineStyle(2);

    // theta = pi/3
    TFile *f2 = new TFile("testFF_WS_Transformations_real_FF2_2D_wedge_pio3_1000ptq.root");
    TH1D *get2 = (TH1D *)f2->Get("FF");
    get2->Scale(197./get2->Integral(), "width");
    get2->SetLineColor(kBlue);
    get2->SetLineStyle(2);


    // theta = pi/9
    TFile *f7 = new TFile("testFF_WS_Transformations_real_FF2_2D_wedge_pio9_1000ptq.root");
    TH1D *get7 = (TH1D *)f7->Get("FF");
    get7->Scale(197./get7->Integral(), "width");
    get7->SetLineColor(kViolet);
    get7->SetLineStyle(2);

    // theta = pi/2
    TFile *f3 = new TFile("testFF_WS_Transformations_real_FF2_2D_wedge_pio2_1000ptq.root");
    TH1D *get3 = (TH1D *)f3->Get("FF");
    get3->Scale(197./get3->Integral(), "width");
    get3->SetLineColor(kBlack);
    get3->SetLineStyle(2);

    // theta = pi/6
    TFile *f4 = new TFile("testFF_WS_Transformations_real_FF2_2D_wedge_pio6_1000ptq.root");
    TH1D *get4 = (TH1D *)f4->Get("FF");
    get4->Scale(197./get4->Integral(), "width");
    get4->SetLineColor(kGreen);
    get4->SetLineStyle(2);

    // theta = pi/4
    TFile *f5 = new TFile("testFF_WS_Transformations_real_FF2_2D_wedge_pio4_1000ptq.root");
    TH1D *get5 = (TH1D *)f5->Get("FF");
    get5->Scale(197./get5->Integral(), "width");
    get5->SetLineColor(kCyan);
    get5->SetLineStyle(2);

    // theta = pi
    TFile *f6 = new TFile("testFF_WS_Transformations_real_FF2_2D_wedge_pi_1000ptq.root");
    TH1D *get6 = (TH1D *)f6->Get("FF");
    get6->Scale(197./get6->Integral(), "width");
    get6->SetLineColor(kPink);
    get6->SetLineStyle(2);

    // theta = 0 (1e-12)
    TFile *f8 = new TFile("testFF_WS_Transformations_real_FF2_2D_wedge_0_1000ptq.root");
    TH1D *get8 = (TH1D *)f8->Get("FF");
    get8->Scale(197./get8->Integral(), "width");
    get8->SetLineColor(kOrange);
    get8->SetLineStyle(2);

    TCanvas *c_f1 = new TCanvas();
    get1->Draw();
    get2->Draw("same");
    get22->Draw("same");
    get7->Draw("same");
    get3->Draw("same");
    get4->Draw("same");
    get5->Draw("same");
    get6->Draw("same");
    get8->Draw("same");
    //gPad->SetLogy(1);

    auto legend = new TLegend(0.63,0.21,0.83,0.5);
	legend->AddEntry(get1,"|F (q)|^{2}","l");
    legend->AddEntry(get22,"|F (q)|^{2} #theta = 2 #pi","l");
    legend->AddEntry(get6,"|F (q)|^{2} #theta = #pi","l");
    legend->AddEntry(get3,"|F (q)|^{2} #theta = #pi/2","l");
    legend->AddEntry(get2,"|F (q)|^{2} #theta = #pi/3","l");
    legend->AddEntry(get5,"|F (q)|^{2} #theta = #pi/4","l");
    legend->AddEntry(get4,"|F (q)|^{2} #theta = #pi/6","l");
    legend->AddEntry(get7,"|F (q)|^{2} #theta = #pi/9","l");
    legend->AddEntry(get8,"|F (q)|^{2} #theta = 0","l");
    legend->Draw();
}

// Draw histogram of |F(t)|^2 with wedge integrated over
void wedge_getFF2t()
{
    // Original parameters
    //int bins = 1000.;
    //double step = 0.0001;
    //double min = 0.0001;
    //double max = 0.1;

    // Parameters from saturation model scan
    int bins = 70.;
    double step = 0.00257;
    double min = 0.;
    double max = 0.18;
    TH1D *FF = new TH1D("FF","Form Factor",bins,min,max);
    FF->Sumw2();
    for(int i=0;i<bins;i++)
    {
        double t = (i+1)*step;
        double q = sqrt(t);
        double val = wedge_FF2_inte(q);
        FF->SetBinContent(i+1,val);
        cout << "t: " << t << endl;
    }
    TFile *FF_WS_Transformations = new TFile("70pt_wedge_ff_pio9.root","recreate");
	FF->SetTitle("Form Factor: |F(t)|^{2}");
    FF->GetYaxis()->SetTitle("|F(t)|^{2}");
	FF->GetXaxis()->SetTitle("t [GeV^{2}]");
	FF->Write();
	FF_WS_Transformations->Close();

    TCanvas *c1 = new TCanvas("c1", "True Form Factor", 800, 600);
    FF->Draw();
}

// Compare histograms |F(t)|^2 with different theta max angles integrated over
void draw_FF2_wedget()
{
    // Real Form Factor: |Ftrue|^2
    TFile *f1 = new TFile("FF_WS_Transformations_real_FF2_1D_1000pt.root");
    TH1D *get1 = (TH1D *)f1->Get("FF");
    get1->SetTitle("Wedge Integration of Form Factor: |F(t)|^{2}");
    get1->GetYaxis()->SetTitle("|F(t)|^{2}");
    get1->GetXaxis()->SetTitle("t [GeV^{2}]");
    get1->Scale(197./get1->Integral(), "width");
    get1->GetYaxis()->SetRangeUser(1e-5,1e5);
    get1->SetLineColor(kRed);

    // theta = pi
    TFile *f = new TFile("testFF_WS_Transformations_real_FF2_2D_wedge_pi_1000ptqt.root");
    TH1D *td1 = (TH1D *)f->Get("FF");
    td1->Scale(197./td1->Integral(), "width");
    td1->SetLineColor(kOrange);
    td1->SetLineStyle(2);

    // theta = pi/3
    TFile *f2 = new TFile("testFF_WS_Transformations_real_FF2_2D_wedge_pio3_1000ptqt.root");
    TH1D *get2 = (TH1D *)f2->Get("FF");
    get2->Scale(197./get2->Integral(), "width");
    get2->SetLineColor(kBlue);
    get2->SetLineStyle(2);

    // theta = pi/9
    TFile *f7 = new TFile("testFF_WS_Transformations_real_FF2_2D_wedge_pio9_1000ptqt.root");
    TH1D *get7 = (TH1D *)f7->Get("FF");
    get7->Scale(197./get7->Integral(), "width");
    get7->SetLineColor(kViolet);
    get7->SetLineStyle(2);

    // theta = pi/2
    TFile *f3 = new TFile("testFF_WS_Transformations_real_FF2_2D_wedge_pio2_1000ptqt.root");
    TH1D *get3 = (TH1D *)f3->Get("FF");
    get3->Scale(197./get3->Integral(), "width");
    get3->SetLineColor(kBlack);
    get3->SetLineStyle(2);

    // theta = pi/6
    TFile *f4 = new TFile("testFF_WS_Transformations_real_FF2_2D_wedge_pio6_1000ptqt.root");
    TH1D *get4 = (TH1D *)f4->Get("FF");
    get4->Scale(197./get4->Integral(), "width");
    get4->SetLineColor(kGreen);
    get4->SetLineStyle(2);

    // theta = pi/4
    TFile *f5 = new TFile("testFF_WS_Transformations_real_FF2_2D_wedge_pio4_1000ptqt.root");
    TH1D *get5 = (TH1D *)f5->Get("FF");
    get5->Scale(197./get5->Integral(), "width");
    get5->SetLineColor(kCyan);
    get5->SetLineStyle(2);

    // theta = 2pi
    TFile *f6 = new TFile("testFF_WS_Transformations_real_FF2_2D_wedge_2pi_1000ptqt.root");
    TH1D *get6 = (TH1D *)f6->Get("FF");
    get6->Scale(197./get6->Integral(), "width");
    get6->SetLineColor(kPink);
    get6->SetLineStyle(2);

    // theta = 0 (1e-12)
    TFile *f8 = new TFile("testFF_WS_Transformations_real_FF2_2D_wedge_0_1000ptqt.root");
    TH1D *get8 = (TH1D *)f8->Get("FF");
    get8->Scale(197./get8->Integral(), "width");
    get8->SetLineColor(kOrange);
    get8->SetLineStyle(2);

    TCanvas *c_f1 = new TCanvas();
    get1->Draw();
    get2->Draw("same");
    get7->Draw("same");
    get3->Draw("same");
    get4->Draw("same");
    get5->Draw("same");
    td1->Draw("same");
    get6->Draw("same");
    get8->Draw("same");
    gPad->SetLogy(1);

    auto legend = new TLegend(0.63,0.21,0.83,0.5);
	legend->AddEntry(get1,"|F (t)|^{2}","l");
    legend->AddEntry(get6,"|F (t)|^{2} #theta = 2 #pi","l");
    legend->AddEntry(td1,"|F (t)|^{2} #theta = #pi","l");
    legend->AddEntry(get3,"|F (t)|^{2} #theta = #pi/2","l");
    legend->AddEntry(get2,"|F (t)|^{2} #theta = #pi/3","l");
    legend->AddEntry(get5,"|F (t)|^{2} #theta = #pi/4","l");
    legend->AddEntry(get4,"|F (t)|^{2} #theta = #pi/6","l");
    legend->AddEntry(get7,"|F (t)|^{2} #theta = #pi/9","l");
    legend->AddEntry(get8,"|F (t)|^{2} #theta = 0","l");
    legend->Draw();
}

// Draw histogram to compare |Ftrue(t)|^2 with |F(t)|^2 with wedge cut, theta = 2pi
void compare()
{
    // Real Form Factor: |Ftrue|^2
    TFile *f1 = new TFile("FF_WS_Transformations_real_FF2_1D_1000pt.root");
    TH1D *get1 = (TH1D *)f1->Get("FF");
    get1->SetTitle("Wedge Integration of Form Factor: |F(t)|^{2}");
    get1->GetYaxis()->SetTitle("|F(t)|^{2}");
    get1->GetXaxis()->SetTitle("t [GeV^{2}]");
    get1->Scale(197./get1->Integral(), "width");
    get1->GetYaxis()->SetRangeUser(1e-5,1e5);
    get1->SetLineColor(kRed);

    // |F(t)|^2, theta = 2pi
    TFile *f5 = new TFile("testFF_WS_Transformations_real_FF2_2D_wedge_2pi_1000ptqt.root");
    TH1D *get5 = (TH1D *)f5->Get("FF");
    get5->Scale(197./get5->Integral(), "width");
    get5->SetLineColor(kCyan);
    get5->SetLineStyle(2);
 
    TCanvas *c_f1 = new TCanvas();
    get1->Draw();
    get5->Draw("same");
    gPad->SetLogy(1);

    auto legend = new TLegend(0.63,0.21,0.83,0.5);
	legend->AddEntry(get1,"|F (t)|^{2}","l");
    legend->AddEntry(get5,"|F (t)|^{2} (wedge) ","l");
    legend->Draw();
}

*/


/********************************************** 
 
 * Add Resolutiuon to form factor
 
 * Parameters:

 * Variables:
    
**********************************************/
double calc_smearing(double *var, double *par)
{
    double sigma = par[0], qx_prime = par[1];;
    double qx = var[0];
    return exp(-(qx-qx_prime)*(qx-qx_prime)/(2*sigma*sigma));
}

TF1 *smear(double qx_min, double qx_max)
{
    double sigma = 0.2;
    TF1 *smearFF = new TF1("Smearing Function", "calc_smearing", qx_min, qx_max, 2);
    smearFF->SetParameters(sigma, 0);
    return smearFF;
}

TF1 *add_res(double qxPrime_min, double qxPrime_max)
{
    double q_min = 0, q_max = 0.5;
    TF1 *wedgeFF = qwedge_inte(q_min, q_max);
    TF1 *smearFF = smear(qx_min, qx_max);
    TF1Convolution *conv = new TF1Convolution(wedge, smearFF, qxPrime_min, qxPrime_max, true);
    TF1 *resFF = new TF1("resFF", "conv", qxPrime_min, qxPrime_max, conv->GetNpar());
    return resFF;
}

void draw_FFres()
{
    double qxPrime_min = 0, qxPrime_max = 0.5;
    TCanvas *c1 = new TCanvas("c1", "Form Factor", 800, 600);
    TF1 *FF2 = add_res(qxPrime_min, qxPrime_max);
    FF2->SetTitle("Form Factor with Resolution");
    FF2->GetYaxis()->SetTitle("F(q)");
    FF2->GetXaxis()->SetTitle("q [GeV]");
    FF2->Draw();

    TFile *ff = new TFile("res_q.root", "recreate");
    FF2->Write("FF2");
    ff->Close();
}


/*
// Draw histogram of |F(q)|^2 with wedge cut and resolution
void testwedge_getFF2()
{
    int bins = 1000.;
    double step = 0.0001;
    double min = 0.0001;
    double max = 0.1;
    TH1D *FF = new TH1D("FF","Form Factor",bins,min,max);
    FF->Sumw2();
    for(int i=0;i<bins;i++)
    {
        double q = (i+1)*step;
        double val = testwedge_FF2_inte(q);
        FF->SetBinContent(i+1,val);
        cout << "q: " << q << endl;
    }
    TFile *FF_WS_Transformations = new TFile("ttestFF_WS_Transformations_real_FF2_2D_wedge_0_1000ptq_10.root","recreate");
	FF->SetTitle("Form Factor: |F(q)|^{2}");
    FF->GetYaxis()->SetTitle("|F(q)|^{2}");
	FF->GetXaxis()->SetTitle("q [GeV]");
	FF->Write();
	FF_WS_Transformations->Close();

    TCanvas *c1 = new TCanvas("c1", "True Form Factor", 800, 600);
    FF->Draw();
}

// Draw histogram of |F(t)|^2 with wedge cut and resolution
void testwedge_getFF2t()
{
    // Original parameters
    //int bins = 1000.;
    //double step = 0.0001;
    //double min = 0.0001;
    //double max = 0.1;

    // Parameters from the saturation model histogram scan
    int bins = 70.;
    double step = 0.00257;
    double min = 0.;
    double max = 0.18;

    TH1D *FF = new TH1D("FF","Form Factor",bins,min,max);
    FF->Sumw2();
    for(int i=0;i<bins;i++)
    {
        double t = (i+1)*step;
        double q = sqrt(t);
        double val = testwedge_FF2_inte(q);
        FF->SetBinContent(i+1,val);
        cout << "t: " << t << endl;
    }
    TFile *FF_WS_Transformations = new TFile("70pt_res_wedge_ff_pio9_200.root","recreate");
	FF->SetTitle("Form Factor: |F(t)|^{2}");
    FF->GetYaxis()->SetTitle("|F(t)|^{2}");
	FF->GetXaxis()->SetTitle("t [GeV^{2}]");
	FF->Write();
	FF_WS_Transformations->Close();

    TCanvas *c1 = new TCanvas("c1", "True Form Factor", 800, 600);
    FF->Draw();
}




                //======================= Histograms of each angle with different resolutions ===========================//

// |F(t)|^2, theta = pi/2
void pio2_draw_FF2_wedget()
{
    // theta = pi/2 no resolution
    TFile *f1 = new TFile("wedge_ff_pio2.root");
    TH1D *get1 = (TH1D *)f1->Get("FF");
    get1->Scale(197./get1->Integral(), "width");
    get1->SetTitle("Wedge Integration of |F(t)|^{2} with Resolution (#theta_{max} = #pi/2)");
    get1->GetYaxis()->SetTitle("|F(t)|^{2}");
    get1->GetXaxis()->SetTitle("t [GeV^{2}]");
    get1->Scale(197./get1->Integral(), "width");
    get1->GetYaxis()->SetRangeUser(1e-5,1e5);
    get1->SetLineColor(kBlack);

    // theta = pi/2 with 100 MeV resolution
    TFile *f3 = new TFile("res_wedge_ff_pio2_100.root");
    TH1D *get3 = (TH1D *)f3->Get("FF");
    get3->Scale(197./get3->Integral(), "width");
    get3->SetLineColor(kBlue);
    //get3->SetLineStyle(2);

    // theta = pi/2 with 50 MeV resolution
    TFile *f4 = new TFile("res_wedge_ff_pio2_50.root");
    TH1D *get4 = (TH1D *)f4->Get("FF");
    get4->Scale(197./get4->Integral(), "width");
    get4->SetLineColor(kGreen);
    //get4->SetLineStyle(2);

    // theta = pi/2 with 10 MeV resolution
    TFile *f5 = new TFile("res_wedge_ff_pio2_10.root");
    TH1D *get5 = (TH1D *)f5->Get("FF");
    get5->Scale(197./get5->Integral(), "width");
    get5->SetLineColor(kCyan);
    //get5->SetLineStyle(2);

    // theta = pi/2 with 150 MeV resolution
    TFile *f6 = new TFile("res_wedge_ff_pio2_150.root");
    TH1D *get6 = (TH1D *)f6->Get("FF");
    get6->Scale(197./get6->Integral(), "width");
    get6->SetLineColor(kViolet);

    // theta = pi/2 with 5 MeV resolution
    TFile *f7 = new TFile("res_wedge_ff_pio2_5.root");
    TH1D *get7 = (TH1D *)f7->Get("FF");
    get7->Scale(197./get7->Integral(), "width");
    get7->SetLineColor(kOrange);

    // theta = pi/2 with 200 MeV resolution
    TFile *f8 = new TFile("res_wedge_ff_pio2_200.root");
    TH1D *get8 = (TH1D *)f8->Get("FF");
    get8->Scale(197./get8->Integral(), "width");
    get8->SetLineColor(kRed);

    TCanvas *c_f1 = new TCanvas();
    get1->Draw();
    get3->Draw("same");
    get4->Draw("same");
    get5->Draw("same");
    get6->Draw("same");
    get7->Draw("same");
    get8->Draw("same");
    gPad->SetLogy(1);
    gStyle->SetOptStat(0);
    

    auto legend = new TLegend(0.63,0.21,0.83,0.5);
	legend->AddEntry(get1,"|F (t)|^{2} no resolution","l");
    legend->AddEntry(get7,"|F (t)|^{2} with 5 MeV resolution","l");
    legend->AddEntry(get5,"|F (t)|^{2} with 10 MeV resolution","l");
    legend->AddEntry(get4,"|F (t)|^{2} with 50 MeV resolution","l");
    legend->AddEntry(get3,"|F (t)|^{2} with 100 MeV resolution","l");
    legend->AddEntry(get6,"|F (t)|^{2} with 150 MeV resolution","l");
    legend->AddEntry(get8,"|F (t)|^{2} with 200 MeV resolution","l");
    legend->Draw();
}

// |F(t)|^2, theta = pi/4
void pio4_draw_FF2_wedget()
{
    // theta = pi/4 no resolution
    TFile *f1 = new TFile("wedge_ff_pio4.root");
    TH1D *get1 = (TH1D *)f1->Get("FF");
    get1->Scale(197./get1->Integral(), "width");
    get1->SetTitle("Wedge Integration of |F(t)|^{2} with Resolution (#theta_{max} = #pi/4)");
    get1->GetYaxis()->SetTitle("|F(t)|^{2}");
    get1->GetXaxis()->SetTitle("t [GeV^{2}]");
    get1->Scale(197./get1->Integral(), "width");
    get1->GetYaxis()->SetRangeUser(1e-5,1e5);
    get1->SetLineColor(kBlack);

    // theta = pi/4 with 100 MeV resolution
    TFile *f3 = new TFile("res_wedge_ff_pio4_100.root");
    TH1D *get3 = (TH1D *)f3->Get("FF");
    get3->Scale(197./get3->Integral(), "width");
    get3->SetLineColor(kBlue);
    //get3->SetLineStyle(2);

    // theta = pi/4 with 50 MeV resolution
    TFile *f4 = new TFile("res_wedge_ff_pio4_50.root");
    TH1D *get4 = (TH1D *)f4->Get("FF");
    get4->Scale(197./get4->Integral(), "width");
    get4->SetLineColor(kGreen);
    //get4->SetLineStyle(2);

    // theta = pi/4 with 10 MeV resolution
    TFile *f5 = new TFile("res_wedge_ff_pio4_10.root");
    TH1D *get5 = (TH1D *)f5->Get("FF");
    get5->Scale(197./get5->Integral(), "width");
    get5->SetLineColor(kCyan);
    //get5->SetLineStyle(2);

    // theta = pi/4 with 150 MeV resolution
    TFile *f6 = new TFile("res_wedge_ff_pio4_150.root");
    TH1D *get6 = (TH1D *)f6->Get("FF");
    get6->Scale(197./get6->Integral(), "width");
    get6->SetLineColor(kViolet);

    // theta = pi/4 with 5 MeV resolution
    TFile *f7 = new TFile("res_wedge_ff_pio4_5.root");
    TH1D *get7 = (TH1D *)f7->Get("FF");
    get7->Scale(197./get7->Integral(), "width");
    get7->SetLineColor(kOrange);
    //get7->SetLineStyle(2);

    // theta = pi/4 with 200 MeV resolution
    TFile *f8 = new TFile("res_wedge_ff_pio4_200.root");
    TH1D *get8 = (TH1D *)f8->Get("FF");
    get8->Scale(197./get8->Integral(), "width");
    get8->SetLineColor(kRed);

    TCanvas *c_f1 = new TCanvas();
    get1->Draw();
    get3->Draw("same");
    get4->Draw("same");
    get5->Draw("same");
    get6->Draw("same");
    get7->Draw("same");
    get8->Draw("same");
    gPad->SetLogy(1);
    gStyle->SetOptStat(0);

    auto legend = new TLegend(0.63,0.21,0.83,0.5);
	legend->AddEntry(get1,"|F (t)|^{2} no resolution","l");
    legend->AddEntry(get7,"|F (t)|^{2} with 5 MeV resolution","l");
    legend->AddEntry(get5,"|F (t)|^{2} with 10 MeV resolution","l");
    legend->AddEntry(get4,"|F (t)|^{2} with 50 MeV resolution","l");
    legend->AddEntry(get3,"|F (t)|^{2} with 100 MeV resolution","l");
    legend->AddEntry(get6,"|F (t)|^{2} with 150 MeV resolution","l");
    legend->AddEntry(get8,"|F (t)|^{2} with 200 MeV resolution","l");
    legend->Draw();
}

// |F(t)|^2, theta = pi/3
void pio3_draw_FF2_wedget()
{
    // theta = pi/3 no resolution
    TFile *f1 = new TFile("wedge_ff_pio3.root");
    TH1D *get1 = (TH1D *)f1->Get("FF");
    get1->Scale(197./get1->Integral(), "width");
    get1->SetTitle("Wedge Integration of |F(t)|^{2} with Resolution (#theta_{max} = #pi/3)");
    get1->GetYaxis()->SetTitle("|F(t)|^{2}");
    get1->GetXaxis()->SetTitle("t [GeV^{2}]");
    get1->Scale(197./get1->Integral(), "width");
    get1->GetYaxis()->SetRangeUser(1e-5,1e5);
    get1->SetLineColor(kBlack);

    // theta = pi/3 with 100 MeV resolution
    TFile *f3 = new TFile("res_wedge_ff_pio3_100.root");
    TH1D *get3 = (TH1D *)f3->Get("FF");
    get3->Scale(197./get3->Integral(), "width");
    get3->SetLineColor(kBlue);
    //get3->SetLineStyle(2);

    // theta = pi/3 with 50 MeV resolution
    TFile *f4 = new TFile("res_wedge_ff_pio3_50.root");
    TH1D *get4 = (TH1D *)f4->Get("FF");
    get4->Scale(197./get4->Integral(), "width");
    get4->SetLineColor(kGreen);
    //get4->SetLineStyle(2);

    // theta = pi/3 with 10 MeV resolution
    TFile *f5 = new TFile("res_wedge_ff_pio3_10.root");
    TH1D *get5 = (TH1D *)f5->Get("FF");
    get5->Scale(197./get5->Integral(), "width");
    get5->SetLineColor(kCyan);
    //get5->SetLineStyle(2);

    // theta = pi/3 with 150 MeV resolution
    TFile *f6 = new TFile("res_wedge_ff_pio3_150.root");
    TH1D *get6 = (TH1D *)f6->Get("FF");
    get6->Scale(197./get6->Integral(), "width");
    get6->SetLineColor(kViolet);

    // theta = pi/3 with 5 MeV resolution
    TFile *f7 = new TFile("res_wedge_ff_pio3_5.root");
    TH1D *get7 = (TH1D *)f7->Get("FF");
    get7->Scale(197./get7->Integral(), "width");
    get7->SetLineColor(kOrange);

    // theta = pi/3 with 200 MeV resolution
    TFile *f8 = new TFile("res_wedge_ff_pio3_200.root");
    TH1D *get8 = (TH1D *)f8->Get("FF");
    get8->Scale(197./get8->Integral(), "width");
    get8->SetLineColor(kRed);

    TCanvas *c_f1 = new TCanvas();
    get1->Draw();
    get3->Draw("same");
    get4->Draw("same");
    get5->Draw("same");
    get6->Draw("same");
    get7->Draw("same");
    get8->Draw("same");
    gPad->SetLogy(1);
    gStyle->SetOptStat(0);

    auto legend = new TLegend(0.63,0.21,0.83,0.5);
	legend->AddEntry(get1,"|F (t)|^{2} no resolution","l");
    legend->AddEntry(get7,"|F (t)|^{2} with 5 MeV resolution","l");
    legend->AddEntry(get5,"|F (t)|^{2} with 10 MeV resolution","l");
    legend->AddEntry(get4,"|F (t)|^{2} with 50 MeV resolution","l");
    legend->AddEntry(get3,"|F (t)|^{2} with 100 MeV resolution","l");
    legend->AddEntry(get6,"|F (t)|^{2} with 150 MeV resolution","l");
    legend->AddEntry(get8,"|F (t)|^{2} with 200 MeV resolution","l");
    legend->Draw();
}

// |F(t)|^2, theta = pi/6
void pio6_draw_FF2_wedget()
{
    // theta = pi/6 no resolution
    TFile *f1 = new TFile("wedge_ff_pio6.root");
    TH1D *get1 = (TH1D *)f1->Get("FF");
    get1->Scale(197./get1->Integral(), "width");
    get1->SetTitle("Wedge Integration of |F(t)|^{2} with Resolution (#theta_{max} = #pi/6)");
    get1->GetYaxis()->SetTitle("|F(t)|^{2}");
    get1->GetXaxis()->SetTitle("t [GeV^{2}]");
    get1->Scale(197./get1->Integral(), "width");
    get1->GetYaxis()->SetRangeUser(1e-5,1e5);
    get1->SetLineColor(kBlack);

    // theta = pi/6 with 100 MeV resolution
    TFile *f3 = new TFile("res_wedge_ff_pio6_100.root");
    TH1D *get3 = (TH1D *)f3->Get("FF");
    get3->Scale(197./get3->Integral(), "width");
    get3->SetLineColor(kBlue);
    //get3->SetLineStyle(2);

    // theta = pi/6 with 50 MeV resolution
    TFile *f4 = new TFile("res_wedge_ff_pio6_50.root");
    TH1D *get4 = (TH1D *)f4->Get("FF");
    get4->Scale(197./get4->Integral(), "width");
    get4->SetLineColor(kGreen);
    //get4->SetLineStyle(2);

    // theta = pi/6 with 10 MeV resolution
    TFile *f5 = new TFile("res_wedge_ff_pio6_10.root");
    TH1D *get5 = (TH1D *)f5->Get("FF");
    get5->Scale(197./get5->Integral(), "width");
    get5->SetLineColor(kCyan);
    //get5->SetLineStyle(2);

    // theta = pi/6 with 150 MeV resolution
    TFile *f6 = new TFile("res_wedge_ff_pio6_150.root");
    TH1D *get6 = (TH1D *)f6->Get("FF");
    get6->Scale(197./get6->Integral(), "width");
    get6->SetLineColor(kViolet);

    // theta = pi/6 with 5 MeV resolution
    TFile *f7 = new TFile("res_wedge_ff_pio6_5.root");
    TH1D *get7 = (TH1D *)f7->Get("FF");
    get7->Scale(197./get7->Integral(), "width");
    get7->SetLineColor(kOrange);
    //get7->SetLineStyle(2);

    // theta = pi/6 with 200 MeV resolution
    TFile *f8 = new TFile("res_wedge_ff_pio6_200.root");
    TH1D *get8 = (TH1D *)f8->Get("FF");
    get8->Scale(197./get8->Integral(), "width");
    get8->SetLineColor(kRed);
    //get8->SetLineStyle(2);

    TCanvas *c_f1 = new TCanvas();
    get1->Draw();
    get3->Draw("same");
    get4->Draw("same");
    get5->Draw("same");
    get6->Draw("same");
    get7->Draw("same");
    get8->Draw("same");
    gPad->SetLogy(1);
    gStyle->SetOptStat(0);

    auto legend = new TLegend(0.63,0.21,0.83,0.5);
	legend->AddEntry(get1,"|F (t)|^{2} no resolution","l");
    legend->AddEntry(get7,"|F (t)|^{2} with 5 MeV resolution","l");
    legend->AddEntry(get5,"|F (t)|^{2} with 10 MeV resolution","l");
    legend->AddEntry(get4,"|F (t)|^{2} with 50 MeV resolution","l");
    legend->AddEntry(get3,"|F (t)|^{2} with 100 MeV resolution","l");
    legend->AddEntry(get6,"|F (t)|^{2} with 150 MeV resolution","l");
    legend->AddEntry(get8,"|F (t)|^{2} with 200 MeV resolution","l");
    legend->Draw();
}

// |F(t)|^2, theta = pi/9
void pio9_draw_FF2_wedget()
{
    // theta = pi/9 no resolution
    TFile *f1 = new TFile("wedge_ff_pio9.root");
    TH1D *get1 = (TH1D *)f1->Get("FF");
    get1->Scale(197./get1->Integral(), "width");
    get1->SetTitle("Wedge Integration of |F(t)|^{2} with Resolution (#theta_{max} = #pi/9)");
    get1->GetYaxis()->SetTitle("|F(t)|^{2}");
    get1->GetXaxis()->SetTitle("t [GeV^{2}]");
    get1->Scale(197./get1->Integral(), "width");
    get1->GetYaxis()->SetRangeUser(1e-5,1e5);
    get1->SetLineColor(kBlack);

    // theta = pi/9 with 100 MeV resolution
    TFile *f3 = new TFile("res_wedge_ff_pio9_100.root");
    TH1D *get3 = (TH1D *)f3->Get("FF");
    get3->Scale(197./get3->Integral(), "width");
    get3->SetLineColor(kBlue);
    //get3->SetLineStyle(2);

    // theta = pi/9 with 50 MeV resolution
    TFile *f4 = new TFile("res_wedge_ff_pio9_50.root");
    TH1D *get4 = (TH1D *)f4->Get("FF");
    get4->Scale(197./get4->Integral(), "width");
    get4->SetLineColor(kGreen);
    //get4->SetLineStyle(2);

    // theta = pi/9 with 10 MeV resolution
    TFile *f5 = new TFile("res_wedge_ff_pio9_10.root");
    TH1D *get5 = (TH1D *)f5->Get("FF");
    get5->Scale(197./get5->Integral(), "width");
    get5->SetLineColor(kCyan);
    //get5->SetLineStyle(2);

    // theta = pi/9 with 150 MeV resolution
    TFile *f6 = new TFile("res_wedge_ff_pio9_150.root");
    TH1D *get6 = (TH1D *)f6->Get("FF");
    get6->Scale(197./get6->Integral(), "width");
    get6->SetLineColor(kViolet);

    // theta = pi/9 with 5 MeV resolution
    TFile *f7 = new TFile("res_wedge_ff_pio9_5.root");
    TH1D *get7 = (TH1D *)f7->Get("FF");
    get7->Scale(197./get7->Integral(), "width");
    get7->SetLineColor(kOrange);
    //get7->SetLineStyle(2);

    // theta = pi/9 with 200 MeV resolution
    TFile *f8 = new TFile("res_wedge_ff_pio9_200.root");
    TH1D *get8 = (TH1D *)f8->Get("FF");
    get8->Scale(197./get8->Integral(), "width");
    get8->SetLineColor(kRed);
    //get8->SetLineStyle(2);

    TCanvas *c_f1 = new TCanvas();
    get1->Draw();
    get3->Draw("same");
    get4->Draw("same");
    get5->Draw("same");
    get6->Draw("same");
    get7->Draw("same");
    get8->Draw("same");
    gPad->SetLogy(1);
    gStyle->SetOptStat(0);

    auto legend = new TLegend(0.63,0.21,0.83,0.5);
	legend->AddEntry(get1,"|F (t)|^{2} no resolution","l");
    legend->AddEntry(get7,"|F (t)|^{2} with 5 MeV resolution","l");
    legend->AddEntry(get5,"|F (t)|^{2} with 10 MeV resolution","l");
    legend->AddEntry(get4,"|F (t)|^{2} with 50 MeV resolution","l");
    legend->AddEntry(get3,"|F (t)|^{2} with 100 MeV resolution","l");
    legend->AddEntry(get6,"|F (t)|^{2} with 150 MeV resolution","l");
    legend->AddEntry(get8,"|F (t)|^{2} with 200 MeV resolution","l");
    legend->Draw();
}

// |F(t)|^2, theta = 0 >> (1e-12)
void zero_draw_FF2_wedget()
{
    // theta = 0 (1e-12) no resolution
    TFile *f1 = new TFile("wedge_ff_0.root");
    TH1D *get1 = (TH1D *)f1->Get("FF");
    get1->Scale(197./get1->Integral(), "width");
    get1->SetTitle("Wedge Integration of |F(t)|^{2} with Resolution (#theta_{max} = 0)");
    get1->GetYaxis()->SetTitle("|F(t)|^{2}");
    get1->GetXaxis()->SetTitle("t [GeV^{2}]");
    get1->Scale(197./get1->Integral(), "width");
    get1->GetYaxis()->SetRangeUser(1e-5,1e5);
    get1->SetLineColor(kBlack);

    // theta = 0 with 100 MeV resolution
    TFile *f3 = new TFile("res_wedge_ff_0_100.root");
    TH1D *get3 = (TH1D *)f3->Get("FF");
    get3->SetTitle("Wedge Integration of |F(t)|^{2} with Resolution (#theta_{max} = 0)");
    get3->Scale(197./get3->Integral(), "width");
    get3->SetLineColor(kBlue);
    //get3->SetLineStyle(2);

    // theta = 0 with 50 MeV resolution
    TFile *f4 = new TFile("res_wedge_ff_0_50.root");
    TH1D *get4 = (TH1D *)f4->Get("FF");
    get4->Scale(197./get4->Integral(), "width");
    get4->SetLineColor(kGreen);
    //get4->SetLineStyle(2);

    // theta = 0 with 10 MeV resolution
    TFile *f5 = new TFile("res_wedge_ff_0_10.root");
    TH1D *get5 = (TH1D *)f5->Get("FF");
    get5->Scale(197./get5->Integral(), "width");
    get5->SetLineColor(kCyan);
    //get5->SetLineStyle(2);

    // theta = 0 with 150 MeV resolution
    TFile *f6 = new TFile("res_wedge_ff_0_150.root");
    TH1D *get6 = (TH1D *)f6->Get("FF");
    get6->Scale(197./get6->Integral(), "width");
    get6->SetLineColor(kViolet);

    // theta = 0 with 5 MeV resolution
    TFile *f7 = new TFile("res_wedge_ff_0_5.root");
    TH1D *get7 = (TH1D *)f7->Get("FF");
    get7->Scale(197./get7->Integral(), "width");
    get7->SetLineColor(kOrange);
    //get7->SetLineStyle(2);

    // theta = 0 with 200 MeV resolution
    TFile *f8 = new TFile("res_wedge_ff_0_200.root");
    TH1D *get8 = (TH1D *)f8->Get("FF");
    get8->Scale(197./get8->Integral(), "width");
    get8->SetLineColor(kRed);
    //get8->SetLineStyle(2);

    TCanvas *c_f1 = new TCanvas();
    get1->Draw();
    get3->Draw();
    get4->Draw("same");
    get5->Draw("same");
    get6->Draw("same");
    get7->Draw("same");
    get8->Draw("same");
    gPad->SetLogy(1);
    gStyle->SetOptStat(0);

    auto legend = new TLegend(0.63,0.21,0.83,0.5);
	legend->AddEntry(get1,"|F (t)|^{2} no resolution","l");
    legend->AddEntry(get7,"|F (t)|^{2} with 5 MeV resolution","l");
    legend->AddEntry(get5,"|F (t)|^{2} with 10 MeV resolution","l");
    legend->AddEntry(get4,"|F (t)|^{2} with 50 MeV resolution","l");
    legend->AddEntry(get3,"|F (t)|^{2} with 100 MeV resolution","l");
    legend->AddEntry(get6,"|F (t)|^{2} with 150 MeV resolution","l");
    legend->AddEntry(get8,"|F (t)|^{2} with 200 MeV resolution","l");
    legend->Draw();
}

// |F(t)|^2, theta = 2pi
void twopi_draw_FF2_wedget()
{
    // theta = 2pi no resolution
    TFile *f1 = new TFile("wedge_ff_2pi.root");
    TH1D *get1 = (TH1D *)f1->Get("FF");
    get1->Scale(197./get1->Integral(), "width");
    get1->SetTitle("Wedge Integration of |F(t)|^{2} with Resolution (#theta_{max} = 2 #pi)");
    get1->GetYaxis()->SetTitle("|F(t)|^{2}");
    get1->GetXaxis()->SetTitle("t [GeV^{2}]");
    get1->Scale(197./get1->Integral(), "width");
    get1->GetYaxis()->SetRangeUser(1e-5,1e5);
    get1->SetLineColor(kBlack);

    // theta = 2pi with 100 MeV resolution
    TFile *f3 = new TFile("res_wedge_ff_2pi_100.root");
    TH1D *get3 = (TH1D *)f3->Get("FF");
    get3->Scale(197./get3->Integral(), "width");
    get3->SetLineColor(kBlue);
    //get3->SetLineStyle(2);

    // theta = 2pi with 10 MeV resolution
    TFile *f4 = new TFile("res_wedge_ff_2pi_10.root");
    TH1D *get4 = (TH1D *)f4->Get("FF");
    get4->Scale(197./get4->Integral(), "width");
    get4->SetLineColor(kGreen);
    //get4->SetLineStyle(2);

    // theta = 2pi with 50 MeV resolution
    TFile *f5 = new TFile("res_wedge_ff_2pi_50.root");
    TH1D *get5 = (TH1D *)f5->Get("FF");
    get5->Scale(197./get5->Integral(), "width");
    get5->SetLineColor(kCyan);
    //get5->SetLineStyle(2);

    // theta = 2pi with 150 MeV resolution
    TFile *f6 = new TFile("res_wedge_ff_2pi_150.root");
    TH1D *get6 = (TH1D *)f6->Get("FF");
    get6->Scale(197./get6->Integral(), "width");
    get6->SetLineColor(kViolet);

    // theta = 2pi with 5 MeV resolution
    TFile *f7 = new TFile("res_wedge_ff_2pi_5.root");
    TH1D *get7 = (TH1D *)f7->Get("FF");
    get7->Scale(197./get7->Integral(), "width");
    get7->SetLineColor(kOrange);
    //get7->SetLineStyle(2);

    // theta = 2pi with 200 MeV resolution
    TFile *f8 = new TFile("res_wedge_ff_2pi_200.root");
    TH1D *get8 = (TH1D *)f8->Get("FF");
    get8->Scale(197./get8->Integral(), "width");
    get8->SetLineColor(kRed);
    //get8->SetLineStyle(2);

    TCanvas *c_f1 = new TCanvas();
    get1->Draw();
    get3->Draw("same");
    get4->Draw("same");
    get5->Draw("same");
    get6->Draw("same");
    get7->Draw("same");
    get8->Draw("same");
    gPad->SetLogy(1);
    gStyle->SetOptStat(0);

    auto legend = new TLegend(0.63,0.21,0.83,0.5);
	legend->AddEntry(get1,"|F (t)|^{2} no resolution","l");
    legend->AddEntry(get7,"|F (t)|^{2} with 5 MeV resolution","l");
    legend->AddEntry(get4,"|F (t)|^{2} with 10 MeV resolution","l");
    legend->AddEntry(get5,"|F (t)|^{2} with 50 MeV resolution","l");
    legend->AddEntry(get3,"|F (t)|^{2} with 100 MeV resolution","l");
    legend->AddEntry(get6,"|F (t)|^{2} with 150 MeV resolution","l");
    legend->AddEntry(get8,"|F (t)|^{2} with 200 MeV resolution","l");
    legend->Draw();
}

// |F(t)|^2, theta = pi
void pi_draw_FF2_wedget()
{
    // theta = pi no resolution
    TFile *f1 = new TFile("wedge_ff_pi.root");
    TH1D *get1 = (TH1D *)f1->Get("FF");
    get1->Scale(197./get1->Integral(), "width");
    get1->SetTitle("Wedge Integration of |F(t)|^{2} with Resolution (#theta_{max} = #pi)");
    get1->GetYaxis()->SetTitle("|F(t)|^{2}");
    get1->GetXaxis()->SetTitle("t [GeV^{2}]");
    get1->Scale(197./get1->Integral(), "width");
    get1->GetYaxis()->SetRangeUser(1e-5,1e5);
    get1->SetLineColor(kBlack);

    // theta = pi with 100 MeV resolution
    TFile *f3 = new TFile("res_wedge_ff_pi_100.root");
    TH1D *get3 = (TH1D *)f3->Get("FF");
    get3->Scale(197./get3->Integral(), "width");
    get3->SetLineColor(kBlue);
    //get3->SetLineStyle(2);

    // theta = pi with 10 MeV resolution
    TFile *f4 = new TFile("res_wedge_ff_pi_10.root");
    TH1D *get4 = (TH1D *)f4->Get("FF");
    get4->Scale(197./get4->Integral(), "width");
    get4->SetLineColor(kGreen);
    //get4->SetLineStyle(2);

    // theta = pi with 50 MeV resolution
    TFile *f5 = new TFile("res_wedge_ff_pi_50.root");
    TH1D *get5 = (TH1D *)f5->Get("FF");
    get5->Scale(197./get5->Integral(), "width");
    get5->SetLineColor(kCyan);
    //get5->SetLineStyle(2);

    // theta = pi with 150 MeV resolution
    TFile *f6 = new TFile("res_wedge_ff_pi_150.root");
    TH1D *get6 = (TH1D *)f6->Get("FF");
    get6->Scale(197./get6->Integral(), "width");
    get6->SetLineColor(kViolet);

    // theta = pi with 5 MeV resolution
    TFile *f7 = new TFile("res_wedge_ff_pi_5.root");
    TH1D *get7 = (TH1D *)f7->Get("FF");
    get7->Scale(197./get7->Integral(), "width");
    get7->SetLineColor(kOrange);
    //get7->SetLineStyle(2);

    // theta = pi with 100 MeV resolution
    TFile *f8 = new TFile("res_wedge_ff_pi_200.root");
    TH1D *get8 = (TH1D *)f8->Get("FF");
    get8->Scale(197./get8->Integral(), "width");
    get8->SetLineColor(kRed);
    //get8->SetLineStyle(2);

    TCanvas *c_f1 = new TCanvas();
    get1->Draw();
    get3->Draw("same");
    get4->Draw("same");
    get5->Draw("same");
    get6->Draw("same");
    get7->Draw("same");
    get8->Draw("same");
    gPad->SetLogy(1);
    gStyle->SetOptStat(0);

    auto legend = new TLegend(0.63,0.21,0.83,0.5);
	legend->AddEntry(get1,"|F (t)|^{2} no resolution","l");
    legend->AddEntry(get7,"|F (t)|^{2} with 5 MeV resolution","l");
    legend->AddEntry(get4,"|F (t)|^{2} with 10 MeV resolution","l");
    legend->AddEntry(get5,"|F (t)|^{2} with 50 MeV resolution","l");
    legend->AddEntry(get3,"|F (t)|^{2} with 100 MeV resolution","l");
    legend->AddEntry(get6,"|F (t)|^{2} with 150 MeV resolution","l");
    legend->AddEntry(get8,"|F (t)|^{2} with 200 MeV resolution","l");
    legend->Draw();
}




                //=============================== Phi saturation model =======================================//

// Draw histogram of saturation plot with wedge cut
void wedge_getFF2t_phi()
{
    int bins = 72.;
    double step = 0.0025;
    double min = 0.;
    double max = 0.18;
    TH1D *FF = new TH1D("FF","Form Factor",bins,min,max);
    FF->Sumw2();
    for(int i=0;i<bins;i++)
    {
        double t = (i+1)*step;
        double q = sqrt(t);
        double val = wedge_FF2_inte(q);
        FF->SetBinContent(i+1,val);
        cout << "t: " << t << endl;
    }
    TFile *FF_WS_Transformations = new TFile("phisat_pio2.root","recreate");
	FF->SetTitle("d#sigma/dt for #phi");
    FF->GetYaxis()->SetTitle("d#sigma/dt [nb/GeV^{2}]");
	FF->GetXaxis()->SetTitle("t [GeV^{2}]");
	FF->Write();
	FF_WS_Transformations->Close();

    TCanvas *c1 = new TCanvas("c1", "True Form Factor", 800, 600);
    FF->Draw();
}

// Draw histogram of saturation plot from scan
void get_phisat() 
{
    const int bins = 72;
    double min = 0;
    double max = .18;
    double x_vals[bins] = {0.00149180044251802,0.00414660644107207,0.00658017860641328,0.00901375077175449,0.0116685567703085,0.0141021289356497,0.0167569349342038,0.019190507099545,
0.0216240792648862,0.0240576514302274,0.0264912235955686,0.0291460295941226,0.031358367926251,0.0340131739248051,0.0364467460901463,0.0388803182554875,0.0413138904208287,
0.0437474625861699,0.0464022685847239,0.0488358407500651,0.0514906467486192,0.0539242189139604,0.0565790249125144,0.0587913632446428,0.0614461692431968,0.063879741408538,
0.0663133135738793,0.0687468857392204,0.0714016917377745,0.0738352639031157,0.0762688360684569,0.078923642067011,0.0813572142323522,0.0837907863976934,0.0864455923962474,
0.0888791645615886,0.0913127367269298,0.093746308892271,0.0961798810576122,0.0988346870561663,0.101268259221507,0.103701831386849,0.10613540355219,0.108347741884318,
0.111223781716085,0.113657353881426,0.115869692213555,0.118745732045322,0.121400538043876,0.123391642542791,0.126267682374558,0.128480020706686,0.130913592872028,0.133568398870582,
0.13578073720271,0.138656777034477,0.140869115366605,0.143523921365159,0.145736259697288,0.148391065695842,0.151045871694396,0.153479443859737,0.155913016025078,0.158346588190419,
0.161001394188973,0.163213732521102,0.165868538519656,0.168302110684997,0.170735682850338,0.173390488848892,0.175824061014234,0.178257633179575};
    double y_vals[bins] = {558853.176136302,224087.526703047,77159.7780282615,21685.4626196854,3859.28092747389,469.329022940005,289.742995521137,760.224527180162,
1470.8807418139,1802.07568265643,1669.9354709361,1434.01276162091,1170.46235898057,841.471980674726,519.487609780666,282.480517526656,122.231121530435,44.2796290689089,
12.7645685962943,8.2906879074936,17.31009647318,34.352592563701,55.6447229185968,73.5686667137878,87.8747301398829,92.4512729862438,94.8281637591125,90.133959553294,
79.3900768010247,68.1741212520722,54.2499745313981,42.0876893247695,31.8336312913207,22.312279620763,14.1287559597991,8.94672119409074,4.99001275499674,2.71340713225333,
1.43847998384123,1.00823457225525,1.11598759598992,1.47546281479756,2.21472290755771,3.15981484859246,4.07292198289417,4.74299606783303,5.523310241139,6.11360281380503,
6.59736568024521,6.59736568024521,6.43200111987633,5.96036388624306,5.523310241139,5.11830428512649,4.50820725397009,3.87130332084099,3.08061332807508,2.64539492870511,
2.15921034425629,1.55230536986868,1.14467925737936,0.888053964582182,0.638442170129982,0.384266614281411,0.29811799949135,0.198608151462336,0.12899755970887,
0.139204998947095,0.150220141959229,0.162106901482951,0.198608151462336,0.237229223974273,};

    TH1D *FF = new TH1D("FF", "Coherent Saturation;|t| [GeV^{2}];d#sigma/dt [nb/GeV^{2}]",bins,min,max);
    FF->Sumw2();
    for (int i=0;i<bins;i++) 
    {
        FF->SetBinContent(i+1,y_vals[i]);
    }
    TFile *FF_WS_Transformations = new TFile("test_phisat.root","recreate");
	FF->SetTitle("Saturation Model: d#sigma/dt");
    FF->GetYaxis()->SetTitle("d#sigma/dt [nb/GeV^{2}]");
	FF->GetXaxis()->SetTitle("t [GeV^{2}]");
	FF->Write();
	FF_WS_Transformations->Close();

    TCanvas *c1 = new TCanvas("c1", "True Form Factor", 800, 600);
    FF->Draw();    
}

// Histograms to compare saturation model with |F(t)|^2 with and without wedge cut
void compare_phisat()
{
    // no wedge: no saturation
    TFile *f2 = new TFile("FF_WS_Transformations_real_FF2_1D_72pt.root");
    TH1D *get2 = (TH1D *)f2->Get("FF");
    get2->Scale(197./get2->Integral(), "width");
    get2->SetLineColor(kBlack);
    get2->SetLineStyle(2);

    // saturation scan
    TFile *f3 = new TFile("test_phisat.root");
    TH1D *get3 = (TH1D *)f3->Get("FF");
    get3->Scale(197./get3->Integral(), "width");
    get3->SetMarkerStyle(4);
    get3->SetMarkerColorAlpha(kBlue,0.35);

    TCanvas *c_f1 = new TCanvas();
    get2->Draw();
    get3->Draw("P HIST same");
    gPad->SetLogy(1);
    gStyle->SetOptStat(0);

    auto legend = new TLegend(0.63,0.71,0.83,0.9);
    legend->AddEntry(get2,"no wedge : no sat","l");
    legend->AddEntry(get3,"with sat","p");
    legend->Draw();
}


    //======= Add resolution to saturation model ======//

// Draw histogram of saturation model with wedge cuts, with resolution
void draw_phisat_wedge()
{
    // Parameters from saturation model scan
    const int bins = 72;
    double step = 0.0025;
    double min = 0.;
    double max = 0.18;
    double x_vals[bins] = {0.00149180044251802,0.00414660644107207,0.00658017860641328,0.00901375077175449,0.0116685567703085,0.0141021289356497,0.0167569349342038,0.019190507099545,
0.0216240792648862,0.0240576514302274,0.0264912235955686,0.0291460295941226,0.031358367926251,0.0340131739248051,0.0364467460901463,0.0388803182554875,0.0413138904208287,
0.0437474625861699,0.0464022685847239,0.0488358407500651,0.0514906467486192,0.0539242189139604,0.0565790249125144,0.0587913632446428,0.0614461692431968,0.063879741408538,
0.0663133135738793,0.0687468857392204,0.0714016917377745,0.0738352639031157,0.0762688360684569,0.078923642067011,0.0813572142323522,0.0837907863976934,0.0864455923962474,
0.0888791645615886,0.0913127367269298,0.093746308892271,0.0961798810576122,0.0988346870561663,0.101268259221507,0.103701831386849,0.10613540355219,0.108347741884318,
0.111223781716085,0.113657353881426,0.115869692213555,0.118745732045322,0.121400538043876,0.123391642542791,0.126267682374558,0.128480020706686,0.130913592872028,0.133568398870582,
0.13578073720271,0.138656777034477,0.140869115366605,0.143523921365159,0.145736259697288,0.148391065695842,0.151045871694396,0.153479443859737,0.155913016025078,0.158346588190419,
0.161001394188973,0.163213732521102,0.165868538519656,0.168302110684997,0.170735682850338,0.173390488848892,0.175824061014234,0.178257633179575};
    double y_vals[bins] = {558853.176136302,224087.526703047,77159.7780282615,21685.4626196854,3859.28092747389,469.329022940005,289.742995521137,760.224527180162,
1470.8807418139,1802.07568265643,1669.9354709361,1434.01276162091,1170.46235898057,841.471980674726,519.487609780666,282.480517526656,122.231121530435,44.2796290689089,
12.7645685962943,8.2906879074936,17.31009647318,34.352592563701,55.6447229185968,73.5686667137878,87.8747301398829,92.4512729862438,94.8281637591125,90.133959553294,
79.3900768010247,68.1741212520722,54.2499745313981,42.0876893247695,31.8336312913207,22.312279620763,14.1287559597991,8.94672119409074,4.99001275499674,2.71340713225333,
1.43847998384123,1.00823457225525,1.11598759598992,1.47546281479756,2.21472290755771,3.15981484859246,4.07292198289417,4.74299606783303,5.523310241139,6.11360281380503,
6.59736568024521,6.59736568024521,6.43200111987633,5.96036388624306,5.523310241139,5.11830428512649,4.50820725397009,3.87130332084099,3.08061332807508,2.64539492870511,
2.15921034425629,1.55230536986868,1.14467925737936,0.888053964582182,0.638442170129982,0.384266614281411,0.29811799949135,0.198608151462336,0.12899755970887,
0.139204998947095,0.150220141959229,0.162106901482951,0.198608151462336,0.237229223974273,};
    TH1D *FF = new TH1D("FF","Form Factor",bins,min,max);
    FF->Sumw2();
    for(int i=0;i<bins;i++)
    {
        double t = x_vals[i];
        double q = sqrt(t);
        double val = testwedge_FF2_inte(q);
        FF->SetBinContent(i+1,val);
        cout << "q: " << q << endl;
    }
    TFile *FF_WS_Transformations = new TFile("phisat_pio9_wRes_200.root","recreate");
	FF->SetTitle("d#sigma/dt for #phi with Resolutiuon");
    FF->GetYaxis()->SetTitle("d#sigma/dt [nb/GeV^{2}]");
	FF->GetXaxis()->SetTitle("t [GeV^{2}]");
	FF->Write();
	FF_WS_Transformations->Close();

    TCanvas *c1 = new TCanvas("c1", "True Form Factor", 800, 600);
    FF->Draw();
}

// Saturation model, theta = pi/2
void compare_phisat_wRes_pio2()
{
    // no resolution
    TFile *f1 = new TFile("phisat_pio2.root");
    TH1D *get1 = (TH1D *)f1->Get("FF");
    get1->Scale(197./get1->Integral(), "width");
    get1->SetTitle("Saturation Model (#phi) with Resolution and Wedge (#theta = #pi/2)");
    get1->GetYaxis()->SetTitle("d#sigma/dt [nb/GeV^{2}]");
    get1->GetXaxis()->SetTitle("|t| [GeV^{2}]");
    get1->GetYaxis()->SetRangeUser(1e-5,1e5);
    get1->SetLineColor(kBlack);

    // with 5 MeV resolution
    TFile *f2 = new TFile("phisat_pio2_wRes_5.root");
    TH1D *get2 = (TH1D *)f2->Get("FF");
    get2->Scale(197./get2->Integral(), "width");
    get2->SetLineColor(kRed);
    get2->SetLineStyle(2);

    // with 10 MeV resolution
    TFile *f3 = new TFile("phisat_pio2_wRes_10.root");
    TH1D *get3 = (TH1D *)f3->Get("FF");
    get3->Scale(197./get3->Integral(), "width");
    get3->SetLineColor(kCyan);
    //get3->SetLineStyle(2);

    // with 50 MeV resolution
    TFile *f4 = new TFile("phisat_pio2_wRes_50.root");
    TH1D *get4 = (TH1D *)f4->Get("FF");
    get4->Scale(197./get4->Integral(), "width");
    get4->SetLineColor(kGreen);
    //get4->SetLineStyle(2);

    // with 100 MeV resolution
    TFile *f5 = new TFile("phisat_pio2_wRes_100.root");
    TH1D *get5 = (TH1D *)f5->Get("FF");
    get5->Scale(197./get5->Integral(), "width");
    get5->SetLineColor(kBlue);
    //get5->SetLineStyle(2);

    // with 150 MeV resolution
    TFile *f6 = new TFile("phisat_pio2_wRes_150.root");
    TH1D *get6 = (TH1D *)f6->Get("FF");
    get6->Scale(197./get6->Integral(), "width");
    get6->SetLineColor(kViolet);
    //get6->SetLineStyle(2);

    // with 200 MeV resolution
    TFile *f7 = new TFile("phisat_pio2_wRes_200.root");
    TH1D *get7 = (TH1D *)f7->Get("FF");
    get7->Scale(197./get7->Integral(), "width");
    get7->SetLineColor(kOrange);
    //get7->SetLineStyle(2);

    TCanvas *c_f1 = new TCanvas();
    get1->Draw();
    get2->Draw("same");
    get3->Draw("same");
    get4->Draw("same");
    get5->Draw("same");
    get6->Draw("same");
    get7->Draw("same");
    gPad->SetLogy(1);
    gStyle->SetOptStat(0);

    auto legend = new TLegend(0.43,0.71,0.63,0.9);
	legend->AddEntry(get1,"no resolution","l");
    legend->AddEntry(get2,"5 MeV resolution","l");
    legend->AddEntry(get3,"10 MeV resolution","l");
    legend->AddEntry(get4,"50 MeV resolution","l");
    legend->AddEntry(get5,"100 MeV resolution","l");
    legend->AddEntry(get6,"150 MeV resolution","l");
    legend->AddEntry(get7,"200 MeV resolution","l");
    legend->Draw();
}

// Saturation model, theta = pi/3
void compare_phisat_wRes_pio3()
{
    // no resolution
    TFile *f1 = new TFile("phisat_pio3.root");
    TH1D *get1 = (TH1D *)f1->Get("FF");
    get1->Scale(197./get1->Integral(), "width");
    get1->SetTitle("Saturation Model (#phi) with Resolution and Wedge Cut (#theta = #pi/3)");
    get1->GetYaxis()->SetTitle("d#sigma/dt [nb/GeV^{2}]");
    get1->GetXaxis()->SetTitle("|t| [GeV^{2}]");
    get1->GetYaxis()->SetRangeUser(1e-5,1e5);
    get1->SetLineColor(kBlack);

    // with 5 MeV resolution
    TFile *f2 = new TFile("phisat_pio3_wRes_5.root");
    TH1D *get2 = (TH1D *)f2->Get("FF");
    get2->Scale(197./get2->Integral(), "width");
    get2->SetLineColor(kRed);
    get2->SetLineStyle(2);

    // with 10 MeV resolution
    TFile *f3 = new TFile("phisat_pio3_wRes_10.root");
    TH1D *get3 = (TH1D *)f3->Get("FF");
    get3->Scale(197./get3->Integral(), "width");
    get3->SetLineColor(kCyan);
    //get3->SetLineStyle(2);

    // with 50 MeV resolution
    TFile *f4 = new TFile("phisat_pio3_wRes_50.root");
    TH1D *get4 = (TH1D *)f4->Get("FF");
    get4->Scale(197./get4->Integral(), "width");
    get4->SetLineColor(kGreen);
    //get4->SetLineStyle(2);

    // with 100 MeV resolution
    TFile *f5 = new TFile("phisat_pio3_wRes_100.root");
    TH1D *get5 = (TH1D *)f5->Get("FF");
    get5->Scale(197./get5->Integral(), "width");
    get5->SetLineColor(kBlue);
    //get5->SetLineStyle(2);

    // with 150 MeV resolution
    TFile *f6 = new TFile("phisat_pio3_wRes_150.root");
    TH1D *get6 = (TH1D *)f6->Get("FF");
    get6->Scale(197./get6->Integral(), "width");
    get6->SetLineColor(kViolet);
    //get6->SetLineStyle(2);

    // with 200 MeV resolution
    TFile *f7 = new TFile("phisat_pio3_wRes_200.root");
    TH1D *get7 = (TH1D *)f7->Get("FF");
    get7->Scale(197./get7->Integral(), "width");
    get7->SetLineColor(kOrange);
    //get7->SetLineStyle(2);

    TCanvas *c_f1 = new TCanvas();
    get1->Draw();
    get2->Draw("same");
    get3->Draw("same");
    get4->Draw("same");
    get5->Draw("same");
    get6->Draw("same");
    get7->Draw("same");
    gPad->SetLogy(1);
    gStyle->SetOptStat(0);

    auto legend = new TLegend(0.43,0.71,0.63,0.9);
	legend->AddEntry(get1,"no resolution","l");
    legend->AddEntry(get2,"5 MeV resolution","l");
    legend->AddEntry(get3,"10 MeV resolution","l");
    legend->AddEntry(get4,"50 MeV resolution","l");
    legend->AddEntry(get5,"100 MeV resolution","l");
    legend->AddEntry(get6,"150 MeV resolution","l");
    legend->AddEntry(get7,"200 MeV resolution","l");
    legend->Draw();
}

// Saturation model, theta = pi/4
void compare_phisat_wRes_pio4()
{
    // no resolution
    TFile *f1 = new TFile("phisat_pio4.root");
    TH1D *get1 = (TH1D *)f1->Get("FF");
    get1->Scale(197./get1->Integral(), "width");
    get1->SetTitle("Saturation Model (#phi) with Resolution and Wedge Cut (#theta = #pi/4)");
    get1->GetYaxis()->SetTitle("d#sigma/dt [nb/GeV^{2}]");
    get1->GetXaxis()->SetTitle("|t| [GeV^{2}]");
    get1->GetYaxis()->SetRangeUser(1e-5,1e5);
    get1->SetLineColor(kBlack);

    // with 5 MeV resolution
    TFile *f2 = new TFile("phisat_pio4_wRes_5.root");
    TH1D *get2 = (TH1D *)f2->Get("FF");
    get2->Scale(197./get2->Integral(), "width");
    get2->SetLineColor(kRed);
    get2->SetLineStyle(2);

    // with 10 MeV resolution
    TFile *f3 = new TFile("phisat_pio4_wRes_10.root");
    TH1D *get3 = (TH1D *)f3->Get("FF");
    get3->Scale(197./get3->Integral(), "width");
    get3->SetLineColor(kCyan);
    //get3->SetLineStyle(2);

    // with 50 MeV resolution
    TFile *f4 = new TFile("phisat_pio4_wRes_50.root");
    TH1D *get4 = (TH1D *)f4->Get("FF");
    get4->Scale(197./get4->Integral(), "width");
    get4->SetLineColor(kGreen);
    //get4->SetLineStyle(2);

    // with 100 MeV resolution
    TFile *f5 = new TFile("phisat_pio4_wRes_100.root");
    TH1D *get5 = (TH1D *)f5->Get("FF");
    get5->Scale(197./get5->Integral(), "width");
    get5->SetLineColor(kBlue);
    //get5->SetLineStyle(2);

    // with 150 MeV resolution
    TFile *f6 = new TFile("phisat_pio4_wRes_150.root");
    TH1D *get6 = (TH1D *)f6->Get("FF");
    get6->Scale(197./get6->Integral(), "width");
    get6->SetLineColor(kViolet);
    //get6->SetLineStyle(2);

    // with 200 MeV resolution
    TFile *f7 = new TFile("phisat_pio4_wRes_200.root");
    TH1D *get7 = (TH1D *)f7->Get("FF");
    get7->Scale(197./get7->Integral(), "width");
    get7->SetLineColor(kOrange);
    //get7->SetLineStyle(2);

    TCanvas *c_f1 = new TCanvas();
    get1->Draw();
    get2->Draw("same");
    get3->Draw("same");
    get4->Draw("same");
    get5->Draw("same");
    get6->Draw("same");
    get7->Draw("same");
    gPad->SetLogy(1);
    gStyle->SetOptStat(0);

    auto legend = new TLegend(0.43,0.71,0.63,0.9);
	legend->AddEntry(get1,"no resolution","l");
    legend->AddEntry(get2,"5 MeV resolution","l");
    legend->AddEntry(get3,"10 MeV resolution","l");
    legend->AddEntry(get4,"50 MeV resolution","l");
    legend->AddEntry(get5,"100 MeV resolution","l");
    legend->AddEntry(get6,"150 MeV resolution","l");
    legend->AddEntry(get7,"200 MeV resolution","l");
    legend->Draw();
}

// Saturation model, theta = pi/6
void compare_phisat_wRes_pio6()
{
    // no resolution
    TFile *f1 = new TFile("phisat_pio6.root");
    TH1D *get1 = (TH1D *)f1->Get("FF");
    get1->Scale(197./get1->Integral(), "width");
    get1->SetTitle("Saturation Model (#phi) with Resolution and Wedge Cut (#theta = #pi/6)");
    get1->GetYaxis()->SetTitle("d#sigma/dt [nb/GeV^{2}]");
    get1->GetXaxis()->SetTitle("|t| [GeV^{2}]");
    get1->GetYaxis()->SetRangeUser(1e-5,1e5);
    get1->SetLineColor(kBlack);

    // with 5 MeV resolution
    TFile *f2 = new TFile("phisat_pio6_wRes_5.root");
    TH1D *get2 = (TH1D *)f2->Get("FF");
    get2->Scale(197./get2->Integral(), "width");
    get2->SetLineColor(kRed);
    get2->SetLineStyle(2);

    // with 10 MeV resolution
    TFile *f3 = new TFile("phisat_pio6_wRes_10.root");
    TH1D *get3 = (TH1D *)f3->Get("FF");
    get3->Scale(197./get3->Integral(), "width");
    get3->SetLineColor(kCyan);
    //get3->SetLineStyle(2);

    // with 50 MeV resolution
    TFile *f4 = new TFile("phisat_pio6_wRes_50.root");
    TH1D *get4 = (TH1D *)f4->Get("FF");
    get4->Scale(197./get4->Integral(), "width");
    get4->SetLineColor(kGreen);
    //get4->SetLineStyle(2);

    // with 100 MeV resolution
    TFile *f5 = new TFile("phisat_pio6_wRes_100.root");
    TH1D *get5 = (TH1D *)f5->Get("FF");
    get5->Scale(197./get5->Integral(), "width");
    get5->SetLineColor(kBlue);
    //get5->SetLineStyle(2);

    // with 150 MeV resolution
    TFile *f6 = new TFile("phisat_pio6_wRes_150.root");
    TH1D *get6 = (TH1D *)f6->Get("FF");
    get6->Scale(197./get6->Integral(), "width");
    get6->SetLineColor(kViolet);
    //get6->SetLineStyle(2);

    // with 200 MeV resolution
    TFile *f7 = new TFile("phisat_pio6_wRes_200.root");
    TH1D *get7 = (TH1D *)f7->Get("FF");
    get7->Scale(197./get7->Integral(), "width");
    get7->SetLineColor(kOrange);
    //get7->SetLineStyle(2);

    TCanvas *c_f1 = new TCanvas();
    get1->Draw();
    get2->Draw("same");
    get3->Draw("same");
    get4->Draw("same");
    get5->Draw("same");
    get6->Draw("same");
    get7->Draw("same");
    gPad->SetLogy(1);
    gStyle->SetOptStat(0);

    auto legend = new TLegend(0.43,0.71,0.63,0.9);
	legend->AddEntry(get1,"no resolution","l");
    legend->AddEntry(get2,"5 MeV resolution","l");
    legend->AddEntry(get3,"10 MeV resolution","l");
    legend->AddEntry(get4,"50 MeV resolution","l");
    legend->AddEntry(get5,"100 MeV resolution","l");
    legend->AddEntry(get6,"150 MeV resolution","l");
    legend->AddEntry(get7,"200 MeV resolution","l");
    legend->Draw();
}

// Saturation model, theta = pi/9
void compare_phisat_wRes_pio9()
{
    // no resolution
    TFile *f1 = new TFile("phisat_pio9.root");
    TH1D *get1 = (TH1D *)f1->Get("FF");
    get1->Scale(197./get1->Integral(), "width");
    get1->SetTitle("Saturation Model (#phi) with Resolution and Wedge Cut (#theta = #pi/9)");
    get1->GetYaxis()->SetTitle("d#sigma/dt [nb/GeV^{2}]");
    get1->GetXaxis()->SetTitle("|t| [GeV^{2}]");
    get1->GetYaxis()->SetRangeUser(1e-5,1e5);
    get1->SetLineColor(kBlack);

    // with 5 MeV resolution
    TFile *f2 = new TFile("phisat_pio9_wRes_5.root");
    TH1D *get2 = (TH1D *)f2->Get("FF");
    get2->Scale(197./get2->Integral(), "width");
    get2->SetLineColor(kRed);
    get2->SetLineStyle(2);

    // with 10 MeV resolution
    TFile *f3 = new TFile("phisat_pio9_wRes_10.root");
    TH1D *get3 = (TH1D *)f3->Get("FF");
    get3->Scale(197./get3->Integral(), "width");
    get3->SetLineColor(kCyan);
    //get3->SetLineStyle(2);

    // with 50 MeV resolution
    TFile *f4 = new TFile("phisat_pio9_wRes_50.root");
    TH1D *get4 = (TH1D *)f4->Get("FF");
    get4->Scale(197./get4->Integral(), "width");
    get4->SetLineColor(kGreen);
    //get4->SetLineStyle(2);

    // with 100 MeV resolution
    TFile *f5 = new TFile("phisat_pio9_wRes_100.root");
    TH1D *get5 = (TH1D *)f5->Get("FF");
    get5->Scale(197./get5->Integral(), "width");
    get5->SetLineColor(kBlue);
    //get5->SetLineStyle(2);

    // with 150 MeV resolution
    TFile *f6 = new TFile("phisat_pio9_wRes_150.root");
    TH1D *get6 = (TH1D *)f6->Get("FF");
    get6->Scale(197./get6->Integral(), "width");
    get6->SetLineColor(kViolet);
    //get6->SetLineStyle(2);

    // with 200 MeV resolution
    TFile *f7 = new TFile("phisat_pio9_wRes_200.root");
    TH1D *get7 = (TH1D *)f7->Get("FF");
    get7->Scale(197./get7->Integral(), "width");
    get7->SetLineColor(kOrange);
    //get7->SetLineStyle(2);

    TCanvas *c_f1 = new TCanvas();
    get1->Draw();
    get2->Draw("same");
    get3->Draw("same");
    get4->Draw("same");
    get5->Draw("same");
    get6->Draw("same");
    get7->Draw("same");
    gPad->SetLogy(1);
    gStyle->SetOptStat(0);

    auto legend = new TLegend(0.43,0.71,0.63,0.9);
	legend->AddEntry(get1,"no resolution","l");
    legend->AddEntry(get2,"5 MeV resolution","l");
    legend->AddEntry(get3,"10 MeV resolution","l");
    legend->AddEntry(get4,"50 MeV resolution","l");
    legend->AddEntry(get5,"100 MeV resolution","l");
    legend->AddEntry(get6,"150 MeV resolution","l");
    legend->AddEntry(get7,"200 MeV resolution","l");
    legend->Draw();
}

// Saturation model, theta = pi
void compare_phisat_wRes_pi()
{
    // no resolution
    TFile *f1 = new TFile("phisat_pi.root");
    TH1D *get1 = (TH1D *)f1->Get("FF");
    get1->Scale(197./get1->Integral(), "width");
    get1->SetTitle("Saturation Model (#phi) with Resolution and Wedge Cut (#theta = #pi)");
    get1->GetYaxis()->SetTitle("d#sigma/dt [nb/GeV^{2}]");
    get1->GetXaxis()->SetTitle("|t| [GeV^{2}]");
    get1->GetYaxis()->SetRangeUser(1e-5,1e5);
    get1->SetLineColor(kBlack);

    // with 5 MeV resolution
    TFile *f2 = new TFile("phisat_pi_wRes_5.root");
    TH1D *get2 = (TH1D *)f2->Get("FF");
    get2->Scale(197./get2->Integral(), "width");
    get2->SetLineColor(kRed);
    get2->SetLineStyle(2);

    // with 10 MeV resolution
    TFile *f3 = new TFile("phisat_pi_wRes_10.root");
    TH1D *get3 = (TH1D *)f3->Get("FF");
    get3->Scale(197./get3->Integral(), "width");
    get3->SetLineColor(kCyan);
    //get3->SetLineStyle(2);

    // with 50 MeV resolution
    TFile *f4 = new TFile("phisat_pi_wRes_50.root");
    TH1D *get4 = (TH1D *)f4->Get("FF");
    get4->Scale(197./get4->Integral(), "width");
    get4->SetLineColor(kGreen);
    //get4->SetLineStyle(2);

    // with 100 MeV resolution
    TFile *f5 = new TFile("phisat_pi_wRes_100.root");
    TH1D *get5 = (TH1D *)f5->Get("FF");
    get5->Scale(197./get5->Integral(), "width");
    get5->SetLineColor(kBlue);
    //get5->SetLineStyle(2);

    // with 150 MeV resolution
    TFile *f6 = new TFile("phisat_pi_wRes_150.root");
    TH1D *get6 = (TH1D *)f6->Get("FF");
    get6->Scale(197./get6->Integral(), "width");
    get6->SetLineColor(kViolet);
    //get6->SetLineStyle(2);

    // with 200 MeV resolution
    TFile *f7 = new TFile("phisat_pi_wRes_200.root");
    TH1D *get7 = (TH1D *)f7->Get("FF");
    get7->Scale(197./get7->Integral(), "width");
    get7->SetLineColor(kOrange);
    //get7->SetLineStyle(2);

    TCanvas *c_f1 = new TCanvas();
    get1->Draw();
    get2->Draw("same");
    get3->Draw("same");
    get4->Draw("same");
    get5->Draw("same");
    get6->Draw("same");
    get7->Draw("same");
    gPad->SetLogy(1);
    gStyle->SetOptStat(0);

    auto legend = new TLegend(0.43,0.71,0.63,0.9);
	legend->AddEntry(get1,"no resolution","l");
    legend->AddEntry(get2,"5 MeV resolution","l");
    legend->AddEntry(get3,"10 MeV resolution","l");
    legend->AddEntry(get4,"50 MeV resolution","l");
    legend->AddEntry(get5,"100 MeV resolution","l");
    legend->AddEntry(get6,"150 MeV resolution","l");
    legend->AddEntry(get7,"200 MeV resolution","l");
    legend->Draw();
}

// Saturation model, theta = 2pi
void compare_phisat_wRes_2pi()
{
    // no resolution
    TFile *f1 = new TFile("phisat_2pi.root");
    TH1D *get1 = (TH1D *)f1->Get("FF");
    get1->Scale(197./get1->Integral(), "width");
    get1->SetTitle("Saturation Model (#phi) with Resolution and Wedge Cut (#theta = 2#pi)");
    get1->GetYaxis()->SetTitle("d#sigma/dt [nb/GeV^{2}]");
    get1->GetXaxis()->SetTitle("|t| [GeV^{2}]");
    get1->GetYaxis()->SetRangeUser(1e-5,1e5);
    get1->SetLineColor(kBlack);

    // with 5 MeV resolution
    TFile *f2 = new TFile("phisat_2pi_wRes_5.root");
    TH1D *get2 = (TH1D *)f2->Get("FF");
    get2->Scale(197./get2->Integral(), "width");
    get2->SetLineColor(kRed);
    get2->SetLineStyle(2);

    // with 10 MeV resolution
    TFile *f3 = new TFile("phisat_2pi_wRes_10.root");
    TH1D *get3 = (TH1D *)f3->Get("FF");
    get3->Scale(197./get3->Integral(), "width");
    get3->SetLineColor(kCyan);
    //get3->SetLineStyle(2);

    // with 50 MeV resolution
    TFile *f4 = new TFile("phisat_2pi_wRes_50.root");
    TH1D *get4 = (TH1D *)f4->Get("FF");
    get4->Scale(197./get4->Integral(), "width");
    get4->SetLineColor(kGreen);
    //get4->SetLineStyle(2);

    // with 100 MeV resolution
    TFile *f5 = new TFile("phisat_2pi_wRes_100.root");
    TH1D *get5 = (TH1D *)f5->Get("FF");
    get5->Scale(197./get5->Integral(), "width");
    get5->SetLineColor(kBlue);
    //get5->SetLineStyle(2);

    // with 150 MeV resolution
    TFile *f6 = new TFile("phisat_2pi_wRes_150.root");
    TH1D *get6 = (TH1D *)f6->Get("FF");
    get6->Scale(197./get6->Integral(), "width");
    get6->SetLineColor(kViolet);
    //get6->SetLineStyle(2);

    // with 200 MeV resolution
    TFile *f7 = new TFile("phisat_2pi_wRes_200.root");
    TH1D *get7 = (TH1D *)f7->Get("FF");
    get7->Scale(197./get7->Integral(), "width");
    get7->SetLineColor(kOrange);
    //get7->SetLineStyle(2);

    TCanvas *c_f1 = new TCanvas();
    get1->Draw();
    get2->Draw("same");
    get3->Draw("same");
    get4->Draw("same");
    get5->Draw("same");
    get6->Draw("same");
    get7->Draw("same");
    gPad->SetLogy(1);
    gStyle->SetOptStat(0);

    auto legend = new TLegend(0.43,0.71,0.63,0.9);
	legend->AddEntry(get1,"no resolution","l");
    legend->AddEntry(get2,"5 MeV resolution","l");
    legend->AddEntry(get3,"10 MeV resolution","l");
    legend->AddEntry(get4,"50 MeV resolution","l");
    legend->AddEntry(get5,"100 MeV resolution","l");
    legend->AddEntry(get6,"150 MeV resolution","l");
    legend->AddEntry(get7,"200 MeV resolution","l");
    legend->Draw();
}

// Saturation model, theta = 0 >> (1e-12)
void compare_phisat_wRes_zero()
{
    // no resolution
    TFile *f1 = new TFile("phisat_0.root");
    TH1D *get1 = (TH1D *)f1->Get("FF");
    get1->Scale(197./get1->Integral(), "width");
    get1->SetTitle("Saturation Model (#phi) with Resolution and Wedge Cut (#theta = 0)");
    get1->GetYaxis()->SetTitle("d#sigma/dt [nb/GeV^{2}]");
    get1->GetXaxis()->SetTitle("|t| [GeV^{2}]");
    get1->GetYaxis()->SetRangeUser(1e-5,1e5);
    get1->SetLineColor(kBlack);

    // with 5 MeV resolution
    TFile *f2 = new TFile("phisat_0_wRes_5.root");
    TH1D *get2 = (TH1D *)f2->Get("FF");
    get2->Scale(197./get2->Integral(), "width");
    get2->SetLineColor(kRed);
    get2->SetLineStyle(2);

    // with 10 MeV resolution
    TFile *f3 = new TFile("phisat_0_wRes_10.root");
    TH1D *get3 = (TH1D *)f3->Get("FF");
    get3->Scale(197./get3->Integral(), "width");
    get3->SetLineColor(kCyan);
    //get3->SetLineStyle(2);

    // with 50 MeV resolution
    TFile *f4 = new TFile("phisat_0_wRes_50.root");
    TH1D *get4 = (TH1D *)f4->Get("FF");
    get4->Scale(197./get4->Integral(), "width");
    get4->SetLineColor(kGreen);
    //get4->SetLineStyle(2);

    // with 100 MeV resolution
    TFile *f5 = new TFile("phisat_0_wRes_100.root");
    TH1D *get5 = (TH1D *)f5->Get("FF");
    get5->Scale(197./get5->Integral(), "width");
    get5->SetLineColor(kBlue);
    //get5->SetLineStyle(2);

    // with 150 MeV resolution
    TFile *f6 = new TFile("phisat_0_wRes_150.root");
    TH1D *get6 = (TH1D *)f6->Get("FF");
    get6->Scale(197./get6->Integral(), "width");
    get6->SetLineColor(kViolet);
    //get6->SetLineStyle(2);

    // with 200 MeV resolution
    TFile *f7 = new TFile("phisat_0_wRes_200.root");
    TH1D *get7 = (TH1D *)f7->Get("FF");
    get7->Scale(197./get7->Integral(), "width");
    get7->SetLineColor(kOrange);
    //get7->SetLineStyle(2);

    TCanvas *c_f1 = new TCanvas();
    get1->Draw();
    get2->Draw("same");
    get3->Draw("same");
    get4->Draw("same");
    get5->Draw("same");
    get6->Draw("same");
    get7->Draw("same");
    gPad->SetLogy(1);
    gStyle->SetOptStat(0);

    auto legend = new TLegend(0.43,0.71,0.63,0.9);
	legend->AddEntry(get1,"no resolution","l");
    legend->AddEntry(get2,"5 MeV resolution","l");
    legend->AddEntry(get3,"10 MeV resolution","l");
    legend->AddEntry(get4,"50 MeV resolution","l");
    legend->AddEntry(get5,"100 MeV resolution","l");
    legend->AddEntry(get6,"150 MeV resolution","l");
    legend->AddEntry(get7,"200 MeV resolution","l");
    legend->Draw();
}




    //=== Histograms of each |F(t)|^2 with each angle, each resolution, and phi saturation ===//

void phipio2_Wsat()
{
    // theta = pi/2
    TFile *f1 = new TFile("72pt_wedge_ff_pio2.root");
    TH1D *get1 = (TH1D *)f1->Get("FF");
    get1->Scale(197./get1->Integral(), "width");
    get1->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/2)");
    get1->GetYaxis()->SetTitle("");
    get1->GetXaxis()->SetTitle("t [GeV^{2}]");
    get1->Scale(197./get1->Integral(), "width");
    get1->GetYaxis()->SetRangeUser(1e-5,1e5);
    get1->SetLineColor(kBlack);

    // theta = pi/2 with 100 MeV resolution
    TFile *f3 = new TFile("72pt_res_wedge_ff_pio2_100.root");
    TH1D *get3 = (TH1D *)f3->Get("FF");
    get3->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/2)");
    get3->GetYaxis()->SetTitle("");
    get3->Scale(197./get3->Integral(), "width");
    get3->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get3->SetLineColor(kBlue);
    //get3->SetLineStyle(2);
    get3->SetMarkerStyle(3);
    get3->SetMarkerSize(0.7);
    get3->SetMarkerColorAlpha(kBlack,0.35);

    // theta = pi/2 with 50 MeV resolution
    TFile *f4 = new TFile("72pt_res_wedge_ff_pio2_50.root");
    TH1D *get4 = (TH1D *)f4->Get("FF");
    get4->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/2)");
    get4->GetYaxis()->SetTitle("");
    get4->Scale(197./get4->Integral(), "width");
    get4->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get4->SetLineColor(kGreen);
    //get4->SetLineStyle(2);
    get4->SetMarkerStyle(3);
    get4->SetMarkerSize(0.7);
    get4->SetMarkerColorAlpha(kBlack,0.35);

    // theta = pi/2 with 10 MeV resolution
    TFile *f5 = new TFile("72pt_res_wedge_ff_pio2_10.root");
    TH1D *get5 = (TH1D *)f5->Get("FF");
    get5->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/2)");
    get5->GetYaxis()->SetTitle("");
    get5->Scale(197./get5->Integral(), "width");
    get5->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get5->SetLineColor(kBlack);
    //get5->SetLineStyle(1);
    get5->SetMarkerStyle(3);
    get5->SetMarkerSize(0.7);
    get5->SetMarkerColorAlpha(kBlack,0.35);

    // theta = pi/2 with 150 MeV resolution
    TFile *f6 = new TFile("72pt_res_wedge_ff_pio2_150.root");
    TH1D *get6 = (TH1D *)f6->Get("FF");
    get6->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/2)");
    get6->GetYaxis()->SetTitle("");
    get6->Scale(197./get6->Integral(), "width");
    get6->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get6->SetLineColor(kViolet);
    get6->SetMarkerStyle(3);
    get6->SetMarkerSize(0.7);
    get6->SetMarkerColorAlpha(kBlack,0.35);

    // theta = pi/2 with 5 MeV resolution
    TFile *f7 = new TFile("72pt_res_wedge_ff_pio2_5.root");
    TH1D *get7 = (TH1D *)f7->Get("FF");
    get7->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/2)");
    get7->GetYaxis()->SetTitle("");
    get7->Scale(197./get7->Integral(), "width");
    get7->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get7->SetLineColor(kOrange);
    //get7->SetLineStyle(2);
    get7->SetMarkerStyle(3);
    get7->SetMarkerSize(0.7);
    get7->SetMarkerColorAlpha(kBlack,0.35);

    // theta = pi/2 with 200 MeV resolution
    TFile *f71 = new TFile("72pt_res_wedge_ff_pio2_200.root");
    TH1D *get71 = (TH1D *)f71->Get("FF");
    get71->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/2)");
    get71->GetYaxis()->SetTitle("");
    get71->Scale(197./get71->Integral(), "width");
    get71->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get71->SetLineColor(kOrange);
    //get71->SetLineStyle(2);
    get71->SetMarkerStyle(3);
    get71->SetMarkerSize(0.7);
    get71->SetMarkerColorAlpha(kBlack,0.35);

    // with saturation
    TFile *f8 = new TFile("phisat_pio2.root");
    TH1D *get8 = (TH1D *)f8->Get("FF");
    get8->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/2)");
    get8->GetYaxis()->SetTitle("");
    get8->Scale(197./get8->Integral(), "width");
    get8->GetYaxis()->SetRangeUser(1e-5,1e5);
    get8->SetLineColor(kRed);
    get8->SetLineStyle(1);
    //get8->SetMarkerStyle(4);
    //get8->SetMarkerSize(0.6);
    //get8->SetMarkerColorAlpha(kRed,0.35);

    // Saturation with 5 MeV resolution
    TFile *f22 = new TFile("phisat_pio2_wRes_5.root");
    TH1D *get22 = (TH1D *)f22->Get("FF");
    get22->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/2)");
    get22->GetYaxis()->SetTitle("");
    get22->Scale(197./get22->Integral(), "width");
    get22->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get22->SetLineColor(kOrange);
    //get22->SetLineStyle(2);
    get22->SetMarkerStyle(4);
    get22->SetMarkerSize(0.7);
    get22->SetMarkerColorAlpha(kRed,0.35);

    // Saturation with 10 MeV resolution
    TFile *f33 = new TFile("phisat_pio2_wRes_10.root");
    TH1D *get33 = (TH1D *)f33->Get("FF");
    get33->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/2)");
    get33->GetYaxis()->SetTitle("");
    get33->Scale(197./get33->Integral(), "width");
    get33->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get33->SetLineColor(kCyan);
    //get33->SetLineStyle(2);
    get33->SetMarkerStyle(4);
    get33->SetMarkerSize(0.7);
    get33->SetMarkerColorAlpha(kRed,0.35);

    // Saturation with 50 MeV resolution
    TFile *f44 = new TFile("phisat_pio2_wRes_50.root");
    TH1D *get44 = (TH1D *)f44->Get("FF");
    get44->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/2)");
    get44->GetYaxis()->SetTitle("");
    get44->Scale(197./get44->Integral(), "width");
    get44->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get44->SetLineColor(kGreen);
    //get44->SetLineStyle(2);
    get44->SetMarkerStyle(4);
    get44->SetMarkerSize(0.7);
    get44->SetMarkerColorAlpha(kRed,0.35);

    // Saturation with 100 MeV resolution
    TFile *f55 = new TFile("phisat_pio2_wRes_100.root");
    TH1D *get55 = (TH1D *)f55->Get("FF");
    get55->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/2)");
    get55->GetYaxis()->SetTitle("");
    get55->Scale(197./get55->Integral(), "width");
    get55->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get55->SetLineColor(kBlue);
    //get55->SetLineStyle(2);
    get55->SetMarkerStyle(4);
    get55->SetMarkerSize(0.7);
    get55->SetMarkerColorAlpha(kRed,0.35);

    // Saturation with 150 MeV resolution
    TFile *f66 = new TFile("phisat_pio2_wRes_150.root");
    TH1D *get66 = (TH1D *)f66->Get("FF");
    get66->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/2)");
    get66->GetYaxis()->SetTitle("");
    get66->Scale(197./get66->Integral(), "width");
    get66->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get66->SetLineColor(kViolet);
    //get66->SetLineStyle(2);
    get66->SetMarkerStyle(4);
    get66->SetMarkerSize(0.7);
    get66->SetMarkerColorAlpha(kRed,0.35);

    // Saturation with 200 MeV resolution
    TFile *f61 = new TFile("phisat_pio2_wRes_200.root");
    TH1D *get61 = (TH1D *)f61->Get("FF");
    get61->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/2)");
    get61->GetYaxis()->SetTitle("");
    get61->Scale(197./get61->Integral(), "width");
    get61->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get61->SetLineColor(kViolet);
    //get61->SetLineStyle(2);
    get61->SetMarkerStyle(4);
    get61->SetMarkerSize(0.7);
    get61->SetMarkerColorAlpha(kRed,0.35);

    TCanvas *c_f1 = new TCanvas();
    //get1->Draw();
    //get3->Draw("P HIST same");
    //get4->Draw("P HIST same");
    //get5->Draw("P HIST same");
    //get6->Draw("P HIST same");
    get7->Draw("P HIST same");
    //get71->Draw("P HIST same");
    //get8->Draw("same");
    get22->Draw("P HIST same");
    //get33->Draw("P HIST same");
    //get44->Draw("P HIST same");
    //get55->Draw("P HIST same");
    //get66->Draw("P HIST same");
    //get61->Draw("P HIST same");
    gPad->SetLogy(1);
    gStyle->SetOptStat(0);

    auto legend = new TLegend(0.43,0.81,0.73,0.9);
    //legend->AddEntry(get8,"d#sigma/dt [nb/GeV^{2}] with sat: no res","l");
    legend->AddEntry(get22,"d#sigma/dt [nb/GeV^{2}] with sat: 5 MeV res","p");
    //legend->AddEntry(get33,"d#sigma/dt [nb/GeV^{2}] with sat: 10 MeV res","p");
    //legend->AddEntry(get44,"d#sigma/dt [nb/GeV^{2}] with sat: 50 MeV res","p");
    //legend->AddEntry(get55,"d#sigma/dt [nb/GeV^{2}] with sat: 100 MeV res","p");
    //legend->AddEntry(get66,"d#sigma/dt [nb/GeV^{2}] with sat: 150 MeV res","p");
    //legend->AddEntry(get61,"d#sigma/dt [nb/GeV^{2}] with sat: 200 MeV res","p");
	//legend->AddEntry(get1,"|F(t)|^{2} no res","l");
    legend->AddEntry(get7,"|F(t)|^{2} with 5 MeV resolution","p");
    //legend->AddEntry(get5,"|F(t)|^{2} with 10 MeV res","p");
    //legend->AddEntry(get4,"|F(t)|^{2} with 50 MeV resolution","p");
    //legend->AddEntry(get3,"|F(t)|^{2} with 100 MeV resolution","p");
    //legend->AddEntry(get6,"|F(t)|^{2} with 150 MeV resolution","p");
    //legend->AddEntry(get71,"|F(t)|^{2} with 200 MeV resolution","p");
    legend->Draw();
}

void phipio4_Wsat()
{
    // theta = pi/4
    TFile *f1 = new TFile("72pt_wedge_ff_pio4.root");
    TH1D *get1 = (TH1D *)f1->Get("FF");
    get1->Scale(197./get1->Integral(), "width");
    get1->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/4)");
    get1->GetYaxis()->SetTitle("");
    get1->GetXaxis()->SetTitle("t [GeV^{2}]");
    get1->Scale(197./get1->Integral(), "width");
    get1->GetYaxis()->SetRangeUser(1e-5,1e5);
    get1->SetLineColor(kBlack);

    // theta = pi/4 with 100 MeV resolution
    TFile *f3 = new TFile("72pt_res_wedge_ff_pio4_100.root");
    TH1D *get3 = (TH1D *)f3->Get("FF");
    get3->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/4)");
    get3->GetYaxis()->SetTitle("");
    get3->Scale(197./get3->Integral(), "width");
    get3->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get3->SetLineColor(kBlue);
    //get3->SetLineStyle(2);
    get3->SetMarkerStyle(3);
    get3->SetMarkerSize(0.7);
    get3->SetMarkerColorAlpha(kBlack,0.35);

    // theta = pi/4 with 50 MeV resolution
    TFile *f4 = new TFile("72pt_res_wedge_ff_pio4_50.root");
    TH1D *get4 = (TH1D *)f4->Get("FF");
    get4->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/4)");
    get4->GetYaxis()->SetTitle("");
    get4->Scale(197./get4->Integral(), "width");
    get4->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get4->SetLineColor(kGreen);
    //get4->SetLineStyle(2);
    get4->SetMarkerStyle(3);
    get4->SetMarkerSize(0.7);
    get4->SetMarkerColorAlpha(kBlack,0.35);

    // theta = pi/4 with 10 MeV resolution
    TFile *f5 = new TFile("72pt_res_wedge_ff_pio4_10.root");
    TH1D *get5 = (TH1D *)f5->Get("FF");
    get5->SetTitle("Wedge Integration ans Saturation Model (#phi) with Resolution (#theta_{max} = #pi/4)");
    get5->GetYaxis()->SetTitle("");
    get5->Scale(197./get5->Integral(), "width");
    get5->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get5->SetLineColor(kCyan);
    //get5->SetLineStyle(2);
    get5->SetMarkerStyle(3);
    get5->SetMarkerSize(0.7);
    get5->SetMarkerColorAlpha(kBlack,0.35);

    // theta = pi/4 with 150 MeV resolution
    TFile *f6 = new TFile("72pt_res_wedge_ff_pio4_150.root");
    TH1D *get6 = (TH1D *)f6->Get("FF");
    get6->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/4)");
    get6->GetYaxis()->SetTitle("");
    get6->Scale(197./get6->Integral(), "width");
    get6->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get6->SetLineColor(kViolet);
    get6->SetMarkerStyle(3);
    get6->SetMarkerSize(0.7);
    get6->SetMarkerColorAlpha(kBlack,0.35);

    // theta = pi/4 with 5 MeV resolution
    TFile *f7 = new TFile("72pt_res_wedge_ff_pio4_5.root");
    TH1D *get7 = (TH1D *)f7->Get("FF");
    get7->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/4)");
    get7->GetYaxis()->SetTitle("");
    get7->Scale(197./get7->Integral(), "width");
    get7->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get7->SetLineColor(kOrange);
    //get7->SetLineStyle(2);
    get7->SetMarkerStyle(3);
    get7->SetMarkerSize(0.7);
    get7->SetMarkerColorAlpha(kBlack,0.35);

    // theta = pi/4 with 200 MeV resolution
    TFile *f71 = new TFile("72pt_res_wedge_ff_pio4_200.root");
    TH1D *get71 = (TH1D *)f71->Get("FF");
    get71->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/4)");
    get71->GetYaxis()->SetTitle("");
    get71->Scale(197./get71->Integral(), "width");
    get71->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get71->SetLineColor(kOrange);
    //get71->SetLineStyle(2);
    get71->SetMarkerStyle(3);
    get71->SetMarkerSize(0.7);
    get71->SetMarkerColorAlpha(kBlack,0.35);

    // with saturation
    TFile *f8 = new TFile("phisat_pio4.root");
    TH1D *get8 = (TH1D *)f8->Get("FF");
    get8->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/4)");
    get8->GetYaxis()->SetTitle("");
    get8->Scale(197./get8->Integral(), "width");
    get8->GetYaxis()->SetRangeUser(1e-5,1e5);
    get8->SetLineColor(kRed);
    get8->SetLineStyle(1);
    //get8->SetMarkerStyle(4);
    //get8->SetMarkerSize(0.6);
    //get8->SetMarkerColorAlpha(kRed,0.35);

    // Saturation with 5 MeV resolution
    TFile *f22 = new TFile("phisat_pio4_wRes_5.root");
    TH1D *get22 = (TH1D *)f22->Get("FF");
    get22->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/4)");
    get22->GetYaxis()->SetTitle("");
    get22->Scale(197./get22->Integral(), "width");
    get22->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get22->SetLineColor(kOrange);
    //get22->SetLineStyle(2);
    get22->SetMarkerStyle(4);
    get22->SetMarkerSize(0.7);
    get22->SetMarkerColorAlpha(kRed,0.35);

    // Saturation with 10 MeV resolution
    TFile *f33 = new TFile("phisat_pio4_wRes_10.root");
    TH1D *get33 = (TH1D *)f33->Get("FF");
    get33->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/4)");
    get33->GetYaxis()->SetTitle("");
    get33->Scale(197./get33->Integral(), "width");
    get33->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get33->SetLineColor(kCyan);
    //get33->SetLineStyle(2);
    get33->SetMarkerStyle(4);
    get33->SetMarkerSize(0.7);
    get33->SetMarkerColorAlpha(kRed,0.35);

    // Saturation with 50 MeV resolution
    TFile *f44 = new TFile("phisat_pio4_wRes_50.root");
    TH1D *get44 = (TH1D *)f44->Get("FF");
    get44->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/4)");
    get44->GetYaxis()->SetTitle("");
    get44->Scale(197./get44->Integral(), "width");
    get44->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get44->SetLineColor(kGreen);
    //get44->SetLineStyle(2);
    get44->SetMarkerStyle(4);
    get44->SetMarkerSize(0.7);
    get44->SetMarkerColorAlpha(kRed,0.35);

    // Saturation with 100 MeV resolution
    TFile *f55 = new TFile("phisat_pio4_wRes_100.root");
    TH1D *get55 = (TH1D *)f55->Get("FF");
    get55->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/4)");
    get55->GetYaxis()->SetTitle("");
    get55->Scale(197./get55->Integral(), "width");
    get55->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get55->SetLineColor(kBlue);
    //get55->SetLineStyle(2);
    get55->SetMarkerStyle(4);
    get55->SetMarkerSize(0.7);
    get55->SetMarkerColorAlpha(kRed,0.35);

    // Saturation with 150 MeV resolution
    TFile *f66 = new TFile("phisat_pio4_wRes_150.root");
    TH1D *get66 = (TH1D *)f66->Get("FF");
    get66->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/4)");
    get66->GetYaxis()->SetTitle("");
    get66->Scale(197./get66->Integral(), "width");
    get66->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get66->SetLineColor(kViolet);
    //get66->SetLineStyle(2);
    get66->SetMarkerStyle(4);
    get66->SetMarkerSize(0.7);
    get66->SetMarkerColorAlpha(kRed,0.35);

    // Saturation with 200 MeV resolution
    TFile *f61 = new TFile("phisat_pio4_wRes_200.root");
    TH1D *get61 = (TH1D *)f61->Get("FF");
    get61->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/4)");
    get61->GetYaxis()->SetTitle("");
    get61->Scale(197./get61->Integral(), "width");
    get61->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get61->SetLineColor(kViolet);
    //get61->SetLineStyle(2);
    get61->SetMarkerStyle(4);
    get61->SetMarkerSize(0.7);
    get61->SetMarkerColorAlpha(kRed,0.35);

    TCanvas *c_f1 = new TCanvas();
    get1->Draw();
    //get3->Draw("P HIST same");
    //get4->Draw("P HIST same");
    //get5->Draw("P HIST same");
    //get6->Draw("P HIST same");
    get7->Draw("P HIST same");
    //get71->Draw("P HIST same");
    get8->Draw("same");
    get22->Draw("P HIST same");
    //get33->Draw("P HIST same");
    //get44->Draw("P HIST same");
    //get55->Draw("P HIST same");
    //get66->Draw("P HIST same");
    //get61->Draw("P HIST same");
    gPad->SetLogy(1);
    gStyle->SetOptStat(0);

    auto legend = new TLegend(0.43,0.71,0.73,0.9);
    legend->AddEntry(get8,"d#sigma/dt [nb/GeV^{2}] with sat: no res","l");
    legend->AddEntry(get22,"d#sigma/dt [nb/GeV^{2}] with sat: 5 MeV res","p");
    //legend->AddEntry(get33,"d#sigma/dt [nb/GeV^{2}] with sat: 10 MeV res","p");
    //legend->AddEntry(get44,"d#sigma/dt [nb/GeV^{2}] with sat: 50 MeV res","p");
    //legend->AddEntry(get55,"d#sigma/dt [nb/GeV^{2}] with sat: 100 MeV res","p");
    //legend->AddEntry(get66,"d#sigma/dt [nb/GeV^{2}] with sat: 150 MeV res","p");
    //legend->AddEntry(get61,"d#sigma/dt [nb/GeV^{2}] with sat: 200 MeV res","p");
    legend->AddEntry(get1,"|F(t)|^{2} (no resolution)","l");
    legend->AddEntry(get7,"|F(t)|^{2} with 5 MeV resolution","p");
    //legend->AddEntry(get5,"|F(t)|^{2} with 10 MeV resolution","p");
    //legend->AddEntry(get4,"|F(t)|^{2} with 50 MeV resolution","p");
    //legend->AddEntry(get3,"|F(t)|^{2} with 100 MeV resolution","p");
    //legend->AddEntry(get6,"|F(t)|^{2} with 150 MeV resolution","p");
    //legend->AddEntry(get71,"|F(t)|^{2} with 200 MeV resolution","p");
    legend->Draw();
}

void phipio3_Wsat()
{
    // theta = pi/3
    TFile *f1 = new TFile("72pt_wedge_ff_pio3.root");
    TH1D *get1 = (TH1D *)f1->Get("FF");
    get1->Scale(197./get1->Integral(), "width");
    get1->SetTitle("Wedge Integration and Saturartion Model (#phi) with Resolution (#theta_{max} = #pi/3)");
    get1->GetYaxis()->SetTitle("");
    get1->GetXaxis()->SetTitle("t [GeV^{2}]");
    get1->Scale(197./get1->Integral(), "width");
    get1->GetYaxis()->SetRangeUser(1e-5,1e5);
    get1->SetLineColor(kBlack);

    // theta = pi/3 with 100 MeV resolution
    TFile *f3 = new TFile("72pt_res_wedge_ff_pio3_100.root");
    TH1D *get3 = (TH1D *)f3->Get("FF");
    get3->SetTitle("Wedge Integration and Saturartion Model (#phi) with Resolution (#theta_{max} = #pi/3)");
    get3->GetYaxis()->SetTitle("");
    get3->Scale(197./get3->Integral(), "width");
    get3->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get3->SetLineColor(kBlue);
    //get3->SetLineStyle(2);
    get3->SetMarkerStyle(3);
    get3->SetMarkerSize(0.7);
    get3->SetMarkerColorAlpha(kBlack,0.35);

    // theta = pi/3 with 50 MeV resolution
    TFile *f4 = new TFile("72pt_res_wedge_ff_pio3_50.root");
    TH1D *get4 = (TH1D *)f4->Get("FF");
    get4->SetTitle("Wedge Integration and Saturartion Model (#phi) with Resolution (#theta_{max} = #pi/3)");
    get4->GetYaxis()->SetTitle("");
    get4->Scale(197./get4->Integral(), "width");
    get4->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get4->SetLineColor(kGreen);
    //get4->SetLineStyle(2);
    get4->SetMarkerStyle(3);
    get4->SetMarkerSize(0.7);
    get4->SetMarkerColorAlpha(kBlack,0.35);

    // theta = pi/3 with 10 MeV resolution
    TFile *f5 = new TFile("72pt_res_wedge_ff_pio3_10.root");
    TH1D *get5 = (TH1D *)f5->Get("FF");
    get5->SetTitle("Wedge Integration and Saturartion Model (#phi) with Resolution (#theta_{max} = #pi/3)");
    get5->GetYaxis()->SetTitle("");
    get5->Scale(197./get5->Integral(), "width");
    get5->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get5->SetLineColor(kCyan);
    //get5->SetLineStyle(2);
    get5->SetMarkerStyle(3);
    get5->SetMarkerSize(0.7);
    get5->SetMarkerColorAlpha(kBlack,0.35);

    // theta = pi/3 with 150 MeV resolution
    TFile *f6 = new TFile("72pt_res_wedge_ff_pio3_150.root");
    TH1D *get6 = (TH1D *)f6->Get("FF");
    get6->SetTitle("Wedge Integration and Saturartion Model (#phi) with Resolution (#theta_{max} = #pi/3)");
    get6->GetYaxis()->SetTitle("");
    get6->Scale(197./get6->Integral(), "width");
    get6->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get6->SetLineColor(kViolet);
    get6->SetMarkerStyle(3);
    get6->SetMarkerSize(0.7);
    get6->SetMarkerColorAlpha(kBlack,0.35);

    // theta = pi/3 with 5 MeV resolution
    TFile *f7 = new TFile("72pt_res_wedge_ff_pio3_5.root");
    TH1D *get7 = (TH1D *)f7->Get("FF");
    get7->SetTitle("Wedge Integration and Saturartion Model (#phi) with Resolution (#theta_{max} = #pi/3)");
    get7->GetYaxis()->SetTitle("");
    get7->Scale(197./get7->Integral(), "width");
    get7->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get7->SetLineColor(kOrange);
    //get7->SetLineStyle(2);
    get7->SetMarkerStyle(3);
    get7->SetMarkerSize(0.7);
    get7->SetMarkerColorAlpha(kBlack,0.35);

    // theta = pi/3 with 200 MeV resolution
    TFile *f71 = new TFile("72pt_res_wedge_ff_pio3_200.root");
    TH1D *get71 = (TH1D *)f71->Get("FF");
    get71->SetTitle("Wedge Integration and Saturartion Model (#phi) with Resolution (#theta_{max} = #pi/3)");
    get71->GetYaxis()->SetTitle("");
    get71->Scale(197./get71->Integral(), "width");
    get71->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get71->SetLineColor(kOrange);
    //get71->SetLineStyle(2);
    get71->SetMarkerStyle(3);
    get71->SetMarkerSize(0.7);
    get71->SetMarkerColorAlpha(kBlack,0.35);

    // with saturation
    TFile *f8 = new TFile("phisat_pio3.root");
    TH1D *get8 = (TH1D *)f8->Get("FF");
    get8->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/3)");
    get8->GetYaxis()->SetTitle("");
    get8->Scale(197./get8->Integral(), "width");
    get8->GetYaxis()->SetRangeUser(1e-5,1e5);
    get8->SetLineColor(kRed);
    get8->SetLineStyle(1);
    //get8->SetMarkerStyle(4);
    //get8->SetMarkerSize(0.6);
    //get8->SetMarkerColorAlpha(kRed,0.35);

    // Saturation with 5 MeV resolution
    TFile *f22 = new TFile("phisat_pio3_wRes_5.root");
    TH1D *get22 = (TH1D *)f22->Get("FF");
    get22->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/3)");
    get22->GetYaxis()->SetTitle("");
    get22->Scale(197./get22->Integral(), "width");
    get22->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get22->SetLineColor(kOrange);
    //get22->SetLineStyle(2);
    get22->SetMarkerStyle(4);
    get22->SetMarkerSize(0.7);
    get22->SetMarkerColorAlpha(kRed,0.35);

    // Saturation with 10 MeV resolution
    TFile *f33 = new TFile("phisat_pio3_wRes_10.root");
    TH1D *get33 = (TH1D *)f33->Get("FF");
    get33->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/3)");
    get33->GetYaxis()->SetTitle("");
    get33->Scale(197./get33->Integral(), "width");
    get33->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get33->SetLineColor(kCyan);
    //get33->SetLineStyle(2);
    get33->SetMarkerStyle(4);
    get33->SetMarkerSize(0.7);
    get33->SetMarkerColorAlpha(kRed,0.35);

    // Saturation with 50 MeV resolution
    TFile *f44 = new TFile("phisat_pio3_wRes_50.root");
    TH1D *get44 = (TH1D *)f44->Get("FF");
    get44->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/3)");
    get44->GetYaxis()->SetTitle("");
    get44->Scale(197./get44->Integral(), "width");
    get44->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get44->SetLineColor(kGreen);
    //get44->SetLineStyle(2);
    get44->SetMarkerStyle(4);
    get44->SetMarkerSize(0.7);
    get44->SetMarkerColorAlpha(kRed,0.35);

    // Saturation with 100 MeV resolution
    TFile *f55 = new TFile("phisat_pio3_wRes_100.root");
    TH1D *get55 = (TH1D *)f55->Get("FF");
    get55->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/3)");
    get55->GetYaxis()->SetTitle("");
    get55->Scale(197./get55->Integral(), "width");
    get55->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get55->SetLineColor(kBlue);
    //get55->SetLineStyle(2);
    get55->SetMarkerStyle(4);
    get55->SetMarkerSize(0.7);
    get55->SetMarkerColorAlpha(kRed,0.35);

    // Saturation with 150 MeV resolution
    TFile *f66 = new TFile("phisat_pio3_wRes_150.root");
    TH1D *get66 = (TH1D *)f66->Get("FF");
    get66->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/3)");
    get66->GetYaxis()->SetTitle("");
    get66->Scale(197./get66->Integral(), "width");
    get66->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get66->SetLineColor(kViolet);
    //get66->SetLineStyle(2);
    get66->SetMarkerStyle(4);
    get66->SetMarkerSize(0.7);
    get66->SetMarkerColorAlpha(kRed,0.35);

    // Saturation with 200 MeV resolution
    TFile *f61 = new TFile("phisat_pio3_wRes_200.root");
    TH1D *get61 = (TH1D *)f61->Get("FF");
    get61->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/3)");
    get61->GetYaxis()->SetTitle("");
    get61->Scale(197./get61->Integral(), "width");
    get61->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get61->SetLineColor(kViolet);
    //get61->SetLineStyle(2);
    get61->SetMarkerStyle(4);
    get61->SetMarkerSize(0.7);
    get61->SetMarkerColorAlpha(kRed,0.35);

    TCanvas *c_f1 = new TCanvas();
    //get1->Draw();
    //get3->Draw("P HIST same");
    //get4->Draw("P HIST same");
    //get5->Draw("P HIST same");
    //get6->Draw("P HIST same");
    get7->Draw("P HIST same");
    //get71->Draw("P HIST same");
    //get8->Draw("same");
    get22->Draw("P HIST same");
    //get33->Draw("P HIST same");
    //get44->Draw("P HIST same");
    //get55->Draw("P HIST same");
    //get66->Draw("P HIST same");
    //get61->Draw("P HIST same");
    gPad->SetLogy(1);
    gStyle->SetOptStat(0);

    auto legend = new TLegend(0.43,0.71,0.73,0.9);
    //legend->AddEntry(get8,"d#sigma/dt [nb/GeV^{2}] with sat: no res","l");
    legend->AddEntry(get22,"d#sigma/dt [nb/GeV^{2}] with sat: 5 MeV res","p");
    //legend->AddEntry(get33,"d#sigma/dt [nb/GeV^{2}] with sat: 10 MeV res","p");
    //legend->AddEntry(get44,"d#sigma/dt [nb/GeV^{2}] with sat: 50 MeV res","p");
    //legend->AddEntry(get55,"d#sigma/dt [nb/GeV^{2}] with sat: 100 MeV res","p");
    //legend->AddEntry(get66,"d#sigma/dt [nb/GeV^{2}] with sat: 150 MeV res","p");
    //legend->AddEntry(get61,"d#sigma/dt [nb/GeV^{2}] with sat: 200 MeV res","p");
    //legend->AddEntry(get1,"|F(t)|^{2} (no resolution)","l");
    legend->AddEntry(get7,"|F(t)|^{2} with 5 MeV resolution","p");
    //legend->AddEntry(get5,"|F(t)|^{2} with 10 MeV resolution","p");
    //legend->AddEntry(get4,"|F(t)|^{2} with 50 MeV resolution","p");
    //legend->AddEntry(get3,"|F(t)|^{2} with 100 MeV resolution","p");
    //legend->AddEntry(get6,"|F(t)|^{2} with 150 MeV resolution","p");
    //legend->AddEntry(get71,"|F(t)|^{2} with 200 MeV resolution","p");
    legend->Draw();
}

void phipio6_Wsat()
{
    // theta = pi/6
    TFile *f1 = new TFile("72pt_wedge_ff_pio6.root");
    TH1D *get1 = (TH1D *)f1->Get("FF");
    get1->Scale(197./get1->Integral(), "width");
    get1->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/6)");
    get1->GetYaxis()->SetTitle("");
    get1->GetXaxis()->SetTitle("t [GeV^{2}]");
    get1->Scale(197./get1->Integral(), "width");
    get1->GetYaxis()->SetRangeUser(1e-5,1e5);
    get1->SetLineColor(kBlack);

    // theta = pi/6 with 100 MeV resolution
    TFile *f3 = new TFile("72pt_res_wedge_ff_pio6_100.root");
    TH1D *get3 = (TH1D *)f3->Get("FF");
    get3->Scale(197./get3->Integral(), "width");
    get3->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/6)");
    get3->GetYaxis()->SetTitle("");
    get3->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get3->SetLineColor(kBlue);
    //get3->SetLineStyle(2);
    get3->SetMarkerStyle(3);
    get3->SetMarkerSize(0.7);
    get3->SetMarkerColorAlpha(kBlack,0.35);

    // theta = pi/6 with 50 MeV resolution
    TFile *f4 = new TFile("72pt_res_wedge_ff_pio6_50.root");
    TH1D *get4 = (TH1D *)f4->Get("FF");
    get4->Scale(197./get4->Integral(), "width");
    get4->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/6)");
    get4->GetYaxis()->SetTitle("");
    get4->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get4->SetLineColor(kGreen);
    //get4->SetLineStyle(2);
    get4->SetMarkerStyle(3);
    get4->SetMarkerSize(0.7);
    get4->SetMarkerColorAlpha(kBlack,0.35);

    // theta = pi/6 with 10 MeV resolution
    TFile *f5 = new TFile("72pt_res_wedge_ff_pio6_10.root");
    TH1D *get5 = (TH1D *)f5->Get("FF");
    get5->Scale(197./get5->Integral(), "width");
    get5->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/6)");
    get5->GetYaxis()->SetTitle("");
    get5->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get5->SetLineColor(kCyan);
    //get5->SetLineStyle(2);
    get5->SetMarkerStyle(3);
    get5->SetMarkerSize(0.7);
    get5->SetMarkerColorAlpha(kBlack,0.35);

    // theta = pi/6 with 150 MeV resolution
    TFile *f6 = new TFile("72pt_res_wedge_ff_pio6_150.root");
    TH1D *get6 = (TH1D *)f6->Get("FF");
    get6->Scale(197./get6->Integral(), "width");
    get6->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/6)");
    get6->GetYaxis()->SetTitle("");
    get6->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get6->SetLineColor(kViolet);
    get6->SetMarkerStyle(3);
    get6->SetMarkerSize(0.7);
    get6->SetMarkerColorAlpha(kBlack,0.35);

    // theta = pi/6 with 5 MeV resolution
    TFile *f7 = new TFile("72pt_res_wedge_ff_pio6_5.root");
    TH1D *get7 = (TH1D *)f7->Get("FF");
    get7->Scale(197./get7->Integral(), "width");
    get7->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/6)");
    get7->GetYaxis()->SetTitle("");
    get7->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get7->SetLineColor(kOrange);
    //get7->SetLineStyle(2);
    get7->SetMarkerStyle(3);
    get7->SetMarkerSize(0.7);
    get7->SetMarkerColorAlpha(kBlack,0.35);

    // theta = pi/6 with 200 MeV resolution
    TFile *f71 = new TFile("72pt_res_wedge_ff_pio6_200.root");
    TH1D *get71 = (TH1D *)f71->Get("FF");
    get71->Scale(197./get71->Integral(), "width");
    get71->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/6)");
    get71->GetYaxis()->SetTitle("");
    get71->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get71->SetLineColor(kOrange);
    //get71->SetLineStyle(2);
    get71->SetMarkerStyle(3);
    get71->SetMarkerSize(0.7);
    get71->SetMarkerColorAlpha(kBlack,0.35);

    // with saturation
    TFile *f8 = new TFile("phisat_pio6.root");
    TH1D *get8 = (TH1D *)f8->Get("FF");
    get8->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/6)");
    get8->GetYaxis()->SetTitle("");
    get8->Scale(197./get8->Integral(), "width");
    get8->GetYaxis()->SetRangeUser(1e-5,1e5);
    get8->SetLineColor(kRed);
    get8->SetLineStyle(1);
    //get8->SetMarkerStyle(4);
    //get8->SetMarkerSize(0.6);
    //get8->SetMarkerColorAlpha(kRed,0.35);

    // Saturation with 5 MeV resolution
    TFile *f22 = new TFile("phisat_pio6_wRes_5.root");
    TH1D *get22 = (TH1D *)f22->Get("FF");
    get22->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/6)");
    get22->GetYaxis()->SetTitle("");
    get22->Scale(197./get22->Integral(), "width");
    get22->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get22->SetLineColor(kOrange);
    //get22->SetLineStyle(2);
    get22->SetMarkerStyle(4);
    get22->SetMarkerSize(0.7);
    get22->SetMarkerColorAlpha(kRed,0.35);

    // Saturation with 10 MeV resolution
    TFile *f33 = new TFile("phisat_pio6_wRes_10.root");
    TH1D *get33 = (TH1D *)f33->Get("FF");
    get33->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/6)");
    get33->GetYaxis()->SetTitle("");
    get33->Scale(197./get33->Integral(), "width");
    get33->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get33->SetLineColor(kCyan);
    //get33->SetLineStyle(2);
    get33->SetMarkerStyle(4);
    get33->SetMarkerSize(0.7);
    get33->SetMarkerColorAlpha(kRed,0.35);

    // Saturation with 50 MeV resolution
    TFile *f44 = new TFile("phisat_pio6_wRes_50.root");
    TH1D *get44 = (TH1D *)f44->Get("FF");
    get44->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/6)");
    get44->GetYaxis()->SetTitle("");
    get44->Scale(197./get44->Integral(), "width");
    get44->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get44->SetLineColor(kGreen);
    //get44->SetLineStyle(2);
    get44->SetMarkerStyle(4);
    get44->SetMarkerSize(0.7);
    get44->SetMarkerColorAlpha(kRed,0.35);

    // Saturation with 100 MeV resolution
    TFile *f55 = new TFile("phisat_pio6_wRes_100.root");
    TH1D *get55 = (TH1D *)f55->Get("FF");
    get55->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/6)");
    get55->GetYaxis()->SetTitle("");
    get55->Scale(197./get55->Integral(), "width");
    get55->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get55->SetLineColor(kBlue);
    //get55->SetLineStyle(2);
    get55->SetMarkerStyle(4);
    get55->SetMarkerSize(0.7);
    get55->SetMarkerColorAlpha(kRed,0.35);

    // Saturation with 150 MeV resolution
    TFile *f66 = new TFile("phisat_pio6_wRes_150.root");
    TH1D *get66 = (TH1D *)f66->Get("FF");
    get66->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/6)");
    get66->GetYaxis()->SetTitle("");
    get66->Scale(197./get66->Integral(), "width");
    get66->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get66->SetLineColor(kViolet);
    //get66->SetLineStyle(2);
    get66->SetMarkerStyle(4);
    get66->SetMarkerSize(0.7);
    get66->SetMarkerColorAlpha(kRed,0.35);

    // Saturation with 200 MeV resolution
    TFile *f61 = new TFile("phisat_pio6_wRes_200.root");
    TH1D *get61 = (TH1D *)f61->Get("FF");
    get61->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/6)");
    get61->GetYaxis()->SetTitle("");
    get61->Scale(197./get61->Integral(), "width");
    get61->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get61->SetLineColor(kViolet);
    //get61->SetLineStyle(2);
    get61->SetMarkerStyle(4);
    get61->SetMarkerSize(0.7);
    get61->SetMarkerColorAlpha(kRed,0.35);

    TCanvas *c_f1 = new TCanvas();
    get1->Draw();
    //get3->Draw("P HIST same");
    //get4->Draw("P HIST same");
    //get5->Draw("P HIST same");
    //get6->Draw("P HIST same");
    get7->Draw("P HIST same");
    //get71->Draw("P HIST same");
    get8->Draw("same");
    get22->Draw("P HIST same");
    //get33->Draw("P HIST same");
    //get44->Draw("P HIST same");
    //get55->Draw("P HIST same");
    //get66->Draw("P HIST same");
    //get61->Draw("P HIST same");
    gPad->SetLogy(1);
    gStyle->SetOptStat(0);

    auto legend = new TLegend(0.43,0.71,0.73,0.9);
    legend->AddEntry(get8,"d#sigma/dt [nb/GeV^{2}] with sat: no res","l");
    legend->AddEntry(get22,"d#sigma/dt [nb/GeV^{2}] with sat: 5 MeV res","p");
    //legend->AddEntry(get33,"d#sigma/dt [nb/GeV^{2}] with sat: 10 MeV res","p");
    //legend->AddEntry(get44,"d#sigma/dt [nb/GeV^{2}] with sat: 50 MeV res","p");
    //legend->AddEntry(get55,"d#sigma/dt [nb/GeV^{2}] with sat: 100 MeV res","p");
    //legend->AddEntry(get66,"d#sigma/dt [nb/GeV^{2}] with sat: 150 MeV res","p");
    //legend->AddEntry(get61,"d#sigma/dt [nb/GeV^{2}] with sat: 200 MeV res","p");
    legend->AddEntry(get1,"|F(t)|^{2} (no resolutiuon)","l");
    legend->AddEntry(get7,"|F(t)|^{2} with 5 MeV resolution","p");
    //legend->AddEntry(get5,"|F(t)|^{2} with 10 MeV resolution","p");
    //legend->AddEntry(get4,"|F(t)|^{2} with 50 MeV resolution","p");
    //legend->AddEntry(get3,"|F(t)|^{2} with 100 MeV resolution","p");
    //legend->AddEntry(get6,"|F(t)|^{2} with 150 MeV resolution","p");
    //legend->AddEntry(get71,"|F(t)|^{2} with 200 MeV resolution","p");
    legend->Draw();
}

void phipio9_Wsat()
{
    // theta = pi/9
    TFile *f1 = new TFile("72pt_wedge_ff_pio9.root");
    TH1D *get1 = (TH1D *)f1->Get("FF");
    get1->Scale(197./get1->Integral(), "width");
    get1->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/9)");
    get1->GetYaxis()->SetTitle("");
    get1->GetXaxis()->SetTitle("t [GeV^{2}]");
    get1->Scale(197./get1->Integral(), "width");
    get1->GetYaxis()->SetRangeUser(1e-5,1e5);
    get1->SetLineColor(kBlack);

    // theta = pi/9 with 100 MeV resolution
    TFile *f3 = new TFile("72pt_res_wedge_ff_pio9_100.root");
    TH1D *get3 = (TH1D *)f3->Get("FF");
    get3->Scale(197./get3->Integral(), "width");
    get3->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/9)");
    get3->GetYaxis()->SetTitle("");
    //get3->SetLineColor(kBlue);
    get3->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get3->SetLineStyle(2);
    get3->SetMarkerStyle(3);
    get3->SetMarkerSize(0.7);
    get3->SetMarkerColorAlpha(kBlack,0.35);

    // theta = pi/9 with 50 MeV resolution
    TFile *f4 = new TFile("72pt_res_wedge_ff_pio9_50.root");
    TH1D *get4 = (TH1D *)f4->Get("FF");
    get4->Scale(197./get4->Integral(), "width");
    get4->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/9)");
    get4->GetYaxis()->SetTitle("");
    //get4->SetLineColor(kGreen);
    get4->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get4->SetLineStyle(2);
    get4->SetMarkerStyle(3);
    get4->SetMarkerSize(0.7);
    get4->SetMarkerColorAlpha(kBlack,0.35);

    // theta = pi/9 with 10 MeV resolution
    TFile *f5 = new TFile("72pt_res_wedge_ff_pio9_10.root");
    TH1D *get5 = (TH1D *)f5->Get("FF");
    get5->Scale(197./get5->Integral(), "width");
    get5->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/9)");
    get5->GetYaxis()->SetTitle("");
    //get5->SetLineColor(kCyan);
    get5->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get5->SetLineStyle(2);
    get5->SetMarkerStyle(3);
    get5->SetMarkerSize(0.7);
    get5->SetMarkerColorAlpha(kBlack,0.35);

    // theta = pi/9 with 150 MeV resolution
    TFile *f6 = new TFile("72pt_res_wedge_ff_pio9_150.root");
    TH1D *get6 = (TH1D *)f6->Get("FF");
    get6->Scale(197./get6->Integral(), "width");
    get6->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/9)");
    get6->GetYaxis()->SetTitle("");
    get6->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get6->SetLineColor(kViolet);
    get6->SetMarkerStyle(3);
    get6->SetMarkerSize(0.7);
    get6->SetMarkerColorAlpha(kBlack,0.35);

    // theta = pi/9 with 5 MeV resolution
    TFile *f7 = new TFile("72pt_res_wedge_ff_pio9_5.root");
    TH1D *get7 = (TH1D *)f7->Get("FF");
    get7->Scale(197./get7->Integral(), "width");
    get7->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/9)");
    get7->GetYaxis()->SetTitle("");
    get7->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get7->SetLineColor(kOrange);
    //get7->SetLineStyle(2);
    get7->SetMarkerStyle(3);
    get7->SetMarkerSize(0.7);
    get7->SetMarkerColorAlpha(kBlack,0.35);

    // theta = pi/9 with 200 MeV resolution
    TFile *f71 = new TFile("72pt_res_wedge_ff_pio9_200.root");
    TH1D *get71 = (TH1D *)f71->Get("FF");
    get71->Scale(197./get71->Integral(), "width");
    get71->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/9)");
    get71->GetYaxis()->SetTitle("");
    get71->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get71->SetLineColor(kOrange);
    //get71->SetLineStyle(2);
    get71->SetMarkerStyle(3);
    get71->SetMarkerSize(0.7);
    get71->SetMarkerColorAlpha(kBlack,0.35);

    // with saturation
    TFile *f8 = new TFile("phisat_pio9.root");
    TH1D *get8 = (TH1D *)f8->Get("FF");
    get8->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/9)");
    get8->GetYaxis()->SetTitle("");
    get8->Scale(197./get8->Integral(), "width");
    get8->GetYaxis()->SetRangeUser(1e-5,1e5);
    get8->SetLineColor(kRed);
    get8->SetLineStyle(1);
    //get8->SetMarkerStyle(4);
    //get8->SetMarkerSize(0.6);
    //get8->SetMarkerColorAlpha(kRed,0.35);

    // Saturation with 5 MeV resolution
    TFile *f22 = new TFile("phisat_pio9_wRes_5.root");
    TH1D *get22 = (TH1D *)f22->Get("FF");
    get22->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/9)");
    get22->GetYaxis()->SetTitle("");
    get22->Scale(197./get22->Integral(), "width");
    get22->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get22->SetLineColor(kOrange);
    //get22->SetLineStyle(2);
    get22->SetMarkerStyle(4);
    get22->SetMarkerSize(0.7);
    get22->SetMarkerColorAlpha(kRed,0.35);

    // Saturation with 10 MeV resolution
    TFile *f33 = new TFile("phisat_pio9_wRes_10.root");
    TH1D *get33 = (TH1D *)f33->Get("FF");
    get33->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/9)");
    get33->GetYaxis()->SetTitle("");
    get33->Scale(197./get33->Integral(), "width");
    get33->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get33->SetLineColor(kCyan);
    //get33->SetLineStyle(2);
    get33->SetMarkerStyle(4);
    get33->SetMarkerSize(0.7);
    get33->SetMarkerColorAlpha(kRed,0.35);

    // Saturation with 50 MeV resolution
    TFile *f44 = new TFile("phisat_pio9_wRes_10.root");
    TH1D *get44 = (TH1D *)f44->Get("FF");
    get44->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/9)");
    get44->GetYaxis()->SetTitle("");
    get44->Scale(197./get44->Integral(), "width");
    get44->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get44->SetLineColor(kGreen);
    //get44->SetLineStyle(2);
    get44->SetMarkerStyle(4);
    get44->SetMarkerSize(0.7);
    get44->SetMarkerColorAlpha(kRed,0.35);

    // Saturation with 100 MeV resolution
    TFile *f55 = new TFile("phisat_pio9_wRes_100.root");
    TH1D *get55 = (TH1D *)f55->Get("FF");
    get55->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/9)");
    get55->GetYaxis()->SetTitle("");
    get55->Scale(197./get55->Integral(), "width");
    get55->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get55->SetLineColor(kBlue);
    //get55->SetLineStyle(2);
    get55->SetMarkerStyle(4);
    get55->SetMarkerSize(0.7);
    get55->SetMarkerColorAlpha(kRed,0.35);

    // Saturation with 150 MeV resolution
    TFile *f66 = new TFile("phisat_pio9_wRes_150.root");
    TH1D *get66 = (TH1D *)f66->Get("FF");
    get66->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/9)");
    get66->GetYaxis()->SetTitle("");
    get66->Scale(197./get66->Integral(), "width");
    get66->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get66->SetLineColor(kViolet);
    //get66->SetLineStyle(2);
    get66->SetMarkerStyle(4);
    get66->SetMarkerSize(0.7);
    get66->SetMarkerColorAlpha(kRed,0.35);

    // Saturation with 200 MeV resolution
    TFile *f61 = new TFile("phisat_pio9_wRes_200.root");
    TH1D *get61 = (TH1D *)f61->Get("FF");
    get61->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi/9)");
    get61->GetYaxis()->SetTitle("");
    get61->Scale(197./get61->Integral(), "width");
    get61->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get61->SetLineColor(kViolet);
    //get61->SetLineStyle(2);
    get61->SetMarkerStyle(4);
    get61->SetMarkerSize(0.7);
    get61->SetMarkerColorAlpha(kRed,0.35);

    TCanvas *c_f1 = new TCanvas();
    //get1->Draw();
    //get3->Draw("P HIST same");
    //get4->Draw("P HIST same");
    //get5->Draw("P HIST same");
    //get6->Draw("P HIST same");
    get7->Draw("P HIST same");
    //get71->Draw("P HIST same");
    //get8->Draw("same");
    get22->Draw("P HIST same");
    //get33->Draw("P HIST same");
    //get44->Draw("P HIST same");
    //get55->Draw("P HIST same");
    //get66->Draw("P HIST same");
    //get61->Draw("P HIST same");
    gPad->SetLogy(1);
    gStyle->SetOptStat(0);

    auto legend = new TLegend(0.43,0.71,0.73,0.9);
    //legend->AddEntry(get8,"d#sigma/dt [nb/GeV^{2}] with sat: no res","l");
    legend->AddEntry(get22,"d#sigma/dt [nb/GeV^{2}] with sat: 5 MeV res","p");
    //legend->AddEntry(get33,"d#sigma/dt [nb/GeV^{2}] with sat: 10 MeV res","p");
    //legend->AddEntry(get44,"d#sigma/dt [nb/GeV^{2}] with sat: 50 MeV res","p");
    //legend->AddEntry(get55,"d#sigma/dt [nb/GeV^{2}] with sat: 100 MeV res","p");
    //legend->AddEntry(get66,"d#sigma/dt [nb/GeV^{2}] with sat: 150 MeV res","p");
    //legend->AddEntry(get61,"d#sigma/dt [nb/GeV^{2}] with sat: 200 MeV res","p");
    //legend->AddEntry(get1,"|F(t)|^{2} (no resolution)","l");
    legend->AddEntry(get7,"|F(t)|^{2} with 5 MeV resolution","p");
    //legend->AddEntry(get5,"|F(t)|^{2} with 10 MeV resolution","p");
    //legend->AddEntry(get4,"|F(t)|^{2} with 50 MeV resolution","p");
    //legend->AddEntry(get3,"|F(t)|^{2} with 100 MeV resolution","p");
    //legend->AddEntry(get6,"|F(t)|^{2} with 150 MeV resolution","p");
    //legend->AddEntry(get71,"|F(t)|^{2} with 200 MeV resolution","p");
    legend->Draw();
}

void phizero_Wsat()
{
    // theta = 0  (1e-12)
    TFile *f1 = new TFile("72pt_wedge_ff_0.root");
    TH1D *get1 = (TH1D *)f1->Get("FF");
    get1->Scale(197./get1->Integral(), "width");
    get1->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = 0)");
    get1->GetYaxis()->SetTitle("");
    get1->GetYaxis()->SetRangeUser(1e-5,1e5);
    get1->GetXaxis()->SetTitle("t [GeV^{2}]");
    get1->Scale(197./get1->Integral(), "width");
    get1->SetLineColor(kBlack);

    // theta = 0 with 100 MeV resolution
    TFile *f3 = new TFile("72pt_res_wedge_ff_0_100.root");
    TH1D *get3 = (TH1D *)f3->Get("FF");
    get3->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = 0)");
    get3->GetYaxis()->SetTitle("");
    get3->GetYaxis()->SetRangeUser(1e-5,1e5);
    get3->Scale(197./get3->Integral(), "width");
    //get3->SetLineColor(kBlue);
    //get3->SetLineStyle(2);
    get3->SetMarkerStyle(3);
    get3->SetMarkerSize(0.7);
    get3->SetMarkerColorAlpha(kBlack,0.35);

    // theta = 0 with 50 MeV resolution
    TFile *f4 = new TFile("72pt_res_wedge_ff_0_50.root");
    TH1D *get4 = (TH1D *)f4->Get("FF");
    get4->Scale(197./get4->Integral(), "width");
    get4->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = 0)");
    get4->GetYaxis()->SetTitle("");
    get4->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get4->SetLineColor(kGreen);
    //get4->SetLineStyle(2);
    get4->SetMarkerStyle(3);
    get4->SetMarkerSize(0.7);
    get4->SetMarkerColorAlpha(kBlack,0.35);

    // theta = 0 with 10 MeV resolution
    TFile *f5 = new TFile("72pt_res_wedge_ff_0_10.root");
    TH1D *get5 = (TH1D *)f5->Get("FF");
    get5->Scale(197./get5->Integral(), "width");
    get5->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = 0)");
    get5->GetYaxis()->SetTitle("");
    get5->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get5->SetLineColor(kCyan);
    //get5->SetLineStyle(2);
    get5->SetMarkerStyle(3);
    get5->SetMarkerSize(0.7);
    get5->SetMarkerColorAlpha(kBlack,0.35);

    // theta = 0 with 150 MeV resolution
    TFile *f6 = new TFile("72pt_res_wedge_ff_0_150.root");
    TH1D *get6 = (TH1D *)f6->Get("FF");
    get6->Scale(197./get6->Integral(), "width");
    get6->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = 0)");
    get6->GetYaxis()->SetTitle("");
    get6->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get6->SetLineColor(kViolet);
    get6->SetMarkerStyle(3);
    get6->SetMarkerSize(0.7);
    get6->SetMarkerColorAlpha(kBlack,0.35);

    // theta = 0 with 5 MeV resolution
    TFile *f7 = new TFile("72pt_res_wedge_ff_0_5.root");
    TH1D *get7 = (TH1D *)f7->Get("FF");
    get7->Scale(197./get7->Integral(), "width");
    get7->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = 0)");
    get7->GetYaxis()->SetTitle("");
    get7->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get7->SetLineColor(kOrange);
    //get7->SetLineStyle(2);
    get7->SetMarkerStyle(3);
    get7->SetMarkerSize(0.7);
    get7->SetMarkerColorAlpha(kBlack,0.35);

    // theta = 0 with 200 MeV resolution
    TFile *f71 = new TFile("72pt_res_wedge_ff_0_200.root");
    TH1D *get71 = (TH1D *)f71->Get("FF");
    get71->Scale(197./get71->Integral(), "width");
    get71->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = 0)");
    get71->GetYaxis()->SetTitle("");
    get71->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get71->SetLineColor(kOrange);
    //get71->SetLineStyle(2);
    get71->SetMarkerStyle(3);
    get71->SetMarkerSize(0.7);
    get71->SetMarkerColorAlpha(kBlack,0.35);

    // with saturation
    TFile *f8 = new TFile("phisat_0.root");
    TH1D *get8 = (TH1D *)f8->Get("FF");
    get8->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = 0)");
    get8->GetYaxis()->SetTitle("");
    get8->Scale(197./get8->Integral(), "width");
    get8->GetYaxis()->SetRangeUser(1e-5,1e5);
    get8->SetLineColor(kRed);
    get8->SetLineStyle(1);
    //get8->SetMarkerStyle(4);
    //get8->SetMarkerSize(0.6);
    //get8->SetMarkerColorAlpha(kRed,0.35);

    // Saturation with 5 MeV resolution
    TFile *f22 = new TFile("phisat_0_wRes_5.root");
    TH1D *get22 = (TH1D *)f22->Get("FF");
    get22->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = 0)");
    get22->GetYaxis()->SetTitle("");
    get22->Scale(197./get22->Integral(), "width");
    get22->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get22->SetLineColor(kOrange);
    //get22->SetLineStyle(2);
    get22->SetMarkerStyle(4);
    get22->SetMarkerSize(0.7);
    get22->SetMarkerColorAlpha(kRed,0.35);

    // Saturation with 10 MeV resolution
    TFile *f33 = new TFile("phisat_0_wRes_10.root");
    TH1D *get33 = (TH1D *)f33->Get("FF");
    get33->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = 0)");
    get33->GetYaxis()->SetTitle("");
    get33->Scale(197./get33->Integral(), "width");
    get33->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get33->SetLineColor(kCyan);
    //get33->SetLineStyle(2);
    get33->SetMarkerStyle(4);
    get33->SetMarkerSize(0.7);
    get33->SetMarkerColorAlpha(kRed,0.35);

    // Saturation with 50 MeV resolution
    TFile *f44 = new TFile("phisat_0_wRes_50.root");
    TH1D *get44 = (TH1D *)f44->Get("FF");
    get44->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = 0)");
    get44->GetYaxis()->SetTitle("");
    get44->Scale(197./get44->Integral(), "width");
    get44->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get44->SetLineColor(kGreen);
    //get44->SetLineStyle(2);
    get44->SetMarkerStyle(4);
    get44->SetMarkerSize(0.7);
    get44->SetMarkerColorAlpha(kRed,0.35);

    // Saturation with 100 MeV resolution
    TFile *f55 = new TFile("phisat_0_wRes_100.root");
    TH1D *get55 = (TH1D *)f55->Get("FF");
    get55->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = 0)");
    get55->GetYaxis()->SetTitle("");
    get55->Scale(197./get55->Integral(), "width");
    get55->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get55->SetLineColor(kBlue);
    //get55->SetLineStyle(2);
    get55->SetMarkerStyle(4);
    get55->SetMarkerSize(0.7);
    get55->SetMarkerColorAlpha(kRed,0.35);

    // Saturation with 150 MeV resolution
    TFile *f66 = new TFile("phisat_0_wRes_150.root");
    TH1D *get66 = (TH1D *)f66->Get("FF");
    get66->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = 0)");
    get66->GetYaxis()->SetTitle("");
    get66->Scale(197./get66->Integral(), "width");
    get66->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get66->SetLineColor(kViolet);
    //get66->SetLineStyle(2);
    get66->SetMarkerStyle(4);
    get66->SetMarkerSize(0.7);
    get66->SetMarkerColorAlpha(kRed,0.35);

    // Saturation with 200 MeV resolution
    TFile *f61 = new TFile("phisat_0_wRes_200.root");
    TH1D *get61 = (TH1D *)f61->Get("FF");
    get61->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = 0)");
    get61->GetYaxis()->SetTitle("");
    get61->Scale(197./get61->Integral(), "width");
    get61->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get61->SetLineColor(kViolet);
    //get61->SetLineStyle(2);
    get61->SetMarkerStyle(4);
    get61->SetMarkerSize(0.7);
    get61->SetMarkerColorAlpha(kRed,0.35);

    TCanvas *c_f1 = new TCanvas();
    //get1->Draw();
    //get3->Draw("P HIST same");
    //get4->Draw("P HIST same");
    get5->Draw("P HIST same");
    //get6->Draw("P HIST same");
    //get7->Draw("P HIST same");
    //get71->Draw("P HIST same");
    //get8->Draw("same");
    //get22->Draw("P HIST same");
    get33->Draw("P HIST same");
    //get44->Draw("P HIST same");
    //get55->Draw("P HIST same");
    //get66->Draw("P HIST same");
    //get61->Draw("P HIST same");
    gPad->SetLogy(1);
    gStyle->SetOptStat(0);

    auto legend = new TLegend(0.43,0.71,0.73,0.9);
    //legend->AddEntry(get8,"d#sigma/dt [nb/GeV^{2}] with sat: no res","l");
    legend->AddEntry(get22,"d#sigma/dt [nb/GeV^{2}] with sat: 5 MeV res","p");
    //legend->AddEntry(get33,"d#sigma/dt [nb/GeV^{2}] with sat: 10 MeV res","p");
    //legend->AddEntry(get44,"d#sigma/dt [nb/GeV^{2}] with sat: 50 MeV res","p");
    //legend->AddEntry(get55,"d#sigma/dt [nb/GeV^{2}] with sat: 100 MeV res","p");
    //legend->AddEntry(get66,"d#sigma/dt [nb/GeV^{2}] with sat: 150 MeV res","p");
    //legend->AddEntry(get61,"d#sigma/dt [nb/GeV^{2}] with sat: 200 MeV res","p");
    //legend->AddEntry(get1,"|F(t)|^{2} (no resolution)","l");
    legend->AddEntry(get7,"|F(t)|^{2} with 5 MeV resolution","p");
    //legend->AddEntry(get5,"|F(t)|^{2} with 10 MeV resolution","p");
    //legend->AddEntry(get4,"|F(t)|^{2} with 50 MeV resolution","p");
    //legend->AddEntry(get3,"|F(t)|^{2} with 100 MeV resolution","p");
    //legend->AddEntry(get6,"|F(t)|^{2} with 150 MeV resolution","p");
    //legend->AddEntry(get71,"|F(t)|^{2} with 200 MeV resolution","p");
    legend->Draw();
}

void phitwopi_Wsat()
{
    // theta = 2pi
    TFile *f1 = new TFile("72pt_wedge_ff_2pi.root");
    TH1D *get1 = (TH1D *)f1->Get("FF");
    get1->Scale(197./get1->Integral(), "width");
    get1->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = 2#pi)");
    get1->GetYaxis()->SetTitle("");
    get1->GetYaxis()->SetRangeUser(1e-5,1e5);
    get1->GetXaxis()->SetTitle("t [GeV^{2}]");
    get1->Scale(197./get1->Integral(), "width");
    //get1->GetYaxis()->SetRangeUser(1e-5,1e5);
    get1->SetLineColor(kBlack);

    // theta = 2pi with 100 MeV resolution
    TFile *f3 = new TFile("72pt_res_wedge_ff_2pi_100.root");
    TH1D *get3 = (TH1D *)f3->Get("FF");
    get3->Scale(197./get3->Integral(), "width");
    get3->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = 2#pi)");
    get3->GetYaxis()->SetTitle("");
    get3->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get3->SetLineColor(kBlue);
    //get3->SetLineStyle(2);
    get3->SetMarkerStyle(3);
    get3->SetMarkerSize(0.7);
    get3->SetMarkerColorAlpha(kBlack,0.35);

    // theta = 2pi with 10 MeV resolution
    TFile *f4 = new TFile("72pt_res_wedge_ff_2pi_10.root");
    TH1D *get4 = (TH1D *)f4->Get("FF");
    get4->Scale(197./get4->Integral(), "width");
    get4->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = 2#pi)");
    get4->GetYaxis()->SetTitle("");
    get4->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get4->SetLineColor(kGreen);
    //get4->SetLineStyle(2);
    get4->SetMarkerStyle(3);
    get4->SetMarkerSize(0.7);
    get4->SetMarkerColorAlpha(kBlack,0.35);

    // theta = 2pi with 50 MeV resolution
    TFile *f5 = new TFile("72pt_res_wedge_ff_2pi_50.root");
    TH1D *get5 = (TH1D *)f5->Get("FF");
    get5->Scale(197./get5->Integral(), "width");
    get5->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = 2#pi)");
    get5->GetYaxis()->SetTitle("");
    get5->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get5->SetLineColor(kCyan);
    //get5->SetLineStyle(2);
    get5->SetMarkerStyle(3);
    get5->SetMarkerSize(0.7);
    get5->SetMarkerColorAlpha(kBlack,0.35);

    // theta = 2pi with 150 MeV resolution
    TFile *f6 = new TFile("72pt_res_wedge_ff_2pi_150.root");
    TH1D *get6 = (TH1D *)f6->Get("FF");
    get6->Scale(197./get6->Integral(), "width");
    get6->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = 2#pi)");
    get6->GetYaxis()->SetTitle("");
    get6->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get6->SetLineColor(kViolet);
    get6->SetMarkerStyle(3);
    get6->SetMarkerSize(0.7);
    get6->SetMarkerColorAlpha(kBlack,0.35);

    // theta = 2pi with 5 MeV resolution
    TFile *f7 = new TFile("72pt_res_wedge_ff_2pi_5.root");
    TH1D *get7 = (TH1D *)f7->Get("FF");
    get7->Scale(197./get7->Integral(), "width");
    get7->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = 2#pi)");
    get7->GetYaxis()->SetTitle("");
    get7->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get7->SetLineColor(kOrange);
    //get7->SetLineStyle(2);
    get7->SetMarkerStyle(3);
    get7->SetMarkerSize(0.7);
    get7->SetMarkerColorAlpha(kBlack,0.35);

    // theta = 2pi with 200 MeV resolution
    TFile *f71 = new TFile("72pt_res_wedge_ff_2pi_200.root");
    TH1D *get71 = (TH1D *)f71->Get("FF");
    get71->Scale(197./get71->Integral(), "width");
    get71->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = 2#pi)");
    get71->GetYaxis()->SetTitle("");
    get71->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get71->SetLineColor(kOrange);
    //get71->SetLineStyle(2);
    get71->SetMarkerStyle(3);
    get71->SetMarkerSize(0.7);
    get71->SetMarkerColorAlpha(kBlack,0.35);

    // with saturation
    TFile *f8 = new TFile("phisat_2pi.root");
    TH1D *get8 = (TH1D *)f8->Get("FF");
    get8->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = 2#pi)");
    get8->GetYaxis()->SetTitle("");
    get8->Scale(197./get8->Integral(), "width");
    get8->GetYaxis()->SetRangeUser(1e-5,1e5);
    get8->SetLineColor(kRed);
    get8->SetLineStyle(1);
    //get8->SetMarkerStyle(4);
    //get8->SetMarkerSize(0.6);
    //get8->SetMarkerColorAlpha(kRed,0.35);

    // Saturation with 5 MeV resolution
    TFile *f22 = new TFile("phisat_2pi_wRes_5.root");
    TH1D *get22 = (TH1D *)f22->Get("FF");
    get22->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = 2#pi)");
    get22->GetYaxis()->SetTitle("");
    get22->Scale(197./get22->Integral(), "width");
    get22->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get22->SetLineColor(kOrange);
    //get22->SetLineStyle(2);
    get22->SetMarkerStyle(4);
    get22->SetMarkerSize(0.7);
    get22->SetMarkerColorAlpha(kRed,0.35);

    // Saturation with 10 MeV resolution
    TFile *f33 = new TFile("phisat_2pi_wRes_10.root");
    TH1D *get33 = (TH1D *)f33->Get("FF");
    get33->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = 2#pi)");
    get33->GetYaxis()->SetTitle("");
    get33->Scale(197./get33->Integral(), "width");
    get33->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get33->SetLineColor(kCyan);
    //get33->SetLineStyle(2);
    get33->SetMarkerStyle(4);
    get33->SetMarkerSize(0.7);
    get33->SetMarkerColorAlpha(kRed,0.35);

    // Saturation with 50 MeV resolution
    TFile *f44 = new TFile("phisat_2pi_wRes_50.root");
    TH1D *get44 = (TH1D *)f44->Get("FF");
    get44->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = 2#pi)");
    get44->GetYaxis()->SetTitle("");
    get44->Scale(197./get44->Integral(), "width");
    get44->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get44->SetLineColor(kGreen);
    //get44->SetLineStyle(2);
    get44->SetMarkerStyle(4);
    get44->SetMarkerSize(0.7);
    get44->SetMarkerColorAlpha(kRed,0.35);

    // Saturation with 100 MeV resolution
    TFile *f55 = new TFile("phisat_2pi_wRes_100.root");
    TH1D *get55 = (TH1D *)f55->Get("FF");
    get55->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = 2#pi)");
    get55->GetYaxis()->SetTitle("");
    get55->Scale(197./get55->Integral(), "width");
    get55->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get55->SetLineColor(kBlue);
    //get55->SetLineStyle(2);
    get55->SetMarkerStyle(4);
    get55->SetMarkerSize(0.7);
    get55->SetMarkerColorAlpha(kRed,0.35);

    // Saturation with 150 MeV resolution
    TFile *f66 = new TFile("phisat_2pi_wRes_150.root");
    TH1D *get66 = (TH1D *)f66->Get("FF");
    get66->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = 2#pi)");
    get66->GetYaxis()->SetTitle("");
    get66->Scale(197./get66->Integral(), "width");
    get66->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get66->SetLineColor(kViolet);
    //get66->SetLineStyle(2);
    get66->SetMarkerStyle(4);
    get66->SetMarkerSize(0.7);
    get66->SetMarkerColorAlpha(kRed,0.35);

    // Saturation with 200 MeV resolution
    TFile *f61 = new TFile("phisat_2pi_wRes_200.root");
    TH1D *get61 = (TH1D *)f61->Get("FF");
    get61->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = 2#pi)");
    get61->GetYaxis()->SetTitle("");
    get61->Scale(197./get61->Integral(), "width");
    get61->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get61->SetLineColor(kViolet);
    //get61->SetLineStyle(2);
    get61->SetMarkerStyle(4);
    get61->SetMarkerSize(0.7);
    get61->SetMarkerColorAlpha(kRed,0.35);

    TCanvas *c_f1 = new TCanvas();
    get1->Draw();
    //get3->Draw("P HIST same");
    //get4->Draw("P HIST same");
    //get5->Draw("P HIST same");
    //get6->Draw("P HIST same");
    get7->Draw("P HIST same");
    //get71->Draw("P HIST same");
    get8->Draw("same");
    get22->Draw("P HIST same");
    //get33->Draw("P HIST same");
    //get44->Draw("P HIST same");
    //get55->Draw("P HIST same");
    //get66->Draw("P HIST same");
    //get61->Draw("P HIST same");
    gPad->SetLogy(1);
    gStyle->SetOptStat(0);

    auto legend = new TLegend(0.43,0.71,0.73,0.9);
    legend->AddEntry(get8,"d#sigma/dt [nb/GeV^{2}] with sat: no res","l");
    legend->AddEntry(get22,"d#sigma/dt [nb/GeV^{2}] with sat: 5 MeV res","p");
    //legend->AddEntry(get33,"d#sigma/dt [nb/GeV^{2}] with sat: 10 MeV res","p");
    //legend->AddEntry(get44,"d#sigma/dt [nb/GeV^{2}] with sat: 50 MeV res","p");
    //legend->AddEntry(get55,"d#sigma/dt [nb/GeV^{2}] with sat: 100 MeV res","p");
    //legend->AddEntry(get66,"d#sigma/dt [nb/GeV^{2}] with sat: 150 MeV res","p");
    //legend->AddEntry(get61,"d#sigma/dt [nb/GeV^{2}] with sat: 200 MeV res","p");
    legend->AddEntry(get1,"|F(t)|^{2} (no resolution)","l");
    legend->AddEntry(get7,"|F(t)|^{2} with 5 MeV resolution","p");
    //legend->AddEntry(get4,"|F(t)|^{2} with 10 MeV resolution","p");
    //legend->AddEntry(get5,"|F(t)|^{2} with 50 MeV resolution","p");
    //legend->AddEntry(get3,"|F(t)|^{2} with 100 MeV resolution","p");
    //legend->AddEntry(get6,"|F(t)|^{2} with 150 MeV resolution","p");
    //legend->AddEntry(get71,"|F(t)|^{2} with 200 MeV resolution","p");
    legend->Draw();
}

void phipi_Wsat()
{
    // theta = pi
    TFile *f1 = new TFile("72pt_wedge_ff_pi.root");
    TH1D *get1 = (TH1D *)f1->Get("FF");
    get1->Scale(197./get1->Integral(), "width");
    get1->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi)");
    get1->GetYaxis()->SetTitle("");
    get1->GetYaxis()->SetRangeUser(1e-5,1e5);
    get1->GetXaxis()->SetTitle("t [GeV^{2}]");
    get1->Scale(197./get1->Integral(), "width");
    get1->SetLineColor(kBlack);

    // theta = pi with 100 MeV resolution
    TFile *f3 = new TFile("72pt_res_wedge_ff_pi_100.root");
    TH1D *get3 = (TH1D *)f3->Get("FF");
    get3->Scale(197./get3->Integral(), "width");
    get3->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi)");
    get3->GetYaxis()->SetTitle("");
    get3->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get3->SetLineColor(kBlue);
    //get3->SetLineStyle(2);
    get3->SetMarkerStyle(3);
    get3->SetMarkerSize(0.7);
    get3->SetMarkerColorAlpha(kBlack,0.35);

    // theta = pi with 10 MeV resolution
    TFile *f4 = new TFile("72pt_res_wedge_ff_pi_10.root");
    TH1D *get4 = (TH1D *)f4->Get("FF");
    get4->Scale(197./get4->Integral(), "width");
    get4->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi)");
    get4->GetYaxis()->SetTitle("");
    get4->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get4->SetLineColor(kGreen);
    //get4->SetLineStyle(2);
    get4->SetMarkerStyle(3);
    get4->SetMarkerSize(0.7);
    get4->SetMarkerColorAlpha(kBlack,0.35);

    // theta = pi with 50 MeV resolution
    TFile *f5 = new TFile("72pt_res_wedge_ff_pi_50.root");
    TH1D *get5 = (TH1D *)f5->Get("FF");
    get5->Scale(197./get5->Integral(), "width");
    get5->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi)");
    get5->GetYaxis()->SetTitle("");
    get5->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get5->SetLineColor(kCyan);
    //get5->SetLineStyle(2);
    get5->SetMarkerStyle(3);
    get5->SetMarkerSize(0.7);
    get5->SetMarkerColorAlpha(kBlack,0.35);

    // theta = pi with 150 MeV resolution
    TFile *f6 = new TFile("72pt_res_wedge_ff_pi_150.root");
    TH1D *get6 = (TH1D *)f6->Get("FF");
    get6->Scale(197./get6->Integral(), "width");
    get6->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi)");
    get6->GetYaxis()->SetTitle("");
    get6->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get6->SetLineColor(kViolet);
    get6->SetMarkerStyle(3);
    get6->SetMarkerSize(0.7);
    get6->SetMarkerColorAlpha(kBlack,0.35);

    // theta = pi with 5 MeV resolution
    TFile *f7 = new TFile("72pt_res_wedge_ff_pi_5.root");
    TH1D *get7 = (TH1D *)f7->Get("FF");
    get7->Scale(197./get7->Integral(), "width");
    get7->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi)");
    get7->GetYaxis()->SetTitle("");
    get7->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get7->SetLineColor(kOrange);
    //get7->SetLineStyle(2);
    get7->SetMarkerStyle(3);
    get7->SetMarkerSize(0.7);
    get7->SetMarkerColorAlpha(kBlack,0.35);

    // theta = pi with 200 MeV resolution
    TFile *f71 = new TFile("72pt_res_wedge_ff_pi_200.root");
    TH1D *get71 = (TH1D *)f71->Get("FF");
    get71->Scale(197./get71->Integral(), "width");
    get71->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi)");
    get71->GetYaxis()->SetTitle("");
    get71->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get71->SetLineColor(kOrange);
    //get71->SetLineStyle(2);
    get71->SetMarkerStyle(3);
    get71->SetMarkerSize(0.7);
    get71->SetMarkerColorAlpha(kBlack,0.35);

    // with saturation
    TFile *f8 = new TFile("phisat_pi.root");
    TH1D *get8 = (TH1D *)f8->Get("FF");
    get8->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi)");
    get8->GetYaxis()->SetTitle("");
    get8->Scale(197./get8->Integral(), "width");
    get8->GetYaxis()->SetRangeUser(1e-5,1e5);
    get8->SetLineColor(kRed);
    get8->SetLineStyle(1);
    //get8->SetMarkerStyle(4);
    //get8->SetMarkerSize(0.6);
    //get8->SetMarkerColorAlpha(kRed,0.35);

    // Saturation with 5 MeV resolution
    TFile *f22 = new TFile("phisat_pi_wRes_5.root");
    TH1D *get22 = (TH1D *)f22->Get("FF");
    get22->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi)");
    get22->GetYaxis()->SetTitle("");
    get22->Scale(197./get22->Integral(), "width");
    get22->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get22->SetLineColor(kOrange);
    //get22->SetLineStyle(2);
    get22->SetMarkerStyle(4);
    get22->SetMarkerSize(0.7);
    get22->SetMarkerColorAlpha(kRed,0.35);

    // Saturation with 10 MeV resolution
    TFile *f33 = new TFile("phisat_pi_wRes_10.root");
    TH1D *get33 = (TH1D *)f33->Get("FF");
    get33->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi)");
    get33->GetYaxis()->SetTitle("");
    get33->Scale(197./get33->Integral(), "width");
    get33->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get33->SetLineColor(kCyan);
    //get33->SetLineStyle(2);
    get33->SetMarkerStyle(4);
    get33->SetMarkerSize(0.7);
    get33->SetMarkerColorAlpha(kRed,0.35);

    // Saturation with 50 MeV resolution
    TFile *f44 = new TFile("phisat_pi_wRes_50.root");
    TH1D *get44 = (TH1D *)f44->Get("FF");
    get44->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi)");
    get44->GetYaxis()->SetTitle("");
    get44->Scale(197./get44->Integral(), "width");
    get44->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get44->SetLineColor(kGreen);
    //get44->SetLineStyle(2);
    get44->SetMarkerStyle(4);
    get44->SetMarkerSize(0.7);
    get44->SetMarkerColorAlpha(kRed,0.35);

    // Saturation with 100 MeV resolution
    TFile *f55 = new TFile("phisat_pi_wRes_100.root");
    TH1D *get55 = (TH1D *)f55->Get("FF");
    get55->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi)");
    get55->GetYaxis()->SetTitle("");
    get55->Scale(197./get55->Integral(), "width");
    get55->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get55->SetLineColor(kBlue);
    //get55->SetLineStyle(2);
    get55->SetMarkerStyle(4);
    get55->SetMarkerSize(0.7);
    get55->SetMarkerColorAlpha(kRed,0.35);

    // Saturation with 150 MeV resolution
    TFile *f66 = new TFile("phisat_pi_wRes_150.root");
    TH1D *get66 = (TH1D *)f66->Get("FF");
    get66->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi)");
    get66->GetYaxis()->SetTitle("");
    get66->Scale(197./get66->Integral(), "width");
    get66->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get66->SetLineColor(kViolet);
    //get66->SetLineStyle(2);
    get66->SetMarkerStyle(4);
    get66->SetMarkerSize(0.7);
    get66->SetMarkerColorAlpha(kRed,0.35);

    // Saturation with 200 MeV resolution
    TFile *f61 = new TFile("phisat_pi_wRes_200.root");
    TH1D *get61 = (TH1D *)f61->Get("FF");
    get61->SetTitle("Wedge Integration and Saturation Model (#phi) with Resolution (#theta_{max} = #pi)");
    get61->GetYaxis()->SetTitle("");
    get61->Scale(197./get61->Integral(), "width");
    get61->GetYaxis()->SetRangeUser(1e-5,1e5);
    //get61->SetLineColor(kViolet);
    //get61->SetLineStyle(2);
    get61->SetMarkerStyle(4);
    get61->SetMarkerSize(0.7);
    get61->SetMarkerColorAlpha(kRed,0.35);

    TCanvas *c_f1 = new TCanvas();
    //get1->Draw();
    //get3->Draw("P HIST same");
    //get4->Draw("P HIST same");
    //get5->Draw("P HIST same");
    //get6->Draw("P HIST same");
    get7->Draw("P HIST same");
    //get71->Draw("P HIST same");
    //get8->Draw("same");
    get22->Draw("P HIST same");
    //get33->Draw("P HIST same");
    //get44->Draw("P HIST same");
    //get55->Draw("P HIST same");
    //get66->Draw("P HIST same");
    //get61->Draw("P HIST same");
    gPad->SetLogy(1);
    gStyle->SetOptStat(0);

    auto legend = new TLegend(0.43,0.71,0.73,0.9);
    //legend->AddEntry(get8,"d#sigma/dt [nb/GeV^{2}] with sat: no res","l");
    legend->AddEntry(get22,"d#sigma/dt [nb/GeV^{2}] with sat: 5 MeV res","p");
    //legend->AddEntry(get33,"d#sigma/dt [nb/GeV^{2}] with sat: 10 MeV res","p");
    //legend->AddEntry(get44,"d#sigma/dt [nb/GeV^{2}] with sat: 50 MeV res","p");
    //legend->AddEntry(get55,"d#sigma/dt [nb/GeV^{2}] with sat: 100 MeV res","p");
    //legend->AddEntry(get66,"d#sigma/dt [nb/GeV^{2}] with sat: 150 MeV res","p");
    //legend->AddEntry(get61,"d#sigma/dt [nb/GeV^{2}] with sat: 200 MeV res","p");
    //legend->AddEntry(get1,"|F(t)|^{2} (no resolution)","l");
    legend->AddEntry(get7,"|F(t)|^{2} with 5 MeV resolution","p");
    //legend->AddEntry(get4,"|F(t)|^{2} with 10 MeV resolution","p");
    //legend->AddEntry(get5,"|F(t)|^{2} with 50 MeV resolution","p");
    //legend->AddEntry(get3,"|F(t)|^{2} with 100 MeV resolution","p");
    //legend->AddEntry(get6,"|F(t)|^{2} with 150 MeV resolution","p");
    //legend->AddEntry(get71,"|F(t)|^{2} with 200 MeV resolution","p");
    legend->Draw();
}



                
                //========================================= Transformations: |F(t)|^2 -> G(r) =============================================//

// Integrand: |F(t)|^2 (wedge, no resolution)->G(x,y) 
double G_integrand(double *x, double *par) 
{
    double xx = par[0];
    double y = par[1];
    double r = sqrt(xx*xx + y*y);
    // integration variable
    double q = x[0];
    return sqrt(2*pi/(r))*sqrt(2*hbarc/(pi*q*r))*sin(q*r/hbarc)*wedge_FF2_inte(q)*sqrt(q)*1/hbarc*1/hbarc*1/sqrt(hbarc); // [GeV^-1 fm^-3]
}

// Integral for G(x,y)
double G_integral_dq(double x, double y)
{
    double q_min = 0.;
    double q_max = 1.;
    TF1 *f = new TF1("f",G_integrand,q_min,q_max,2);  
    f->SetParameters(x,y);
    return f->Integral(q_min,q_max,1e-12); //[1/fm^3]
}

// Draw 2D G(x,y) histogram (from wedge with no resolution)
void FF_transform_to_G()
{
    int bins = 70;
    double step = 0.3;
    int min = -10;
    int max = 10;

    TH2D *hankelTransformFF = new TH2D("hankelTransformFF","Transform Form Factor",bins,min,max,bins,min,max);
    hankelTransformFF->Sumw2();
    for(int i=min;i<bins;i++)
    {
        double x = (i+min)*step;
        for (int j=min;j<bins;j++)
        {
            double y = (j+min)*step;
            double val = G_integral_dq(x,y);
            hankelTransformFF->SetBinContent(i+1,j+1,val);
        }
        cout << "x: " << x << endl;
    }
    TFile *eAu_5_41_dsigmadydt = new TFile("transform_FF.root","recreate");
	hankelTransformFF->SetTitle("Transformed WS: F(t) #rightarrow G(r)");
    hankelTransformFF->GetYaxis()->SetTitle("y [fm]");
	hankelTransformFF->GetXaxis()->SetTitle("x [fm]");
	hankelTransformFF->Write();
	eAu_5_41_dsigmadydt->Close();

    TCanvas *c1 = new TCanvas("c1", "Hankel Transform X-proj", 800, 600);
    hankelTransformFF->Draw();
}

// Integrand: |F(t)|^2 (wedge with resolution)-> G(x,y)
double G_integrand_wedge_wRes(double *x, double *par) 
{
    double xx = par[0];
    double y = par[1];
    double r = sqrt(xx*xx + y*y);
    // integration variable
    double q = x[0];
    //double q = qy;
    return sqrt(2*pi/(r))*sqrt(2*hbarc/(pi*q*r))*sin(q*r/hbarc)*testwedge_FF2_inte(q)*sqrt(q)*1/hbarc*1/hbarc*1/sqrt(hbarc); // [GeV^-1 fm^-3]
}

// Integral for G(x,y)
double G_integral_dq_wedge_wRes(double x, double y)
{
    double q_min = 0.;
    double q_max = 1.;
    TF1 *f = new TF1("f",G_integrand_wedge_wRes,q_min,q_max,2);  
    f->SetParameters(x,y);
    return f->Integral(q_min,q_max,1e-12); //[1/fm^3]
}

// Draw 2D G(x,y) histogram (from wedge with resolution)
void FF_transform_to_G_wedge_wRes()
{
    int bins = 70;
    double step = 0.286;
    double min = -10.;
    double max = 10.;

    TH2D *hankelTransformFF = new TH2D("hankelTransformFF","Transform Form Factor",bins,min,max,bins,min,max);
    hankelTransformFF->Sumw2();
    for(int i=0;i<bins;i++)
    {
        double x = min+(i*step);
        for (int j=0; j<bins; j++)
        {
            double y = min+(j*step);
            double val = G_integral_dq_wedge_wRes(x,y);
            hankelTransformFF->SetBinContent(i+1,j+1,val);
        }
        cout << "x: " << x << endl;
    }
    TFile *eAu_5_41_dsigmadydt = new TFile("transform_wedge_pio3_150.root","recreate");
	hankelTransformFF->SetTitle("Transformed WS: F(t) #rightarrow G(r)");
    hankelTransformFF->GetYaxis()->SetTitle("y [fm]");
	hankelTransformFF->GetXaxis()->SetTitle("x [fm]");
	hankelTransformFF->Write();
	eAu_5_41_dsigmadydt->Close();

    TCanvas *c1 = new TCanvas("c1", "Hankel Transform X-proj", 800, 600);
    hankelTransformFF->Draw();
}

// 1D histogram of G(x,y) 
void draw_trans()
{
    TFile *f = new TFile("transform_FF.root");
    TH2D *td = (TH2D *)f->Get("hankelTransformFF");
    TH1D *td1 = (TH1D *)td->ProjectionY();
    td1->SetTitle("Transformed WS: F(t) #rightarrow G(r)");
    td1->GetXaxis()->SetTitle("y [fm]");
    td1->GetYaxis()->SetTitle("G(r) [fm^{-3}]"); 
    td1->Draw();
}

// G(r) no wedge saturation model
void draw_satTransform()
{
    const int bins = 70;
    int min = -10;
    int max = 10;
    double x_vals[bins] = {-8.85701049264268,-8.21604516146429,-7.86642770809425,-7.41969429545476,-6.87584492354581,-6.58449704573745,-6.35141874349076,
-6.11834044124407,-5.88526213899738,-5.76872298787404,-5.67160702860458,-5.53564468562734,-5.419105534504,-5.26371999967287,-5.14718084854953,
-5.03064169742618,-4.91410254630284,-4.75871701147171,-4.68102424405615,-4.5644850929328,-4.44794594180945,-4.33140679068611,-4.23429083141665,
-4.09832848843942,-3.90409656990051,-3.76813422692327,-3.4962095409688,-3.26313123872211,-3.01062974462153,-2.79697463422873,-2.56389633198204,
-2.33081802973535,-2.11716291934255,-1.88408461709586,-1.61215993114139,-1.3790816288947,-1.02946417552466,-0.641000338446842,0.0582345682932299,
0.5438143646405,1.04881735284166,1.35958842250392,1.70920587587395,2.07824652109788,2.40844078261403,2.6803654685685,2.91344377081519,3.14652207306188,
3.37960037530857,3.61267867755526,3.84575697980195,4.07883528204864,4.31191358429533,4.5838382702498,4.8169165724965,5.04999487474318,5.28307317698987,
5.49672828738268,5.72980658962936,5.96288489187606,6.11827042670718,6.39019511266165,6.60385022305445,6.83692852530114,7.05058363569395,7.24481555423285,
7.49731704833344,7.80808811799569,8.06058961209627,8.48789983288187};
    double y_vals[bins] = {-0.000906383912261291,-0.000166321921966014,0.000388724570755436,0.00297894153678888,0.00889943745915102,0.0153749798742346,
0.0214804912941706,0.0279560337092542,0.0353566536122069,0.0386869325685355,0.0425722580175857,0.0455325059787668,0.0488627849350955,0.0521930638914242,
0.0557083583453267,0.058483590808934,0.0618138697652627,0.0636640247410008,0.0664392572046081,0.0684744276779201,0.0706946136488059,
0.0725447686245441,0.0747649545954299,0.0764300940735942,0.0793903420347753,0.0823505899959564,0.0849408069619898,0.0877160394255971,0.089751209898909,
0.0921564113673687,0.093821550845533,0.0951166593285497,0.0965967833091403,0.0977068762945832,0.0991870002751737,0.100112077763043,0.101037155250912,
0.101962232738781,0.102332263733929,0.101777217241207,0.101037155250912,0.100112077763043,0.0991870002751737,0.0965967833091403,0.0949316438309759,
0.0936365353479592,0.0906762873867781,0.0893811789037614,0.0869759774353018,0.0832756674838255,0.0810554815129397,0.0779102180541847,0.074579939097856,
0.0705095981512321,0.065514179716739,0.061443838770115,0.0568184513307696,0.0501578934181122,0.0436823510030286,0.0381318860758141,0.0316563436607305,
0.0244407392553516,0.0185202433329895,0.0131547939033488,0.00871442196157719,0.00519912750767468,0.00297894153678888,0.00168383305377216,
0.000388724570755436,-0.000536352917113653};

    TH1D *hist = new TH1D("hist", "Coherent Saturation;|t| [GeV^{2}];d#sigma/dt [nb/GeV^{2}] #rightarrow G(r)", bins, min, max);
    for (int i = 0; i < bins; ++i) 
    {
        hist->SetBinContent(i+1, y_vals[i]);
    }
    TCanvas *c1 = new TCanvas("c1", "True Form Factor", 800, 600);
    hist->Draw();
    TFile *file = new TFile("test_sat_transform.root", "RECREATE");
    hist->Write();
    file->Close();
}



*/

