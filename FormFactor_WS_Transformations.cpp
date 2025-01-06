#include <iostream>
#include "TComplex.h"
#include "TVector2.h"
#include "TF2.h"
#include "TVirtualFFT.h"
#include "EICvalueconst.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"
#include <Riostream.h>
#include "TLegend.h"
#include "TLegendEntry.h"
#include "Math/IFunction.h"
#include <cmath>
#include "TSystem.h"
#include "TAxis.h"
#include "TPaveLabel.h"
#include "TFormula.h"



using namespace std;

            //======================= Original Woods-Saxon: G(r) and Form Factor: F(q) =====================================//

double calcWS(double r)
{
    double Vo = 2.12;  
    double A = 197.;
    double a = depth_Au;
    double R = 6.38;
    return Vo/(1. + exp((r-R)/a)); // [1/fm^3]
}

void getWS()
{
    double bins = 1000;
    double step = 0.015;
    double min = 0.;
    double max = 15;
    TH1D *ws = new TH1D("ws","Woods-Saxon",bins,min,max);
    ws->Sumw2();
    for(int i=0;i<bins;i++)
    {
        double r = (i+1)*step;
        double a = calcWS(r);
        ws->SetBinContent(i+1,a);
        cout << "r: " << r << endl;
    }
	TFile *FF_WS_Transformations = new TFile("FF_WS_Transformations_real_WS.root","recreate");
	ws->SetTitle("Woods-Saxon");
    ws->GetYaxis()->SetTitle("G(r) [fm^{-3}]");
	ws->GetXaxis()->SetTitle("r [fm]");
	ws->Write();
	FF_WS_Transformations->Close();
}

void drawWS() 
{
    TFile *f = new TFile("FF_WS_Transformations_real_WS.root");
    TH1D *td = (TH1D *)f->Get("ws");
    TCanvas *drawws = new TCanvas();
	td->Draw();   
}

double Gr_ws_integrand(double *x, double *par) 
{
    double r = x[0];
    double q = par[0];
    return calcWS(r)*r*r;
}

void Gr_ws_integral()
{
    double r_min = 0.;
    double r_max = 15;
    TF1 *f = new TF1("f",Gr_ws_integrand,r_min,r_max);  
    double q = 0;
    f->SetParameters(q);
    cout << "Integral Value: " << f->Integral(r_min,r_max,1e-6) << endl; //[-]
}

double calcFF(double q)
{
	if(q==0){return 0;}
    else{
        double A = 197;
        double ro = 1.25;
        double rho0 = 2.12; // 0.181
        double R = 6.38; //ro*5.81865;  // A^(1/3) = 5.81865 >> R = 7.27331 
	    const double arg1 = q * R / hbarc;  
	    const double arg2 = hbarc /q;  // RWS_Au = 6.38 ( from .h file)
	    const double sph  = (sin(arg1) - arg1 * cos(arg1)) * 4*pi*rho0 * arg2 * arg2 *arg2/ double(A); 
	    const double a0   = .70;
	    return sph * 1/(1. + (a0 * a0 * q*q/(hbarc*hbarc)))*1/0.00163359;  // [-]
    }
}

void getFF()
{
    double bins = 1000.;
    double step = 0.0005;
    double min = 0.0005;
    double max = 0.5;
    TH1D *FF = new TH1D("FF","Form Factor",bins,min,max);
    FF->Sumw2();

    for(int i=0;i<bins;i++)
    {
        double q = (i+1)*step;
        double a = calcFF(q);
        FF->SetBinContent(i+1,a);
        cout << "q: " << q << endl;
    }
	TFile *FF_WS_Transformations = new TFile("FF_WS_Transformations_real_FF.root","recreate");
	FF->SetTitle("Form Factor: F(q)");
    FF->GetYaxis()->SetTitle("F(q)");
	FF->GetXaxis()->SetTitle("q [GeV]");
	FF->Write();
	FF_WS_Transformations->Close();
}

void drawFF() 
{
    TFile *f = new TFile("FF_WS_Transformations_real_FF.root");
    TH1D *td = (TH1D *)f->Get("FF");
    TCanvas *drawws = new TCanvas();
	td->Draw();   
}

double FF_integrand(double *x, double *par) 
{
    double q = x[0];
    double r = par[0];
    return calcFF(q)*q*q;
}

void FF_integral()
{
    double q_min = 0.0005;
    double q_max = 5;
    TF1 *f = new TF1("f",FF_integrand,q_min,q_max);  
    double r = 0;
    f->SetParameters(r);
    cout << "Integral Value: " << f->Integral(q_min,q_max,1e-6) << endl; //[-]
}


                        //================================== Transformations =====================================//

//***Transforming Form Factor: F(q)->G(r)->F(q)***//
// F(q) -> G(r) (uses real Form Factor)
double ff_integrand_dq(double *x, double *par) 
{
    double q = x[0];
    double r = par[0];
    return  2 * pi * sqrt(2 * pi/(r)) *sqrt(2*hbarc/(pi*q*r)) * sin(q *r/hbarc) *calcFF(q) * q *sqrt(q) * 1/hbarc * 1/hbarc * 1/sqrt(hbarc); //[1/fm^3Gev] 
}

double ff_integral_dq(double r)
{
    double q_min = 0.0005;
    double q_max = 5;
    TF1 *f = new TF1("f",ff_integrand_dq,q_min,q_max,1);  
    f->SetParameters(r);
    return f->Integral(q_min,q_max,1e-6)*196.253/231.072; //[1/fm^3]
}

void ff_transform_to_r()
{
    int bins = 1000.;
    double step = 0.015;
    double min = 0.;
    double max = 15;

    TH1D *hankelTransformFF = new TH1D("hankelTransformFF","Transform Form Factor",bins,min,max);
    hankelTransformFF->Sumw2();
    double b[1000];
    for(int i=0;i<bins;i++)
    {
        double r = (i+1)*step;
	    double a = ff_integral_dq(r);
        b[i] = a;
        hankelTransformFF->SetBinContent(i+1,a);
        cout << "r: " << r << endl;
    }
    TFile *eAu_5_41_dsigmadydt = new TFile("FF_WS_Transformations_Fq_to_Gr.root","recreate");
	hankelTransformFF->SetTitle("Transformed WS: F(q) #rightarrow G(r)");
    hankelTransformFF->GetYaxis()->SetTitle("G(r) [fm^{-3}]");
	hankelTransformFF->GetXaxis()->SetTitle("r [fm]");
	hankelTransformFF->Write();
	eAu_5_41_dsigmadydt->Close();

    TCanvas *c1 = new TCanvas("c1", "Hankel Transform", 800, 600);
    hankelTransformFF->Draw();
    c1->Update();
}

double Gr_ff_integrand_dq(double *x, double *par) 
{
    double r = x[0];
    double q = par[0];
    return ff_integral_dq(r)*r*r; //[-]
}

void Gr_ff_integral_dq()
{
    double r_min = 0.;
    double r_max = 15;
    TF1 *f = new TF1("f",Gr_ff_integrand_dq,r_min,r_max);  
    double q = 0;
    f->SetParameters(q);
    cout << "Integral Value: " << f->Integral(r_min,r_max,1e-6)  << endl;
}


//G(r) -> F(q) (uses transformed form factor) ie: F(q) -> G(r) -> F(q)
double G_integrand_dr(double *x, double *par) 
{
    double r = x[0];
    double q = par[0];
    return 2 * pi * sqrt(2 * pi/q) * sqrt(2*hbarc/(pi*q*r)) * sqrt(hbarc) * ff_integral_dq(r) * r * sqrt(r) * sin(q *r/hbarc);
}

double G_integral_dr(double q)
{
    double r_min = 0.;
    double r_max = 30;
    TF1 *f = new TF1("f",G_integrand_dr,r_min,r_max,1);  
    f->SetParameters(q);
    return f->Integral(r_min,r_max,1e-3); // [-]
}

void G_transform_to_F() 
{
    int bins = 1000.;
    double step = 0.0005;
    double min = 0.0005;
    double max = 0.5;
    TH1D *hankelTransform = new TH1D("hankelTransform","Transform Woods-Saxon: G(r) #rightarrow F(q)",bins,min,max);
    hankelTransform->Sumw2();
    int b[bins];
    for(int i=0;i<bins;i++)
    {
        double q = (i+1)*step;
	    double a = G_integral_dr(q);
        b[i] = a;
        hankelTransform->SetBinContent(i+1,a);
        cout << "q:" << q << endl;
    }
    TFile *FF_WS_Transformations = new TFile("FF_WS_Transformations_Gr_to_Fq.root","recreate");
	hankelTransform->SetTitle("Transformed WS: G(r) #rightarrow F(q)");
    hankelTransform->GetYaxis()->SetTitle("F(q)");
	hankelTransform->GetXaxis()->SetTitle("q [GeV]");
	hankelTransform->Write();
	FF_WS_Transformations->Close();

    TCanvas *c1 = new TCanvas("c1", "Hankel Transform", 800, 600);
    hankelTransform->Draw();
    c1->Update();
}

double ff_G_integrand_dq(double *x, double *par) 
{
    double q = x[0];
    double r = par[0];
    return G_integral_dr(q)*q*q; //[-]
}

void ff_G_integral_dq()
{
    double q_min = 0.0005;
    double q_max = 5;
    TF1 *f = new TF1("f",ff_G_integrand_dq,q_min,q_max);  
    double r = 0;
    f->SetParameters(r);
    cout << "thinking..." << endl;
    cout << "Integral Value: " << f->Integral(q_min,q_max,1e-6)  << endl;
    cout << "done" << endl;
}


//***Transforming Woods-Saxon: rho(r)->F(q)->G(r)***//
// rho(r) -> F(q) (uses real Woods-Saxon)
double ws_integrand_dr(double *x, double *par) 
{
    double r = x[0];
    double q = par[0];
    return 2 * pi * sqrt(2 * pi/q) * sqrt(2*hbarc/(pi*q*r)) * calcWS(r) * sqrt(hbarc) * r * sqrt(r) * sin(q *r/hbarc); //[1/fm]
}

double ws_integral_dr(double q)
{
    double r_min = 0.;
    double r_max = 30;
    TF1 *f = new TF1("f",ws_integrand_dr,r_min,r_max,1);  
    f->SetParameters(q);
    return f->Integral(r_min,r_max,1e-3)*1/0.316509;// [-]
}

void ws_transform_to_q()
{
    int bins = 1000.;
    double step = 0.0005;
    double min = 0.0005;
    double max = 0.5;

    TH1D *hankelTransformWS = new TH1D("hankelTransformWS","Transform Woods-Saxon",bins,min,max);
    hankelTransformWS->Sumw2();
    int b[bins];
    for(int i=0;i<bins;i++)
    {
        double q = (i+1)*step;
	    double a = ws_integral_dr(q);
        b[i] = a;
        hankelTransformWS->SetBinContent(i+1,a);
        cout << "q:" << q << endl;
    }
    TFile *eAu_5_41_dsigmadydt = new TFile("transform_WS_to_F_of_q.root","recreate");
	hankelTransformWS->SetTitle("Transformed WS: #rho(r) #rightarrow F(q)");
    hankelTransformWS->GetYaxis()->SetTitle("F(q)");
	hankelTransformWS->GetXaxis()->SetTitle("q [GeV]");
	hankelTransformWS->Write();
	eAu_5_41_dsigmadydt->Close();

    TCanvas *c1 = new TCanvas("c1", "Hankel Transform", 800, 600);
    hankelTransformWS->Draw();
    c1->Update();
}

double ff_ws_integrand_dq(double *x, double *par) 
{
    double q = x[0];
    double r = par[0];
    return ws_integral_dr(q)*q*q; // [-]
}

void ff_ws_integral_dq()
{
    double q_min = 0.0005;
    double q_max = 5;
    TF1 *f = new TF1("f",ff_ws_integrand_dq,q_min,q_max);  
    double r = 0;
    f->SetParameters(r);
    cout << "Integral Value: " << f->Integral(q_min,q_max,1e-6) << endl;
}


// F(q) -> G(r) (uses transformed Woods-Saxon) ie: rho(r) -> F(q) -> G(r)
double ws_integrand_dq(double *x, double *par)
{
    double q = x[0];
    double r = par[0];
    return 2 * pi * sqrt(2*pi/r) * sqrt(2*hbarc/(pi*q*r)) * ws_integral_dr(q) * sin(q*r/hbarc) * q * sqrt(q) * 1/hbarc * 1/hbarc * 1/sqrt(hbarc);  // [1/fm^3GeV]
}

double ws_integral_dq(double r)
{
    double q_min = 0.0005;
    double q_max = 5;
    TF1 *f = new TF1("f",ws_integrand_dq,q_min,q_max,1);
    f->SetParameters(r);
    return f->Integral(q_min,q_max,1e-3)*196.253/48680.4; // [1/fm^3]
}

void ws_transform_to_r()
{
    int bins = 1000.;
    double step = 0.015;
    double min = 0.;
    double max = 15;

    TH1D *hankelTransformWS = new TH1D("hankelTransformWS","Transform Woods-Saxon",bins,min,max);
    hankelTransformWS->Sumw2();
    double b[bins];
    for(int i=0;i<bins;i++)
    {
        double r = (i+1)*step;
	    double a = ws_integral_dq(r);
        b[i] = a;
        hankelTransformWS->SetBinContent(i+1,a);
        cout << "r: " << r << endl;
    }
    TFile *eAu_5_41_dsigmadydt = new TFile("transform_transform_WS_of_q_to_WS_of_r.root","recreate");
	hankelTransformWS->SetTitle("Transformed-Transformed WS: F(t) #rightarrow #rho(r)");
    hankelTransformWS->GetYaxis()->SetTitle("G(r) [fm^{-3}]");
	hankelTransformWS->GetXaxis()->SetTitle("r [fm]");
	hankelTransformWS->Write();
	eAu_5_41_dsigmadydt->Close();

    TCanvas *c1 = new TCanvas("c1", "Hankel Transform", 800, 600);
    hankelTransformWS->Draw();
    c1->Update();
}

double Gr_ws_integrand_dq(double *x, double *par) 
{
    double r = x[0];
    double q = par[0];
    return ws_integral_dq(r)*r*r; // [-]
}

void Gr_ws_integral_dq()
{
    double r_min = 0.;
    double r_max = 15;
    TF1 *f = new TF1("f",Gr_ws_integrand_dq,r_min,r_max);  
    double q = 0;
    f->SetParameters(q);
    cout << "Integral Value: " << f->Integral(r_min,r_max,1e-6) << endl;
}


//***Compare all 3 G(r)***//

void stack_all_Gr()
{
    // Real Woods-Saxon: rho(r)
    TFile *f1 = new TFile("FF_WS_Transformations_real_WS.root");
    TH1D *get1 = (TH1D *)f1->Get("ws");
    get1->SetTitle("Woods-Saxon G(r)");
    get1->GetYaxis()->SetTitle("G(r) [fm^{-3}]");
    get1->GetXaxis()->SetTitle("r [fm]");
    get1->GetYaxis()->SetRangeUser(0,2.5);
    get1->SetLineColor(kBlack);
    //TH1*h0 = (TH1*)(get1->Clone("h0"));
    //h0->Scale(1./h0->Integral(), "width");
    //h0->SetLineColor(kBlack);

    // G(r) from real form factor: F(q)->G(r)
    TFile *f2 = new TFile("FF_WS_Transformations_Fq_to_Gr.root");
    TH1D *get2 = (TH1D *)f2->Get("hankelTransformFF");
    get2->SetLineColor(kRed);
    //TH1*h1 = (TH1*)(get2->Clone("h1"));
    //h1->Scale(1./h1->Integral(), "width");
    //h1->SetLineColor(kRed);

    // G(r) from double transformed Woods-Saxon: rho(r)->F(q)->G(r)
    TFile *f3 = new TFile("transform_transform_WS_of_q_to_WS_of_r.root");
    TH1D *get3 = (TH1D *)f3->Get("hankelTransformWS");
    get3->SetLineColor(kMagenta);
    get3->SetLineStyle(2);
    //TH1*h2 = (TH1*)(get3->Clone("h2"));
    //h2->Scale(1./h2->Integral(), "width");
    //h2->SetLineColor(kBlue);

    TCanvas *c_f1 = new TCanvas();
	get1->Draw();
    get2->Draw("same");
    get3->Draw("same");
    //h1->Draw("same");
    //h2->Draw("same");

    auto legend = new TLegend(0.43,0.78,0.63,0.9);
	legend->AddEntry(get1,"Real Woods-Saxon: G(r)","l");
    legend->AddEntry(get2,"F(q) #rightarrow G(r)","l");
    legend->AddEntry(get3,"G(r) #rightarrow F(q) #rightarrow G(r)","l");
    //legend->AddEntry(h1,"F(q) #rightarrow G(r)","l");
    //legend->AddEntry(h2,"G(r) #rightarrow F(q) #rightarrow G(r)","l");
    legend->Draw();
}

//***Compare all 3 F(q)***//

void stack_all_Fq()
{
    // Real Form Factor: F(q)
    TFile *f1 = new TFile("FF_WS_Transformations_real_FF.root");
    TH1D *get1 = (TH1D *)f1->Get("FF");
    get1->SetTitle("Form Factor F(q)");
    get1->GetYaxis()->SetTitle("F(q)");
    get1->GetXaxis()->SetTitle("q [GeV]");
    //get1->GetYaxis()->SetRangeUser(1e-9,15);
    get1->SetLineColor(kBlack);

    // F(q) from real Woods-Saxon: rho(r)->F(q)
    TFile *f3 = new TFile("transform_WS_to_F_of_q.root");
    TH1D *get3 = (TH1D *)f3->Get("hankelTransformWS");
    get3->SetLineColor(kRed);
    //TH1*h1 = (TH1*)(get3->Clone("h1"));
    //h1->Scale(1./h1->Integral(), "width");
    //h1->SetLineColor(kRed);

    // F(q) from double transformed Form Factor: F(q)->G(r)->F(q)
    TFile *f2 = new TFile("FF_WS_Transformations_Gr_to_Fq.root");
    TH1D *get2 = (TH1D *)f2->Get("hankelTransform");
    get2->SetLineColor(kBlue);
    //TH1*h2 = (TH1*)(get2->Clone("h2"));
    //h2->Scale(1./h2->Integral(), "width");
    //h2->SetLineColor(kBlue);

    TCanvas *c_f1 = new TCanvas();
	get1->Draw();
    get3->Draw("same");
    get2->Draw("same");
    //h1->Draw("same");
    //h2->Draw("same");

    auto legend = new TLegend(0.43,0.78,0.63,0.9);
	legend->AddEntry(get1,"Real F(q)","l");
    legend->AddEntry(get3,"G(r) #rightarrow F(q)","l");
    legend->AddEntry(get2,"F(q) #rightarrow G(r) #rightarrow F(q)","l");
    //legend->AddEntry(h1,"G(r) #rightarrow F(q)","l");
    //legend->AddEntry(h2,"F(q) #rightarrow G(r) #rightarrow F(q)","l");
    legend->Draw();
}