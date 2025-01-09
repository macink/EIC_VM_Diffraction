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

    TCanvas *c1 = new TCanvas("c1", "True Woods-Saxon", 800, 600);
    ws->Draw();
    c1->Update();
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
	    return sph * 1/(1. + (a0 * a0 * q*q/(hbarc*hbarc)));  // [-]
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
        double t = (i+1)*step;
        double q = sqrt(t);
        double a = calcFF(q)*calcFF(q);
        FF->SetBinContent(i+1,a);
        cout << "q: " << q << endl;
    }
	TFile *FF_WS_Transformations = new TFile("FF_WS_Transformations_real_FF.root","recreate");
	FF->SetTitle("Form Factor: F(t)");
    FF->GetYaxis()->SetTitle("F(t)");
	FF->GetXaxis()->SetTitle("t [GeV^2]");
	FF->Write();
	FF_WS_Transformations->Close();

    TCanvas *c1 = new TCanvas("c1", "True Form Factor", 800, 600);
    FF->Draw();
    c1->Update();
}




                        //================================== Transformations =====================================//

//***Transforming Form Factor: Ftrue(q)->G1(r)->F1(q)***//
// Ftrue(q) -> G1(r) (uses real Form Factor)    **** G1 should == Gtrue
double G1_integrand_dq(double *x, double *par) 
{
    double q = x[0];
    double r = par[0];
    return  2 * pi * sqrt(2 * pi/(r)) *sqrt(2*hbarc/(pi*q*r)) * sin(q *r/hbarc) *calcFF(q) * q *sqrt(q) * 1/hbarc * 1/hbarc * 1/sqrt(hbarc); //[1/fm^3Gev] 
}

double G1_integral_dq(double r)
{
    double q_min = 0.0005;
    double q_max = 5;
    TF1 *f = new TF1("f",G1_integrand_dq,q_min,q_max,1);  
    f->SetParameters(r);
    return f->Integral(q_min,q_max,1e-6); //[1/fm^3]
}

// G1(r)
void F_transform_to_G1()
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
	    double a = G1_integral_dq(r);
        b[i] = a;
        hankelTransformFF->SetBinContent(i+1,a);
        cout << "r: " << r << endl;
    }
    TFile *eAu_5_41_dsigmadydt = new TFile("FF_WS_Transformations_Fq_to_Gr.root","recreate");
    //hankelTransformFF->Scale(197/hankelTransformFF->Integral());
	hankelTransformFF->SetTitle("Transformed WS: F(q) #rightarrow G(r)");
    hankelTransformFF->GetYaxis()->SetTitle("G(r) [fm^{-3}]");
	hankelTransformFF->GetXaxis()->SetTitle("r [fm]");
	hankelTransformFF->Write();
	eAu_5_41_dsigmadydt->Close();

    TCanvas *c1 = new TCanvas("c1", "Hankel Transform", 800, 600);
    hankelTransformFF->Draw();
    c1->Update();
}


//G1(r) -> F1(q)        ****F1 should == Ftrue
double F1_integrand_dr(double *x, double *par) 
{
    double r = x[0];
    double q = par[0];
    return 2 * pi * sqrt(2 * pi/q) * sqrt(2*hbarc/(pi*q*r))  * sin(q *r/hbarc) * sqrt(hbarc) * G1_integral_dq(r) * r * sqrt(r);
}

double F1_integral_dr(double q)
{
    double r_min = 0.;
    double r_max = 30;
    TF1 *f = new TF1("f",F1_integrand_dr,r_min,r_max,1);  
    f->SetParameters(q);
    return f->Integral(r_min,r_max,1e-6); // [-]
}

// F1(q) == Ftrue
void G1_transform_to_F1() 
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
	    double a = F1_integral_dr(q);
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


//***Transforming Woods-Saxon: Gtrue(r)->F2(q)->G2(r)***//
// Gtrue(r) -> F2(q) (uses real Woods-Saxon)      ****F2 should == Ftrue
double F2_integrand_dr(double *x, double *par) 
{
    double r = x[0];
    double q = par[0];
    return 2 * pi * sqrt(2 * pi/q) * sqrt(2*hbarc/(pi*q*r))  * sin(q *r/hbarc) * calcWS(r) * sqrt(hbarc) * r * sqrt(r); //[1/fm]
}

double F2_integral_dr(double q)
{
    double r_min = 0.;
    double r_max = 30;
    TF1 *f = new TF1("f",F2_integrand_dr,r_min,r_max,1);  
    f->SetParameters(q);
    return f->Integral(r_min,r_max,1e-6);// [-]
}

// F2(q)
void G_transform_to_F2()
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
	    double a = F2_integral_dr(q);
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


// F2(q) -> G2(r)   **G2 should == Gtrue
double G2_integrand_dq(double *x, double *par)
{
    double q = x[0];
    double r = par[0];
    return 2 * pi * sqrt(2*pi/r) * sqrt(2*hbarc/(pi*q*r))  * sin(q *r/hbarc) * F2_integral_dr(q) * q * sqrt(q) * 1/hbarc * 1/hbarc * 1/sqrt(hbarc);  // [1/fm^3GeV]
}

double G2_integral_dq(double r)
{
    double q_min = 0.0005;
    double q_max = 5;
    TF1 *f = new TF1("f",G2_integrand_dq,q_min,q_max,1);
    f->SetParameters(r);
    return f->Integral(q_min,q_max,1e-6); // [1/fm^3]
}

//G2(r)
void F2_transform_to_G2()
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
	    double a = G2_integral_dq(r);
        b[i] = a;
        hankelTransformWS->SetBinContent(i+1,a);
        cout << "r: " << r << endl;
    }
    TFile *eAu_5_41_dsigmadydt = new TFile("transform_transform_WS_of_q_to_WS_of_r.root","recreate");
	hankelTransformWS->SetTitle("Transformed-Transformed WS: F(q) #rightarrow #rho(r)");
    hankelTransformWS->GetYaxis()->SetTitle("G(r) [fm^{-3}]");
	hankelTransformWS->GetXaxis()->SetTitle("r [fm]");
	hankelTransformWS->Write();
	eAu_5_41_dsigmadydt->Close();

    TCanvas *c1 = new TCanvas("c1", "Hankel Transform", 800, 600);
    hankelTransformWS->Draw();
    c1->Update();
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
    get1->GetYaxis()->SetRangeUser(0,0.012);
    get1->Scale(197./get1->Integral(), "width");
    get1->SetLineColor(kBlack);
    //TH1*h0 = (TH1*)(get1->Clone("h0"));
    //h0->Scale(1./h0->Integral(), "width");
    //h0->SetLineColor(kBlack);

    // Ftrue(q)->G1(r)
    TFile *f2 = new TFile("FF_WS_Transformations_Fq_to_Gr.root");
    TH1D *get2 = (TH1D *)f2->Get("hankelTransformFF");
    get2->Scale(197./get2->Integral(), "width");
    get2->SetLineColor(kRed);
    //TH1*h1 = (TH1*)(get2->Clone("h1"));
    //h1->Scale(1./h1->Integral(), "width");
    //h1->SetLineColor(kRed);

    // Gtrue(r)->F2(q)->G2(r)
    TFile *f3 = new TFile("transform_transform_WS_of_q_to_WS_of_r.root");
    TH1D *get3 = (TH1D *)f3->Get("hankelTransformWS");
    get3->Scale(197./get3->Integral(), "width");
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
    // Real Form Factor: Ftrue(q)
    TFile *f1 = new TFile("FF_WS_Transformations_real_FF.root");
    TH1D *get1 = (TH1D *)f1->Get("FF");
    get1->SetTitle("Form Factor F(q)");
    get1->GetYaxis()->SetTitle("F(q)");
    get1->GetXaxis()->SetTitle("q [GeV]");
    get1->Scale(197./get1->Integral(), "width");
    //get1->GetYaxis()->SetRangeUser(1e-9,15);
    get1->SetLineColor(kBlack);

    // Gtrue(r)->F1(q)
    TFile *f3 = new TFile("transform_WS_to_F_of_q.root");
    TH1D *get3 = (TH1D *)f3->Get("hankelTransformWS");
    get3->Scale(197./get3->Integral(), "width");
    get3->SetLineColor(kRed);
    //TH1*h1 = (TH1*)(get3->Clone("h1"));
    //h1->Scale(1./h1->Integral(), "width");
    //h1->SetLineColor(kRed);

    // Ftrue(q)->G1(r)->F2(q)
    TFile *f2 = new TFile("FF_WS_Transformations_Gr_to_Fq.root");
    TH1D *get2 = (TH1D *)f2->Get("hankelTransform");
    get2->Scale(197./get2->Integral(), "width");
    get2->SetLineColor(kMagenta);
    get2->SetLineStyle(2);
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


//=========================================== 2D ========================================//

/*
double calcWS_2D(double x, double y)
{
    double Vo = 2.12;  
    double A = 197.;
    double a = depth_Au;
    double R = 6.38;
    double r = sqrt(x*x + y*y);
    return Vo/(1. + exp((r-R)/a)); // [1/fm^3]
}

void getWS_2D()
{
    double bins = 1000;
    double step = .015;
    double min = 0.;
    double max = 15;
    TH2D *ws = new TH2D("ws","Woods-Saxon",bins,min,max,bins,min,max);
    ws->Sumw2();
    for(int i=0;i<bins;i++)
    {
        double x = (i+1)*step;
        for(int j=0; j<bins; j++)
        {
            double y = (j+1)*step;
            double a = calcWS_2D(x,y);
            ws->SetBinContent(i+1,j+1,a);
        }
        cout << "r: " << x << endl;
    }
	TFile *FF_WS_Transformations = new TFile("FF_WS_Transformations_real_WS_2D.root","recreate");
	ws->SetTitle("Woods-Saxon");
    ws->GetYaxis()->SetTitle("G(r) [fm^{-3}]");
	ws->GetXaxis()->SetTitle("r [fm]");
	ws->Write();
	FF_WS_Transformations->Close();

    TCanvas *c1 = new TCanvas("c1", "True Woods-Saxon", 800, 600);
    //ws->ProjectionY()->Draw();
    ws->Draw();
}

double calcFF_2D(double qx, double qy)
{
    double q = sqrt(qx*qx + qy*qy);
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
	    return sph * 1/(1. + (a0 * a0 * q*q/(hbarc*hbarc)));  // [-]
    }
}

void getFF_2D()
{
    double bins = 1000.;
    double step = 0.0001;
    double min = 0.0001;
    double max = 0.1;
    TH2D *FF = new TH2D("FF","Form Factor",bins,min,max,bins,min,max);
    FF->Sumw2();

    for(int i=0;i<bins;i++)
    {
        double tx = (i+1)*step;
        double qx = sqrt(tx);
        //double qx = (i+1)*step;
        for(int j=0; j<bins; j++)
        {
            double ty = (j+1)*step;
            double qy = sqrt(ty);
            //double qy = (j+1)*step;
            double a = calcFF_2D(qx,qy);
            FF->SetBinContent(i+1,j+1,a);
        }
        cout << "qx: " << qx << endl;
    }
    TFile *FF_WS_Transformations = new TFile("FF_WS_Transformations_real_FF_2D_test.root","recreate");
	//TFile *FF_WS_Transformations = new TFile("FF_WS_Transformations_real_FF_2D.root","recreate");
	FF->SetTitle("Form Factor: F(q)");
    FF->GetYaxis()->SetTitle("F(q)");
	FF->GetXaxis()->SetTitle("q [GeV]");
	FF->Write();
	FF_WS_Transformations->Close();

    TCanvas *c1 = new TCanvas("c1", "True Form Factor", 800, 600);
    //FF->ProjectionY()->Draw();
    FF->Draw();
}
*/

                        //================================== Transformations =====================================//

//***Transforming Form Factor: Ftrue(q)->G1(r)->F1(q)***//
// Ftrue(q) -> G1(r) (uses real Form Factor)    **** G1 should == Gtrue
/*
double G1_integrand_dq_2D(double *x, double *par) 
{
    double qx = x[0];
    double qy = x[1];
    double xx = par[0];
    double y = par[1];
    double r = sqrt(xx*xx + y*y);
    double q = sqrt(qx*qx + qy*qy);
    return  4 * pi * pi * sqrt(2 * pi/(r)) * sqrt(2*hbarc/(pi*q*r)) * sin(q *r/hbarc) * calcFF_2D(qx,qy) *sqrt(q) * 1/hbarc * 1/hbarc * 1/sqrt(hbarc);
    //return  4 * pi * pi * sqrt(2 * pi/(r)) *sqrt(2*hbarc/(pi*q*r)) * sin(q *r/hbarc) *calcFF_2D(qx,qy) * qx*qy/(qx*qx+qy*qy) *sqrt(q) * 1/hbarc * 1/hbarc * 1/sqrt(hbarc); //[1/fm^3Gev] 
}

double G1_integral_dq_2D(double x, double y)
{
    double r = sqrt(x*x + y*y);
    double q_min = 0.0005;
    double q_max = 1;
    TF2 *f = new TF2("f",G1_integrand_dq_2D,q_min,q_max,q_min,q_max,2);  
    f->SetParameters(x,y);
    return f->Integral(q_min,q_max,q_min,q_max,1e-6); //[1/fm^3]
}

// G1(r)
void F_transform_to_G1_2D()
{
    int bins = 1000.;
    double step = .015;
    double min = 0.;
    double max = 15;

    TH2D *hankelTransformFF = new TH2D("hankelTransformFF","Transform Form Factor",bins,min,max,bins,min,max);
    hankelTransformFF->Sumw2();
    for(int i=0;i<bins;i++)
    {
        double x = (i+1)*step;
        for (int j =0; j<bins; j++)
        {
            double y = (j+1)*step;
            double a = G1_integral_dq_2D(x,y);
            hankelTransformFF->SetBinContent(i+1,j+1,a);
        }
        cout << "x: " << x << endl;
    }
    TFile *eAu_5_41_dsigmadydt = new TFile("FF_WS_Transformations_Fq_to_Gr_2D.root","recreate");
	hankelTransformFF->SetTitle("Transformed WS: F(q) #rightarrow G(r)");
    hankelTransformFF->GetYaxis()->SetTitle("G(r) [fm^{-3}]");
	hankelTransformFF->GetXaxis()->SetTitle("r [fm]");
	hankelTransformFF->Write();
	eAu_5_41_dsigmadydt->Close();

    TCanvas *c1 = new TCanvas("c1", "Hankel Transform X-proj", 800, 600);
    //hankelTransformFF->ProjectionX()->Draw();

    //TCanvas *c2 = new TCanvas("c2", "Hankel Transform Y-proj", 800, 600);
    //hankelTransformFF->ProjectionY()->Draw();
    hankelTransformFF->Draw();
}


//G1(r) -> F1(q)        ****F1 should == Ftrue
double F1_integrand_dr_2D(double *x, double *par) 
{
    double xx = x[0];
    double yy = x[1];
    double r = sqrt(xx*xx + yy*yy);
    double qx = par[0];
    double qy = par[1];
    double q = sqrt(qx*qx + qy*qy);
    //cout << "r: " << r << endl;
    //return 4 * pi * pi * sqrt(2 * pi/q) * sqrt(2*hbarc/(pi*q*r)) * sqrt(hbarc) * G1_integral_dq_2D(xx,yy) * xx*yy/(xx*xx+yy*yy) * sqrt(r) * sin(q *r/hbarc);
    return 4 * pi * pi * sqrt(2*pi/q) * sqrt(2*hbarc/(pi*q*r)) * sin(q *r/hbarc) * sqrt(hbarc) * G1_integral_dq_2D(xx,yy) * sqrt(r);
}

double F1_integral_dr_2D(double qx, double qy)
{
    double r_min = 0.;
    double r_max = 30;
    TF2 *f = new TF2("f",F1_integrand_dr_2D,r_min,r_max,r_min,r_max,2);  
    f->SetParameters(qx,qy);
    cout << "hi" << endl;
    return f->Integral(r_min,r_max,r_min,r_max,1e-6); // [-]
}

double F1_integrand_dr_2D_test(double *x, double *par) 
{
    double xx = x[0];
    double yy = x[1];
    double r = sqrt(xx*xx + yy*yy);
    double qx = par[0];
    double qy = par[1];
    double q = sqrt(qx*qx + qy*qy);
    //return  4 * pi * pi * sqrt(2 * pi/(r)) *sqrt(2*hbarc/(pi*q*r)) * sin(q *r/hbarc) *calcFF_2D(qx,qy) * qx*qy/(qx*qx+qy*qy) *sqrt(q) * 1/hbarc * 1/hbarc * 1/sqrt(hbarc); //[1/fm^3Gev] 
    return  4 * pi * pi * sqrt(2 * pi/(r)) *sqrt(2*hbarc/(pi*q*r)) * sin(q*r/hbarc) *calcFF_2D(qx,qy) *sqrt(q) * 1/hbarc * 1/hbarc * 1/sqrt(hbarc) *4 * pi * pi * sqrt(2 * pi/q) * sqrt(2*hbarc/(pi*q*r)) * sqrt(hbarc) * sqrt(r) * sin(q *r/hbarc);
    //cout << "r: " << r << endl;
    //return 4 * pi * sqrt(2 * pi/q) * sqrt(2*hbarc/(pi*q*r)) * sqrt(hbarc) * sqrt(r) * sin(q *r/hbarc);
}

double F1_integral_dr_2D_test(double qx, double qy)
{
    double r_min = 0.;
    double r_max = 30;
    TF2 *f = new TF2("f",F1_integrand_dr_2D_test,r_min,r_max,r_min,r_max,2);  
    f->SetParameters(qx,qy);
    return f->Integral(r_min,r_max,r_min,r_max,1e-6); // [-]
}

// F1(q) == Ftrue
void G1_transform_to_F1_2D() 
{
    int bins = 1000.;
    double step_r = .015;
    double step = 0.0005;
    double min = 0.0005;
    double max = 0.5;
    TH2D *hankelTransform = new TH2D("hankelTransform","Transform Woods-Saxon: G(r) #rightarrow F(q)",bins,min,max,bins,min,max);
    hankelTransform->Sumw2();

    for(int i=0;i<bins;i++)
    {
        double qx = (i+1)*step;
        //double x = (i+1)*step_r;
        for (int j =0; j<bins; j++)
        {
            double qy = (j+1)*step;
            //double y = (j+1)*step_r;
            //double a = G1_integral_dq_2D(x,y);
            double b = F1_integral_dr_2D_test(qx,qy);
            //double val = a*b;
            hankelTransform->SetBinContent(i+1,j+1,b);
        }
        cout << "qx: " << qx << endl;
    }
    TFile *FF_WS_Transformations = new TFile("FF_WS_Transformations_Gr_to_Fq_2D.root","recreate");
	hankelTransform->SetTitle("Transformed WS: G(r) #rightarrow F(q)");
    hankelTransform->GetYaxis()->SetTitle("F(q)");
	hankelTransform->GetXaxis()->SetTitle("q [GeV]");
	hankelTransform->Write();
	FF_WS_Transformations->Close();

    TCanvas *c1 = new TCanvas("c1", "Hankel Transform", 800, 600);
    //hankelTransform->ProjectionY()->Draw();
    hankelTransform->Draw();
}
*/

//***Transforming Woods-Saxon: Gtrue(r)->F2(q)->G2(r)***//
// Gtrue(r) -> F2(q) (uses real Woods-Saxon)      ****F2 should == Ftrue
/*
double F2_integrand_dr_2D(double *x, double *par) 
{
    double xx = x[0];
    double yy = x[1];
    double r = sqrt(xx*xx + yy*yy);
    double qx = par[0];
    double qy = par[1];
    double q = sqrt(qx*qx + qy*qy);
    return 4 * pi * pi * sqrt(2 * pi/q) * sqrt(2*hbarc/(pi*q*r)) * sin(q *r/hbarc) * calcWS_2D(xx,yy) * sqrt(hbarc) * sqrt(r); //[1/fm]
    //double q = sqrt(qx*qx + qy*qy);
    //return 4 * pi * pi * sqrt(2 * pi/q) * sqrt(2*hbarc/(pi*q*r)) * calcWS_2D(xx,yy) * sqrt(hbarc) * xx*yy/(xx*xx+yy*yy) * sqrt(r) * sin(q *r/hbarc); //[1/fm]
}

double F2_integral_dr_2D(double qx, double qy)
{
    double r_min = 0.;
    double r_max = 30;
    TF2 *f = new TF2("f",F2_integrand_dr_2D,r_min,r_max,r_min,r_max,2);  
    f->SetParameters(qx,qy);
    return f->Integral(r_min,r_max,r_min,r_max,1e-6);// [-]
}

// F2(q)
void G_transform_to_F2_2D()
{
    int bins = 1000.;
    double step = 0.0001;
    double min = 0.0001;
    double max = 0.1;

    TH2D *hankelTransformWS = new TH2D("hankelTransformWS","Transform Woods-Saxon",bins,min,max,bins,min,max);
    hankelTransformWS->Sumw2();
    for(int i=0;i<bins;i++)
    {
        double tx = (i+1)*step;
        double qx = sqrt(tx);
        //double qx = (i+1)*step;
        for(int j=0; j<bins; j++)
        {
            double ty = (j+1)*step;
            double qy = sqrt(ty);
            //double qy = (j+1)*step;
	        double a = F2_integral_dr_2D(qx,qy);
            hankelTransformWS->SetBinContent(i+1,j+1,a);
        }
        cout << "qx:" << qx << endl;
    }
    TFile *eAu_5_41_dsigmadydt = new TFile("transform_WS_to_F_of_t_2D_test.root","recreate");
    //TFile *eAu_5_41_dsigmadydt = new TFile("transform_WS_to_F_of_q_2D.root","recreate");
	hankelTransformWS->SetTitle("Transformed WS: #rho(r) #rightarrow F(q)");
    hankelTransformWS->GetYaxis()->SetTitle("F(q)");
	hankelTransformWS->GetXaxis()->SetTitle("q [GeV]");
	hankelTransformWS->Write();
	eAu_5_41_dsigmadydt->Close();

    TCanvas *c1 = new TCanvas("c1", "Hankel Transform", 800, 600);
    hankelTransformWS->Draw();
    //hankelTransformWS->ProjectionY()->Draw();
}

// F2(q) -> G2(r)   **G2 should == Gtrue
double G2_integrand_dq_2D(double *x, double *par)
{
    double qx = x[0];
    double qy = x[1];
    double q = sqrt(qx*qx + qy*qy);
    double xx = par[0];
    double y = par[1];
    double r = sqrt(xx*xx + y*y);
    return 4 * pi * pi * sqrt(2*pi/r) * sqrt(2*hbarc/(pi*q*r)) * sin(q *r/hbarc) * F2_integral_dr_2D(qx,qy) * sqrt(q) * 1/hbarc * 1/hbarc * 1/sqrt(hbarc);  // [1/fm^3GeV]
}

double G2_integral_dq_2D(double x, double y)
{
    double q_min = 0.0005;
    double q_max = 1;
    TF2 *f = new TF2("f",G2_integrand_dq_2D,q_min,q_max,q_min,q_max,2);
    f->SetParameters(x,y);
    return f->Integral(q_min,q_max,q_min,q_max,1e-6); // [1/fm^3]
}

double G2_integrand_dq_2D_test(double *x, double *par)
{
    double qx = x[0];
    double qy = x[1];
    double q = sqrt(qx*qx + qy*qy);
    double xx = par[0];
    double y = par[1];
    double r = sqrt(xx*xx + y*y);
    return 2 * pi * sqrt(2*pi/r) * sqrt(2*hbarc/(pi*q*r)) * sin(q*r/hbarc) * sqrt(q) * 1/hbarc * 1/hbarc * 1/sqrt(hbarc);  // [1/fm^3GeV]
}

double G2_integral_dq_2D_test(double x, double y)
{
    double q_min = 0.0005;
    double q_max = 1;
    TF2 *f = new TF2("f",G2_integrand_dq_2D_test,q_min,q_max,q_min,q_max,2);
    f->SetParameters(x,y);
    return f->Integral(q_min,q_max,q_min,q_max,1e-3); // [1/fm^3]
}

//G2(r)
void F2_transform_to_G2_2D()
{
    int bins = 1000.;
    double step = .015;
    double min = 0.;
    double max = 15;

    TH2D *hankelTransformWS = new TH2D("hankelTransformWS","Transform Woods-Saxon",bins,min,max,bins,min,max);
    hankelTransformWS->Sumw2();

     // Gtrue(r)->F1(q)
    TFile *f3 = new TFile("transform_WS_to_F_of_q_2D.root");
    TH2D *get3_3 = (TH2D *)f3->Get("hankelTransformWS");
    double inte = get3_3->Integral();

    for(int i=0;i<bins;i++)
    {
        double x = (i+1)*step;
        for(int j=0; j<bins; j++)
        {
            double y = (j+1)*step;
	        double a = G2_integral_dq_2D_test(x,y)*inte;
            hankelTransformWS->SetBinContent(i+1,j+1,a);
            cout << "y: " << y << endl;
        }
        cout << "x: " << x << endl;
    }
    TFile *eAu_5_41_dsigmadydt = new TFile("transform_transform_WS_of_q_to_WS_of_r_2D.root","recreate");
	hankelTransformWS->SetTitle("Transformed-Transformed WS: F(q) #rightarrow #rho(r)");
    hankelTransformWS->GetYaxis()->SetTitle("G(r) [fm^{-3}]");
	hankelTransformWS->GetXaxis()->SetTitle("r [fm]");
	hankelTransformWS->Write();
	eAu_5_41_dsigmadydt->Close();

    TCanvas *c1 = new TCanvas("c1", "Hankel Transform", 800, 600);
    hankelTransformWS->Draw();
    //hankelTransformWS->ProjectionY()->Draw();
}
*/

//***Compare all 3 G(r)***//
/*
void stack_all_Gr_2D()
{
    // Real Woods-Saxon: rho(r)
    TFile *f1 = new TFile("FF_WS_Transformations_real_WS_2D.root");
    TH2D *get1_1 = (TH2D *)f1->Get("ws");
    TH1D *get1 = (TH1D *)get1_1->ProjectionY();
    get1->SetTitle("Woods-Saxon G(r)");
    get1->GetYaxis()->SetTitle("G(r) [fm^{-3}]");
    get1->GetXaxis()->SetTitle("r [fm]");
    get1->GetYaxis()->SetRangeUser(0,0.012);
    get1->Scale(197./get1->Integral(), "width");
    get1->SetLineColor(kBlack);
    
    // Ftrue(q)->G1(r)
    TFile *f2 = new TFile("FF_WS_Transformations_Fq_to_Gr_2D.root");
    TH2D *get2_2 = (TH2D *)f2->Get("hankelTransformFF");
    TH1D *get2 = (TH1D *)get2_2->ProjectionY();
    get2->Scale(197./get2->Integral(), "width");
    get2->SetLineColor(kRed);

    // Gtrue(r)->F2(q)->G2(r)
    TFile *f3 = new TFile("transform_transform_WS_of_q_to_WS_of_r_2D.root");
    TH2D *get3_3 = (TH2D *)f3->Get("hankelTransformWS");
    TH1D *get3 = (TH1D *)get3_3->ProjectionY();
    get3->Scale(197./get3->Integral(), "width");
    get3->SetLineColor(kMagenta);
    get3->SetLineStyle(2);

    TCanvas *c_f1 = new TCanvas();
	get1->Draw();
    get2->Draw("same");
    //get3->Draw("same");

    auto legend = new TLegend(0.43,0.78,0.63,0.9);
	legend->AddEntry(get1,"Real Woods-Saxon: G(r)","l");
    legend->AddEntry(get2,"F(q) #rightarrow G(r)","l");
    //legend->AddEntry(get3,"G(r) #rightarrow F(q) #rightarrow G(r)","l");
    legend->Draw();
}
*/
//***Compare all 3 F(q)***//
/*
void stack_all_Fq_2D()
{
    // Real Form Factor: Ftrue(q)
    TFile *f1 = new TFile("FF_WS_Transformations_real_FF_2D.root");
    TH2D *get1_1 = (TH2D *)f1->Get("FF");
    TH1D *get1 = (TH1D *)get1_1->ProjectionY();
    get1->SetTitle("Form Factor F(q)");
    get1->GetYaxis()->SetTitle("F(q)");
    get1->GetXaxis()->SetTitle("q [GeV]");
    get1->Scale(197./get1->Integral(), "width");
    //get1->GetYaxis()->SetRangeUser(1e-9,15);
    get1->SetLineColor(kBlack);

    // Gtrue(r)->F1(q)
    TFile *f3 = new TFile("transform_WS_to_F_of_q_2D.root");
    TH2D *get3_3 = (TH2D *)f3->Get("hankelTransformWS");
    TH1D *get3 = (TH1D *)get3_3->ProjectionY();
    get3->Scale(197./get3->Integral(), "width");
    get3->SetLineColor(kRed);

    // Ftrue(q)->G1(r)->F2(q)
    TFile *f2 = new TFile("FF_WS_Transformations_Gr_to_Fq_2D.root");
    TH2D *get2_2 = (TH2D *)f2->Get("hankelTransform");
    TH1D *get2 = (TH1D *)get2_2->ProjectionY();
    get2->Scale(197./get2->Integral(), "width");
    get2->SetLineColor(kMagenta);
    get2->SetLineStyle(2);

    TCanvas *c_f1 = new TCanvas();
	get1->Draw();
    get3->Draw("same");
    //get2->Draw("same");

    auto legend = new TLegend(0.43,0.78,0.63,0.9);
	legend->AddEntry(get1,"Real F(q)","l");
    legend->AddEntry(get3,"G(r) #rightarrow F(q)","l");
    //legend->AddEntry(get2,"F(q) #rightarrow G(r) #rightarrow F(q)","l");
    legend->Draw();
}

*/






//=========================================== 2D of t = (tx+ty) ========================================//


double calcWS_2D(double x, double y)
{
    double Vo = 2.12;  
    double A = 197.;
    double a = depth_Au;
    double R = 6.38;
    double r = sqrt(x*x + y*y);
    return Vo/(1. + exp((r-R)/a)); // [1/fm^3]
}

void getWS_2D()
{
    double bins = 1000;
    double step = .015;
    double min = 0.;
    double max = 15;
    TH2D *ws = new TH2D("ws","Woods-Saxon",bins,min,max,bins,min,max);
    ws->Sumw2();
    for(int i=0;i<bins;i++)
    {
        double x = (i+1)*step;
        for(int j=0; j<bins; j++)
        {
            double y = (j+1)*step;
            double a = calcWS_2D(x,y);
            ws->SetBinContent(i+1,j+1,a);
        }
        cout << "r: " << x << endl;
    }
	TFile *FF_WS_Transformations = new TFile("FF_WS_Transformations_real_WS_2D.root","recreate");
	ws->SetTitle("Woods-Saxon: G(r) [fm^{-3}]");
    ws->GetXaxis()->SetTitle("x [fm]");
    ws->GetYaxis()->SetTitle("y [fm]");
	ws->Write();
	FF_WS_Transformations->Close();

    TCanvas *c1 = new TCanvas("c1", "True Woods-Saxon", 800, 600);
    //ws->ProjectionY()->Draw();
    ws->Draw();
}

double calcFF_2D(double qx, double qy)
{
    double q = sqrt(qx*qx + qy*qy);
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
	    return sph * 1/(1. + (a0 * a0 * q*q/(hbarc*hbarc)));  // [-]
    }
}

void tgetFF_2D()
{
    double bins = 1000.;
    double step = 0.0001;
    double min = 0.0001;
    double max = 0.1;
    TH2D *FF = new TH2D("FF","Form Factor",bins,min,max,bins,min,max);
    FF->Sumw2();

    for(int i=0;i<bins;i++)
    {
        double tx = (i+1)*step;
        double qx = sqrt(tx);
        for(int j=0; j<bins; j++)
        {
            double ty = (j+1)*step;
            double qy = sqrt(ty);
            double a = calcFF_2D(qx,qy);
            FF->SetBinContent(i+1,j+1,a);
        }
        cout << "tx: " << tx << endl;
    }
    TFile *FF_WS_Transformations = new TFile("FF_WS_Transformations_real_FF_2D_test.root","recreate");
	FF->SetTitle("Form Factor: F(t)");
    FF->GetYaxis()->SetTitle("ty [GeV^{2}]");
	FF->GetXaxis()->SetTitle("tx [GeV^{2}]");
	FF->Write();
	FF_WS_Transformations->Close();

    TCanvas *c1 = new TCanvas("c1", "True Form Factor", 800, 600);
    FF->ProjectionY()->Draw();
    //FF->Draw();
}


                //================================== Transformations =====================================//

//***Transforming Form Factor: Ftrue(t)->G1(r)->F1(t)***//
// Ftrue(t) -> G1(r) (uses real Form Factor)    **** G1 should == Gtrue
double tG1_integrand_dq_2D(double *x, double *par) 
{
    double qx = x[0];
    double qy = x[1];
    double xx = par[0];
    double y = par[1];
    double r = sqrt(xx*xx + y*y);
    double q = sqrt(qx*qx + qy*qy);
    return  4 * pi * pi * sqrt(2 * pi/(r)) * sqrt(2*hbarc/(pi*q*r)) * sin(q *r/hbarc) * calcFF_2D(qx,qy) *sqrt(q) * 1/hbarc * 1/hbarc * 1/sqrt(hbarc);
    //return  4 * pi * pi * sqrt(2 * pi/(r)) *sqrt(2*hbarc/(pi*q*r)) * sin(q *r/hbarc) *tcalcFF_2D(qx,qy) * qx*qy/(qx*qx+qy*qy) *sqrt(q) * 1/hbarc * 1/hbarc * 1/sqrt(hbarc); //[1/fm^3Gev] 
}

double tG1_integral_dq_2D(double x, double y)
{
    double r = sqrt(x*x + y*y);
    double q_min = 0.0005;
    double q_max = 1;
    TF2 *f = new TF2("f",tG1_integrand_dq_2D,q_min,q_max,q_min,q_max,2);  
    f->SetParameters(x,y);
    return f->Integral(q_min,q_max,q_min,q_max,1e-3); //[1/fm^3]
}

// G1(r)
void tF_transform_to_G1_2D()
{
    int bins = 1000.;
    double step = .015;
    double min = 0.;
    double max = 15;

    TH2D *hankelTransformFF = new TH2D("hankelTransformFF","Transform Form Factor",bins,min,max,bins,min,max);
    hankelTransformFF->Sumw2();
    for(int i=0;i<bins;i++)
    {
        double x = (i+1)*step;
        for (int j =0; j<bins; j++)
        {
            double y = (j+1)*step;
            double a = tG1_integral_dq_2D(x,y);
            hankelTransformFF->SetBinContent(i+1,j+1,a);
        }
        cout << "x: " << x << endl;
    }
    TFile *eAu_5_41_dsigmadydt = new TFile("FF_WS_Transformations_Fq_to_Gr_2D_test.root","recreate");
	hankelTransformFF->SetTitle("Transformed WS: F(tx,ty) #rightarrow G(r)");
    hankelTransformFF->GetYaxis()->SetTitle("y [fm]");
	hankelTransformFF->GetXaxis()->SetTitle("x [fm]");
	hankelTransformFF->Write();
	eAu_5_41_dsigmadydt->Close();

    TCanvas *c1 = new TCanvas("c1", "Hankel Transform X-proj", 800, 600);
    //hankelTransformFF->ProjectionX()->Draw();

    //TCanvas *c2 = new TCanvas("c2", "Hankel Transform Y-proj", 800, 600);
    //hankelTransformFF->ProjectionY()->Draw();
    hankelTransformFF->Draw();
}


//G1(r) -> F1(t)        ****F1 should == Ftrue
double tF1_integrand_dr_2D(double *x, double *par) 
{
    double xx = x[0];
    double yy = x[1];
    double r = sqrt(xx*xx + yy*yy);
    double qx = par[0];
    double qy = par[1];
    double q = sqrt(qx*qx + qy*qy);
    //cout << "r: " << r << endl;
    //return 4 * pi * pi * sqrt(2 * pi/q) * sqrt(2*hbarc/(pi*q*r)) * sqrt(hbarc) * tG1_integral_dq_2D(xx,yy) * xx*yy/(xx*xx+yy*yy) * sqrt(r) * sin(q *r/hbarc);
    return 4 * pi * pi * sqrt(2*pi/q) * sqrt(2*hbarc/(pi*q*r)) * sin(q *r/hbarc) * sqrt(hbarc) * tG1_integral_dq_2D(xx,yy) * sqrt(r);
}

double tF1_integral_dr_2D(double qx, double qy)
{
    double r_min = 0.;
    double r_max = 30;
    TF2 *f = new TF2("f",tF1_integrand_dr_2D,r_min,r_max,r_min,r_max,2);  
    f->SetParameters(qx,qy);
    cout << "hi" << endl;
    return f->Integral(r_min,r_max,r_min,r_max,1e-3); // [-]
}

double tF1_integrand_dr_2D_test(double *x, double *par) 
{
    double xx = x[0];
    double yy = x[1];
    double r = sqrt(xx*xx + yy*yy);
    double qx = par[0];
    double qy = par[1];
    double q = sqrt(qx*qx + qy*qy);
    //return  4 * pi * pi * sqrt(2 * pi/(r)) *sqrt(2*hbarc/(pi*q*r)) * sin(q *r/hbarc) *calcFF_2D(qx,qy) * qx*qy/(qx*qx+qy*qy) *sqrt(q) * 1/hbarc * 1/hbarc * 1/sqrt(hbarc); //[1/fm^3Gev] 
    return  4 * pi * pi * sqrt(2 * pi/q) * sqrt(2*hbarc/(pi*q*r)) * sin(q *r/hbarc) * calcWS_2D(xx,yy) * sqrt(hbarc) * sqrt(r)*4 * pi * pi * sqrt(2 * pi/(r)) *sqrt(2*hbarc/(pi*q*r)) * sin(q*r/hbarc) *calcFF_2D(qx,qy) *sqrt(q) * 1/hbarc * 1/hbarc * 1/sqrt(hbarc) *4 * pi * pi * sqrt(2 * pi/q) * sqrt(2*hbarc/(pi*q*r)) * sqrt(hbarc) * sqrt(r) * sin(q *r/hbarc);
    //cout << "r: " << r << endl;
    //return 4 * pi * sqrt(2 * pi/q) * sqrt(2*hbarc/(pi*q*r)) * sqrt(hbarc) * sqrt(r) * sin(q *r/hbarc);
}

double tF1_integral_dr_2D_test(double qx, double qy)
{
    double r_min = 0.;
    double r_max = 30;
    TF2 *f = new TF2("f",tF1_integrand_dr_2D_test,r_min,r_max,r_min,r_max,2);  
    f->SetParameters(qx,qy);
    return f->Integral(r_min,r_max,r_min,r_max,1e-3); // [-]
}

// F1(t) == Ftrue
void tG1_transform_to_F1_2D() 
{
    int bins = 1000.;
    double step = 0.0001;
    double min = 0.0001;
    double max = 0.1;
    TH2D *hankelTransform = new TH2D("hankelTransform","Transform Woods-Saxon: G(r) #rightarrow F(t)",bins,min,max,bins,min,max);
    hankelTransform->Sumw2();

    for(int i=0;i<bins;i++)
    {
        double tx = (i+1)*step;
        double qx = sqrt(tx);
        for (int j =0; j<bins; j++)
        {
            double ty = (j+1)*step;
            double qy = sqrt(ty);
            //double a = tG1_integral_dq_2D(x,y);
            double b = tF1_integral_dr_2D_test(qx,qy);
            hankelTransform->SetBinContent(i+1,j+1,b);
        }
        cout << "tx: " << tx << endl;
    }
    TFile *FF_WS_Transformations = new TFile("FF_WS_Transformations_Gr_to_Fq_2D_test.root","recreate");
	hankelTransform->SetTitle("Transformed WS: G(r) #rightarrow F(t)");
    hankelTransform->GetYaxis()->SetTitle("ty [GeV^{2}]");
	hankelTransform->GetXaxis()->SetTitle("tx [GeV^{2}]");
	hankelTransform->Write();
	FF_WS_Transformations->Close();

    TCanvas *c1 = new TCanvas("c1", "Hankel Transform", 800, 600);
    //hankelTransform->ProjectionY()->Draw();
    hankelTransform->Draw();
}


//***Transforming Woods-Saxon: Gtrue(r)->F2(t)->G2(r)***//
// Gtrue(r) -> F2(t) (uses real Woods-Saxon)      ****F2 should == Ftrue
double tF2_integrand_dr_2D(double *x, double *par) 
{
    double xx = x[0];
    double yy = x[1];
    double r = sqrt(xx*xx + yy*yy);
    double qx = par[0];
    double qy = par[1];
    double q = sqrt(qx*qx + qy*qy);
    return 4 * pi * pi * sqrt(2 * pi/q) * sqrt(2*hbarc/(pi*q*r)) * sin(q *r/hbarc) * calcWS_2D(xx,yy) * sqrt(hbarc) * sqrt(r); //[1/fm]
    //double q = sqrt(qx*qx + qy*qy);
    //return 4 * pi * pi * sqrt(2 * pi/q) * sqrt(2*hbarc/(pi*q*r)) * calcWS_2D(xx,yy) * sqrt(hbarc) * xx*yy/(xx*xx+yy*yy) * sqrt(r) * sin(q *r/hbarc); //[1/fm]
}

double tF2_integral_dr_2D(double qx, double qy)
{
    double r_min = 0.;
    double r_max = 30;
    TF2 *f = new TF2("f",tF2_integrand_dr_2D,r_min,r_max,r_min,r_max,2);  
    f->SetParameters(qx,qy);
    return f->Integral(r_min,r_max,r_min,r_max,1e-3);// [-]
}

// F2(t)
void tG_transform_to_F2_2D()
{
    int bins = 1000.;
    double step = 0.0001;
    double min = 0.0001;
    double max = 0.1;

    TH2D *hankelTransformWS = new TH2D("hankelTransformWS","Transform Woods-Saxon",bins,min,max,bins,min,max);
    hankelTransformWS->Sumw2();
    for(int i=0;i<bins;i++)
    {
        double tx = (i+1)*step;
        double qx = sqrt(tx);
        //double qx = (i+1)*step;
        for(int j=0; j<bins; j++)
        {
            double ty = (j+1)*step;
            double qy = sqrt(ty);
            //double qy = (j+1)*step;
	        double a = tF2_integral_dr_2D(qx,qy);
            hankelTransformWS->SetBinContent(i+1,j+1,a);
        }
        cout << "tx:" << tx << endl;
    }
    TFile *eAu_5_41_dsigmadydt = new TFile("transform_WS_to_F_of_q_2D_test.root","recreate");
	hankelTransformWS->SetTitle("Transformed WS: G_{true} (r) #rightarrow F(t)");
    hankelTransformWS->GetYaxis()->SetTitle("ty [GeV^{2}]");
	hankelTransformWS->GetXaxis()->SetTitle("tx [GeV^{2}]");
	hankelTransformWS->Write();
	eAu_5_41_dsigmadydt->Close();

    TCanvas *c1 = new TCanvas("c1", "Hankel Transform", 800, 600);
    hankelTransformWS->Draw();
    //hankelTransformWS->ProjectionY()->Draw();
}


// F2(t) -> G2(r)   **G2 should == Gtrue
double tG2_integrand_dq_2D(double *x, double *par)
{
    double qx = x[0];
    double qy = x[1];
    double q = sqrt(qx*qx + qy*qy);
    double xx = par[0];
    double y = par[1];
    double r = sqrt(xx*xx + y*y);
    return 4 * pi * pi * sqrt(2*pi/r) * sqrt(2*hbarc/(pi*q*r)) * sin(q *r/hbarc) * tF2_integral_dr_2D(qx,qy) * sqrt(q) * 1/hbarc * 1/hbarc * 1/sqrt(hbarc);  // [1/fm^3GeV]
}

double tG2_integral_dq_2D(double x, double y)
{
    double q_min = 0.0005;
    double q_max = 1;
    TF2 *f = new TF2("f",tG2_integrand_dq_2D,q_min,q_max,q_min,q_max,2);
    f->SetParameters(x,y);
    return f->Integral(q_min,q_max,q_min,q_max,1e-3); // [1/fm^3]
}

double tG2_integrand_dq_2D_test(double *x, double *par)
{
    double qx = x[0];
    double qy = x[1];
    double q = sqrt(qx*qx + qy*qy);
    double xx = par[0];
    double y = par[1];
    double r = sqrt(xx*xx + y*y);
    return 4 * pi * pi * sqrt(2 * pi/(r)) * sqrt(2*hbarc/(pi*q*r)) * sin(q *r/hbarc) * calcFF_2D(qx,qy) *sqrt(q) * 1/hbarc * 1/hbarc * 1/sqrt(hbarc)*2 * pi * sqrt(2*pi/r) * sqrt(2*hbarc/(pi*q*r)) * sin(q*r/hbarc) * sqrt(q) * 1/hbarc * 1/hbarc * 1/sqrt(hbarc);  // [1/fm^3GeV]
}

double tG2_integral_dq_2D_test(double x, double y)
{
    double q_min = 0.0005;
    double q_max = 1;
    TF2 *f = new TF2("f",tG2_integrand_dq_2D_test,q_min,q_max,q_min,q_max,2);
    f->SetParameters(x,y);
    return f->Integral(q_min,q_max,q_min,q_max,1e-3); // [1/fm^3]
}

//G2(r)
void tF2_transform_to_G2_2D()
{
    int bins = 1000.;
    double step = .015;
    double min = 0.;
    double max = 15;

    TH2D *hankelTransformWS = new TH2D("hankelTransformWS","Transform Woods-Saxon",bins,min,max,bins,min,max);
    hankelTransformWS->Sumw2();

    for(int i=0;i<bins;i++)
    {
        double x = (i+1)*step;
        for(int j=0; j<bins; j++)
        {
            double y = (j+1)*step;
	        double a = tG2_integral_dq_2D_test(x,y);
            hankelTransformWS->SetBinContent(i+1,j+1,a);
            cout << "y: " << y << endl;
        }
        cout << "x: " << x << endl;
    }
    TFile *eAu_5_41_dsigmadydt = new TFile("transform_transform_WS_of_q_to_WS_of_r_2D_test.root","recreate");
	hankelTransformWS->SetTitle("Transformed-Transformed WS: F(t) #rightarrow G(r)");
    hankelTransformWS->GetYaxis()->SetTitle("y [fm]");
	hankelTransformWS->GetXaxis()->SetTitle("x [fm]");
	hankelTransformWS->Write();
	eAu_5_41_dsigmadydt->Close();

    TCanvas *c1 = new TCanvas("c1", "Hankel Transform", 800, 600);
    hankelTransformWS->Draw();
    //hankelTransformWS->ProjectionY()->Draw();
}


//***Compare all 3 G(r)***//

void tstack_all_Gr_2D()
{
    // Real Woods-Saxon: Gtrue(r)
    TFile *f1 = new TFile("FF_WS_Transformations_real_WS_2D.root");
    TH2D *get1_1 = (TH2D *)f1->Get("ws");
    TH1D *get1 = (TH1D *)get1_1->ProjectionY();
    get1->SetTitle("Woods-Saxon: G(r) [fm^{-3}]");
    get1->GetYaxis()->SetTitle("y [fm]");
    get1->GetXaxis()->SetTitle("x [fm]");
    get1->GetYaxis()->SetRangeUser(0,0.012);
    get1->Scale(197./get1->Integral(), "width");
    get1->SetLineColor(kBlack);
    
    // Ftrue(t)->G1(r)
    TFile *f2 = new TFile("FF_WS_Transformations_Fq_to_Gr_2D_test.root");
    TH2D *get2_2 = (TH2D *)f2->Get("hankelTransformFF");
    TH1D *get2 = (TH1D *)get2_2->ProjectionY();
    get2->Scale(197./get2->Integral(), "width");
    get2->SetLineColor(kRed);
/*
    // Gtrue(r)->F2(t)->G2(r)
    TFile *f3 = new TFile("transform_transform_WS_of_t_to_WS_of_r_2D_test.root");
    TH2D *get3_3 = (TH2D *)f3->Get("hankelTransformWS");
    TH1D *get3 = (TH1D *)get3_3->ProjectionY();
    get3->Scale(197./get3->Integral(), "width");
    get3->SetLineColor(kMagenta);
    get3->SetLineStyle(2);
*/
    TCanvas *c_f1 = new TCanvas();
	get1->Draw();
    get2->Draw("same");
    //get3->Draw("same");

    auto legend = new TLegend(0.43,0.78,0.63,0.9);
	legend->AddEntry(get1,"G_{true} (r)","l");
    legend->AddEntry(get2,"F_{true} (t) #rightarrow G(r)","l");
    //legend->AddEntry(get3,"G_{true} (r) #rightarrow F(t) #rightarrow G(r)","l");
    legend->Draw();
}

//***Compare all 3 F(q)***//

void tstack_all_Ft_2D()
{
    // Real Form Factor: Ftrue(t)
    TFile *f1 = new TFile("FF_WS_Transformations_real_FF_2D_test.root");
    TH2D *get1_1 = (TH2D *)f1->Get("FF");
    TH1D *get1 = (TH1D *)get1_1->ProjectionY();
    get1->SetTitle("Form Factor F(t)");
    get1->GetYaxis()->SetTitle("ty [GeV^{2}]");
    get1->GetXaxis()->SetTitle("tx [GeV^{2}]");
    get1->Scale(197./get1->Integral(), "width");
    //get1->GetYaxis()->SetRangeUser(1e-9,15);
    get1->SetLineColor(kBlack);

    // Gtrue(r)->F1(t)
    TFile *f3 = new TFile("transform_WS_to_F_of_q_2D_test.root");
    TH2D *get3_3 = (TH2D *)f3->Get("hankelTransformWS");
    TH1D *get3 = (TH1D *)get3_3->ProjectionY();
    get3->Scale(197./get3->Integral(), "width");
    get3->SetLineColor(kRed);

    // Ftrue(t)->G1(r)->F2(t)
    TFile *f2 = new TFile("FF_WS_Transformations_Gr_to_Fq_2D_test.root");
    TH2D *get2_2 = (TH2D *)f2->Get("hankelTransform");
    TH1D *get2 = (TH1D *)get2_2->ProjectionY();
    get2->Scale(197./get2->Integral(), "width");
    get2->SetLineColor(kMagenta);
    get2->SetLineStyle(2);

    TCanvas *c_f1 = new TCanvas();
	get1->Draw();
    get3->Draw("same");
    get2->Draw("same");

    auto legend = new TLegend(0.43,0.78,0.63,0.9);
	legend->AddEntry(get1,"F_{true} (t)","l");
    legend->AddEntry(get3,"G_{true} (r) #rightarrow F(t)","l");
    legend->AddEntry(get2,"F_{true} (t) #rightarrow G(r) #rightarrow F(t)","l");
    legend->Draw();
}
















 //======================= Original Form Factor: F(t) =====================================//

double calcFF_oft(double t)
{
	if(t==0){return 0;}
    else{
        double q = sqrt(t);
        double A = 197;
        double ro = 1.25;
        double rho0 = 2.12; // 0.181
        double R = 6.38; //ro*5.81865;  // A^(1/3) = 5.81865 >> R = 7.27331 
	    const double arg1 = q * R / hbarc;  
	    const double arg2 = hbarc /q;  // RWS_Au = 6.38 ( from .h file)
	    const double sph  = (sin(arg1) - arg1 * cos(arg1)) * 4*pi*rho0 * arg2 * arg2 *arg2/ double(A); 
	    const double a0   = .70;
	    return sph * 1/(1. + (a0 * a0 * q*q/(hbarc*hbarc)));  // [-]
    }
}

void getFF_oft()
{
    double bins = 1000.;
    double step = 0.0001;
    double min = 0.0001;
    double max = 0.1;
    TH1D *FF = new TH1D("FF","Form Factor",bins,min,max);
    FF->Sumw2();

    for(int i=0;i<bins;i++)
    {
        double t = (i+1)*step;
        double a = calcFF_oft(t)*calcFF_oft(t);
        FF->SetBinContent(i+1,a);
        cout << "t: " << t << endl;
    }
	TFile *FF_WS_Transformations = new TFile("FF_WS_Transformations_real_FF_t_v2.root","recreate");
	FF->SetTitle("Form Factor: F(t)");
    FF->GetYaxis()->SetTitle("F(t)");
	FF->GetXaxis()->SetTitle("t [GeV^2]");
	FF->Write();
	FF_WS_Transformations->Close();

    TCanvas *c1 = new TCanvas("c1", "True Form Factor", 800, 600);
    FF->Draw();
    c1->Update();
}


                        //================================== Transformations =====================================//

//***Transforming Form Factor: Ftrue(t)->G1(r)->F1(t)***//
// Ftrue(t) -> G1(r) (uses real Form Factor)    **** G1 should == Gtrue
double G1_integrand_dq_t(double *x, double *par) 
{
    double q = x[0];
    double r = par[0];
    //double q = sqrt(t);
    return sqrt(2*pi/r)*sqrt(2*hbarc/(pi*q*r))*sin(q*r/hbarc)*calcFF_oft(q)*sqrt(q)*1/hbarc*1/hbarc*1/sqrt(hbarc);
    //return  2 * pi * sqrt(2 * pi/(r)) *sqrt(2*hbarc/(pi*q*r)) * sin(q*r/hbarc) *calcFF_oft(t) * 1/2 * 1/sqrt(t) * q *sqrt(q) * 1/hbarc * 1/hbarc * 1/sqrt(hbarc) ; //[1/fm^3Gev^2] 
}

double G1_integral_dq_t(double r)
{
    double t_min = 0.0005;
    double t_max = 5;
    TF1 *f = new TF1("f",G1_integrand_dq_t,t_min,t_max,1);  
    f->SetParameters(r);
    return f->Integral(t_min,t_max,1e-6); //[1/fm^3]
}

// G1(r)
void F_transform_to_G1_t()
{
    int bins = 1000.;
    double step = 0.015;
    double min = 0.;
    double max = 15;

    TH1D *hankelTransformFF = new TH1D("hankelTransformFF","Transform Form Factor",bins,min,max);
    hankelTransformFF->Sumw2();
    for(int i=0;i<bins;i++)
    {
        double r = (i+1)*step;
	    double a = G1_integral_dq_t(r);
        hankelTransformFF->SetBinContent(i+1,a);
        cout << "r: " << r << endl;
    }
    TFile *eAu_5_41_dsigmadydt = new TFile("FF_WS_Transformations_Ft_to_Gr_v2.root","recreate");
    //hankelTransformFF->Scale(197/hankelTransformFF->Integral());
	hankelTransformFF->SetTitle("Transformed WS: F(t) #rightarrow G(r)");
    hankelTransformFF->GetYaxis()->SetTitle("G(r) [fm^{-3}]");
	hankelTransformFF->GetXaxis()->SetTitle("r [fm]");
	hankelTransformFF->Write();
	eAu_5_41_dsigmadydt->Close();

    TCanvas *c1 = new TCanvas("c1", "Hankel Transform", 800, 600);
    hankelTransformFF->Draw();
    c1->Update();
}


//G1(r) -> F1(t)        ****F1 should == Ftrue
double F1_integrand_dr_t(double *x, double *par) 
{
    double r = x[0];
    double t = par[0];
    return 2 * pi * sqrt(2 * pi/sqrt(t)) * sqrt(2*hbarc/(pi*sqrt(t)*r)) * sqrt(hbarc) * G1_integral_dq_t(r) * r * sqrt(r) * sin(sqrt(t) *r/hbarc);
}

double F1_integral_dr_t(double t)
{
    double r_min = 0.;
    double r_max = 30;
    TF1 *f = new TF1("f",F1_integrand_dr_t,r_min,r_max,1);  
    f->SetParameters(t);
    return f->Integral(r_min,r_max,1e-6); // [-]
}

// F1(t) == Ftrue
void G1_transform_to_F1_t() 
{
    int bins = 1000.;
    double step = 0.0001;
    double min = 0.0001;
    double max = 0.1;
    TH1D *hankelTransform = new TH1D("hankelTransform","Transform Woods-Saxon: G(r) #rightarrow F(t)",bins,min,max);
    hankelTransform->Sumw2();
    int b[bins];
    for(int i=0;i<bins;i++)
    {
        double t = (i+1)*step;
	    double a = F1_integral_dr_t(t);
        b[i] = a;
        hankelTransform->SetBinContent(i+1,a);
        cout << "t:" << t << endl;
    }
    TFile *FF_WS_Transformations = new TFile("FF_WS_Transformations_Gr_to_Ft.root","recreate");
	hankelTransform->SetTitle("Transformed WS: G(r) #rightarrow F(t)");
    hankelTransform->GetYaxis()->SetTitle("F(t)");
	hankelTransform->GetXaxis()->SetTitle("t [GeV^2]");
	hankelTransform->Write();
	FF_WS_Transformations->Close();

    TCanvas *c1 = new TCanvas("c1", "Hankel Transform", 800, 600);
    hankelTransform->Draw();
    c1->Update();
}


//***Transforming Woods-Saxon: Gtrue(r)->F2(t)->G2(r)***//
// Gtrue(r) -> F2(t) (uses real Woods-Saxon)      ****F2 should == Ftrue

double F2_integrand_dr_t(double *x, double *par) 
{
    double r = x[0];
    double t = par[0];
    return 2 * pi * sqrt(2 * pi/sqrt(t)) * sqrt(2*hbarc/(pi*sqrt(t)*r)) * calcWS(r) * sqrt(hbarc) * r * sqrt(r) * sin(sqrt(t) *r/hbarc); //[1/fm]
}

double F2_integral_dr_t(double t)
{
    double r_min = 0.;
    double r_max = 30;
    TF1 *f = new TF1("f",F2_integrand_dr_t,r_min,r_max,1);  
    f->SetParameters(t);
    return f->Integral(r_min,r_max,1e-6);// [-]
}

// F2(t)
void G_transform_to_F2_t()
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
        double t = (i+1)*step;
	    double a = F2_integral_dr_t(t);
        b[i] = a;
        hankelTransformWS->SetBinContent(i+1,a);
        cout << "t:" << t << endl;
    }
    TFile *eAu_5_41_dsigmadydt = new TFile("transform_WS_to_F_of_t.root","recreate");
	hankelTransformWS->SetTitle("Transformed WS: #rho(r) #rightarrow F(t)");
    hankelTransformWS->GetYaxis()->SetTitle("F(t)");
	hankelTransformWS->GetXaxis()->SetTitle("t [GeV^2]");
	hankelTransformWS->Write();
	eAu_5_41_dsigmadydt->Close();

    TCanvas *c1 = new TCanvas("c1", "Hankel Transform", 800, 600);
    hankelTransformWS->Draw();
    c1->Update();
}


// F2(t) -> G2(r)   **G2 should == Gtrue
double G2_integrand_dq_t(double *x, double *par)
{
    double t = x[0];
    double r = par[0];
    double q = sqrt(t);
    return 2 * pi * sqrt(2*pi/r) * sqrt(2*hbarc/(pi*q*r)) * sin(q *r/hbarc) * F2_integral_dr_t(t) * 1/2 * 1/sqrt(t) * q * sqrt(q) * 1/hbarc * 1/hbarc * 1/sqrt(hbarc);  // [1/fm^3GeV]
}

double G2_integral_dq_t(double r)
{
    double t_min = 0.0005;
    double t_max = 5;
    TF1 *f = new TF1("f",G2_integrand_dq_t,t_min,t_max,1);
    f->SetParameters(r);
    return f->Integral(t_min,t_max,1e-6); // [1/fm^3]
}

//G2(r)
void F2_transform_to_G2_t()
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
	    double a = G2_integral_dq_t(r);
        b[i] = a;
        hankelTransformWS->SetBinContent(i+1,a);
        cout << "r: " << r << endl;
    }
    TFile *eAu_5_41_dsigmadydt = new TFile("transform_transform_WS_of_t_to_WS_of_r.root","recreate");
	hankelTransformWS->SetTitle("Transformed-Transformed WS: F(t) #rightarrow #rho(r)");
    hankelTransformWS->GetYaxis()->SetTitle("G(r) [fm^{-3}]");
	hankelTransformWS->GetXaxis()->SetTitle("r [fm]");
	hankelTransformWS->Write();
	eAu_5_41_dsigmadydt->Close();

    TCanvas *c1 = new TCanvas("c1", "Hankel Transform", 800, 600);
    hankelTransformWS->Draw();
    c1->Update();
}


//***Compare all 3 G(r)***//

void stack_all_Gr_t()
{
    // Real Woods-Saxon: rho(r)
    TFile *f1 = new TFile("FF_WS_Transformations_real_WS.root");
    TH1D *get1 = (TH1D *)f1->Get("ws");
    get1->SetTitle("Woods-Saxon G(r)");
    get1->GetYaxis()->SetTitle("G(r) [fm^{-3}]");
    get1->GetXaxis()->SetTitle("r [fm]");
    get1->GetYaxis()->SetRangeUser(0,0.012);
    get1->Scale(197./get1->Integral(), "width");
    get1->SetLineColor(kBlack);
    //TH1*h0 = (TH1*)(get1->Clone("h0"));
    //h0->Scale(1./h0->Integral(), "width");
    //h0->SetLineColor(kBlack);

    // Ftrue(t)->G1(r)
    TFile *f2 = new TFile("FF_WS_Transformations_Ft_to_Gr.root");
    TH1D *get2 = (TH1D *)f2->Get("hankelTransformFF");
    get2->Scale(197./get2->Integral(), "width");
    get2->SetLineColor(kRed);
    //TH1*h1 = (TH1*)(get2->Clone("h1"));
    //h1->Scale(1./h1->Integral(), "width");
    //h1->SetLineColor(kRed);

    // Gtrue(r)->F2(t)->G2(r)
    TFile *f3 = new TFile("transform_transform_WS_of_t_to_WS_of_r.root");
    TH1D *get3 = (TH1D *)f3->Get("hankelTransformWS");
    get3->Scale(197./get3->Integral(), "width");
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
    legend->AddEntry(get2,"F(t) #rightarrow G(r)","l");
    legend->AddEntry(get3,"G(r) #rightarrow F(t) #rightarrow G(r)","l");
    //legend->AddEntry(h1,"F(t) #rightarrow G(r)","l");
    //legend->AddEntry(h2,"G(r) #rightarrow F(t) #rightarrow G(r)","l");
    legend->Draw();
}


//***Compare all 3 F(t)***//

void stack_all_Ft()
{
    // Real Form Factor: Ftrue(t)
    TFile *f1 = new TFile("FF_WS_Transformations_real_FF_t.root");
    TH1D *get1 = (TH1D *)f1->Get("FF");
    get1->SetTitle("Form Factor F(t)");
    get1->GetYaxis()->SetTitle("F(t)");
    get1->GetXaxis()->SetTitle("t [GeV^2]");
    get1->Scale(197./get1->Integral(), "width");
    //get1->GetYaxis()->SetRangeUser(1e-9,15);
    get1->SetLineColor(kBlack);

    // Gtrue(r)->F1(t)
    TFile *f3 = new TFile("transform_WS_to_F_of_t.root");
    TH1D *get3 = (TH1D *)f3->Get("hankelTransformWS");
    get3->Scale(197./get3->Integral(), "width");
    get3->SetLineColor(kRed);
    //TH1*h1 = (TH1*)(get3->Clone("h1"));
    //h1->Scale(1./h1->Integral(), "width");
    //h1->SetLineColor(kRed);

    // Ftrue(t)->G1(r)->F2(t)
    TFile *f2 = new TFile("FF_WS_Transformations_Gr_to_Ft.root");
    TH1D *get2 = (TH1D *)f2->Get("hankelTransform");
    get2->Scale(197./get2->Integral(), "width");
    get2->SetLineColor(kMagenta);
    get2->SetLineStyle(2);
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
	legend->AddEntry(get1,"Real F(t)","l");
    legend->AddEntry(get3,"G(r) #rightarrow F(t)","l");
    legend->AddEntry(get2,"F(t) #rightarrow G(r) #rightarrow F(t)","l");
    //legend->AddEntry(h1,"G(r) #rightarrow F(t)","l");
    //legend->AddEntry(h2,"F(t) #rightarrow G(r) #rightarrow F(t)","l");
    legend->Draw();
}

void stack_all_Ft_v2()
{
    // Real Form Factor: Ftrue(t)
    TFile *f1 = new TFile("FF_WS_Transformations_real_FF_t_v2.root");
    TH1D *get1 = (TH1D *)f1->Get("FF");
    get1->SetTitle("Form Factor F(t)");
    get1->GetYaxis()->SetTitle("F(t)");
    get1->GetXaxis()->SetTitle("t [GeV^2]");
    get1->Scale(197./get1->Integral(), "width");
    //get1->GetYaxis()->SetRangeUser(1e-9,15);
    get1->SetLineColor(kBlack);

    // Gtrue(r)->F1(t)
    TFile *f3 = new TFile("transform_WS_to_F_of_t.root");
    TH1D *get3 = (TH1D *)f3->Get("hankelTransformWS");
    get3->Scale(197./get3->Integral(), "width");
    get3->SetLineColor(kRed);
    //TH1*h1 = (TH1*)(get3->Clone("h1"));
    //h1->Scale(1./h1->Integral(), "width");
    //h1->SetLineColor(kRed);

    // Ftrue(t)->G1(r)->F2(t)
    TFile *f2 = new TFile("FF_WS_Transformations_Gr_to_Ft.root");
    TH1D *get2 = (TH1D *)f2->Get("hankelTransform");
    get2->Scale(197./get2->Integral(), "width");
    get2->SetLineColor(kMagenta);
    get2->SetLineStyle(2);
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
	legend->AddEntry(get1,"Real F(t)","l");
    legend->AddEntry(get3,"G(r) #rightarrow F(t)","l");
    legend->AddEntry(get2,"F(t) #rightarrow G(r) #rightarrow F(t)","l");
    //legend->AddEntry(h1,"G(r) #rightarrow F(t)","l");
    //legend->AddEntry(h2,"F(t) #rightarrow G(r) #rightarrow F(t)","l");
    legend->Draw();
}


//=========================================== 2D ========================================//


double calcFF_2D_oft(double tx, double ty)
{
    double q = sqrt(tx + ty);
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
	    return sph * 1/(1. + (a0 * a0 * q*q/(hbarc*hbarc)));  // [-]
    }
}

void getFF_2D_oft()
{
    double bins = 1000.;
    double step = 0.0001;
    double min = 0.0001;
    double max = 0.1;
    TH2D *FF = new TH2D("FF","Form Factor",bins,min,max,bins,min,max);
    FF->Sumw2();

    for(int i=0;i<bins;i++)
    {
        double tx = (i+1)*step;
        for(int j=0; j<bins; j++)
        {
            double ty = (j+1)*step;
            double a = calcFF_2D_oft(tx,ty);
            FF->SetBinContent(i+1,j+1,a);
        }
        cout << "tx: " << tx << endl;
    }
	TFile *FF_WS_Transformations = new TFile("FF_WS_Transformations_real_FF_2D_t.root","recreate");
	FF->SetTitle("Form Factor: F(t)");
    FF->GetYaxis()->SetTitle("t_{y} [GeV^{2}]");
	FF->GetXaxis()->SetTitle("t_{x} [GeV^{2}]");
	FF->Write();
	FF_WS_Transformations->Close();

    TCanvas *c1 = new TCanvas("c1", "True Form Factor", 800, 600);
    //FF->ProjectionY()->Draw();
    FF->Draw();
}


                        //================================== Transformations =====================================//

//***Transforming Form Factor: Ftrue(tx,ty)->G1(x,y)->F1(tx,ty)***//
// Ftrue(tx,ty) -> G1(x,y) (uses real Form Factor)    **** G1 should == Gtrue
double G1_integrand_dq_2D_oft(double *x, double *par) 
{
    double qx = x[0];
    double qy = x[1];
    double xx = par[0];
    double y = par[1];
    double r = sqrt(xx*xx + y*y);
    double q = sqrt(qx*qx + qy*qy);
    return  4 * pi * pi * sqrt(2 * pi/(r)) * sqrt(2*hbarc/(pi*q*r)) * sin(q *r/hbarc) * calcFF_2D(qx,qy) *sqrt(q) * 1/hbarc * 1/hbarc * 1/sqrt(hbarc);
    //return  4 * pi * pi * sqrt(2 * pi/(r)) *sqrt(2*hbarc/(pi*q*r)) * sin(q *r/hbarc) *calcFF_2D(qx,qy) * qx*qy/(qx*qx+qy*qy) *sqrt(q) * 1/hbarc * 1/hbarc * 1/sqrt(hbarc); //[1/fm^3Gev] 
/*
    double tx = x[0];
    double ty = x[1];
    double xx = par[0];
    double y = par[1];
    double r = sqrt(xx*xx + y*y);
    double qx = sqrt(tx);
    double qy = sqrt(ty);
    double q = sqrt(tx + ty);
    return  pi * pi * sqrt(2 * pi/(r)) *sqrt(2*hbarc/(pi*q*r)) * sin(q *r/hbarc) *calcFF_2D_oft(tx,ty) * 1/(tx+ty) *sqrt(q) * 1/hbarc * 1/hbarc * 1/sqrt(hbarc); //[1/fm^3Gev]
    */ 
}

double G1_integral_dq_2D_oft(double x, double y)
{
    double r = sqrt(x*x + y*y);
    double q_min = 0.0005;
    double q_max = 1;
    TF2 *f = new TF2("f",G1_integrand_dq_2D_oft,q_min,q_max,q_min,q_max,2);  
    f->SetParameters(x,y);
    return f->Integral(q_min,q_max,q_min,q_max,1e-6); //[1/fm^3]
/*
    double r = sqrt(x*x + y*y);
    double t_min = 0.0005;
    double t_max = 1;
    TF2 *f = new TF2("f",G1_integrand_dq_2D_oft,t_min,t_max,t_min,t_max,2);  
    f->SetParameters(x,y);
    return f->Integral(t_min,t_max,t_min,t_max,1e-3); //[1/fm^3]
    */
}

// G1(x,y)
void F_transform_to_G1_2D_oft()
{
    int bins = 1000.;
    double step = .015;
    double min = 0.;
    double max = 15;

    TH2D *hankelTransformFF = new TH2D("hankelTransformFF","Transform Form Factor",bins,min,max,bins,min,max);
    hankelTransformFF->Sumw2();
    for(int i=0;i<bins;i++)
    {
        double x = (i+1)*step;
        for (int j =0; j<bins; j++)
        {
            double y = (j+1)*step;
            double a = G1_integral_dq_2D_oft(x,y);
            hankelTransformFF->SetBinContent(i+1,j+1,a);
        }
        cout << "x: " << x << endl;
    }
    TFile *eAu_5_41_dsigmadydt = new TFile("FF_WS_Transformations_Ft_to_Gr_2D.root","recreate");
	hankelTransformFF->SetTitle("Transformed WS: F(tx,ty) #rightarrow G(x,y)");
    hankelTransformFF->GetYaxis()->SetTitle("y [fm]");
	hankelTransformFF->GetXaxis()->SetTitle("x [fm]");
	hankelTransformFF->Write();
	eAu_5_41_dsigmadydt->Close();

    TCanvas *c1 = new TCanvas("c1", "Hankel Transform", 800, 600);
    //hankelTransformFF->ProjectionX()->Draw();

    //TCanvas *c2 = new TCanvas("c2", "Hankel Transform Y-proj", 800, 600);
    //hankelTransformFF->ProjectionY()->Draw();
    hankelTransformFF->Draw();
}


//G1(x,y) -> F1(tx,ty)        ****F1 should == Ftrue
double F1_integrand_dr_2D_oft(double *x, double *par) 
{
    double xx = x[0];
    double yy = x[1];
    double r = sqrt(xx*xx + yy*yy);
    double tx = par[0];
    double ty = par[1];
    double qx = sqrt(tx);
    double qy = sqrt(ty);
    double q = sqrt(qx*qx + qy*qy);
    //cout << "r: " << r << endl;
    return 4 * pi * pi * sqrt(2 * pi/q) * sqrt(2*hbarc/(pi*q*r)) * sin(q *r/hbarc) * sqrt(hbarc) * G1_integral_dq_2D_oft(xx,yy) * sqrt(r);
}

double F1_integral_dr_2D_oft(double tx, double ty)
{
    double r_min = 0.;
    double r_max = 30;
    TF2 *f = new TF2("f",F1_integrand_dr_2D_oft,r_min,r_max,r_min,r_max,2);  
    f->SetParameters(tx,ty);
    cout << "hi" << endl;
    return f->Integral(r_min,r_max,r_min,r_max,1e-3); // [-]
}

double F1_integrand_dr_2D_test_oft(double *x, double *par) 
{
    double xx = x[0];
    double yy = x[1];
    double r = sqrt(xx*xx + yy*yy);
    double tx = par[0];
    double ty = par[1];
    double qx = sqrt(tx);
    double qy = sqrt(ty);
    double q = sqrt(qx*qx + qy*qy);
    //cout << "r: " << r << endl;
    return 4 * pi * sqrt(2 * pi/q) * sqrt(2*hbarc/(pi*q*r)) * sqrt(hbarc) * sqrt(r) * sin(q *r/hbarc);
}

double F1_integral_dr_2D_test_oft(double tx, double ty)
{
    double r_min = 0.;
    double r_max = 30;
    TF2 *f = new TF2("f",F1_integrand_dr_2D_test_oft,r_min,r_max,r_min,r_max,2);  
    f->SetParameters(tx,ty);
    return f->Integral(r_min,r_max,r_min,r_max,1e-3); // [-]
}

// F1(tx,ty) == Ftrue
void G1_transform_to_F1_2D_oft() 
{
    int bins = 1000.;
    double step = 0.0001;
    double min = 0.0001;
    double max = 0.1;
    TH2D *hankelTransform = new TH2D("hankelTransform","Transform Woods-Saxon: G(x,y) #rightarrow F(tx,ty)",bins,min,max,bins,min,max);
    hankelTransform->Sumw2();

     // Ftrue(tx,ty)->G1(x,y)
    TFile *f2 = new TFile("FF_WS_Transformations_Ft_to_Gr_2D.root");
    TH2D *get2_2 = (TH2D *)f2->Get("hankelTransformFF");
    double inte = get2_2->Integral();

    for(int i=0;i<bins;i++)
    {
        double tx = (i+1)*step;
        for(int j=0; j<bins; j++)
        {
            double ty = (j+1)*step;
	        double a = F1_integral_dr_2D_oft(tx,ty);
            hankelTransform->SetBinContent(i+1,j+1,a);
        }
        cout << "tx:" << tx << endl;
    }
    TFile *FF_WS_Transformations = new TFile("FF_WS_Transformations_Gr_to_Ft_2D.root","recreate");
	hankelTransform->SetTitle("Transformed WS: G(x,y) #rightarrow F(tx,ty)");
    hankelTransform->GetYaxis()->SetTitle("F(tx,ty)");
	hankelTransform->GetXaxis()->SetTitle("t [GeV^2]");
	hankelTransform->Write();
	FF_WS_Transformations->Close();

    TCanvas *c1 = new TCanvas("c1", "Hankel Transform", 800, 600);
    //hankelTransform->ProjectionY()->Draw();
    hankelTransform->Draw();
}


//***Transforming Woods-Saxon: Gtrue(x,y)->F2(tx,ty)->G2(x,y)***//
// Gtrue(x,y) -> F2(tx,ty) (uses real Woods-Saxon)      ****F2 should == Ftrue
double F2_integrand_dr_2D_oft(double *x, double *par) 
{
    double xx = x[0];
    double yy = x[1];
    double r = sqrt(xx*xx + yy*yy);
    double tx = par[0];
    double ty = par[1];
    double qx = sqrt(tx);
    double qy = sqrt(ty);
    double q = sqrt(qx*qx + qy*qy);
    return 4 * pi * pi * sqrt(2 * pi/q) * sqrt(2*hbarc/(pi*q*r)) * sin(q *r/hbarc) * calcWS_2D(xx,yy) * sqrt(hbarc) * sqrt(r); //[1/fm]
}

double F2_integral_dr_2D_oft(double tx, double ty)
{
    double r_min = 0.;
    double r_max = 30;
    TF2 *f = new TF2("f",F2_integrand_dr_2D_oft,r_min,r_max,r_min,r_max,2);  
    f->SetParameters(tx,ty);
    return f->Integral(r_min,r_max,r_min,r_max,1e-3);// [-]
}

// F2(tx,ty)
void G_transform_to_F2_2D_oft()
{
    int bins = 1000.;
    double step = 0.0001;
    double min = 0.0001;
    double max = 0.1;

    TH2D *hankelTransformWS = new TH2D("hankelTransformWS","Transform Woods-Saxon",bins,min,max,bins,min,max);
    hankelTransformWS->Sumw2();
    for(int i=0;i<bins;i++)
    {
        double tx = (i+1)*step;
        for(int j=0; j<bins; j++)
        {
            double ty = (j+1)*step;
	        double a = F2_integral_dr_2D_oft(tx,ty);
            hankelTransformWS->SetBinContent(i+1,j+1,a);
        }
        cout << "tx:" << tx << endl;
    }
    TFile *eAu_5_41_dsigmadydt = new TFile("transform_WS_to_F_of_t_2D.root","recreate");
	hankelTransformWS->SetTitle("Transformed WS: #rho(x,y) #rightarrow F(tx,ty)");
    hankelTransformWS->GetYaxis()->SetTitle("F(tx,ty)");
	hankelTransformWS->GetXaxis()->SetTitle("tx [GeV^2]");
	hankelTransformWS->Write();
	eAu_5_41_dsigmadydt->Close();

    TCanvas *c1 = new TCanvas("c1", "Hankel Transform", 800, 600);
    hankelTransformWS->Draw();
    //hankelTransformWS->ProjectionY()->Draw();
}


// F2(tx,ty) -> G2(x,y)   **G2 should == Gtrue
double G2_integrand_dq_2D_oft(double *x, double *par)
{
    double tx = x[0];
    double ty = x[1];
    double qx = sqrt(tx);
    double qy = sqrt(ty);
    double q = sqrt(qx*qx + qy*qy);
    double xx = par[0];
    double y = par[1];
    double r = sqrt(xx*xx + y*y);
    return pi * pi * sqrt(2*pi/r) * sqrt(2*hbarc/(pi*q*r)) * sin(q *r/hbarc) * F2_integral_dr_2D_oft(tx,ty) * 1/sqrt(tx*ty) * sqrt(q) * 1/hbarc * 1/hbarc * 1/sqrt(hbarc);  // [1/fm^3GeV]
}

double G2_integral_dq_2D_oft(double x, double y)
{
    double t_min = 0.0005;
    double t_max = 1;
    TF2 *f = new TF2("f",G2_integrand_dq_2D_oft,t_min,t_max,t_min,t_max,2);
    f->SetParameters(x,y);
    return f->Integral(t_min,t_max,t_min,t_max,1e-3); // [1/fm^3]
}

double G2_integrand_dq_2D_test_oft(double *x, double *par)
{
    double tx = x[0];
    double ty = x[1];
    double qx = sqrt(tx);
    double qy = sqrt(ty);
    double q = sqrt(qx*qx + qy*qy);
    double xx = par[0];
    double y = par[1];
    double r = sqrt(xx*xx + y*y);
    return 2 * pi * sqrt(2*pi/r) * sqrt(2*hbarc/(pi*q*r)) * sin(q*r/hbarc) * sqrt(q) * 1/hbarc * 1/hbarc * 1/sqrt(hbarc);  // [1/fm^3GeV]
}

double G2_integral_dq_2D_test_oft(double x, double y)
{
    double t_min = 0.0005;
    double t_max = 1;
    TF2 *f = new TF2("f",G2_integrand_dq_2D_test_oft,t_min,t_max,t_min,t_max,2);
    f->SetParameters(x,y);
    return f->Integral(t_min,t_max,t_min,t_max,1e-3); // [1/fm^3]
}

//G2(x,y)
void F2_transform_to_G2_2D_oft()
{
    int bins = 1000.;
    double step = .015;
    double min = 0.;
    double max = 15;

    TH2D *hankelTransformWS = new TH2D("hankelTransformWS","Transform Woods-Saxon",bins,min,max,bins,min,max);
    hankelTransformWS->Sumw2();

     // Gtrue(x,y)->F1(tx,ty)
    TFile *f3 = new TFile("transform_WS_to_F_of_t_2D.root");
    TH2D *get3_3 = (TH2D *)f3->Get("hankelTransformWS");
    double inte = get3_3->Integral();

    for(int i=0;i<bins;i++)
    {
        double x = (i+1)*step;
        for(int j=0; j<bins; j++)
        {
            double y = (j+1)*step;
	        double a = G2_integral_dq_2D_test_oft(x,y)*inte;
            hankelTransformWS->SetBinContent(i+1,j+1,a);
            cout << "y: " << y << endl;
        }
        cout << "x: " << x << endl;
    }
    TFile *eAu_5_41_dsigmadydt = new TFile("transform_transform_WS_of_t_to_WS_of_r_2D.root","recreate");
	hankelTransformWS->SetTitle("Transformed-Transformed WS: F(t) #rightarrow #rho(r)");
    hankelTransformWS->GetYaxis()->SetTitle("G(r) [fm^{-3}]");
	hankelTransformWS->GetXaxis()->SetTitle("r [fm]");
	hankelTransformWS->Write();
	eAu_5_41_dsigmadydt->Close();

    TCanvas *c1 = new TCanvas("c1", "Hankel Transform", 800, 600);
    hankelTransformWS->Draw();
    //hankelTransformWS->ProjectionY()->Draw();
}


//***Compare all 3 G(r)***//

void stack_all_Gr_2D_oft()
{
    // Real Woods-Saxon: rho(r)
    TFile *f1 = new TFile("FF_WS_Transformations_real_WS_2D.root");
    TH2D *get1_1 = (TH2D *)f1->Get("ws");
    TH1D *get1 = (TH1D *)get1_1->ProjectionY();
    get1->SetTitle("Woods-Saxon G(x,y)");
    get1->GetYaxis()->SetTitle("G(r) [fm^{-3}]");
    get1->GetXaxis()->SetTitle("r [fm]");
    get1->Scale(197./get1->Integral(), "width");
    get1->GetYaxis()->SetRangeUser(0,70);
    get1->SetLineColor(kBlack);
   
    // Ftrue(q)->G1(r)
    TFile *f2 = new TFile("FF_WS_Transformations_Ft_to_Gr_2D.root");
    TH2D *get2_2 = (TH2D *)f2->Get("hankelTransformFF");
    TH1D *get2 = (TH1D *)get2_2->ProjectionY();
    get2->Scale(197./get2->Integral(), "width");
    get2->SetLineColor(kRed);
/*
    // Gtrue(r)->F2(q)->G2(r)
    TFile *f3 = new TFile("transform_transform_WS_of_q_to_WS_of_r_2D.root");
    TH2D *get3_3 = (TH2D *)f3->Get("hankelTransformWS");
    TH1D *get3 = (TH1D *)get3_3->ProjectionY();
    get3->Scale(197./get3->Integral(), "width");
    get3->SetLineColor(kMagenta);
    get3->SetLineStyle(2);
*/
    TCanvas *c_f1 = new TCanvas();
	get1->Draw();
    get2->Draw("same");
    //get3->Draw("same");

    auto legend = new TLegend(0.43,0.78,0.63,0.9);
	legend->AddEntry(get1,"Real Woods-Saxon: G(x,y)","l");
    legend->AddEntry(get2,"F(q) #rightarrow G(r)","l");
    //legend->AddEntry(get3,"G(r) #rightarrow F(q) #rightarrow G(r)","l");
    legend->Draw();
}

//***Compare all 3 F(q)***//
void stack_all_Ft_2D()
{
    // Real Form Factor: Ftrue(q)
    TFile *f1 = new TFile("FF_WS_Transformations_real_FF_2D_t.root");
    TH2D *get1_1 = (TH2D *)f1->Get("FF");
    TH1D *get1 = (TH1D *)get1_1->ProjectionY();
    get1->SetTitle("Form Factor F(q)");
    get1->GetYaxis()->SetTitle("F(q)");
    get1->GetXaxis()->SetTitle("q [GeV]");
    get1->Scale(197./get1->Integral(), "width");
    //get1->GetYaxis()->SetRangeUser(1e-9,15);
    get1->SetLineColor(kBlack);

    // Gtrue(r)->F1(q)
    TFile *f3 = new TFile("transform_WS_to_F_of_q_2D.root");
    TH2D *get3_3 = (TH2D *)f3->Get("hankelTransformWS");
    TH1D *get3 = (TH1D *)get3_3->ProjectionY();
    get3->Scale(197./get3->Integral(), "width");
    get3->SetLineColor(kRed);

    // Ftrue(q)->G1(r)->F2(q)
    TFile *f2 = new TFile("FF_WS_Transformations_Gr_to_Fq_2D.root");
    TH2D *get2_2 = (TH2D *)f2->Get("hankelTransform");
    TH1D *get2 = (TH1D *)get2_2->ProjectionY();
    get2->Scale(197./get2->Integral(), "width");
    get2->SetLineColor(kMagenta);
    get2->SetLineStyle(2);

    TCanvas *c_f1 = new TCanvas();
	get1->Draw();
    //get3->Draw("same");
    //get2->Draw("same");

    auto legend = new TLegend(0.43,0.78,0.63,0.9);
	legend->AddEntry(get1,"Real F(q)","l");
    //legend->AddEntry(get3,"G(r) #rightarrow F(q)","l");
    //legend->AddEntry(get2,"F(q) #rightarrow G(r) #rightarrow F(q)","l");
    legend->Draw();
}

