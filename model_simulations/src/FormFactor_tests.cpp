#include "FormFactor_t_1D.h"
#include "FormFactor_t_2D.h"
#include "FormFactor_q_1D.h"
#include "FormFactor_q_2D.h"
#include "FormFactor_Wedge.h"
#include "FormFactor_resolution.h"
#include "FormFactor_saturation_data.h"
#include "FF_WS_Transforms_1D.h"
#include "FF_WS_Transforms_2D.h"
#include "WoodsSaxon_1D.h"
#include "WoodsSaxon_2D.h"
#include "dsigma_resolution_add_wedge_1D.h"
#include "WedgeResolution.h"
#include "FormFactor_resolution_add_wedge_1D.h"
#include "FormFactor_resolution_add_wedge_1D2.h"
#include "FormFactor_resolution_add_wedge_2D.h"
//#include "FormFactor_transform_resolution_add_wedge_1D.h"
#include "FormFactor_transform_resolution_add_wedge_2D.h"
#include "FormFactor_q_2DwResCut.h"
//#include "RooUnfoldResponse.h"
//#include "RooUnfoldBayes.h"
//#include "RooUnfoldSvd.h"
//#include "RooUnfoldTUnfold.h"
//#include "RooUnfoldIds.h"
#include "TSVDUnfold.h"
#include <TMatrixD.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TPad.h>
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
#include "TLegendEntry.h"
#include "Math/IFunction.h"
#include <cmath>
#include "TSystem.h"
#include "TAxis.h"
#include "TPaveLabel.h"
#include "TFormula.h"


using namespace std;


// Woods-Saxon 1D tests
void WS_1D_func()
{
    double Vo = 2.12, R = 6.38, a = 0.535, r_min = 0, r_max = 10;
    double bins = 1000, t_min = 0, t_max = 0.25;
    WoodsSaxon_1D ws(Vo,R,a,r_min,r_max,bins,t_min,t_max);
    TF1 *test = ws.getWoodsSaxon1D_fun();
    test->Draw();
}

void WS_1D_hist()
{
    double Vo = 2.12, R = 6.38, a = 0.535, r_min = 0, r_max = 10;
    double bins = 1000, t_min = 0, t_max = 0.25;
    WoodsSaxon_1D ws(Vo,R,a,r_min,r_max,bins,t_min,t_max);
    TH1D *test = ws.getWoodsSaxon1D_hist();
    test->Draw();
}

void WS_transform_1D_func()
{
    double Vo = 2.12, R = 6.38, a = 0.535, r_min = 0, r_max = 10;
    double bins = 1000, t_min = 0.0001, t_max = 0.1;
    WoodsSaxon_1D ws(Vo,R,a,r_min,r_max,bins,t_min,t_max);
    //TF1 *test = ws.get_transformed_WoodsSaxon1D_fun();
    //test->Draw();
}

void WS_transform_1D_hist()
{
    double Vo = 2.12, R = 6.38, a = 0.535, r_min = 0, r_max = 10;
    double bins = 1000, t_min = 0.0001, t_max = 0.1, q_min = 0.01, q_max = 0.31623;
    WoodsSaxon_1D ws(Vo,R,a,r_min,r_max,bins,t_min,t_max);
    TH1D *test = ws.get_transformed_WoodsSaxon1D_hist();
    test->Draw();
}


// Woods-Saon 2D tests
void WS_2D_func()
{
    double Vo = 2.12, R =6.38, a = 0.535;
    double x_min = 0, x_max = 10, y_min = 0, y_max = 10, bins = 100;
    double tx_min = 0, tx_max = 0.1, ty_min = 0, ty_max = 0.1;
    WoodsSaxon_2D ws(Vo,R,a,x_min,x_max,y_min,y_max,bins,tx_min,tx_max,ty_min,ty_max);
    TF2 *test = ws.getWoodsSaxon2D_fun();
    test->Draw();
}

void WS_2D_hist()
{
    double Vo = 2.12, R =6.38, a = 0.535;
    double x_min = 0, x_max = 10, y_min = 0, y_max = 10, bins = 100;
    double tx_min = 0, tx_max = 0.1, ty_min = 0, ty_max = 0.1;
    WoodsSaxon_2D ws(Vo,R,a,x_min,x_max,y_min,y_max,bins,tx_min,tx_max,ty_min,ty_max);
    TH2D *test = ws.getWoodsSaxon2D_hist();
    test->ProjectionY()->Draw();
}

void WS_transform_2D_fun()
{
    double Vo = 2.12, R =6.38, a = 0.535;
    double x_min = 0, x_max = 10, y_min = 0, y_max = 10, bins = 100;
    double tx_min = 0, tx_max = 0.1, ty_min = 0, ty_max = 0.1;
    WoodsSaxon_2D ws(Vo,R,a,x_min,x_max,y_min,y_max,bins,tx_min,tx_max,ty_min,ty_max);
    //TF2 *test = ws.get_transformed_WoodsSaxon2D_fun();
    //test->Draw();
}

void WS_transform_2D_hist()
{
    double Vo = 2.12, R =6.38, a = 0.535;
    double x_min = 0, x_max = 10, y_min = 0, y_max = 10, bins = 100;
    double tx_min = 0, tx_max = 0.1, ty_min = 0, ty_max = 0.1;
    WoodsSaxon_2D ws(Vo,R,a,x_min,x_max,y_min,y_max,bins,tx_min,tx_max,ty_min,ty_max);
    TH2D *test = ws.get_transformed_WoodsSaxon2D_hist();
    test->ProjectionY()->Draw();
}


// FormFactor_q_1D tests
void FFq_squared_fun()
{
    double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7;
    double q_min = 0.001, q_max = 0.3;
    FormFactor_q_1D ff(A,Vo,R,a0,q_min,q_max);
    TF1 *test = ff.getFormFactorq_1D();
    test->Draw();
}


// FormFactor_t_1D tests
void FFt_squared_fun()
{
    double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7;
    double t_min = 0.0001, t_max = 0.1;
    double bins = 1000, r_min = 0, r_max = 15;
    FormFactor_t_1D ff(A,Vo,R,a0,t_min,t_max,bins,r_min,r_max);
    TF1 *test = ff.getFormFactort_1D();
    test->Draw();
}

void FFt_hist()
{
    double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7;
    double t_min = 0, t_max = 0.2;
    double bins = 1000, r_min = 0, r_max = 15;
    FormFactor_t_1D ff(A,Vo,R,a0,t_min,t_max,bins,r_min,r_max);
    TH1D *test = ff.getFF_hist();
    test->Draw();
}

void FFt_transform_fun()
{
    double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7;
    double t_min = 0.0001, t_max = 0.1;
    double bins = 1000, r_min = 0, r_max = 15;
    FormFactor_t_1D ff(A,Vo,R,a0,t_min,t_max,bins,r_min,r_max);
    TF1 *test = ff.getTransformedFF_TF1();
    test->Draw();
}

void FFt_transform_hist()
{
    double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7;
    double t_min = 0.0001, t_max = 0.1;
    double bins = 1000, r_min = 0, r_max = 15;
    FormFactor_t_1D ff(A,Vo,R,a0,t_min,t_max,bins,r_min,r_max);
    TH1D *test = ff.getTransformedFF_hist();
    test->Draw();
}


// FormFactor_q_2D tests
void FFq2_2D()
{
    double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7;
    double qy_min = 0, qy_max = 0.5, qx_min = 0, qx_max = 0.5, bins = 1000;
    double x_min= 0.015, x_max = 15, y_min = 0.015, y_max = 15;
    FormFactor_q_2D ff(A,Vo,R,a0,qy_min,qy_max,qx_min,qx_max,bins);
    TF2* test = ff.getFormFactorq2_2D();
    test->Draw("AP");
}

void FFq_2D_hist()
{
    double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7;
    double qy_min = 0, qy_max = 0.5, qx_min = 0, qx_max = 0.5, bins = 1000;
    FormFactor_q_2D ff(A,Vo,R,a0,qy_min,qy_max,qx_min,qx_max,bins);
    TH2D *ff_histz = ff.getFormFactorq_2D_hist();
    TCanvas *c1 = new TCanvas("c1", "c1", 1000, 1000);
    c1->SetLeftMargin(0.13); 
    c1->SetBottomMargin(0.13);

    ff_histz->SetTitle("|F(q_{x},q_{y})|^{2}");
    ff_histz->GetYaxis()->SetTitleSize(0.04);
    ff_histz->GetXaxis()->SetTitleSize(0.04);
    ff_histz->GetXaxis()->SetTitle("q_{x} = #sqrt{t_{x}} [GeV/c]");
    ff_histz->GetYaxis()->SetTitle("q_{y} = #sqrt{t_{y}} [GeV/c]");
    ff_histz->GetYaxis()->SetTitleOffset(1.4);
    ff_histz->GetXaxis()->SetTitleOffset(1.1);
    ff_histz->Draw("col");

    TGaxis *xAxis = new TGaxis(ff_histz->GetXaxis()->GetXmin(), 0, ff_histz->GetXaxis()->GetXmax(), 0, ff_histz->GetXaxis()->GetXmin(), ff_histz->GetXaxis()->GetXmax(), 510, "U");
    xAxis->SetLabelSize(0.02);
    xAxis->SetLabelOffset(0.01);
    xAxis->Draw();

TGaxis *yAxis = new TGaxis(0, ff_histz->GetYaxis()->GetXmin(), 0, ff_histz->GetYaxis()->GetXmax(), ff_histz->GetYaxis()->GetXmin(), ff_histz->GetYaxis()->GetXmax(), 510, "U");
    yAxis->SetLabelSize(0.01);
    yAxis->SetLabelOffset(0.02);
    yAxis->Draw();

double length = .5; 
double theta = pi/12; 
double x_end = length * sin(theta);
double y_end = length * cos(theta);

TArrow *arrow = new TArrow(0, 0, x_end, y_end, 0.02, "|>");
    arrow->SetLineColor(kBlack);
    arrow->SetLineWidth(2);
    arrow->Draw();

double theta2 = -pi/12; 
double x_end2 = length * sin(theta2);
double y_end2 = length * cos(theta2);

TArrow *arrow2 = new TArrow(0, 0, x_end2, y_end2, 0.02, "|>");
    arrow2->SetLineColor(kBlack);
    arrow2->SetLineWidth(2);
    arrow2->Draw();

double theta3 = -13*pi/12; 
double x_end3 = length * sin(theta3);
double y_end3 = length * cos(theta3);

TArrow *arrow3 = new TArrow(0, 0, x_end3, y_end3, 0.02, "|>");
    arrow3->SetLineColor(kBlack);
    arrow3->SetLineWidth(2);
    arrow3->Draw();

double theta4 = 13*pi/12; 
double x_end4 = length * sin(theta4);
double y_end4 = length * cos(theta4);

TArrow *arrow4 = new TArrow(0, 0, x_end4, y_end4, 0.02, "|>");
    arrow4->SetLineColor(kBlack);
    arrow4->SetLineWidth(2);
    arrow4->Draw();

    auto el3 = new TEllipse(0, 0, 0.3, 0.3, -15, 15);
   el3->SetFillStyle(0);
   el3->SetLineColor(kBlack);
   el3->SetLineWidth(2);
   el3->SetTheta(270);
   el3->Draw();

    auto el = new TEllipse(0, 0, -0.3, -0.3, -15, 15);
    el->SetFillStyle(0);
    el->SetLineColor(kBlack);
    el->SetLineWidth(2);
    el->SetTheta(270);
    el->Draw();

    gPad->SetLogz(1);
    gStyle->SetOptStat(0);
    c1->Draw();
}

void FFq_2D_hist_projection()
{
    double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7;
    double qy_min = 0, qy_max = 0.5, qx_min = 0, qx_max = 0.5, bins = 1000;
    FormFactor_q_2D ff(A,Vo,R,a0,qy_min,qy_max,qx_min,qx_max,bins);
    TH2D *test = ff.getFormFactorq_2D_hist();
    test->ProjectionY()->Draw();
}


// FormFactor_t_2D tests
void FFt2_2D()
{
    double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7;
    double ty_min = 0, ty_max = 0.25, tx_min = 0, tx_max = 0.25, bins = 1000;
    double x_min = 0, x_max = 15, y_min = 0, y_max = 15;
    FormFactor_t_2D ff(A,Vo,R,a0,tx_min,tx_max,ty_min,ty_max,bins,x_min,x_max,y_min,y_max);
    TF2 *test = ff.getFormFactort2_2D();
    test->Draw();
}

void FFt_2D_hist()
{
    double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7;
    double ty_min = 0, ty_max = 0.25, tx_min = 0, tx_max = 0.25, bins = 1000;
    double x_min = 0, x_max = 15, y_min = 0, y_max = 15;
    FormFactor_t_2D ff(A,Vo,R,a0,tx_min,tx_max,ty_min,ty_max,bins,x_min,x_max,y_min,y_max);
    TH2D *test = ff.getFormFactort_hist();
    test->Draw();
}

void FFt_2D_hist_projection()
{
    double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7;
    double ty_min = 0, ty_max = 0.25, tx_min = 0, tx_max = 0.25, bins = 1000;
    double x_min = 0, x_max = 15, y_min = 0, y_max = 15;
    FormFactor_t_2D ff(A,Vo,R,a0,tx_min,tx_max,ty_min,ty_max,bins,x_min,x_max,y_min,y_max);
    TH2D *test = ff.getFormFactort_hist();
    test->ProjectionY()->Draw();
}

void FFt_2D_transform_fun()
{
    double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7;
    double ty_min = 0, ty_max = 0.25, tx_min = 0, tx_max = 0.25, bins = 1000;
    double x_min = 0, x_max = 15, y_min = 0, y_max = 15;
    FormFactor_t_2D ff(A,Vo,R,a0,tx_min,tx_max,ty_min,ty_max,bins,x_min,x_max,y_min,y_max);
    TF2 *test = ff.getTransformed_TF2();
    test->Draw();
}

void FFt_transform_hist_2d()
{
    double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7;
    double ty_min = 0, ty_max = 0.25, tx_min = 0, tx_max = 0.25, bins = 1000;
    double x_min = 0, x_max = 15, y_min = 0, y_max = 15;
    FormFactor_t_2D ff(A,Vo,R,a0,tx_min,tx_max,ty_min,ty_max,bins,x_min,x_max,y_min,y_max);
    TH2D *test = ff.getTransformed_hist();
    test->Draw();
}

void FFt_transform_hist_proj()
{
    double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7;
    double ty_min = 0, ty_max = 0.25, tx_min = 0, tx_max = 0.25, bins = 1000;
    double x_min = 0, x_max = 15, y_min = 0, y_max = 15;
    FormFactor_t_2D ff(A,Vo,R,a0,tx_min,tx_max,ty_min,ty_max,bins,x_min,x_max,y_min,y_max);
    TH2D *test = ff.getTransformed_hist();
    test->ProjectionY()->Draw();
}


// Compare transforms
void trueFF_transformed_WS_1d()
{
    double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7, a = 0.535;
    double t_min = 0.0001, t_max = 0.1;
    double bins = 1000, r_min = 0, r_max = 15;

    FormFactor_t_1D ff(A,Vo,R,a0,t_min,t_max,bins,r_min,r_max);
    TF1 *test1 = ff.getFormFactort_1D();

    WoodsSaxon_1D ws(Vo,R,a,r_min,r_max,bins,t_min,t_max);
    TF1 *test2 = ws.getWoodsSaxon1D_fun();
    test2->SetLineStyle(3);
    test2->SetLineColor(kBlack);
    test2->Draw();

   
}

void trueFF_transformed_WS_1d_hist()
{
    double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7, a = 0.535;
    double t_min = 0.0001, t_max = 0.1;
    double bins = 1000, r_min = 0, r_max = 15;

    FormFactor_t_1D ff(A,Vo,R,a0,t_min,t_max,bins,r_min,r_max);
    TH1D *test1 = ff.getFF_hist();
    test1->SetLineStyle(2);
    test1->SetLineColor(kRed);
    test1->Draw();

    WoodsSaxon_1D ws(Vo,R,a,r_min,r_max,bins,t_min,t_max);
    TH1D *test2 = ws.get_transformed_WoodsSaxon1D_hist();
    test2->SetLineStyle(3);
    test2->SetLineColor(kBlack);
    test2->Draw("same");
}

void trueFF_transformed_WS_2d()
{
    double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7, a = 0.535;
    double ty_min = 0, ty_max = 0.25, tx_min = 0, tx_max = 0.25, bins = 1000;
    double x_min = 0, x_max = 15, y_min = 0, y_max = 15;
    double r_min = 0, r_max =15;


double tx_prime_min = 0, tx_prime_max = 0.25; 

double phi_min = 0, phi_max = pi/9, sigma = 0.1;


    FormFactor_t_2D ff(A,Vo,R,a0,tx_min,tx_max,ty_min,ty_max,bins,x_min,x_max,y_min,y_max);
    TF2 *test1 = ff.getFormFactort2_2D();
    test1->SetLineStyle(2);
    test1->SetLineColor(kRed);
    test1->Draw();

    WoodsSaxon_2D ws(Vo,R,a,x_min,x_max,y_min,y_max,bins,tx_min,tx_max,ty_min,ty_max);
    //TF2 *test2 = ws.get_transformed_WoodsSaxon2D_fun();
    //test2->SetLineStyle(3);
    //test2->SetLineColor(kBlack);
    //test2->Draw("same");

   

    
}

void trueFF_transformed_WS_2d_hist()
{
    double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7, a = 0.535;
    double ty_min = 0, ty_max = 0.25, tx_min = 0, tx_max = 0.25, bins = 1000;
    double x_min = 0, x_max = 15, y_min = 0, y_max = 15;

    FormFactor_t_2D ff(A,Vo,R,a0,tx_min,tx_max,ty_min,ty_max,bins,x_min,x_max,y_min,y_max);
    TH2D *test1 = ff.getFormFactort_hist();
    test1->SetLineStyle(2);
    test1->SetLineColor(kRed);
    test1->Draw();

    WoodsSaxon_2D ws(Vo,R,a,x_min,x_max,y_min,y_max,bins,tx_min,tx_max,ty_min,ty_max);
    TH2D *test2 = ws.get_transformed_WoodsSaxon2D_hist();
    test2->SetLineStyle(3);
    test2->SetLineColor(kBlack);
    test2->Draw("same");
}

void trueFF_transformed_WS_2d_hist_projection()
{
    double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7, a = 0.535;
    double ty_min = 0, ty_max = 0.25, tx_min = 0, tx_max = 0.25, bins = 1000;
    double x_min = 0, x_max = 15, y_min = 0, y_max = 15;

    FormFactor_t_2D ff(A,Vo,R,a0,tx_min,tx_max,ty_min,ty_max,bins,x_min,x_max,y_min,y_max);
    TH2D *test1 = ff.getFormFactort_hist();
    test1->SetLineStyle(2);
    test1->SetLineColor(kRed);
    test1->ProjectionY()->Draw();

    WoodsSaxon_2D ws(Vo,R,a,x_min,x_max,y_min,y_max,bins,tx_min,tx_max,ty_min,ty_max);
    TH2D *test2 = ws.get_transformed_WoodsSaxon2D_hist();
    test2->SetLineStyle(3);
    test2->SetLineColor(kBlack);
    test2->ProjectionY()->Draw("same");
}

void trueWS_transformed_FF_1d()
{
    double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7, a = 0.535;
    double t_min = 0.0001, t_max = 0.1;
    double bins = 1000, r_min = 0, r_max = 15;
    FormFactor_t_1D ff(A,Vo,R,a0,t_min,t_max,bins,r_min,r_max);
    TF1 *test1 = ff.getTransformedFF_TF1();
    test1->SetLineStyle(2);
    test1->SetLineColor(kRed);
    test1->Draw();

    WoodsSaxon_1D ws(Vo,R,a,r_min,r_max,bins,t_min,t_max);
    TF1 *test2 = ws.getWoodsSaxon1D_fun();
    test2->SetLineStyle(3);
    test2->SetLineColor(kBlack);
    test2->Draw("same");
}

void trueWS_transformed_FF_1d_hist()
{
    double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7, a = 0.535;
    double t_min = 0.0001, t_max = 0.1;
    double bins = 1000, r_min = 0, r_max = 15;

    FormFactor_t_1D ff(A,Vo,R,a0,t_min,t_max,bins,r_min,r_max);
    TH1D *test1 = ff.getTransformedFF_hist();
    test1->SetLineStyle(2);
    test1->SetLineColor(kRed);
    test1->Draw();

    WoodsSaxon_1D ws(Vo,R,a,r_min,r_max,bins,t_min,t_max);
    TH1D *test2 = ws.getWoodsSaxon1D_hist();
    test2->SetLineStyle(3);
    test2->SetLineColor(kBlack);
    test2->Draw("same");
}

void trueWS_transformed_FF_2d()
{
    double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7, a = 0.535;
    double ty_min = 0, ty_max = 0.25, tx_min = 0, tx_max = 0.25, bins = 1000;
    double x_min = 0, x_max = 15, y_min = 0, y_max = 15;


    FormFactor_t_2D ff(A,Vo,R,a0,tx_min,tx_max,ty_min,ty_max,bins,x_min,x_max,y_min,y_max);
    TF2 *test1 = ff.getTransformed_TF2();
    test1->SetLineStyle(2);
    test1->SetLineColor(kRed);
    test1->Draw();

    WoodsSaxon_2D ws(Vo,R,a,x_min,x_max,y_min,y_max,bins,tx_min,tx_max,ty_min,ty_max);
    TF2 *test2 = ws.getWoodsSaxon2D_fun();
    test2->SetLineStyle(3);
    test2->SetLineColor(kBlack);
    test2->Draw("same");
    
}

void trueWS_transformed_FF_2d_hist()
{
    double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7, a = 0.535;
    double ty_min = 0, ty_max = 0.25, tx_min = 0, tx_max = 0.25, bins = 1000;
    double x_min = 0, x_max = 15, y_min = 0, y_max = 15;

    FormFactor_t_2D ff(A,Vo,R,a0,tx_min,tx_max,ty_min,ty_max,bins,x_min,x_max,y_min,y_max);
    TH2D *test1 = ff.getTransformed_hist();
    test1->SetLineStyle(2);
    test1->SetLineColor(kRed);
    test1->Draw();

    WoodsSaxon_2D ws(Vo,R,a,x_min,x_max,y_min,y_max,bins,tx_min,tx_max,ty_min,ty_max);
    TH2D *test2 = ws.getWoodsSaxon2D_hist();
    test2->SetLineStyle(3);
    test2->SetLineColor(kBlack);
    test2->Draw("same");
    
}

void trueWS_transformed_FF_2d_hist_projection()
{
    double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7, a = 0.535;
    double ty_min = 0, ty_max = 0.25, tx_min = 0, tx_max = 0.25, bins = 1000;
    double x_min = 0, x_max = 15, y_min = 0, y_max = 15;

    FormFactor_t_2D ff(A,Vo,R,a0,tx_min,tx_max,ty_min,ty_max,bins,x_min,x_max,y_min,y_max);
    TH2D *test1 = ff.getTransformed_hist();
    test1->SetLineStyle(2);
    test1->SetLineColor(kRed);
    test1->ProjectionY()->Draw();

    WoodsSaxon_2D ws(Vo,R,a,x_min,x_max,y_min,y_max,bins,tx_min,tx_max,ty_min,ty_max);
    TH2D *test2 = ws.getWoodsSaxon2D_hist();
    test2->SetLineStyle(3);
    test2->SetLineColor(kBlack);
    test2->ProjectionY()->Draw("same");
}


// FormFactor_Wedge tests
void compare_true_FFt_cut_1D()
{
    double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7;
    double t_min = 0.0001, t_max = 0.1;
    double tx_min = 0, tx_max = 0.1, ty_min = 0, ty_max = 0.1;
    double phi_min = 0, phi_max = 2*pi, bins = 1000;
    double r_min = 0, r_max = 15;
    FormFactor_Wedge ff(A,Vo,R,a0,t_min,t_max,phi_min,phi_max, bins);
    FormFactor_t_1D ff2(A,Vo,R,a0,t_min,t_max,bins,r_min,r_max);

    TF1 *trueFF = ff2.getFormFactort_1D();
    double trueFFIntegral = trueFF->Integral(trueFF->GetXmin(), trueFF->GetXmax());
    TF1 *normalizedTrueFF = new TF1("", [trueFF, trueFFIntegral](double *x, double *par) 
    {
        return trueFF->Eval(x[0]) / trueFFIntegral;
    }, trueFF->GetXmin(), trueFF->GetXmax(), 0);
    normalizedTrueFF->SetLineColor(kBlack);
    normalizedTrueFF->SetLineStyle(3);
    normalizedTrueFF->Draw();

    TF1 *wedgeFF = ff.getWedgeCut1D_fun();
    double wedgeFFIntegral = wedgeFF->Integral(wedgeFF->GetXmin(), wedgeFF->GetXmax());
    TF1 *normalizedWedgeFF = new TF1("", [wedgeFF, wedgeFFIntegral](double *x, double *par) 
    {
        return wedgeFF->Eval(x[0]) / wedgeFFIntegral;
    }, wedgeFF->GetXmin(), wedgeFF->GetXmax(), 0);
    normalizedWedgeFF->SetLineColor(kRed);
    normalizedWedgeFF->SetLineStyle(2);
    normalizedWedgeFF->Draw("same");

    TLegend *legend = new TLegend(0.75, 0.8, 0.9, 0.9);
    legend->AddEntry(normalizedTrueFF, "True Form Factor", "l");
    legend->AddEntry(normalizedWedgeFF, "Wedge Form Factor", "l");
    legend->Draw();
}

void compare_true_FFt_cut_log_1D()
{
    double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7;
    double t_min = 0.0001, t_max = 0.1;
    double tx_min = 0, tx_max = 0.1, ty_min = 0, ty_max = 0.1;
    double phi_min = 0, phi_max = 2*pi, bins = 1000;
    double r_min = 0, r_max = 15;
    FormFactor_Wedge ff(A,Vo,R,a0,t_min,t_max,phi_min,phi_max,bins);
    FormFactor_t_1D ff2(A,Vo,R,a0,t_min,t_max,bins,r_min,r_max);
    
    TF1 *trueFF = ff2.getFormFactort_1D();
    double trueFFIntegral = trueFF->Integral(trueFF->GetXmin(), trueFF->GetXmax());
    TF1 *normalizedTrueFF = new TF1("", [trueFF, trueFFIntegral](double *x, double *par) 
    {
        return trueFF->Eval(x[0]) / trueFFIntegral;
    }, trueFF->GetXmin(), trueFF->GetXmax(), 0);
    normalizedTrueFF->SetLineColor(kBlack);
    normalizedTrueFF->SetLineStyle(3);
    normalizedTrueFF->Draw();

    TF1 *wedgeFF = ff.getWedgeCut1D_fun();
    double wedgeFFIntegral = wedgeFF->Integral(wedgeFF->GetXmin(), wedgeFF->GetXmax());
    TF1 *normalizedWedgeFF = new TF1("", [wedgeFF, wedgeFFIntegral](double *x, double *par) 
    {
        return wedgeFF->Eval(x[0]) / wedgeFFIntegral;
    }, wedgeFF->GetXmin(), wedgeFF->GetXmax(), 0);
    normalizedWedgeFF->SetLineColor(kRed);
    normalizedWedgeFF->SetLineStyle(2);
    normalizedWedgeFF->Draw("same");

    TLegend *legend = new TLegend(0.75, 0.8, 0.9, 0.9);
    legend->AddEntry(normalizedTrueFF, "True Form Factor", "l");
    legend->AddEntry(normalizedWedgeFF, "Wedge Form Factor", "l");
    legend->Draw();

    gPad->SetLogy(1);
}


// FormFactor_resolution tests
void smeared_FFt()
{
    double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7;
    double ty_min = 0, ty_max = 0.1, tx_prime_min = 0, tx_prime_max = 0.1, bins = 1000; 
    double sigma = 0.01;
    FormFactor_resolution ff(A,Vo,R,a0,ty_min,ty_max,tx_prime_min,tx_prime_max,bins,sigma);
    TF2 *test = ff.getSmearedFormFactor();
    test->Draw("colz");
    //ff.testwedge_getFF2t();
}

void smeared_FFt_hist()
{
    double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7;
    double ty_min = 0, ty_max = 0.1, tx_prime_min = 0, tx_prime_max = 0.1, bins = 1000; 
    double sigma = 0.01;
    FormFactor_resolution ff(A,Vo,R,a0,ty_min,ty_max,tx_prime_min,tx_prime_max,bins,sigma);
    TH2D *test = ff.getSmeared_hist();
    test->Draw();
}

void smeared_FFt_projection()
{
    double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7;
    double ty_min = 0, ty_max = 0.1, tx_prime_min = 0, tx_prime_max = 0.1, bins = 100; 
    double sigma = 0.01;
    FormFactor_resolution ff(A,Vo,R,a0,ty_min,ty_max,tx_prime_min,tx_prime_max,bins,sigma);
    TH2D *test = ff.getSmeared_hist();
    test->ProjectionY()->Draw();
}

void test_smeared_TF2()
{
    double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7;
    double ty_min = 0, ty_max = 0.1, tx_prime_min = 0, tx_prime_max = 0.1, bins = 100; 
    double sigma = 0.01;
    double x_min = 0, x_max = 10, y_min = 0, y_max = 10, tx_min = 0, tx_max = 0.1;
    FormFactor_resolution ff(A,Vo,R,a0,ty_min,ty_max,tx_prime_min,tx_prime_max,bins,sigma);
    TF2 *test = ff.getSmearedFormFactor();
    test->SetLineStyle(2);
    test->SetLineColor(kBlack);
    test->Draw();

    FormFactor_t_2D ff_2D(A,Vo,R,a0,tx_min,tx_max,ty_min,ty_max,bins,x_min,x_max,y_min,y_max);
    TF2 *check = ff_2D.getFormFactort2_2D();
    check->SetLineStyle(10);
    check->SetLineColor(kRed);
    check->Draw("same");
}


// FormFactor_saturation_data tests
void FF_phi_saturation_graph()
{
    // Data from phi saturation model plot scan
    const int bins = 72;
    double min = 0, max = .18;
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
    
    vector<double> x_vals_vec(x_vals, x_vals + bins); 
    vector<double> y_vals_vec(y_vals, y_vals + bins);

    FormFactor_saturation_data saturation(bins,x_vals,y_vals);
    TGraph *test = saturation.getGraph();
    test->Draw();
    gPad->SetLogy(1);
}

void FF_phi_saturation_TH1D()
{
    // Data from phi saturation model plot scan
    const int bins = 72;
    double min = 0, max = .18;
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
    
    vector<double> x_vals_vec(x_vals, x_vals + bins); 
    vector<double> y_vals_vec(y_vals, y_vals + bins);

    FormFactor_saturation_data saturation(bins,x_vals,y_vals);
    TH1D *test = saturation.getHist();
    test->Draw();
    gPad->SetLogy(1);
}

void FF_phi_saturation_compare()
{
    // Data from phi saturation model plot scan
    const int bins = 72;
    double min = 0, max = .18;
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
    
    vector<double> x_vals_vec(x_vals, x_vals + bins); 
    vector<double> y_vals_vec(y_vals, y_vals + bins);

    FormFactor_saturation_data saturation(bins,x_vals,y_vals);
    TH1D *test = saturation.getHist();
    test->SetMarkerColor(kBlack);
    test->SetMarkerStyle(1);
    test->Draw();
    gPad->SetLogy(1);

    double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7;
    double t_min = 0.0001, t_max = 0.1;
    double r_min = 0, r_max = 15;
    FormFactor_t_1D ff(A,Vo,R,a0,t_min,t_max,bins,r_min,r_max);
    TH1D *test2 = ff.getFF_hist();
    test2->Scale(197./test2->Integral(), "width");
    test2->SetLineColor(kRed);
    test2->Draw("same");
}

void FF_Jpsi_saturation_graph()
{
    // data from J/psi saturation model plot scan
    const int bins = 70;
    double min = 0, max = .18;
    double x_vals[bins] = {0.00136281334124255, 0.00417373659242085, 0.00679726496018727, 0.0090460035611299, 0.0114821370454844, 0.0141056654132508, 0.0165417988976054, 0.0189779323819599,
0.0216014607497263, 0.0240375942340808, 0.0266611226018472, 0.0339695230549108, 0.0364056565392653, 0.0390291849070318, 0.0414653183913863, 0.0439014518757408,
0.0465249802435072, 0.0487737188444499, 0.0513972472122163, 0.0538333806965708, 0.0564569090643372, 0.0588930425486917, 0.0613291760330463, 0.0637653095174008,
0.0662014430017553, 0.0686375764861098, 0.0712611048538762, 0.0736972383382308, 0.0763207667059972, 0.0787569001903517, 0.0811930336747062, 0.0838165620424727,
0.0862526955268272, 0.0886888290111817, 0.0913123573789481, 0.0937484908633026, 0.0961846243476572, 0.0986207578320117, 0.101056891316366, 0.103680419684133,
0.106116553168487, 0.108552686652842, 0.110988820137196, 0.113612348504963, 0.116048481989317, 0.118484615473672, 0.120920748958026, 0.123356882442381,
0.125980410810147, 0.12822914941109, 0.131040072662268, 0.133288811263211, 0.135912339630977, 0.138348473115332, 0.140784606599686, 0.143408134967453,
0.145844268451807, 0.148467796819574, 0.150903930303928, 0.153340063788283, 0.155776197272637, 0.158399725640404, 0.160835859124758, 0.163271992609113,
0.165895520976879, 0.168331654461233, 0.170767787945588, 0.173391316313354, 0.175640054914297, 0.178076188398652};
    double y_vals[bins] = {13448.9218102926, 6937.46196167787, 3204.7878500576, 1305.08427031405, 446.872520224438, 103.18239774925, 8.55278296449268,
11.005917498039, 31.6399824551623, 54.9298751813913, 66.3663164066068, 60.3779220946169, 46.9201466951133, 32.653225248372, 21.3360051124315,
12.4849168326106, 6.24035001776384, 2.58162492032996, 0.912279863030455, 0.596094495215495, 1.1375147582384, 2.17069581349244, 3.37486862233193,
4.3428584331333, 5.16499432428701, 5.58848994755982, 5.7674564552788, 5.67726804222898, 5.41507684299759, 4.84942389709097, 4.14229378615772,
3.48294584122589, 2.83767526873991, 2.17069581349244, 1.60896053571582, 1.15558518403178, 0.779253636732963, 0.485658601659348, 0.288701178756081,
0.163693431856085, 0.110384508079367, 0.100424245161287, 0.148922974953952, 0.197774527183426, 0.27106213627279, 0.343354491666976, 0.408354089359272,
0.470588418718807, 0.525479418352291, 0.533827122649325, 0.586773086901781, 0.577597441807907, 0.568565280566255, 0.533827122649325, 0.525479418352291,
0.493373717141273, 0.434927240554336, 0.377409027059443, 0.322376245262033, 0.27974266435082, 0.238951332128076, 0.188640778332054, 0.151288747835883,
0.123260124544862, 0.0885276988118201, 0.06258796389094, 0.0486375871383818, 0.032285236309424, 0.0250891049438687, 0.0191920533212076};

    vector<double> x_vals_vec(x_vals, x_vals + bins); 
    vector<double> y_vals_vec(y_vals, y_vals + bins);

    FormFactor_saturation_data saturation(bins,x_vals,y_vals);
    TGraph *test = saturation.getGraph();
    test->Draw();

    gPad->SetLogy(1);
}

void FF_Jpsi_saturation_hist()
{
    // data from J/psi saturation model plot scan
    const int bins = 70;
    double min = 0, max = .18;
    double x_vals[bins] = {0.00136281334124255, 0.00417373659242085, 0.00679726496018727, 0.0090460035611299, 0.0114821370454844, 0.0141056654132508, 0.0165417988976054, 0.0189779323819599,
0.0216014607497263, 0.0240375942340808, 0.0266611226018472, 0.0339695230549108, 0.0364056565392653, 0.0390291849070318, 0.0414653183913863, 0.0439014518757408,
0.0465249802435072, 0.0487737188444499, 0.0513972472122163, 0.0538333806965708, 0.0564569090643372, 0.0588930425486917, 0.0613291760330463, 0.0637653095174008,
0.0662014430017553, 0.0686375764861098, 0.0712611048538762, 0.0736972383382308, 0.0763207667059972, 0.0787569001903517, 0.0811930336747062, 0.0838165620424727,
0.0862526955268272, 0.0886888290111817, 0.0913123573789481, 0.0937484908633026, 0.0961846243476572, 0.0986207578320117, 0.101056891316366, 0.103680419684133,
0.106116553168487, 0.108552686652842, 0.110988820137196, 0.113612348504963, 0.116048481989317, 0.118484615473672, 0.120920748958026, 0.123356882442381,
0.125980410810147, 0.12822914941109, 0.131040072662268, 0.133288811263211, 0.135912339630977, 0.138348473115332, 0.140784606599686, 0.143408134967453,
0.145844268451807, 0.148467796819574, 0.150903930303928, 0.153340063788283, 0.155776197272637, 0.158399725640404, 0.160835859124758, 0.163271992609113,
0.165895520976879, 0.168331654461233, 0.170767787945588, 0.173391316313354, 0.175640054914297, 0.178076188398652};
    double y_vals[bins] = {13448.9218102926, 6937.46196167787, 3204.7878500576, 1305.08427031405, 446.872520224438, 103.18239774925, 8.55278296449268,
11.005917498039, 31.6399824551623, 54.9298751813913, 66.3663164066068, 60.3779220946169, 46.9201466951133, 32.653225248372, 21.3360051124315,
12.4849168326106, 6.24035001776384, 2.58162492032996, 0.912279863030455, 0.596094495215495, 1.1375147582384, 2.17069581349244, 3.37486862233193,
4.3428584331333, 5.16499432428701, 5.58848994755982, 5.7674564552788, 5.67726804222898, 5.41507684299759, 4.84942389709097, 4.14229378615772,
3.48294584122589, 2.83767526873991, 2.17069581349244, 1.60896053571582, 1.15558518403178, 0.779253636732963, 0.485658601659348, 0.288701178756081,
0.163693431856085, 0.110384508079367, 0.100424245161287, 0.148922974953952, 0.197774527183426, 0.27106213627279, 0.343354491666976, 0.408354089359272,
0.470588418718807, 0.525479418352291, 0.533827122649325, 0.586773086901781, 0.577597441807907, 0.568565280566255, 0.533827122649325, 0.525479418352291,
0.493373717141273, 0.434927240554336, 0.377409027059443, 0.322376245262033, 0.27974266435082, 0.238951332128076, 0.188640778332054, 0.151288747835883,
0.123260124544862, 0.0885276988118201, 0.06258796389094, 0.0486375871383818, 0.032285236309424, 0.0250891049438687, 0.0191920533212076};

    vector<double> x_vals_vec(x_vals, x_vals + bins); 
    vector<double> y_vals_vec(y_vals, y_vals + bins);

    FormFactor_saturation_data saturation(bins,x_vals,y_vals);
    TH1D *test = saturation.getHist();
    test->Draw();

    gPad->SetLogy(1);
}

void FF_Jpsi_saturation_compare()
{
    // data from J/psi saturation model plot scan
    const int bins = 70;
    double min = 0, max = .18;
    double x_vals[bins] = {0.00136281334124255, 0.00417373659242085, 0.00679726496018727, 0.0090460035611299, 0.0114821370454844, 0.0141056654132508, 0.0165417988976054, 0.0189779323819599,
0.0216014607497263, 0.0240375942340808, 0.0266611226018472, 0.0339695230549108, 0.0364056565392653, 0.0390291849070318, 0.0414653183913863, 0.0439014518757408,
0.0465249802435072, 0.0487737188444499, 0.0513972472122163, 0.0538333806965708, 0.0564569090643372, 0.0588930425486917, 0.0613291760330463, 0.0637653095174008,
0.0662014430017553, 0.0686375764861098, 0.0712611048538762, 0.0736972383382308, 0.0763207667059972, 0.0787569001903517, 0.0811930336747062, 0.0838165620424727,
0.0862526955268272, 0.0886888290111817, 0.0913123573789481, 0.0937484908633026, 0.0961846243476572, 0.0986207578320117, 0.101056891316366, 0.103680419684133,
0.106116553168487, 0.108552686652842, 0.110988820137196, 0.113612348504963, 0.116048481989317, 0.118484615473672, 0.120920748958026, 0.123356882442381,
0.125980410810147, 0.12822914941109, 0.131040072662268, 0.133288811263211, 0.135912339630977, 0.138348473115332, 0.140784606599686, 0.143408134967453,
0.145844268451807, 0.148467796819574, 0.150903930303928, 0.153340063788283, 0.155776197272637, 0.158399725640404, 0.160835859124758, 0.163271992609113,
0.165895520976879, 0.168331654461233, 0.170767787945588, 0.173391316313354, 0.175640054914297, 0.178076188398652};
    double y_vals[bins] = {13448.9218102926, 6937.46196167787, 3204.7878500576, 1305.08427031405, 446.872520224438, 103.18239774925, 8.55278296449268,
11.005917498039, 31.6399824551623, 54.9298751813913, 66.3663164066068, 60.3779220946169, 46.9201466951133, 32.653225248372, 21.3360051124315,
12.4849168326106, 6.24035001776384, 2.58162492032996, 0.912279863030455, 0.596094495215495, 1.1375147582384, 2.17069581349244, 3.37486862233193,
4.3428584331333, 5.16499432428701, 5.58848994755982, 5.7674564552788, 5.67726804222898, 5.41507684299759, 4.84942389709097, 4.14229378615772,
3.48294584122589, 2.83767526873991, 2.17069581349244, 1.60896053571582, 1.15558518403178, 0.779253636732963, 0.485658601659348, 0.288701178756081,
0.163693431856085, 0.110384508079367, 0.100424245161287, 0.148922974953952, 0.197774527183426, 0.27106213627279, 0.343354491666976, 0.408354089359272,
0.470588418718807, 0.525479418352291, 0.533827122649325, 0.586773086901781, 0.577597441807907, 0.568565280566255, 0.533827122649325, 0.525479418352291,
0.493373717141273, 0.434927240554336, 0.377409027059443, 0.322376245262033, 0.27974266435082, 0.238951332128076, 0.188640778332054, 0.151288747835883,
0.123260124544862, 0.0885276988118201, 0.06258796389094, 0.0486375871383818, 0.032285236309424, 0.0250891049438687, 0.0191920533212076};


    FormFactor_saturation_data saturation(bins,x_vals,y_vals);
    TH1D *test = saturation.getHist();
    test->SetMarkerColor(kBlack);
    test->SetMarkerStyle(1);
    test->Draw();
    gPad->SetLogy(1);

    double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7;
    double t_min = 0.0001, t_max = 0.1;
    double r_min = 0, r_max = 15;
    FormFactor_t_1D ff(A,Vo,R,a0,t_min,t_max,bins,r_min,r_max);
    TH1D *test2 = ff.getFF_hist();
    test2->Scale(197./test2->Integral(), "width");
    test2->SetLineColor(kRed);
    test2->Draw("same");
}


// FF_WS_Transforms 
void compare_transforms_F_to_G()
{
double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7, a = 0.535;
double ty_min = 0, ty_max = 0.25, tx_min = 0, tx_max = 0.25, bins = 1000;
double x_min = 0, x_max = 15, y_min = 0, y_max = 15;

// Initialize constructors
FF_WS_Transforms_2D transformer_2D;
FormFactor_t_2D ff_2d(A,Vo,R,a0,tx_min,tx_max,ty_min,ty_max,bins,x_min,x_max,y_min,y_max);
WoodsSaxon_2D ws_2d(Vo,R,a,x_min,x_max,y_min,y_max,bins,tx_min,tx_max,ty_min,ty_max);

TF2 *ws_fun = ws_2d.getWoodsSaxon2D_fun();
    double trueWSIntegral = ws_fun->Integral(ws_fun->GetXmin(), ws_fun->GetXmax());
    TF2 *normalizedTrueWS = new TF2("", [ws_fun, trueWSIntegral](double *x, double *par) 
    {
        return ws_fun->Eval(x[0], x[1]) / trueWSIntegral;
    }, ws_fun->GetXmin(), ws_fun->GetXmax(), ws_fun->GetYmin(), ws_fun->GetYmax(), 0);



TF2 *transformed_ff = ff_2d.getTransformed_TF2();
    double transformedFFIntegral = transformed_ff->Integral(transformed_ff->GetXmin(), transformed_ff->GetXmax());
    TF2 *normalizedTransformedFF = new TF2("", [transformed_ff, transformedFFIntegral](double *x, double *par) 
    {
        return transformed_ff->Eval(x[0], x[1]) / transformedFFIntegral;
    }, transformed_ff->GetXmin(), transformed_ff->GetXmax(), transformed_ff->GetYmin(), transformed_ff->GetYmax(), 0);



//TF2 *transformer_fun_FtoG = transformer_2D.transformTF2_FtoG(transformed_ff, x_min, x_max, y_min, y_max);
//    double transFFIntegral = transformer_fun_FtoG->Integral(transformer_fun_FtoG->GetXmin(), transformer_fun_FtoG->GetXmax());
//    TF2 *normalizedTransFF = new TF2("", [transformer_fun_FtoG, transFFIntegral](double *x, double *par) 
//    {
//        return transformer_fun_FtoG->Eval(x[0], x[1]) / transFFIntegral;
//    }, transformer_fun_FtoG->GetXmin(), transformer_fun_FtoG->GetXmax(), transformer_fun_FtoG->GetYmin(), transformer_fun_FtoG->GetYmax(), 0);

    normalizedTrueWS->SetTitle("Fourier-Bessel Transformation from |F(t_{x},t_{y})|^{2} to G(x,y)");
normalizedTrueWS->GetZaxis()->SetTitle("G(x,y) [fm^{-3}");
normalizedTrueWS->GetXaxis()->SetTitle("x [fm]");
normalizedTrueWS->GetYaxis()->SetTitle("y [fm]");
normalizedTrueWS->SetLineStyle(2);
normalizedTrueWS->SetLineColor(kRed);
normalizedTrueWS->Draw();

normalizedTransformedFF->SetLineStyle(10);
normalizedTransformedFF->SetLineColor(kOrange);
normalizedTransformedFF->Draw("same");

//normalizedTransFF->SetLineStyle(3);
//normalizedTransFF->SetLineColor(kBlue);
//normalizedTransFF->Draw("same");

auto legend = new TLegend(0.63,0.75,0.9,0.9);
	legend->AddEntry(normalizedTrueWS,"True Woods-Saxon","l");
    legend->AddEntry(normalizedTransformedFF,"Transformed Form Factor","l");
//    legend->AddEntry(normalizedTransFF,"Transformed FF with Transformer Func","l");
    legend->Draw();


}

void compare_transforms_G_to_F()
{

double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7, a = 0.535;
double ty_min = 0, ty_max = 0.25, tx_min = 0, tx_max = 0.25, bins = 1000;
double x_min = 0, x_max = 15, y_min = 0, y_max = 15;

// Initialize constructors
FF_WS_Transforms_2D transformer_2D;
FormFactor_t_2D ff_2d(A,Vo,R,a0,tx_min,tx_max,ty_min,ty_max,bins,x_min,x_max,y_min,y_max);
WoodsSaxon_2D ws_2d(Vo,R,a,x_min,x_max,y_min,y_max,bins,tx_min,tx_max,ty_min,ty_max);

    TF2 *ws_fun = ws_2d.getWoodsSaxon2D_fun();
    //TF2 *transformer_fun_GtoF = transformer_2D.transformTF2_GtoF(ws_fun, tx_min, tx_max, ty_min, ty_max);
    //double transGIntegral = transformer_fun_GtoF->Integral(transformer_fun_GtoF->GetXmin(), transformer_fun_GtoF->GetXmax());
    //TF2 *normalizedTransformedG = new TF2("", [transformer_fun_GtoF, transGIntegral](double *x, double *par) 
    //{
    //    return transformer_fun_GtoF->Eval(x[0],x[1]) / transGIntegral;
    //}, transformer_fun_GtoF->GetXmin(), transformer_fun_GtoF->GetXmax(), transformer_fun_GtoF->GetYmin(), transformer_fun_GtoF->GetYmax(), 0);

    //TF2 *transformed_ws = ws_2d.get_transformed_WoodsSaxon2D_fun();
    //double transformedWSIntegral = transformed_ws->Integral(transformed_ws->GetXmin(), transformed_ws->GetXmax());
    //TF2 *normalizedTransformedWS = new TF2("", [transformed_ws, transformedWSIntegral](double *x, double *par) 
    //{
    //    return transformed_ws->Eval(x[0], x[1]) / transformedWSIntegral;
    //}, transformed_ws->GetXmin(), transformed_ws->GetXmax(), transformed_ws->GetYmin(), transformed_ws->GetYmax(), 0);

    TF2 *ff_fun = ff_2d.getFormFactort2_2D();
    double trueFFIntegral = ff_fun->Integral(ff_fun->GetXmin(), ff_fun->GetXmax());
    TF2 *normalizedTrueFF = new TF2("", [ff_fun, trueFFIntegral](double *x, double *par) 
    {
        return ff_fun->Eval(x[0], x[1]) / trueFFIntegral;
    }, ff_fun->GetXmin(), ff_fun->GetXmax(), ff_fun->GetYmin(), ff_fun->GetYmax(), 0);

    //normalizedTransformedWS->SetTitle("Fourier-Bessel Transformation from G(x,y) to |F(t_{x},t_{y})|^{2}");
//normalizedTransformedWS->GetYaxis()->SetTitle("|F(t_{x},t_{y})|^{2}");
//normalizedTransformedWS->GetXaxis()->SetTitle("t_{x} [GeV^{2}]");
//normalizedTransformedWS->GetYaxis()->SetTitle("t_{y} [GeV^{2}]");
//normalizedTransformedWS->SetLineColor(kBlack);
//normalizedTransformedWS->SetLineStyle(3);
//normalizedTransformedWS->Draw();

normalizedTrueFF->SetLineColor(kBlue);
normalizedTrueFF->SetLineStyle(10);
normalizedTrueFF->Draw("same");

//normalizedTransformedG->SetLineColor(kRed);
//normalizedTransformedG->SetLineStyle(2);
//normalizedTransformedG->Draw("same");

auto legend = new TLegend(0.63,0.75,0.9,0.9);
	//legend->AddEntry(normalizedTransformedWS,"Transformed Woods-Saxon","l");
    legend->AddEntry(normalizedTrueFF,"True Form Factor","l");
//    legend->AddEntry(normalizedTransformedG,"Transformed WS with Transformer Func","l");
    legend->Draw();


//gPad->SetLogy(1);
}

void compare_transform_F_to_G_hist()
{
    double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7, a = 0.535;
double ty_min = 0, ty_max = 0.25, tx_min = 0, tx_max = 0.25, bins = 1000;
double x_min = 0, x_max = 15, y_min = 0, y_max = 15;

// Initialize constructors
FF_WS_Transforms_2D transformer_2D;
FormFactor_t_2D ff_2d(A,Vo,R,a0,tx_min,tx_max,ty_min,ty_max,bins,x_min,x_max,y_min,y_max);
WoodsSaxon_2D ws_2d(Vo,R,a,x_min,x_max,y_min,y_max,bins,tx_min,tx_max,ty_min,ty_max);

    TH2D *ff_transform_hist = ff_2d.getTransformed_hist();
TH2D *ws_hist = ws_2d.getWoodsSaxon2D_hist();

ff_transform_hist->Scale(197./ff_transform_hist->Integral(), "width");
ff_transform_hist->SetTitle("Fourier-Bessel Transformation Histogram from |F(t_{x},t_{y})|^{2} to G(x,y)");
ff_transform_hist->GetZaxis()->SetTitle("G(x,y) [fm^{-3}]");
ff_transform_hist->GetXaxis()->SetTitle("x [fm]");
ff_transform_hist->GetYaxis()->SetTitle("y [fm]");
ff_transform_hist->SetLineColor(kBlack);
ff_transform_hist->SetLineStyle(3);
ff_transform_hist->Draw();

ws_hist->Scale(197./ws_hist->Integral(), "width");
ws_hist->SetLineColor(kRed);
ws_hist->SetLineStyle(10);
ws_hist->Draw("same");

auto legend = new TLegend(0.63,0.75,0.9,0.9);
	legend->AddEntry(ws_hist,"True Woods-Saxon","l");
    legend->AddEntry(ff_transform_hist,"Transformed Form Factor","l");
    legend->Draw();

gStyle->SetOptStat(0);
gPad->SetLogy(0);

}

void compare_G_to_F_hist()
{
    double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7, a = 0.535;
double ty_min = 0, ty_max = 0.25, tx_min = 0, tx_max = 0.25, bins = 1000;
double x_min = 0, x_max = 15, y_min = 0, y_max = 15;

// Initialize constructors
FF_WS_Transforms_2D transformer_2D;
FormFactor_t_2D ff_2d(A,Vo,R,a0,tx_min,tx_max,ty_min,ty_max,bins,x_min,x_max,y_min,y_max);
WoodsSaxon_2D ws_2d(Vo,R,a,x_min,x_max,y_min,y_max,bins,tx_min,tx_max,ty_min,ty_max);

    TH2D *ws_transform_hist = ws_2d.get_transformed_WoodsSaxon2D_hist();
TH2D *ff_hist = ff_2d.getFormFactort_hist();

ff_hist->Scale(197./ff_hist->Integral(), "width");
ff_hist->SetTitle("Fourier-Bessel Transformation Histogram from G(x,y) to |F(t_{x},t_{y)|^{2}");
ff_hist->GetZaxis()->SetTitle("|F(t_{x},t_{y})|^{2}");
ff_hist->GetXaxis()->SetTitle("t_{x} [GeV^{2}]");
ff_hist->GetYaxis()->SetTitle("t_{y} [GeV^{2}]");
ff_hist->SetLineColor(kBlack);
ff_hist->SetLineStyle(3);
ff_hist->ProjectionY()->Draw();

ws_transform_hist->Scale(197./ws_transform_hist->Integral(), "width");
ws_transform_hist->SetLineColor(kRed);
ws_transform_hist->SetLineStyle(2);
ws_transform_hist->ProjectionY()->Draw("same");

auto legend = new TLegend(0.63,0.75,0.9,0.9);
	legend->AddEntry(ff_hist,"True Form Factor","l");
    legend->AddEntry(ws_transform_hist,"Transformed Woods-Saxon","l");
    legend->Draw();

gStyle->SetOptStat(0);
gPad->SetLogy(1);

}


// FormFactor_resolution_add_wedge tests
void compare_1D_true_wResCut()
{
    double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7;
double q_min = 0, q_max = 0.5, t_min = 0, t_max = 0.18;
double qy_min = 0, qy_max = 0.5, qx_prime_min = 0, qx_prime_max = 0.5;
double ty_min = 0, ty_max = 0.18, tx_min = 0, tx_max = 0.18, tx_prime_min = 0, tx_prime_max = 0.18, bins = 1000; 
double x_min = 0, x_max = 15, y_min = 0, y_max = 15, r_min = 0, r_max = 15;
double phi_min = 0, phi_max = pi/12, phi_max_whole = pi/2, sigma = 0.025;

    FormFactor_resolution_add_wedge_1D ff(A,Vo,R,a0,t_min,t_max,bins,phi_min,phi_max,sigma,r_min,r_max);
    FormFactor_t_1D trueff(A,Vo,R,a0,t_min,t_max,bins,r_min,r_max);
    
    TF1 *test = trueff.getFormFactort_1D();
    TF1 *check = ff.getWedgeRes_fun_1D();
    
    check->SetLineColor(kRed);
    check->SetLineStyle(2);
    check->Draw();

    test->SetLineColor(kBlack);
    test->SetLineStyle(10);
    test->Draw("same");
}

void compare_1D_true_wResCut_hist()
{
    double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7;
double q_min = 0, q_max = 0.5, t_min = 0, t_max = 0.25;
double qy_min = 0, qy_max = 0.5, qx_prime_min = 0, qx_prime_max = 0.5;
double ty_min = 0, ty_max = 0.25, tx_min = 0, tx_max = 0.25, tx_prime_min = 0, tx_prime_max = 0.25, bins = 10; 
double x_min = 0, x_max = 15, y_min = 0, y_max = 15, r_min = 0, r_max = 15;
double phi_min = 0, phi_max = pi/9, sigma = 0.1;

    FormFactor_resolution_add_wedge_1D ff(A,Vo,R,a0,t_min,t_max,bins,phi_min,phi_max,sigma,r_min,r_max);
    //FormFactor_t_1D trueff(A,Vo,R,a0,t_min,t_max,bins,r_min,r_max);
    
    //TH1D *test = trueff.getFF_hist();
    TH1D *check = ff.getWedgeRes_hist_1D();
    
    check->SetLineColor(kRed);
    check->SetLineStyle(2);
    check->Draw();

    //test->SetLineColor(kBlack);
    //test->SetLineStyle(10);
    //test->Draw("same");
}

void compare_2D_true_wResCut()
{
    double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7;
double ty_min = 0, ty_max = 0.25, tx_prime_min = 0, tx_prime_max = 0.25, bins = 10; 
double x_min = 0, x_max = 15, y_min = 0, y_max = 15;
double phi_min = 0, phi_max = pi/9, sigma = 0.1;
double tx_min = 0, tx_max = 0.25;

    FormFactor_resolution_add_wedge_2D test(A,Vo,R,a0,ty_min,ty_max,tx_prime_min,tx_prime_max,bins,phi_min,phi_max,sigma,x_min,x_max,y_min,y_max);
    FormFactor_t_2D ff(A,Vo,R,a0,tx_min,tx_max,ty_min,ty_max,bins,x_min,x_max,y_min,y_max);

    TF2 *check = test.getWedgeRes_fun_2D();
    TF2 *trueFF = ff.getFormFactort2_2D();

    check->SetLineColor(kRed);
    check->SetLineStyle(2);
    check->Draw();

    trueFF->SetLineColor(kBlack);
    trueFF->SetLineStyle(10);
    trueFF->Draw("same");
}

void compare_2D_true_wResCut_hist()
{
    double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7;
double ty_min = 0, ty_max = 0.25, tx_prime_min = 0, tx_prime_max = 0.25, bins = 10; 
double x_min = 0, x_max = 15, y_min = 0, y_max = 15;
double phi_min = 0, phi_max = pi/9, sigma = 0.1;
double tx_min = 0, tx_max = 0.25;

    FormFactor_resolution_add_wedge_2D test(A,Vo,R,a0,ty_min,ty_max,tx_prime_min,tx_prime_max,bins,phi_min,phi_max,sigma,x_min,x_max,y_min,y_max);
    FormFactor_t_2D ff(A,Vo,R,a0,tx_min,tx_max,ty_min,ty_max,bins,x_min,x_max,y_min,y_max);

    TH2D *check = test.getWedgeRes_hist_2D();
    TH2D *trueFF = ff.getFormFactort_hist();

    check->SetLineColor(kRed);
    check->SetLineStyle(2);
    check->ProjectionY()->Draw();

    trueFF->SetLineColor(kBlack);
    trueFF->SetLineStyle(10);
    trueFF->ProjectionY()->Draw("same");
}

void compare_more_hists()
{
    double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7;
double ty_min = 0, ty_max = 0.25, tx_prime_min = 0, tx_prime_max = 0.25, bins_2d = 10; 
double x_min = 0, x_max = 15, y_min = 0, y_max = 15;
double phi_min = 0, phi_max = pi/9, sigma = 0.1;
double tx_min = 0, tx_max = 0.25;
double sigma25 = 0.025, sigma50 = 0.05, sigma100 = 0.1, sigma150 = 0.15;
FormFactor_t_2D ff_2d(A,Vo,R,a0,tx_min,tx_max,ty_min,ty_max,bins_2d,x_min,x_max,y_min,y_max);
FormFactor_resolution_add_wedge_2D ff_wResCut_2d_25(A,Vo,R,a0,ty_min,ty_max,tx_prime_min,tx_prime_max,bins_2d,phi_min,phi_max,sigma25,x_min,x_max,y_min,y_max);
FormFactor_resolution_add_wedge_2D ff_wResCut_2d_50(A,Vo,R,a0,ty_min,ty_max,tx_prime_min,tx_prime_max,bins_2d,phi_min,phi_max,sigma50,x_min,x_max,y_min,y_max);
FormFactor_resolution_add_wedge_2D ff_wResCut_2d_100(A,Vo,R,a0,ty_min,ty_max,tx_prime_min,tx_prime_max,bins_2d,phi_min,phi_max,sigma100,x_min,x_max,y_min,y_max);
FormFactor_resolution_add_wedge_2D ff_wResCut_2d_150(A,Vo,R,a0,ty_min,ty_max,tx_prime_min,tx_prime_max,bins_2d,phi_min,phi_max,sigma150,x_min,x_max,y_min,y_max);

    TCanvas *c1 = new TCanvas("", "", 1200, 800);

TH2D *ff_cutwRes_hist2 = ff_wResCut_2d_100.getWedgeRes_hist_2D();
    ff_cutwRes_hist2->Scale(197./ff_cutwRes_hist2->Integral(), "width");
    ff_cutwRes_hist2->SetTitle("|F(t_{x},t_{y})|^{2} y-Projection");
    ff_cutwRes_hist2->GetZaxis()->SetTitle("|F(t_{x},t_{y})|^{2}");
    ff_cutwRes_hist2->GetYaxis()->SetTitle("t_{y} [GeV^{2}]");
    ff_cutwRes_hist2->GetXaxis()->SetTitle("t_{x} [GeV^{2}]");
    ff_cutwRes_hist2->SetLineColor(kRed);
    ff_cutwRes_hist2->SetLineStyle(2);
    ff_cutwRes_hist2->ProjectionY()->Draw();
    
TH2D *trueFF_hist = ff_2d.getFormFactort_hist();
    trueFF_hist->Scale(197./trueFF_hist->Integral(), "width");
    trueFF_hist->SetLineColor(kBlack);
    trueFF_hist->SetLineStyle(10);
    trueFF_hist->ProjectionY()->Draw("same");

auto legend = new TLegend(0.63,0.75,0.9,0.9);
	legend->AddEntry(ff_cutwRes_hist2,"Form Factor with 100 MeV resolution","l");
    legend->AddEntry(trueFF_hist,"True Form Factor","l");
    legend->Draw();

gStyle->SetOptStat(0);
gPad->SetLogy(1);
c1->Draw();

gPad->SetLogy(1);
gStyle->SetOptStat(0);
c1->Update();
c1->Draw();
}


// Plots for paper
void result_plot()
{
     double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7, a = 0.535;
double bins = 1000, phi_min = 0, phi_max_cut = pi/12, phi_max_whole = pi/2;

//sigma = 25 MeV
double sigma = 0.025;

// 1D params
double t_min = 0, t_max = 0.2, q_min = 0, q_max = 0.5;
double r_min = 0, r_max = 15;

// 2D params
double tx_min = 0, tx_max = 0.2, ty_min = 0, ty_max = 0.2;
double x_min = 0, x_max = 15, y_min = 0, y_max = 15;
double tx_prime_min = 0, tx_prime_max = 0.2;

FormFactor_t_1D ff_1d(A,Vo,R,a0,t_min,t_max,bins,r_min,r_max);

FormFactor_resolution_add_wedge_1D ff_wResCut_1d_cut(A,Vo,R,a0,t_min,t_max,bins,phi_min,phi_max_cut,sigma,r_min,r_max);

FormFactor_resolution_add_wedge_1D ff_wResCut_1d_whole(A,Vo,R,a0,t_min,t_max,bins,phi_min,phi_max_whole,sigma,r_min,r_max);

TCanvas *c = new TCanvas("c", "Canvas", 800, 600);
c->SetLeftMargin(0.12); 
c->SetBottomMargin(0.13);

TH1D *ff_true_hist = ff_1d.getFF_hist();
    ff_true_hist->GetYaxis()->SetTitleSize(0.047);
    ff_true_hist->GetXaxis()->SetTitleSize(0.047);
    ff_true_hist->SetTitle("");
    ff_true_hist->GetYaxis()->SetTitleOffset(1.03);
    ff_true_hist->GetXaxis()->SetTitleOffset(0.97);
    ff_true_hist->GetYaxis()->SetTitle("|F_{#hat{n}}(t)|^{2}");
    ff_true_hist->GetXaxis()->SetTitle("|t| [GeV^{2}/c^{2}]");
    ff_true_hist->GetYaxis()->SetRangeUser(1e-5,1e5);
    ff_true_hist->SetLineStyle(1);
    ff_true_hist->SetLineWidth(2);
    ff_true_hist->SetLineColor(kBlack);
    ff_true_hist->Scale(197./ff_true_hist->Integral(), "width");
    ff_true_hist->Draw();
      
TH1D *ff_hist_25 = ff_wResCut_1d_cut.getWedgeRes_hist_1D();
    ff_hist_25->GetYaxis()->SetRangeUser(1e-5,1e5);
    ff_hist_25->Scale(197./ff_hist_25->Integral(), "width");
    //ff_hist_25->SetMarkerStyle(3);
    //ff_hist_25->SetMarkerSize(0.9);
    //ff_hist_25->SetMarkerColor(kMagenta);
    ff_hist_25->SetLineStyle(2);
    ff_hist_25->SetLineColor(kMagenta);
    ff_hist_25->SetLineWidth(2);
    ff_hist_25->Draw("same");

TH1D *ff_hist_50 = ff_wResCut_1d_whole.getWedgeRes_hist_1D();
    ff_hist_50->GetYaxis()->SetRangeUser(1e-5,1e5);
    ff_hist_50->Scale(197./ff_hist_50->Integral(), "width");
    //ff_hist_50->SetMarkerStyle(27);
    //ff_hist_50->SetMarkerSize(0.9);
    //ff_hist_50->SetMarkerColor(kBlue);
    ff_hist_50->SetLineColor(kBlue);
    ff_hist_50->SetLineWidth(2);
    ff_hist_50->SetLineStyle(9);
    ff_hist_50->Draw("same");

auto legend = new TLegend(0.45,0.65,0.86,0.86);
    legend->SetHeader("Coherent J/#psi","C");
    legend->AddEntry(ff_true_hist,"Truth","l");
	legend->AddEntry(ff_hist_25,"25 MeV/c res. w. #theta_{max} = #pi/12","l");
    legend->AddEntry(ff_hist_50,"25 MeV/c res. w. #theta_{max} = #pi/2","l");
    legend->SetTextSize(0.04);
    legend->SetBorderSize(0);
    legend->Draw();

gPad->SetLogy(1);
gStyle->SetOptStat(0);
c->Draw();
}

// FormFactor_transform_resolution_add_wedge_1D tests
/*void ff_transform_wResCut_fun()
{
    double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7;
    double bins = 10, r_min = 0, r_max = 15;
    double phi_min = 0, phi_max = pi/6, sigma = 0.1;

    FormFactor_transform_resolution_add_wedge_1D ff(A,Vo,R,a0,bins,phi_min,phi_max,sigma,r_min,r_max);
    TF1 *test = ff.getWedgeResTransform_fun_1D();
    test->Draw();*/
/*
    FormFactor_transform_resolution_add_wedge_1D ff;
    TF1 *test = ff.transform_FF_to_G_fun(phi_min,phi_max,sigma,r_min,r_max);
    test->Draw();
    */
//}

/*void ff_transform_wResCut_hist()
{
    double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7, a = 0.535;
    double bins = 10, r_min = 0, r_max = 15;
    double t_min = 0, t_max = 0.1;
    double phi_min = 0, phi_max = pi/6, sigma = 0.1;

     FormFactor_transform_resolution_add_wedge_1D ff(A,Vo,R,a0,bins,phi_min,phi_max,sigma,r_min,r_max);
    TH1D *test = ff.getWedgeResTransform_hist_1D();
    test->Draw();
    */
    // Initialize constructor
/*FormFactor_transform_resolution_add_wedge_1D ff_trans_wResCut(A,Vo,R,a0,bins,phi_min,phi_max,sigma,r_min,r_max);
TF1 *ff_transform_wResCut = ff_wResCut.getWedgeResTransform_fun_1D();
ff_transform_wResCut->Draw();


   FormFactor_transform_resolution_add_wedge_1D ff;
    TH1D *test = ff.transform_FF_to_G_hist(phi_min,phi_max,sigma,r_min,r_max,bins);
    test->Draw();

    FormFactor_transform_resolution_add_wedge_1D ff1(A,Vo,R,a0,bins,phi_min,phi_max,sigma,r_min,r_max);
    TH1D *test1 = ff1.getWedgeResTransform_hist_1D();
    test1->SetTitle("Form Factor with Resolution and #pi/9 Wedge Cut Transformed");
    test1->GetXaxis()->SetTitle("r [fm]");
    test1->GetYaxis()->SetTitle("G(r) [fm^{-3}]");
    test1->SetLineStyle(2);
    test1->SetLineColor(kBlue);
    test1->Draw();

double sigma2 = .15;
    FormFactor_transform_resolution_add_wedge_1D ff2(A,Vo,R,a0,bins,phi_min,phi_max,sigma2,r_min,r_max);
    TH1D *test2 = ff2.getWedgeResTransform_hist_1D();
    test2->Draw("same");

double sigma3 = 0.05;
    FormFactor_transform_resolution_add_wedge_1D ff3(A,Vo,R,a0,bins,phi_min,phi_max,sigma3,r_min,r_max);
    TH1D *test3 = ff3.getWedgeResTransform_hist_1D();
    test3->Draw("same");

    FormFactor_t_1D ff4(A,Vo,R,a0,t_min,t_max,bins,r_min,r_max);
    TH1D *test4 = ff4.getFF_hist();
    test4->Draw();

    WoodsSaxon_1D ws(Vo,R,a,r_min,r_max,bins,t_min,t_max);
    TH1D *test5 = ws.getWS_hist();
    test5->Draw("same");

    auto legend = new TLegend(0.43,0.71,0.73,0.9);
    legend->AddEntry(test1,"|F(t)|^{2} -> G(r) with 50 MeV res","l");
    legend->AddEntry(test2,"|F(t)|^{2} -> G(r) with 50 MeV res","l");
    legend->AddEntry(test3,"|F(t)|^{2} -> G(r) with 50 MeV res","l");
    legend->AddEntry(test4,"|F(t)|^{2} -> G(r) with 50 MeV res","l");
    legend->AddEntry(test5,"|F(t)|^{2} -> G(r) with 50 MeV res","p");
    legend->Draw();
*/
    
//}



void a()
{
    double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7, a = 0.535;
double bins = 1000, r_min = 0, r_max = 15;
double phi_min = 0, phi_max = pi/9, sigma = 0.1;
double t_min = 0, t_max = 0.1;

// Initialize constructor
//FormFactor_transform_resolution_add_wedge_1D ff_trans_wResCut;
WoodsSaxon_1D ws_1d(Vo,R,a,r_min,r_max,bins,t_min,t_max);
FormFactor_t_1D ff_1d(A,Vo,R,a0,t_min,t_max,bins,r_min,r_max);
/*
TH1D *transformed_ff_wResCut = ff_trans_wResCut.transform_FF_to_G_hist(phi_min,phi_max,sigma,r_min,r_max,bins);
    transformed_ff_wResCut->SetTitle("Histogram of |F(t)|^{2} with Resolution and #pi/9 Wedge Cut");
    transformed_ff_wResCut->GetYaxis()->SetTitle("|F(t)|^{2}");
    transformed_ff_wResCut->GetXaxis()->SetTitle("t [GeV^{2}]");
    transformed_ff_wResCut->Scale(197./transformed_ff_wResCut->Integral(), "width");
    transformed_ff_wResCut->SetLineStyle(10);
    transformed_ff_wResCut->SetLineColor(kRed);
    transformed_ff_wResCut->Draw();
*/
TH1D *ff_hist_5 = ff_1d.getTransformedFF_hist();
    ff_hist_5->SetLineStyle(4);
    ff_hist_5->SetLineColor(kOrange);
    ff_hist_5->Scale(197./ff_hist_5->Integral(), "width");
    ff_hist_5->Draw("same");

TH1D *ff_hist_50 = ws_1d.getWoodsSaxon1D_hist();
    ff_hist_50->SetLineStyle(2);
    ff_hist_50->SetLineColor(kBlack);
    ff_hist_50->Scale(197./ff_hist_50->Integral(), "width");
    ff_hist_50->Draw("same");

auto legend = new TLegend(0.53,0.65,0.9,0.9);
	legend->AddEntry(ff_hist_5,"True Woods-Saxon","l");
    legend->AddEntry(ff_hist_50,"Transformed True Form Factor","l");
    //legend->AddEntry(transformed_ff_wResCut,"Transformed FF with Added Resolution and Wedge Cut","l");
    legend->Draw();

gPad->SetLogy(0);
gStyle->SetOptStat(0);
}

void b()
{

    double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7, a = 0.535;
double bins = 1000, phi_min = 0, phi_max = pi/9;

// 1D params
double t_min = 0, t_max = 0.1, q_min = 0, q_max = 0.316;
double r_min = 0, r_max = 15;

WoodsSaxon_1D ws_1d(Vo,R,a,r_min,r_max,bins,t_min,t_max);

FormFactor_t_1D ff_1d(A,Vo,R,a0,t_min,t_max,bins,r_min,r_max);

// sigma = 5% MeV
double sigma5 = 0.005;
    FormFactor_resolution_add_wedge_1D ff_wResCut_1d_5(A,Vo,R,a0,t_min,t_max,bins,phi_min,phi_max,sigma5,r_min,r_max);

// sigma = 50% MeV
double sigma50 = 0.05;
    FormFactor_resolution_add_wedge_1D ff_wResCut_1d_50(A,Vo,R,a0,t_min,t_max,bins,phi_min,phi_max,sigma50,r_min,r_max);

// sigma = 100MeV
double sigma100 = 0.1;
    FormFactor_resolution_add_wedge_1D ff_wResCut_1d_100(A,Vo,R,a0,t_min,t_max,bins,phi_min,phi_max,sigma100,r_min,r_max);

// sigma = 150 MeV
double sigma150 = 0.15;
    FormFactor_resolution_add_wedge_1D ff_wResCut_1d_150(A,Vo,R,a0,t_min,t_max,bins,phi_min,phi_max,sigma150,r_min,r_max);

   TF1 *trueff_1d = ff_1d.getFormFactort_1D();
    double trueFFIntegral = trueff_1d->Integral(trueff_1d->GetXmin(), trueff_1d->GetXmax());
    TF1 *trueFF = new TF1("", [trueff_1d, trueFFIntegral](double *x, double *par) 
    {
        return trueff_1d->Eval(x[0]) / trueFFIntegral;
    }, trueff_1d->GetXmin(), trueff_1d->GetXmax(), 0);
    trueFF->SetTitle("Plot |F(t)|^{2} with Resolution and #pi/9 Wedge Cut");
    trueFF->GetYaxis()->SetTitle("|F(t)|^{2}");
    trueFF->GetXaxis()->SetTitle("t [GeV^{2}]");
    trueFF->SetLineColor(kRed);
    trueFF->SetLineStyle(2);
    trueFF->Draw();

TF1 *ff_fun_5 = ff_wResCut_1d_5.getWedgeRes_fun_1D();
    double FFwRes5Integral = ff_fun_5->Integral(ff_fun_5->GetXmin(), ff_fun_5->GetXmax());
    TF1 *FF_wCutRes_5 = new TF1("", [ff_fun_5, FFwRes5Integral](double *x, double *par) 
    {
        return ff_fun_5->Eval(x[0]) / FFwRes5Integral;
    }, ff_fun_5->GetXmin(), ff_fun_5->GetXmax(), 0);
    FF_wCutRes_5->SetLineColor(kBlue);
    FF_wCutRes_5->SetLineStyle(10);
    FF_wCutRes_5->Draw("same");

TF1 *ff_fun_50 = ff_wResCut_1d_50.getWedgeRes_fun_1D();
    double FFwRes50Integral = ff_fun_50->Integral(ff_fun_50->GetXmin(), ff_fun_50->GetXmax());
    TF1 *FF_wCutRes_50 = new TF1("", [ff_fun_50, FFwRes50Integral](double *x, double *par) 
    {
        return ff_fun_50->Eval(x[0]) / FFwRes50Integral;
    }, ff_fun_50->GetXmin(), ff_fun_50->GetXmax(), 0);
    FF_wCutRes_50->SetLineColor(kMagenta);
    FF_wCutRes_50->SetLineStyle(5);
    FF_wCutRes_50->Draw("same");

TF1 *ff_fun_100 = ff_wResCut_1d_100.getWedgeRes_fun_1D();
    double FFwRes100Integral = ff_fun_100->Integral(ff_fun_100->GetXmin(), ff_fun_100->GetXmax());
    TF1 *FF_wCutRes_100 = new TF1("", [ff_fun_100, FFwRes100Integral](double *x, double *par) 
    {
        return ff_fun_100->Eval(x[0]) / FFwRes100Integral;
    }, ff_fun_100->GetXmin(), ff_fun_100->GetXmax(), 0);
    FF_wCutRes_100->SetLineColor(kBlack);
    FF_wCutRes_100->SetLineStyle(4);
    FF_wCutRes_100->Draw("same");

TF1 *ff_fun_150 = ff_wResCut_1d_150.getWedgeRes_fun_1D();
    double FFwRes150Integral = ff_fun_150->Integral(ff_fun_150->GetXmin(), ff_fun_150->GetXmax());
    TF1 *FF_wCutRes_150 = new TF1("", [ff_fun_150, FFwRes150Integral](double *x, double *par) 
    {
        return ff_fun_150->Eval(x[0]) / FFwRes150Integral;
    }, ff_fun_150->GetXmin(), ff_fun_150->GetXmax(), 0);
    FF_wCutRes_150->SetLineColor(kOrange);
    FF_wCutRes_150->SetLineStyle(3);
    FF_wCutRes_150->Draw("same");

auto legend = new TLegend(0.53,0.65,0.9,0.9);
	legend->AddEntry(FF_wCutRes_5,"Form Factor with 5 MeV resolution","l");
    legend->AddEntry(FF_wCutRes_50,"Form Factor with 50 MeV resolution","l");
    legend->AddEntry(FF_wCutRes_100,"Form Factor with 100 MeV resolution","l");
    legend->AddEntry(FF_wCutRes_150,"Form Factor with 150 MeV resolution","l");
    legend->AddEntry(trueFF,"True Form Factor","l");
    legend->Draw();

gPad->SetLogy(1);
}

void c()
{
    double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7, a = 0.535;
double trans_bins = 1000, phi_min = 0, phi_max = pi/9, bins = 1000;

// 1D params
double t_min = 0, t_max = 0.1, q_min = 0, q_max = 0.316;
double r_min = 0, r_max = 15;

WoodsSaxon_1D ws_1d(Vo,R,a,r_min,r_max,bins,t_min,t_max);

FormFactor_t_1D ff_1d(A,Vo,R,a0,t_min,t_max,bins,r_min,r_max);
double sigma50 = 0.05, sigma100 = 0.1, sigma150 = 0.15, sigma5 = 0.005;
//FormFactor_transform_resolution_add_wedge_1D ff_transform;

/*

    TH1D *ff_transform_5 = ff_transform.transform_FF_to_G_hist(phi_min,phi_max,sigma5,r_min,r_max,trans_bins);
    ff_transform_5->SetTitle("Histogram of |F(t)|^{2} -> G(r) Transformation with Resolution and #phi = #pi/9 Cut");
    ff_transform_5->GetYaxis()->SetTitle("G(r) [fm^{-3}]");
    ff_transform_5->GetXaxis()->SetTitle("r [fm]");
    ff_transform_5->SetLineStyle(2);
    ff_transform_5->SetLineColor(kBlue);
    ff_transform_5->Scale(197./ff_transform_5->Integral(), "width");
    ff_transform_5->Draw();

TH1D *ff_transform_50 = ff_transform.transform_FF_to_G_hist(phi_min,phi_max,sigma50,r_min,r_max,trans_bins);
    ff_transform_50->SetLineStyle(3);
    ff_transform_50->SetLineColor(kMagenta);
    ff_transform_50->Scale(197./ff_transform_50->Integral(), "width");
    ff_transform_50->Draw("same");

TH1D *ff_transform_100 = ff_transform.transform_FF_to_G_hist(phi_min,phi_max,sigma100,r_min,r_max,trans_bins);
    ff_transform_100->SetLineStyle(4);
    ff_transform_100->SetLineColor(kBlack);
    ff_transform_100->Scale(197./ff_transform_100->Integral(), "width");
    ff_transform_100->Draw("same");

TH1D *ff_transform_150 = ff_transform.transform_FF_to_G_hist(phi_min,phi_max,sigma150,r_min,r_max,trans_bins);
    ff_transform_150->SetLineStyle(2);
    ff_transform_150->SetLineColor(kOrange);
    ff_transform_150->Scale(197./ff_transform_150->Integral(), "width");
    ff_transform_150->Draw("same");*/

TH1D *ff_true_transform = ff_1d.getTransformedFF_hist();
    ff_true_transform->SetLineStyle(3);
    ff_true_transform->SetLineColor(kRed);
    ff_true_transform->Scale(197./ff_true_transform->Integral(), "width");
    ff_true_transform->Draw("same");

TH1D *ws_true = ws_1d.getWoodsSaxon1D_hist();
    ws_true->SetLineStyle(3);
    ws_true->SetLineColor(kCyan);
    ws_true->Scale(197./ws_true->Integral(), "width");
    ws_true->Draw("same");

auto legend = new TLegend(0.53,0.65,0.9,0.9);
	//legend->AddEntry(ff_transform_5,"FF with 5 MeV resolution","l");
    //legend->AddEntry(ff_transform_50,"FF with 50 MeV resolution","l");
    //legend->AddEntry(ff_transform_100,"FF with 100 MeV resolution","l");
    //legend->AddEntry(ff_transform_150,"FF with 150 MeV resolution","l");
    legend->AddEntry(ff_true_transform,"True FF to G","l");
    legend->AddEntry(ws_true,"True Woods-Saxon","l");
    legend->Draw();
}
/*
void unfold_test()
{
    double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7;
    double t_min = 0.0001, t_max = 0.1;
    double bins = 1000, r_min = 0, r_max = 15;
    double r_min = 0, r_max = 15;
    double phi_min = 0, phi_max = pi/9, sigma = 0.1;
    
    // generate truth distribution
    FormFactor_t_1D ff(A,Vo,R,a0,t_min,t_max,bins,r_min,r_max);
    TH1D* histTrue = ff.getFF_hist();
    // Fill the truth histogram using your form factor model
    for (int i = 1; i <= bins; ++i) 
    {
        double t = histTrue->GetXaxis()->GetBinCenter(i);
        double truthValue = calcFF(sqrt(t));
        histTrue->SetBinContent(i, truthValue);
    }

    // generate measured distribution
    FormFactor_resolution_add_wedge_1D ff_wResCut(A,Vo,R,a0,t_min,t_max,bins,phi_min,phi_max,sigma,r_min,r_max);
    TH1D* histMeasured = ff_wResCut.getWedgeRes_hist_1D();
    // Apply the resolution and cut to the truth distribution
    for (int i = 1; i <= bins; ++i) 
    {
        double t = histMeasured->GetXaxis()->GetBinCenter(i);
        double measuredValue = cutFFt2->Eval(t); // Using your TF1 object for the cut form factor
        histMeasured->SetBinContent(i, measuredValue);
    }

    // set up response matrix
    RooUnfoldResponse response(histMeasured, histTrue);
    // Fill the response matrix using your resolution model
    for (int i = 1; i <= bins; ++i) 
    {
        double t_true = histTrue->GetXaxis()->GetBinCenter(i);
        double t_measured = histMeasured->GetXaxis()->GetBinCenter(i);
        double weight = cutFFt2->Eval(t_true); // Using your TF1 object for the cut form factor
        response.Fill(t_measured, t_true, weight);
    }

    // perform unfolding
    RooUnfoldBayes unfold(&response, histMeasured, 4); // Using Bayesian unfolding with 4 iterations
    TH1D* histUnfolded = (TH1D*)unfold.Hreco();
    histUnfolded->SetTitle("Unfolded Distribution");

    // compare unfolded with truth
    TCanvas* c = new TCanvas("c", "Unfolding Example", 800, 600);
    histTrue->SetLineColor(kRed);
    histMeasured->SetLineColor(kBlue);
    histUnfolded->SetLineColor(kGreen);
    histTrue->Draw("HIST");
    histMeasured->Draw("HIST SAME");
    histUnfolded->Draw("HIST SAME");


}
*/


// FormFactor_q_2DwResCut tests
void hist_2d_wResCut_q()
{
    double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7;
    double qy_min = -.5, qy_max = .5, bins = 1000;
    double qx_prime_min = -.5, qx_prime_max = .5, phi_min = 0, sigma = 0.025;

    FormFactor_q_2DwResCut ff_wResCut25(A,Vo,R,a0,qy_min,qy_max,qx_prime_min,qx_prime_max,bins,sigma);
    TCanvas *c1 = new TCanvas("", "", 1000, 1000);
    c1->SetLeftMargin(0.13); 
    c1->SetBottomMargin(0.13);

TH2D *ff_wRes25 = ff_wResCut25.getSmeared_hist();
    ff_wRes25->SetTitle("|F(q_{x},q_{y})|^{2} with 25 MeV/c Resolution");
    ff_wRes25->GetYaxis()->SetTitleSize(0.04);
    ff_wRes25->GetXaxis()->SetTitleSize(0.04);
    ff_wRes25->GetXaxis()->SetTitle("q_{x} = #sqrt{t_{x}} [GeV/c]");
    ff_wRes25->GetYaxis()->SetTitle("q_{y} = #sqrt{t_{y}} [GeV/c]");
    ff_wRes25->GetYaxis()->SetTitleOffset(1.4);
    ff_wRes25->GetXaxis()->SetTitleOffset(1.1);
    ff_wRes25->Draw("col");

TGaxis *xAxis = new TGaxis(ff_wRes25->GetXaxis()->GetXmin(), 0, ff_wRes25->GetXaxis()->GetXmax(), 0, ff_wRes25->GetXaxis()->GetXmin(), ff_wRes25->GetXaxis()->GetXmax(), 510, "U");
    xAxis->SetLabelSize(0.02);
    xAxis->SetLabelOffset(0.01);
    xAxis->Draw();

TGaxis *yAxis = new TGaxis(0, ff_wRes25->GetYaxis()->GetXmin(), 0, ff_wRes25->GetYaxis()->GetXmax(), ff_wRes25->GetYaxis()->GetXmin(), ff_wRes25->GetYaxis()->GetXmax(), 510, "U");
    yAxis->SetLabelSize(0.01);
    yAxis->SetLabelOffset(0.02);
    yAxis->Draw();

double length = .5; 
double theta = pi/12; 
double x_end = length * sin(theta);
double y_end = length * cos(theta);

TArrow *arrow = new TArrow(0, 0, x_end, y_end, 0.02, "|>");
    arrow->SetLineColor(kBlack);
    arrow->SetLineWidth(2);
    arrow->Draw();

double theta2 = -pi/12; 
double x_end2 = length * sin(theta2);
double y_end2 = length * cos(theta2);

TArrow *arrow2 = new TArrow(0, 0, x_end2, y_end2, 0.02, "|>");
    arrow2->SetLineColor(kBlack);
    arrow2->SetLineWidth(2);
    arrow2->Draw();

double theta3 = -13*pi/12; 
double x_end3 = length * sin(theta3);
double y_end3 = length * cos(theta3);

TArrow *arrow3 = new TArrow(0, 0, x_end3, y_end3, 0.02, "|>");
    arrow3->SetLineColor(kBlack);
    arrow3->SetLineWidth(2);
    arrow3->Draw();

double theta4 = 13*pi/12; 
double x_end4 = length * sin(theta4);
double y_end4 = length * cos(theta4);

TArrow *arrow4 = new TArrow(0, 0, x_end4, y_end4, 0.02, "|>");
    arrow4->SetLineColor(kBlack);
    arrow4->SetLineWidth(2);
    arrow4->Draw();

    auto el3 = new TEllipse(0, 0, 0.3, 0.3, -15, 15);
   el3->SetFillStyle(0);
   el3->SetLineColor(kBlack);
   el3->SetLineWidth(2);
   el3->SetTheta(270);
   el3->Draw();

    auto el = new TEllipse(0, 0, -0.3, -0.3, -15, 15);
    el->SetFillStyle(0);
    el->SetLineColor(kBlack);
    el->SetLineWidth(2);
    el->SetTheta(270);
    el->Draw();


gStyle->SetOptStat(0);
gPad->SetLogz(1);
c1->Draw();
}

// ATHENA plot scan
void ATHENA_data()
{
    double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7, a = 0.535;
    double bins = 1000, phi_min = 0, phi_max = pi/2;

    //sigma = 25 MeV, 50 MeV, 100 MeV
    double sigma25 = 0.025, sigma50 = 0.05, sigma100 = 0.1;

    // 1D params
    double t_min = 0, t_max = 0.17;
    double r_min = 0, r_max = 15;

    FormFactor_t_1D ff(A,Vo,R,a0,t_min,t_max,bins,r_min,r_max);
    TH1D *ff_true = ff.getFF_hist();
    ff_true->GetYaxis()->SetTitle("|F(t)|^{2}");
    ff_true->GetXaxis()->SetTitle("|t| [GeV^{2}/c^{2}]");
    ff_true->GetYaxis()->SetTitleSize(0.047);
    ff_true->GetXaxis()->SetTitleSize(0.047);
    ff_true->Scale(1/(pi)*197./ff_true->Integral(), "width");
    ff_true->GetYaxis()->SetTitleOffset(0.9);
    ff_true->GetXaxis()->SetTitleOffset(0.9);
    ff_true->SetLineStyle(1);
    ff_true->SetLineWidth(2);
    ff_true->SetLineColor(kBlack);
    ff_true->Draw();

    const int incoherent_bins = 91;
    double incoherent_min = 0.001, incoherent_max = 0.18; 
    double incoherent_xvals[incoherent_bins] = {0.0010133167966835457,0.0030211660950511393,0.0030211660950511393,0.003221951024887905,0.0052298003232554985,0.005631370182929016,
        0.008241574270806892,0.008241574270806892,0.008241574270806892,0.010450208499011251,0.011855703007868562,0.013863552306236163,0.016473756394114032,0.018481605692481633,
        0.020690239920685992,0.022898874148890344,0.025308293306931462,0.02912320697382989,0.031733411061707774,0.034945969939095926,0.037355389097137044,0.03936323839550463,
        0.04117030276403547,0.042977367132566315,0.04538678629060742,0.04759542051881178,0.04960326981717938,0.051410334185710206,0.05361896841391458,0.055626817712282166,
        0.057835451940486525, 0.05964251630901737,0.06185115053722171,0.06405978476542606, 0.0658668491339569,0.06827626829199804,0.0712880422395494,0.0712880422395494,
        0.07309510660808025,0.07450060111693757,0.07691002027497867,0.07911865450318303,0.08092571887171388,0.08293356817008146,0.08735083662649018,0.08735083662649018,
        0.08554377225795934,0.0893586859248578,0.09196889001273567,0.09477987903045029,0.09638615846914436,0.0981932228376752,0.10000028720620605,0.1022089214344104,
        0.10421677073277802,0.10602383510130883,0.10863403918918671,0.11044110355771755,0.11305130764559543,0.11485837201412627,0.11706700624233063,0.11927564047053499,
        0.1210827048390658,0.12349212399710696,0.12529918836563775,0.12811017738335242,0.12991724175188324,0.1321258759800876,0.1341337252784552,0.13634235950665954,0.13835020880502713,
        0.14035805810339474,0.14236590740176236,0.14417297177029317,0.14618082106866076,0.14818867036702837,0.15039730459523273,0.15200358403392678,0.15461378812180468,
        0.1566216374201723,0.1584287017887031,0.1604365510870707,0.16264518531527508,0.16465303461364264,0.16706245377168377,0.1690703030700514,0.17067658250874546,
        0.1728852167369498,0.17509385096515417,0.17710170026352176,0.17910954956188935};
    double incoherent_yvals[incoherent_bins] = {26.4259758563496,26.4259758563496,32.1883134302757,39.9882250899722,39.9882250899722,50.6678823310276,50.6678823310276,
        50.6678823310276,62.9457858323086,61.7163084604715,69.4704726439587,75.1739081002721,75.1739081002721,82.9661089645002,82.9661089645002,82.9661089645002,82.9661089645002,
        82.9661089645002,82.9661089645002,84.6189128366981,91.5660155319823,88.0239551889403,95.2506074330171,91.5660155319823,91.5660155319823,91.5660155319823,
        95.2506074330171,93.3901418761231,89.7775186113819,91.5660155319823,91.5660155319823,89.7775186113819,89.7775186113819,89.7775186113819,89.7775186113819,91.5660155319823,
        88.0239551889403,88.0239551889403,88.0239551889403,89.7775186113819,88.0239551889403,84.6189128366981,84.6189128366981,84.6189128366981,84.6189128366981,84.6189128366981,
        84.6189128366981,84.6189128366981,84.6189128366981,84.6189128366981,84.6189128366981,84.6189128366981,82.9661089645002,79.7567199496749,78.1988859956862,
        81.3455881901156,79.7567199496749,78.1988859956862,78.1988859956862,78.1988859956862,76.671480156466,78.1988859956862,76.671480156466,76.671480156466,
        75.1739081002721,73.7055871040409,76.671480156466,75.1739081002721,75.1739081002721,75.1739081002721,72.2659458266437,72.2659458266437,69.4704726439587,
        68.1135529868728,69.4704726439587,69.4704726439587,69.4704726439587,68.1135529868728,68.1135529868728,64.1997561565852,66.783137121768,66.783137121768,66.783137121768,
        64.1997561565852,64.1997561565852,64.1997561565852,61.7163084604715,61.7163084604715,65.4787073680391,62.9457858323086,60.5108456368344};

    TGraph *incoherent_hist = new TGraph(incoherent_bins, incoherent_xvals, incoherent_yvals);
    incoherent_hist->SetMarkerColor(kBlack);
    incoherent_hist->GetYaxis()->SetRangeUser(1e-2,1e5);
    incoherent_hist->SetMarkerSize(0.3);
    incoherent_hist->SetMarkerStyle(25);
    incoherent_hist->Draw("P same");

    const int detector_bins = 36;
    double detector_min = 0.003, detector_max = 0.18;
    double detector_xvals[detector_bins] = {0.0030211660950511393,0.008040789340970134,0.013261197516725887,0.018280820762644874,0.02309965907872711,0.02791849739480933,0.03333969050040186,
         0.03815852881648408,0.042977367132566315,0.047796205448648536,0.05301661362440429,0.058036236870323284,0.06305586011624228,0.06807548336216127,0.0728943216782435,
         0.07811472985399925,0.08273278324024472,0.0877524064861637,0.09317359959175622,0.0981932228376752,0.10301206115375744,0.10803168439967645,0.11285052271575868,
         0.11807093089151444,0.12288976920759664,0.12831096231318914,0.13312980062927138,0.1377478540155169,0.14296826219127262,0.1479878854371916,0.15320829361294738,
         0.1580271319290296,0.16304675517494857,0.1678655934910308,0.17308600166678656,0.17810562491270557};
    double detector_yvals[detector_bins] = {6.385976005504609,6.910256330465505,6.775282977133568,6.642945966831217,5.347204743387561,3.9776451066954124,2.6809681618761285,2.074551939390839,
         1.9943017711225621,2.0340311224297665,1.8797094233711418,1.6053028194928454,1.344177122426765,1.217931619705514,1.0608546254453768,0.9424440169171399,0.8372501789771866,
         0.7150254482006665,0.6106435142231171,0.4915343150578518,0.4632908104997826,0.48193350842950883,0.4915343150578518,0.45424166510334535,0.4115791570074048,
         0.4035400602899766,0.45424166510334535,0.43667017442776856,0.3803526953736092,0.35149538166571154,0.35849767374771574,0.365639461529973,0.3803526953736092,
         0.372923523966886,0.365639461529973,0.3879298665289226};
    
     TGraph *detector_hist = new TGraph(detector_bins, detector_xvals, detector_yvals);
     detector_hist->SetMarkerColor(kRed);
     detector_hist->GetYaxis()->SetRangeUser(1e-2,1e5);
     detector_hist->SetMarkerSize(.5);
     detector_hist->SetMarkerStyle(27);
     detector_hist->Draw("P SAME ");

    // Coherent from plot scan
    const int coherent_bins = 149;
    double coherent_min = 0.0006, coherent_max = 0.18; 
    double coherent_xvals[coherent_bins] = {0.000611746937010025,0.00261959623537762,0.00282038116521438,0.0042258756740717,0.0042258756740717,0.00663529483211282,
        0.00663529483211282,0.00683607976194958,0.00844235920064365,0.00804078934097013,0.0100486386393377,0.010650993428848,0.010650993428848,0.0122572728675421,
        0.0122572728675421,0.0124580577973788,0.0146666920255832,0.0146666920255832,0.0146666920255832,0.0144659070957464,0.0144659070957464,0.0162729714642773,0.016473756394114,
        0.0158714016046038,0.0162729714642773,0.0160721865344405,0.0182808207626449,0.0182808207626449,0.0188831755521552,0.0220957344295433,0.0224973042892168,
        0.0222965193593801,0.0222965193593801,0.0243043686577477,0.0243043686577477,0.026513002885952,0.0269145727456255,0.0293239919036667,0.0323357658512181,
        0.0343436151495857,0.0353475397987694,0.0381585288164841,0.0387608836059944,0.0403671630446884,0.0421742274132193,0.0431781520624031,0.0443828616414236,0.045788356150281,
        0.0467922807994647,0.0477962054486485,0.0483985602381588,0.0479969903784853,0.0504064095365264,0.0506071944663632,0.0504064095365264,0.0512095492558735,
        0.0522134739050572,0.0522134739050572,0.0522134739050572,0.0522134739050572,0.0532173985542411,0.0544221081332616,0.0544221081332616,0.0556268177122822,
        0.0562291725017924,0.0572330971509763,0.0584378067299968,0.0584378067299968,0.0586385916598336,0.0596425163090174,0.0602448710985276,0.0602448710985276,
        0.0620519354670585,0.0628550751864055,0.0642605696952628,0.0664692039234672,0.0686778381516716,0.0706856874500391,0.073295891537917,0.0757053106959582,
        0.0779139449241625,0.0815280736612242,0.0827327832402447,0.0849414174684491,0.0873508366264902,0.0889571160651843,0.0913665352232254,0.0927720297320827,
        0.0943783091707768,0.0963861584691444,0.0965869433989811,0.0985947926973488,0.100201072136043,0.10040185706588,0.102610491294084,0.102610491294084,
        0.103212846083594,0.104216770732778,0.104216770732778,0.106626189890819,0.106626189890819,0.106626189890819,0.107027759750493,0.107027759750493,0.108232469329513,
        0.108232469329513,0.108232469329513,0.112649737785922,0.112248167926248,0.112248167926248,0.114256017224616,0.114858372014126,0.114858372014126,
        0.114858372014126,0.116263866522984,0.11646465145282,0.117870145961678,0.118472500751188,0.121082704839066,0.124496048646291,0.126905467804332,
        0.127708607523679,0.130318811611557,0.132929015699435,0.135940789646986,0.139354133454211,0.142365907401762,0.14437375670013,0.147385530647681,
        0.149995734735559,0.153007508683111,0.155818497700825,0.157625562069356,0.160034981227397,0.161038905876581,0.164251464753969,0.166259314052337,
        0.16666088391201,0.16666088391201,0.167464023631357,0.170475797578909,0.170676582508745,0.173286786596623,0.174692281105481,0.174692281105481,0.176900915333685,
        0.178507194772379,0.178306409842542,0.1797119043514};
    double coherent_yvals[coherent_bins] = {18825.5808273721,18825.5808273721,11725.9255897844,9254.36253303648,6617.77945297608,6749.61511023115,5020.8613291811,
        3734.88681579103,3006.37768200003,2149.855744337,1836.01222879312,1339.08476545981,1036.19316957847,902.555035945294,596.448017766563,394.159052611198,
        336.618325480646,236.013484786255,175.564229896955,120.689191843757,89.7775186113819,81.3455881901156,58.1700965564148,41.5973405401449,29.7461899231103,
        23.4763578718976,22.5682188015885,17.8113426454577,12.4880813915112,14.914070425225,18.5280662953482,25.4037363217018,37.6905063285368,43.2712063643629,
        58.1700965564148,64.1997561565852,89.7775186113819,95.2506074330171,95.2506074330171,91.5660155319823,72.2659458266437,61.7163084604715,47.7565117069936,
        41.5973405401449,34.1506031120067,26.9524191928772,20.0491964420711,16.4599981998138,12.7368618819768,10.0522160478075,7.62654349327295,5.24276141731504,
        4.65757419665132,3.33062366092517,2.24486995561805,1.7370960932444,1.2921801402308,0.999898042087226,0.68736602379934,0.481933508429509,0.380352695373609,
        0.300182847524764,0.223298132738852,0.176231876749734,0.246444247486369,0.31848283412103,0.411579157007405,0.587021910040612,0.888291293349121,1.12552729274239,
        1.39826644923732,1.91715593077436,2.07455193939084,2.95886567919123,3.74909008162436,4.75035975407174,5.7862032986717,6.13894644732955,6.38597600550461,
        6.38597600550461,6.38597600550461,6.13894644732955,4.94151295621955,4.38995052517442,3.97764510669541,3.0779296403613,2.68096816187613,2.03403112242977,
        1.70316657260834,1.45453231604845,1.08198837095495,0.888291293349121,0.758615400877701,0.511313524034245,0.491534315057852,0.344629860608996,0.251353770321779,
        0.194499307790862,0.156561204380984,0.131094262486242,0.0883586983008236,0.0607409601406383,0.0417555295602655,0.0298593106316738,0.025002262348988,0.0165225933309593,
        0.0115844995193824,0.0111363747855147,0.0155732070404027,0.0217777133218168,0.0255003432707026,0.0370948554716042,0.0529071808931029,0.0784962657076955,0.0975175772611055,
        0.150504926513317,0.186975465229049,0.256361097892814,0.324827468912678,0.351495381665712,0.419778404060304,0.553291699383871,0.660776552586587,0.789141881098933,
        0.758615400877701,0.673940169133192,0.701059340744442,0.587021910040612,0.68736602379934,0.531888643073776,0.50132638352274,0.445369270530689,0.33789843912183,
        0.300182847524764,0.214660261405754,0.202325907758978,0.169414675649648,0.123561601415304,0.0901189302861088,0.0683725795511678,0.063185160082014,0.045183604233252,
        0.0409399471020499,0.0342804732040046,0.0255003432707026,0.0255003432707026,0.0205263684479667,0.0161998687139596,0.0132997754119192};
    
    const int methodL_bins = 46;
    double methodL_min = 0.0015, methodL_max = 0.18;
    double methodL_xvals[methodL_bins] ={0.0016156715861938217, 0.004426660603908457,0.006635294832112816,0.009245498919990693,0.011855703007868562,0.01426512216590968,
        0.016674541323950798,0.019083960481991916,0.021694164569869785,0.02430436865774767,0.02651300288595202,0.029323991903666656,0.03153262613187101,0.03434361514958565,
        0.03675303430762677,0.03916245346566787,0.04177265755354575,0.04418207671158688,0.04679228079946474,0.04940248488734262,0.05161111911554698,0.054221323203424855,
        0.05683152729130273,0.05924094644934385,0.06144958067754821,0.06426056969526284,0.0664692039234672,0.06928019294118182,0.07148882716938618,0.07409903125726405,
        0.07751237506448896,0.08152807366122417,0.0853429873281226,0.0893586859248578,0.09357516945142974,0.09739008311832817,0.10140578171506337,0.1052206953819618,
        0.10943717890853377,0.1162638665229836,0.12650389794465833,0.1365431444364963, 0.14638160599849753,0.15642085249033552,0.16646009898217348,0.17649934547401147};
    double methodL_yvals[methodL_bins] = {11959.522843672647,6884.077122844144,3661.9358030159865,1986.7464450059745,1036.193169578465,529.8736059242736,287.4775964137653,
        179.06172113301182,128.0467427873466,105.12387176520791,95.25060743301714,84.61891283669812,78.19888599568615,66.78313712176801,58.17009655641482,45.909140220349855,
        36.23251946518077,28.036977310803344,22.127409079820325,16.13849621108551,12.005003280046312,9.855872924845716,7.778475176878387,7.047918546560998,
        6.775282977133568,6.261243064418155,6.01903856717346,5.901472766370377,5.786203298671697,5.901472766370377,5.673185311359277,4.47740463526015,4.389950525174419,
        3.3306236609251663,2.577259914313217,2.0340311224297665,1.5739475275473165,1.2921801402307969,1.0608546254453768,0.9240358966958786,1.0608546254453768,0.8882912933491209,
        0.905987329810412,0.6739401691331922,0.6106435142231171,0.4725202278498894};

    TGraph *methodL_hist = new TGraph(methodL_bins, methodL_xvals, methodL_yvals);
    methodL_hist->SetMarkerColor(kBlue);
    methodL_hist->GetYaxis()->SetRangeUser(1e-2,1e5);
    methodL_hist->SetMarkerSize(0.9);
    methodL_hist->SetMarkerStyle(4);
    methodL_hist->SetLineStyle(1);
    methodL_hist->SetLineWidth(1);
    methodL_hist->Draw("P SAME");
    methodL_hist->Draw("L SAME");

    FormFactor_resolution_add_wedge_1D ff_wResCut_25(A,Vo,R,a0,t_min,t_max,bins,phi_min,phi_max,sigma25,r_min,r_max);
    TH1D *ff_25 = ff_wResCut_25.getWedgeRes_hist_1D();
    ff_25->Scale((1/pi)*197./ff_25->Integral(), "width");
    ff_25->GetYaxis()->SetRangeUser(1e-2,1e5);
    ff_25->SetLineStyle(2);
    ff_25->SetLineWidth(3);
    ff_25->SetLineColor(kGreen+2);
    ff_25->Draw("same");

    FormFactor_resolution_add_wedge_1D ff_wResCut_50(A,Vo,R,a0,t_min,t_max,bins,phi_min,phi_max,sigma50,r_min,r_max);
    TH1D *ff_50 = ff_wResCut_50.getWedgeRes_hist_1D();
    ff_50->Scale((1/pi)*197./ff_50->Integral(), "width");
    ff_50->GetYaxis()->SetRangeUser(1e-2,1e5);
    ff_50->SetLineStyle(3);
    ff_50->SetLineWidth(3);
    ff_50->SetLineColor(kOrange+2);
    ff_50->Draw("same");

    FormFactor_resolution_add_wedge_1D ff_wResCut_100(A,Vo,R,a0,t_min,t_max,bins,phi_min,phi_max,sigma100,r_min,r_max);
    TH1D *ff_100 = ff_wResCut_100.getWedgeRes_hist_1D();
    ff_100->Scale((1/pi)*197./ff_100->Integral(), "width");
    ff_100->GetYaxis()->SetRangeUser(1e-2,1e5);
    ff_100->SetLineStyle(9);
    ff_100->SetLineColor(kMagenta);
    ff_100->Draw("same");
   
    
    // Upper right legend
    auto legend = new TLegend(0.55,0.74,0.77,0.84);
    legend->AddEntry(incoherent_hist,"Truth: incoherent","p");
    legend->AddEntry(detector_hist,"Reco.: incoherent","p");
    legend->AddEntry(methodL_hist,"Reco.: coherent","p");
    legend->SetTextSize(0.04);
    legend->SetBorderSize(0);
    legend->Draw("same");

    // Lower legend
    auto legend3 = new TLegend(0.15,0.14,0.27,0.3);
    legend3->SetHeader("Applied res.");
    legend3->AddEntry(ff_25,"   25 MeV/c","l");
    legend3->AddEntry(ff_50,"   50 MeV/c","l");
    legend3->AddEntry(ff_100,"   100 MeV/c","l");
    legend3->AddEntry(ff_true,"Truth: coherent","l");
    legend3->SetTextSize(0.04);
    legend3->SetBorderSize(0);
    legend3->Draw("same");
   
gPad->SetLogy(1);
gStyle->SetOptStat(0);
}




Double_t Reconstruct( Double_t xt, TRandom3& R )
{
   // apply some Gaussian smearing + bias and efficiency corrections to fake reconstruction
   const Double_t cutdummy = -99999.0;
   Double_t xeff = 0.3 + (1.0 - 0.3)/20.0*(xt + 10.0);  // efficiency
   Double_t x    = R.Rndm();
   if (x > xeff) return cutdummy;
   else {
     Double_t xsmear= R.Gaus(-2.5,0.2); // bias and smear
     return xt+xsmear;
   }
}

void x()
{
    

gROOT->SetStyle("Plain");
gStyle->SetOptStat(0);

TRandom3 R;

const Double_t cutdummy= -99999.0;

Int_t nbins = 40;
TH1D *xini = new TH1D("xini", "MC truth", nbins, -10.0, 10.0);
TH1D *bini = new TH1D("bini", "MC reco", nbins, -10.0, 10.0);
TH2D *Adet = new TH2D("Adet", "detector response", nbins, -10.0, 10.0, nbins, -10.0, 10.0);

TH1D *data = new TH1D("data", "data", nbins, -10.0, 10.0);

TH1D *datatrue = new TH1D("datatrue", "data truth", nbins, -10.0, 10.0);

TH2D *statcov = new TH2D("statcov", "covariance matrix", nbins, -10.0, 10.0, nbins, -10.0, 10.0);

for (Int_t i= 0; i<100000; i++) {
   Double_t xt = R.BreitWigner(0.3, 2.5);
   xini->Fill(xt);
   Double_t x = Reconstruct( xt, R );
   if (x != cutdummy) {
      Adet->Fill(x, xt);
      bini->Fill(x);
   }
}

for (Int_t i=0; i<10000; i++) {
   Double_t xt = R.Gaus(0.0, 2.0);
   datatrue->Fill(xt);
   Double_t x = Reconstruct( xt, R );
   if (x != cutdummy)
   data->Fill(x);
}

cout << "Created toy distributions and errors for: " << endl;
cout << "... \"true MC\"   and \"reconstructed (smeared) MC\"" << endl;
cout << "... \"true data\" and \"reconstructed (smeared) data\"" << endl;
cout << "... the \"detector response matrix\"" << endl;

for (int i=1; i<=data->GetNbinsX(); i++) {
    statcov->SetBinContent(i,i,data->GetBinError(i)*data->GetBinError(i));
}

TSVDUnfold *tsvdunf = new TSVDUnfold( data, statcov, bini, xini, Adet );

tsvdunf->SetNormalize( kFALSE ); // no normalisation here

TH1D* unfres = tsvdunf->Unfold( 13 );

TH1D* ddist = tsvdunf->GetD();

TH1D* svdist = tsvdunf->GetSV();

TH2D* ustatcov = tsvdunf->GetUnfoldCovMatrix( statcov, 100 );

TH2D* uadetcov = tsvdunf->GetAdetCovMatrix( 100 );

ustatcov->Add( uadetcov );

TH2D* utaucov = tsvdunf->GetXtau();
utaucov->Add( uadetcov );

TH2D* uinvcov = tsvdunf->GetXinv();

for (int i=1; i<=unfres->GetNbinsX(); i++) {
   unfres->SetBinError(i, TMath::Sqrt(utaucov->GetBinContent(i,i)));
}

xini->Scale(0.7*datatrue->Integral()/xini->Integral());

TLegend *leg = new TLegend(0.58,0.60,0.99,0.88);
leg->SetBorderSize(0);
leg->SetFillColor(0);
leg->SetFillStyle(0);
leg->AddEntry(unfres,"Unfolded Data","p");
leg->AddEntry(datatrue,"True Data","l");
leg->AddEntry(data,"Reconstructed Data","l");
leg->AddEntry(xini,"True MC","l");

TCanvas *c1 = new TCanvas( "c1", "Unfolding toy example with TSVDUnfold", 1000, 900 );

c1->Divide(1,2);
TVirtualPad * c11 = c1->cd(1);

TH1D* frame = new TH1D( *unfres );
frame->SetTitle( "Unfolding toy example with TSVDUnfold" );
frame->GetXaxis()->SetTitle( "x variable" );
frame->GetYaxis()->SetTitle( "Events" );
frame->GetXaxis()->SetTitleOffset( 1.25 );
frame->GetYaxis()->SetTitleOffset( 1.29 );
frame->Draw();

data->SetLineStyle(kDashed);
data->SetLineColor(4);
data->SetLineWidth(2);
unfres->SetMarkerStyle(20);
datatrue->SetLineColor(2);
datatrue->SetLineWidth(2);
xini->SetLineStyle(kDashed);
xini->SetLineColor(8);
xini->SetLineWidth(2);

unfres->Draw("same");
datatrue->Draw("same");
data->Draw("same");
xini->Draw("same");

leg->Draw();

TVirtualPad * c12 = c1->cd(2);
c12->Divide(2,1);
TVirtualPad * c2 = c12->cd(1);
c2->SetRightMargin   ( 0.15         );

TH2D* covframe = new TH2D( *ustatcov );
covframe->SetTitle( "TSVDUnfold covariance matrix" );
covframe->GetXaxis()->SetTitle( "x variable" );
covframe->GetYaxis()->SetTitle( "x variable" );
covframe->GetXaxis()->SetTitleOffset( 1.25 );
covframe->GetYaxis()->SetTitleOffset( 1.29 );
covframe->Draw();

ustatcov->SetLineWidth( 2 );
ustatcov->Draw( "colzsame" );

TVirtualPad * c3 = c12->cd(2);
c3->SetLogy();

TLine *line = new TLine( 0.,1.,40.,1. );
line->SetLineStyle(kDashed);

TH1D* dframe = new TH1D( *ddist );
dframe->SetTitle( "TSVDUnfold |d_{i}|" );
dframe->GetXaxis()->SetTitle( "i" );
dframe->GetYaxis()->SetTitle( "|d_{i}|" );
dframe->GetXaxis()->SetTitleOffset( 1.25 );
dframe->GetYaxis()->SetTitleOffset( 1.29 );
dframe->SetMinimum( 0.001 );
dframe->Draw();

ddist->SetLineWidth( 2 );
ddist->Draw( "same" );
line->Draw();

gROOT->GetListOfCanvases()->Draw();
}

void f()
{
double bins = 100;

double A = 197, Vo = 2.12, RR = 6.38, a0 = 0.7, a = 0.535;
double phi_min = 0, phi_max = pi/9;
double tx_min = 0, tx_max = 0.1, ty_min = 0, ty_max = 0.1;
double x_min = 0, x_max = 15, y_min = 0, y_max = 15;

//sigma = 5 MeV, 50 MeV, 100 MeV, 150 MeV respectively
double sigma5 = 0.005, sigma50 = 0.05, sigma100 = 0.1, sigma150 = 0.15;

// 1D params
double t_min = 0, t_max = 0.1, q_min = 0, q_max = 0.316;
double r_min = 0, r_max = 15;

FormFactor_resolution_add_wedge_1D ff_wResCut_100(A,Vo,RR,a0,t_min,t_max,bins,phi_min,phi_max,sigma100,r_min,r_max);
FormFactor_t_1D ff(A,Vo,RR,a0,t_min,t_max,bins,r_min,r_max);

ROOT::Math::IntegratorOneDimOptions::SetDefaultAbsTolerance(1e-4);
ROOT::Math::IntegratorOneDimOptions::SetDefaultRelTolerance(1e-2);
ROOT::Math::IntegratorOneDimOptions::SetDefaultIntegrator("Gauss");

// truth
TH1D* histTrue = ff.getFF_hist();

// measured (100 MeV resolution with pi/9 wedge cut)
TH1D* histMeasured = ff_wResCut_100.getWedgeRes_hist_1D();

TH2D* Adet = new TH2D("", "Response Matrix", bins, tx_min, tx_max, bins, ty_min, ty_max);
for (int i = 1; i <= bins; ++i) 
{
    double t_true = histTrue->GetXaxis()->GetBinCenter(i);
    for (int j = 1; j <= bins; ++j) 
    {
        double t_measured = histMeasured->GetXaxis()->GetBinCenter(j);
        double response_value = ff_wResCut_100.getWedgeRes_fun_1D()->Eval(t_true); 
        Adet->SetBinContent(i, j, response_value);
    }
}

// guesses
TH1D* bini = (TH1D*)histMeasured->Clone("bini");
TH1D* xini = (TH1D*)histTrue->Clone("xini");

TH2D* Bcov = new TH2D("", "Covariance Matrix", bins, tx_min, tx_max, bins, ty_min, ty_max);
for (int i = 1; i <= bins; ++i) 
{
    for (int j = 1; j <= bins; ++j) 
    {
        double covariance = 0.0;
        if (i == j) 
        {
            // Diagonal elements: variance of each bin
            covariance = histMeasured->GetBinError(i) * histMeasured->GetBinError(i);
        } else 
        {
            // Off-diagonal elements: covariance between different bins
            covariance = 0.0;
        }
        Bcov->SetBinContent(i, j, covariance);
    }
}


TSVDUnfold unfold(histMeasured, Bcov, bini, xini, Adet);
TH1D* histUnfolded = (TH1D*)unfold.Unfold(3); 
histUnfolded->SetTitle("Unfolded Distribution");



// compare
TCanvas* c1 = new TCanvas("", "", 800, 600);
histTrue->SetLineColor(kRed);
histTrue->SetLineStyle(2);
histTrue->Scale(197./histTrue->Integral(), "width");
//histTrue->GetYaxis()->SetRangeUser(-1e5,1e8);
//histMeasured->GetYaxis()->SetRangeUser(-1e5,1e8);
histMeasured->SetLineColor(kBlue);
histMeasured->SetLineStyle(3);
histMeasured->Scale(197./histMeasured->Integral(), "width");
//histUnfolded->SetLineColor(kGreen);
//histUnfolded->SetMarkerStyle(4);
//histUnfolded->SetMarkerSize(0.7);
//histUnfolded->SetMarkerColorAlpha(kGreen,0.35);
histUnfolded->SetLineStyle(10);
histUnfolded->SetLineColor(kGreen);
//histUnfolded->GetYaxis()->SetRangeUser(-1e5,1e8);
histUnfolded->Scale(197./histUnfolded->Integral(), "width");
histTrue->Draw();
histMeasured->Draw("same");
histUnfolded->Draw("same");

auto legend = new TLegend(0.53,0.65,0.9,0.9);
    legend->AddEntry(histUnfolded,"Unfolded FF","l");
    legend->AddEntry(histMeasured,"Measured FF (with res and wedge)","l");
    legend->AddEntry(histTrue,"True FF","l");
    legend->Draw();

gPad->SetLogy(0);
gStyle->SetOptStat(0);
c1->Draw();
}

Double_t Reconstruct2( Double_t xt, TRandom3& R )
{
   // apply some Gaussian smearing + bias and efficiency corrections to fake reconstruction
   const Double_t cutdummy = -99999.0;
   Double_t xeff = 0.3 + (1.0 - 0.3)/20.0*(xt + 10.0);  // efficiency
   Double_t x    = R.Rndm();
   if (x > xeff) return cutdummy;
   else {
     Double_t xsmear= R.Gaus(-2.5,0.2); // bias and smear
     return xt+xsmear;
   }
}

void y()
{
    double bins = 40;

double A = 197, Vo = 2.12, RR = 6.38, a0 = 0.7, a = 0.535;
double phi_min = 0, phi_max = pi/9;
double tx_min = 0, tx_max = 0.1, ty_min = 0, ty_max = 0.1;
double x_min = 0, x_max = 15, y_min = 0, y_max = 15;

//sigma = 5 MeV, 50 MeV, 100 MeV, 150 MeV respectively
double sigma5 = 0.005, sigma50 = 0.05, sigma100 = 0.1, sigma150 = 0.15;

// 1D params
double t_min = 0, t_max = 0.1, q_min = 0, q_max = 0.316;
double r_min = 0, r_max = 15;

FormFactor_resolution_add_wedge_1D ff_wResCut_100(A,Vo,RR,a0,t_min,t_max,bins,phi_min,phi_max,sigma100,r_min,r_max);
FormFactor_t_1D ff(A,Vo,RR,a0,t_min,t_max,bins,r_min,r_max);

ROOT::Math::IntegratorOneDimOptions::SetDefaultAbsTolerance(1e-4);
ROOT::Math::IntegratorOneDimOptions::SetDefaultRelTolerance(1e-2);
ROOT::Math::IntegratorOneDimOptions::SetDefaultIntegrator("Gauss");


gROOT->SetStyle("Plain");
gStyle->SetOptStat(0);

TRandom3 R;

const Double_t cutdummy= -99999.0;

Int_t nbins = 40;
// truth
//TH1D* xini = ff.getFF_hist();
TH1D* datatrue = ff.getFF_hist();

// measured (100 MeV resolution with pi/9 wedge cut)
//TH1D* bini = ff_wResCut_100.getWedgeRes_hist_1D();
TH1D* data = ff_wResCut_100.getWedgeRes_hist_1D();
TH1D *xini = new TH1D("", "Truth", nbins, t_min, t_max);
TH1D *bini = new TH1D("", "Measured", nbins, t_min, t_max);
TH2D *Adet = new TH2D("", "Response Matrix", nbins, t_min, t_max, nbins, t_min, t_max);

//TH1D *data = new TH1D("data", "data", nbins, t_min, t_max);

//TH1D *datatrue = new TH1D("datatrue", "data truth", nbins, t_min, t_max);

TH2D *statcov = new TH2D("", "Covariance Matrix", nbins, t_min, t_max, nbins, t_min, t_max);

for (Int_t i= 0; i<100000; i++) 
{
   Double_t xt = R.BreitWigner(0.3, 2.5);
   xini->Fill(xt);
   Double_t x = Reconstruct2( xt, R );
   if (x != cutdummy) 
   {
      Adet->Fill(x, xt);
      bini->Fill(x);
   }
}

for (Int_t i=0; i<10000; i++) 
{
   Double_t xt = R.Gaus(0.0, 2.0);
   datatrue->Fill(xt);
   Double_t x = Reconstruct2( xt, R );
   if (x != cutdummy)
   data->Fill(x);
}

for (int i=1; i<=data->GetNbinsX(); i++) 
{
    statcov->SetBinContent(i,i,data->GetBinError(i)*data->GetBinError(i));
}

TSVDUnfold *tsvdunf = new TSVDUnfold( data, statcov, bini, xini, Adet );

tsvdunf->SetNormalize( kFALSE ); // no normalisation here

TH1D* unfres = tsvdunf->Unfold( 13 );

TH1D* ddist = tsvdunf->GetD();

TH1D* svdist = tsvdunf->GetSV();

TH2D* ustatcov = tsvdunf->GetUnfoldCovMatrix( statcov, 100 );

TH2D* uadetcov = tsvdunf->GetAdetCovMatrix( 100 );

ustatcov->Add( uadetcov );

TH2D* utaucov = tsvdunf->GetXtau();
utaucov->Add( uadetcov );

TH2D* uinvcov = tsvdunf->GetXinv();

for (int i=1; i<=unfres->GetNbinsX(); i++) 
{
   unfres->SetBinError(i, TMath::Sqrt(utaucov->GetBinContent(i,i)));
}

xini->Scale(0.7*datatrue->Integral()/xini->Integral());

TLegend *leg = new TLegend(0.58,0.60,0.99,0.88);
leg->SetBorderSize(0);
leg->SetFillColor(0);
leg->SetFillStyle(0);
leg->AddEntry(unfres,"Unfolded Data","p");
leg->AddEntry(datatrue,"True Data","l");
leg->AddEntry(data,"Reconstructed Data","l");
leg->AddEntry(xini,"True MC","l");

TCanvas *c1 = new TCanvas( "c1", "Unfolding toy example with TSVDUnfold", 1000, 900 );

c1->Divide(1,2);
TVirtualPad * c11 = c1->cd(1);

TH1D* frame = new TH1D( *unfres );
frame->SetTitle( "Unfolding toy example with TSVDUnfold" );
frame->GetXaxis()->SetTitle( "x variable" );
frame->GetYaxis()->SetTitle( "Events" );
frame->GetXaxis()->SetTitleOffset( 1.25 );
frame->GetYaxis()->SetTitleOffset( 1.29 );
frame->Draw();

data->SetLineStyle(kDashed);
data->SetLineColor(4);
data->SetLineWidth(2);
unfres->SetMarkerStyle(20);
datatrue->SetLineColor(2);
datatrue->SetLineWidth(2);
xini->SetLineStyle(kDashed);
xini->SetLineColor(8);
xini->SetLineWidth(2);

unfres->Draw("same");
datatrue->Draw("same");
data->Draw("same");
xini->Draw("same");

leg->Draw();

TVirtualPad * c12 = c1->cd(2);
c12->Divide(2,1);
TVirtualPad * c2 = c12->cd(1);
c2->SetRightMargin   ( 0.15         );

TH2D* covframe = new TH2D( *ustatcov );
covframe->SetTitle( "TSVDUnfold covariance matrix" );
covframe->GetXaxis()->SetTitle( "x variable" );
covframe->GetYaxis()->SetTitle( "x variable" );
covframe->GetXaxis()->SetTitleOffset( 1.25 );
covframe->GetYaxis()->SetTitleOffset( 1.29 );
covframe->Draw();

ustatcov->SetLineWidth( 2 );
ustatcov->Draw( "colzsame" );

TVirtualPad * c3 = c12->cd(2);
c3->SetLogy();

TLine *line = new TLine( 0.,1.,40.,1. );
line->SetLineStyle(kDashed);

TH1D* dframe = new TH1D( *ddist );
dframe->SetTitle( "TSVDUnfold |d_{i}|" );
dframe->GetXaxis()->SetTitle( "i" );
dframe->GetYaxis()->SetTitle( "|d_{i}|" );
dframe->GetXaxis()->SetTitleOffset( 1.25 );
dframe->GetYaxis()->SetTitleOffset( 1.29 );
dframe->SetMinimum( 0.001 );
dframe->Draw();

ddist->SetLineWidth( 2 );
ddist->Draw( "same" );
line->Draw();

gROOT->GetListOfCanvases()->Draw();
}

void z()
{

    int bins = 50;  // Start with a smaller number of bins

    // Initialize parameters
    double A = 197, Vo = 2.12, RR = 6.38, a0 = 0.7, phi_min = 0, phi_max = pi/9;
    double tx_min = 0, tx_max = 0.1, ty_min = 0, ty_max = 0.1;
    double x_min = 0, x_max = 15, y_min = 0, y_max = 15;
    double t_min = 0, t_max = 0.1, q_min = 0, q_max = 0.316;
    double r_min = 0, r_max = 15;
    double sigma100 = 0.1; // 100 MeV

    // Adjust integration settings
    ROOT::Math::IntegratorOneDimOptions::SetDefaultAbsTolerance(1e-3); // Increase tolerance
    ROOT::Math::IntegratorOneDimOptions::SetDefaultRelTolerance(1e-1); // Increase tolerance

    // Create instances of your classes
    FormFactor_resolution_add_wedge_1D ff_wResCut_100(A, Vo, RR, a0, t_min, t_max, bins, phi_min, phi_max, sigma100, r_min, r_max);
    FormFactor_t_1D ff(A, Vo, RR, a0, t_min, t_max, bins, r_min, r_max);

    // Get histograms
    TH1D* histTrue = ff.getFF_hist();
    TH1D* histMeasured = ff_wResCut_100.getWedgeRes_hist_1D();

    // Construct response matrix
    TH2D* Adet = new TH2D("", "Response Matrix", bins, 0, bins, bins, 0, bins);
    for (int i = 1; i <= bins; ++i) {
        double t_true = histTrue->GetXaxis()->GetBinCenter(i);
        for (int j = 1; j <= bins; ++j) {
            double t_measured = histMeasured->GetXaxis()->GetBinCenter(j);
            double response_value = ff_wResCut_100.getWedgeRes_fun_1D()->Eval(t_true);
            Adet->SetBinContent(i, j, response_value);
        }
    }

    // Check the condition number of the response matrix
    TMatrixD AdetMatrix(bins, bins);
    for (int i = 0; i < bins; ++i) {
        for (int j = 0; j < bins; ++j) {
            AdetMatrix(i, j) = Adet->GetBinContent(i + 1, j + 1);
        }
    }
    
    TDecompSVD svd(AdetMatrix);
    double conditionNumber = svd.Condition();
    std::cout << "Condition number of the response matrix: " << conditionNumber << std::endl;

    // Regularize the response matrix if needed
    if (conditionNumber > 1e10) {  // Arbitrary threshold for regularization
        for (int i = 0; i < bins; ++i) {
            AdetMatrix(i, i) += 1e-6;  // Add a small value to diagonal elements
        }
        for (int i = 0; i < bins; ++i) {
            for (int j = 0; j < bins; ++j) {
                Adet->SetBinContent(i + 1, j + 1, AdetMatrix(i, j));
            }
        }
    }

    // Create initial guesses
    TH1D* bini = (TH1D*)histMeasured->Clone("bini");
    TH1D* xini = (TH1D*)histTrue->Clone("xini");

    // Construct covariance matrix
    TH2D* Bcov = new TH2D("", "Covariance Matrix", bins, 0, bins, bins, 0, bins);
    for (int i = 1; i <= bins; ++i) {
        for (int j = 1; j <= bins; ++j) {
            double covariance = 0.0;
            if (i == j) {
                covariance = histMeasured->GetBinError(i) * histMeasured->GetBinError(i);
            }
            Bcov->SetBinContent(i, j, covariance);
        }
    }

    // Unfold the measured data
    double regularizationParameter = 3;  // Increase regularization parameter if needed
    TSVDUnfold unfold(histMeasured, Bcov, bini, xini, Adet);
    TH1D* histUnfolded = (TH1D*)unfold.Unfold(regularizationParameter);
    histUnfolded->SetTitle("Unfolded Distribution");

    // Compare histograms
    TCanvas* c1 = new TCanvas("", "", 800, 600);
    histTrue->SetLineColor(kRed);
    histTrue->SetLineStyle(1);
    histTrue->Scale(197. / histTrue->Integral(), "width");
    histMeasured->SetLineColor(kBlue);
    histMeasured->SetLineStyle(3);
    histMeasured->Scale(197. / histMeasured->Integral(), "width");
    histUnfolded->SetLineColor(kGreen);
    histUnfolded->Scale(197. / histUnfolded->Integral(), "width");

    histTrue->Draw();
    histMeasured->Draw("same");
    histUnfolded->Draw("same");

    auto legend = new TLegend(0.53, 0.65, 0.9, 0.9);
    legend->AddEntry(histUnfolded, "Unfolded FF", "l");
    legend->AddEntry(histMeasured, "Measured FF (with res and wedge)", "l");
    legend->AddEntry(histTrue, "True FF", "l");
    legend->Draw();

    gPad->SetLogy(0);
    gStyle->SetOptStat(0);
    c1->Draw();

    
}



/* differential cross-section */

/*
void params()
{
double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7, a = 0.535;
    double bins = 1000, phi_min = 0, phi_max = pi;

    //sigma = 25 MeV, 50 MeV
    double sigma25 = 0.025, sigma50 = 0.05;

    // 1D params
    double t_min = 0, t_max = 0.17;
    double r_min = 0, r_max = 15;


    FormFactor_t_1D ff(A,Vo,R,a0,t_min,t_max,bins,r_min,r_max);
    TH1D *ff_true = ff.getFF_hist();

    FormFactor_resolution_add_wedge_1D ff_wResCut_25(A,Vo,R,a0,t_min,t_max,bins,phi_min,phi_max,sigma25,r_min,r_max);

    TH1D *ff_25 = ff_wResCut_25.getWedgeRes_hist_1D();
    ff_25->Scale((1/pi)*197./ff_25->Integral(), "width");

    ff_true->GetYaxis()->SetTitle("|F(t)|^{2}");
    ff_true->GetXaxis()->SetTitle("|t| [GeV^{2}/c^{2}]");
    ff_true->GetYaxis()->SetTitleSize(0.047);
    ff_true->GetXaxis()->SetTitleSize(0.047);
    ff_true->Scale(1/(pi)*197./ff_true->Integral(), "width");
    ff_true->GetYaxis()->SetTitleOffset(0.9);
    ff_true->GetXaxis()->SetTitleOffset(0.9);
    ff_true->SetLineStyle(1);
    ff_true->SetLineWidth(2);
    ff_true->SetLineColor(kBlack);
    ff_true->Draw();
    TH1D *ff_25 = ff_wResCut_25.getWedgeRes_hist_1D();
    ff_25->Scale((1/pi)*197./ff_25->Integral(), "width");
    ff_25->GetYaxis()->SetRangeUser(1e-2,1e5);
    //ff_25->SetMarkerStyle(26);
    //ff_25->SetMarkerSize(1);
    //ff_25->SetMarkerColor(kGreen+2);
    ff_25->SetLineStyle(2);
    ff_25->SetLineWidth(3);
    ff_25->SetLineColor(kGreen+2);
    ff_25->Draw("same");

    FormFactor_resolution_add_wedge_1D ff_wResCut_50(A,Vo,R,a0,t_min,t_max,bins,phi_min,phi_max,sigma50,r_min,r_max);

    TH1D *ff_50 = ff_wResCut_50.getWedgeRes_hist_1D();
    ff_50->Scale((1/pi)*197./ff_50->Integral(), "width");
    ff_50->GetYaxis()->SetRangeUser(1e-2,1e5);
    //ff_50->SetMarkerStyle(27);
    //ff_50->SetMarkerSize(0.9);
    //ff_50->SetMarkerColor(kOrange+7);
    ff_50->SetLineStyle(3);
    ff_50->SetLineWidth(3);
    ff_50->SetLineColor(kOrange+2);
    ff_50->Draw("same");

    double sigma100 = 0.1;
    FormFactor_resolution_add_wedge_1D ff_wResCut_100(A,Vo,R,a0,t_min,t_max,bins,phi_min,phi_max,sigma100,r_min,r_max);

    TH1D *ff_100 = ff_wResCut_100.getWedgeRes_hist_1D();
    ff_100->Scale((1/pi)*197./ff_100->Integral(), "width");
    ff_100->GetYaxis()->SetRangeUser(1e-2,1e5);
    //ff_100->SetMarkerStyle(26);
    //ff_100->SetMarkerSize(1);
    //ff_100->SetMarkerColor(kMagenta);
    ff_100->SetLineStyle(9);
    //ff_100->SetMarkerSize(1);
    ff_100->SetLineColor(kMagenta);
    ff_100->Draw("same");
}
*/

double qmax2_new(double E,double omega)
{
	/* 
    Eqn (4)
	maximum momentum exchange allowed for a photon emitted by an electron
	determined by electron energy loss
    returns [GeV^2]
	*/
    if((E-mass_electric)<=omega)
    {
        return 0;
    }
    else{
        double C1 = sqrt(E*E-mass_electric*mass_electric);  // [GeV]
        double C2 = sqrt((E-omega)*(E-omega)-mass_electric*mass_electric);  // [GeV]
        double pgamma = C1+C2;  // [GeV]
        return -(omega*omega-pgamma*pgamma);  // [GeV^2]
    }
    
}

double qmin2_new(double E,double omega)
{
	/*
    Eqn (3)
	photon flux small at Q^2 min
	depends on final state
    returns [GeV^2]
	*/
    if((E-mass_electric)<=omega)
    {
        return 0;
    }
    else{
       return -(2*E*omega-2*E*E+2*mass_electric*mass_electric+2*sqrt((E*E-mass_electric*mass_electric)*((E-omega)*(E-omega)-mass_electric*mass_electric)));  // [GeV^2]
    }
    
}

double getpz(double Egamma) // this function is never used
{ 
    /* 
    function call for dsigma_dydtdQ2 but param is then never used
    */
	double a = Egamma;
    double b = 0.938;
    double c = 3.1;
	double temp1 = 8.0*(2.0*a*b+b*b);
	double temp2 = 8.0*a*a*b+8.0*a*b*b-4.0*a*c*c;
	double temp3 = 4.0*b*(8.0*a+4.0*b)*c*c*(-4.0*a*b-4.0*b*b+c*c)+a*a*pow(-8.0*a*b-8.0*b*b+4.0*c*c,2);
	double temp4 = 4.0*b*(8.0*a+4.0*b)*c*c*(-4.0*a*b-4.0*b*b+c*c);
	double temp5 = a*a*pow((-8.0*a*b-8.0*b*b+4.0*c*c),2);

	//cout<<temp4<<"  "<<temp5<<"  "<<temp3<<endl;
	return (temp2-pow(temp3,0.5))/temp1;
}

double N_omega_Q2(double omega,double Q2)
{
	/*
    Eqn (2)
	photon flux, k corresponds to omega, rest frame
    returns [dimensionless]
	*/
    double E = Gamma_electron*mass_electric;  // [GeV]
    if((E-10*mass_electric)>omega)
    {
        double qmin2 = qmin2_new(E,omega);  // [GeV^2]
        double C1 = alpha/pi;  // [dimensionless]
        double C2 = 1-omega/E+omega*omega/(2*E*E);  // [dimensionless]
        double C3 = (1-omega/E)*abs(qmin2/Q2);  // [dimensionless]
        return C1*(C2-C3);  // [dimensionless]
    }
    else{return 0;}
    
}

double dsigma_dt_origin(const double Wgp)
{
    /*
    used in dsigma_gammaA_VA_overdt_origin function below
    returns [fm^2/GeV^2]
    */
	double sigmagp_r = 0;  // [fm^2] >> assuming since in .h file there's a value of cross section in [fm^2]
    if(Wgp>(mass_Jpsi+mass_proton))
    {
       sigmagp_r = bv_Jpsi*1.E-7*80.2*pow((1-(mass_proton+mass_Jpsi)*(mass_proton+mass_Jpsi)/(Wgp*Wgp)),1.5)*pow((Wgp*Wgp/10000.),0.321); // [GeV^-2]
        return sigmagp_r;  // [fm^2/GeV^2]
    }
    else{return 0;}
}

double rws(double r)
{
    /*
    A: scaling for nuclear effects. If no nuclear effect, A=1
    Woods-Saxon potential (nuclear density distribution)
    used for nuclear thickness calculation
    returns [dimensionless]
    */
    double A = 197;
	return 1.0 / (1.+ exp((r - pow(A, (1./3.))*1.16 * (1. - 1.16 * pow(A, (-2. / 3.)))) / depth_Au));  // [dimensionless]
}

double thickness(double r)
{
    /*
    T_AA in paper, used to calculate next function sigma_total_VA
    nuclear thickness = 2 R_A
    returns [fm]
    */
	const unsigned int nmbPoints         = 5;  // [dimensionless]
	const double       xg[nmbPoints + 1] = {0., 0.1488743390, 0.4333953941, 0.6794095683, 0.8650633667, 0.9739065285};  // [dimensionless]
	const double       ag[nmbPoints + 1] = {0., 0.2955242247, 0.2692667193, 0.2190863625, 0.1494513492, 0.0666713443};  // [dimensionless]
    const double zMin   = 0;  // [fm]
	const double zMax   = 15;  // [fm]
	const double zRange = 0.5 * (zMax - zMin);  // [fm]
	const double zMean  = 0.5 * (zMax + zMin);  // [fm]
	double       sum    = 0;  // [dimensionless]
	for(unsigned int i = 1; i <= nmbPoints; ++i){
		double zsp    = zRange * xg[i] + zMean;  // [fm]
		double radius = sqrt(r * r + zsp * zsp);  // [fm]
		sum          += ag[i] * rws(radius);  // [dimensionless]
		zsp           = zRange * (-xg[i]) + zMean;  // [fm]
		radius        = sqrt(r * r + zsp * zsp);  // [fm]
		sum          += ag[i] * rws(radius);  // [dimensionless]
	}
    return 2. * zRange * sum;  // [fm]
}

double sigma_total_VA(const double cross_Vp)
{
    /*
    Eqn (9)
    quantum Glauber calculation
    returns [fm^2]
    */
    double r;  
	double arg,Pint;
    Pint = 0;
    for(int i=0;i<1000;i++)
    {
        r = i*0.1;  // [fm]
        arg = - cross_Vp * thicknessnorm * thickness(r);  // [fm]
	    Pint += 2*pi*r*0.1*(1.0 - exp(arg * GT));  // [fm]
    }
	return Pint*.75/ GT;  // [fm^2]
}

double dsigma_gammaA_VA_overdt_origin(double Egamma)
{
	/*
	used in Eqn (13) and the next function
    returns [fm^2/GeV^2]
	*/
    double Wgp;
    double Mp = mass_proton;  // [GeV]
    Wgp = sqrt(2. * Egamma * (Ep+ sqrt(Ep * Ep - Mp * Mp))+ Mp * Mp);  // rest frame  // [GeV]
    double sigma_total_Vp = sqrt(16.*pi*vmphotoncoupling_Jpsi/(alpha)*NQ*NQ*dsigma_dt_origin(Wgp)*hbarc*hbarc);  // total VM-proton cross-section  // [fm^2]
    double sigma_tot_VA = sigma_total_VA(sigma_total_Vp);  // total VM nucleon cross-section  // [fm^2]
    double sigma_gA_VA_overdt = alpha/(16*pi*vmphotoncoupling_Jpsi*NQ*NQ)*sigma_tot_VA*sigma_tot_VA;  // [fm^4]
    return sigma_gA_VA_overdt/(hbarc*hbarc);  // [fm^2/GeV^2]
}

double formFactor(const double t)
{
    /*
    Eqn (14) , F(q^2)
    Diffractive minima when t = multiple of pi/R_A
    A: scaling for nuclear effects. If no nuclear effect, A=1
    returns [dimensionless]
    */
	if(t==0){return 0;}
    else{
        double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7;
        double q = sqrt(t);
	    const double arg1 = q*R / hbarc;  
	    const double arg2 = hbarc / q;
	    const double arg3  = (sin(arg1)-arg1*cos(arg1)) * 4*pi*Vo * arg2*arg2*arg2 / A; 
	    double result = 1/(2*pi)*arg3 / (1. + (a0*a0 * q*q/(hbarc*hbarc))); 
        return result; // [-]
        //double A = 197;
	    //const double R    = pow(A, (1./3.))*1.16 * (1. - 1.16 * pow(A, -2. / 3.));  // [fm] >> assuming since they put that a0 below is in [fm]
	    //const double q    = sqrt(t);  // [GeV]
	    //const double arg1 = q * R / hbarc;  // [dimensionless]
	    //const double arg2 = hbarc / (q * 1.16 * (1. - 1.16 * pow(A, -2. / 3.)));  // [dimensionless]
	    //const double sph  = (sin(arg1) - arg1 * cos(arg1)) * 3. * arg2 * arg2 * arg2 / double(A);  // [dimensionless]
	    //const double a0   = 0.70;  // [fm]
	    //return sph / (1. + (a0 * a0 * t) / (hbarc * hbarc));  // [dimensionless]
    }
    
}

double dsigma_dydtdQ2(double *x,double *par)
{
    /*
    Eqn (10) >> and (13)??
    returns [fm^2/GeV^4]
    */
    double Q2 = x[0];  // [GeV^2]
    double y = par[0];  // [dimensionless] >> rapidity
    double t = par[1];
    double Egamma = 0.5*mass_Jpsi*exp(-y);  // Eqn (25)  // [GeV]
    double Egamma_target = Egamma*2.*Gamma_proton;  // [GeV]
    double pz = getpz(Egamma_target);  // doesn't get used
    double tmin = (mass_Jpsi*mass_Jpsi/(4*Gamma_proton*Egamma))*(mass_Jpsi*mass_Jpsi/(4*Gamma_proton*Egamma)); 
    if(t>tmin)
    {
        double t_true = t+tmin;  
        double n_Jpsi = c1_Jpsi+c2_Jpsi*(Q2+mass_Jpsi*mass_Jpsi);  // [dimensionless]
        double C = pow(mass_Jpsi*mass_Jpsi/(mass_Jpsi*mass_Jpsi+Q2),n_Jpsi);  // [dimensionless]
        double temp = N_omega_Q2(Egamma,Q2)/Q2*dsigma_gammaA_VA_overdt_origin(Egamma)*formFactor(t_true)*formFactor(t_true)*C;  // [fm^2/GeV^4]
        return temp;
    }
    else{return 0;}
}

double inte_Q2(double y, double t)
{
	/*
    Eqn (11)
	single differential photon flux
    returns [fm^2/GeV^2]
	*/
    double E = Gamma_electron*mass_electric;  // [GeV]
    double omega = 0.5*mass_Jpsi*exp(-y);  // [GeV]
    if(omega>(E-10*mass_electric)){return 0;}
    //returns a=0 until y=-1.2 for function below
    else{
        double qmax2 = 1;  // [GeV^2]
        double qmin2 = qmin2_new(E*2*Gamma_proton,omega*2*Gamma_proton);  // [GeV^2]
        TF1 *f = new TF1("",dsigma_dydtdQ2,qmin2,qmax2,2);  // [fm^2/GeV^4]
        f->SetParameters(y,t);
        return f->Integral(qmin2,qmax2,1e-3);  // [fm^2/GeV^2]
    }
}

void drawdsigmadydt()
{
    TH2D *d2sigmadydt = new TH2D("d2sigmadydt_10pts","",50,-2,3,10,0,0.1);
    for(int i=0;i<50;i++)
    {
        double y = -2+i*0.1;  // rapidity  // [dimensionless]
        for(int j=0;j<50;j++)
        {
            double t = (j+1)*0.01;  // [GeV^2]
	        double a = inte_Q2(y,t);  // [fm^2/GeV^2]
            d2sigmadydt->SetBinContent(i+1,j+1,a);
        }
        cout<<"y:"<<y<<endl;
    }
	TFile *eAu_5_41_dsigmadydt = new TFile("d2sigmadydt_10pts.root","recreate");
	//d2sigmadydt->SetTitle("Coherent");
	//d2sigmadydt->GetXaxis()->SetTitle("Rapidity");
	//d2sigmadydt->GetYaxis()->SetTitle("d^{2}#sigma/dydt");
	d2sigmadydt->Write();
	eAu_5_41_dsigmadydt->Close();
}

void drawdisgmadt()
{
	/*
	differential cross-section as a function of t
	y-projection
	*/
    TFile *f = new TFile("d2sigmadydt_10pts.root");
    TH2D *td = (TH2D *)f->Get("d2sigmadydt_10pts");
    TH1D *td1 = (TH1D *)td->ProjectionY();
    td1->Scale(197./td->Integral(), "width");
    td1->SetTitle("Coherent y-Projection");
    td1->GetXaxis()->SetTitle("t_{#perp_x} [(GeV/c)^2]");
    td1->GetYaxis()->SetTitle("d#sigma/dt [fm^2/(GeV/c)^2]");
    td1->SetFillStyle(0);
    TCanvas *c1 = new TCanvas();
	td1->Draw();
    c1->Update();
	//td12->Draw("same");
	gPad->SetLogy(1);
    //td1->Draw();
}

double dSigma_dt(double *x, double *par) 
{
    double t = x[0];
    TFile *file = new TFile("d2sigmadydt_10pts.root");
    TH2D *d2sigmadydt = (TH2D *)file->Get("d2sigmadydt_10pts");

    int tBin = d2sigmadydt->GetYaxis()->FindBin(t); 
    double integratedValue = 0.0;

    if (tBin >= 1 && tBin <= d2sigmadydt->GetNbinsY()) 
    {
        TH1D *yProjection = d2sigmadydt->ProjectionX("_py", tBin, tBin);
        integratedValue = yProjection->Integral();  // Integrate over the y-axis
        delete yProjection;  // Clean up
    }

    file->Close();
    return integratedValue;  // Return the integrated value
}

void createTF1() {
    // Scale the function to be in nanobarns (nb)
    auto scaled_dSigma_dt = [](double *x, double *p) {
        return 1e7 * dSigma_dt(x, p); // Multiply the original function by 10^7
    };

    // Create the TF1 object with the scaled function
    TF1 *f_dSigma_dt = new TF1("f_dSigma_dt", scaled_dSigma_dt, 0, 0.1);

    // Update the title to reflect the new units
    f_dSigma_dt->SetTitle("; |t| [GeV^{2}/c^{2}]; d#sigma/dt [nb/(GeV/c)^{2}]");

    // Create a canvas and draw the function
    TCanvas *c = new TCanvas();
    f_dSigma_dt->Draw();
}

    
double raw_dsigma_dt(double t)
{
    const double y_min = -2.0;  // Rapid range: [y_min, y_max]
    const double y_max = 3.0;

    // Integrate inte_Q2 over y using ROOT's TF1 integration
    TF1 *dsigma_dy = new TF1("dsigma_dy", [](double *y, double *params) {
        double y_val = y[0];
        double t = params[0];  // Transfer momentum
        return inte_Q2(y_val, t);  // Your d^2σ/dy dt function
    }, y_min, y_max, 1);

    dsigma_dy->SetParameter(0, t);  // Pass t to the function
    double result = dsigma_dy->Integral(y_min, y_max, 1e-6);  // Perform integration
    delete dsigma_dy;

    return result;  // [fm^2/GeV^2]
}

double smear_dsigma_dt(double t, double sigma)
{
    TF1 *smearing_func = new TF1("smearing_func", [](double *x, double *params) {
        double t_prime = x[0];
        double t = params[0];
        double sigma = params[1];
        return exp(-(t_prime - t) * (t_prime - t) / (2 * sigma * sigma)) * raw_dsigma_dt(t_prime);
    }, t - 5 * sigma, t + 5 * sigma, 2);

    smearing_func->SetParameters(t, sigma);
    double smeared_result = smearing_func->Integral(t - 5 * sigma, t + 5 * sigma, 1e-6);
    delete smearing_func;

    return smeared_result;  // [fm^2/GeV^2]
}

double wedge_cut_dsigma_dt(double t, double phi_min, double phi_max, double sigma)
{
    TF1 *wedge_func = new TF1("wedge_func", [](double *x, double *params) {
        double phi = x[0];
        double t = params[0];
        double sigma = params[1];
        return smear_dsigma_dt(t, sigma);
    }, phi_min, phi_max, 2);

    wedge_func->SetParameters(t, sigma);
    double wedge_result = wedge_func->Integral(phi_min, phi_max, 1e-6);
    delete wedge_func;

    return wedge_result;  // [fm^2/GeV^2]
}

void process_dsigma_dt_with_smearing_and_wedge()
{
    double phi_min = 0.0;
    double phi_max = M_PI / 12;  // Wedge range
    double sigma = 0.05;        // Smearing width
    double t_min = 0.0;
    double t_max = .1;
    int bins = 10;

    // Initialize histogram
    TH1D *hist = new TH1D("Smearing and Wedge Cut on dσ/dt", "", bins, t_min, t_max);

    // Fill histogram
    for (int i = 1; i <= hist->GetNbinsX(); i++) {
        double t = hist->GetBinCenter(i);
        cout << "t:" << t << endl;
        double value = wedge_cut_dsigma_dt(t, phi_min, phi_max, sigma);
        hist->SetBinContent(i, value);
    }

    // Draw and save the histogram
    TCanvas *canvas = new TCanvas("canvas", "Smeared and Wedge-Cut dσ/dt", 800, 600);
    hist->Draw();
    canvas->SaveAs("smearing_wedge_dsigma_dt_final.png");

    delete canvas;
    delete hist;
}

 
void dsigma_test_fun()
{
    double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7;
    double t_min = 0, t_max = 0.1, bins = 100;
    double phi_min = 0,  phi_max = pi/12, sigma = 0.05;
    double r_min = 0, r_max = 15;

    dsigma_resolution_add_wedge_1D test(A,Vo,R,a0,t_min,t_max,bins,phi_min,phi_max,sigma,r_min,r_max);
    TF1 *test_fun = test.get_fun_dsigma();
    test_fun->Draw();
    gPad->SetLogy(1);

}

void dsigma_test_hist()
{
    double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7;
    double t_min = 0, t_max = 0.2, bins = 1000;
    double phi_min = 0,  phi_max = pi/12, sigma = 0.05;
    double r_min = 0, r_max = 15;

    dsigma_resolution_add_wedge_1D test(A,Vo,R,a0,t_min,t_max,bins,phi_min,phi_max,sigma,r_min,r_max);
    
    TH1D *test_hist = test.get_hist_dsigma();
    test_hist->Draw();
    gPad->SetLogy(1);
}

void dsigma_test_fun_cut()
{
    double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7;
    double t_min = 0, t_max = 0.1, bins = 100;
    double phi_min = 0,  phi_max = pi/12, sigma = 0.1;
    double r_min = 0, r_max = 15;

    dsigma_resolution_add_wedge_1D test(A,Vo,R,a0,t_min,t_max,bins,phi_min,phi_max,sigma,r_min,r_max);
    TF1 *test_fun = test.get_fun_dsigma();
    //test_fun->Draw();

    WedgeResolution cut(test_fun,A,Vo,R,a0,t_min,t_max,bins,phi_min,phi_max,sigma,r_min,r_max);
    TF1 *test_cut = cut.getWedgeRes_fun_1D();
    test_cut->Draw();

    gPad->SetLogy(1);

}

void dsigma_test_hist_cut()
{
    double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7;
    double t_min = 0, t_max = 0.1, bins = 10;
    double phi_min = 0,  phi_max = pi/12, sigma = 0.05;
    double r_min = 0, r_max = 15;

    dsigma_resolution_add_wedge_1D test(A,Vo,R,a0,t_min,t_max,bins,phi_min,phi_max,sigma,r_min,r_max);
    
    //TH1D *test_hist = test.getWedgeRes_hist_1D();
    //test_hist->Draw();
    gPad->SetLogy(1);
}




    // Integrand to calcuate smearing
    static double guassian(double *x, double *par)
    {
        double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7;
    double t_min = 0, t_max = 0.2, bins = 1000;
    double phi_min = 0,  phi_max = pi/12, sigma = 0.05;
    double r_min = 0, r_max = 15;
        //dsigma_resolution_add_wedge_1D test(A,Vo,R,a0,t_min,t_max,bins,phi_min,phi_max,sigma,r_min,r_max);

        double qy = par[0], qx_prime = par[1];
        //double qy = sqrt(ty), qx_prime = sqrt(tx_prime);
        double qx = x[0];  // integration variable
        double q = sqrt(qx*qx+qy*qy);
        double t = q*q;
        TFile *f = new TFile("dsigma_dt_hist_1000pt.root");
        TH1D *td = (TH1D *)f->Get("dsigma_dt_hist");
        //TH1D *td1 = (TH1D *)td->ProjectionY();
        //td1->Scale(197./td->Integral(), "width");
        double dsig = td->GetXaxis()->GetBinCenter(t);
        //cout << "d_sigma: " << dsig << endl;
	    return dsig*exp(-(qx-qx_prime)*(qx-qx_prime)/(2*sigma*sigma));
    }
// Integral over qx -> new form factor with smearing
static double FF_wRes(double qy, double qx_prime, double sigma) 
{
    double qqx_min = 0, qqx_max = 1;
    TF1 *f = new TF1("f", guassian, qqx_min, qqx_max, 3);  
    f->SetParameters(qy, qx_prime, sigma);
    return f->Integral(qqx_min, qqx_max, 1e-6); // F(qx',qy) [GeV]
}


// Integral over qx -> new form factor with smearing
void FF_wRes_hist() 
{
    double sigma = 0.05, bins = 10, qx_prime_min = 0, qx_prime_max = 0.3162, qy_min = 0, qy_max = 0.3162;
    TH2D *testHist = new TH2D("testHist","",bins, qx_prime_min,qx_prime_max,bins,qy_min,qy_max);
    double step = qy_max/bins;
    for(int i=0; i<bins; i++)
    {
        double qy = (i+1)*step;
        for(int j=0;j<bins;j++)
        {
            double qx_prime = (j+1)*step;
            double result = FF_wRes(qy,qx_prime,sigma);
            testHist->SetBinContent(i+1,j+1,result);
        }
    }
    TFile *dsigmadt = new TFile("dsigma_smeared.root","recreate");
    if (!dsigmadt || dsigmadt->IsZombie()) {
        cerr << "Error: File could not be created." << endl;
        return;
    }
	testHist->Write();
	dsigmadt->Close();

    testHist->ProjectionY()->Draw();
    gPad->SetLogy(1);
}

static double FF_cut_wRes(double *x, double *par)
    {
        double t = par[0], sigma = par[1];
        double phi = x[0];  // integration variable
        double q = sqrt(t), qx_prime = q*sin(phi), qy = q*cos(phi); 
        // call form factor with resolution
	    return FF_wRes(qy, qx_prime, sigma); // [GeV^2]
    }

/*static double FF_cut_wRes(double* x, double* par) {
    static std::map<std::tuple<double, double, double>, double> cache; // Caching values
    double t = par[0], sigma = par[1];
    double phi = x[0];  // Integration variable
    double q = sqrt(t), qx_prime = q * sin(phi), qy = q * cos(phi);

    // Use caching for FF_wRes
    auto key = std::make_tuple(qy, qx_prime, sigma);
    if (cache.find(key) == cache.end()) {
        cache[key] = FF_wRes(qy, qx_prime, sigma);
    }

    return cache[key]; // [GeV^2]
}*/
// Integral to integrate over theta
static double wedge_FF(double t, double sigma, double phi_min, double phi_max) {
    static TF1* f = new TF1("FF_cut_wRes", [](double* x, double* par) {
        return FF_cut_wRes(x, par);
    }, phi_min, phi_max, 3);

    f->SetParameters(t, sigma);
    return f->Integral(phi_min, phi_max, 1e-6); // Use coarser precision
}


void drawthis()
{
    double sigma = 0.05, phi_min = 0, phi_max = pi/12;
    int bins = 1000;
    double t_min = 0, t_max = 0.2;

    TH1D* d2sigmadydt = new TH1D("d2sigmadydt", "", bins, t_min, t_max);

    for (int i = 0; i < bins; i++) {
        double t = d2sigmadydt->GetBinCenter(i + 1); 
        double a = wedge_FF(t, sigma, phi_min, phi_max);  // [fm^2/GeV^2]
        d2sigmadydt->SetBinContent(i + 1, a);
        std::cout << "t_hist: " << t << ", a: " << a << std::endl; 
    }
    TFile *dsigmadt = new TFile("dsigma_method_test3.root","recreate");
    if (!dsigmadt || dsigmadt->IsZombie()) {
        cerr << "Error: File could not be created." << endl;
        return;
    }
	d2sigmadydt->Write();
	dsigmadt->Close();

    d2sigmadydt->Draw();
    gPad->SetLogy(1);
}
/*
void drawthis2()
{
    double sigma = 0.05, phi_min = 0, phi_max = pi / 4;
    int bins = 10;
    double t_min = 0, t_max = 0.1;

    // Initialize |F(t)|^2 wedge cut TF1
    TF1 *cutFFt2 = new TF1("Cut Form Factor: |F(t)|^{2}", [this, phi_min, phi_max, sigma] (double *var, double *par)
    {
        double t = var[0];  // t is variable
        TF1 fft2("", FF_cut_wRes, phi_min, phi_max, 2);
        fft2.SetParameters(t, sigma);
        double integral = fft2.Integral(phi_min, phi_max, 1e-12);
        return integral;
    }, t_min, t_max, 0);

    TH1D* d2sigmadydt = new TH1D("d2sigmadydt", "", bins, t_min, t_max);
    double step3 = t_max/bins;
    for(int i=0;i<bins;i++)
        {
            double t = hist1D->GetXaxis()->GetBinCenter(i+1);
            hist1D->SetBinContent(i+1, cutFFt2->Eval(t)); 
        }
    TFile *dsigmadt = new TFile("dsigma_method_test2.root","recreate");
    if (!dsigmadt || dsigmadt->IsZombie()) {
        cerr << "Error: File could not be created." << endl;
        return;
    }
	d2sigmadydt->Write();
	dsigmadt->Close();

    d2sigmadydt->Draw();
    gPad->SetLogy(1);
}

*/

