#include "FormFactor_t_1D.h"
#include "FormFactor_resolution_add_wedge_1D.h"
#include "FormFactor_q_2D.h"
#include "FormFactor_q_2DwResCut.h"
#include "FormFactor_Wedge.h"
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
    ff_true_hist->GetYaxis()->SetTitle("|F_{#hat{n}}(t)|^{2}, |F(t)|^{2}");
    ff_true_hist->GetXaxis()->SetTitle("|t| [GeV^{2}/c^{2}]");
    ff_true_hist->GetYaxis()->SetRangeUser(1e-9,1e3);
    ff_true_hist->SetLineStyle(1);
    ff_true_hist->SetLineWidth(2);
    ff_true_hist->SetLineColor(kBlack);
    //ff_true_hist->Scale(197./ff_true_hist->Integral(), "width");
    ff_true_hist->Draw();
      

    TH1D *ff_hist_pi2 = ff_wResCut_1d_whole.getWedgeRes_hist_1D();
    ff_hist_pi2->GetYaxis()->SetRangeUser(1e-9,1e3);
    double pi2_integral = ff_hist_pi2->Integral();
    //ff_hist_pi2->Scale(197./ff_hist_pi2->Integral(), "width");
    ff_hist_pi2->Scale(ff_true_hist->Integral()/pi2_integral);
    //ff_hist_pi2->SetMarkerStyle(27);
    //ff_hist_pi2->SetMarkerSize(0.9);
    //ff_hist_pi2->SetMarkerColor(kBlue);->SetLineStyle(3);
    ff_hist_pi2->SetLineColor(kRed);
    ff_hist_pi2->SetLineWidth(3);
    ff_hist_pi2->SetLineStyle(3);
    ff_hist_pi2->Draw("same");

     TH1D *ff_hist_pi12 = ff_wResCut_1d_cut.getWedgeRes_hist_1D();
    ff_hist_pi12->GetYaxis()->SetRangeUser(1e-9,1e3);
    ff_hist_pi12->Scale(ff_true_hist->Integral()/pi2_integral);
    //ff_hist_pi12->Scale(197./ff_hist_pi12->Integral(), "width");
    //ff_hist_pi12->SetMarkerStyle(3);
    //ff_hist_pi12->SetMarkerSize(0.9);
    //ff_hist_pi12->SetMarkerColor(kMagenta);
    ff_hist_pi12->SetLineStyle(2);
    ff_hist_pi12->SetLineColor(kBlue);
    ff_hist_pi12->SetLineWidth(2);
    ff_hist_pi12->Draw("same");

    auto legend = new TLegend(0.45,0.65,0.86,0.86);
    legend->SetHeader("Coherent production","C");
    legend->AddEntry(ff_true_hist,"Truth","l");
    legend->AddEntry(ff_hist_pi2,"25 MeV/c res. w. #theta_{max} = #pi/2","l");
    legend->AddEntry(ff_hist_pi12,"25 MeV/c res. w. #theta_{max} = #pi/12","l");
    legend->SetTextSize(0.04);
    legend->SetBorderSize(0);
    legend->Draw();

    gPad->SetLogy(1);
    gStyle->SetOptStat(0);
    c->Draw();
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

    TCanvas *c = new TCanvas("c", "Canvas", 800, 600);

    FormFactor_t_1D ff(A,Vo,R,a0,t_min,t_max,bins,r_min,r_max);
    TH1D *ff_true = ff.getFF_hist();
    ff_true->GetYaxis()->SetTitle("|F(t)|^{2}");
    ff_true->GetXaxis()->SetTitle("|t| [GeV^{2}/c^{2}]");
    ff_true->GetYaxis()->SetTitleSize(0.047);
    ff_true->GetXaxis()->SetTitleSize(0.047);
    ff_true->GetYaxis()->SetRangeUser(1e-9,1e5);
    //ff_true->Scale(1/(pi)*197./ff_true->Integral(), "width");
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
    incoherent_hist->SetLineColor(kGreen+2);
    incoherent_hist->GetYaxis()->SetRangeUser(1e-9,1e5);
    //for (int i=0;i<incoherent_his1/(1/(pi)*197./incoherent_hist->Integral())t->GetN();i++) incoherent_hist->GetY()[i] *= ;
    double scaleFactor = 1/(1/(pi)*197./incoherent_hist->Integral());
    for(int i = 0; i < incoherent_hist->GetN(); ++i) 
    {
        double x, y;
        incoherent_hist->GetPoint(i, x, y);
        incoherent_hist->SetPoint(i, x, y * scaleFactor);
    }

    //incoherent_hist->Scale(1/(1/(pi)*197./incoherent_hist->Integral()), "width");
    //incoherent_hist->Scale(ff_true->Integral()/incoherent_hist->Integral(),"width");
    incoherent_hist->SetLineStyle(2);
    incoherent_hist->SetLineWidth(3);
    incoherent_hist->Draw("same");

    FormFactor_resolution_add_wedge_1D ff_wResCut_25(A,Vo,R,a0,t_min,t_max,bins,phi_min,phi_max,sigma25,r_min,r_max);
    TH1D *ff_25 = ff_wResCut_25.getWedgeRes_hist_1D();
    //ff_25->Scale(1/(pi)*197./ff_25->Integral(), "width");
    ff_25->Scale(ff_true->Integral()/ff_25->Integral());
    ff_25->GetYaxis()->SetRangeUser(1e-9,1e5);
    ff_25->SetLineStyle(3);
    ff_25->SetLineWidth(3);
    ff_25->SetLineColor(kRed);
    ff_25->Draw("same");

    // Upper right legend
    auto legend = new TLegend(0.47,0.72,0.77,0.88);
    legend->AddEntry(ff_true,"Truth: coherent","l");
    legend->AddEntry(incoherent_hist,"Truth: incoherent","l");
    legend->AddEntry(ff_25,"Coherent with 25 MeV/c res.","l");
    legend->SetTextSize(0.04);
    legend->SetBorderSize(0);
    legend->Draw("same");
   
    gPad->SetLogy(1);
    gStyle->SetOptStat(0);
    c->Draw();
}

void FFq_2D_hist()
{
    double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7;
    double qy_min = 0, qy_max = 0.5, qx_min = 0, qx_max = 0.5, bins = 1000;
    FormFactor_q_2D ff(A,Vo,R,a0,qy_min,qy_max,qx_min,qx_max,bins);
    TH2D *ff_histz = ff.getFormFactorq_2D_hist();
    TCanvas *c1 = new TCanvas("c1", "c1",800, 800);
    c1->SetLeftMargin(0.14); 
    c1->SetRightMargin(0.14);
    c1->SetBottomMargin(0.13);

    ff_histz->SetTitle("|F(t)|^{2}");
    ff_histz->GetYaxis()->SetTitleSize(0.04);
    ff_histz->GetXaxis()->SetTitleSize(0.04);
    ff_histz->GetXaxis()->SetTitle("q_{x} = #sqrt{t_{x}} [GeV/c]");
    ff_histz->GetYaxis()->SetTitle("q_{y} = #sqrt{t_{y}} [GeV/c]");
    ff_histz->GetYaxis()->SetTitleOffset(1.6);
    ff_histz->GetXaxis()->SetTitleOffset(1.1);
    ff_histz->Draw("colz");

    TGaxis *xAxis = new TGaxis(ff_histz->GetXaxis()->GetXmin(), 0, ff_histz->GetXaxis()->GetXmax(), 0, ff_histz->GetXaxis()->GetXmin(), ff_histz->GetXaxis()->GetXmax(), 510, "U");
    xAxis->SetLabelSize(0.02);
    xAxis->SetLabelOffset(0.02);
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

    /*   double theta2 = -pi/12; 
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

    auto el3 = new TEllipse(0, 0, -0.3, -0.3, -15, 15);
    el3->SetFillStyle(0);
    el3->SetLineColor(kBlack);
    el3->SetLineWidth(2);
    el3->SetTheta(270);
    el3->Draw();
    */
    auto el = new TEllipse(0, 0, -0.3, 0.3, 0, 15);
    el->SetFillStyle(0);
    el->SetLineColor(kBlack);
    el->SetLineWidth(2);
    el->SetTheta(270);
    el->Draw();

    gPad->SetLogz(1);
    gStyle->SetOptStat(0);
    c1->Draw();
}

void hist_2d_wResCut_q()
{
    double A = 197, Vo = 2.12, R = 6.38, a0 = 0.7;
    double qy_min = 0, qy_max = .5, bins = 1000;
    double qx_prime_min = 0, qx_prime_max = .5, phi_min = 0, sigma = 0.025;

    FormFactor_q_2DwResCut ff_wResCut25(A,Vo,R,a0,qy_min,qy_max,qx_prime_min,qx_prime_max,bins,sigma);
    TCanvas *c1 = new TCanvas("", "", 800, 800);
    c1->SetLeftMargin(0.14); 
    c1->SetRightMargin(0.14);
    c1->SetBottomMargin(0.13);

    TH2D *ff_wRes25 = ff_wResCut25.getSmeared_hist();
    ff_wRes25->SetTitle("|F(t)|^{2} with 25 MeV/c Resolution");
    ff_wRes25->GetYaxis()->SetTitleSize(0.04);
    ff_wRes25->GetXaxis()->SetTitleSize(0.04);
    ff_wRes25->GetXaxis()->SetTitle("q_{x} = #sqrt{t_{x}} [GeV/c]");
    ff_wRes25->GetYaxis()->SetTitle("q_{y} = #sqrt{t_{y}} [GeV/c]");
    ff_wRes25->GetYaxis()->SetTitleOffset(1.6);
    ff_wRes25->GetXaxis()->SetTitleOffset(1.1);
    ff_wRes25->Draw("colz");

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

    /*    double theta2 = -pi/12; 
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

    auto el3 = new TEllipse(0, 0, -0.3, -0.3, -15, 15);
    el3->SetFillStyle(0);
    el3->SetLineColor(kBlack);
    el3->SetLineWidth(2);
    el3->SetTheta(270);
    el3->Draw();
    */
    auto el = new TEllipse(0, 0, -0.3, 0.3, 0, 15);
    el->SetFillStyle(0);
    el->SetLineColor(kBlack);
    el->SetLineWidth(2);
    el->SetTheta(270);
    el->Draw();


    gStyle->SetOptStat(0);
    gPad->SetLogz(1);
    c1->Draw();
}

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
    double phi_min = 0, phi_max = pi/2, bins = 1000;
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
    normalizedTrueFF->SetLineStyle(1);
    //normalizedTrueFF->Draw();
    trueFF->SetLineColor(kBlack);
    trueFF->SetLineStyle(1);
    trueFF->Draw();

    TF1 *wedgeFF = ff.getWedgeCut1D_fun();
    double wedgeFFIntegral = wedgeFF->Integral(wedgeFF->GetXmin(), wedgeFF->GetXmax());
    TF1 *normalizedWedgeFF = new TF1("", [wedgeFF, wedgeFFIntegral](double *x, double *par) 
    {
        return wedgeFF->Eval(x[0]) / wedgeFFIntegral;
    }, wedgeFF->GetXmin(), wedgeFF->GetXmax(), 0);
     normalizedWedgeFF->SetLineColor(kRed);
    normalizedWedgeFF->SetLineStyle(2);
    //normalizedWedgeFF->Draw("same");
    wedgeFF->SetLineColor(kRed);
    wedgeFF->SetLineStyle(2);
    wedgeFF->Draw("same");

    TLegend *legend = new TLegend(0.75, 0.8, 0.9, 0.9);
    legend->AddEntry(normalizedTrueFF, "True Form Factor", "l");
    legend->AddEntry(normalizedWedgeFF, "Wedge Form Factor", "l");
    legend->Draw();

    gPad->SetLogy(1);
}
