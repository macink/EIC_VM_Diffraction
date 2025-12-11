#include "RiceStyle.h"
#include "ePIC_style.C"
#include <TFile.h>
#include <TH1D.h>
#include <TRandom3.h>
#include <TMath.h>
#include <vector>

using namespace std;

void plot_t_resolution()
{
	TFile* file = TFile::Open("phi_sartre_25_10_2_eAu_10x100v5.root","READ");	
	TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";
	TH2D* h_t_res = (TH2D*) file->Get("h_t_res_proj_percent_pi12");

	TCanvas* c2 = new TCanvas("c2","c2",1,1,1000,800);
	gPad->SetTicks();
	gPad->SetLogy(1);
	gPad->SetLeftMargin(0.13);
	gPad->SetRightMargin(0.1);
	gPad->SetTopMargin(0.1);
	gPad->SetBottomMargin(0.13);

	TH1D* base1 = makeHist("base1", "", "|t| [GeV/c]^{2}", " #delta t/|t| (resolution) ", 100,0,0.2,kBlack);
	base1->GetYaxis()->SetRangeUser(1e-3, 1000);
	base1->GetXaxis()->SetTitleColor(kBlack);
	TGaxis::SetMaxDigits(3);
	fixedFontHist1D(base1,1.2,1.6);
	base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize()*1.5);
	base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize()*1.5);
	base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize()*1.8);
	base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize()*1.7);
	base1->GetXaxis()->SetNdivisions(5,5,0);
	base1->GetYaxis()->SetNdivisions(5,5,0);
	base1->Draw();
	TH2D* h_res = (TH2D*) h_t_res;

	TH1D* h_res_1D = new TH1D("h_res_1D","",100,0,0.2);
	for(int ibin=0;ibin<h_res->GetNbinsX();ibin++)
	{
		TH1D* tmp=h_res->ProjectionY("tmp",ibin+1,ibin+1);
		double sigma = tmp->GetStdDev();
		double sigmaerror = tmp->GetStdDevError();
		h_res_1D->SetBinContent(ibin+1, sigma);
		h_res_1D->SetBinError(ibin+1, sigmaerror);
	}

	h_res_1D->SetMarkerSize(1.6);
	h_res_1D->SetMarkerColor(kBlack);
	h_res_1D->SetLineColor(kBlack);
	h_res_1D->SetMarkerStyle(20);

	h_res_1D->Fit("pol0","RMS0","",0.011,0.019);//first  dip
	h_res_1D->Fit("pol0","RMS0","",0.045,0.053);//second dip
	h_res_1D->Fit("pol0","RMS0","",0.095,0.102);//third  dip
    h_res_1D->Fit("pol0","RMS0","",0.167,0.175);//fourth  dip
	h_res_1D->Draw("Psame");

	TLatex* r42 = new TLatex(0.15, 0.91, "eAu 10x100 GeV");
	r42->SetNDC();
	r42->SetTextSize(25);
	r42->SetTextFont(43);
	r42->SetTextColor(kBlack);
	r42->Draw("same");

	TLatex* r43 = new TLatex(0.8,0.91, "ePIC");
	r43->SetNDC();
	r43->SetTextSize(0.04);
	r43->Draw("same");

	TLatex* r44_2 = new TLatex(0.17, 0.18, vm_label+" #rightarrow "+daug_label );
	r44_2->SetNDC();
	r44_2->SetTextSize(30);
	r44_2->SetTextFont(43);
	r44_2->SetTextColor(kBlack);
	r44_2->Draw("same");

	TPad* drawPad = new TPad("pad_etalab_11","pad_etalab_11",0.16,0.53,0.47,0.83);
	drawPad->SetLeftMargin(0.08);
	drawPad->SetRightMargin(0.08);
	drawPad->SetTopMargin(0.0);
	drawPad->SetBottomMargin(0.08);
	drawPad->Draw("same");
	drawPad->SetTicks();
	drawPad->SetLogz(1);
	drawPad->cd();
	TH1D* base2 = makeHist("base2", "", "MC", " resolution ", 100,0,0.2,kBlack);
	base2->GetYaxis()->SetRangeUser(-10, 1.5);
	base2->GetXaxis()->SetTitleColor(kBlack);
	base2->GetXaxis()->SetLabelColor(kBlack);
	base2->GetYaxis()->SetLabelColor(kBlack);
	base2->GetXaxis()->SetTitle("Resolution");
	base2->GetYaxis()->SetTitle("MC");
	TGaxis::SetMaxDigits(3);
	fixedFontHist1D(base2,3,3);
	base2->GetYaxis()->SetTitleSize(base2->GetYaxis()->GetTitleSize()*1);
	base2->GetXaxis()->SetTitleSize(base2->GetXaxis()->GetTitleSize()*1);
	base2->GetYaxis()->SetLabelSize(base2->GetYaxis()->GetLabelSize()*1);
	base2->GetXaxis()->SetLabelSize(base2->GetXaxis()->GetLabelSize()*1);
	base2->GetXaxis()->SetNdivisions(5,5,0);
	base2->GetYaxis()->SetNdivisions(5,5,0);
	base2->Draw();

	h_res->Draw("colzsame");
	gPad->Update(); 
    TPaletteAxis *palette = (TPaletteAxis*)h_res->GetListOfFunctions()->FindObject("palette");
    palette->SetX1NDC(0.9); 
    palette->SetX2NDC(0.95); 
    palette->SetY1NDC(0.08); 
    palette->SetY2NDC(0.99); 
    gPad->Modified();
    gPad->Update();

	c2->Print("./figures/plot_t_resolution.pdf");

	//TFile* fout = new TFile("plot_t_resolution.root","RECREATE");
    //c2->Write();
    //h_res->Write();
    //h_res_1D->Write();
    //fout->Close();
}

void plot_tdist_preTDR()
{
    TFile* file = TFile::Open("phi_sartre_25_10_2_eAu_10x100v5.root","READ");
	TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

	//t distribution
	TH1D* h_t_MC = dynamic_cast<TH1D*>(file->Get("h_t_MC"));
	TH1D* h_t_REC_L = dynamic_cast<TH1D*>(file->Get("h_t_REC_EEMC_cut"));
	TH1D* h_t_REC_proj12 = dynamic_cast<TH1D*>(file->Get("h_t_REC_wRES_cut_pi12"));
    TH1D* h_t_REC_proj2 = dynamic_cast<TH1D*>(file->Get("h_t_REC_wRES_cut_pi2"));
    TH1D* h_t_REC_proj3 = dynamic_cast<TH1D*>(file->Get("h_t_REC_wRES_cut_pi3"));
    TH1D* h_t_REC_proj4 = dynamic_cast<TH1D*>(file->Get("h_t_REC_wRES_cut_pi4"));
    TH1D* h_t_REC_proj6 = dynamic_cast<TH1D*>(file->Get("h_t_REC_wRES_cut_pi6"));
    TH1D* h_t_REC_proj9 = dynamic_cast<TH1D*>(file->Get("h_t_REC_wRES_cut_pi9"));
    TH1D* h_t_REC_proj16 = dynamic_cast<TH1D*>(file->Get("h_t_REC_wRES_cut_pi16"));
    TH1D* h_t_REC_proj20 = dynamic_cast<TH1D*>(file->Get("h_t_REC_wRES_cut_pi20"));
    TH1D* h_t_REC_proj24 = dynamic_cast<TH1D*>(file->Get("h_t_REC_wRES_cut_pi24"));
    TH1F* h_phi_sartre_events = dynamic_cast<TH1F*>(file->Get("h_Nevents"));

	TCanvas* c1 = new TCanvas("c1","c1",1,1,1000,800);
	gPad->SetLogy(1);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	gPad->SetRightMargin(0.1);
	TH1D* base1 = makeHist("base1", "", "|t| [GeV/c]^{2}", "d#sigma/d|t| [nb/(GeV/c)^{2}] ", 100,0,0.18,kBlack);
	base1->GetYaxis()->SetRangeUser(1e-2, 1e6);
	base1->GetXaxis()->SetTitleColor(kBlack);
	fixedFontHist1D(base1,1.,1.2);
	base1->GetYaxis()->SetTitleOffset(1.5);
	base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize()*1.5);
	base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize()*1.5);
	base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize()*1.5);
	base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize()*1.5);
	base1->GetXaxis()->SetNdivisions(4,4,0);
	base1->GetYaxis()->SetNdivisions(5,5,0);
	base1->Draw();

	double dN_dt = 0.002;
    double sigma_coherent = 459.05; // nb
    double sigma_incoherent = 261.95; // nb
    double sigma_DIS = 50184.8; // nb
    double nEvents_coherent = h_phi_sartre_events->GetEntries(); // 6.36679 M
    double preTDR_lumi = 50761.4; // nb^-1
    int rebin_width = 2;   

    h_t_MC->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    double simu_coherent_lumi = nEvents_coherent/sigma_coherent; // nb^-1 
    double ratio_coherent_phi_preTDR = preTDR_lumi/simu_coherent_lumi;
    for (int i = 1; i <= h_t_MC->GetNbinsX(); ++i) 
    {
        double binError = h_t_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_MC->SetBinError(i, newError);
    }
	h_t_MC->SetLineStyle(1);   
	h_t_MC->SetLineWidth(1);   
	h_t_MC->SetLineColor(kBlack);
	h_t_MC->Draw("HISTsame");
	h_t_MC->Draw("same");

    h_t_REC_L->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    for (int i = 1; i <= h_t_REC_L->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_L->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_L->SetBinError(i, newError);
    }
	h_t_REC_L->SetMarkerStyle(20); // method L RECO
	h_t_REC_L->SetMarkerColor(kP8Blue);
	h_t_REC_L->Draw("PEsame");

	h_t_REC_proj2->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_proj2->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/2)));
    for (int i = 1; i <= h_t_REC_proj2->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_proj2->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_proj2->SetBinError(i, newError);
    }
	h_t_REC_proj2->SetMarkerStyle(42);
	h_t_REC_proj2->SetMarkerColor(kP10Violet);
	h_t_REC_proj2->SetLineColor(kP10Violet);
	//h_t_REC_proj2->Draw("PEsame");

    h_t_REC_proj3->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_proj3->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/3)));
    for (int i = 1; i <= h_t_REC_proj3->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_proj3->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_proj3->SetBinError(i, newError);
    }
	h_t_REC_proj3->SetMarkerStyle(4);
	h_t_REC_proj3->SetMarkerColor(kP8Cyan);
	h_t_REC_proj3->SetLineColor(kP8Cyan);
	//h_t_REC_proj3->Draw("PEsame");

    h_t_REC_proj4->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_proj4->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/4)));
    for (int i = 1; i <= h_t_REC_proj4->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_proj4->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_proj4->SetBinError(i, newError);
    }
	h_t_REC_proj4->SetMarkerStyle(24);
	h_t_REC_proj4->SetMarkerColor(kP8Orange);
	h_t_REC_proj4->SetLineColor(kP8Orange);
	//h_t_REC_proj4->Draw("PEsame");

    h_t_REC_proj6->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_proj6->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/6)));
    for (int i = 1; i <= h_t_REC_proj6->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_proj6->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_proj6->SetBinError(i, newError);
    }
	h_t_REC_proj6->SetMarkerStyle(31);
	h_t_REC_proj6->SetMarkerColor(kP8Gray);
	h_t_REC_proj6->SetLineColor(kP8Gray);
	//h_t_REC_proj6->Draw("PEsame");

    h_t_REC_proj9->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_proj9->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/9)));
    for (int i = 1; i <= h_t_REC_proj9->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_proj9->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_proj9->SetBinError(i, newError);
    }
	h_t_REC_proj9->SetMarkerStyle(33);
	h_t_REC_proj9->SetMarkerColor(kP8Green);
	h_t_REC_proj9->SetLineColor(kP8Green);
	//h_t_REC_proj9->Draw("PEsame");

    h_t_REC_proj12->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_proj12->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_t_REC_proj12->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_proj12->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_proj12->SetBinError(i, newError);
    }
	h_t_REC_proj12->SetMarkerStyle(30);
	h_t_REC_proj12->SetMarkerColor(kP8Pink);
	h_t_REC_proj12->SetLineColor(kP8Pink);
	h_t_REC_proj12->Draw("PEsame");

    h_t_REC_proj16->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_proj16->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/16)));
    for (int i = 1; i <= h_t_REC_proj16->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_proj16->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_proj16->SetBinError(i, newError);
    }
	h_t_REC_proj16->SetMarkerStyle(28);
	h_t_REC_proj16->SetMarkerColor(kP8Azure);
	h_t_REC_proj16->SetLineColor(kP8Azure);
	//h_t_REC_proj16->Draw("PEsame");

    h_t_REC_proj20->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_proj20->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/20)));
    for (int i = 1; i <= h_t_REC_proj20->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_proj20->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_proj20->SetBinError(i, newError);
    }
	h_t_REC_proj20->SetMarkerStyle(13);
	h_t_REC_proj20->SetMarkerColor(kP8Red);
	h_t_REC_proj20->SetLineColor(kP8Red);
	//h_t_REC_proj20->Draw("PEsame");

    h_t_REC_proj24->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_proj24->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/24)));
    for (int i = 1; i <= h_t_REC_proj24->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_proj24->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_proj24->SetBinError(i, newError);
    }
	h_t_REC_proj24->SetMarkerStyle(38);
	h_t_REC_proj24->SetMarkerColor(kP10Yellow);
	h_t_REC_proj24->SetLineColor(kP10Yellow);
	//h_t_REC_proj24->Draw("PEsame");
	
	TLatex* r44 = new TLatex(0.2, 0.84, "1<Q^{2}<10 GeV^{2}, 0.01 < y < 0.85");
	r44->SetNDC();
	r44->SetTextSize(20);
	r44->SetTextFont(43);
	r44->SetTextColor(kBlack);
	r44->Draw("same");
	
	TLatex* r44_0 = new TLatex(0.2, 0.79, "  |y_{"+vm_label+"}|<3.5, |M_{inv} #minus M_{"+vm_label+"}| < 0.02 GeV");
	r44_0->SetNDC();
	r44_0->SetTextSize(20);
	r44_0->SetTextFont(43);
	r44_0->SetTextColor(kBlack);
	r44_0->Draw("same");
	
    // Add labels
    TLatex* ep = new TLatex(0.18, 0.33, "ePIC");
    ep->SetNDC();
    ep->SetTextFont(62);
    ep->SetTextSize(0.030);
    ep->SetTextColor(kBlack);
    ep->Draw("same");

    TLatex* r421 = new TLatex(0.23, 0.33, " Simulation 25.10.2");
	r421->SetNDC();
	r421->SetTextSize(25);
	r421->SetTextFont(43);
	r421->SetTextColor(kBlack);
	r421->Draw("same");

	TLatex* label_5 = new TLatex(0.18, 0.23, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->SetTextSize(25);
	label_5->SetTextFont(43);
	label_5->SetTextColor(kBlack);
	label_5->Draw("same");

    TLatex* r4211 = new TLatex(0.18, 0.28, "eAu #rightarrow e'Au'#phi, 10x100 GeV");
	r4211->SetNDC();
	r4211->SetTextSize(25);
	r4211->SetTextFont(43);
	r4211->SetTextColor(kBlack);
	r4211->Draw("same");

    TLatex* r142111 = new TLatex(0.18, 0.18, "L_{int} = 10 fb^{-1}/A"); // preTDR lumi
	r142111->SetNDC();
	r142111->SetTextSize(25);
	r142111->SetTextFont(43);
	r142111->SetTextColor(kBlack);
	r142111->Draw("same");
	
	TLegend *w7 = new TLegend(0.58,0.7,0.73,0.85);
	w7->SetLineColor(kWhite);
	w7->SetFillColor(0);
	w7->SetTextSize(20);
	w7->SetTextFont(45);
	w7->AddEntry(h_t_MC, "Sartre "+vm_label+" MC ", "L");
	w7->AddEntry(h_t_REC_L, "Sartre "+vm_label+" Method L RECO", "P");
	//w7->AddEntry(h_t_REC_proj2, "Sartre "+vm_label+" RECO #theta_{max} = #pi/2", "P");
    //w7->AddEntry(h_t_REC_proj3, "Sartre "+vm_label+" RECO #theta_{max} = #pi/3", "P");
    //w7->AddEntry(h_t_REC_proj4, "Sartre "+vm_label+" RECO #theta_{max} = #pi/4", "P");
    //w7->AddEntry(h_t_REC_proj6, "Sartre "+vm_label+" RECO #theta_{max} = #pi/6", "P");
    //w7->AddEntry(h_t_REC_proj9, "Sartre "+vm_label+" RECO #theta_{max} = #pi/9", "P");
    w7->AddEntry(h_t_REC_proj12, "Sartre "+vm_label+" RECO #theta_{max} = #pi/12", "P");
    //w7->AddEntry(h_t_REC_proj16, "Sartre "+vm_label+" RECO #theta_{max} = #pi/16", "P");
    //w7->AddEntry(h_t_REC_proj20, "Sartre "+vm_label+" RECO #theta_{max} = #pi/20", "P");
    //w7->AddEntry(h_t_REC_proj24, "Sartre "+vm_label+" RECO #theta_{max} = #pi/24", "P");
	w7->Draw("same");

	c1->Print("./figures/plot_t_dist_preTDR_proj.pdf");

	//TFile* fout_preTDR = new TFile("plot_t_dist_preTDR.root","RECREATE");
    //c1->Write();
    //h_t_MC->Write();
    //h_t_REC_L->Write();
    //h_t_REC_proj12->Write();
    //fout_preTDR->Close();
}

void plot_tdist_ES()
{
    TFile* file = TFile::Open("phi_sartre_25_10_2_eAu_10x100v5.root","READ");
	TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

	//t distribution
	TH1D* h_t_MC = dynamic_cast<TH1D*>(file->Get("h_t_MC"));
	TH1D* h_t_REC_L = dynamic_cast<TH1D*>(file->Get("h_t_REC_EEMC_cut"));
	TH1D* h_t_REC_proj12 = dynamic_cast<TH1D*>(file->Get("h_t_REC_wRES_cut_pi12"));
    TH1D* h_t_REC_proj2 = dynamic_cast<TH1D*>(file->Get("h_t_REC_wRES_cut_pi2"));
    TH1D* h_t_REC_proj3 = dynamic_cast<TH1D*>(file->Get("h_t_REC_wRES_cut_pi3"));
    TH1D* h_t_REC_proj4 = dynamic_cast<TH1D*>(file->Get("h_t_REC_wRES_cut_pi4"));
    TH1D* h_t_REC_proj6 = dynamic_cast<TH1D*>(file->Get("h_t_REC_wRES_cut_pi6"));
    TH1D* h_t_REC_proj9 = dynamic_cast<TH1D*>(file->Get("h_t_REC_wRES_cut_pi9"));
    TH1D* h_t_REC_proj16 = dynamic_cast<TH1D*>(file->Get("h_t_REC_wRES_cut_pi16"));
    TH1D* h_t_REC_proj20 = dynamic_cast<TH1D*>(file->Get("h_t_REC_wRES_cut_pi20"));
    TH1D* h_t_REC_proj24 = dynamic_cast<TH1D*>(file->Get("h_t_REC_wRES_cut_pi24"));
    TH1F* h_phi_sartre_events = dynamic_cast<TH1F*>(file->Get("h_Nevents"));

	TCanvas* c1 = new TCanvas("c1","c1",1,1,1000,800);
	gPad->SetLogy(1);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	gPad->SetRightMargin(0.1);
	TH1D* base1 = makeHist("base1", "", "|t| [GeV/c]^{2}", "d#sigma/d|t| [nb/(GeV/c)^{2}] ", 100,0,0.18,kBlack);
	base1->GetYaxis()->SetRangeUser(1e-2, 1e6);
	base1->GetXaxis()->SetTitleColor(kBlack);
	fixedFontHist1D(base1,1.,1.2);
	base1->GetYaxis()->SetTitleOffset(1.5);
	base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize()*1.5);
	base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize()*1.5);
	base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize()*1.5);
	base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize()*1.5);
	base1->GetXaxis()->SetNdivisions(4,4,0);
	base1->GetYaxis()->SetNdivisions(5,5,0);
	base1->Draw();

	double dN_dt = 0.002;
    double sigma_coherent = 459.05; // nb
    double sigma_incoherent = 261.95; // nb
    double sigma_DIS = 50184.8; // nb
    double nEvents_coherent = h_phi_sartre_events->GetEntries(); // 6.36679 M
    double early_science_lumi = 8121.83; // nb^-1
    int rebin_width = 2;   

    h_t_MC->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    double simu_coherent_lumi = nEvents_coherent/sigma_coherent; // nb^-1 
    double ratio_coherent_phi_ES = early_science_lumi/simu_coherent_lumi;
    for (int i = 1; i <= h_t_MC->GetNbinsX(); ++i) 
    {
        double binError = h_t_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_MC->SetBinError(i, newError);
    }
	h_t_MC->SetLineStyle(1);   
	h_t_MC->SetLineWidth(1);   
	h_t_MC->SetLineColor(kBlack);
	h_t_MC->Draw("HISTsame");
	h_t_MC->Draw("same");

    h_t_REC_L->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    for (int i = 1; i <= h_t_REC_L->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_L->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_REC_L->SetBinError(i, newError);
    }
	h_t_REC_L->SetMarkerStyle(20); // method L RECO
	h_t_REC_L->SetMarkerColor(kP8Blue);
	h_t_REC_L->Draw("PEsame");

	h_t_REC_proj2->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_proj2->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/2)));
    for (int i = 1; i <= h_t_REC_proj2->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_proj2->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_REC_proj2->SetBinError(i, newError);
    }
	h_t_REC_proj2->SetMarkerStyle(42);
	h_t_REC_proj2->SetMarkerColor(kP10Violet);
	h_t_REC_proj2->SetLineColor(kP10Violet);
	//h_t_REC_proj2->Draw("PEsame");

    h_t_REC_proj3->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_proj3->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/3)));
    for (int i = 1; i <= h_t_REC_proj3->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_proj3->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_REC_proj3->SetBinError(i, newError);
    }
	h_t_REC_proj3->SetMarkerStyle(4);
	h_t_REC_proj3->SetMarkerColor(kP8Cyan);
	h_t_REC_proj3->SetLineColor(kP8Cyan);
	//h_t_REC_proj3->Draw("PEsame");

    h_t_REC_proj4->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_proj4->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/4)));
    for (int i = 1; i <= h_t_REC_proj4->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_proj4->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_REC_proj4->SetBinError(i, newError);
    }
	h_t_REC_proj4->SetMarkerStyle(24);
	h_t_REC_proj4->SetMarkerColor(kP8Orange);
	h_t_REC_proj4->SetLineColor(kP8Orange);
	//h_t_REC_proj4->Draw("PEsame");

    h_t_REC_proj6->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_proj6->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/6)));
    for (int i = 1; i <= h_t_REC_proj6->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_proj6->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_REC_proj6->SetBinError(i, newError);
    }
	h_t_REC_proj6->SetMarkerStyle(31);
	h_t_REC_proj6->SetMarkerColor(kP8Gray);
	h_t_REC_proj6->SetLineColor(kP8Gray);
	//h_t_REC_proj6->Draw("PEsame");

    h_t_REC_proj9->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_proj9->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/9)));
    for (int i = 1; i <= h_t_REC_proj9->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_proj9->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_REC_proj9->SetBinError(i, newError);
    }
	h_t_REC_proj9->SetMarkerStyle(33);
	h_t_REC_proj9->SetMarkerColor(kP8Green);
	h_t_REC_proj9->SetLineColor(kP8Green);
	//h_t_REC_proj9->Draw("PEsame");

    h_t_REC_proj12->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_proj12->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_t_REC_proj12->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_proj12->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_REC_proj12->SetBinError(i, newError);
    }
	h_t_REC_proj12->SetMarkerStyle(30);
	h_t_REC_proj12->SetMarkerColor(kP8Pink);
	h_t_REC_proj12->SetLineColor(kP8Pink);
	h_t_REC_proj12->Draw("PEsame");

    h_t_REC_proj16->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_proj16->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/16)));
    for (int i = 1; i <= h_t_REC_proj16->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_proj16->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_REC_proj16->SetBinError(i, newError);
    }
	h_t_REC_proj16->SetMarkerStyle(28);
	h_t_REC_proj16->SetMarkerColor(kP8Azure);
	h_t_REC_proj16->SetLineColor(kP8Azure);
	//h_t_REC_proj16->Draw("PEsame");

    h_t_REC_proj20->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_proj20->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/20)));
    for (int i = 1; i <= h_t_REC_proj20->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_proj20->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_REC_proj20->SetBinError(i, newError);
    }
	h_t_REC_proj20->SetMarkerStyle(13);
	h_t_REC_proj20->SetMarkerColor(kP8Red);
	h_t_REC_proj20->SetLineColor(kP8Red);
	//h_t_REC_proj20->Draw("PEsame");

    h_t_REC_proj24->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_proj24->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/24)));
    for (int i = 1; i <= h_t_REC_proj24->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_proj24->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_REC_proj24->SetBinError(i, newError);
    }
	h_t_REC_proj24->SetMarkerStyle(38);
	h_t_REC_proj24->SetMarkerColor(kP10Yellow);
	h_t_REC_proj24->SetLineColor(kP10Yellow);
	//h_t_REC_proj24->Draw("PEsame");
	
	TLatex* r44 = new TLatex(0.2, 0.84, "1<Q^{2}<10 GeV^{2}, 0.01 < y < 0.85");
	r44->SetNDC();
	r44->SetTextSize(20);
	r44->SetTextFont(43);
	r44->SetTextColor(kBlack);
	r44->Draw("same");
	
	TLatex* r44_0 = new TLatex(0.2, 0.79, "  |y_{"+vm_label+"}|<3.5, |M_{inv} #minus M_{"+vm_label+"}| < 0.02 GeV");
	r44_0->SetNDC();
	r44_0->SetTextSize(20);
	r44_0->SetTextFont(43);
	r44_0->SetTextColor(kBlack);
	r44_0->Draw("same");
	
    // Add labels
    TLatex* ep = new TLatex(0.18, 0.33, "ePIC");
    ep->SetNDC();
    ep->SetTextFont(62);
    ep->SetTextSize(0.030);
    ep->SetTextColor(kBlack);
    ep->Draw("same");

    TLatex* r421 = new TLatex(0.23, 0.33, " Simulation 25.10.2");
	r421->SetNDC();
	r421->SetTextSize(25);
	r421->SetTextFont(43);
	r421->SetTextColor(kBlack);
	r421->Draw("same");

	TLatex* label_5 = new TLatex(0.18, 0.23, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->SetTextSize(25);
	label_5->SetTextFont(43);
	label_5->SetTextColor(kBlack);
	label_5->Draw("same");

    TLatex* r4211 = new TLatex(0.18, 0.28, "eAu #rightarrow e'Au'#phi, 10x100 GeV");
	r4211->SetNDC();
	r4211->SetTextSize(25);
	r4211->SetTextFont(43);
	r4211->SetTextColor(kBlack);
	r4211->Draw("same");

    TLatex* r42111 = new TLatex(0.18, 0.18, "L_{int} = 1.6 fb^{-1}/A"); // ES lumi
	r42111->SetNDC();
	r42111->SetTextSize(25);
	r42111->SetTextFont(43);
	r42111->SetTextColor(kBlack);
	r42111->Draw("same");
	
	TLegend *w7 = new TLegend(0.58,0.65,0.73,0.85);
	w7->SetLineColor(kWhite);
	w7->SetFillColor(0);
	w7->SetTextSize(20);
	w7->SetTextFont(45);
	w7->AddEntry(h_t_MC, "Sartre "+vm_label+" MC ", "L");
	w7->AddEntry(h_t_REC_L, "Sartre "+vm_label+" Method L RECO", "P");
	//w7->AddEntry(h_t_REC_proj2, "Sartre "+vm_label+" RECO #theta_{max} = #pi/2", "P");
    //w7->AddEntry(h_t_REC_proj3, "Sartre "+vm_label+" RECO #theta_{max} = #pi/3", "P");
    //w7->AddEntry(h_t_REC_proj4, "Sartre "+vm_label+" RECO #theta_{max} = #pi/4", "P");
    //w7->AddEntry(h_t_REC_proj6, "Sartre "+vm_label+" RECO #theta_{max} = #pi/6", "P");
    //w7->AddEntry(h_t_REC_proj9, "Sartre "+vm_label+" RECO #theta_{max} = #pi/9", "P");
    w7->AddEntry(h_t_REC_proj12, "Sartre "+vm_label+" RECO #theta_{max} = #pi/12", "P");
    //w7->AddEntry(h_t_REC_proj16, "Sartre "+vm_label+" RECO #theta_{max} = #pi/16", "P");
    //w7->AddEntry(h_t_REC_proj20, "Sartre "+vm_label+" RECO #theta_{max} = #pi/20", "P");
    //w7->AddEntry(h_t_REC_proj24, "Sartre "+vm_label+" RECO #theta_{max} = #pi/24", "P");
	w7->Draw("same");

	c1->Print("./figures/plot_t_dist_ES_allMethods.pdf");

	//TFile* fout_ES = new TFile("plot_t_dist_ES.root","RECREATE");
    //c1->Write();
    //h_t_MC->Write();
    //h_t_REC_L->Write();
    //h_t_REC_proj12->Write();
    //fout_ES->Close();
}

void plot_DIS_noVetoes_preTDR()
{
	TFile* phi_t_file = TFile::Open("phi_sartre_25_10_2_eAu_10x100v5.root","READ");
	TH1D* h_t_MC = (TH1D*) phi_t_file->Get("h_t_MC");
	TH1D* h_t_REC_EEMC_cut = (TH1D*) phi_t_file->Get("h_t_REC_EEMC_cut");
	TH1D* h_t_REC_wRES_cut_pi12 = (TH1D*) phi_t_file->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_phi_sartre_events = (TH1D*) phi_t_file->Get("h_Nevents");

	TFile* DIS_noVetoes = TFile::Open("DIS_beagle_25_10_2_eAu_noVetoes_10x100v4.root","READ");
    TH1D* h_DIS_noVetoes_MC = (TH1D*) DIS_noVetoes->Get("h_t_MC");
    TH1D* h_DIS_noVetoes_L = (TH1D*) DIS_noVetoes->Get("h_t_REC_EEMC_cut");
    TH1D* h_DIS_noVetoes_proj = (TH1D*) DIS_noVetoes->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_DIS_events = (TH1D*) DIS_noVetoes->Get("h_Nevents");

    double dN_dt = 0.002;
    double sigma_coherent = 459.05; // nb
    double sigma_DIS = 50184.8; // nb
    double nEvents_coherent = h_phi_sartre_events->GetEntries();
    double nEvents_DIS = h_DIS_events->GetEntries();
    double preTDR_lumi = 50761.4; //nb^-1
    int rebin_width = 2;   

    TCanvas* c14_213 = new TCanvas("c14_213","c14_213",1,1,1200,800);
	gStyle->SetOptStat(0);
    gPad->SetLogy(1);
    h_t_MC->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
	h_t_MC->GetXaxis()->SetTitle("|t| [GeV/c]^{2}");
    h_t_MC->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    double simu_coherent_lumi = nEvents_coherent/sigma_coherent; // nb^-1 
    double ratio_coherent_phi_preTDR = preTDR_lumi/simu_coherent_lumi;
    for (int i = 1; i <= h_t_MC->GetNbinsX(); ++i) 
    {
        double binError = h_t_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_MC->SetBinError(i, newError);
    }
    h_t_MC->SetLineColor(kBlack);
	h_t_MC->Draw("HISTsame");

    h_t_REC_EEMC_cut->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    for (int i = 1; i <= h_t_REC_EEMC_cut->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_EEMC_cut->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_EEMC_cut->SetBinError(i, newError);
    }
	h_t_REC_EEMC_cut->SetMarkerStyle(20); // method L RECO
	h_t_REC_EEMC_cut->SetMarkerColor(kP8Blue);
	//h_t_REC_EEMC_cut->Draw("PEsame");

    h_t_REC_wRES_cut_pi12->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_wRES_cut_pi12->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_t_REC_wRES_cut_pi12->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_wRES_cut_pi12->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_wRES_cut_pi12->SetBinError(i, newError);
    }
	h_t_REC_wRES_cut_pi12->SetMarkerStyle(30);
	h_t_REC_wRES_cut_pi12->SetMarkerColor(kP8Pink);
    h_t_REC_wRES_cut_pi12->SetLineColor(kP8Pink);
	h_t_REC_wRES_cut_pi12->Draw("PEsame");

    h_DIS_noVetoes_MC->Rebin(rebin_width);
    h_DIS_noVetoes_MC->Scale(sigma_DIS*(1/nEvents_DIS)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_DIS_lumi = nEvents_DIS/sigma_DIS; // nb^-1 
    double ratio_DIS_preTDR = preTDR_lumi/simu_DIS_lumi;
    for (int i = 1; i <= h_DIS_noVetoes_MC->GetNbinsX(); ++i) 
    {
        double binError = h_DIS_noVetoes_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_DIS_preTDR);
        h_DIS_noVetoes_MC->SetBinError(i, newError);
    }
    h_DIS_noVetoes_MC->SetMarkerStyle(2);
    h_DIS_noVetoes_MC->SetMarkerColor(kP8Cyan);
    h_DIS_noVetoes_MC->SetLineColor(kP8Cyan);
    //h_DIS_noVetoes_MC->Draw("PEsame");

    h_DIS_noVetoes_L->Rebin(rebin_width);
    h_DIS_noVetoes_L->Scale(sigma_DIS*(1/nEvents_DIS)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_DIS_noVetoes_L->GetNbinsX(); ++i) 
    {
        double binError = h_DIS_noVetoes_L->GetBinError(i);
        double newError = binError / sqrt(ratio_DIS_preTDR);
        h_DIS_noVetoes_L->SetBinError(i, newError);
    }
    h_DIS_noVetoes_L->SetMarkerStyle(2);
    h_DIS_noVetoes_L->SetMarkerColor(kP8Cyan);
    h_DIS_noVetoes_L->SetLineColor(kP8Cyan);
    //h_DIS_noVetoes_L->Draw("PEsame");

    h_DIS_noVetoes_proj->Rebin(rebin_width);
    h_DIS_noVetoes_proj->Scale(sigma_DIS*(1/nEvents_DIS)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_DIS_noVetoes_proj->GetNbinsX(); ++i) 
    {
        double binError = h_DIS_noVetoes_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_DIS_preTDR);
        h_DIS_noVetoes_proj->SetBinError(i, newError);
    }
    h_DIS_noVetoes_proj->SetMarkerStyle(2);
    h_DIS_noVetoes_proj->SetMarkerColor(kP8Cyan);
    h_DIS_noVetoes_proj->SetLineColor(kP8Cyan);
    h_DIS_noVetoes_proj->Draw("PEsame");
        
	TLegend *w14_213 = new TLegend(0.6,0.68,0.72,0.9);
	w14_213->AddEntry(h_t_MC,"#phi: MC", "L");
    w14_213->AddEntry(h_t_REC_wRES_cut_pi12,"#phi: #theta_{max}= #pi/12", "P");    
    //w14_213->AddEntry(h_t_REC_EEMC_cut,"#phi: L", "P");
    //w14_213->AddEntry(h_DIS_noVetoes_MC,"DIS: MC", "P");
    //w14_213->AddEntry(h_DIS_noVetoes_L,"DIS: L", "P");
    w14_213->AddEntry(h_DIS_noVetoes_proj,"DIS", "P");
	w14_213->SetBorderSize(0);   
    w14_213->SetFillStyle(0);
    w14_213->SetTextSize(30);
	w14_213->SetTextFont(45);
    w14_213->Draw("same");	

    c14_213->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title14_213 = new TLatex();
    title14_213->SetNDC(); 
    title14_213->SetTextSize(0.05);
    title14_213->SetTextAlign(22);  
    //title14_213->DrawLatex(0.5, 0.97, "|t| Distribution with DIS No Vetoes");  

    TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

    // Add labels
    TLatex* ep = new TLatex(0.13, 0.85, "ePIC");
    ep->SetNDC();
    ep->SetTextFont(62);
    ep->SetTextSize(.038);
    ep->SetTextColor(kBlack);
    ep->Draw("same");

    TLatex* r421 = new TLatex(0.18, 0.85, " Simulation 25.10.2, 10x100 GeV");
	r421->SetNDC();
	r421->SetTextSize(30);
	r421->SetTextFont(43);
	r421->SetTextColor(kBlack);
	r421->Draw("same");

	TLatex* label_5 = new TLatex(0.17, 0.75, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->SetTextSize(30);
	label_5->SetTextFont(43);
	label_5->SetTextColor(kBlack);
	label_5->Draw("same");

    TLatex* r4211 = new TLatex(0.15, 0.8, "eAu #rightarrow e'Au'#phi");
	r4211->SetNDC();
	r4211->SetTextSize(30);
	r4211->SetTextFont(43);
	r4211->SetTextColor(kBlack);
	r4211->Draw("same"); 

    TLatex* r421111 = new TLatex(0.15, 0.15, "L_{int} = 10 fb^{-1}/A"); // preTDR lumi
	r421111->SetNDC();
	r421111->SetTextSize(25);
	r421111->SetTextFont(43);
	r421111->SetTextColor(kBlack);
	r421111->Draw("same");

    c14_213->Print("./figures/plot_t_DIS_noVetoes_preTDR_proj.pdf");

	//TFile* fout_preTDR = new TFile("plot_t_DIS_noVetoes_preTDR.root","RECREATE");
    //c14_213->Write();
    //h_t_MC->Write();
    //h_t_REC_EEMC_cut->Write();
    //h_t_REC_wRES_cut_pi12->Write();
    //h_DIS_noVetoes_MC->Write();
    //h_DIS_noVetoes_L->Write();
    //h_DIS_noVetoes_proj->Write();
    //fout_preTDR->Close();
}

void plot_DIS_noVetoes_ES()
{
	TFile* phi_t_file = TFile::Open("phi_sartre_25_10_2_eAu_10x100v5.root","READ");
	TH1D* h_t_MC = (TH1D*) phi_t_file->Get("h_t_MC");
	TH1D* h_t_REC_EEMC_cut = (TH1D*) phi_t_file->Get("h_t_REC_EEMC_cut");
	TH1D* h_t_REC_wRES_cut_pi12 = (TH1D*) phi_t_file->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_phi_sartre_events = (TH1D*) phi_t_file->Get("h_Nevents");

	TFile* DIS_noVetoes = TFile::Open("DIS_beagle_25_10_2_eAu_noVetoes_10x100v4.root","READ");
    TH1D* h_DIS_noVetoes_MC = (TH1D*) DIS_noVetoes->Get("h_t_MC");
    TH1D* h_DIS_noVetoes_L = (TH1D*) DIS_noVetoes->Get("h_t_REC_EEMC_cut");
    TH1D* h_DIS_noVetoes_proj = (TH1D*) DIS_noVetoes->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_DIS_events = (TH1D*) DIS_noVetoes->Get("h_Nevents");

    double dN_dt = 0.002;
    double sigma_coherent = 459.05; // nb
    double sigma_DIS = 50184.8; // nb
    double nEvents_coherent = h_phi_sartre_events->GetEntries();
    double nEvents_DIS = h_DIS_events->GetEntries();
    double early_science_lumi = 8121.83; // nb^-1
    double preTDR_lumi = 50761.4; //nb^-1
    int rebin_width = 2;   

    TCanvas* c14_213 = new TCanvas("c14_213","c14_213",1,1,1200,800);
	gStyle->SetOptStat(0);
    gPad->SetLogy(1);
    h_t_MC->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
	h_t_MC->GetXaxis()->SetTitle("|t| [GeV/c]^{2}");
    h_t_MC->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    double simu_coherent_lumi = nEvents_coherent/sigma_coherent; // nb^-1 
    double ratio_coherent_phi_ES = early_science_lumi/simu_coherent_lumi;
    for (int i = 1; i <= h_t_MC->GetNbinsX(); ++i) 
    {
        double binError = h_t_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_MC->SetBinError(i, newError);
    }
    h_t_MC->SetLineColor(kBlack);
	h_t_MC->Draw("HISTsame");

    h_t_REC_EEMC_cut->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    for (int i = 1; i <= h_t_REC_EEMC_cut->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_EEMC_cut->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_REC_EEMC_cut->SetBinError(i, newError);
    }
	h_t_REC_EEMC_cut->SetMarkerStyle(20); // method L RECO
	h_t_REC_EEMC_cut->SetMarkerColor(kP8Blue);
	//h_t_REC_EEMC_cut->Draw("PEsame");

    h_t_REC_wRES_cut_pi12->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_wRES_cut_pi12->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_t_REC_wRES_cut_pi12->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_wRES_cut_pi12->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_REC_wRES_cut_pi12->SetBinError(i, newError);
    }
	h_t_REC_wRES_cut_pi12->SetMarkerStyle(30);
	h_t_REC_wRES_cut_pi12->SetMarkerColor(kP8Pink);
    h_t_REC_wRES_cut_pi12->SetLineColor(kP8Pink);
	h_t_REC_wRES_cut_pi12->Draw("PEsame");

    h_DIS_noVetoes_MC->Rebin(rebin_width);
    h_DIS_noVetoes_MC->Scale(sigma_DIS*(1/nEvents_DIS)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_DIS_lumi = nEvents_DIS/sigma_DIS; // nb^-1 
    double ratio_DIS_ES = early_science_lumi/simu_DIS_lumi;
    for (int i = 1; i <= h_DIS_noVetoes_MC->GetNbinsX(); ++i) 
    {
        double binError = h_DIS_noVetoes_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_DIS_ES);
        h_DIS_noVetoes_MC->SetBinError(i, newError);
    }
    h_DIS_noVetoes_MC->SetMarkerStyle(2);
    h_DIS_noVetoes_MC->SetMarkerColor(kP8Cyan);
    h_DIS_noVetoes_MC->SetLineColor(kP8Cyan);
    //h_DIS_noVetoes_MC->Draw("PEsame");

    h_DIS_noVetoes_L->Rebin(rebin_width);
    h_DIS_noVetoes_L->Scale(sigma_DIS*(1/nEvents_DIS)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_DIS_noVetoes_L->GetNbinsX(); ++i) 
    {
        double binError = h_DIS_noVetoes_L->GetBinError(i);
        double newError = binError / sqrt(ratio_DIS_ES);
        h_DIS_noVetoes_L->SetBinError(i, newError);
    }
    h_DIS_noVetoes_L->SetMarkerStyle(2);
    h_DIS_noVetoes_L->SetMarkerColor(kP8Cyan);
    h_DIS_noVetoes_L->SetLineColor(kP8Cyan);
    //h_DIS_noVetoes_L->Draw("PEsame");

    h_DIS_noVetoes_proj->Rebin(rebin_width);
    h_DIS_noVetoes_proj->Scale(sigma_DIS*(1/nEvents_DIS)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_DIS_noVetoes_proj->GetNbinsX(); ++i) 
    {
        double binError = h_DIS_noVetoes_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_DIS_ES);
        h_DIS_noVetoes_proj->SetBinError(i, newError);
    }
    h_DIS_noVetoes_proj->SetMarkerStyle(2);
    h_DIS_noVetoes_proj->SetMarkerColor(kP8Cyan);
    h_DIS_noVetoes_proj->SetLineColor(kP8Cyan);
    h_DIS_noVetoes_proj->Draw("PEsame");

	TLegend *w14_213 = new TLegend(0.6,0.68,0.72,0.9);
	w14_213->AddEntry(h_t_MC,"#phi: MC", "L");
    w14_213->AddEntry(h_t_REC_wRES_cut_pi12,"#phi: #theta_{max}= #pi/12", "P");    
    //w14_213->AddEntry(h_t_REC_EEMC_cut,"#phi: L", "P");
    //w14_213->AddEntry(h_DIS_noVetoes_MC,"DIS: MC", "P");
    //w14_213->AddEntry(h_DIS_noVetoes_L,"DIS: L", "P");
    w14_213->AddEntry(h_DIS_noVetoes_proj,"DIS", "P");
	w14_213->SetBorderSize(0);   
    w14_213->SetFillStyle(0);
    w14_213->SetTextSize(30);
	w14_213->SetTextFont(45);
    w14_213->Draw("same");	

    c14_213->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title14_213 = new TLatex();
    title14_213->SetNDC(); 
    title14_213->SetTextSize(0.05);
    title14_213->SetTextAlign(22);  
    //title14_213->DrawLatex(0.5, 0.97, "|t| Distribution with DIS No Vetoes");  

    TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

    // Add labels
    TLatex* ep = new TLatex(0.13, 0.85, "ePIC");
    ep->SetNDC();
    ep->SetTextFont(62);
    ep->SetTextSize(.038);
    ep->SetTextColor(kBlack);
    ep->Draw("same");

    TLatex* r421 = new TLatex(0.18, 0.85, " Simulation 25.10.2, 10x100 GeV");
	r421->SetNDC();
	r421->SetTextSize(30);
	r421->SetTextFont(43);
	r421->SetTextColor(kBlack);
	r421->Draw("same");

	TLatex* label_5 = new TLatex(0.17, 0.75, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->SetTextSize(30);
	label_5->SetTextFont(43);
	label_5->SetTextColor(kBlack);
	label_5->Draw("same");

    TLatex* r4211 = new TLatex(0.15, 0.8, "eAu #rightarrow e'Au'#phi");
	r4211->SetNDC();
	r4211->SetTextSize(30);
	r4211->SetTextFont(43);
	r4211->SetTextColor(kBlack);
	r4211->Draw("same"); 

    TLatex* r42111 = new TLatex(0.15, 0.15, "L_{int} = 1.6 fb^{-1}/A"); // ES lumi
	r42111->SetNDC();
	r42111->SetTextSize(30);
	r42111->SetTextFont(43);
	r42111->SetTextColor(kBlack);
	//r42111->Draw("same");

    c14_213->Print("./figures/plot_t_DIS_noVetoes_ES_proj.pdf");

	//TFile* fout_ES = new TFile("plot_t_DIS_noVetoes_ES.root","RECREATE");
    //c14_213->Write();
    //h_t_MC->Write();
    //h_t_REC_EEMC_cut->Write();
    //h_t_REC_wRES_cut_pi12->Write();
    //h_DIS_noVetoes_MC->Write();
    //h_DIS_noVetoes_L->Write();
    //h_DIS_noVetoes_proj->Write();
    //fout_ES->Close();
}

void plot_DIS_allVetoes_preTDR()
{
	TFile* phi_t_file = TFile::Open("phi_sartre_25_10_2_eAu_10x100v5.root","READ");
	TH1D* h_t_MC = (TH1D*) phi_t_file->Get("h_t_MC");
	TH1D* h_t_REC_EEMC_cut = (TH1D*) phi_t_file->Get("h_t_REC_EEMC_cut");
	TH1D* h_t_REC_wRES_cut_pi12 = (TH1D*) phi_t_file->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_phi_sartre_events = (TH1D*) phi_t_file->Get("h_Nevents");

	TFile* DIS_allVetoes = TFile::Open("DIS_beagle_25_10_2_eAu_allVetoes_10x100v4.root","READ");
    TH1D* h_DIS_allVetoes_MC = (TH1D*) DIS_allVetoes->Get("h_t_MC");
    TH1D* h_DIS_allVetoes_L = (TH1D*) DIS_allVetoes->Get("h_t_REC_EEMC_cut");
    TH1D* h_DIS_allVetoes_proj = (TH1D*) DIS_allVetoes->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_DIS_events = (TH1D*) DIS_allVetoes->Get("h_Nevents");

	double dN_dt = 0.002;
    double sigma_coherent = 459.05; // nb
    double sigma_DIS = 50184.8; // nb
    double nEvents_coherent = h_phi_sartre_events->GetEntries();//6.2e6;
    double nEvents_DIS = h_DIS_events->GetEntries();//1.1e6;
    double preTDR_lumi = 50761.4; //nb^-1
    int rebin_width = 2;   

	TCanvas* c14_2131 = new TCanvas("c14_2131","c14_2131",1,1,1200,800);
    gPad->SetLogy(1);
    gStyle->SetOptStat(0);
    h_t_MC->GetYaxis()->SetTitle("dN/d|t| [GeV/c]^{-2}");
	h_t_MC->GetXaxis()->SetTitle("|t| [GeV/c]^{2}");
    h_t_MC->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    double simu_coherent_lumi = nEvents_coherent/sigma_coherent; // nb^-1 
    double ratio_coherent_phi_preTDR = preTDR_lumi/simu_coherent_lumi;
    for (int i = 1; i <= h_t_MC->GetNbinsX(); ++i) 
    {
        double binError = h_t_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_MC->SetBinError(i, newError);
    }
    h_t_MC->SetLineColor(kBlack);
	h_t_MC->Draw("HISTsame");

    h_t_REC_EEMC_cut->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    for (int i = 1; i <= h_t_REC_EEMC_cut->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_EEMC_cut->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_EEMC_cut->SetBinError(i, newError);
    }
	h_t_REC_EEMC_cut->SetMarkerStyle(20); // method L RECO
	h_t_REC_EEMC_cut->SetMarkerColor(kP8Blue);
	//h_t_REC_EEMC_cut->Draw("PEsame");

    h_t_REC_wRES_cut_pi12->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_wRES_cut_pi12->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_t_REC_wRES_cut_pi12->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_wRES_cut_pi12->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_wRES_cut_pi12->SetBinError(i, newError);
    }
	h_t_REC_wRES_cut_pi12->SetMarkerStyle(30);
	h_t_REC_wRES_cut_pi12->SetMarkerColor(kP8Pink);
    h_t_REC_wRES_cut_pi12->SetLineColor(kP8Pink);
	h_t_REC_wRES_cut_pi12->Draw("PEsame");

    h_DIS_allVetoes_MC->Rebin(rebin_width);
    h_DIS_allVetoes_MC->Scale(sigma_DIS*(1/nEvents_DIS)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_DIS_lumi = nEvents_DIS/sigma_DIS; // nb^-1 
    double ratio_DIS_preTDR = preTDR_lumi/simu_DIS_lumi;
    for (int i = 1; i <= h_DIS_allVetoes_MC->GetNbinsX(); ++i) 
    {
        double binError = h_DIS_allVetoes_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_DIS_preTDR);
        h_DIS_allVetoes_MC->SetBinError(i, newError);
    }
    h_DIS_allVetoes_MC->SetMarkerStyle(2);
    h_DIS_allVetoes_MC->SetMarkerColor(kP8Cyan);
    h_DIS_allVetoes_MC->SetLineColor(kP8Cyan);
    //h_DIS_allVetoes_MC->Draw("PEsame");

    h_DIS_allVetoes_L->Rebin(rebin_width);
    h_DIS_allVetoes_L->Scale(sigma_DIS*(1/nEvents_DIS)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_DIS_allVetoes_L->GetNbinsX(); ++i) 
    {
        double binError = h_DIS_allVetoes_L->GetBinError(i);
        double newError = binError / sqrt(ratio_DIS_preTDR);
        h_DIS_allVetoes_L->SetBinError(i, newError);
    }
    h_DIS_allVetoes_L->SetMarkerStyle(2);
    h_DIS_allVetoes_L->SetMarkerColor(kP8Cyan);
    h_DIS_allVetoes_L->SetLineColor(kP8Cyan);
    //h_DIS_allVetoes_L->Draw("PEsame");

    h_DIS_allVetoes_proj->Rebin(rebin_width);
    h_DIS_allVetoes_proj->Scale(sigma_DIS*(1/nEvents_DIS)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_DIS_allVetoes_proj->GetNbinsX(); ++i) 
    {
        double binError = h_DIS_allVetoes_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_DIS_preTDR);
        h_DIS_allVetoes_proj->SetBinError(i, newError);
    }
    h_DIS_allVetoes_proj->SetMarkerStyle(2);
    h_DIS_allVetoes_proj->SetMarkerColor(kP8Cyan);
    h_DIS_allVetoes_proj->SetLineColor(kP8Cyan);
    h_DIS_allVetoes_proj->Draw("PEsame");

	TLegend *w14_2131 = new TLegend(0.6,0.68,0.72,0.9);
	w14_2131->AddEntry(h_t_MC,"#phi: MC", "L");
    w14_2131->AddEntry(h_t_REC_wRES_cut_pi12,"#phi: #theta_{max}= #pi/12", "P");    
    //w14_2131->AddEntry(h_t_REC_EEMC_cut,"#phi: L", "P");
    //w14_2131->AddEntry(h_DIS_allVetoes_MC,"DIS: MC", "P");
    //w14_2131->AddEntry(h_DIS_allVetoes_L,"DIS: L", "P");
    w14_2131->AddEntry(h_DIS_allVetoes_proj,"DIS: #theta_{max}= #pi/12", "P");
    w14_2131->SetBorderSize(0);   
    w14_2131->SetFillStyle(0);
	w14_2131->SetTextSize(30);
	w14_2131->SetTextFont(45);
    w14_2131->Draw("same");	

    c14_2131->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title14_2131 = new TLatex();
    title14_2131->SetNDC(); 
    title14_2131->SetTextSize(0.05);
    title14_2131->SetTextAlign(22);  
    //title14_2131->DrawLatex(0.5, 0.97, "|t| Distribution with DIS All Vetoes");  

    TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

    // Add labels
    TLatex* ep = new TLatex(0.13, 0.85, "ePIC");
    ep->SetNDC();
    ep->SetTextFont(62);
    ep->SetTextSize(.038);
    ep->SetTextColor(kBlack);
    ep->Draw("same");

    TLatex* r421 = new TLatex(0.18, 0.85, " Simulation 25.10.2, 10x100 GeV");
	r421->SetNDC();
	r421->SetTextSize(30);
	r421->SetTextFont(43);
	r421->SetTextColor(kBlack);
	r421->Draw("same");

	TLatex* label_5 = new TLatex(0.17, 0.75, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->SetTextSize(30);
	label_5->SetTextFont(43);
	label_5->SetTextColor(kBlack);
	label_5->Draw("same");

    TLatex* r4211 = new TLatex(0.15, 0.8, "eAu #rightarrow e'Au'#phi");
	r4211->SetNDC();
	r4211->SetTextSize(30);
	r4211->SetTextFont(43);
	r4211->SetTextColor(kBlack);
	r4211->Draw("same"); 

    TLatex* r421111 = new TLatex(0.15, 0.15, "L_{int} = 10 fb^{-1}/A"); // preTDR lumi
	r421111->SetNDC();
	r421111->SetTextSize(25);
	r421111->SetTextFont(43);
	r421111->SetTextColor(kBlack);
	r421111->Draw("same");

	
    c14_2131->Print("./figures/plot_t_DIS_allvetoes_preTDR_proj.pdf");

	//TFile* fout_preTDR = new TFile("plot_t_DIS_allvetoes_L_preTDR.root","RECREATE");
    //c14_2131->Write();
    //h_t_MC->Write();
    //h_t_REC_EEMC_cut->Write();
    //h_t_REC_wRES_cut_pi12->Write();
    //h_DIS_allVetoes_MC->Write();
    //h_DIS_allVetoes_L->Write();
    //h_DIS_allVetoes_proj->Write();
    //fout_preTDR->Close();

}

void plot_DIS_allVetoes_ES()
{
	TFile* phi_t_file = TFile::Open("phi_sartre_25_10_2_eAu_10x100v5.root","READ");
	TH1D* h_t_MC = (TH1D*) phi_t_file->Get("h_t_MC");
	TH1D* h_t_REC_EEMC_cut = (TH1D*) phi_t_file->Get("h_t_REC_EEMC_cut");
	TH1D* h_t_REC_wRES_cut_pi12 = (TH1D*) phi_t_file->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_phi_sartre_events = (TH1D*) phi_t_file->Get("h_Nevents");

	TFile* DIS_allVetoes = TFile::Open("DIS_beagle_25_10_2_eAu_allVetoes_10x100v4.root","READ");
    TH1D* h_DIS_allVetoes_MC = (TH1D*) DIS_allVetoes->Get("h_t_MC");
    TH1D* h_DIS_allVetoes_L = (TH1D*) DIS_allVetoes->Get("h_t_REC_EEMC_cut");
    TH1D* h_DIS_allVetoes_proj = (TH1D*) DIS_allVetoes->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_DIS_events = (TH1D*) DIS_allVetoes->Get("h_Nevents");

	double dN_dt = 0.002;
    double sigma_coherent = 459.05; // nb
    double sigma_DIS = 50184.8; // nb
    double nEvents_coherent = h_phi_sartre_events->GetEntries();//6.2e6;
    double nEvents_DIS = h_DIS_events->GetEntries();//1.1e6;
    double early_science_lumi = 8121.83; // nb^-1
    int rebin_width = 2;   

	TCanvas* c14_2131 = new TCanvas("c14_2131","c14_2131",1,1,1200,800);
    gPad->SetLogy(1);
    gStyle->SetOptStat(0);
    h_t_MC->GetYaxis()->SetTitle("dN/d|t| [GeV/c]^{-2}");
	h_t_MC->GetXaxis()->SetTitle("|t| [GeV/c]^{2}");
    h_t_MC->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    double simu_coherent_lumi = nEvents_coherent/sigma_coherent; // nb^-1 
    double ratio_coherent_phi_ES = early_science_lumi/simu_coherent_lumi;
    for (int i = 1; i <= h_t_MC->GetNbinsX(); ++i) 
    {
        double binError = h_t_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_MC->SetBinError(i, newError);
    }
    h_t_MC->SetLineColor(kBlack);
	h_t_MC->Draw("HISTsame");

    h_t_REC_EEMC_cut->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    for (int i = 1; i <= h_t_REC_EEMC_cut->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_EEMC_cut->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_REC_EEMC_cut->SetBinError(i, newError);
    }
	h_t_REC_EEMC_cut->SetMarkerStyle(20); // method L RECO
	h_t_REC_EEMC_cut->SetMarkerColor(kP8Blue);
	//h_t_REC_EEMC_cut->Draw("PEsame");

    h_t_REC_wRES_cut_pi12->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_wRES_cut_pi12->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_t_REC_wRES_cut_pi12->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_wRES_cut_pi12->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_REC_wRES_cut_pi12->SetBinError(i, newError);
    }
	h_t_REC_wRES_cut_pi12->SetMarkerStyle(30);
	h_t_REC_wRES_cut_pi12->SetMarkerColor(kP8Pink);
    h_t_REC_wRES_cut_pi12->SetLineColor(kP8Pink);
	h_t_REC_wRES_cut_pi12->Draw("PEsame");

    h_DIS_allVetoes_MC->Rebin(rebin_width);
    h_DIS_allVetoes_MC->Scale(sigma_DIS*(1/nEvents_DIS)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_DIS_lumi = nEvents_DIS/sigma_DIS; // nb^-1 
    double ratio_DIS_ES = early_science_lumi/simu_DIS_lumi;
    for (int i = 1; i <= h_DIS_allVetoes_MC->GetNbinsX(); ++i) 
    {
        double binError = h_DIS_allVetoes_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_DIS_ES);
        h_DIS_allVetoes_MC->SetBinError(i, newError);
    }
    h_DIS_allVetoes_MC->SetMarkerStyle(2);
    h_DIS_allVetoes_MC->SetMarkerColor(kP8Cyan);
    h_DIS_allVetoes_MC->SetLineColor(kP8Cyan);
    //h_DIS_allVetoes_MC->Draw("PEsame");

    h_DIS_allVetoes_L->Rebin(rebin_width);
    h_DIS_allVetoes_L->Scale(sigma_DIS*(1/nEvents_DIS)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_DIS_allVetoes_L->GetNbinsX(); ++i) 
    {
        double binError = h_DIS_allVetoes_L->GetBinError(i);
        double newError = binError / sqrt(ratio_DIS_ES);
        h_DIS_allVetoes_L->SetBinError(i, newError);
    }
    h_DIS_allVetoes_L->SetMarkerStyle(2);
    h_DIS_allVetoes_L->SetMarkerColor(kP8Cyan);
    h_DIS_allVetoes_L->SetLineColor(kP8Cyan);
    //h_DIS_allVetoes_L->Draw("PEsame");

    h_DIS_allVetoes_proj->Rebin(rebin_width);
    h_DIS_allVetoes_proj->Scale(sigma_DIS*(1/nEvents_DIS)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_DIS_allVetoes_proj->GetNbinsX(); ++i) 
    {
        double binError = h_DIS_allVetoes_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_DIS_ES);
        h_DIS_allVetoes_proj->SetBinError(i, newError);
    }
    h_DIS_allVetoes_proj->SetMarkerStyle(2);
    h_DIS_allVetoes_proj->SetMarkerColor(kP8Cyan);
    h_DIS_allVetoes_proj->SetLineColor(kP8Cyan);
    h_DIS_allVetoes_proj->Draw("PEsame");
        
	TLegend *w14_2131 = new TLegend(0.6,0.68,0.72,0.9);
	w14_2131->AddEntry(h_t_MC,"#phi: MC", "L");
    w14_2131->AddEntry(h_t_REC_wRES_cut_pi12,"#phi: #theta_{max}= #pi/12", "P");    
    //w14_2131->AddEntry(h_t_REC_EEMC_cut,"#phi: L", "P");
    //w14_2131->AddEntry(h_DIS_allVetoes_MC,"DIS: MC", "P");
    //w14_2131->AddEntry(h_DIS_allVetoes_L,"DIS: L", "P");
    w14_2131->AddEntry(h_DIS_allVetoes_proj,"DIS: #theta_{max}= #pi/12", "P");
    w14_2131->SetBorderSize(0);   
    w14_2131->SetFillStyle(0);
	w14_2131->SetTextSize(30);
	w14_2131->SetTextFont(45);
    w14_2131->Draw("same");	

    c14_2131->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title14_2131 = new TLatex();
    title14_2131->SetNDC(); 
    title14_2131->SetTextSize(0.05);
    title14_2131->SetTextAlign(22);  
    //title14_2131->DrawLatex(0.5, 0.97, "|t| Distribution with DIS All Vetoes");  

    TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

    // Add labels
    TLatex* ep = new TLatex(0.13, 0.85, "ePIC");
    ep->SetNDC();
    ep->SetTextFont(62);
    ep->SetTextSize(.038);
    ep->SetTextColor(kBlack);
    ep->Draw("same");

    TLatex* r421 = new TLatex(0.18, 0.85, " Simulation 25.10.2, 10x100 GeV");
	r421->SetNDC();
	r421->SetTextSize(30);
	r421->SetTextFont(43);
	r421->SetTextColor(kBlack);
	r421->Draw("same");

	TLatex* label_5 = new TLatex(0.17, 0.75, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->SetTextSize(30);
	label_5->SetTextFont(43);
	label_5->SetTextColor(kBlack);
	label_5->Draw("same");

    TLatex* r4211 = new TLatex(0.15, 0.8, "eAu #rightarrow e'Au'#phi");
	r4211->SetNDC();
	r4211->SetTextSize(30);
	r4211->SetTextFont(43);
	r4211->SetTextColor(kBlack);
	r4211->Draw("same"); 

    TLatex* r42111 = new TLatex(0.15, 0.15, "L_{int} = 1.6 fb^{-1}/A"); // ES lumi
	r42111->SetNDC();
	r42111->SetTextSize(30);
	r42111->SetTextFont(43);
	r42111->SetTextColor(kBlack);
	r42111->Draw("same");
	
    c14_2131->Print("./figures/plot_t_DIS_allvetoes_ES_proj.pdf");

	//TFile* fout_ES = new TFile("plot_t_DIS_allvetoes_L_ES.root","RECREATE");
    //c14_2131->Write();
    //h_t_MC->Write();
    //h_t_REC_EEMC_cut->Write();
    //h_t_REC_wRES_cut_pi12->Write();
    //h_DIS_allVetoes_MC->Write();
    //h_DIS_allVetoes_L->Write();
    //h_DIS_allVetoes_proj->Write();
    //fout_ES->Close();

}

void plot_DIS_allDetectorVetoes_preTDR()
{
	TFile* phi_t_file = TFile::Open("phi_sartre_25_10_2_eAu_10x100v5.root","READ");
	TH1D* h_t_MC = (TH1D*) phi_t_file->Get("h_t_MC");
	TH1D* h_t_REC_EEMC_cut = (TH1D*) phi_t_file->Get("h_t_REC_EEMC_cut");
	TH1D* h_t_REC_wRES_cut_pi12 = (TH1D*) phi_t_file->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_phi_sartre_events = (TH1D*) phi_t_file->Get("h_Nevents");

	TFile* DIS_onlyDetectors = TFile::Open("DIS_beagle_25_10_2_eAu_allDetectorVetoes_10x100v4.root");
    TH1D* h_DIS_onlyDetectors_MC = (TH1D*) DIS_onlyDetectors->Get("h_t_MC");
    TH1D* h_DIS_onlyDetectors_L = (TH1D*) DIS_onlyDetectors->Get("h_t_REC_EEMC_cut");
    TH1D* h_DIS_onlyDetectors_proj = (TH1D*) DIS_onlyDetectors->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_DIS_events = (TH1D*) DIS_onlyDetectors->Get("h_Nevents");

	double dN_dt = 0.002;
    double sigma_coherent = 459.05; // nb
    double sigma_DIS = 50184.8; // nb
    double nEvents_coherent = h_phi_sartre_events->GetEntries();//6.2e6;
    double nEvents_DIS = h_DIS_events->GetEntries();//1.1e6;
    double preTDR_lumi = 50761.4; // nb^-1
    int rebin_width = 2; 

	TCanvas* c14_2132 = new TCanvas("c14_2132","c14_2132",1,1,1200,800);
    gPad->SetLogy(1);
    gStyle->SetOptStat(0);
    h_t_MC->GetYaxis()->SetTitle("dN/d|t| [GeV/c]^{-2}");
	h_t_MC->GetXaxis()->SetTitle("|t| [GeV/c]^{2}");
    h_t_MC->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    double simu_coherent_lumi = nEvents_coherent/sigma_coherent; // nb^-1 
    double ratio_coherent_phi_preTDR = preTDR_lumi/simu_coherent_lumi;
    for (int i = 1; i <= h_t_MC->GetNbinsX(); ++i) 
    {
        double binError = h_t_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_MC->SetBinError(i, newError);
    }
    h_t_MC->SetLineColor(kBlack);
	h_t_MC->Draw("HISTsame");

    h_t_REC_EEMC_cut->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    for (int i = 1; i <= h_t_REC_EEMC_cut->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_EEMC_cut->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_EEMC_cut->SetBinError(i, newError);
    }
	h_t_REC_EEMC_cut->SetMarkerStyle(20); // method L RECO
	h_t_REC_EEMC_cut->SetMarkerColor(kP8Blue);
	//h_t_REC_EEMC_cut->Draw("PEsame");
    
    h_t_REC_wRES_cut_pi12->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_wRES_cut_pi12->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_t_REC_wRES_cut_pi12->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_wRES_cut_pi12->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_wRES_cut_pi12->SetBinError(i, newError);
    }
	h_t_REC_wRES_cut_pi12->SetMarkerStyle(30);
	h_t_REC_wRES_cut_pi12->SetMarkerColor(kP8Pink);
    h_t_REC_wRES_cut_pi12->SetLineColor(kP8Pink);
	h_t_REC_wRES_cut_pi12->Draw("PEsame");

    h_DIS_onlyDetectors_MC->Rebin(rebin_width);
    h_DIS_onlyDetectors_MC->Scale(sigma_DIS*(1/nEvents_DIS)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_DIS_lumi = nEvents_DIS/sigma_DIS; // nb^-1 
    double ratio_DIS_preTDR = preTDR_lumi/simu_DIS_lumi;
    for (int i = 1; i <= h_DIS_onlyDetectors_MC->GetNbinsX(); ++i) 
    {
        double binError = h_DIS_onlyDetectors_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_DIS_preTDR);
        h_DIS_onlyDetectors_MC->SetBinError(i, newError);
    }
    h_DIS_onlyDetectors_MC->SetMarkerStyle(2);
    h_DIS_onlyDetectors_MC->SetMarkerColor(kP8Cyan);
    h_DIS_onlyDetectors_MC->SetLineColor(kP8Cyan);
    //h_DIS_onlyDetectors_MC->Draw("PEsame");

    h_DIS_onlyDetectors_L->Rebin(rebin_width);
    h_DIS_onlyDetectors_L->Scale(sigma_DIS*(1/nEvents_DIS)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_DIS_onlyDetectors_L->GetNbinsX(); ++i) 
    {
        double binError = h_DIS_onlyDetectors_L->GetBinError(i);
        double newError = binError / sqrt(ratio_DIS_preTDR);
        h_DIS_onlyDetectors_L->SetBinError(i, newError);
    }
    h_DIS_onlyDetectors_L->SetMarkerStyle(2);
    h_DIS_onlyDetectors_L->SetMarkerColor(kP8Cyan);
    h_DIS_onlyDetectors_L->SetLineColor(kP8Cyan);
    //h_DIS_onlyDetectors_L->Draw("PEsame");

    h_DIS_onlyDetectors_proj->Rebin(rebin_width);
    h_DIS_onlyDetectors_proj->Scale(sigma_DIS*(1/nEvents_DIS)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_DIS_onlyDetectors_proj->GetNbinsX(); ++i) 
    {
        double binError = h_DIS_onlyDetectors_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_DIS_preTDR);
        h_DIS_onlyDetectors_proj->SetBinError(i, newError);
    }
    h_DIS_onlyDetectors_proj->SetMarkerStyle(2);
    h_DIS_onlyDetectors_proj->SetMarkerColor(kP8Cyan);
    h_DIS_onlyDetectors_proj->SetLineColor(kP8Cyan);
    h_DIS_onlyDetectors_proj->Draw("PEsame");

	TLegend *w14_2132 = new TLegend(0.6,0.68,0.72,0.9);
	w14_2132->AddEntry(h_t_MC,"#phi: MC", "L");
    w14_2132->AddEntry(h_t_REC_wRES_cut_pi12,"#phi: #theta_{max}= #pi/12", "P");    
    //w14_2132->AddEntry(h_t_REC_EEMC_cut,"#phi: L", "P");
    //w14_2132->AddEntry(h_DIS_onlyDetectors_MC,"DIS: MC", "P");
    //w14_2132->AddEntry(h_DIS_onlyDetectors_L,"DIS: L", "P");
    w14_2132->AddEntry(h_DIS_onlyDetectors_proj,"DIS: #theta_{max}= #pi/12", "P");
    w14_2132->SetBorderSize(0);   
    w14_2132->SetFillStyle(0);
	w14_2132->SetTextSize(30);
	w14_2132->SetTextFont(45);
    w14_2132->Draw("same");	

    c14_2132->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title14_2132 = new TLatex();
    title14_2132->SetNDC(); 
    title14_2132->SetTextSize(0.05);
    title14_2132->SetTextAlign(22);  
    //title14_2132->DrawLatex(0.5, 0.97, "|t| Distribution with DIS Only Detector Vetoes");  

    TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

    // Add labels
    TLatex* ep = new TLatex(0.13, 0.85, "ePIC");
    ep->SetNDC();
    ep->SetTextFont(62);
    ep->SetTextSize(.038);
    ep->SetTextColor(kBlack);
    ep->Draw("same");

    TLatex* r421 = new TLatex(0.18, 0.85, " Simulation 25.10.2, 10x100 GeV");
	r421->SetNDC();
	r421->SetTextSize(30);
	r421->SetTextFont(43);
	r421->SetTextColor(kBlack);
	r421->Draw("same");

	TLatex* label_5 = new TLatex(0.17, 0.75, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->SetTextSize(30);
	label_5->SetTextFont(43);
	label_5->SetTextColor(kBlack);
	label_5->Draw("same");

    TLatex* r4211 = new TLatex(0.15, 0.8, "eAu #rightarrow e'Au'#phi");
	r4211->SetNDC();
	r4211->SetTextSize(30);
	r4211->SetTextFont(43);
	r4211->SetTextColor(kBlack);
	r4211->Draw("same"); 

    TLatex* r421111 = new TLatex(0.15, 0.15, "L_{int} = 10 fb^{-1}/A"); // preTDR lumi
	r421111->SetNDC();
	r421111->SetTextSize(25);
	r421111->SetTextFont(43);
	r421111->SetTextColor(kBlack);
	r421111->Draw("same");

    c14_2132->Print("./figures/plot_t_DIS_detectorsOnly_preTDR_proj.pdf");

	//TFile* fout_preTDR = new TFile("plot_t_DIS_detectorsOnly_preTDR.root","RECREATE");
    //c14_2132->Write();
    //h_t_MC->Write();
    //h_t_REC_EEMC_cut->Write();
    //h_t_REC_wRES_cut_pi12->Write();
    //h_DIS_onlyDetectors_MC->Write();
    //h_DIS_onlyDetectors_L->Write();
    //h_DIS_onlyDetectors_proj->Write();
    //fout_preTDR->Close();
}

void plot_DIS_allDetectorVetoes_ES()
{
	TFile* phi_t_file = TFile::Open("phi_sartre_25_10_2_eAu_10x100v3.root","READ");
	TH1D* h_t_MC = (TH1D*) phi_t_file->Get("h_t_MC");
	TH1D* h_t_REC_EEMC_cut = (TH1D*) phi_t_file->Get("h_t_REC_EEMC_cut");
	TH1D* h_t_REC_wRES_cut_pi12 = (TH1D*) phi_t_file->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_phi_sartre_events = (TH1D*) phi_t_file->Get("h_Nevents");

	TFile* DIS_onlyDetectors = TFile::Open("DIS_beagle_25_10_2_eAu_allDetectorVetoes_10x100v3.root");
    TH1D* h_DIS_onlyDetectors_MC = (TH1D*) DIS_onlyDetectors->Get("h_t_MC");
    TH1D* h_DIS_onlyDetectors_L = (TH1D*) DIS_onlyDetectors->Get("h_t_REC_EEMC_cut");
    TH1D* h_DIS_onlyDetectors_proj = (TH1D*) DIS_onlyDetectors->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_DIS_events = (TH1D*) DIS_onlyDetectors->Get("h_Nevents");

	double dN_dt = 0.002;
    double sigma_coherent = 459.05; // nb
    double sigma_DIS = 50184.8; // nb
    double nEvents_coherent = h_phi_sartre_events->GetEntries();//6.2e6;
    double nEvents_DIS = h_DIS_events->GetEntries();//1.1e6;
    double early_science_lumi = 8121.83; // nb^-1
    double preTDR_lumi = 50761.4; // nb^-1
    int rebin_width = 2; 

	TCanvas* c14_2132 = new TCanvas("c14_2132","c14_2132",1,1,1200,800);
    gPad->SetLogy(1);
    gStyle->SetOptStat(0);
    h_t_MC->GetYaxis()->SetTitle("dN/d|t| [GeV/c]^{-2}");
	h_t_MC->GetXaxis()->SetTitle("|t| [GeV/c]^{2}");
    h_t_MC->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    double simu_coherent_lumi = nEvents_coherent/sigma_coherent; // nb^-1 
    double ratio_coherent_phi_ES = early_science_lumi/simu_coherent_lumi;
    for (int i = 1; i <= h_t_MC->GetNbinsX(); ++i) 
    {
        double binError = h_t_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_MC->SetBinError(i, newError);
    }
    h_t_MC->SetLineColor(kBlack);
	h_t_MC->Draw("HISTsame");

    h_t_REC_EEMC_cut->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    for (int i = 1; i <= h_t_REC_EEMC_cut->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_EEMC_cut->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_REC_EEMC_cut->SetBinError(i, newError);
    }
	h_t_REC_EEMC_cut->SetMarkerStyle(20); // method L RECO
	h_t_REC_EEMC_cut->SetMarkerColor(kP8Blue);
	//h_t_REC_EEMC_cut->Draw("PEsame");
    
    h_t_REC_wRES_cut_pi12->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_wRES_cut_pi12->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_t_REC_wRES_cut_pi12->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_wRES_cut_pi12->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_REC_wRES_cut_pi12->SetBinError(i, newError);
    }
	h_t_REC_wRES_cut_pi12->SetMarkerStyle(30);
	h_t_REC_wRES_cut_pi12->SetMarkerColor(kP8Pink);
    h_t_REC_wRES_cut_pi12->SetLineColor(kP8Pink);
	h_t_REC_wRES_cut_pi12->Draw("PEsame");

    h_DIS_onlyDetectors_MC->Rebin(rebin_width);
    h_DIS_onlyDetectors_MC->Scale(sigma_DIS*(1/nEvents_DIS)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_DIS_lumi = nEvents_DIS/sigma_DIS; // nb^-1 
    double ratio_DIS_ES = early_science_lumi/simu_DIS_lumi;
    for (int i = 1; i <= h_DIS_onlyDetectors_MC->GetNbinsX(); ++i) 
    {
        double binError = h_DIS_onlyDetectors_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_DIS_ES);
        h_DIS_onlyDetectors_MC->SetBinError(i, newError);
    }
    h_DIS_onlyDetectors_MC->SetMarkerStyle(2);
    h_DIS_onlyDetectors_MC->SetMarkerColor(kP8Cyan);
    h_DIS_onlyDetectors_MC->SetLineColor(kP8Cyan);
    //h_DIS_onlyDetectors_MC->Draw("PEsame");

    h_DIS_onlyDetectors_L->Rebin(rebin_width);
    h_DIS_onlyDetectors_L->Scale(sigma_DIS*(1/nEvents_DIS)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_DIS_onlyDetectors_L->GetNbinsX(); ++i) 
    {
        double binError = h_DIS_onlyDetectors_L->GetBinError(i);
        double newError = binError / sqrt(ratio_DIS_ES);
        h_DIS_onlyDetectors_L->SetBinError(i, newError);
    }
    h_DIS_onlyDetectors_L->SetMarkerStyle(2);
    h_DIS_onlyDetectors_L->SetMarkerColor(kP8Cyan);
    h_DIS_onlyDetectors_L->SetLineColor(kP8Cyan);
    //h_DIS_onlyDetectors_L->Draw("PEsame");

    h_DIS_onlyDetectors_proj->Rebin(rebin_width);
    h_DIS_onlyDetectors_proj->Scale(sigma_DIS*(1/nEvents_DIS)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_DIS_onlyDetectors_proj->GetNbinsX(); ++i) 
    {
        double binError = h_DIS_onlyDetectors_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_DIS_ES);
        h_DIS_onlyDetectors_proj->SetBinError(i, newError);
    }
    h_DIS_onlyDetectors_proj->SetMarkerStyle(2);
    h_DIS_onlyDetectors_proj->SetMarkerColor(kP8Cyan);
    h_DIS_onlyDetectors_proj->SetLineColor(kP8Cyan);
    h_DIS_onlyDetectors_proj->Draw("PEsame");
    
	TLegend *w14_2132 = new TLegend(0.6,0.68,0.72,0.9);
	w14_2132->AddEntry(h_t_MC,"#phi: MC", "L");
    w14_2132->AddEntry(h_t_REC_wRES_cut_pi12,"#phi: #theta_{max}= #pi/12", "P");    
    //w14_2132->AddEntry(h_t_REC_EEMC_cut,"#phi: L", "P");
    //w14_2132->AddEntry(h_DIS_onlyDetectors_MC,"DIS: MC", "P");
    //w14_2132->AddEntry(h_DIS_onlyDetectors_L,"DIS: L", "P");
    w14_2132->AddEntry(h_DIS_onlyDetectors_proj,"DIS: #theta_{max}= #pi/12", "P");
    w14_2132->SetBorderSize(0);   
    w14_2132->SetFillStyle(0);
	w14_2132->SetTextSize(30);
	w14_2132->SetTextFont(45);
    w14_2132->Draw("same");	

    c14_2132->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title14_2132 = new TLatex();
    title14_2132->SetNDC(); 
    title14_2132->SetTextSize(0.05);
    title14_2132->SetTextAlign(22);  
    //title14_2132->DrawLatex(0.5, 0.97, "|t| Distribution with DIS Only Detector Vetoes");  

    TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

    // Add labels
    TLatex* ep = new TLatex(0.13, 0.85, "ePIC");
    ep->SetNDC();
    ep->SetTextFont(62);
    ep->SetTextSize(.038);
    ep->SetTextColor(kBlack);
    ep->Draw("same");

    TLatex* r421 = new TLatex(0.18, 0.85, " Simulation 25.10.2, 10x100 GeV");
	r421->SetNDC();
	r421->SetTextSize(30);
	r421->SetTextFont(43);
	r421->SetTextColor(kBlack);
	r421->Draw("same");

	TLatex* label_5 = new TLatex(0.17, 0.75, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->SetTextSize(30);
	label_5->SetTextFont(43);
	label_5->SetTextColor(kBlack);
	label_5->Draw("same");

    TLatex* r4211 = new TLatex(0.15, 0.8, "eAu #rightarrow e'Au'#phi");
	r4211->SetNDC();
	r4211->SetTextSize(30);
	r4211->SetTextFont(43);
	r4211->SetTextColor(kBlack);
	r4211->Draw("same"); 

    TLatex* r421111 = new TLatex(0.15, 0.15, "L_{int} = 10 fb^{-1}/A"); // preTDR lumi
	r421111->SetNDC();
	r421111->SetTextSize(25);
	r421111->SetTextFont(43);
	r421111->SetTextColor(kBlack);
	r421111->Draw("same");

    c14_2132->Print("./figures/plot_t_DIS_detectorsOnly_ES_proj.pdf");

	//TFile* fout_ES = new TFile("plot_t_DIS_detectorsOnly_ES.root","RECREATE");
    //c14_2132->Write();
    //h_t_MC->Write();
    //h_t_REC_EEMC_cut->Write();
    //h_t_REC_wRES_cut_pi12->Write();
    //h_DIS_onlyDetectors_MC->Write();
    //h_DIS_onlyDetectors_L->Write();
    //h_DIS_onlyDetectors_proj->Write();
    //fout_ES->Close();
}

void plot_DIS_etaVetoes_preTDR()
{
	TFile* phi_t_file = TFile::Open("phi_sartre_25_10_2_eAu_10x100v5.root","READ");
	TH1D* h_t_MC = (TH1D*) phi_t_file->Get("h_t_MC");
	TH1D* h_t_REC_EEMC_cut = (TH1D*) phi_t_file->Get("h_t_REC_EEMC_cut");
	TH1D* h_t_REC_wRES_cut_pi12 = (TH1D*) phi_t_file->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_phi_sartre_events = (TH1D*) phi_t_file->Get("h_Nevents");

	TFile* DIS_etaOnly = TFile::Open("DIS_beagle_25_10_2_eAu_etaVetoes_10x100v4.root");
    TH1D* h_DIS_etaOnly_MC = (TH1D*) DIS_etaOnly->Get("h_t_MC");
    TH1D* h_DIS_etaOnly_L = (TH1D*) DIS_etaOnly->Get("h_t_REC_EEMC_cut");
    TH1D* h_DIS_etaOnly_proj = (TH1D*) DIS_etaOnly->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_DIS_events = (TH1D*) DIS_etaOnly->Get("h_Nevents");

	double dN_dt = 0.002;
    double sigma_coherent = 459.05; // nb
    double sigma_DIS = 50184.8; // nb
    double nEvents_coherent = h_phi_sartre_events->GetEntries();//6.2e6;
    double nEvents_DIS = h_DIS_events->GetEntries();//1.1e6;
    double preTDR_lumi = 50761.4; // nb^-1
    int rebin_width = 2; 

	TCanvas* c14_2133 = new TCanvas("c14_2133","c14_2133",1,1,1200,800);
    gPad->SetLogy(1);
    gStyle->SetOptStat(0);
    h_t_MC->GetYaxis()->SetTitle("dN/d|t| [GeV/c]^{-2}");
	h_t_MC->GetXaxis()->SetTitle("|t| [GeV/c]^{2}");
    h_t_MC->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    double simu_coherent_lumi = nEvents_coherent/sigma_coherent; // nb^-1 
    double ratio_coherent_phi_preTDR = preTDR_lumi/simu_coherent_lumi;
    for (int i = 1; i <= h_t_MC->GetNbinsX(); ++i) 
    {
        double binError = h_t_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_MC->SetBinError(i, newError);
    }
    h_t_MC->SetLineColor(kBlack);
	h_t_MC->Draw("HISTsame");
    
    h_t_REC_EEMC_cut->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    for (int i = 1; i <= h_t_REC_EEMC_cut->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_EEMC_cut->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_EEMC_cut->SetBinError(i, newError);
    }
	h_t_REC_EEMC_cut->SetMarkerStyle(20); // method L RECO
	h_t_REC_EEMC_cut->SetMarkerColor(kP8Blue);
	//h_t_REC_EEMC_cut->Draw("PEsame");

    h_t_REC_wRES_cut_pi12->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_wRES_cut_pi12->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_t_REC_wRES_cut_pi12->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_wRES_cut_pi12->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_wRES_cut_pi12->SetBinError(i, newError);
    }
	h_t_REC_wRES_cut_pi12->SetMarkerStyle(30);
	h_t_REC_wRES_cut_pi12->SetMarkerColor(kP8Pink);
    h_t_REC_wRES_cut_pi12->SetLineColor(kP8Pink);
	h_t_REC_wRES_cut_pi12->Draw("PEsame");

    h_DIS_etaOnly_MC->Rebin(rebin_width);
    h_DIS_etaOnly_MC->Scale(sigma_DIS*(1/nEvents_DIS)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_DIS_lumi = nEvents_DIS/sigma_DIS; // nb^-1 
    double ratio_DIS_preTDR = preTDR_lumi/simu_DIS_lumi;
    for (int i = 1; i <= h_DIS_etaOnly_MC->GetNbinsX(); ++i) 
    {
        double binError = h_DIS_etaOnly_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_DIS_preTDR);
        h_DIS_etaOnly_MC->SetBinError(i, newError);
    }
    h_DIS_etaOnly_MC->SetMarkerStyle(2);
    h_DIS_etaOnly_MC->SetMarkerColor(kP8Cyan);
    h_DIS_etaOnly_MC->SetLineColor(kP8Cyan);
    //h_DIS_etaOnly_MC->Draw("PEsame");

    h_DIS_etaOnly_L->Rebin(rebin_width);
    h_DIS_etaOnly_L->Scale(sigma_DIS*(1/nEvents_DIS)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_DIS_etaOnly_L->GetNbinsX(); ++i) 
    {
        double binError = h_DIS_etaOnly_L->GetBinError(i);
        double newError = binError / sqrt(ratio_DIS_preTDR);
        h_DIS_etaOnly_L->SetBinError(i, newError);
    }
    h_DIS_etaOnly_L->SetMarkerStyle(2);
    h_DIS_etaOnly_L->SetMarkerColor(kP8Cyan);
    h_DIS_etaOnly_L->SetLineColor(kP8Cyan);
    //h_DIS_etaOnly_L->Draw("PEsame");

    h_DIS_etaOnly_proj->Rebin(rebin_width);
    h_DIS_etaOnly_proj->Scale(sigma_DIS*(1/nEvents_DIS)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_DIS_etaOnly_proj->GetNbinsX(); ++i) 
    {
        double binError = h_DIS_etaOnly_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_DIS_preTDR);
        h_DIS_etaOnly_proj->SetBinError(i, newError);
    }
    h_DIS_etaOnly_proj->SetMarkerStyle(2);
    h_DIS_etaOnly_proj->SetMarkerColor(kP8Cyan);
    h_DIS_etaOnly_proj->SetLineColor(kP8Cyan);
    h_DIS_etaOnly_proj->Draw("PEsame");

	TLegend *w14_2133 = new TLegend(0.6,0.68,0.72,0.9);
    w14_2133->AddEntry(h_t_MC,"#phi: MC", "L");
    w14_2133->AddEntry(h_t_REC_wRES_cut_pi12,"#phi: #theta_{max}= #pi/12", "P");    
    //w14_2133->AddEntry(h_t_REC_EEMC_cut,"#phi: L", "P");
    //w14_2133->AddEntry(h_DIS_etaOnly_MC,"DIS: MC", "P");
    //w14_2133->AddEntry(h_DIS_etaOnly_L,"DIS: L", "P");
    w14_2133->AddEntry(h_DIS_etaOnly_proj,"DIS: #theta_{max}= #pi/12", "P");
    w14_2133->SetBorderSize(0);   
    w14_2133->SetFillStyle(0);
	w14_2133->SetTextSize(30);
	w14_2133->SetTextFont(45);
    w14_2133->Draw("same");	

	c14_2133->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title14_2133 = new TLatex();
    title14_2133->SetNDC(); 
    title14_2133->SetTextSize(0.05);
    title14_2133->SetTextAlign(22);  
    //title14_2133->DrawLatex(0.5, 0.97, "|t| Distribution DIS Only #eta Vetoes");  

    TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

    // Add labels
    TLatex* ep = new TLatex(0.13, 0.85, "ePIC");
    ep->SetNDC();
    ep->SetTextFont(62);
    ep->SetTextSize(.038);
    ep->SetTextColor(kBlack);
    ep->Draw("same");

    TLatex* r421 = new TLatex(0.18, 0.85, " Simulation 25.10.2, 10x100 GeV");
	r421->SetNDC();
	r421->SetTextSize(30);
	r421->SetTextFont(43);
	r421->SetTextColor(kBlack);
	r421->Draw("same");

	TLatex* label_5 = new TLatex(0.17, 0.75, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->SetTextSize(30);
	label_5->SetTextFont(43);
	label_5->SetTextColor(kBlack);
	label_5->Draw("same");

    TLatex* r4211 = new TLatex(0.15, 0.8, "eAu #rightarrow e'Au'#phi");
	r4211->SetNDC();
	r4211->SetTextSize(30);
	r4211->SetTextFont(43);
	r4211->SetTextColor(kBlack);
	r4211->Draw("same"); 

    TLatex* r421111 = new TLatex(0.15, 0.15, "L_{int} = 10 fb^{-1}/A"); // preTDR lumi
	r421111->SetNDC();
	r421111->SetTextSize(25);
	r421111->SetTextFont(43);
	r421111->SetTextColor(kBlack);
	r421111->Draw("same");

    c14_2133->Print("./figures/plot_t_DIS_etaOnly_preTDR_proj.pdf");

	//TFile* fout_preTDR = new TFile("plot_t_DIS_etaOnly_MC__preTDR.root","RECREATE");
    //c14_2133->Write();
    //h_t_MC->Write();
    //h_t_REC_EEMC_cut->Write();
    //h_t_REC_wRES_cut_pi12->Write();
    //h_DIS_etaOnly_MC->Write();
    //h_DIS_etaOnly_L->Write();
    //h_DIS_etaOnly_proj->Write();
    //fout_preTDR->Close();
}

void plot_DIS_etaVetoes_ES()
{
	TFile* phi_t_file = TFile::Open("phi_sartre_25_10_2_eAu_10x100v5.root","READ");
	TH1D* h_t_MC = (TH1D*) phi_t_file->Get("h_t_MC");
	TH1D* h_t_REC_EEMC_cut = (TH1D*) phi_t_file->Get("h_t_REC_EEMC_cut");
	TH1D* h_t_REC_wRES_cut_pi12 = (TH1D*) phi_t_file->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_phi_sartre_events = (TH1D*) phi_t_file->Get("h_Nevents");

	TFile* DIS_etaOnly = TFile::Open("DIS_beagle_25_10_2_eAu_etaVetoes_10x100v4.root");
    TH1D* h_DIS_etaOnly_MC = (TH1D*) DIS_etaOnly->Get("h_t_MC");
    TH1D* h_DIS_etaOnly_L = (TH1D*) DIS_etaOnly->Get("h_t_REC_EEMC_cut");
    TH1D* h_DIS_etaOnly_proj = (TH1D*) DIS_etaOnly->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_DIS_events = (TH1D*) DIS_etaOnly->Get("h_Nevents");

	double dN_dt = 0.002;
    double sigma_coherent = 459.05; // nb
    double sigma_DIS = 50184.8; // nb
    double nEvents_coherent = h_phi_sartre_events->GetEntries();//6.2e6;
    double nEvents_DIS = h_DIS_events->GetEntries();//1.1e6;
    double early_science_lumi = 8121.83; // nb^-1
    int rebin_width = 2; 

	TCanvas* c14_2133 = new TCanvas("c14_2133","c14_2133",1,1,1200,800);
    gPad->SetLogy(1);
    gStyle->SetOptStat(0);
    h_t_MC->GetYaxis()->SetTitle("dN/d|t| [GeV/c]^{-2}");
	h_t_MC->GetXaxis()->SetTitle("|t| [GeV/c]^{2}");
    h_t_MC->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    double simu_coherent_lumi = nEvents_coherent/sigma_coherent; // nb^-1 
    double ratio_coherent_phi_ES = early_science_lumi/simu_coherent_lumi;
    for (int i = 1; i <= h_t_MC->GetNbinsX(); ++i) 
    {
        double binError = h_t_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_MC->SetBinError(i, newError);
    }
    h_t_MC->SetLineColor(kBlack);
	h_t_MC->Draw("HISTsame");
    
    h_t_REC_EEMC_cut->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    for (int i = 1; i <= h_t_REC_EEMC_cut->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_EEMC_cut->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_REC_EEMC_cut->SetBinError(i, newError);
    }
	h_t_REC_EEMC_cut->SetMarkerStyle(20); // method L RECO
	h_t_REC_EEMC_cut->SetMarkerColor(kP8Blue);
	//h_t_REC_EEMC_cut->Draw("PEsame");

    h_t_REC_wRES_cut_pi12->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_wRES_cut_pi12->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_t_REC_wRES_cut_pi12->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_wRES_cut_pi12->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_REC_wRES_cut_pi12->SetBinError(i, newError);
    }
	h_t_REC_wRES_cut_pi12->SetMarkerStyle(30);
	h_t_REC_wRES_cut_pi12->SetMarkerColor(kP8Pink);
    h_t_REC_wRES_cut_pi12->SetLineColor(kP8Pink);
	h_t_REC_wRES_cut_pi12->Draw("PEsame");

    h_DIS_etaOnly_MC->Rebin(rebin_width);
    h_DIS_etaOnly_MC->Scale(sigma_DIS*(1/nEvents_DIS)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_DIS_lumi = nEvents_DIS/sigma_DIS; // nb^-1 
    double ratio_DIS_ES = early_science_lumi/simu_DIS_lumi;
    for (int i = 1; i <= h_DIS_etaOnly_MC->GetNbinsX(); ++i) 
    {
        double binError = h_DIS_etaOnly_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_DIS_ES);
        h_DIS_etaOnly_MC->SetBinError(i, newError);
    }
    h_DIS_etaOnly_MC->SetMarkerStyle(2);
    h_DIS_etaOnly_MC->SetMarkerColor(kP8Cyan);
    h_DIS_etaOnly_MC->SetLineColor(kP8Cyan);
    //h_DIS_etaOnly_MC->Draw("PEsame");

    h_DIS_etaOnly_L->Rebin(rebin_width);
    h_DIS_etaOnly_L->Scale(sigma_DIS*(1/nEvents_DIS)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_DIS_etaOnly_L->GetNbinsX(); ++i) 
    {
        double binError = h_DIS_etaOnly_L->GetBinError(i);
        double newError = binError / sqrt(ratio_DIS_ES);
        h_DIS_etaOnly_L->SetBinError(i, newError);
    }
    h_DIS_etaOnly_L->SetMarkerStyle(2);
    h_DIS_etaOnly_L->SetMarkerColor(kP8Cyan);
    h_DIS_etaOnly_L->SetLineColor(kP8Cyan);
    //h_DIS_etaOnly_L->Draw("PEsame");

    h_DIS_etaOnly_proj->Rebin(rebin_width);
    h_DIS_etaOnly_proj->Scale(sigma_DIS*(1/nEvents_DIS)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_DIS_etaOnly_proj->GetNbinsX(); ++i) 
    {
        double binError = h_DIS_etaOnly_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_DIS_ES);
        h_DIS_etaOnly_proj->SetBinError(i, newError);
    }
    h_DIS_etaOnly_proj->SetMarkerStyle(2);
    h_DIS_etaOnly_proj->SetMarkerColor(kP8Cyan);
    h_DIS_etaOnly_proj->SetLineColor(kP8Cyan);
    h_DIS_etaOnly_proj->Draw("PEsame");

	TLegend *w14_2133 = new TLegend(0.6,0.68,0.72,0.9);
    w14_2133->AddEntry(h_t_MC,"#phi: MC", "L");
    w14_2133->AddEntry(h_t_REC_wRES_cut_pi12,"#phi: #theta_{max}= #pi/12", "P");    
    //w14_2133->AddEntry(h_t_REC_EEMC_cut,"#phi: L", "P");
    //w14_2133->AddEntry(h_DIS_etaOnly_MC,"DIS: MC", "P");
    //w14_2133->AddEntry(h_DIS_etaOnly_L,"DIS: L", "P");
    w14_2133->AddEntry(h_DIS_etaOnly_proj,"DIS: #theta_{max}= #pi/12", "P");
    w14_2133->SetBorderSize(0);   
    w14_2133->SetFillStyle(0);
	w14_2133->SetTextSize(30);
	w14_2133->SetTextFont(45);
    w14_2133->Draw("same");	

	c14_2133->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title14_2133 = new TLatex();
    title14_2133->SetNDC(); 
    title14_2133->SetTextSize(0.05);
    title14_2133->SetTextAlign(22);  
    //title14_2133->DrawLatex(0.5, 0.97, "|t| Distribution DIS Only #eta Vetoes");  

    TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

    // Add labels
    TLatex* ep = new TLatex(0.13, 0.85, "ePIC");
    ep->SetNDC();
    ep->SetTextFont(62);
    ep->SetTextSize(.038);
    ep->SetTextColor(kBlack);
    ep->Draw("same");

    TLatex* r421 = new TLatex(0.18, 0.85, " Simulation 25.10.2, 10x100 GeV");
	r421->SetNDC();
	r421->SetTextSize(30);
	r421->SetTextFont(43);
	r421->SetTextColor(kBlack);
	r421->Draw("same");

	TLatex* label_5 = new TLatex(0.17, 0.75, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->SetTextSize(30);
	label_5->SetTextFont(43);
	label_5->SetTextColor(kBlack);
	label_5->Draw("same");

    TLatex* r4211 = new TLatex(0.15, 0.8, "eAu #rightarrow e'Au'#phi");
	r4211->SetNDC();
	r4211->SetTextSize(30);
	r4211->SetTextFont(43);
	r4211->SetTextColor(kBlack);
	r4211->Draw("same"); 

    TLatex* r42111 = new TLatex(0.15, 0.15, "L_{int} = 1.6 fb^{-1}/A"); // ES lumi
	r42111->SetNDC();
	r42111->SetTextSize(30);
	r42111->SetTextFont(43);
	r42111->SetTextColor(kBlack);
	r42111->Draw("same");

    c14_2133->Print("./figures/plot_t_DIS_etaOnly_ES_proj.pdf");

	//TFile* fout_ES = new TFile("plot_t_DIS_etaOnly_MC_ES.root","RECREATE");
    //c14_2133->Write();
    //h_t_MC->Write();
    //h_t_REC_EEMC_cut->Write();
    //h_t_REC_wRES_cut_pi12->Write();
    //h_DIS_etaOnly_MC->Write();
    //h_DIS_etaOnly_L->Write();
    //h_DIS_etaOnly_proj->Write();
    //fout_ES->Close();
}

void plot_DIS_omdVetoes_preTDR()
{
    TFile* phi_t_file = TFile::Open("phi_sartre_25_10_2_eAu_10x100v5.root","READ");
	TH1D* h_t_MC = (TH1D*) phi_t_file->Get("h_t_MC");
	TH1D* h_t_REC_EEMC_cut = (TH1D*) phi_t_file->Get("h_t_REC_EEMC_cut");
	TH1D* h_t_REC_wRES_cut_pi12 = (TH1D*) phi_t_file->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_phi_sartre_events = (TH1D*) phi_t_file->Get("h_Nevents");

    TFile* DIS_OMDonly = TFile::Open("DIS_beagle_25_10_2_eAu_omdVetoes_10x100v4.root");
    TH1D* h_DIS_OMDonly_MC = (TH1D*) DIS_OMDonly->Get("h_t_MC");
    TH1D* h_DIS_OMDonly_L = (TH1D*) DIS_OMDonly->Get("h_t_REC_EEMC_cut");
    TH1D* h_DIS_OMDonly_proj = (TH1D*) DIS_OMDonly->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_DIS_events = (TH1D*) DIS_OMDonly->Get("h_Nevents");

    double dN_dt = 0.002;
    double sigma_coherent = 459.05; // nb
    double sigma_DIS = 50184.8; // nb
    double nEvents_coherent = h_phi_sartre_events->GetEntries();//6.2e6;
    double nEvents_DIS = h_DIS_events->GetEntries();//1.1e6;
    double preTDR_lumi = 50761.4; // nb^-1
    int rebin_width = 2; 

    TCanvas* c14_2134 = new TCanvas("c14_2134","c14_2134",1,1,1200,800);
    gPad->SetLogy(1);
    gStyle->SetOptStat(0);
    h_t_MC->GetYaxis()->SetTitle("dN/d|t| [GeV/c]^{-2}");
	h_t_MC->GetXaxis()->SetTitle("|t| [GeV/c]^{2}");
    h_t_MC->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    double simu_coherent_lumi = nEvents_coherent/sigma_coherent; // nb^-1 
    double ratio_coherent_phi_preTDR = preTDR_lumi/simu_coherent_lumi;
    for (int i = 1; i <= h_t_MC->GetNbinsX(); ++i) 
    {
        double binError = h_t_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_MC->SetBinError(i, newError);
    }
    h_t_MC->SetLineColor(kBlack);
	h_t_MC->Draw("HISTsame");

    h_t_REC_EEMC_cut->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    for (int i = 1; i <= h_t_REC_EEMC_cut->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_EEMC_cut->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_EEMC_cut->SetBinError(i, newError);
    }
	h_t_REC_EEMC_cut->SetMarkerStyle(20); // method L RECO
	h_t_REC_EEMC_cut->SetMarkerColor(kP8Blue);
    h_t_REC_EEMC_cut->SetLineColor(kP8Blue);
	//h_t_REC_EEMC_cut->Draw("PEsame");

    h_t_REC_wRES_cut_pi12->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_wRES_cut_pi12->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_t_REC_wRES_cut_pi12->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_wRES_cut_pi12->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_wRES_cut_pi12->SetBinError(i, newError);
    }
	h_t_REC_wRES_cut_pi12->SetMarkerStyle(30);
	h_t_REC_wRES_cut_pi12->SetMarkerColor(kP8Pink);
    h_t_REC_wRES_cut_pi12->SetLineColor(kP8Pink);
	h_t_REC_wRES_cut_pi12->Draw("PEsame");

    h_DIS_OMDonly_MC->Rebin(rebin_width);
    h_DIS_OMDonly_MC->Scale(sigma_DIS*(1/nEvents_DIS)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_DIS_lumi = nEvents_DIS/sigma_DIS; // nb^-1 
    double ratio_DIS_preTDR = preTDR_lumi/simu_DIS_lumi;
    for (int i = 1; i <= h_DIS_OMDonly_MC->GetNbinsX(); ++i) 
    {
        double binError = h_DIS_OMDonly_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_DIS_preTDR);
        h_DIS_OMDonly_MC->SetBinError(i, newError);
    }
    h_DIS_OMDonly_MC->SetMarkerStyle(2);
    h_DIS_OMDonly_MC->SetMarkerColor(kP8Cyan);
    h_DIS_OMDonly_MC->SetLineColor(kP8Cyan);
    //h_DIS_OMDonly_MC->Draw("PEsame");

    h_DIS_OMDonly_L->Rebin(rebin_width);
    h_DIS_OMDonly_L->Scale(sigma_DIS*(1/nEvents_DIS)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_DIS_OMDonly_L->GetNbinsX(); ++i) 
    {
        double binError = h_DIS_OMDonly_L->GetBinError(i);
        double newError = binError / sqrt(ratio_DIS_preTDR);
        h_DIS_OMDonly_L->SetBinError(i, newError);
    }
    h_DIS_OMDonly_L->SetMarkerStyle(2);
    h_DIS_OMDonly_L->SetMarkerColor(kP8Cyan);
    h_DIS_OMDonly_L->SetLineColor(kP8Cyan);
    //h_DIS_OMDonly_L->Draw("PEsame");

    h_DIS_OMDonly_proj->Rebin(rebin_width);
    h_DIS_OMDonly_proj->Scale(sigma_DIS*(1/nEvents_DIS)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_DIS_OMDonly_proj->GetNbinsX(); ++i) 
    {
        double binError = h_DIS_OMDonly_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_DIS_preTDR);
        h_DIS_OMDonly_proj->SetBinError(i, newError);
    }
    h_DIS_OMDonly_proj->SetMarkerStyle(2);
    h_DIS_OMDonly_proj->SetMarkerColor(kP8Cyan);
    h_DIS_OMDonly_proj->SetLineColor(kP8Cyan);
    h_DIS_OMDonly_proj->Draw("PEsame");

    TLegend *w14_2134 = new TLegend(0.6,0.68,0.72,0.9);
	w14_2134->AddEntry(h_t_MC,"#phi: MC", "L");
    w14_2134->AddEntry(h_t_REC_wRES_cut_pi12,"#phi: #theta_{max}= #pi/12", "P");    
    //w14_2134->AddEntry(h_t_REC_EEMC_cut,"#phi: L", "P");
    //w14_2134->AddEntry(h_DIS_OMDonly_MC,"DIS: MC", "P");
    //w14_2134->AddEntry(h_DIS_OMDonly_L,"DIS: L", "P");
    w14_2134->AddEntry(h_DIS_OMDonly_proj,"DIS: #theta_{max}= #pi/12", "P");
    w14_2134->SetBorderSize(0);   
    w14_2134->SetFillStyle(0);
	w14_2134->SetTextSize(30);
	w14_2134->SetTextFont(45);
    w14_2134->Draw("same");	

    c14_2134->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title14_2134 = new TLatex();
    title14_2134->SetNDC(); 
    title14_2134->SetTextSize(0.05);
    title14_2134->SetTextAlign(22);  
    //title14_2134->DrawLatex(0.5, 0.97, "|t| Distribution DIS OMD Vetoes Only");  

    TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

    // Add labels
    TLatex* ep = new TLatex(0.13, 0.85, "ePIC");
    ep->SetNDC();
    ep->SetTextFont(62);
    ep->SetTextSize(.038);
    ep->SetTextColor(kBlack);
    ep->Draw("same");

    TLatex* r421 = new TLatex(0.18, 0.85, " Simulation 25.10.2, 10x100 GeV");
	r421->SetNDC();
	r421->SetTextSize(30);
	r421->SetTextFont(43);
	r421->SetTextColor(kBlack);
	r421->Draw("same");

	TLatex* label_5 = new TLatex(0.17, 0.75, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->SetTextSize(30);
	label_5->SetTextFont(43);
	label_5->SetTextColor(kBlack);
	label_5->Draw("same");

    TLatex* r4211 = new TLatex(0.15, 0.8, "eAu #rightarrow e'Au'#phi");
	r4211->SetNDC();
	r4211->SetTextSize(30);
	r4211->SetTextFont(43);
	r4211->SetTextColor(kBlack);
	r4211->Draw("same"); 

    TLatex* r421111 = new TLatex(0.15, 0.15, "L_{int} = 10 fb^{-1}/A"); // preTDR lumi
	r421111->SetNDC();
	r421111->SetTextSize(25);
	r421111->SetTextFont(43);
	r421111->SetTextColor(kBlack);
	r421111->Draw("same");

    c14_2134->Print("./figures/plot_t_DIS_OMDonly_preTDR_proj.pdf");
}

void plot_DIS_omdVetoes_ES()
{
    TFile* phi_t_file = TFile::Open("phi_sartre_25_10_2_eAu_10x100v5.root","READ");
	TH1D* h_t_MC = (TH1D*) phi_t_file->Get("h_t_MC");
	TH1D* h_t_REC_EEMC_cut = (TH1D*) phi_t_file->Get("h_t_REC_EEMC_cut");
	TH1D* h_t_REC_wRES_cut_pi12 = (TH1D*) phi_t_file->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_phi_sartre_events = (TH1D*) phi_t_file->Get("h_Nevents");

    TFile* DIS_OMDonly = TFile::Open("DIS_beagle_25_10_2_eAu_omdVetoes_10x100v4.root");
    TH1D* h_DIS_OMDonly_MC = (TH1D*) DIS_OMDonly->Get("h_t_MC");
    TH1D* h_DIS_OMDonly_L = (TH1D*) DIS_OMDonly->Get("h_t_REC_EEMC_cut");
    TH1D* h_DIS_OMDonly_proj = (TH1D*) DIS_OMDonly->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_DIS_events = (TH1D*) DIS_OMDonly->Get("h_Nevents");

    double dN_dt = 0.002;
    double sigma_coherent = 459.05; // nb
    double sigma_DIS = 50184.8; // nb
    double nEvents_coherent = h_phi_sartre_events->GetEntries();//6.2e6;
    double nEvents_DIS = h_DIS_events->GetEntries();//1.1e6;
    double early_science_lumi = 8121.83; // nb^-1
    int rebin_width = 2; 

    TCanvas* c14_2134 = new TCanvas("c14_2134","c14_2134",1,1,1200,800);
    gPad->SetLogy(1);
    gStyle->SetOptStat(0);
    h_t_MC->GetYaxis()->SetTitle("dN/d|t| [GeV/c]^{-2}");
	h_t_MC->GetXaxis()->SetTitle("|t| [GeV/c]^{2}");
    h_t_MC->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    double simu_coherent_lumi = nEvents_coherent/sigma_coherent; // nb^-1 
    double ratio_coherent_phi_ES = early_science_lumi/simu_coherent_lumi;
    for (int i = 1; i <= h_t_MC->GetNbinsX(); ++i) 
    {
        double binError = h_t_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_MC->SetBinError(i, newError);
    }
    h_t_MC->SetLineColor(kBlack);
	h_t_MC->Draw("HISTsame");

    h_t_REC_EEMC_cut->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    for (int i = 1; i <= h_t_REC_EEMC_cut->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_EEMC_cut->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_REC_EEMC_cut->SetBinError(i, newError);
    }
	h_t_REC_EEMC_cut->SetMarkerStyle(20); // method L RECO
	h_t_REC_EEMC_cut->SetMarkerColor(kP8Blue);
    h_t_REC_EEMC_cut->SetLineColor(kP8Blue);
	//h_t_REC_EEMC_cut->Draw("PEsame");

    h_t_REC_wRES_cut_pi12->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_wRES_cut_pi12->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_t_REC_wRES_cut_pi12->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_wRES_cut_pi12->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_REC_wRES_cut_pi12->SetBinError(i, newError);
    }
	h_t_REC_wRES_cut_pi12->SetMarkerStyle(30);
	h_t_REC_wRES_cut_pi12->SetMarkerColor(kP8Pink);
    h_t_REC_wRES_cut_pi12->SetLineColor(kP8Pink);
	h_t_REC_wRES_cut_pi12->Draw("PEsame");

    h_DIS_OMDonly_MC->Rebin(rebin_width);
    h_DIS_OMDonly_MC->Scale(sigma_DIS*(1/nEvents_DIS)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_DIS_lumi = nEvents_DIS/sigma_DIS; // nb^-1 
    double ratio_DIS_ES = early_science_lumi/simu_DIS_lumi;
    for (int i = 1; i <= h_DIS_OMDonly_MC->GetNbinsX(); ++i) 
    {
        double binError = h_DIS_OMDonly_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_DIS_ES);
        h_DIS_OMDonly_MC->SetBinError(i, newError);
    }
    h_DIS_OMDonly_MC->SetMarkerStyle(2);
    h_DIS_OMDonly_MC->SetMarkerColor(kP8Cyan);
    h_DIS_OMDonly_MC->SetLineColor(kP8Cyan);
    //h_DIS_OMDonly_MC->Draw("PEsame");

    h_DIS_OMDonly_L->Rebin(rebin_width);
    h_DIS_OMDonly_L->Scale(sigma_DIS*(1/nEvents_DIS)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_DIS_OMDonly_L->GetNbinsX(); ++i) 
    {
        double binError = h_DIS_OMDonly_L->GetBinError(i);
        double newError = binError / sqrt(ratio_DIS_ES);
        h_DIS_OMDonly_L->SetBinError(i, newError);
    }
    h_DIS_OMDonly_L->SetMarkerStyle(2);
    h_DIS_OMDonly_L->SetMarkerColor(kP8Cyan);
    h_DIS_OMDonly_L->SetLineColor(kP8Cyan);
    //h_DIS_OMDonly_L->Draw("PEsame");

    h_DIS_OMDonly_proj->Rebin(rebin_width);
    h_DIS_OMDonly_proj->Scale(sigma_DIS*(1/nEvents_DIS)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_DIS_OMDonly_proj->GetNbinsX(); ++i) 
    {
        double binError = h_DIS_OMDonly_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_DIS_ES);
        h_DIS_OMDonly_proj->SetBinError(i, newError);
    }
    h_DIS_OMDonly_proj->SetMarkerStyle(2);
    h_DIS_OMDonly_proj->SetMarkerColor(kP8Cyan);
    h_DIS_OMDonly_proj->SetLineColor(kP8Cyan);
    h_DIS_OMDonly_proj->Draw("PEsame");

    TLegend *w14_2134 = new TLegend(0.6,0.68,0.72,0.9);
	w14_2134->AddEntry(h_t_MC,"#phi: MC", "L");
    w14_2134->AddEntry(h_t_REC_wRES_cut_pi12,"#phi: #theta_{max}= #pi/12", "P");    
    //w14_2134->AddEntry(h_t_REC_EEMC_cut,"#phi: L", "P");
    //w14_2134->AddEntry(h_DIS_OMDonly_MC,"DIS: MC", "P");
    //w14_2134->AddEntry(h_DIS_OMDonly_L,"DIS: L", "P");
    w14_2134->AddEntry(h_DIS_OMDonly_proj,"DIS: #theta_{max}= #pi/12", "P");
    w14_2134->SetBorderSize(0);   
    w14_2134->SetFillStyle(0);
	w14_2134->SetTextSize(30);
	w14_2134->SetTextFont(45);
    w14_2134->Draw("same");	

    c14_2134->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title14_2134 = new TLatex();
    title14_2134->SetNDC(); 
    title14_2134->SetTextSize(0.05);
    title14_2134->SetTextAlign(22);  
    //title14_2134->DrawLatex(0.5, 0.97, "|t| Distribution DIS OMD Vetoes Only");  

    TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

    // Add labels
    TLatex* ep = new TLatex(0.13, 0.85, "ePIC");
    ep->SetNDC();
    ep->SetTextFont(62);
    ep->SetTextSize(.038);
    ep->SetTextColor(kBlack);
    ep->Draw("same");

    TLatex* r421 = new TLatex(0.18, 0.85, " Simulation 25.10.2, 10x100 GeV");
	r421->SetNDC();
	r421->SetTextSize(30);
	r421->SetTextFont(43);
	r421->SetTextColor(kBlack);
	r421->Draw("same");

	TLatex* label_5 = new TLatex(0.17, 0.75, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->SetTextSize(30);
	label_5->SetTextFont(43);
	label_5->SetTextColor(kBlack);
	label_5->Draw("same");

    TLatex* r4211 = new TLatex(0.15, 0.8, "eAu #rightarrow e'Au'#phi");
	r4211->SetNDC();
	r4211->SetTextSize(30);
	r4211->SetTextFont(43);
	r4211->SetTextColor(kBlack);
	r4211->Draw("same"); 

    TLatex* r42111 = new TLatex(0.15, 0.15, "L_{int} = 1.6 fb^{-1}/A"); // ES lumi
	r42111->SetNDC();
	r42111->SetTextSize(30);
	r42111->SetTextFont(43);
	r42111->SetTextColor(kBlack);
	r42111->Draw("same");

    c14_2134->Print("./figures/plot_t_DIS_OMDonly_ES_proj.pdf");
}

void plot_DIS_rpVetoes_preTDR()
{
    TFile* phi_t_file = TFile::Open("phi_sartre_25_10_2_eAu_10x100v5.root","READ");
	TH1D* h_t_MC = (TH1D*) phi_t_file->Get("h_t_MC");
	TH1D* h_t_REC_EEMC_cut = (TH1D*) phi_t_file->Get("h_t_REC_EEMC_cut");
	TH1D* h_t_REC_wRES_cut_pi12 = (TH1D*) phi_t_file->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_phi_sartre_events = (TH1D*) phi_t_file->Get("h_Nevents");

    TFile* DIS_RPonly = TFile::Open("DIS_beagle_25_10_2_eAu_rpVetoes_10x100v4.root");
    TH1D* h_DIS_RPonly_MC = (TH1D*) DIS_RPonly->Get("h_t_MC");
    TH1D* h_DIS_RPonly_L = (TH1D*) DIS_RPonly->Get("h_t_REC_EEMC_cut");
    TH1D* h_DIS_RPonly_proj = (TH1D*) DIS_RPonly->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_DIS_events = (TH1D*) DIS_RPonly->Get("h_Nevents");

    double dN_dt = 0.002;
    double sigma_coherent = 459.05; // nb
    double sigma_DIS = 50184.8; // nb
    double nEvents_coherent = h_phi_sartre_events->GetEntries();//6.2e6;
    double nEvents_DIS = h_DIS_events->GetEntries();//1.1e6;
    double preTDR_lumi = 50761.4; // nb^-1
    int rebin_width = 2; 

    TCanvas* c14_2135 = new TCanvas("c14_2135","c14_2135",1,1,1200,800);
    gPad->SetLogy(1);
    gStyle->SetOptStat(0);
    h_t_MC->GetYaxis()->SetTitle("dN/d|t| [GeV/c]^{-2}");
	h_t_MC->GetXaxis()->SetTitle("|t| [GeV/c]^{2}");
    h_t_MC->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    double simu_coherent_lumi = nEvents_coherent/sigma_coherent; // nb^-1 
    double ratio_coherent_phi_preTDR = preTDR_lumi/simu_coherent_lumi;
    for (int i = 1; i <= h_t_MC->GetNbinsX(); ++i) 
    {
        double binError = h_t_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_MC->SetBinError(i, newError);
    }
    h_t_MC->SetLineColor(kBlack);
	h_t_MC->Draw("HISTsame");

    h_t_REC_EEMC_cut->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    for (int i = 1; i <= h_t_REC_EEMC_cut->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_EEMC_cut->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_EEMC_cut->SetBinError(i, newError);
    }
	h_t_REC_EEMC_cut->SetMarkerStyle(20); // method L RECO
	h_t_REC_EEMC_cut->SetMarkerColor(kP8Blue);
	//h_t_REC_EEMC_cut->Draw("PEsame");

    h_t_REC_wRES_cut_pi12->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_wRES_cut_pi12->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_t_REC_wRES_cut_pi12->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_wRES_cut_pi12->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_wRES_cut_pi12->SetBinError(i, newError);
    }
	h_t_REC_wRES_cut_pi12->SetMarkerStyle(30);
	h_t_REC_wRES_cut_pi12->SetMarkerColor(kP8Pink);
    h_t_REC_wRES_cut_pi12->SetLineColor(kP8Pink);
	h_t_REC_wRES_cut_pi12->Draw("PEsame");

    h_DIS_RPonly_MC->Rebin(rebin_width);
    h_DIS_RPonly_MC->Scale(sigma_DIS*(1/nEvents_DIS)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_DIS_lumi = nEvents_DIS/sigma_DIS; // nb^-1 
    double ratio_DIS_preTDR = preTDR_lumi/simu_DIS_lumi;
    for (int i = 1; i <= h_DIS_RPonly_MC->GetNbinsX(); ++i) 
    {
        double binError = h_DIS_RPonly_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_DIS_preTDR);
        h_DIS_RPonly_MC->SetBinError(i, newError);
    }
    h_DIS_RPonly_MC->SetMarkerStyle(2);
    h_DIS_RPonly_MC->SetMarkerColor(kP8Cyan);
    h_DIS_RPonly_MC->SetLineColor(kP8Cyan);
    //h_DIS_RPonly_MC->Draw("PEsame");

    h_DIS_RPonly_L->Rebin(rebin_width);
    h_DIS_RPonly_L->Scale(sigma_DIS*(1/nEvents_DIS)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_DIS_RPonly_L->GetNbinsX(); ++i) 
    {
        double binError = h_DIS_RPonly_L->GetBinError(i);
        double newError = binError / sqrt(ratio_DIS_preTDR);
        h_DIS_RPonly_L->SetBinError(i, newError);
    }
    h_DIS_RPonly_L->SetMarkerStyle(2);
    h_DIS_RPonly_L->SetMarkerColor(kP8Cyan);
    h_DIS_RPonly_L->SetLineColor(kP8Cyan);
    //h_DIS_RPonly_L->Draw("PEsame");

    h_DIS_RPonly_proj->Rebin(rebin_width);
    h_DIS_RPonly_proj->Scale(sigma_DIS*(1/nEvents_DIS)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_DIS_RPonly_proj->GetNbinsX(); ++i) 
    {
        double binError = h_DIS_RPonly_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_DIS_preTDR);
        h_DIS_RPonly_proj->SetBinError(i, newError);
    }
    h_DIS_RPonly_proj->SetMarkerStyle(2);
    h_DIS_RPonly_proj->SetMarkerColor(kP8Cyan);
    h_DIS_RPonly_proj->SetLineColor(kP8Cyan);
    h_DIS_RPonly_proj->Draw("PEsame");
    
    TLegend *w14_2135 = new TLegend(0.62,0.68,0.72,0.9);
	w14_2135->AddEntry(h_t_MC,"#phi: MC", "L");
    w14_2135->AddEntry(h_t_REC_wRES_cut_pi12,"#phi: #theta_{max}= #pi/12", "P");    
    //w14_2135->AddEntry(h_t_REC_EEMC_cut,"#phi: L", "P");
    //w14_2135->AddEntry(h_DIS_RPonly_MC,"DIS: MC", "P");
    //w14_2135->AddEntry(h_DIS_RPonly_L,"DIS: L", "P");
    w14_2135->AddEntry(h_DIS_RPonly_proj,"DIS: #theta_{max}= #pi/12", "P");
    w14_2135->SetBorderSize(0);   
    w14_2135->SetFillStyle(0);
	w14_2135->SetTextSize(30);
	w14_2135->SetTextFont(45);
    w14_2135->Draw("same");	

    c14_2135->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title14_2135 = new TLatex();
    title14_2135->SetNDC(); 
    title14_2135->SetTextSize(0.05);
    title14_2135->SetTextAlign(22);  
    //title14_2135->DrawLatex(0.5, 0.97, "|t| Distribution DIS RP Vetoes Only");  

    TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

    // Add labels
    TLatex* ep = new TLatex(0.13, 0.85, "ePIC");
    ep->SetNDC();
    ep->SetTextFont(62);
    ep->SetTextSize(.038);
    ep->SetTextColor(kBlack);
    ep->Draw("same");

    TLatex* r421 = new TLatex(0.18, 0.85, " Simulation 25.10.2, 10x100 GeV");
	r421->SetNDC();
	r421->SetTextSize(30);
	r421->SetTextFont(43);
	r421->SetTextColor(kBlack);
	r421->Draw("same");

	TLatex* label_5 = new TLatex(0.17, 0.75, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->SetTextSize(30);
	label_5->SetTextFont(43);
	label_5->SetTextColor(kBlack);
	label_5->Draw("same");

    TLatex* r4211 = new TLatex(0.15, 0.8, "eAu #rightarrow e'Au'#phi");
	r4211->SetNDC();
	r4211->SetTextSize(30);
	r4211->SetTextFont(43);
	r4211->SetTextColor(kBlack);
	r4211->Draw("same"); 

    TLatex* r421111 = new TLatex(0.15, 0.15, "L_{int} = 10 fb^{-1}/A"); // preTDR lumi
	r421111->SetNDC();
	r421111->SetTextSize(25);
	r421111->SetTextFont(43);
	r421111->SetTextColor(kBlack);
	r421111->Draw("same");

    c14_2135->Print("./figures/plot_t_DIS_RPonly_preTDR_proj.pdf");
}

void plot_DIS_rpVetoes_ES()
{
    TFile* phi_t_file = TFile::Open("phi_sartre_25_10_2_eAu_10x100v5.root","READ");
	TH1D* h_t_MC = (TH1D*) phi_t_file->Get("h_t_MC");
	TH1D* h_t_REC_EEMC_cut = (TH1D*) phi_t_file->Get("h_t_REC_EEMC_cut");
	TH1D* h_t_REC_wRES_cut_pi12 = (TH1D*) phi_t_file->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_phi_sartre_events = (TH1D*) phi_t_file->Get("h_Nevents");

    TFile* DIS_RPonly = TFile::Open("DIS_beagle_25_10_2_eAu_rpVetoes_10x100v4.root");
    TH1D* h_DIS_RPonly_MC = (TH1D*) DIS_RPonly->Get("h_t_MC");
    TH1D* h_DIS_RPonly_L = (TH1D*) DIS_RPonly->Get("h_t_REC_EEMC_cut");
    TH1D* h_DIS_RPonly_proj = (TH1D*) DIS_RPonly->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_DIS_events = (TH1D*) DIS_RPonly->Get("h_Nevents");

    double dN_dt = 0.002;
    double sigma_coherent = 459.05; // nb
    double sigma_DIS = 50184.8; // nb
    double nEvents_coherent = h_phi_sartre_events->GetEntries();//6.2e6;
    double nEvents_DIS = h_DIS_events->GetEntries();//1.1e6;
    double early_science_lumi = 8121.83; // nb^-1
    int rebin_width = 2; 

    TCanvas* c14_2135 = new TCanvas("c14_2135","c14_2135",1,1,1200,800);
    gPad->SetLogy(1);
    gStyle->SetOptStat(0);
    h_t_MC->GetYaxis()->SetTitle("dN/d|t| [GeV/c]^{-2}");
	h_t_MC->GetXaxis()->SetTitle("|t| [GeV/c]^{2}");
    h_t_MC->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    double simu_coherent_lumi = nEvents_coherent/sigma_coherent; // nb^-1 
    double ratio_coherent_phi_ES = early_science_lumi/simu_coherent_lumi;
    for (int i = 1; i <= h_t_MC->GetNbinsX(); ++i) 
    {
        double binError = h_t_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_MC->SetBinError(i, newError);
    }
    h_t_MC->SetLineColor(kBlack);
	h_t_MC->Draw("HISTsame");

    h_t_REC_EEMC_cut->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    for (int i = 1; i <= h_t_REC_EEMC_cut->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_EEMC_cut->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_REC_EEMC_cut->SetBinError(i, newError);
    }
	h_t_REC_EEMC_cut->SetMarkerStyle(20); // method L RECO
	h_t_REC_EEMC_cut->SetMarkerColor(kP8Blue);
	//h_t_REC_EEMC_cut->Draw("PEsame");

    h_t_REC_wRES_cut_pi12->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_wRES_cut_pi12->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_t_REC_wRES_cut_pi12->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_wRES_cut_pi12->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_REC_wRES_cut_pi12->SetBinError(i, newError);
    }
	h_t_REC_wRES_cut_pi12->SetMarkerStyle(30);
	h_t_REC_wRES_cut_pi12->SetMarkerColor(kP8Pink);
    h_t_REC_wRES_cut_pi12->SetLineColor(kP8Pink);
	h_t_REC_wRES_cut_pi12->Draw("PEsame");

    h_DIS_RPonly_MC->Rebin(rebin_width);
    h_DIS_RPonly_MC->Scale(sigma_DIS*(1/nEvents_DIS)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_DIS_lumi = nEvents_DIS/sigma_DIS; // nb^-1 
    double ratio_DIS_ES = early_science_lumi/simu_DIS_lumi;
    for (int i = 1; i <= h_DIS_RPonly_MC->GetNbinsX(); ++i) 
    {
        double binError = h_DIS_RPonly_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_DIS_ES);
        h_DIS_RPonly_MC->SetBinError(i, newError);
    }
    h_DIS_RPonly_MC->SetMarkerStyle(2);
    h_DIS_RPonly_MC->SetMarkerColor(kP8Cyan);
    h_DIS_RPonly_MC->SetLineColor(kP8Cyan);
    //h_DIS_RPonly_MC->Draw("PEsame");

    h_DIS_RPonly_L->Rebin(rebin_width);
    h_DIS_RPonly_L->Scale(sigma_DIS*(1/nEvents_DIS)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_DIS_RPonly_L->GetNbinsX(); ++i) 
    {
        double binError = h_DIS_RPonly_L->GetBinError(i);
        double newError = binError / sqrt(ratio_DIS_ES);
        h_DIS_RPonly_L->SetBinError(i, newError);
    }
    h_DIS_RPonly_L->SetMarkerStyle(2);
    h_DIS_RPonly_L->SetMarkerColor(kP8Cyan);
    h_DIS_RPonly_L->SetLineColor(kP8Cyan);
    //h_DIS_RPonly_L->Draw("PEsame");

    h_DIS_RPonly_proj->Rebin(rebin_width);
    h_DIS_RPonly_proj->Scale(sigma_DIS*(1/nEvents_DIS)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_DIS_RPonly_proj->GetNbinsX(); ++i) 
    {
        double binError = h_DIS_RPonly_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_DIS_ES);
        h_DIS_RPonly_proj->SetBinError(i, newError);
    }
    h_DIS_RPonly_proj->SetMarkerStyle(2);
    h_DIS_RPonly_proj->SetMarkerColor(kP8Cyan);
    h_DIS_RPonly_proj->SetLineColor(kP8Cyan);
    h_DIS_RPonly_proj->Draw("PEsame");
    
    TLegend *w14_2135 = new TLegend(0.62,0.68,0.72,0.9);
	w14_2135->AddEntry(h_t_MC,"#phi: MC", "L");
    w14_2135->AddEntry(h_t_REC_wRES_cut_pi12,"#phi: #theta_{max}= #pi/12", "P");    
    //w14_2135->AddEntry(h_t_REC_EEMC_cut,"#phi: L", "P");
    //w14_2135->AddEntry(h_DIS_RPonly_MC,"DIS: MC", "P");
    //w14_2135->AddEntry(h_DIS_RPonly_L,"DIS: L", "P");
    w14_2135->AddEntry(h_DIS_RPonly_proj,"DIS: #theta_{max}= #pi/12", "P");
    w14_2135->SetBorderSize(0);   
    w14_2135->SetFillStyle(0);
	w14_2135->SetTextSize(30);
	w14_2135->SetTextFont(45);
    w14_2135->Draw("same");	

    c14_2135->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title14_2135 = new TLatex();
    title14_2135->SetNDC(); 
    title14_2135->SetTextSize(0.05);
    title14_2135->SetTextAlign(22);  
    //title14_2135->DrawLatex(0.5, 0.97, "|t| Distribution DIS RP Vetoes Only");  

    TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

    // Add labels
    TLatex* ep = new TLatex(0.13, 0.85, "ePIC");
    ep->SetNDC();
    ep->SetTextFont(62);
    ep->SetTextSize(.038);
    ep->SetTextColor(kBlack);
    ep->Draw("same");

    TLatex* r421 = new TLatex(0.18, 0.85, " Simulation 25.10.2, 10x100 GeV");
	r421->SetNDC();
	r421->SetTextSize(30);
	r421->SetTextFont(43);
	r421->SetTextColor(kBlack);
	r421->Draw("same");

	TLatex* label_5 = new TLatex(0.17, 0.75, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->SetTextSize(30);
	label_5->SetTextFont(43);
	label_5->SetTextColor(kBlack);
	label_5->Draw("same");

    TLatex* r4211 = new TLatex(0.15, 0.8, "eAu #rightarrow e'Au'#phi");
	r4211->SetNDC();
	r4211->SetTextSize(30);
	r4211->SetTextFont(43);
	r4211->SetTextColor(kBlack);
	r4211->Draw("same"); 

    TLatex* r42111 = new TLatex(0.15, 0.15, "L_{int} = 1.6 fb^{-1}/A"); // ES lumi
	r42111->SetNDC();
	r42111->SetTextSize(30);
	r42111->SetTextFont(43);
	r42111->SetTextColor(kBlack);
	r42111->Draw("same");

    c14_2135->Print("./figures/plot_t_DIS_RPonly_ES_proj.pdf");
}

void plot_DIS_together_preTDR()
{
    TFile* phi_t_file = TFile::Open("phi_sartre_25_10_2_eAu_10x100v5.root","READ");
	TH1D* h_t_MC = (TH1D*) phi_t_file->Get("h_t_MC");
	TH1D* h_t_REC_EEMC_cut = (TH1D*) phi_t_file->Get("h_t_REC_EEMC_cut");
	TH1D* h_t_REC_wRES_cut_pi12 = (TH1D*) phi_t_file->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_phi_sartre_events = (TH1D*) phi_t_file->Get("h_Nevents");

	TFile* DIS_noVetoes = TFile::Open("DIS_beagle_25_10_2_eAu_noVetoes_10x100v4.root","READ");
    TH1D* h_DIS_noVetoes_MC = (TH1D*) DIS_noVetoes->Get("h_t_MC");
    TH1D* h_DIS_noVetoes_L = (TH1D*) DIS_noVetoes->Get("h_t_REC_EEMC_cut");
    TH1D* h_DIS_noVetoes_proj = (TH1D*) DIS_noVetoes->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_DIS_events = (TH1D*) DIS_noVetoes->Get("h_Nevents");

    TFile* DIS_allVetoes = TFile::Open("DIS_beagle_25_10_2_eAu_allVetoes_10x100v4.root","READ");
    TH1D* h_DIS_allVetoes_MC = (TH1D*) DIS_allVetoes->Get("h_t_MC");
    TH1D* h_DIS_allVetoes_L = (TH1D*) DIS_allVetoes->Get("h_t_REC_EEMC_cut");
    TH1D* h_DIS_allVetoes_proj = (TH1D*) DIS_allVetoes->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_DIS_events2 = (TH1D*) DIS_allVetoes->Get("h_Nevents");

    TFile* DIS_onlyDetectors = TFile::Open("DIS_beagle_25_10_2_eAu_allDetectorVetoes_10x100v4.root");
    TH1D* h_DIS_onlyDetectors_MC = (TH1D*) DIS_onlyDetectors->Get("h_t_MC");
    TH1D* h_DIS_onlyDetectors_L = (TH1D*) DIS_onlyDetectors->Get("h_t_REC_EEMC_cut");
    TH1D* h_DIS_onlyDetectors_proj = (TH1D*) DIS_onlyDetectors->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_DIS_events3 = (TH1D*) DIS_onlyDetectors->Get("h_Nevents");

    TFile* DIS_etaOnly = TFile::Open("DIS_beagle_25_10_2_eAu_etaVetoes_10x100v4.root");
    TH1D* h_DIS_etaOnly_MC = (TH1D*) DIS_etaOnly->Get("h_t_MC");
    TH1D* h_DIS_etaOnly_L = (TH1D*) DIS_etaOnly->Get("h_t_REC_EEMC_cut");
    TH1D* h_DIS_etaOnly_proj = (TH1D*) DIS_etaOnly->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_DIS_events4 = (TH1D*) DIS_etaOnly->Get("h_Nevents");

    TFile* DIS_OMDonly = TFile::Open("DIS_beagle_25_10_2_eAu_omdVetoes_10x100v4.root");
    TH1D* h_DIS_OMDonly_MC = (TH1D*) DIS_OMDonly->Get("h_t_MC");
    TH1D* h_DIS_OMDonly_L = (TH1D*) DIS_OMDonly->Get("h_t_REC_EEMC_cut");
    TH1D* h_DIS_OMDonly_proj = (TH1D*) DIS_OMDonly->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_DIS_events5 = (TH1D*) DIS_OMDonly->Get("h_Nevents");

    TFile* DIS_RPonly = TFile::Open("DIS_beagle_25_10_2_eAu_rpVetoes_10x100v4.root");
    TH1D* h_DIS_RPonly_MC = (TH1D*) DIS_RPonly->Get("h_t_MC");
    TH1D* h_DIS_RPonly_L = (TH1D*) DIS_RPonly->Get("h_t_REC_EEMC_cut");
    TH1D* h_DIS_RPonly_proj = (TH1D*) DIS_RPonly->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_DIS_events6 = (TH1D*) DIS_RPonly->Get("h_Nevents");

    double dN_dt = 0.002;
    double sigma_coherent = 459.05; // nb
    double sigma_DIS = 50184.8; // nb
    double nEvents_coherent = h_phi_sartre_events->GetEntries();//6.2e6;
    double nEvents_DIS = h_DIS_events->GetEntries();//1.1e6;
    double nEvents_DIS2 = h_DIS_events2->GetEntries();//1.1e6;
    double nEvents_DIS3 = h_DIS_events3->GetEntries();//1.1e6;
    double nEvents_DIS4 = h_DIS_events4->GetEntries();//1.1e6;
    double nEvents_DIS5 = h_DIS_events5->GetEntries();//1.1e6;
    double nEvents_DIS6 = h_DIS_events6->GetEntries();//1.1e6;
    double preTDR_lumi = 50761.4; //nb^-1
    int rebin_width = 2;   

    TCanvas* c14_213 = new TCanvas("c14_213","c14_213",1,1,1200,800);
	gStyle->SetOptStat(0);
    gPad->SetLogy(1);
    h_t_MC->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
	h_t_MC->GetXaxis()->SetTitle("|t| [GeV/c]^{2}");
    h_t_MC->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    double simu_coherent_lumi = nEvents_coherent/sigma_coherent; // nb^-1 
    double ratio_coherent_phi_preTDR = preTDR_lumi/simu_coherent_lumi;
    for (int i = 1; i <= h_t_MC->GetNbinsX(); ++i) 
    {
        double binError = h_t_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_MC->SetBinError(i, newError);
    }
    h_t_MC->SetLineColor(kBlack);
	h_t_MC->Draw("HISTsame");

    h_t_REC_wRES_cut_pi12->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_wRES_cut_pi12->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_t_REC_wRES_cut_pi12->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_wRES_cut_pi12->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_wRES_cut_pi12->SetBinError(i, newError);
    }
	h_t_REC_wRES_cut_pi12->SetMarkerStyle(30);
	h_t_REC_wRES_cut_pi12->SetMarkerColor(kP8Pink);
    h_t_REC_wRES_cut_pi12->SetLineColor(kP8Pink);
	h_t_REC_wRES_cut_pi12->Draw("PEsame");

    h_DIS_noVetoes_proj->Rebin(rebin_width);
    h_DIS_noVetoes_proj->Scale(sigma_DIS*(1/nEvents_DIS)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_coherent_lumi1 = nEvents_DIS/sigma_DIS; // nb^-1 
    double ratio_coherent_phi_preTDR1 = preTDR_lumi/simu_coherent_lumi1;
    for (int i = 1; i <= h_DIS_noVetoes_proj->GetNbinsX(); ++i) 
    {
        double binError = h_DIS_noVetoes_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR1);
        h_DIS_noVetoes_proj->SetBinError(i, newError);
    }
    h_DIS_noVetoes_proj->SetMarkerStyle(2);
    h_DIS_noVetoes_proj->SetMarkerColor(kP8Cyan);
    h_DIS_noVetoes_proj->SetLineColor(kP8Cyan);
    h_DIS_noVetoes_proj->Draw("PEsame");

    h_DIS_allVetoes_proj->Rebin(rebin_width);
    h_DIS_allVetoes_proj->Scale(sigma_DIS*(1/nEvents_DIS2)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_coherent_lumi2 = nEvents_DIS2/sigma_DIS; // nb^-1 
    double ratio_coherent_phi_preTDR2 = preTDR_lumi/simu_coherent_lumi2;
    for (int i = 1; i <= h_DIS_allVetoes_proj->GetNbinsX(); ++i) 
    {
        double binError = h_DIS_allVetoes_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR2);
        h_DIS_allVetoes_proj->SetBinError(i, newError);
    }
    h_DIS_allVetoes_proj->SetMarkerStyle(42);
    h_DIS_allVetoes_proj->SetMarkerColor(kP8Gray);
    h_DIS_allVetoes_proj->SetLineColor(kP8Gray);
    h_DIS_allVetoes_proj->Draw("PEsame");

    h_DIS_onlyDetectors_proj->Rebin(rebin_width);
    h_DIS_onlyDetectors_proj->Scale(sigma_DIS*(1/nEvents_DIS3)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_coherent_lumi3 = nEvents_DIS3/sigma_DIS; // nb^-1 
    double ratio_coherent_phi_preTDR3 = preTDR_lumi/simu_coherent_lumi3;
    for (int i = 1; i <= h_DIS_onlyDetectors_proj->GetNbinsX(); ++i) 
    {
        double binError = h_DIS_onlyDetectors_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR3);
        h_DIS_onlyDetectors_proj->SetBinError(i, newError);
    }
    h_DIS_onlyDetectors_proj->SetMarkerStyle(31);
    h_DIS_onlyDetectors_proj->SetMarkerColor(kP8Orange);
    h_DIS_onlyDetectors_proj->SetLineColor(kP8Orange);
    h_DIS_onlyDetectors_proj->Draw("PEsame");

    h_DIS_etaOnly_proj->Rebin(rebin_width);
    h_DIS_etaOnly_proj->Scale(sigma_DIS*(1/nEvents_DIS4)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_coherent_lumi4 = nEvents_DIS4/sigma_DIS; // nb^-1 
    double ratio_coherent_phi_preTDR4 = preTDR_lumi/simu_coherent_lumi4;
    for (int i = 1; i <= h_DIS_etaOnly_proj->GetNbinsX(); ++i) 
    {
        double binError = h_DIS_etaOnly_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR4);
        h_DIS_etaOnly_proj->SetBinError(i, newError);
    }
    h_DIS_etaOnly_proj->SetMarkerStyle(4);
    h_DIS_etaOnly_proj->SetMarkerColor(kP10Violet);
    h_DIS_etaOnly_proj->SetLineColor(kP10Violet);
    h_DIS_etaOnly_proj->Draw("PEsame");

    h_DIS_OMDonly_proj->Rebin(rebin_width);
    h_DIS_OMDonly_proj->Scale(sigma_DIS*(1/nEvents_DIS5)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_coherent_lumi5 = nEvents_DIS5/sigma_DIS; // nb^-1 
    double ratio_coherent_phi_preTDR5 = preTDR_lumi/simu_coherent_lumi5;
    for (int i = 1; i <= h_DIS_OMDonly_proj->GetNbinsX(); ++i) 
    {
        double binError = h_DIS_OMDonly_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR5);
        h_DIS_OMDonly_proj->SetBinError(i, newError);
    }
    h_DIS_OMDonly_proj->SetMarkerStyle(21);
    h_DIS_OMDonly_proj->SetMarkerColor(kP8Green);
    h_DIS_OMDonly_proj->SetLineColor(kP8Green);
    h_DIS_OMDonly_proj->Draw("PEsame");

    h_DIS_RPonly_proj->Rebin(rebin_width);
    h_DIS_RPonly_proj->Scale(sigma_DIS*(1/nEvents_DIS6)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_coherent_lumi6 = nEvents_DIS6/sigma_DIS; // nb^-1 
    double ratio_coherent_phi_preTDR6 = preTDR_lumi/simu_coherent_lumi6;
    for (int i = 1; i <= h_DIS_RPonly_proj->GetNbinsX(); ++i) 
    {
        double binError = h_DIS_RPonly_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR6);
        h_DIS_RPonly_proj->SetBinError(i, newError);
    }
    h_DIS_RPonly_proj->SetMarkerStyle(25);
    h_DIS_RPonly_proj->SetMarkerColor(kP8Red);
    h_DIS_RPonly_proj->SetLineColor(kP8Red);
    h_DIS_RPonly_proj->Draw("PEsame");
    
	TLegend *w14_213 = new TLegend(0.6,0.55,0.72,0.9);
	w14_213->AddEntry(h_t_MC,"#phi: MC", "L");
    w14_213->AddEntry(h_t_REC_wRES_cut_pi12,"#phi: #theta_{max}= #pi/12", "P");    
    w14_213->AddEntry(h_DIS_noVetoes_proj,"DIS: no vetoes", "P");
    w14_213->AddEntry(h_DIS_OMDonly_proj,"DIS: OMD only", "P");
    w14_213->AddEntry(h_DIS_RPonly_proj,"DIS: RP only", "P");
    w14_213->AddEntry(h_DIS_onlyDetectors_proj,"DIS: only detectors", "P");
    w14_213->AddEntry(h_DIS_etaOnly_proj,"DIS: #eta only", "P");
    w14_213->AddEntry(h_DIS_allVetoes_proj,"DIS: all vetoes", "P");
	w14_213->SetBorderSize(0);   
    w14_213->SetFillStyle(0);
    w14_213->SetTextSize(30);
	w14_213->SetTextFont(45);
    w14_213->Draw("same");	

    TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

    // Add labels
    TLatex* ep = new TLatex(0.13, 0.85, "ePIC");
    ep->SetNDC();
    ep->SetTextFont(62);
    ep->SetTextSize(.038);
    ep->SetTextColor(kBlack);
    ep->Draw("same");

    TLatex* r421 = new TLatex(0.18, 0.85, " Simulation 25.10.2, 10x100 GeV");
	r421->SetNDC();
	r421->SetTextSize(30);
	r421->SetTextFont(43);
	r421->SetTextColor(kBlack);
	r421->Draw("same");

	TLatex* label_5 = new TLatex(0.17, 0.75, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->SetTextSize(30);
	label_5->SetTextFont(43);
	label_5->SetTextColor(kBlack);
	label_5->Draw("same");

    TLatex* r4211 = new TLatex(0.15, 0.8, "eAu #rightarrow e'Au'#phi");
	r4211->SetNDC();
	r4211->SetTextSize(30);
	r4211->SetTextFont(43);
	r4211->SetTextColor(kBlack);
	r4211->Draw("same"); 

    TLatex* r421111 = new TLatex(0.15, 0.15, "L_{int} = 10 fb^{-1}/A"); // preTDR lumi
	r421111->SetNDC();
	r421111->SetTextSize(25);
	r421111->SetTextFont(43);
	r421111->SetTextColor(kBlack);
	r421111->Draw("same");

    c14_213->Print("./figures/plot_t_DIS_compareALL_preTDR_proj.pdf");
}

void plot_DIS_together_ES()
{
    TFile* phi_t_file = TFile::Open("phi_sartre_25_10_2_eAu_10x100v5.root","READ");
	TH1D* h_t_MC = (TH1D*) phi_t_file->Get("h_t_MC");
	TH1D* h_t_REC_EEMC_cut = (TH1D*) phi_t_file->Get("h_t_REC_EEMC_cut");
	TH1D* h_t_REC_wRES_cut_pi12 = (TH1D*) phi_t_file->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_phi_sartre_events = (TH1D*) phi_t_file->Get("h_Nevents");

	TFile* DIS_noVetoes = TFile::Open("DIS_beagle_25_10_2_eAu_noVetoes_10x100v4.root","READ");
    TH1D* h_DIS_noVetoes_MC = (TH1D*) DIS_noVetoes->Get("h_t_MC");
    TH1D* h_DIS_noVetoes_L = (TH1D*) DIS_noVetoes->Get("h_t_REC_EEMC_cut");
    TH1D* h_DIS_noVetoes_proj = (TH1D*) DIS_noVetoes->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_DIS_events = (TH1D*) DIS_noVetoes->Get("h_Nevents");

    TFile* DIS_allVetoes = TFile::Open("DIS_beagle_25_10_2_eAu_allVetoes_10x100v4.root","READ");
    TH1D* h_DIS_allVetoes_MC = (TH1D*) DIS_allVetoes->Get("h_t_MC");
    TH1D* h_DIS_allVetoes_L = (TH1D*) DIS_allVetoes->Get("h_t_REC_EEMC_cut");
    TH1D* h_DIS_allVetoes_proj = (TH1D*) DIS_allVetoes->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_DIS_events2 = (TH1D*) DIS_allVetoes->Get("h_Nevents");

    TFile* DIS_onlyDetectors = TFile::Open("DIS_beagle_25_10_2_eAu_allDetectorVetoes_10x100v4.root");
    TH1D* h_DIS_onlyDetectors_MC = (TH1D*) DIS_onlyDetectors->Get("h_t_MC");
    TH1D* h_DIS_onlyDetectors_L = (TH1D*) DIS_onlyDetectors->Get("h_t_REC_EEMC_cut");
    TH1D* h_DIS_onlyDetectors_proj = (TH1D*) DIS_onlyDetectors->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_DIS_events3 = (TH1D*) DIS_onlyDetectors->Get("h_Nevents");

    TFile* DIS_etaOnly = TFile::Open("DIS_beagle_25_10_2_eAu_etaVetoes_10x100v4.root");
    TH1D* h_DIS_etaOnly_MC = (TH1D*) DIS_etaOnly->Get("h_t_MC");
    TH1D* h_DIS_etaOnly_L = (TH1D*) DIS_etaOnly->Get("h_t_REC_EEMC_cut");
    TH1D* h_DIS_etaOnly_proj = (TH1D*) DIS_etaOnly->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_DIS_events4 = (TH1D*) DIS_etaOnly->Get("h_Nevents");

    TFile* DIS_OMDonly = TFile::Open("DIS_beagle_25_10_2_eAu_omdVetoes_10x100v4.root");
    TH1D* h_DIS_OMDonly_MC = (TH1D*) DIS_OMDonly->Get("h_t_MC");
    TH1D* h_DIS_OMDonly_L = (TH1D*) DIS_OMDonly->Get("h_t_REC_EEMC_cut");
    TH1D* h_DIS_OMDonly_proj = (TH1D*) DIS_OMDonly->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_DIS_events5 = (TH1D*) DIS_OMDonly->Get("h_Nevents");

    TFile* DIS_RPonly = TFile::Open("DIS_beagle_25_10_2_eAu_rpVetoes_10x100v4.root");
    TH1D* h_DIS_RPonly_MC = (TH1D*) DIS_RPonly->Get("h_t_MC");
    TH1D* h_DIS_RPonly_L = (TH1D*) DIS_RPonly->Get("h_t_REC_EEMC_cut");
    TH1D* h_DIS_RPonly_proj = (TH1D*) DIS_RPonly->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_DIS_events6 = (TH1D*) DIS_RPonly->Get("h_Nevents");

    double dN_dt = 0.002;
    double sigma_coherent = 459.05; // nb
    double sigma_DIS = 50184.8; // nb
    double nEvents_coherent = h_phi_sartre_events->GetEntries();//6.2e6;
    double nEvents_DIS = h_DIS_events->GetEntries();//1.1e6;
    double nEvents_DIS2 = h_DIS_events2->GetEntries();//1.1e6;
    double nEvents_DIS3 = h_DIS_events3->GetEntries();//1.1e6;
    double nEvents_DIS4 = h_DIS_events4->GetEntries();//1.1e6;
    double nEvents_DIS5 = h_DIS_events5->GetEntries();//1.1e6;
    double nEvents_DIS6 = h_DIS_events6->GetEntries();//1.1e6;
    double early_science_lumi = 8121.83; // nb^-1
    int rebin_width = 2;   

    TCanvas* c14_213 = new TCanvas("c14_213","c14_213",1,1,1200,800);
	gStyle->SetOptStat(0);
    gPad->SetLogy(1);
    h_t_MC->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
	h_t_MC->GetXaxis()->SetTitle("|t| [GeV/c]^{2}");
    h_t_MC->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    double simu_coherent_lumi = nEvents_coherent/sigma_coherent; // nb^-1 
    double ratio_coherent_phi_ES = early_science_lumi/simu_coherent_lumi;
    for (int i = 1; i <= h_t_MC->GetNbinsX(); ++i) 
    {
        double binError = h_t_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_MC->SetBinError(i, newError);
    }
    h_t_MC->SetLineColor(kBlack);
	h_t_MC->Draw("HISTsame");

    h_t_REC_wRES_cut_pi12->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_wRES_cut_pi12->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_t_REC_wRES_cut_pi12->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_wRES_cut_pi12->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_REC_wRES_cut_pi12->SetBinError(i, newError);
    }
	h_t_REC_wRES_cut_pi12->SetMarkerStyle(30);
	h_t_REC_wRES_cut_pi12->SetMarkerColor(kP8Pink);
    h_t_REC_wRES_cut_pi12->SetLineColor(kP8Pink);
	h_t_REC_wRES_cut_pi12->Draw("PEsame");

    h_DIS_noVetoes_proj->Rebin(rebin_width);
    h_DIS_noVetoes_proj->Scale(sigma_DIS*(1/nEvents_DIS)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_coherent_lumi1 = nEvents_DIS/sigma_DIS; // nb^-1 
    double ratio_DIS_ES1 = early_science_lumi/simu_coherent_lumi1;
    for (int i = 1; i <= h_DIS_noVetoes_proj->GetNbinsX(); ++i) 
    {
        double binError = h_DIS_noVetoes_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_DIS_ES1);
        h_DIS_noVetoes_proj->SetBinError(i, newError);
    }
    h_DIS_noVetoes_proj->SetMarkerStyle(2);
    h_DIS_noVetoes_proj->SetMarkerColor(kP8Cyan);
    h_DIS_noVetoes_proj->SetLineColor(kP8Cyan);
    h_DIS_noVetoes_proj->Draw("PEsame");

    h_DIS_allVetoes_proj->Rebin(rebin_width);
    h_DIS_allVetoes_proj->Scale(sigma_DIS*(1/nEvents_DIS2)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_coherent_lumi2 = nEvents_DIS2/sigma_DIS; // nb^-1 
    double ratio_DIS_ES2 = early_science_lumi/simu_coherent_lumi2;
    for (int i = 1; i <= h_DIS_allVetoes_proj->GetNbinsX(); ++i) 
    {
        double binError = h_DIS_allVetoes_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_DIS_ES2);
        h_DIS_allVetoes_proj->SetBinError(i, newError);
    }
    h_DIS_allVetoes_proj->SetMarkerStyle(42);
    h_DIS_allVetoes_proj->SetMarkerColor(kP8Gray);
    h_DIS_allVetoes_proj->SetLineColor(kP8Gray);
    h_DIS_allVetoes_proj->Draw("PEsame");

    h_DIS_onlyDetectors_proj->Rebin(rebin_width);
    h_DIS_onlyDetectors_proj->Scale(sigma_DIS*(1/nEvents_DIS3)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_coherent_lumi3 = nEvents_DIS3/sigma_DIS; // nb^-1 
    double ratio_DIS_ES3 = early_science_lumi/simu_coherent_lumi3;
    for (int i = 1; i <= h_DIS_onlyDetectors_proj->GetNbinsX(); ++i) 
    {
        double binError = h_DIS_onlyDetectors_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_DIS_ES3);
        h_DIS_onlyDetectors_proj->SetBinError(i, newError);
    }
    h_DIS_onlyDetectors_proj->SetMarkerStyle(31);
    h_DIS_onlyDetectors_proj->SetMarkerColor(kP8Orange);
    h_DIS_onlyDetectors_proj->SetLineColor(kP8Orange);
    h_DIS_onlyDetectors_proj->Draw("PEsame");

    h_DIS_etaOnly_proj->Rebin(rebin_width);
    h_DIS_etaOnly_proj->Scale(sigma_DIS*(1/nEvents_DIS4)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_coherent_lumi4 = nEvents_DIS4/sigma_DIS; // nb^-1 
    double ratio_DIS_ES4 = early_science_lumi/simu_coherent_lumi4;
    for (int i = 1; i <= h_DIS_etaOnly_proj->GetNbinsX(); ++i) 
    {
        double binError = h_DIS_etaOnly_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_DIS_ES4);
        h_DIS_etaOnly_proj->SetBinError(i, newError);
    }
    h_DIS_etaOnly_proj->SetMarkerStyle(4);
    h_DIS_etaOnly_proj->SetMarkerColor(kP10Violet);
    h_DIS_etaOnly_proj->SetLineColor(kP10Violet);
    h_DIS_etaOnly_proj->Draw("PEsame");

    h_DIS_OMDonly_proj->Rebin(rebin_width);
    h_DIS_OMDonly_proj->Scale(sigma_DIS*(1/nEvents_DIS5)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_coherent_lumi5 = nEvents_DIS5/sigma_DIS; // nb^-1 
    double ratio_DIS_ES5 = early_science_lumi/simu_coherent_lumi5;
    for (int i = 1; i <= h_DIS_OMDonly_proj->GetNbinsX(); ++i) 
    {
        double binError = h_DIS_OMDonly_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_DIS_ES5);
        h_DIS_OMDonly_proj->SetBinError(i, newError);
    }
    h_DIS_OMDonly_proj->SetMarkerStyle(21);
    h_DIS_OMDonly_proj->SetMarkerColor(kP8Green);
    h_DIS_OMDonly_proj->SetLineColor(kP8Green);
    h_DIS_OMDonly_proj->Draw("PEsame");

    h_DIS_RPonly_proj->Rebin(rebin_width);
    h_DIS_RPonly_proj->Scale(sigma_DIS*(1/nEvents_DIS6)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_coherent_lumi6 = nEvents_DIS6/sigma_DIS; // nb^-1 
    double ratio_DIS_ES6 = early_science_lumi/simu_coherent_lumi6;
    for (int i = 1; i <= h_DIS_RPonly_proj->GetNbinsX(); ++i) 
    {
        double binError = h_DIS_RPonly_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_DIS_ES6);
        h_DIS_RPonly_proj->SetBinError(i, newError);
    }
    h_DIS_RPonly_proj->SetMarkerStyle(25);
    h_DIS_RPonly_proj->SetMarkerColor(kP8Red);
    h_DIS_RPonly_proj->SetLineColor(kP8Red);
    h_DIS_RPonly_proj->Draw("PEsame");
    
	TLegend *w14_213 = new TLegend(0.6,0.55,0.72,0.9);
	w14_213->AddEntry(h_t_MC,"#phi: MC", "L");
    w14_213->AddEntry(h_t_REC_wRES_cut_pi12,"#phi: #theta_{max}= #pi/12", "P");    
    w14_213->AddEntry(h_DIS_noVetoes_proj,"DIS: no vetoes", "P");
    w14_213->AddEntry(h_DIS_OMDonly_proj,"DIS: OMD only", "P");
    w14_213->AddEntry(h_DIS_RPonly_proj,"DIS: RP only", "P");
    w14_213->AddEntry(h_DIS_onlyDetectors_proj,"DIS: only detectors", "P");
    w14_213->AddEntry(h_DIS_etaOnly_proj,"DIS: #eta only", "P");
    w14_213->AddEntry(h_DIS_allVetoes_proj,"DIS: all vetoes", "P");
	w14_213->SetBorderSize(0);   
    w14_213->SetFillStyle(0);
    w14_213->SetTextSize(30);
	w14_213->SetTextFont(45);
    w14_213->Draw("same");	

    TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

    // Add labels
    TLatex* ep = new TLatex(0.13, 0.85, "ePIC");
    ep->SetNDC();
    ep->SetTextFont(62);
    ep->SetTextSize(.038);
    ep->SetTextColor(kBlack);
    ep->Draw("same");

    TLatex* r421 = new TLatex(0.18, 0.85, " Simulation 25.10.2, 10x100 GeV");
	r421->SetNDC();
	r421->SetTextSize(30);
	r421->SetTextFont(43);
	r421->SetTextColor(kBlack);
	r421->Draw("same");

	TLatex* label_5 = new TLatex(0.17, 0.75, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->SetTextSize(30);
	label_5->SetTextFont(43);
	label_5->SetTextColor(kBlack);
	label_5->Draw("same");

    TLatex* r4211 = new TLatex(0.15, 0.8, "eAu #rightarrow e'Au'#phi");
	r4211->SetNDC();
	r4211->SetTextSize(30);
	r4211->SetTextFont(43);
	r4211->SetTextColor(kBlack);
	r4211->Draw("same"); 

    TLatex* r42111 = new TLatex(0.15, 0.15, "L_{int} = 1.6 fb^{-1}/A"); // ES lumi
	r42111->SetNDC();
	r42111->SetTextSize(30);
	r42111->SetTextFont(43);
	r42111->SetTextColor(kBlack);
	r42111->Draw("same");

    c14_213->Print("./figures/plot_t_DIS_compareALL_ES_proj.pdf");
}

void plot_incoherent_phi_allVetoes_preTDR()
{
	TFile* phi_t_file = TFile::Open("phi_sartre_25_10_2_eAu_10x100v5.root","READ");
	TH1D* h_t_MC = (TH1D*) phi_t_file->Get("h_t_MC");
	TH1D* h_t_REC_EEMC_cut = (TH1D*) phi_t_file->Get("h_t_REC_EEMC_cut");
	TH1D* h_t_REC_wRES_cut_pi12 = (TH1D*) phi_t_file->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_phi_sartre_events = (TH1D*) phi_t_file->Get("h_Nevents");

	TFile* incoherent_phi_allVetoes = TFile::Open("phi_beagle_25_10_3_eAu_allVetoes_10x100v4.root","READ");
    TH1D* h_incoherent_phi_allVetoes_MC = (TH1D*) incoherent_phi_allVetoes->Get("h_t_MC");
    TH1D* h_incoherent_phi_allVetoes_L = (TH1D*) incoherent_phi_allVetoes->Get("h_t_REC_EEMC_cut");
    TH1D* h_incoherent_phi_allVetoes_proj = (TH1D*) incoherent_phi_allVetoes->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_incoherent_events = (TH1D*) incoherent_phi_allVetoes->Get("h_Nevents");

	double dN_dt = 0.002;
    double sigma_coherent = 459.05; // nb
    double sigma_incoherent = 261.95; // nb
    double nEvents_coherent = h_phi_sartre_events->GetEntries();
    double nEvents_incoherent = h_incoherent_events->GetEntries();
    double preTDR_lumi = 50761.4; // nb^-1
    int rebin_width = 2;   
	
	TCanvas* c14_2131 = new TCanvas("c14_2131","c14_2131",1,1,1200,800);
    gPad->SetLogy(1);
	gStyle->SetOptStat(0);
    h_t_MC->GetYaxis()->SetTitle("dN/d|t| [GeV/c]^{-2}");
	h_t_MC->GetXaxis()->SetTitle("|t| [GeV/c]^{2}");
    h_t_MC->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    double simu_coherent_lumi = nEvents_coherent/sigma_coherent; // nb^-1 
    double ratio_coherent_phi_preTDR = preTDR_lumi/simu_coherent_lumi;
    for (int i = 1; i <= h_t_MC->GetNbinsX(); ++i) 
    {
        double binError = h_t_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_MC->SetBinError(i, newError);
    }
    h_t_MC->SetLineColor(kBlack);
	h_t_MC->Draw("HISTsame");

    h_t_REC_EEMC_cut->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    for (int i = 1; i <= h_t_REC_EEMC_cut->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_EEMC_cut->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_EEMC_cut->SetBinError(i, newError);
    }
	h_t_REC_EEMC_cut->SetMarkerStyle(20); // method L RECO
	h_t_REC_EEMC_cut->SetMarkerColor(kP8Blue);
	//h_t_REC_EEMC_cut->Draw("PEsame");

    h_t_REC_wRES_cut_pi12->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_wRES_cut_pi12->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_t_REC_wRES_cut_pi12->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_wRES_cut_pi12->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_wRES_cut_pi12->SetBinError(i, newError);
    }
	h_t_REC_wRES_cut_pi12->SetMarkerStyle(30);
	h_t_REC_wRES_cut_pi12->SetMarkerColor(kP8Pink);
    h_t_REC_wRES_cut_pi12->SetLineColor(kP8Pink);
	h_t_REC_wRES_cut_pi12->Draw("PEsame");

    h_incoherent_phi_allVetoes_MC->Rebin(rebin_width);
    h_incoherent_phi_allVetoes_MC->Scale(sigma_incoherent*(1/nEvents_incoherent)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_incoherent_lumi = nEvents_incoherent/sigma_incoherent; // nb^-1 
    double ratio_incoherent_preTDR = preTDR_lumi/simu_incoherent_lumi;
    for (int i = 1; i <= h_incoherent_phi_allVetoes_MC->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_allVetoes_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_preTDR);
        h_incoherent_phi_allVetoes_MC->SetBinError(i, newError);
    }
    h_incoherent_phi_allVetoes_MC->SetMarkerStyle(21);
    h_incoherent_phi_allVetoes_MC->SetMarkerColor(kP8Orange);
    h_incoherent_phi_allVetoes_MC->SetLineColor(kP8Orange);
    //h_incoherent_phi_allVetoes_MC->Draw("PEsame");

    h_incoherent_phi_allVetoes_L->Rebin(rebin_width);
    h_incoherent_phi_allVetoes_L->Scale(sigma_incoherent*(1/nEvents_incoherent)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_incoherent_phi_allVetoes_L->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_allVetoes_L->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_preTDR);
        h_incoherent_phi_allVetoes_L->SetBinError(i, newError);
    }
    h_incoherent_phi_allVetoes_L->SetMarkerStyle(21);
    h_incoherent_phi_allVetoes_L->SetMarkerColor(kP8Orange);
    h_incoherent_phi_allVetoes_L->SetLineColor(kP8Orange);
    //h_incoherent_phi_allVetoes_L->Draw("PEsame");

    h_incoherent_phi_allVetoes_proj->Rebin(rebin_width);
    h_incoherent_phi_allVetoes_proj->Scale(sigma_incoherent*(1/nEvents_incoherent)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_incoherent_phi_allVetoes_proj->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_allVetoes_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_preTDR);
        h_incoherent_phi_allVetoes_proj->SetBinError(i, newError);
    }
    h_incoherent_phi_allVetoes_proj->SetMarkerStyle(21);
    h_incoherent_phi_allVetoes_proj->SetMarkerColor(kP8Orange);
    h_incoherent_phi_allVetoes_proj->SetLineColor(kP8Orange);
    h_incoherent_phi_allVetoes_proj->Draw("PEsame");

	TLegend *w14_2131 = new TLegend(0.6,0.7,0.72,0.9);
    w14_2131->AddEntry(h_t_MC,"Coherent #phi: MC", "L");
    w14_2131->AddEntry(h_t_REC_wRES_cut_pi12,"Coherent #phi: #theta_{max}= #pi/12", "P");    
    //w14_2131->AddEntry(h_t_REC_EEMC_cut,"Coherent #phi: L", "P");
    //w14_2131->AddEntry(h_incoherent_phi_allVetoes_MC,"Incoherent #phi: MC", "P");
    //w14_2131->AddEntry(h_incoherent_phi_allVetoes_L,"Incoherent #phi: L", "P");
    w14_2131->AddEntry(h_incoherent_phi_allVetoes_proj,"Incoherent: #theta_{max}= #pi/12", "P");
	w14_2131->SetBorderSize(0);   
    w14_2131->SetFillStyle(0);
    w14_2131->SetTextSize(30);
	w14_2131->SetTextFont(45);
    w14_2131->Draw("same");	

    c14_2131->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title14_2131 = new TLatex();
    title14_2131->SetNDC(); 
    title14_2131->SetTextSize(0.05);
    title14_2131->SetTextAlign(22);  
    //title14_2131->DrawLatex(0.5, 0.97, "|t| Distribution with Incoherent #phi All Vetoes");  

    TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

    // Add labels
    TLatex* ep = new TLatex(0.13, 0.85, "ePIC");
    ep->SetNDC();
    ep->SetTextFont(62);
    ep->SetTextSize(.038);
    ep->SetTextColor(kBlack);
    ep->Draw("same");

    TLatex* r421 = new TLatex(0.18, 0.85, " Simulation 25.10.2, 10x100 GeV");
	r421->SetNDC();
	r421->SetTextSize(30);
	r421->SetTextFont(43);
	r421->SetTextColor(kBlack);
	r421->Draw("same");

	TLatex* label_5 = new TLatex(0.17, 0.75, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->SetTextSize(30);
	label_5->SetTextFont(43);
	label_5->SetTextColor(kBlack);
	label_5->Draw("same");

    TLatex* r4211 = new TLatex(0.15, 0.8, "eAu #rightarrow e'Au'#phi");
	r4211->SetNDC();
	r4211->SetTextSize(30);
	r4211->SetTextFont(43);
	r4211->SetTextColor(kBlack);
	r4211->Draw("same"); 

    TLatex* r421111 = new TLatex(0.18, 0.13, "L_{int} = 10 fb^{-1}/A"); // preTDR lumi
	r421111->SetNDC();
	r421111->SetTextSize(25);
	r421111->SetTextFont(43);
	r421111->SetTextColor(kBlack);
	r421111->Draw("same");

    c14_2131->Print("./figures/plot_t_incoherent_allvetoes_preTDR_proj.pdf");

	//TFile* fout_preTDR = new TFile("plot_t_incoherent_allvetoes_preTDR.root","RECREATE");
    //c14_2131->Write();
    //h_t_MC->Write();
    //h_t_REC_EEMC_cut->Write();
    //h_t_REC_wRES_cut_pi12->Write();
    //h_incoherent_phi_allVetoes_MC->Write();
    //h_incoherent_phi_allVetoes_L->Write();
    //h_incoherent_phi_allVetoes_proj->Write();
    //fout_preTDR->Close();
}

void plot_incoherent_phi_allVetoes_ES()
{
	TFile* phi_t_file = TFile::Open("phi_sartre_25_10_2_eAu_10x100v5.root","READ");
	TH1D* h_t_MC = (TH1D*) phi_t_file->Get("h_t_MC");
	TH1D* h_t_REC_EEMC_cut = (TH1D*) phi_t_file->Get("h_t_REC_EEMC_cut");
	TH1D* h_t_REC_wRES_cut_pi12 = (TH1D*) phi_t_file->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_phi_sartre_events = (TH1D*) phi_t_file->Get("h_Nevents");

	TFile* incoherent_phi_allVetoes = TFile::Open("phi_beagle_25_10_3_eAu_allVetoes_10x100v4.root","READ");
    TH1D* h_incoherent_phi_allVetoes_MC = (TH1D*) incoherent_phi_allVetoes->Get("h_t_MC");
    TH1D* h_incoherent_phi_allVetoes_L = (TH1D*) incoherent_phi_allVetoes->Get("h_t_REC_EEMC_cut");
    TH1D* h_incoherent_phi_allVetoes_proj = (TH1D*) incoherent_phi_allVetoes->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_incoherent_events = (TH1D*) incoherent_phi_allVetoes->Get("h_Nevents");

	double dN_dt = 0.002;
    double sigma_coherent = 459.05; // nb
    double sigma_incoherent = 261.95; // nb
    double nEvents_coherent = h_phi_sartre_events->GetEntries();
    double nEvents_incoherent = h_incoherent_events->GetEntries();
    double early_science_lumi = 8121.83; // nb^-1
    int rebin_width = 2;   
	
	TCanvas* c14_2131 = new TCanvas("c14_2131","c14_2131",1,1,1200,800);
    gPad->SetLogy(1);
	gStyle->SetOptStat(0);
    h_t_MC->GetYaxis()->SetTitle("dN/d|t| [GeV/c]^{-2}");
	h_t_MC->GetXaxis()->SetTitle("|t| [GeV/c]^{2}");
    h_t_MC->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    double simu_coherent_lumi = nEvents_coherent/sigma_coherent; // nb^-1 
    double ratio_coherent_phi_ES = early_science_lumi/simu_coherent_lumi;
    for (int i = 1; i <= h_t_MC->GetNbinsX(); ++i) 
    {
        double binError = h_t_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_MC->SetBinError(i, newError);
    }
    h_t_MC->SetLineColor(kBlack);
	h_t_MC->Draw("HISTsame");

    h_t_REC_EEMC_cut->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    for (int i = 1; i <= h_t_REC_EEMC_cut->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_EEMC_cut->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_REC_EEMC_cut->SetBinError(i, newError);
    }
	h_t_REC_EEMC_cut->SetMarkerStyle(20); // method L RECO
	h_t_REC_EEMC_cut->SetMarkerColor(kP8Blue);
	//h_t_REC_EEMC_cut->Draw("PEsame");

    h_t_REC_wRES_cut_pi12->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_wRES_cut_pi12->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_t_REC_wRES_cut_pi12->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_wRES_cut_pi12->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_REC_wRES_cut_pi12->SetBinError(i, newError);
    }
	h_t_REC_wRES_cut_pi12->SetMarkerStyle(30);
	h_t_REC_wRES_cut_pi12->SetMarkerColor(kP8Pink);
    h_t_REC_wRES_cut_pi12->SetLineColor(kP8Pink);
	h_t_REC_wRES_cut_pi12->Draw("PEsame");

    h_incoherent_phi_allVetoes_MC->Rebin(rebin_width);
    h_incoherent_phi_allVetoes_MC->Scale(sigma_incoherent*(1/nEvents_incoherent)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_incoherent_lumi = nEvents_incoherent/sigma_incoherent; // nb^-1 
    double ratio_incoherent_ES = early_science_lumi/simu_incoherent_lumi;
    for (int i = 1; i <= h_incoherent_phi_allVetoes_MC->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_allVetoes_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_ES);
        h_incoherent_phi_allVetoes_MC->SetBinError(i, newError);
    }
    h_incoherent_phi_allVetoes_MC->SetMarkerStyle(21);
    h_incoherent_phi_allVetoes_MC->SetMarkerColor(kP8Orange);
    h_incoherent_phi_allVetoes_MC->SetLineColor(kP8Orange);
    //h_incoherent_phi_allVetoes_MC->Draw("PEsame");

    h_incoherent_phi_allVetoes_L->Rebin(rebin_width);
    h_incoherent_phi_allVetoes_L->Scale(sigma_incoherent*(1/nEvents_incoherent)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_incoherent_phi_allVetoes_L->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_allVetoes_L->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_ES);
        h_incoherent_phi_allVetoes_L->SetBinError(i, newError);
    }
    h_incoherent_phi_allVetoes_L->SetMarkerStyle(21);
    h_incoherent_phi_allVetoes_L->SetMarkerColor(kP8Orange);
    h_incoherent_phi_allVetoes_L->SetLineColor(kP8Orange);
    //h_incoherent_phi_allVetoes_L->Draw("PEsame");

    h_incoherent_phi_allVetoes_proj->Rebin(rebin_width);
    h_incoherent_phi_allVetoes_proj->Scale(sigma_incoherent*(1/nEvents_incoherent)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_incoherent_phi_allVetoes_proj->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_allVetoes_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_ES);
        h_incoherent_phi_allVetoes_proj->SetBinError(i, newError);
    }
    h_incoherent_phi_allVetoes_proj->SetMarkerStyle(21);
    h_incoherent_phi_allVetoes_proj->SetMarkerColor(kP8Orange);
    h_incoherent_phi_allVetoes_proj->SetLineColor(kP8Orange);
    h_incoherent_phi_allVetoes_proj->Draw("PEsame");

	TLegend *w14_2131 = new TLegend(0.6,0.7,0.72,0.9);
    w14_2131->AddEntry(h_t_MC,"Coherent #phi: MC", "L");
    w14_2131->AddEntry(h_t_REC_wRES_cut_pi12,"Coherent #phi: #theta_{max}= #pi/12", "P");    
    //w14_2131->AddEntry(h_t_REC_EEMC_cut,"Coherent #phi: L", "P");
    //w14_2131->AddEntry(h_incoherent_phi_allVetoes_MC,"Incoherent #phi: MC", "P");
    //w14_2131->AddEntry(h_incoherent_phi_allVetoes_L,"Incoherent #phi: L", "P");
    w14_2131->AddEntry(h_incoherent_phi_allVetoes_proj,"Incoherent: #theta_{max}= #pi/12", "P");
	w14_2131->SetBorderSize(0);   
    w14_2131->SetFillStyle(0);
    w14_2131->SetTextSize(30);
	w14_2131->SetTextFont(45);
    w14_2131->Draw("same");	

    c14_2131->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title14_2131 = new TLatex();
    title14_2131->SetNDC(); 
    title14_2131->SetTextSize(0.05);
    title14_2131->SetTextAlign(22);  
    //title14_2131->DrawLatex(0.5, 0.97, "|t| Distribution with Incoherent #phi All Vetoes");  

    TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

    // Add labels
    TLatex* ep = new TLatex(0.13, 0.85, "ePIC");
    ep->SetNDC();
    ep->SetTextFont(62);
    ep->SetTextSize(.038);
    ep->SetTextColor(kBlack);
    ep->Draw("same");

    TLatex* r421 = new TLatex(0.18, 0.85, " Simulation 25.10.2, 10x100 GeV");
	r421->SetNDC();
	r421->SetTextSize(30);
	r421->SetTextFont(43);
	r421->SetTextColor(kBlack);
	r421->Draw("same");

	TLatex* label_5 = new TLatex(0.17, 0.75, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->SetTextSize(30);
	label_5->SetTextFont(43);
	label_5->SetTextColor(kBlack);
	label_5->Draw("same");

    TLatex* r4211 = new TLatex(0.15, 0.8, "eAu #rightarrow e'Au'#phi");
	r4211->SetNDC();
	r4211->SetTextSize(30);
	r4211->SetTextFont(43);
	r4211->SetTextColor(kBlack);
	r4211->Draw("same"); 

    TLatex* r42111 = new TLatex(0.18, 0.15, "L_{int} = 1.6 fb^{-1}/A"); // ES lumi
	r42111->SetNDC();
	r42111->SetTextSize(30);
	r42111->SetTextFont(43);
	r42111->SetTextColor(kBlack);
	r42111->Draw("same");

    c14_2131->Print("./figures/plot_t_incoherent_allvetoes_ES_proj.pdf");

	//TFile* fout_ES = new TFile("plot_t_incoherent_allvetoes_ES.root","RECREATE");
    //c14_2131->Write();
    //h_t_MC->Write();
    //h_t_REC_EEMC_cut->Write();
    //h_t_REC_wRES_cut_pi12->Write();
    //h_incoherent_phi_allVetoes_MC->Write();
    //h_incoherent_phi_allVetoes_L->Write();
    //h_incoherent_phi_allVetoes_proj->Write();
    //fout_ES->Close();
}

void plot_incoherent_phi_noVetoes_preTDR()
{
	TFile* phi_t_file = TFile::Open("phi_sartre_25_10_2_eAu_10x100v5.root","READ");
	TH1D* h_t_MC = (TH1D*) phi_t_file->Get("h_t_MC");
	TH1D* h_t_REC_EEMC_cut = (TH1D*) phi_t_file->Get("h_t_REC_EEMC_cut");
	TH1D* h_t_REC_wRES_cut_pi12 = (TH1D*) phi_t_file->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_phi_sartre_events = (TH1D*) phi_t_file->Get("h_Nevents");

	TFile* incoherent_phi_noVetoes = TFile::Open("phi_beagle_25_10_3_eAu_noVetoes_10x100v4.root","READ");
    TH1D* h_incoherent_phi_noVetoes_MC = (TH1D*) incoherent_phi_noVetoes->Get("h_t_MC");
    TH1D* h_incoherent_phi_noVetoes_L = (TH1D*) incoherent_phi_noVetoes->Get("h_t_REC_EEMC_cut");
    TH1D* h_incoherent_phi_noVetoes_proj = (TH1D*) incoherent_phi_noVetoes->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_incoherent_phi_noVetoes_projNoScale = (TH1D*) incoherent_phi_noVetoes->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_incoherent_phi_noVetoes_proj_pi2 = (TH1D*) incoherent_phi_noVetoes->Get("h_t_REC_wRES_cut_pi2");
    TH1D* h_incoherent_events = (TH1D*) incoherent_phi_noVetoes->Get("h_Nevents");

	double dN_dt = 0.002;
    double sigma_coherent = 459.05; // nb
    double sigma_incoherent = 261.95; // nb
    double nEvents_coherent = h_phi_sartre_events->GetEntries();
    double nEvents_incoherent = h_incoherent_events->GetEntries();
    double preTDR_lumi = 50761.4; // nb^-1
    int rebin_width = 2;   

	TCanvas* c14_213 = new TCanvas("c14_213","c14_213",1,1,1200,800);
    gPad->SetLogy(1);
	gStyle->SetOptStat(0);
    h_t_MC->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
	h_t_MC->GetXaxis()->SetTitle("|t| [GeV/c]^{2}");
    h_t_MC->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    double simu_coherent_lumi = nEvents_coherent/sigma_coherent; // nb^-1 
    double ratio_coherent_phi_preTDR = preTDR_lumi/simu_coherent_lumi;
    for (int i = 1; i <= h_t_MC->GetNbinsX(); ++i) 
    {
        double binError = h_t_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_MC->SetBinError(i, newError);
    }
    h_t_MC->SetLineColor(kBlack);
	h_t_MC->Draw("HISTsame");

    h_t_REC_EEMC_cut->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    for (int i = 1; i <= h_t_REC_EEMC_cut->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_EEMC_cut->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_EEMC_cut->SetBinError(i, newError);
    }
	h_t_REC_EEMC_cut->SetMarkerStyle(20); // method L RECO
	h_t_REC_EEMC_cut->SetMarkerColor(kP8Blue);
	//h_t_REC_EEMC_cut->Draw("PEsame");

    h_t_REC_wRES_cut_pi12->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_wRES_cut_pi12->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_t_REC_wRES_cut_pi12->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_wRES_cut_pi12->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_wRES_cut_pi12->SetBinError(i, newError);
    }
	h_t_REC_wRES_cut_pi12->SetMarkerStyle(30);
	h_t_REC_wRES_cut_pi12->SetMarkerColor(kP8Pink);
    h_t_REC_wRES_cut_pi12->SetLineColor(kP8Pink);
	h_t_REC_wRES_cut_pi12->Draw("PEsame");

    h_incoherent_phi_noVetoes_MC->Rebin(rebin_width);
    h_incoherent_phi_noVetoes_MC->Scale(sigma_incoherent*(1/nEvents_incoherent)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_incoherent_lumi = nEvents_incoherent/sigma_incoherent; // nb^-1 
    double ratio_incoherent_preTDR = preTDR_lumi/simu_incoherent_lumi;
    for (int i = 1; i <= h_incoherent_phi_noVetoes_MC->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_noVetoes_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_preTDR);
        h_incoherent_phi_noVetoes_MC->SetBinError(i, newError);
    }
    h_incoherent_phi_noVetoes_MC->SetMarkerStyle(21);
    h_incoherent_phi_noVetoes_MC->SetMarkerColor(kP8Orange);
    h_incoherent_phi_noVetoes_MC->SetLineColor(kP8Orange);
    //h_incoherent_phi_noVetoes_MC->Draw("PEsame");

    h_incoherent_phi_noVetoes_L->Rebin(rebin_width);
    h_incoherent_phi_noVetoes_L->Scale(sigma_incoherent*(1/nEvents_incoherent)*(1/dN_dt)/rebin_width);
    for (int i = 1; i <= h_incoherent_phi_noVetoes_L->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_noVetoes_L->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_preTDR);
        h_incoherent_phi_noVetoes_L->SetBinError(i, newError);
    }
    h_incoherent_phi_noVetoes_L->SetMarkerStyle(31);
    h_incoherent_phi_noVetoes_L->SetMarkerColor(kP8Blue);
    h_incoherent_phi_noVetoes_L->SetLineColor(kP8Blue);
    //h_incoherent_phi_noVetoes_L->Draw("PEsame");

    h_incoherent_phi_noVetoes_proj->Rebin(rebin_width);
    h_incoherent_phi_noVetoes_proj->Scale(sigma_incoherent*(1/nEvents_incoherent)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_incoherent_phi_noVetoes_proj->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_noVetoes_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_preTDR);
        h_incoherent_phi_noVetoes_proj->SetBinError(i, newError);
    }
    h_incoherent_phi_noVetoes_proj->SetMarkerStyle(21);
    h_incoherent_phi_noVetoes_proj->SetMarkerColor(kP8Orange);
    h_incoherent_phi_noVetoes_proj->SetLineColor(kP8Orange);
    h_incoherent_phi_noVetoes_proj->Draw("PEsame");

    h_incoherent_phi_noVetoes_proj_pi2->Rebin(rebin_width);
    h_incoherent_phi_noVetoes_proj_pi2->Scale(sigma_incoherent*(1/nEvents_incoherent)*(1/dN_dt)/rebin_width);
    for (int i = 1; i <= h_incoherent_phi_noVetoes_proj_pi2->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_noVetoes_proj_pi2->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_preTDR);
        h_incoherent_phi_noVetoes_proj_pi2->SetBinError(i, newError);
    }
    h_incoherent_phi_noVetoes_proj_pi2->SetMarkerStyle(42);
    h_incoherent_phi_noVetoes_proj_pi2->SetMarkerColor(kP8Gray);
    h_incoherent_phi_noVetoes_proj_pi2->SetLineColor(kP8Gray);
    //h_incoherent_phi_noVetoes_proj_pi2->Draw("PEsame");
    
	TLegend *w14_213 = new TLegend(0.6,0.7,0.72,0.9);
    w14_213->AddEntry(h_t_MC,"Coherent #phi: MC", "L");
    w14_213->AddEntry(h_t_REC_wRES_cut_pi12,"Coherent #phi: #theta_{max}= #pi/12 ", "P");    
    //w14_213->AddEntry(h_t_REC_EEMC_cut,"Coherent #phi: L", "P");
    //w14_213->AddEntry(h_incoherent_phi_noVetoes_MC,"Incoherent #phi: MC", "P");
    //w14_213->AddEntry(h_incoherent_phi_noVetoes_L,"Incoherent #phi: L", "P");
    w14_213->AddEntry(h_incoherent_phi_noVetoes_proj,"Incoherent #phi: #theta_{max}= #pi/12", "P");
    //w14_213->AddEntry(h_incoherent_phi_noVetoes_proj_pi2,"Incoherent #phi: #theta_{max}= #pi/2", "P");
	w14_213->SetBorderSize(0);   
    w14_213->SetFillStyle(0);
	w14_213->SetTextSize(30);
	w14_213->SetTextFont(45);
    w14_213->Draw("same");	

    c14_213->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title14_213 = new TLatex();
    title14_213->SetNDC(); 
    title14_213->SetTextSize(0.05);
    title14_213->SetTextAlign(22);  
    //title14_213->DrawLatex(0.5, 0.97, "|t| Distribution with Incoherent #phi No Vetoes");  

    TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

    // Add labels
    TLatex* ep = new TLatex(0.13, 0.85, "ePIC");
    ep->SetNDC();
    ep->SetTextFont(62);
    ep->SetTextSize(.038);
    ep->SetTextColor(kBlack);
    ep->Draw("same");

    TLatex* r421 = new TLatex(0.18, 0.85, " Simulation 25.10.2, 10x100 GeV");
	r421->SetNDC();
	r421->SetTextSize(30);
	r421->SetTextFont(43);
	r421->SetTextColor(kBlack);
	r421->Draw("same");

	TLatex* label_5 = new TLatex(0.17, 0.75, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->SetTextSize(30);
	label_5->SetTextFont(43);
	label_5->SetTextColor(kBlack);
	label_5->Draw("same");

    TLatex* r4211 = new TLatex(0.15, 0.8, "eAu #rightarrow e'Au'#phi");
	r4211->SetNDC();
	r4211->SetTextSize(30);
	r4211->SetTextFont(43);
	r4211->SetTextColor(kBlack);
	r4211->Draw("same"); 

    TLatex* r421111 = new TLatex(0.18, 0.15, "L_{int} = 10 fb^{-1}/A"); // preTDR lumi
	r421111->SetNDC();
	r421111->SetTextSize(25);
	r421111->SetTextFont(43);
	r421111->SetTextColor(kBlack);
	r421111->Draw("same");

    c14_213->Print("./figures/plot_t_incoherent_noVetoes_preTDR_proj.pdf");

	//TFile* fout_preTDR = new TFile("plot_t_incoherent_noVetoes_preTDR.root","RECREATE");
    //c14_213->Write();
    //h_t_MC->Write();
    //h_t_REC_EEMC_cut->Write();
    //h_t_REC_wRES_cut_pi12->Write();
    //h_incoherent_phi_noVetoes_MC->Write();
    //h_incoherent_phi_noVetoes_L->Write();
    //h_incoherent_phi_noVetoes_proj->Write();
    //fout_preTDR->Close();
}

void plot_incoherent_phi_noVetoes_ES()
{
	TFile* phi_t_file = TFile::Open("phi_sartre_25_10_2_eAu_10x100v5.root","READ");
	TH1D* h_t_MC = (TH1D*) phi_t_file->Get("h_t_MC");
	TH1D* h_t_REC_EEMC_cut = (TH1D*) phi_t_file->Get("h_t_REC_EEMC_cut");
	TH1D* h_t_REC_wRES_cut_pi12 = (TH1D*) phi_t_file->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_phi_sartre_events = (TH1D*) phi_t_file->Get("h_Nevents");

	TFile* incoherent_phi_noVetoes = TFile::Open("phi_beagle_25_10_3_eAu_noVetoes_10x100v4.root","READ");
    TH1D* h_incoherent_phi_noVetoes_MC = (TH1D*) incoherent_phi_noVetoes->Get("h_t_MC");
    TH1D* h_incoherent_phi_noVetoes_L = (TH1D*) incoherent_phi_noVetoes->Get("h_t_REC_EEMC_cut");
    TH1D* h_incoherent_phi_noVetoes_proj = (TH1D*) incoherent_phi_noVetoes->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_incoherent_phi_noVetoes_projNoScale = (TH1D*) incoherent_phi_noVetoes->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_incoherent_phi_noVetoes_proj_pi2 = (TH1D*) incoherent_phi_noVetoes->Get("h_t_REC_wRES_cut_pi2");
    TH1D* h_incoherent_events = (TH1D*) incoherent_phi_noVetoes->Get("h_Nevents");

	double dN_dt = 0.002;
    double sigma_coherent = 459.05; // nb
    double sigma_incoherent = 261.95; // nb
    double nEvents_coherent = h_phi_sartre_events->GetEntries();
    double nEvents_incoherent = h_incoherent_events->GetEntries();
    double early_science_lumi = 8121.83; // nb^-1
    int rebin_width = 2;   

	TCanvas* c14_213 = new TCanvas("c14_213","c14_213",1,1,1200,800);
    gPad->SetLogy(1);
	gStyle->SetOptStat(0);
    h_t_MC->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
	h_t_MC->GetXaxis()->SetTitle("|t| [GeV/c]^{2}");
    h_t_MC->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    double simu_coherent_lumi = nEvents_coherent/sigma_coherent; // nb^-1 
    double ratio_coherent_phi_ES = early_science_lumi/simu_coherent_lumi;
    for (int i = 1; i <= h_t_MC->GetNbinsX(); ++i) 
    {
        double binError = h_t_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_MC->SetBinError(i, newError);
    }
    h_t_MC->SetLineColor(kBlack);
	h_t_MC->Draw("HISTsame");

    h_t_REC_EEMC_cut->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    for (int i = 1; i <= h_t_REC_EEMC_cut->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_EEMC_cut->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_REC_EEMC_cut->SetBinError(i, newError);
    }
	h_t_REC_EEMC_cut->SetMarkerStyle(20); // method L RECO
	h_t_REC_EEMC_cut->SetMarkerColor(kP8Blue);
	//h_t_REC_EEMC_cut->Draw("PEsame");

    h_t_REC_wRES_cut_pi12->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_wRES_cut_pi12->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_t_REC_wRES_cut_pi12->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_wRES_cut_pi12->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_REC_wRES_cut_pi12->SetBinError(i, newError);
    }
	h_t_REC_wRES_cut_pi12->SetMarkerStyle(30);
	h_t_REC_wRES_cut_pi12->SetMarkerColor(kP8Pink);
    h_t_REC_wRES_cut_pi12->SetLineColor(kP8Pink);
	h_t_REC_wRES_cut_pi12->Draw("PEsame");

    h_incoherent_phi_noVetoes_MC->Rebin(rebin_width);
    h_incoherent_phi_noVetoes_MC->Scale(sigma_incoherent*(1/nEvents_incoherent)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_incoherent_lumi = nEvents_incoherent/sigma_incoherent; // nb^-1 
    double ratio_incoherent_ES = early_science_lumi/simu_incoherent_lumi;
    for (int i = 1; i <= h_incoherent_phi_noVetoes_MC->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_noVetoes_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_ES);
        h_incoherent_phi_noVetoes_MC->SetBinError(i, newError);
    }
    h_incoherent_phi_noVetoes_MC->SetMarkerStyle(21);
    h_incoherent_phi_noVetoes_MC->SetMarkerColor(kP8Orange);
    h_incoherent_phi_noVetoes_MC->SetLineColor(kP8Orange);
    //h_incoherent_phi_noVetoes_MC->Draw("PEsame");

    h_incoherent_phi_noVetoes_L->Rebin(rebin_width);
    h_incoherent_phi_noVetoes_L->Scale(sigma_incoherent*(1/nEvents_incoherent)*(1/dN_dt)/rebin_width);
    for (int i = 1; i <= h_incoherent_phi_noVetoes_L->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_noVetoes_L->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_ES);
        h_incoherent_phi_noVetoes_L->SetBinError(i, newError);
    }
    h_incoherent_phi_noVetoes_L->SetMarkerStyle(31);
    h_incoherent_phi_noVetoes_L->SetMarkerColor(kP8Blue);
    h_incoherent_phi_noVetoes_L->SetLineColor(kP8Blue);
    //h_incoherent_phi_noVetoes_L->Draw("PEsame");

    h_incoherent_phi_noVetoes_proj->Rebin(rebin_width);
    h_incoherent_phi_noVetoes_proj->Scale(sigma_incoherent*(1/nEvents_incoherent)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_incoherent_phi_noVetoes_proj->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_noVetoes_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_ES);
        h_incoherent_phi_noVetoes_proj->SetBinError(i, newError);
    }
    h_incoherent_phi_noVetoes_proj->SetMarkerStyle(21);
    h_incoherent_phi_noVetoes_proj->SetMarkerColor(kP8Orange);
    h_incoherent_phi_noVetoes_proj->SetLineColor(kP8Orange);
    h_incoherent_phi_noVetoes_proj->Draw("PEsame");

    h_incoherent_phi_noVetoes_proj_pi2->Rebin(rebin_width);
    h_incoherent_phi_noVetoes_proj_pi2->Scale(sigma_incoherent*(1/nEvents_incoherent)*(1/dN_dt)/rebin_width);
    for (int i = 1; i <= h_incoherent_phi_noVetoes_proj->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_noVetoes_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_ES);
        h_incoherent_phi_noVetoes_proj->SetBinError(i, newError);
    }
    h_incoherent_phi_noVetoes_proj_pi2->SetMarkerStyle(42);
    h_incoherent_phi_noVetoes_proj_pi2->SetMarkerColor(kP8Gray);
    h_incoherent_phi_noVetoes_proj_pi2->SetLineColor(kP8Gray);
    //h_incoherent_phi_noVetoes_proj_pi2->Draw("PEsame");
    
	TLegend *w14_213 = new TLegend(0.6,0.7,0.72,0.9);
    w14_213->AddEntry(h_t_MC,"Coherent #phi: MC", "L");
    w14_213->AddEntry(h_t_REC_wRES_cut_pi12,"Coherent #phi: #theta_{max}= #pi/12 ", "P");    
    //w14_213->AddEntry(h_t_REC_EEMC_cut,"Coherent #phi: L", "P");
    //w14_213->AddEntry(h_incoherent_phi_noVetoes_MC,"Incoherent #phi: MC", "P");
    //w14_213->AddEntry(h_incoherent_phi_noVetoes_L,"Incoherent #phi: L", "P");
    w14_213->AddEntry(h_incoherent_phi_noVetoes_proj,"Incoherent #phi: #theta_{max}= #pi/12", "P");
    //w14_213->AddEntry(h_incoherent_phi_noVetoes_proj_pi2,"Incoherent #phi: #theta_{max}= #pi/2", "P");
	w14_213->SetBorderSize(0);   
    w14_213->SetFillStyle(0);
	w14_213->SetTextSize(30);
	w14_213->SetTextFont(45);
    w14_213->Draw("same");	

    c14_213->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title14_213 = new TLatex();
    title14_213->SetNDC(); 
    title14_213->SetTextSize(0.05);
    title14_213->SetTextAlign(22);  
    //title14_213->DrawLatex(0.5, 0.97, "|t| Distribution with Incoherent #phi No Vetoes");  

    TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

    // Add labels
    TLatex* ep = new TLatex(0.13, 0.85, "ePIC");
    ep->SetNDC();
    ep->SetTextFont(62);
    ep->SetTextSize(.038);
    ep->SetTextColor(kBlack);
    ep->Draw("same");

    TLatex* r421 = new TLatex(0.18, 0.85, " Simulation 25.10.2, 10x100 GeV");
	r421->SetNDC();
	r421->SetTextSize(30);
	r421->SetTextFont(43);
	r421->SetTextColor(kBlack);
	r421->Draw("same");

	TLatex* label_5 = new TLatex(0.17, 0.75, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->SetTextSize(30);
	label_5->SetTextFont(43);
	label_5->SetTextColor(kBlack);
	label_5->Draw("same");

    TLatex* r4211 = new TLatex(0.15, 0.8, "eAu #rightarrow e'Au'#phi");
	r4211->SetNDC();
	r4211->SetTextSize(30);
	r4211->SetTextFont(43);
	r4211->SetTextColor(kBlack);
	r4211->Draw("same"); 

    TLatex* r42111 = new TLatex(0.18, 0.15, "L_{int} = 1.6 fb^{-1}/A"); // ES lumi
	r42111->SetNDC();
	r42111->SetTextSize(30);
	r42111->SetTextFont(43);
	r42111->SetTextColor(kBlack);
	r42111->Draw("same");

    c14_213->Print("./figures/plot_t_incoherent_noVetoes_ES_proj.pdf");

	//TFile* fout_ES = new TFile("plot_t_incoherent_noVetoes_ES.root","RECREATE");
    //c14_213->Write();
    //h_t_MC->Write();
    //h_t_REC_EEMC_cut->Write();
    //h_t_REC_wRES_cut_pi12->Write();
    //h_incoherent_phi_noVetoes_MC->Write();
    //h_incoherent_phi_noVetoes_L->Write();
    //h_incoherent_phi_noVetoes_proj->Write();
    //fout_ES->Close();
}

void plot_incoherent_phi_allDetectorVetoes_preTDR()
{
	TFile* phi_t_file = TFile::Open("phi_sartre_25_10_2_eAu_10x100v5.root","READ");
	TH1D* h_t_MC = (TH1D*) phi_t_file->Get("h_t_MC");
	TH1D* h_t_REC_EEMC_cut = (TH1D*) phi_t_file->Get("h_t_REC_EEMC_cut");
	TH1D* h_t_REC_wRES_cut_pi12 = (TH1D*) phi_t_file->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_phi_sartre_events = (TH1D*) phi_t_file->Get("h_Nevents");

	TFile* incoherent_phi_allDetectorVetoes = TFile::Open("phi_beagle_25_10_3_eAu_allDetectorVetoes_10x100v4.root");
    TH1D* h_incoherent_phi_onlyDetectors_MC = (TH1D*) incoherent_phi_allDetectorVetoes->Get("h_t_MC");
    TH1D* h_incoherent_phi_onlyDetectors_L = (TH1D*) incoherent_phi_allDetectorVetoes->Get("h_t_REC_EEMC_cut");
    TH1D* h_incoherent_phi_onlyDetectors_proj = (TH1D*) incoherent_phi_allDetectorVetoes->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_incoherent_events = (TH1D*) incoherent_phi_allDetectorVetoes->Get("h_Nevents");

	double dN_dt = 0.002;
    double sigma_coherent = 459.05; // nb
    double sigma_incoherent = 261.95; // nb
    double nEvents_coherent = h_phi_sartre_events->GetEntries();
    double nEvents_incoherent = h_incoherent_events->GetEntries();
    double preTDR_lumi = 50761.4; // nb^-1
    int rebin_width = 2;   

	TCanvas* c14_2132 = new TCanvas("c14_2132","c14_2132",1,1,1200,800);
    gPad->SetLogy(1);
    gStyle->SetOptStat(0);
    h_t_MC->GetYaxis()->SetTitle("dN/d|t| [GeV/c]^{-2}");
	h_t_MC->GetXaxis()->SetTitle("|t| [GeV/c]^{2}");
    h_t_MC->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    double simu_coherent_lumi = nEvents_coherent/sigma_coherent; // nb^-1 
    double ratio_coherent_phi_preTDR = preTDR_lumi/simu_coherent_lumi;
    for (int i = 1; i <= h_t_MC->GetNbinsX(); ++i) 
    {
        double binError = h_t_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_MC->SetBinError(i, newError);
    }
    h_t_MC->SetLineColor(kBlack);
	h_t_MC->Draw("HISTsame");

    h_t_REC_EEMC_cut->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    for (int i = 1; i <= h_t_REC_EEMC_cut->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_EEMC_cut->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_EEMC_cut->SetBinError(i, newError);
    }
	h_t_REC_EEMC_cut->SetMarkerStyle(20); // method L RECO
	h_t_REC_EEMC_cut->SetMarkerColor(kP8Blue);
	//h_t_REC_EEMC_cut->Draw("PEsame");

    h_t_REC_wRES_cut_pi12->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_wRES_cut_pi12->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_t_REC_wRES_cut_pi12->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_wRES_cut_pi12->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_wRES_cut_pi12->SetBinError(i, newError);
    }
	h_t_REC_wRES_cut_pi12->SetMarkerStyle(30);
	h_t_REC_wRES_cut_pi12->SetMarkerColor(kP8Pink);
    h_t_REC_wRES_cut_pi12->SetLineColor(kP8Pink);
	h_t_REC_wRES_cut_pi12->Draw("PEsame");

    h_incoherent_phi_onlyDetectors_MC->Rebin(rebin_width);
    h_incoherent_phi_onlyDetectors_MC->Scale(sigma_incoherent*(1/nEvents_incoherent)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_incoherent_lumi = nEvents_incoherent/sigma_incoherent; // nb^-1 
    double ratio_incoherent_preTDR = preTDR_lumi/simu_incoherent_lumi;
    for (int i = 1; i <= h_incoherent_phi_onlyDetectors_MC->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_onlyDetectors_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_preTDR);
        h_incoherent_phi_onlyDetectors_MC->SetBinError(i, newError);
    }
    h_incoherent_phi_onlyDetectors_MC->SetMarkerStyle(21);
    h_incoherent_phi_onlyDetectors_MC->SetMarkerColor(kP8Orange);
    h_incoherent_phi_onlyDetectors_MC->SetLineColor(kP8Orange);
    //h_incoherent_phi_onlyDetectors_MC->Draw("PEsame");

    h_incoherent_phi_onlyDetectors_L->Rebin(rebin_width);
    h_incoherent_phi_onlyDetectors_L->Scale(sigma_incoherent*(1/nEvents_incoherent)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_incoherent_phi_onlyDetectors_L->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_onlyDetectors_L->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_preTDR);
        h_incoherent_phi_onlyDetectors_L->SetBinError(i, newError);
    }
    h_incoherent_phi_onlyDetectors_L->SetMarkerStyle(21);
    h_incoherent_phi_onlyDetectors_L->SetMarkerColor(kP8Orange);
    h_incoherent_phi_onlyDetectors_L->SetLineColor(kP8Orange);
    //h_incoherent_phi_onlyDetectors_L->Draw("PEsame");

    h_incoherent_phi_onlyDetectors_proj->Rebin(rebin_width);
    h_incoherent_phi_onlyDetectors_proj->Scale(sigma_incoherent*(1/nEvents_incoherent)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_incoherent_phi_onlyDetectors_proj->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_onlyDetectors_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_preTDR);
        h_incoherent_phi_onlyDetectors_proj->SetBinError(i, newError);
    }
    h_incoherent_phi_onlyDetectors_proj->SetMarkerStyle(21);
    h_incoherent_phi_onlyDetectors_proj->SetMarkerColor(kP8Orange);
    h_incoherent_phi_onlyDetectors_proj->SetLineColor(kP8Orange);
    h_incoherent_phi_onlyDetectors_proj->Draw("PEsame");

	TLegend *w14_2132 = new TLegend(0.6,0.73,0.72,0.9);
    w14_2132->AddEntry(h_t_MC,"Coherent #phi: MC", "L");
    w14_2132->AddEntry(h_t_REC_wRES_cut_pi12,"Coherent #phi: #theta_{max}= #pi/12", "P");    
    //w14_2132->AddEntry(h_t_REC_EEMC_cut,"Coherent #phi: L", "P");
    //w14_2132->AddEntry(h_incoherent_phi_onlyDetectors_MC,"Incoherent #phi: MC", "P");
    //w14_2132->AddEntry(h_incoherent_phi_onlyDetectors_L,"Incoherent #phi: L", "P");
    w14_2132->AddEntry(h_incoherent_phi_onlyDetectors_proj,"Incoherent #phi: #theta_{max}= #pi/12", "P");
    w14_2132->SetBorderSize(0);   
    w14_2132->SetFillStyle(0);
    w14_2132->SetTextSize(30);
	w14_2132->SetTextFont(45);
    w14_2132->Draw("same");	
        
    c14_2132->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title14_2132 = new TLatex();
    title14_2132->SetNDC(); 
    title14_2132->SetTextSize(0.05);
    title14_2132->SetTextAlign(22);  
    //title14_2132->DrawLatex(0.5, 0.97, "|t| Distribution with Incoherent #phi Only Detector Vetoes");  

    TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

    // Add labels
    TLatex* ep = new TLatex(0.13, 0.85, "ePIC");
    ep->SetNDC();
    ep->SetTextFont(62);
    ep->SetTextSize(.038);
    ep->SetTextColor(kBlack);
    ep->Draw("same");

    TLatex* r421 = new TLatex(0.18, 0.85, " Simulation 25.10.2, 10x100 GeV");
	r421->SetNDC();
	r421->SetTextSize(30);
	r421->SetTextFont(43);
	r421->SetTextColor(kBlack);
	r421->Draw("same");

	TLatex* label_5 = new TLatex(0.17, 0.75, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->SetTextSize(30);
	label_5->SetTextFont(43);
	label_5->SetTextColor(kBlack);
	label_5->Draw("same");

    TLatex* r4211 = new TLatex(0.15, 0.8, "eAu #rightarrow e'Au'#phi");
	r4211->SetNDC();
	r4211->SetTextSize(30);
	r4211->SetTextFont(43);
	r4211->SetTextColor(kBlack);
	r4211->Draw("same"); 

    TLatex* r421111 = new TLatex(0.18, 0.13, "L_{int} = 10 fb^{-1}/A"); // preTDR lumi
	r421111->SetNDC();
	r421111->SetTextSize(25);
	r421111->SetTextFont(43);
	r421111->SetTextColor(kBlack);
	r421111->Draw("same");

    c14_2132->Print("./figures/plot_t_incoherent_detectorsOnly_preTDR_proj.pdf");

	//TFile* fout_preTDR = new TFile("plot_t_incoherent_detectorsOnly_MC_preTDR.root","RECREATE");
    //c14_2132->Write();
    //h_t_MC->Write();
    //h_t_REC_EEMC_cut->Write();
    //h_t_REC_wRES_cut_pi12->Write();
    //h_incoherent_phi_onlyDetectors_MC->Write();
    //h_incoherent_phi_onlyDetectors_L->Write();
    //h_incoherent_phi_onlyDetectors_proj->Write();
    //fout_preTDR->Close();
}

void plot_incoherent_phi_allDetectorVetoes_ES()
{
	TFile* phi_t_file = TFile::Open("phi_sartre_25_10_2_eAu_10x100v5.root","READ");
	TH1D* h_t_MC = (TH1D*) phi_t_file->Get("h_t_MC");
	TH1D* h_t_REC_EEMC_cut = (TH1D*) phi_t_file->Get("h_t_REC_EEMC_cut");
	TH1D* h_t_REC_wRES_cut_pi12 = (TH1D*) phi_t_file->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_phi_sartre_events = (TH1D*) phi_t_file->Get("h_Nevents");

	TFile* incoherent_phi_allDetectorVetoes = TFile::Open("phi_beagle_25_10_3_eAu_allDetectorVetoes_10x100v4.root");
    TH1D* h_incoherent_phi_onlyDetectors_MC = (TH1D*) incoherent_phi_allDetectorVetoes->Get("h_t_MC");
    TH1D* h_incoherent_phi_onlyDetectors_L = (TH1D*) incoherent_phi_allDetectorVetoes->Get("h_t_REC_EEMC_cut");
    TH1D* h_incoherent_phi_onlyDetectors_proj = (TH1D*) incoherent_phi_allDetectorVetoes->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_incoherent_events = (TH1D*) incoherent_phi_allDetectorVetoes->Get("h_Nevents");

	double dN_dt = 0.002;
    double sigma_coherent = 459.05; // nb
    double sigma_incoherent = 261.95; // nb
    double nEvents_coherent = h_phi_sartre_events->GetEntries();
    double nEvents_incoherent = h_incoherent_events->GetEntries();
    double early_science_lumi = 8121.83; // nb^-1
    int rebin_width = 2;   

	TCanvas* c14_2132 = new TCanvas("c14_2132","c14_2132",1,1,1200,800);
    gPad->SetLogy(1);
    gStyle->SetOptStat(0);
    h_t_MC->GetYaxis()->SetTitle("dN/d|t| [GeV/c]^{-2}");
	h_t_MC->GetXaxis()->SetTitle("|t| [GeV/c]^{2}");
    h_t_MC->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    double simu_coherent_lumi = nEvents_coherent/sigma_coherent; // nb^-1 
    double ratio_coherent_phi_ES = early_science_lumi/simu_coherent_lumi;
    for (int i = 1; i <= h_t_MC->GetNbinsX(); ++i) 
    {
        double binError = h_t_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_MC->SetBinError(i, newError);
    }
    h_t_MC->SetLineColor(kBlack);
	h_t_MC->Draw("HISTsame");

    h_t_REC_EEMC_cut->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    for (int i = 1; i <= h_t_REC_EEMC_cut->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_EEMC_cut->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_REC_EEMC_cut->SetBinError(i, newError);
    }
	h_t_REC_EEMC_cut->SetMarkerStyle(20); // method L RECO
	h_t_REC_EEMC_cut->SetMarkerColor(kP8Blue);
	//h_t_REC_EEMC_cut->Draw("PEsame");

    h_t_REC_wRES_cut_pi12->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_wRES_cut_pi12->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_t_REC_wRES_cut_pi12->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_wRES_cut_pi12->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_REC_wRES_cut_pi12->SetBinError(i, newError);
    }
	h_t_REC_wRES_cut_pi12->SetMarkerStyle(30);
	h_t_REC_wRES_cut_pi12->SetMarkerColor(kP8Pink);
    h_t_REC_wRES_cut_pi12->SetLineColor(kP8Pink);
	h_t_REC_wRES_cut_pi12->Draw("PEsame");

    h_incoherent_phi_onlyDetectors_MC->Rebin(rebin_width);
    h_incoherent_phi_onlyDetectors_MC->Scale(sigma_incoherent*(1/nEvents_incoherent)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_incoherent_lumi = nEvents_incoherent/sigma_incoherent; // nb^-1 
    double ratio_incoherent_ES = early_science_lumi/simu_incoherent_lumi;
    for (int i = 1; i <= h_incoherent_phi_onlyDetectors_MC->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_onlyDetectors_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_ES);
        h_incoherent_phi_onlyDetectors_MC->SetBinError(i, newError);
    }
    h_incoherent_phi_onlyDetectors_MC->SetMarkerStyle(21);
    h_incoherent_phi_onlyDetectors_MC->SetMarkerColor(kP8Orange);
    h_incoherent_phi_onlyDetectors_MC->SetLineColor(kP8Orange);
    //h_incoherent_phi_onlyDetectors_MC->Draw("PEsame");

    h_incoherent_phi_onlyDetectors_L->Rebin(rebin_width);
    h_incoherent_phi_onlyDetectors_L->Scale(sigma_incoherent*(1/nEvents_incoherent)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_incoherent_phi_onlyDetectors_L->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_onlyDetectors_L->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_ES);
        h_incoherent_phi_onlyDetectors_L->SetBinError(i, newError);
    }
    h_incoherent_phi_onlyDetectors_L->SetMarkerStyle(21);
    h_incoherent_phi_onlyDetectors_L->SetMarkerColor(kP8Orange);
    h_incoherent_phi_onlyDetectors_L->SetLineColor(kP8Orange);
    //h_incoherent_phi_onlyDetectors_L->Draw("PEsame");

    h_incoherent_phi_onlyDetectors_proj->Rebin(rebin_width);
    h_incoherent_phi_onlyDetectors_proj->Scale(sigma_incoherent*(1/nEvents_incoherent)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_incoherent_phi_onlyDetectors_proj->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_onlyDetectors_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_ES);
        h_incoherent_phi_onlyDetectors_proj->SetBinError(i, newError);
    }
    h_incoherent_phi_onlyDetectors_proj->SetMarkerStyle(21);
    h_incoherent_phi_onlyDetectors_proj->SetMarkerColor(kP8Orange);
    h_incoherent_phi_onlyDetectors_proj->SetLineColor(kP8Orange);
    h_incoherent_phi_onlyDetectors_proj->Draw("PEsame");

	TLegend *w14_2132 = new TLegend(0.6,0.73,0.72,0.9);
    w14_2132->AddEntry(h_t_MC,"Coherent #phi: MC", "L");
    w14_2132->AddEntry(h_t_REC_wRES_cut_pi12,"Coherent #phi: #theta_{max}= #pi/12", "P");    
    //w14_2132->AddEntry(h_t_REC_EEMC_cut,"Coherent #phi: L", "P");
    //w14_2132->AddEntry(h_incoherent_phi_onlyDetectors_MC,"Incoherent #phi: MC", "P");
    //w14_2132->AddEntry(h_incoherent_phi_onlyDetectors_L,"Incoherent #phi: L", "P");
    w14_2132->AddEntry(h_incoherent_phi_onlyDetectors_proj,"Incoherent #phi: #theta_{max}= #pi/12", "P");
    w14_2132->SetBorderSize(0);   
    w14_2132->SetFillStyle(0);
    w14_2132->SetTextSize(30);
	w14_2132->SetTextFont(45);
    w14_2132->Draw("same");	
        
    c14_2132->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title14_2132 = new TLatex();
    title14_2132->SetNDC(); 
    title14_2132->SetTextSize(0.05);
    title14_2132->SetTextAlign(22);  
    //title14_2132->DrawLatex(0.5, 0.97, "|t| Distribution with Incoherent #phi Only Detector Vetoes");  

    TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

    // Add labels
    TLatex* ep = new TLatex(0.13, 0.85, "ePIC");
    ep->SetNDC();
    ep->SetTextFont(62);
    ep->SetTextSize(.038);
    ep->SetTextColor(kBlack);
    ep->Draw("same");

    TLatex* r421 = new TLatex(0.18, 0.85, " Simulation 25.10.2, 10x100 GeV");
	r421->SetNDC();
	r421->SetTextSize(30);
	r421->SetTextFont(43);
	r421->SetTextColor(kBlack);
	r421->Draw("same");

	TLatex* label_5 = new TLatex(0.17, 0.75, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->SetTextSize(30);
	label_5->SetTextFont(43);
	label_5->SetTextColor(kBlack);
	label_5->Draw("same");

    TLatex* r4211 = new TLatex(0.15, 0.8, "eAu #rightarrow e'Au'#phi");
	r4211->SetNDC();
	r4211->SetTextSize(30);
	r4211->SetTextFont(43);
	r4211->SetTextColor(kBlack);
	r4211->Draw("same"); 

    TLatex* r42111 = new TLatex(0.18, 0.13, "L_{int} = 1.6 fb^{-1}/A"); // ES lumi
	r42111->SetNDC();
	r42111->SetTextSize(30);
	r42111->SetTextFont(43);
	r42111->SetTextColor(kBlack);
	r42111->Draw("same");

    c14_2132->Print("./figures/plot_t_incoherent_detectorsOnly_ES_proj.pdf");

	//TFile* fout_ES = new TFile("plot_t_incoherent_detectorsOnly_MC_ES.root","RECREATE");
    //c14_2132->Write();
    //h_t_MC->Write();
    //h_t_REC_EEMC_cut->Write();
    //h_t_REC_wRES_cut_pi12->Write();
    //h_incoherent_phi_onlyDetectors_MC->Write();
    //h_incoherent_phi_onlyDetectors_L->Write();
    //h_incoherent_phi_onlyDetectors_proj->Write();
    //fout_ES->Close();
}

void plot_incoherent_phi_etaVetoes_preTDR()
{
	TFile* phi_t_file = TFile::Open("phi_sartre_25_10_2_eAu_10x100v5.root","READ");
	TH1D* h_t_MC = (TH1D*) phi_t_file->Get("h_t_MC");
	TH1D* h_t_REC_EEMC_cut = (TH1D*) phi_t_file->Get("h_t_REC_EEMC_cut");
	TH1D* h_t_REC_wRES_cut_pi12 = (TH1D*) phi_t_file->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_phi_sartre_events = (TH1D*) phi_t_file->Get("h_Nevents");

	TFile* incoherent_phi_etaVetoes = TFile::Open("phi_beagle_25_10_3_eAu_etaVetoes_10x100v4.root");
    TH1D* h_incoherent_phi_etaOnly_MC = (TH1D*) incoherent_phi_etaVetoes->Get("h_t_MC");
    TH1D* h_incoherent_phi_etaOnly_L = (TH1D*) incoherent_phi_etaVetoes->Get("h_t_REC_EEMC_cut");
    TH1D* h_incoherent_phi_etaOnly_proj = (TH1D*) incoherent_phi_etaVetoes->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_incoherent_events = (TH1D*) incoherent_phi_etaVetoes->Get("h_Nevents");

	double dN_dt = 0.002;
    double sigma_coherent = 459.05; // nb
    double sigma_incoherent = 261.95; // nb
    double nEvents_coherent = h_phi_sartre_events->GetEntries();
    double nEvents_incoherent = h_incoherent_events->GetEntries();
    double preTDR_lumi = 50761.4; // nb^-1
    int rebin_width = 2;   

	TCanvas* c14_2133 = new TCanvas("c14_2133","c14_2133",1,1,1200,800);
    gPad->SetLogy(1);
    gStyle->SetOptStat(0);
    h_t_MC->GetYaxis()->SetTitle("dN/d|t| [GeV/c]^{-2}");
	h_t_MC->GetXaxis()->SetTitle("|t| [GeV/c]^{2}");
    h_t_MC->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    double simu_coherent_lumi = nEvents_coherent/sigma_coherent; // nb^-1 
    double ratio_coherent_phi_preTDR = preTDR_lumi/simu_coherent_lumi;
    for (int i = 1; i <= h_t_MC->GetNbinsX(); ++i) 
    {
        double binError = h_t_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_MC->SetBinError(i, newError);
    }
    h_t_MC->SetLineColor(kBlack);
	h_t_MC->Draw("HISTsame");

    h_t_REC_EEMC_cut->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    for (int i = 1; i <= h_t_REC_EEMC_cut->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_EEMC_cut->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_EEMC_cut->SetBinError(i, newError);
    }
	h_t_REC_EEMC_cut->SetMarkerStyle(20); // method L RECO
	h_t_REC_EEMC_cut->SetMarkerColor(kP8Blue);
	//h_t_REC_EEMC_cut->Draw("PEsame");

    h_t_REC_wRES_cut_pi12->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_wRES_cut_pi12->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_t_REC_wRES_cut_pi12->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_wRES_cut_pi12->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_wRES_cut_pi12->SetBinError(i, newError);
    }
	h_t_REC_wRES_cut_pi12->SetMarkerStyle(30);
	h_t_REC_wRES_cut_pi12->SetMarkerColor(kP8Pink);
    h_t_REC_wRES_cut_pi12->SetLineColor(kP8Pink);
	h_t_REC_wRES_cut_pi12->Draw("PEsame");

    h_incoherent_phi_etaOnly_MC->Rebin(rebin_width);
    h_incoherent_phi_etaOnly_MC->Scale(sigma_incoherent*(1/nEvents_incoherent)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_incoherent_lumi = nEvents_incoherent/sigma_incoherent; // nb^-1 
    double ratio_incoherent_preTDR = preTDR_lumi/simu_incoherent_lumi;
    for (int i = 1; i <= h_incoherent_phi_etaOnly_MC->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_etaOnly_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_preTDR);
        h_incoherent_phi_etaOnly_MC->SetBinError(i, newError);
    }
    h_incoherent_phi_etaOnly_MC->SetMarkerStyle(21);
    h_incoherent_phi_etaOnly_MC->SetMarkerColor(kP8Orange);
    h_incoherent_phi_etaOnly_MC->SetLineColor(kP8Orange);
    //h_incoherent_phi_etaOnly_MC->Draw("PEsame");

    h_incoherent_phi_etaOnly_L->Rebin(rebin_width);
    h_incoherent_phi_etaOnly_L->Scale(sigma_incoherent*(1/nEvents_incoherent)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_incoherent_phi_etaOnly_L->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_etaOnly_L->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_preTDR);
        h_incoherent_phi_etaOnly_L->SetBinError(i, newError);
    }
    h_incoherent_phi_etaOnly_L->SetMarkerStyle(21);
    h_incoherent_phi_etaOnly_L->SetMarkerColor(kP8Orange);
    h_incoherent_phi_etaOnly_L->SetLineColor(kP8Orange);
    //h_incoherent_phi_etaOnly_L->Draw("PEsame");

    h_incoherent_phi_etaOnly_proj->Rebin(rebin_width);
    h_incoherent_phi_etaOnly_proj->Scale(sigma_incoherent*(1/nEvents_incoherent)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_incoherent_phi_etaOnly_proj->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_etaOnly_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_preTDR);
        h_incoherent_phi_etaOnly_proj->SetBinError(i, newError);
    }
    h_incoherent_phi_etaOnly_proj->SetMarkerStyle(21);
    h_incoherent_phi_etaOnly_proj->SetMarkerColor(kP8Orange);
    h_incoherent_phi_etaOnly_proj->SetLineColor(kP8Orange);
    h_incoherent_phi_etaOnly_proj->Draw("PEsame");
    
	TLegend *w14_2133 = new TLegend(0.6,0.7,0.72,0.9);
    w14_2133->AddEntry(h_t_MC,"Coherent #phi: MC", "L");
    w14_2133->AddEntry(h_t_REC_wRES_cut_pi12,"Coherent #phi: #theta_{max}= #pi/12", "P");    
    //w14_2133->AddEntry(h_t_REC_EEMC_cut,"Coherent #phi: L", "P");
    //w14_2133->AddEntry(h_incoherent_phi_etaOnly_MC,"Incoherent #phi: MC", "P");
    //w14_2133->AddEntry(h_incoherent_phi_etaOnly_L,"Incoherent #phi: L", "P");
    w14_2133->AddEntry(h_incoherent_phi_etaOnly_proj,"Incoherent #phi: #theta_{max}= #pi/12", "P");
    w14_2133->SetBorderSize(0);   
    w14_2133->SetFillStyle(0);
    w14_2133->SetTextSize(30);
	w14_2133->SetTextFont(45);
    w14_2133->Draw("same");	

    c14_2133->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title14_2133 = new TLatex();
    title14_2133->SetNDC(); 
    title14_2133->SetTextSize(0.05);
    title14_2133->SetTextAlign(22);  
    //title14_2133->DrawLatex(0.5, 0.97, "|t| Distribution with Incoherent #phi Only #eta Vetoes");  

    TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

    // Add labels
    TLatex* ep = new TLatex(0.13, 0.85, "ePIC");
    ep->SetNDC();
    ep->SetTextFont(62);
    ep->SetTextSize(.038);
    ep->SetTextColor(kBlack);
    ep->Draw("same");

    TLatex* r421 = new TLatex(0.18, 0.85, " Simulation 25.10.2, 10x100 GeV");
	r421->SetNDC();
	r421->SetTextSize(30);
	r421->SetTextFont(43);
	r421->SetTextColor(kBlack);
	r421->Draw("same");

	TLatex* label_5 = new TLatex(0.17, 0.75, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->SetTextSize(30);
	label_5->SetTextFont(43);
	label_5->SetTextColor(kBlack);
	label_5->Draw("same");

    TLatex* r4211 = new TLatex(0.15, 0.8, "eAu #rightarrow e'Au'#phi");
	r4211->SetNDC();
	r4211->SetTextSize(30);
	r4211->SetTextFont(43);
	r4211->SetTextColor(kBlack);
	r4211->Draw("same"); 

    TLatex* r421111 = new TLatex(0.18, 0.15, "L_{int} = 10 fb^{-1}/A"); // preTDR lumi
	r421111->SetNDC();
	r421111->SetTextSize(25);
	r421111->SetTextFont(43);
	r421111->SetTextColor(kBlack);
	r421111->Draw("same");

    c14_2133->Print("./figures/plot_t_incoherent_etaOnly_preTDR_proj.pdf");

	//TFile* fout_preTDR = new TFile("plot_t_incoherent_etaOnly_L_preTDR.root","RECREATE");
    //c14_2133->Write();
    //h_t_MC->Write();
    //h_t_REC_EEMC_cut->Write();
    //h_t_REC_wRES_cut_pi12->Write();
    //h_incoherent_phi_etaOnly_MC->Write();
    //h_incoherent_phi_etaOnly_L->Write();
    //h_incoherent_phi_etaOnly_proj->Write();
    //fout_preTDR->Close();
}

void plot_incoherent_phi_etaVetoes_ES()
{
	TFile* phi_t_file = TFile::Open("phi_sartre_25_10_2_eAu_10x100v5.root","READ");
	TH1D* h_t_MC = (TH1D*) phi_t_file->Get("h_t_MC");
	TH1D* h_t_REC_EEMC_cut = (TH1D*) phi_t_file->Get("h_t_REC_EEMC_cut");
	TH1D* h_t_REC_wRES_cut_pi12 = (TH1D*) phi_t_file->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_phi_sartre_events = (TH1D*) phi_t_file->Get("h_Nevents");

	TFile* incoherent_phi_etaVetoes = TFile::Open("phi_beagle_25_10_3_eAu_etaVetoes_10x100v4.root");
    TH1D* h_incoherent_phi_etaOnly_MC = (TH1D*) incoherent_phi_etaVetoes->Get("h_t_MC");
    TH1D* h_incoherent_phi_etaOnly_L = (TH1D*) incoherent_phi_etaVetoes->Get("h_t_REC_EEMC_cut");
    TH1D* h_incoherent_phi_etaOnly_proj = (TH1D*) incoherent_phi_etaVetoes->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_incoherent_events = (TH1D*) incoherent_phi_etaVetoes->Get("h_Nevents");

	double dN_dt = 0.002;
    double sigma_coherent = 459.05; // nb
    double sigma_incoherent = 261.95; // nb
    double nEvents_coherent = h_phi_sartre_events->GetEntries();
    double nEvents_incoherent = h_incoherent_events->GetEntries();
    double early_science_lumi = 8121.83; // nb^-1
    int rebin_width = 2;   

	TCanvas* c14_2133 = new TCanvas("c14_2133","c14_2133",1,1,1200,800);
    c14_2133->Divide(1,1,0.01,0.01);
    c14_2133->cd(1);
    gPad->SetLogy(1);
    gStyle->SetOptStat(0);
    h_t_MC->GetYaxis()->SetTitle("dN/d|t| [GeV/c]^{-2}");
	h_t_MC->GetXaxis()->SetTitle("|t| [GeV/c]^{2}");
    h_t_MC->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    double simu_coherent_lumi = nEvents_coherent/sigma_coherent; // nb^-1 
    double ratio_coherent_phi_ES = early_science_lumi/simu_coherent_lumi;
    for (int i = 1; i <= h_t_MC->GetNbinsX(); ++i) 
    {
        double binError = h_t_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_MC->SetBinError(i, newError);
    }
    h_t_MC->SetLineColor(kBlack);
	h_t_MC->Draw("HISTsame");

    h_t_REC_EEMC_cut->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    for (int i = 1; i <= h_t_REC_EEMC_cut->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_EEMC_cut->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_REC_EEMC_cut->SetBinError(i, newError);
    }
	h_t_REC_EEMC_cut->SetMarkerStyle(20); // method L RECO
	h_t_REC_EEMC_cut->SetMarkerColor(kP8Blue);
	//h_t_REC_EEMC_cut->Draw("PEsame");

    h_t_REC_wRES_cut_pi12->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_wRES_cut_pi12->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_t_REC_wRES_cut_pi12->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_wRES_cut_pi12->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_REC_wRES_cut_pi12->SetBinError(i, newError);
    }
	h_t_REC_wRES_cut_pi12->SetMarkerStyle(30);
	h_t_REC_wRES_cut_pi12->SetMarkerColor(kP8Pink);
    h_t_REC_wRES_cut_pi12->SetLineColor(kP8Pink);
	h_t_REC_wRES_cut_pi12->Draw("PEsame");

    h_incoherent_phi_etaOnly_MC->Rebin(rebin_width);
    h_incoherent_phi_etaOnly_MC->Scale(sigma_incoherent*(1/nEvents_incoherent)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_incoherent_lumi = nEvents_incoherent/sigma_incoherent; // nb^-1 
    double ratio_incoherent_ES = early_science_lumi/simu_incoherent_lumi;
    for (int i = 1; i <= h_incoherent_phi_etaOnly_MC->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_etaOnly_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_ES);
        h_incoherent_phi_etaOnly_MC->SetBinError(i, newError);
    }
    h_incoherent_phi_etaOnly_MC->SetMarkerStyle(21);
    h_incoherent_phi_etaOnly_MC->SetMarkerColor(kP8Orange);
    h_incoherent_phi_etaOnly_MC->SetLineColor(kP8Orange);
    //h_incoherent_phi_etaOnly_MC->Draw("PEsame");

    h_incoherent_phi_etaOnly_L->Rebin(rebin_width);
    h_incoherent_phi_etaOnly_L->Scale(sigma_incoherent*(1/nEvents_incoherent)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_incoherent_phi_etaOnly_L->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_etaOnly_L->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_ES);
        h_incoherent_phi_etaOnly_L->SetBinError(i, newError);
    }
    h_incoherent_phi_etaOnly_L->SetMarkerStyle(21);
    h_incoherent_phi_etaOnly_L->SetMarkerColor(kP8Orange);
    h_incoherent_phi_etaOnly_L->SetLineColor(kP8Orange);
    //h_incoherent_phi_etaOnly_L->Draw("PEsame");

    h_incoherent_phi_etaOnly_proj->Rebin(rebin_width);
    h_incoherent_phi_etaOnly_proj->Scale(sigma_incoherent*(1/nEvents_incoherent)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_incoherent_phi_etaOnly_proj->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_etaOnly_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_ES);
        h_incoherent_phi_etaOnly_proj->SetBinError(i, newError);
    }
    h_incoherent_phi_etaOnly_proj->SetMarkerStyle(21);
    h_incoherent_phi_etaOnly_proj->SetMarkerColor(kP8Orange);
    h_incoherent_phi_etaOnly_proj->SetLineColor(kP8Orange);
    h_incoherent_phi_etaOnly_proj->Draw("PEsame");
    
	TLegend *w14_2133 = new TLegend(0.6,0.7,0.72,0.9);
    w14_2133->AddEntry(h_t_MC,"Coherent #phi: MC", "L");
    w14_2133->AddEntry(h_t_REC_wRES_cut_pi12,"Coherent #phi: #theta_{max}= #pi/12", "P");    
    //w14_2133->AddEntry(h_t_REC_EEMC_cut,"Coherent #phi: L", "P");
    //w14_2133->AddEntry(h_incoherent_phi_etaOnly_MC,"Incoherent #phi: MC", "P");
    //w14_2133->AddEntry(h_incoherent_phi_etaOnly_L,"Incoherent #phi: L", "P");
    w14_2133->AddEntry(h_incoherent_phi_etaOnly_proj,"Incoherent #phi: #theta_{max}= #pi/12", "P");
    w14_2133->SetBorderSize(0);   
    w14_2133->SetFillStyle(0);
    w14_2133->SetTextSize(30);
	w14_2133->SetTextFont(45);
    w14_2133->Draw("same");	

    c14_2133->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title14_2133 = new TLatex();
    title14_2133->SetNDC(); 
    title14_2133->SetTextSize(0.05);
    title14_2133->SetTextAlign(22);  
    //title14_2133->DrawLatex(0.5, 0.97, "|t| Distribution with Incoherent #phi Only #eta Vetoes");  

    TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

    // Add labels
    TLatex* ep = new TLatex(0.13, 0.85, "ePIC");
    ep->SetNDC();
    ep->SetTextFont(62);
    ep->SetTextSize(.038);
    ep->SetTextColor(kBlack);
    ep->Draw("same");

    TLatex* r421 = new TLatex(0.18, 0.85, " Simulation 25.10.2, 10x100 GeV");
	r421->SetNDC();
	r421->SetTextSize(30);
	r421->SetTextFont(43);
	r421->SetTextColor(kBlack);
	r421->Draw("same");

	TLatex* label_5 = new TLatex(0.17, 0.75, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->SetTextSize(30);
	label_5->SetTextFont(43);
	label_5->SetTextColor(kBlack);
	label_5->Draw("same");

    TLatex* r4211 = new TLatex(0.15, 0.8, "eAu #rightarrow e'Au'#phi");
	r4211->SetNDC();
	r4211->SetTextSize(30);
	r4211->SetTextFont(43);
	r4211->SetTextColor(kBlack);
	r4211->Draw("same"); 

    TLatex* r42111 = new TLatex(0.18, 0.15, "L_{int} = 1.6 fb^{-1}/A"); // ES lumi
	r42111->SetNDC();
	r42111->SetTextSize(30);
	r42111->SetTextFont(43);
	r42111->SetTextColor(kBlack);
	r42111->Draw("same");

    c14_2133->Print("./figures/plot_t_incoherent_etaOnly_ES_proj.pdf");

	//TFile* fout_ES = new TFile("plot_t_incoherent_etaOnly_L_ES.root","RECREATE");
    //c14_2133->Write();
    //h_t_MC->Write();
    //h_t_REC_EEMC_cut->Write();
    //h_t_REC_wRES_cut_pi12->Write();
    //h_incoherent_phi_etaOnly_MC->Write();
    //h_incoherent_phi_etaOnly_L->Write();
    //h_incoherent_phi_etaOnly_proj->Write();
    //fout_ES->Close();
}

void plot_incoherent_phi_rpVetoes_preTDR()
{
	TFile* phi_t_file = TFile::Open("phi_sartre_25_10_2_eAu_10x100v5.root","READ");
	TH1D* h_t_MC = (TH1D*) phi_t_file->Get("h_t_MC");
	TH1D* h_t_REC_EEMC_cut = (TH1D*) phi_t_file->Get("h_t_REC_EEMC_cut");
	TH1D* h_t_REC_wRES_cut_pi12 = (TH1D*) phi_t_file->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_phi_sartre_events = (TH1D*) phi_t_file->Get("h_Nevents");

	TFile* incoherent_phi_rpVetoes = TFile::Open("phi_beagle_25_10_3_eAu_rpVetoes_10x100v4.root");
    TH1D* h_incoherent_phi_RPonly_MC = (TH1D*) incoherent_phi_rpVetoes->Get("h_t_MC");
    TH1D* h_incoherent_phi_RPonly_L = (TH1D*) incoherent_phi_rpVetoes->Get("h_t_REC_EEMC_cut");
    TH1D* h_incoherent_phi_RPonly_proj = (TH1D*) incoherent_phi_rpVetoes->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_incoherent_events = (TH1D*) incoherent_phi_rpVetoes->Get("h_Nevents");

	double dN_dt = 0.002;
    double sigma_coherent = 459.05; // nb
    double sigma_incoherent = 261.95; // nb
    double nEvents_coherent = h_phi_sartre_events->GetEntries();
    double nEvents_incoherent = h_incoherent_events->GetEntries();
    double preTDR_lumi = 50761.4; // nb^-1
    int rebin_width = 2;   

	TCanvas* c14_2135 = new TCanvas("c14_2135","c14_2135",1,1,1200,800);
    gPad->SetLogy(1);
    gStyle->SetOptStat(0);
    h_t_MC->GetYaxis()->SetTitle("dN/d|t| [GeV/c]^{-2}");
	h_t_MC->GetXaxis()->SetTitle("|t| [GeV/c]^{2}");
    h_t_MC->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    double simu_coherent_lumi = nEvents_coherent/sigma_coherent; // nb^-1 
    double ratio_coherent_phi_preTDR = preTDR_lumi/simu_coherent_lumi;
    for (int i = 1; i <= h_t_MC->GetNbinsX(); ++i) 
    {
        double binError = h_t_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_MC->SetBinError(i, newError);
    }
    h_t_MC->SetLineColor(kBlack);
	h_t_MC->Draw("HISTsame");

    h_t_REC_EEMC_cut->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    for (int i = 1; i <= h_t_REC_EEMC_cut->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_EEMC_cut->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_EEMC_cut->SetBinError(i, newError);
    }
	h_t_REC_EEMC_cut->SetMarkerStyle(20); // method L RECO
	h_t_REC_EEMC_cut->SetMarkerColor(kP8Blue);
	//h_t_REC_EEMC_cut->Draw("PEsame");

    h_t_REC_wRES_cut_pi12->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_wRES_cut_pi12->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_t_REC_wRES_cut_pi12->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_wRES_cut_pi12->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_wRES_cut_pi12->SetBinError(i, newError);
    }
	h_t_REC_wRES_cut_pi12->SetMarkerStyle(30);
	h_t_REC_wRES_cut_pi12->SetMarkerColor(kP8Pink);
    h_t_REC_wRES_cut_pi12->SetLineColor(kP8Pink);
	h_t_REC_wRES_cut_pi12->Draw("PEsame");

    h_incoherent_phi_RPonly_MC->Rebin(rebin_width);
    h_incoherent_phi_RPonly_MC->Scale(sigma_incoherent*(1/nEvents_incoherent)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_incoherent_lumi = nEvents_incoherent/sigma_incoherent; // nb^-1 
    double ratio_incoherent_preTDR = preTDR_lumi/simu_incoherent_lumi;
    for (int i = 1; i <= h_incoherent_phi_RPonly_MC->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_RPonly_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_preTDR);
        h_incoherent_phi_RPonly_MC->SetBinError(i, newError);
    }
    h_incoherent_phi_RPonly_MC->SetMarkerStyle(21);
    h_incoherent_phi_RPonly_MC->SetMarkerColor(kP8Orange);
    h_incoherent_phi_RPonly_MC->SetLineColor(kP8Orange);
    //h_incoherent_phi_RPonly_MC->Draw("PEsame");

    h_incoherent_phi_RPonly_L->Rebin(rebin_width);
    h_incoherent_phi_RPonly_L->Scale(sigma_incoherent*(1/nEvents_incoherent)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_incoherent_phi_RPonly_L->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_RPonly_L->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_preTDR);
        h_incoherent_phi_RPonly_L->SetBinError(i, newError);
    }
    h_incoherent_phi_RPonly_L->SetMarkerStyle(21);
    h_incoherent_phi_RPonly_L->SetMarkerColor(kP8Orange);
    h_incoherent_phi_RPonly_L->SetLineColor(kP8Orange);
    //h_incoherent_phi_RPonly_L->Draw("PEsame");

    h_incoherent_phi_RPonly_proj->Rebin(rebin_width);
    h_incoherent_phi_RPonly_proj->Scale(sigma_incoherent*(1/nEvents_incoherent)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_incoherent_phi_RPonly_proj->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_RPonly_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_preTDR);
        h_incoherent_phi_RPonly_proj->SetBinError(i, newError);
    }
    h_incoherent_phi_RPonly_proj->SetMarkerStyle(21);
    h_incoherent_phi_RPonly_proj->SetMarkerColor(kP8Orange);
    h_incoherent_phi_RPonly_proj->SetLineColor(kP8Orange);
    h_incoherent_phi_RPonly_proj->Draw("PEsame");
    
	TLegend *w14_2135 = new TLegend(0.6,0.7,0.72,0.9);
	w14_2135->AddEntry(h_t_MC,"Coherent #phi: MC", "L");
    w14_2135->AddEntry(h_t_REC_wRES_cut_pi12,"Coherent #phi: #theta_{max}= #pi/12", "P");    
    //w14_2135->AddEntry(h_t_REC_EEMC_cut,"Coherent #phi: L", "P");
    //w14_2135->AddEntry(h_incoherent_phi_RPonly_MC,"Incoherent #phi: MC", "P");
    //w14_2135->AddEntry(h_incoherent_phi_RPonly_L,"Incoherent #phi: L", "P");
    w14_2135->AddEntry(h_incoherent_phi_RPonly_proj,"Incoherent #phi: #theta_{max}= #pi/12", "P");
    w14_2135->SetBorderSize(0);   
    w14_2135->SetFillStyle(0);
    w14_2135->SetTextSize(30);
	w14_2135->SetTextFont(45);
    w14_2135->Draw("same");	

    c14_2135->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title14_2135 = new TLatex();
    title14_2135->SetNDC(); 
    title14_2135->SetTextSize(0.05);
    title14_2135->SetTextAlign(22);  
    //title14_2135->DrawLatex(0.5, 0.97, "|t| Distribution with Incoherent #phi RP Vetoes Only");  

    TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

    // Add labels
    TLatex* ep = new TLatex(0.13, 0.85, "ePIC");
    ep->SetNDC();
    ep->SetTextFont(62);
    ep->SetTextSize(.038);
    ep->SetTextColor(kBlack);
    ep->Draw("same");

    TLatex* r421 = new TLatex(0.18, 0.85, " Simulation 25.10.2, 10x100 GeV");
	r421->SetNDC();
	r421->SetTextSize(30);
	r421->SetTextFont(43);
	r421->SetTextColor(kBlack);
	r421->Draw("same");

	TLatex* label_5 = new TLatex(0.17, 0.75, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->SetTextSize(30);
	label_5->SetTextFont(43);
	label_5->SetTextColor(kBlack);
	label_5->Draw("same");

    TLatex* r4211 = new TLatex(0.15, 0.8, "eAu #rightarrow e'Au'#phi");
	r4211->SetNDC();
	r4211->SetTextSize(30);
	r4211->SetTextFont(43);
	r4211->SetTextColor(kBlack);
	r4211->Draw("same"); 

    TLatex* r421111 = new TLatex(0.18, 0.15, "L_{int} = 10 fb^{-1}/A"); // preTDR lumi
	r421111->SetNDC();
	r421111->SetTextSize(25);
	r421111->SetTextFont(43);
	r421111->SetTextColor(kBlack);
	r421111->Draw("same");

    c14_2135->Print("./figures/plot_t_incoherent_RPonly_preTDR_proj.pdf");

	//TFile* fout_preTDR = new TFile("plot_t_incoherent_RPonly_MC_preTDR.root","RECREATE");
    //c14_2135->Write();
    //h_t_MC->Write();
    //h_t_REC_EEMC_cut->Write();
    //h_t_REC_wRES_cut_pi12->Write();
    //h_incoherent_phi_RPonly_L->Write();
    //h_incoherent_phi_RPonly_MC->Write();
    //h_incoherent_phi_RPonly_proj->Write();
    //fout_preTDR->Close();
}

void plot_incoherent_phi_rpVetoes()
{
	TFile* phi_t_file = TFile::Open("phi_sartre_25_10_2_eAu_10x100v5.root","READ");
	TH1D* h_t_MC = (TH1D*) phi_t_file->Get("h_t_MC");
	TH1D* h_t_REC_EEMC_cut = (TH1D*) phi_t_file->Get("h_t_REC_EEMC_cut");
	TH1D* h_t_REC_wRES_cut_pi12 = (TH1D*) phi_t_file->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_phi_sartre_events = (TH1D*) phi_t_file->Get("h_Nevents");

	TFile* incoherent_phi_rpVetoes = TFile::Open("phi_beagle_25_10_3_eAu_rpVetoes_10x100v4.root");
    TH1D* h_incoherent_phi_RPonly_MC = (TH1D*) incoherent_phi_rpVetoes->Get("h_t_MC");
    TH1D* h_incoherent_phi_RPonly_L = (TH1D*) incoherent_phi_rpVetoes->Get("h_t_REC_EEMC_cut");
    TH1D* h_incoherent_phi_RPonly_proj = (TH1D*) incoherent_phi_rpVetoes->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_incoherent_events = (TH1D*) incoherent_phi_rpVetoes->Get("h_Nevents");

	double dN_dt = 0.002;
    double sigma_coherent = 459.05; // nb
    double sigma_incoherent = 261.95; // nb
    double nEvents_coherent = h_phi_sartre_events->GetEntries();
    double nEvents_incoherent = h_incoherent_events->GetEntries();
    double early_science_lumi = 8121.83; // nb^-1
    int rebin_width = 2;   

	TCanvas* c14_2135 = new TCanvas("c14_2135","c14_2135",1,1,1200,800);
    gPad->SetLogy(1);
    gStyle->SetOptStat(0);
    h_t_MC->GetYaxis()->SetTitle("dN/d|t| [GeV/c]^{-2}");
	h_t_MC->GetXaxis()->SetTitle("|t| [GeV/c]^{2}");
    h_t_MC->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    double simu_coherent_lumi = nEvents_coherent/sigma_coherent; // nb^-1 
    double ratio_coherent_phi_ES = early_science_lumi/simu_coherent_lumi;
    for (int i = 1; i <= h_t_MC->GetNbinsX(); ++i) 
    {
        double binError = h_t_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_MC->SetBinError(i, newError);
    }
    h_t_MC->SetLineColor(kBlack);
	h_t_MC->Draw("HISTsame");

    h_t_REC_EEMC_cut->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    for (int i = 1; i <= h_t_REC_EEMC_cut->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_EEMC_cut->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_REC_EEMC_cut->SetBinError(i, newError);
    }
	h_t_REC_EEMC_cut->SetMarkerStyle(20); // method L RECO
	h_t_REC_EEMC_cut->SetMarkerColor(kP8Blue);
	//h_t_REC_EEMC_cut->Draw("PEsame");

    h_t_REC_wRES_cut_pi12->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_wRES_cut_pi12->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_t_REC_wRES_cut_pi12->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_wRES_cut_pi12->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_REC_wRES_cut_pi12->SetBinError(i, newError);
    }
	h_t_REC_wRES_cut_pi12->SetMarkerStyle(30);
	h_t_REC_wRES_cut_pi12->SetMarkerColor(kP8Pink);
    h_t_REC_wRES_cut_pi12->SetLineColor(kP8Pink);
	h_t_REC_wRES_cut_pi12->Draw("PEsame");

    h_incoherent_phi_RPonly_MC->Rebin(rebin_width);
    h_incoherent_phi_RPonly_MC->Scale(sigma_incoherent*(1/nEvents_incoherent)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_incoherent_lumi = nEvents_incoherent/sigma_incoherent; // nb^-1 
    double ratio_incoherent_ES = early_science_lumi/simu_incoherent_lumi;
    for (int i = 1; i <= h_incoherent_phi_RPonly_MC->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_RPonly_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_ES);
        h_incoherent_phi_RPonly_MC->SetBinError(i, newError);
    }
    h_incoherent_phi_RPonly_MC->SetMarkerStyle(21);
    h_incoherent_phi_RPonly_MC->SetMarkerColor(kP8Orange);
    h_incoherent_phi_RPonly_MC->SetLineColor(kP8Orange);
    //h_incoherent_phi_RPonly_MC->Draw("PEsame");

    h_incoherent_phi_RPonly_L->Rebin(rebin_width);
    h_incoherent_phi_RPonly_L->Scale(sigma_incoherent*(1/nEvents_incoherent)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_incoherent_phi_RPonly_L->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_RPonly_L->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_ES);
        h_incoherent_phi_RPonly_L->SetBinError(i, newError);
    }
    h_incoherent_phi_RPonly_L->SetMarkerStyle(21);
    h_incoherent_phi_RPonly_L->SetMarkerColor(kP8Orange);
    h_incoherent_phi_RPonly_L->SetLineColor(kP8Orange);
    //h_incoherent_phi_RPonly_L->Draw("PEsame");

    h_incoherent_phi_RPonly_proj->Rebin(rebin_width);
    h_incoherent_phi_RPonly_proj->Scale(sigma_incoherent*(1/nEvents_incoherent)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_incoherent_phi_RPonly_proj->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_RPonly_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_ES);
        h_incoherent_phi_RPonly_proj->SetBinError(i, newError);
    }
    h_incoherent_phi_RPonly_proj->SetMarkerStyle(21);
    h_incoherent_phi_RPonly_proj->SetMarkerColor(kP8Orange);
    h_incoherent_phi_RPonly_proj->SetLineColor(kP8Orange);
    h_incoherent_phi_RPonly_proj->Draw("PEsame");
    
	TLegend *w14_2135 = new TLegend(0.6,0.7,0.72,0.9);
	w14_2135->AddEntry(h_t_MC,"Coherent #phi: MC", "L");
    w14_2135->AddEntry(h_t_REC_wRES_cut_pi12,"Coherent #phi: #theta_{max}= #pi/12", "P");    
    //w14_2135->AddEntry(h_t_REC_EEMC_cut,"Coherent #phi: L", "P");
    //w14_2135->AddEntry(h_incoherent_phi_RPonly_MC,"Incoherent #phi: MC", "P");
    //w14_2135->AddEntry(h_incoherent_phi_RPonly_L,"Incoherent #phi: L", "P");
    w14_2135->AddEntry(h_incoherent_phi_RPonly_proj,"Incoherent #phi: #theta_{max}= #pi/12", "P");
    w14_2135->SetBorderSize(0);   
    w14_2135->SetFillStyle(0);
    w14_2135->SetTextSize(30);
	w14_2135->SetTextFont(45);
    w14_2135->Draw("same");	

    c14_2135->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title14_2135 = new TLatex();
    title14_2135->SetNDC(); 
    title14_2135->SetTextSize(0.05);
    title14_2135->SetTextAlign(22);  
    //title14_2135->DrawLatex(0.5, 0.97, "|t| Distribution with Incoherent #phi RP Vetoes Only");  

    TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

    // Add labels
    TLatex* ep = new TLatex(0.13, 0.85, "ePIC");
    ep->SetNDC();
    ep->SetTextFont(62);
    ep->SetTextSize(.038);
    ep->SetTextColor(kBlack);
    ep->Draw("same");

    TLatex* r421 = new TLatex(0.18, 0.85, " Simulation 25.10.2, 10x100 GeV");
	r421->SetNDC();
	r421->SetTextSize(30);
	r421->SetTextFont(43);
	r421->SetTextColor(kBlack);
	r421->Draw("same");

	TLatex* label_5 = new TLatex(0.17, 0.75, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->SetTextSize(30);
	label_5->SetTextFont(43);
	label_5->SetTextColor(kBlack);
	label_5->Draw("same");

    TLatex* r4211 = new TLatex(0.15, 0.8, "eAu #rightarrow e'Au'#phi");
	r4211->SetNDC();
	r4211->SetTextSize(30);
	r4211->SetTextFont(43);
	r4211->SetTextColor(kBlack);
	r4211->Draw("same"); 

    TLatex* r42111 = new TLatex(0.18, 0.15, "L_{int} = 1.6 fb^{-1}/A"); // ES lumi
	r42111->SetNDC();
	r42111->SetTextSize(30);
	r42111->SetTextFont(43);
	r42111->SetTextColor(kBlack);
	r42111->Draw("same");

    c14_2135->Print("./figures/plot_t_incoherent_RPonly_ES_proj.pdf");

	//TFile* fout_ES = new TFile("plot_t_incoherent_RPonly_MC_ES.root","RECREATE");
    //c14_2135->Write();
    //h_t_MC->Write();
    //h_t_REC_EEMC_cut->Write();
    //h_t_REC_wRES_cut_pi12->Write();
    //h_incoherent_phi_RPonly_L->Write();
    //h_incoherent_phi_RPonly_MC->Write();
    //h_incoherent_phi_RPonly_proj->Write();
    //fout_ES->Close();
}

void plot_incoherent_phi_omdVetoes_preTDR()
{
	TFile* phi_t_file = TFile::Open("phi_sartre_25_10_2_eAu_10x100v5.root","READ");
	TH1D* h_t_MC = (TH1D*) phi_t_file->Get("h_t_MC");
	TH1D* h_t_REC_EEMC_cut = (TH1D*) phi_t_file->Get("h_t_REC_EEMC_cut");
	TH1D* h_t_REC_wRES_cut_pi12 = (TH1D*) phi_t_file->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_phi_sartre_events = (TH1D*) phi_t_file->Get("h_Nevents");

	TFile* incoherent_phi_omdVetoes = TFile::Open("phi_beagle_25_10_3_eAu_omdVetoes_10x100v4.root");
    TH1D* h_incoherent_phi_OMDonly_MC = (TH1D*) incoherent_phi_omdVetoes->Get("h_t_MC");
    TH1D* h_incoherent_phi_OMDonly_L = (TH1D*) incoherent_phi_omdVetoes->Get("h_t_REC_EEMC_cut");
    TH1D* h_incoherent_phi_OMDonly_proj = (TH1D*) incoherent_phi_omdVetoes->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_incoherent_events = (TH1D*) incoherent_phi_omdVetoes->Get("h_Nevents");

	double dN_dt = 0.002;
    double sigma_coherent = 459.05; // nb
    double sigma_incoherent = 261.95; // nb
    double nEvents_coherent = h_phi_sartre_events->GetEntries();
    double nEvents_incoherent = h_incoherent_events->GetEntries();
    double preTDR_lumi = 50761.4; // nb^-1
    int rebin_width = 2;   

	TCanvas* c14_2134 = new TCanvas("c14_2134","c14_2134",1,1,1200,800);
    gPad->SetLogy(1);
    gStyle->SetOptStat(0);
    h_t_MC->GetYaxis()->SetTitle("dN/d|t| [GeV/c]^{-2}");
	h_t_MC->GetXaxis()->SetTitle("|t| [GeV/c]^{2}");
    h_t_MC->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    double simu_coherent_lumi = nEvents_coherent/sigma_coherent; // nb^-1 
    double ratio_coherent_phi_preTDR = preTDR_lumi/simu_coherent_lumi;
    for (int i = 1; i <= h_t_MC->GetNbinsX(); ++i) 
    {
        double binError = h_t_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_MC->SetBinError(i, newError);
    }
    h_t_MC->SetLineColor(kBlack);
	h_t_MC->Draw("HISTsame");

    h_t_REC_EEMC_cut->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    for (int i = 1; i <= h_t_REC_EEMC_cut->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_EEMC_cut->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_EEMC_cut->SetBinError(i, newError);
    }
	h_t_REC_EEMC_cut->SetMarkerStyle(20); // method L RECO
	h_t_REC_EEMC_cut->SetMarkerColor(kP8Blue);
	//h_t_REC_EEMC_cut->Draw("PEsame");

    h_t_REC_wRES_cut_pi12->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_wRES_cut_pi12->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_t_REC_wRES_cut_pi12->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_wRES_cut_pi12->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_wRES_cut_pi12->SetBinError(i, newError);
    }
	h_t_REC_wRES_cut_pi12->SetMarkerStyle(30);
	h_t_REC_wRES_cut_pi12->SetMarkerColor(kP8Pink);
    h_t_REC_wRES_cut_pi12->SetLineColor(kP8Pink);
	h_t_REC_wRES_cut_pi12->Draw("PEsame");

    h_incoherent_phi_OMDonly_MC->Rebin(rebin_width);
    h_incoherent_phi_OMDonly_MC->Scale(sigma_incoherent*(1/nEvents_incoherent)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_incoherent_lumi = nEvents_incoherent/sigma_incoherent; // nb^-1 
    double ratio_incoherent_preTDR = preTDR_lumi/simu_incoherent_lumi;
    for (int i = 1; i <= h_incoherent_phi_OMDonly_MC->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_OMDonly_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_preTDR);
        h_incoherent_phi_OMDonly_MC->SetBinError(i, newError);
    }
    h_incoherent_phi_OMDonly_MC->SetMarkerStyle(21);
    h_incoherent_phi_OMDonly_MC->SetMarkerColor(kP8Orange);
    h_incoherent_phi_OMDonly_MC->SetLineColor(kP8Orange);
    //h_incoherent_phi_OMDonly_MC->Draw("PEsame");

    h_incoherent_phi_OMDonly_L->Rebin(rebin_width);
    h_incoherent_phi_OMDonly_L->Scale(sigma_incoherent*(1/nEvents_incoherent)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_incoherent_phi_OMDonly_L->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_OMDonly_L->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_preTDR);
        h_incoherent_phi_OMDonly_L->SetBinError(i, newError);
    }
    h_incoherent_phi_OMDonly_L->SetMarkerStyle(21);
    h_incoherent_phi_OMDonly_L->SetMarkerColor(kP8Orange);
    h_incoherent_phi_OMDonly_L->SetLineColor(kP8Orange);
    //h_incoherent_phi_OMDonly_L->Draw("PEsame");

    h_incoherent_phi_OMDonly_proj->Rebin(rebin_width);
    h_incoherent_phi_OMDonly_proj->Scale(sigma_incoherent*(1/nEvents_incoherent)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_incoherent_phi_OMDonly_proj->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_OMDonly_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_preTDR);
        h_incoherent_phi_OMDonly_proj->SetBinError(i, newError);
    }
    h_incoherent_phi_OMDonly_proj->SetMarkerStyle(21);
    h_incoherent_phi_OMDonly_proj->SetMarkerColor(kP8Orange);
    h_incoherent_phi_OMDonly_proj->SetLineColor(kP8Orange);
    h_incoherent_phi_OMDonly_proj->Draw("PEsame");

	TLegend *w14_2134 = new TLegend(0.6,0.75,0.72,0.9);
    w14_2134->AddEntry(h_t_MC,"Coherent #phi: MC", "L");
    w14_2134->AddEntry(h_t_REC_wRES_cut_pi12,"Coherent #phi: #theta_{max}= #pi/12", "P");    
    //w14_2134->AddEntry(h_t_REC_EEMC_cut,"Coherent #phi: L", "P");
    //w14_2134->AddEntry(h_incoherent_phi_OMDonly_MC,"Incoherent #phi: MC", "P");
    //w14_2134->AddEntry(h_incoherent_phi_OMDonly_L,"Incoherent #phi: L", "P");
    w14_2134->AddEntry(h_incoherent_phi_OMDonly_proj,"Incoherent #phi: #theta_{max}= #pi/12", "P");
    w14_2134->SetBorderSize(0);   
    w14_2134->SetFillStyle(0);
    w14_2134->SetTextSize(30);
	w14_2134->SetTextFont(45);
    w14_2134->Draw("same");	

    c14_2134->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title14_2134 = new TLatex();
    title14_2134->SetNDC(); 
    title14_2134->SetTextSize(0.05);
    title14_2134->SetTextAlign(22);  
    //title14_2134->DrawLatex(0.5, 0.97, "|t| Distribution with Incoherent #phi OMD Vetoes Only");  

    TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

    // Add labels
    TLatex* ep = new TLatex(0.13, 0.85, "ePIC");
    ep->SetNDC();
    ep->SetTextFont(62);
    ep->SetTextSize(.038);
    ep->SetTextColor(kBlack);
    ep->Draw("same");

    TLatex* r421 = new TLatex(0.18, 0.85, " Simulation 25.10.2, 10x100 GeV");
	r421->SetNDC();
	r421->SetTextSize(30);
	r421->SetTextFont(43);
	r421->SetTextColor(kBlack);
	r421->Draw("same");

	TLatex* label_5 = new TLatex(0.17, 0.75, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->SetTextSize(30);
	label_5->SetTextFont(43);
	label_5->SetTextColor(kBlack);
	label_5->Draw("same");

    TLatex* r4211 = new TLatex(0.15, 0.8, "eAu #rightarrow e'Au'#phi");
	r4211->SetNDC();
	r4211->SetTextSize(30);
	r4211->SetTextFont(43);
	r4211->SetTextColor(kBlack);
	r4211->Draw("same"); 

    TLatex* r421111 = new TLatex(0.18, 0.15, "L_{int} = 10 fb^{-1}/A"); // preTDR lumi
	r421111->SetNDC();
	r421111->SetTextSize(25);
	r421111->SetTextFont(43);
	r421111->SetTextColor(kBlack);
	r421111->Draw("same");

    c14_2134->Print("./figures/plot_t_incoherent_OMDonly_preTDR_proj.pdf");

	//TFile* fout_preTDR = new TFile("plot_t_incoherent_OMDonly_MC_preTDR.root","RECREATE");
    //c14_2134->Write();
    //h_t_MC->Write();
    //h_t_REC_EEMC_cut->Write();
    //h_t_REC_wRES_cut_pi12->Write();
    //h_incoherent_phi_OMDonly_MC->Write();
    //h_incoherent_phi_OMDonly_L->Write();
    //h_incoherent_phi_OMDonly_proj->Write();
    //fout_preTDR->Close();
}

void plot_incoherent_phi_omdVetoes_ES()
{
	TFile* phi_t_file = TFile::Open("phi_sartre_25_10_2_eAu_10x100v5.root","READ");
	TH1D* h_t_MC = (TH1D*) phi_t_file->Get("h_t_MC");
	TH1D* h_t_REC_EEMC_cut = (TH1D*) phi_t_file->Get("h_t_REC_EEMC_cut");
	TH1D* h_t_REC_wRES_cut_pi12 = (TH1D*) phi_t_file->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_phi_sartre_events = (TH1D*) phi_t_file->Get("h_Nevents");

	TFile* incoherent_phi_omdVetoes = TFile::Open("phi_beagle_25_10_3_eAu_omdVetoes_10x100v4.root");
    TH1D* h_incoherent_phi_OMDonly_MC = (TH1D*) incoherent_phi_omdVetoes->Get("h_t_MC");
    TH1D* h_incoherent_phi_OMDonly_L = (TH1D*) incoherent_phi_omdVetoes->Get("h_t_REC_EEMC_cut");
    TH1D* h_incoherent_phi_OMDonly_proj = (TH1D*) incoherent_phi_omdVetoes->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_incoherent_events = (TH1D*) incoherent_phi_omdVetoes->Get("h_Nevents");

	double dN_dt = 0.002;
    double sigma_coherent = 459.05; // nb
    double sigma_incoherent = 261.95; // nb
    double nEvents_coherent = h_phi_sartre_events->GetEntries();
    double nEvents_incoherent = h_incoherent_events->GetEntries();
    double early_science_lumi = 8121.83; // nb^-1
    int rebin_width = 2;   

	TCanvas* c14_2134 = new TCanvas("c14_2134","c14_2134",1,1,1200,800);
    gPad->SetLogy(1);
    gStyle->SetOptStat(0);
    h_t_MC->GetYaxis()->SetTitle("dN/d|t| [GeV/c]^{-2}");
	h_t_MC->GetXaxis()->SetTitle("|t| [GeV/c]^{2}");
    h_t_MC->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    double simu_coherent_lumi = nEvents_coherent/sigma_coherent; // nb^-1 
    double ratio_coherent_phi_ES = early_science_lumi/simu_coherent_lumi;
    for (int i = 1; i <= h_t_MC->GetNbinsX(); ++i) 
    {
        double binError = h_t_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_MC->SetBinError(i, newError);
    }
    h_t_MC->SetLineColor(kBlack);
	h_t_MC->Draw("HISTsame");

    h_t_REC_EEMC_cut->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    for (int i = 1; i <= h_t_REC_EEMC_cut->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_EEMC_cut->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_REC_EEMC_cut->SetBinError(i, newError);
    }
	h_t_REC_EEMC_cut->SetMarkerStyle(20); // method L RECO
	h_t_REC_EEMC_cut->SetMarkerColor(kP8Blue);
	//h_t_REC_EEMC_cut->Draw("PEsame");

    h_t_REC_wRES_cut_pi12->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_wRES_cut_pi12->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_t_REC_wRES_cut_pi12->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_wRES_cut_pi12->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_REC_wRES_cut_pi12->SetBinError(i, newError);
    }
	h_t_REC_wRES_cut_pi12->SetMarkerStyle(30);
	h_t_REC_wRES_cut_pi12->SetMarkerColor(kP8Pink);
    h_t_REC_wRES_cut_pi12->SetLineColor(kP8Pink);
	h_t_REC_wRES_cut_pi12->Draw("PEsame");

    h_incoherent_phi_OMDonly_MC->Rebin(rebin_width);
    h_incoherent_phi_OMDonly_MC->Scale(sigma_incoherent*(1/nEvents_incoherent)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_incoherent_lumi = nEvents_incoherent/sigma_incoherent; // nb^-1 
    double ratio_incoherent_ES = early_science_lumi/simu_incoherent_lumi;
    for (int i = 1; i <= h_incoherent_phi_OMDonly_MC->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_OMDonly_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_ES);
        h_incoherent_phi_OMDonly_MC->SetBinError(i, newError);
    }
    h_incoherent_phi_OMDonly_MC->SetMarkerStyle(21);
    h_incoherent_phi_OMDonly_MC->SetMarkerColor(kP8Orange);
    h_incoherent_phi_OMDonly_MC->SetLineColor(kP8Orange);
    //h_incoherent_phi_OMDonly_MC->Draw("PEsame");

    h_incoherent_phi_OMDonly_L->Rebin(rebin_width);
    h_incoherent_phi_OMDonly_L->Scale(sigma_incoherent*(1/nEvents_incoherent)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_incoherent_phi_OMDonly_L->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_OMDonly_L->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_ES);
        h_incoherent_phi_OMDonly_L->SetBinError(i, newError);
    }
    h_incoherent_phi_OMDonly_L->SetMarkerStyle(21);
    h_incoherent_phi_OMDonly_L->SetMarkerColor(kP8Orange);
    h_incoherent_phi_OMDonly_L->SetLineColor(kP8Orange);
    //h_incoherent_phi_OMDonly_L->Draw("PEsame");

    h_incoherent_phi_OMDonly_proj->Rebin(rebin_width);
    h_incoherent_phi_OMDonly_proj->Scale(sigma_incoherent*(1/nEvents_incoherent)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_incoherent_phi_OMDonly_proj->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_OMDonly_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_ES);
        h_incoherent_phi_OMDonly_proj->SetBinError(i, newError);
    }
    h_incoherent_phi_OMDonly_proj->SetMarkerStyle(21);
    h_incoherent_phi_OMDonly_proj->SetMarkerColor(kP8Orange);
    h_incoherent_phi_OMDonly_proj->SetLineColor(kP8Orange);
    h_incoherent_phi_OMDonly_proj->Draw("PEsame");

	TLegend *w14_2134 = new TLegend(0.6,0.75,0.72,0.9);
    w14_2134->AddEntry(h_t_MC,"Coherent #phi: MC", "L");
    w14_2134->AddEntry(h_t_REC_wRES_cut_pi12,"Coherent #phi: #theta_{max}= #pi/12", "P");    
    //w14_2134->AddEntry(h_t_REC_EEMC_cut,"Coherent #phi: L", "P");
    //w14_2134->AddEntry(h_incoherent_phi_OMDonly_MC,"Incoherent #phi: MC", "P");
    //w14_2134->AddEntry(h_incoherent_phi_OMDonly_L,"Incoherent #phi: L", "P");
    w14_2134->AddEntry(h_incoherent_phi_OMDonly_proj,"Incoherent #phi: #theta_{max}= #pi/12", "P");
    w14_2134->SetBorderSize(0);   
    w14_2134->SetFillStyle(0);
    w14_2134->SetTextSize(30);
	w14_2134->SetTextFont(45);
    w14_2134->Draw("same");	

    c14_2134->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title14_2134 = new TLatex();
    title14_2134->SetNDC(); 
    title14_2134->SetTextSize(0.05);
    title14_2134->SetTextAlign(22);  
    //title14_2134->DrawLatex(0.5, 0.97, "|t| Distribution with Incoherent #phi OMD Vetoes Only");  

    TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

    // Add labels
    TLatex* ep = new TLatex(0.13, 0.85, "ePIC");
    ep->SetNDC();
    ep->SetTextFont(62);
    ep->SetTextSize(.038);
    ep->SetTextColor(kBlack);
    ep->Draw("same");

    TLatex* r421 = new TLatex(0.18, 0.85, " Simulation 25.10.2, 10x100 GeV");
	r421->SetNDC();
	r421->SetTextSize(30);
	r421->SetTextFont(43);
	r421->SetTextColor(kBlack);
	r421->Draw("same");

	TLatex* label_5 = new TLatex(0.17, 0.75, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->SetTextSize(30);
	label_5->SetTextFont(43);
	label_5->SetTextColor(kBlack);
	label_5->Draw("same");

    TLatex* r4211 = new TLatex(0.15, 0.8, "eAu #rightarrow e'Au'#phi");
	r4211->SetNDC();
	r4211->SetTextSize(30);
	r4211->SetTextFont(43);
	r4211->SetTextColor(kBlack);
	r4211->Draw("same"); 

    TLatex* r42111 = new TLatex(0.18, 0.15, "L_{int} = 1.6 fb^{-1}/A"); // ES lumi
	r42111->SetNDC();
	r42111->SetTextSize(30);
	r42111->SetTextFont(43);
	r42111->SetTextColor(kBlack);
	r42111->Draw("same");

    c14_2134->Print("./figures/plot_t_incoherent_OMDonly_ES_proj.pdf");

	//TFile* fout_ES = new TFile("plot_t_incoherent_OMDonly_MC_ES.root","RECREATE");
    //c14_2134->Write();
    //h_t_MC->Write();
    //h_t_REC_EEMC_cut->Write();
    //h_t_REC_wRES_cut_pi12->Write();
    //h_incoherent_phi_OMDonly_MC->Write();
    //h_incoherent_phi_OMDonly_L->Write();
    //h_incoherent_phi_OMDonly_proj->Write();
    //fout_ES->Close();
}

void plot_incoherent_phi_zdcVetoes_preTDR()
{
	TFile* phi_t_file = TFile::Open("phi_sartre_25_10_2_eAu_10x100v5.root","READ");
	TH1D* h_t_MC = (TH1D*) phi_t_file->Get("h_t_MC");
	TH1D* h_t_REC_EEMC_cut = (TH1D*) phi_t_file->Get("h_t_REC_EEMC_cut");
	TH1D* h_t_REC_wRES_cut_pi12 = (TH1D*) phi_t_file->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_phi_sartre_events = (TH1D*) phi_t_file->Get("h_Nevents");

	TFile* incoherent_phi_zdcVetoes = TFile::Open("phi_beagle_25_10_3_eAu_zdcVetoes_10x100v4.root");
    TH1D* h_incoherent_phi_ZDConly_MC = (TH1D*) incoherent_phi_zdcVetoes->Get("h_t_MC");
    TH1D* h_incoherent_phi_ZDConly_L = (TH1D*) incoherent_phi_zdcVetoes->Get("h_t_REC_EEMC_cut");
    TH1D* h_incoherent_phi_ZDConly_proj = (TH1D*) incoherent_phi_zdcVetoes->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_incoherent_events = (TH1D*) incoherent_phi_zdcVetoes->Get("h_Nevents");

	double dN_dt = 0.002;
    double sigma_coherent = 459.05; // nb
    double sigma_incoherent = 261.95; // nb
    double nEvents_coherent = h_phi_sartre_events->GetEntries();
    double nEvents_incoherent = h_incoherent_events->GetEntries();
    double preTDR_lumi = 50761.4; // nb^-1
    int rebin_width = 2; 

	TCanvas* c14_2136 = new TCanvas("c14_2136","c14_2136",1,1,1200,800);
    gPad->SetLogy(1);
    gStyle->SetOptStat(0);
    h_t_MC->GetYaxis()->SetTitle("dN/d|t| [GeV/c]^{-2}");
	h_t_MC->GetXaxis()->SetTitle("|t| [GeV/c]^{2}");
    h_t_MC->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    double simu_coherent_lumi = nEvents_coherent/sigma_coherent; // nb^-1 
    double ratio_coherent_phi_preTDR = preTDR_lumi/simu_coherent_lumi;
    for (int i = 1; i <= h_t_MC->GetNbinsX(); ++i) 
    {
        double binError = h_t_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_MC->SetBinError(i, newError);
    }
    h_t_MC->SetLineColor(kBlack);
	h_t_MC->Draw("HISTsame");

    h_t_REC_EEMC_cut->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    for (int i = 1; i <= h_t_REC_EEMC_cut->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_EEMC_cut->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_EEMC_cut->SetBinError(i, newError);
    }
	h_t_REC_EEMC_cut->SetMarkerStyle(20); // method L RECO
	h_t_REC_EEMC_cut->SetMarkerColor(kP8Blue);
	//h_t_REC_EEMC_cut->Draw("Psame");

    h_t_REC_wRES_cut_pi12->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_wRES_cut_pi12->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_t_REC_wRES_cut_pi12->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_wRES_cut_pi12->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_wRES_cut_pi12->SetBinError(i, newError);
    }
	h_t_REC_wRES_cut_pi12->SetMarkerStyle(30);
	h_t_REC_wRES_cut_pi12->SetMarkerColor(kP8Pink);
    h_t_REC_wRES_cut_pi12->SetLineColor(kP8Pink);
	h_t_REC_wRES_cut_pi12->Draw("Psame");

    h_incoherent_phi_ZDConly_MC->Rebin(rebin_width);
    h_incoherent_phi_ZDConly_MC->Scale(sigma_incoherent*(1/nEvents_incoherent)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_incoherent_lumi = nEvents_incoherent/sigma_incoherent; // nb^-1 
    double ratio_incoherent_preTDR = preTDR_lumi/simu_incoherent_lumi;
    for (int i = 1; i <= h_incoherent_phi_ZDConly_MC->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_ZDConly_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_preTDR);
        h_incoherent_phi_ZDConly_MC->SetBinError(i, newError);
    }
    h_incoherent_phi_ZDConly_MC->SetMarkerStyle(21);
    h_incoherent_phi_ZDConly_MC->SetMarkerColor(kP8Orange);
    h_incoherent_phi_ZDConly_MC->SetLineColor(kP8Orange);
    //h_incoherent_phi_ZDConly_MC->Draw("PEsame");

    h_incoherent_phi_ZDConly_L->Rebin(rebin_width);
    h_incoherent_phi_ZDConly_L->Scale(sigma_incoherent*(1/nEvents_incoherent)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_incoherent_phi_ZDConly_L->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_ZDConly_L->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_preTDR);
        h_incoherent_phi_ZDConly_L->SetBinError(i, newError);
    }
    h_incoherent_phi_ZDConly_L->SetMarkerStyle(21);
    h_incoherent_phi_ZDConly_L->SetMarkerColor(kP8Orange);
    h_incoherent_phi_ZDConly_L->SetLineColor(kP8Orange);
    //h_incoherent_phi_ZDConly_L->Draw("PEsame");

    h_incoherent_phi_ZDConly_proj->Rebin(rebin_width);
    h_incoherent_phi_ZDConly_proj->Scale(sigma_incoherent*(1/nEvents_incoherent)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_incoherent_phi_ZDConly_proj->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_ZDConly_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_preTDR);
        h_incoherent_phi_ZDConly_proj->SetBinError(i, newError);
    }
    h_incoherent_phi_ZDConly_proj->SetMarkerStyle(21);
    h_incoherent_phi_ZDConly_proj->SetMarkerColor(kP8Orange);
    h_incoherent_phi_ZDConly_proj->SetLineColor(kP8Orange);
    h_incoherent_phi_ZDConly_proj->Draw("PEsame");
    
	TLegend *w14_2136 = new TLegend(0.45,0.63,0.72,0.8);
    w14_2136->AddEntry(h_t_MC,"Coherent #phi: MC", "L");
    w14_2136->AddEntry(h_t_REC_wRES_cut_pi12,"Coherent #phi: #theta_{max}= #pi/12", "P");    
    //w14_2136->AddEntry(h_t_REC_EEMC_cut,"Coherent #phi: L", "P");
    //w14_2136->AddEntry(h_incoherent_phi_ZDConly_MC,"Incoherent #phi: MC", "P");
    //w14_2136->AddEntry(h_incoherent_phi_ZDConly_L,"Incoherent #phi: L", "P");
    w14_2136->AddEntry(h_incoherent_phi_ZDConly_proj,"Incoh. #phi: #theta_{max}= #pi/12", "P");
    w14_2136->SetBorderSize(0);   
    w14_2136->SetFillStyle(0);
    w14_2136->SetTextSize(30);
	w14_2136->SetTextFont(45);
    w14_2136->Draw("same");	

    c14_2136->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title14_2136 = new TLatex();
    title14_2136->SetNDC(); 
    title14_2136->SetTextSize(0.05);
    title14_2136->SetTextAlign(22);  
    //title14_2136->DrawLatex(0.5, 0.97, "|t| Distribution with Incoherent #phi ZDC Vetoes Only");  

    TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

    // Add labels
    TLatex* ep = new TLatex(0.13, 0.85, "ePIC");
    ep->SetNDC();
    ep->SetTextFont(62);
    ep->SetTextSize(.038);
    ep->SetTextColor(kBlack);
    ep->Draw("same");

    TLatex* r421 = new TLatex(0.18, 0.85, " Simulation 25.10.2, 10x100 GeV");
	r421->SetNDC();
	r421->SetTextSize(30);
	r421->SetTextFont(43);
	r421->SetTextColor(kBlack);
	r421->Draw("same");

	TLatex* label_5 = new TLatex(0.17, 0.75, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->SetTextSize(30);
	label_5->SetTextFont(43);
	label_5->SetTextColor(kBlack);
	label_5->Draw("same");

    TLatex* r4211 = new TLatex(0.15, 0.8, "eAu #rightarrow e'Au'#phi");
	r4211->SetNDC();
	r4211->SetTextSize(30);
	r4211->SetTextFont(43);
	r4211->SetTextColor(kBlack);
	r4211->Draw("same"); 

    TLatex* r421111 = new TLatex(0.18, 0.15, "L_{int} = 10 fb^{-1}/A"); // preTDR lumi
	r421111->SetNDC();
	r421111->SetTextSize(25);
	r421111->SetTextFont(43);
	r421111->SetTextColor(kBlack);
	r421111->Draw("same");

    c14_2136->Print("./figures/plot_t_incoherent_ZDConly_preTDR_proj.pdf");

	//TFile* fout_preTDR = new TFile("plot_t_incoherent_ZDConly_L_preTDR.root","RECREATE");
    //c14_2136->Write();
    //h_t_MC->Write();
    //h_t_REC_EEMC_cut->Write();
    //h_t_REC_wRES_cut_pi12->Write();
    //h_incoherent_phi_ZDConly_MC->Write();
    //h_incoherent_phi_ZDConly_L->Write();
    //h_incoherent_phi_ZDConly_proj->Write();
    //fout_preTDR->Close();
}

void plot_incoherent_phi_zdcVetoes_ES()
{
	TFile* phi_t_file = TFile::Open("phi_sartre_25_10_2_eAu_10x100v5.root","READ");
	TH1D* h_t_MC = (TH1D*) phi_t_file->Get("h_t_MC");
	TH1D* h_t_REC_EEMC_cut = (TH1D*) phi_t_file->Get("h_t_REC_EEMC_cut");
	TH1D* h_t_REC_wRES_cut_pi12 = (TH1D*) phi_t_file->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_phi_sartre_events = (TH1D*) phi_t_file->Get("h_Nevents");

	TFile* incoherent_phi_zdcVetoes = TFile::Open("phi_beagle_25_10_3_eAu_zdcVetoes_10x100v4.root");
    TH1D* h_incoherent_phi_ZDConly_MC = (TH1D*) incoherent_phi_zdcVetoes->Get("h_t_MC");
    TH1D* h_incoherent_phi_ZDConly_L = (TH1D*) incoherent_phi_zdcVetoes->Get("h_t_REC_EEMC_cut");
    TH1D* h_incoherent_phi_ZDConly_proj = (TH1D*) incoherent_phi_zdcVetoes->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_incoherent_events = (TH1D*) incoherent_phi_zdcVetoes->Get("h_Nevents");

	double dN_dt = 0.002;
    double sigma_coherent = 459.05; // nb
    double sigma_incoherent = 261.95; // nb
    double nEvents_coherent = h_phi_sartre_events->GetEntries();
    double nEvents_incoherent = h_incoherent_events->GetEntries();
    double early_science_lumi = 8121.83; // nb^-1
    int rebin_width = 2; 

	TCanvas* c14_2136 = new TCanvas("c14_2136","c14_2136",1,1,1200,800);
    gPad->SetLogy(1);
    gStyle->SetOptStat(0);
    h_t_MC->GetYaxis()->SetTitle("dN/d|t| [GeV/c]^{-2}");
	h_t_MC->GetXaxis()->SetTitle("|t| [GeV/c]^{2}");
    h_t_MC->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    double simu_coherent_lumi = nEvents_coherent/sigma_coherent; // nb^-1 
    double ratio_coherent_phi_ES = early_science_lumi/simu_coherent_lumi;
    for (int i = 1; i <= h_t_MC->GetNbinsX(); ++i) 
    {
        double binError = h_t_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_MC->SetBinError(i, newError);
    }
    h_t_MC->SetLineColor(kBlack);
	h_t_MC->Draw("HISTsame");

    h_t_REC_EEMC_cut->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    for (int i = 1; i <= h_t_REC_EEMC_cut->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_EEMC_cut->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_REC_EEMC_cut->SetBinError(i, newError);
    }
	h_t_REC_EEMC_cut->SetMarkerStyle(20); // method L RECO
	h_t_REC_EEMC_cut->SetMarkerColor(kP8Blue);
	//h_t_REC_EEMC_cut->Draw("Psame");

    h_t_REC_wRES_cut_pi12->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_wRES_cut_pi12->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_t_REC_wRES_cut_pi12->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_wRES_cut_pi12->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_REC_wRES_cut_pi12->SetBinError(i, newError);
    }
	h_t_REC_wRES_cut_pi12->SetMarkerStyle(30);
	h_t_REC_wRES_cut_pi12->SetMarkerColor(kP8Pink);
    h_t_REC_wRES_cut_pi12->SetLineColor(kP8Pink);
	h_t_REC_wRES_cut_pi12->Draw("Psame");

    h_incoherent_phi_ZDConly_MC->Rebin(rebin_width);
    h_incoherent_phi_ZDConly_MC->Scale(sigma_incoherent*(1/nEvents_incoherent)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_incoherent_lumi = nEvents_incoherent/sigma_incoherent; // nb^-1 
    double ratio_incoherent_ES = early_science_lumi/simu_incoherent_lumi;
    for (int i = 1; i <= h_incoherent_phi_ZDConly_MC->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_ZDConly_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_ES);
        h_incoherent_phi_ZDConly_MC->SetBinError(i, newError);
    }
    h_incoherent_phi_ZDConly_MC->SetMarkerStyle(21);
    h_incoherent_phi_ZDConly_MC->SetMarkerColor(kP8Orange);
    h_incoherent_phi_ZDConly_MC->SetLineColor(kP8Orange);
    //h_incoherent_phi_ZDConly_MC->Draw("PEsame");

    h_incoherent_phi_ZDConly_L->Rebin(rebin_width);
    h_incoherent_phi_ZDConly_L->Scale(sigma_incoherent*(1/nEvents_incoherent)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_incoherent_phi_ZDConly_L->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_ZDConly_L->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_ES);
        h_incoherent_phi_ZDConly_L->SetBinError(i, newError);
    }
    h_incoherent_phi_ZDConly_L->SetMarkerStyle(21);
    h_incoherent_phi_ZDConly_L->SetMarkerColor(kP8Orange);
    h_incoherent_phi_ZDConly_L->SetLineColor(kP8Orange);
    //h_incoherent_phi_ZDConly_L->Draw("PEsame");

    h_incoherent_phi_ZDConly_proj->Rebin(rebin_width);
    h_incoherent_phi_ZDConly_proj->Scale(sigma_incoherent*(1/nEvents_incoherent)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_incoherent_phi_ZDConly_proj->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_ZDConly_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_ES);
        h_incoherent_phi_ZDConly_proj->SetBinError(i, newError);
    }
    h_incoherent_phi_ZDConly_proj->SetMarkerStyle(21);
    h_incoherent_phi_ZDConly_proj->SetMarkerColor(kP8Orange);
    h_incoherent_phi_ZDConly_proj->SetLineColor(kP8Orange);
    h_incoherent_phi_ZDConly_proj->Draw("PEsame");
    
	TLegend *w14_2136 = new TLegend(0.45,0.63,0.72,0.8);
    w14_2136->AddEntry(h_t_MC,"Coherent #phi: MC", "L");
    w14_2136->AddEntry(h_t_REC_wRES_cut_pi12,"Coherent #phi: #theta_{max}= #pi/12", "P");    
    //w14_2136->AddEntry(h_t_REC_EEMC_cut,"Coherent #phi: L", "P");
    //w14_2136->AddEntry(h_incoherent_phi_ZDConly_MC,"Incoherent #phi: MC", "P");
    //w14_2136->AddEntry(h_incoherent_phi_ZDConly_L,"Incoherent #phi: L", "P");
    w14_2136->AddEntry(h_incoherent_phi_ZDConly_proj,"Incoh. #phi: #theta_{max}= #pi/12", "P");
    w14_2136->SetBorderSize(0);   
    w14_2136->SetFillStyle(0);
    w14_2136->SetTextSize(30);
	w14_2136->SetTextFont(45);
    w14_2136->Draw("same");	

    c14_2136->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title14_2136 = new TLatex();
    title14_2136->SetNDC(); 
    title14_2136->SetTextSize(0.05);
    title14_2136->SetTextAlign(22);  
    //title14_2136->DrawLatex(0.5, 0.97, "|t| Distribution with Incoherent #phi ZDC Vetoes Only");  

    TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

    // Add labels
    TLatex* ep = new TLatex(0.13, 0.85, "ePIC");
    ep->SetNDC();
    ep->SetTextFont(62);
    ep->SetTextSize(.038);
    ep->SetTextColor(kBlack);
    ep->Draw("same");

    TLatex* r421 = new TLatex(0.18, 0.85, " Simulation 25.10.2, 10x100 GeV");
	r421->SetNDC();
	r421->SetTextSize(30);
	r421->SetTextFont(43);
	r421->SetTextColor(kBlack);
	r421->Draw("same");

	TLatex* label_5 = new TLatex(0.17, 0.75, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->SetTextSize(30);
	label_5->SetTextFont(43);
	label_5->SetTextColor(kBlack);
	label_5->Draw("same");

    TLatex* r4211 = new TLatex(0.15, 0.8, "eAu #rightarrow e'Au'#phi");
	r4211->SetNDC();
	r4211->SetTextSize(30);
	r4211->SetTextFont(43);
	r4211->SetTextColor(kBlack);
	r4211->Draw("same"); 

    TLatex* r42111 = new TLatex(0.18, 0.15, "L_{int} = 1.6 fb^{-1}/A"); // ES lumi
	r42111->SetNDC();
	r42111->SetTextSize(30);
	r42111->SetTextFont(43);
	r42111->SetTextColor(kBlack);
	r42111->Draw("same");

    c14_2136->Print("./figures/plot_t_incoherent_ZDConly_ES_proj.pdf");

	//TFile* fout_ES = new TFile("plot_t_incoherent_ZDConly_L_ES.root","RECREATE");
    //c14_2136->Write();
    //h_t_MC->Write();
    //h_t_REC_EEMC_cut->Write();
    //h_t_REC_wRES_cut_pi12->Write();
    //h_incoherent_phi_ZDConly_MC->Write();
    //h_incoherent_phi_ZDConly_L->Write();
    //h_incoherent_phi_ZDConly_proj->Write();
    //fout_ES->Close();
}

void plot_incoherent_together_preTDR()
{
    TFile* phi_t_file = TFile::Open("phi_sartre_25_10_2_eAu_10x100v5.root","READ");
	TH1D* h_t_MC = (TH1D*) phi_t_file->Get("h_t_MC");
	TH1D* h_t_REC_EEMC_cut = (TH1D*) phi_t_file->Get("h_t_REC_EEMC_cut");
	TH1D* h_t_REC_wRES_cut_pi12 = (TH1D*) phi_t_file->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_phi_sartre_events = (TH1D*) phi_t_file->Get("h_Nevents");

	TFile* incoherent_phi_noVetoes = TFile::Open("phi_beagle_25_10_3_eAu_noVetoes_10x100v4.root","READ");
    TH1D* h_incoherent_phi_noVetoes_MC = (TH1D*) incoherent_phi_noVetoes->Get("h_t_MC");
    TH1D* h_incoherent_phi_noVetoes_L = (TH1D*) incoherent_phi_noVetoes->Get("h_t_REC_EEMC_cut");
    TH1D* h_incoherent_phi_noVetoes_proj = (TH1D*) incoherent_phi_noVetoes->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_incoherent_events = (TH1D*) incoherent_phi_noVetoes->Get("h_Nevents");

    TFile* incoherent_phi_noVetoes_before = TFile::Open("phi_beagle_25_10_3_eAu_noVetoes_10x100v3.root","READ");
    TH1D* h_incoherent_phi_noVetoes_MC_before = (TH1D*) incoherent_phi_noVetoes_before->Get("h_t_MC");
    TH1D* h_incoherent_phi_noVetoes_L_before = (TH1D*) incoherent_phi_noVetoes_before->Get("h_t_REC_EEMC_cut");
    TH1D* h_incoherent_phi_noVetoes_proj_before = (TH1D*) incoherent_phi_noVetoes_before->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_incoherent_events_before = (TH1D*) incoherent_phi_noVetoes_before->Get("h_Nevents");

    TFile* incoherent_phi_allDetectorVetoes = TFile::Open("phi_beagle_25_10_3_eAu_allDetectorVetoes_10x100v4.root");
    TH1D* h_incoherent_phi_onlyDetectors_MC = (TH1D*) incoherent_phi_allDetectorVetoes->Get("h_t_MC");
    TH1D* h_incoherent_phi_onlyDetectors_L = (TH1D*) incoherent_phi_allDetectorVetoes->Get("h_t_REC_EEMC_cut");
    TH1D* h_incoherent_phi_onlyDetectors_proj = (TH1D*) incoherent_phi_allDetectorVetoes->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_incoherent_events2 = (TH1D*) incoherent_phi_allDetectorVetoes->Get("h_Nevents");

    TFile* incoherent_phi_rpVetoes = TFile::Open("phi_beagle_25_10_3_eAu_rpVetoes_10x100v4.root");
    TH1D* h_incoherent_phi_RPonly_MC = (TH1D*) incoherent_phi_rpVetoes->Get("h_t_MC");
    TH1D* h_incoherent_phi_RPonly_L = (TH1D*) incoherent_phi_rpVetoes->Get("h_t_REC_EEMC_cut");
    TH1D* h_incoherent_phi_RPonly_proj = (TH1D*) incoherent_phi_rpVetoes->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_incoherent_events3 = (TH1D*) incoherent_phi_rpVetoes->Get("h_Nevents");

    TFile* incoherent_phi_omdVetoes = TFile::Open("phi_beagle_25_10_3_eAu_omdVetoes_10x100v4.root");
    TH1D* h_incoherent_phi_OMDonly_MC = (TH1D*) incoherent_phi_omdVetoes->Get("h_t_MC");
    TH1D* h_incoherent_phi_OMDonly_L = (TH1D*) incoherent_phi_omdVetoes->Get("h_t_REC_EEMC_cut");
    TH1D* h_incoherent_phi_OMDonly_proj = (TH1D*) incoherent_phi_omdVetoes->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_incoherent_events4 = (TH1D*) incoherent_phi_omdVetoes->Get("h_Nevents");

    TFile* incoherent_phi_zdcVetoes = TFile::Open("phi_beagle_25_10_3_eAu_zdcVetoes_10x100v4.root");
    TH1D* h_incoherent_phi_ZDConly_MC = (TH1D*) incoherent_phi_zdcVetoes->Get("h_t_MC");
    TH1D* h_incoherent_phi_ZDConly_L = (TH1D*) incoherent_phi_zdcVetoes->Get("h_t_REC_EEMC_cut");
    TH1D* h_incoherent_phi_ZDConly_proj = (TH1D*) incoherent_phi_zdcVetoes->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_incoherent_events5 = (TH1D*) incoherent_phi_zdcVetoes->Get("h_Nevents");

    TFile* incoherent_phi_etaVetoes = TFile::Open("phi_beagle_25_10_3_eAu_etaVetoes_10x100v4.root");
    TH1D* h_incoherent_phi_etaOnly_MC = (TH1D*) incoherent_phi_etaVetoes->Get("h_t_MC");
    TH1D* h_incoherent_phi_etaOnly_L = (TH1D*) incoherent_phi_etaVetoes->Get("h_t_REC_EEMC_cut");
    TH1D* h_incoherent_phi_etaOnly_proj = (TH1D*) incoherent_phi_etaVetoes->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_incoherent_events6 = (TH1D*) incoherent_phi_etaVetoes->Get("h_Nevents");

    TFile* incoherent_phi_allVetoes = TFile::Open("phi_beagle_25_10_3_eAu_allVetoes_10x100v4.root","READ");
    TH1D* h_incoherent_phi_allVetoes_MC = (TH1D*) incoherent_phi_allVetoes->Get("h_t_MC");
    TH1D* h_incoherent_phi_allVetoes_L = (TH1D*) incoherent_phi_allVetoes->Get("h_t_REC_EEMC_cut");
    TH1D* h_incoherent_phi_allVetoes_proj = (TH1D*) incoherent_phi_allVetoes->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_incoherent_events7 = (TH1D*) incoherent_phi_allVetoes->Get("h_Nevents");

    TFile* incoherent_phi_allVetoes_before = TFile::Open("phi_beagle_25_10_3_eAu_allVetoes_10x100v3.root","READ");
    TH1D* h_incoherent_phi_allVetoes_MC_before = (TH1D*) incoherent_phi_allVetoes_before->Get("h_t_MC");
    TH1D* h_incoherent_phi_allVetoes_L_before = (TH1D*) incoherent_phi_allVetoes_before->Get("h_t_REC_EEMC_cut");
    TH1D* h_incoherent_phi_allVetoes_proj_before = (TH1D*) incoherent_phi_allVetoes_before->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_incoherent_events7_before = (TH1D*) incoherent_phi_allVetoes_before->Get("h_Nevents");

    double dN_dt = 0.002;
    double sigma_coherent = 459.05; // nb
    double sigma_incoherent = 261.95; // nb
    double sigma_incoherent_before = 274; // nb
    double nEvents_coherent = h_phi_sartre_events->GetEntries();
    double nEvents_incoherent = h_incoherent_events->GetEntries();
    double nEvents_incoherent_before = h_incoherent_events_before->GetEntries();
    double nEvents_incoherent2 = h_incoherent_events2->GetEntries();
    double nEvents_incoherent3 = h_incoherent_events3->GetEntries();
    double nEvents_incoherent4 = h_incoherent_events4->GetEntries();
    double nEvents_incoherent5 = h_incoherent_events5->GetEntries();
    double nEvents_incoherent6 = h_incoherent_events6->GetEntries();
    double nEvents_incoherent7 = h_incoherent_events7->GetEntries();
    double nEvents_incoherent7_before = h_incoherent_events7_before->GetEntries();
    double preTDR_lumi = 50761.4; // nb^-1
    int rebin_width = 2;   

	TCanvas* c14_213 = new TCanvas("c14_213","c14_213",1,1,1200,800);
    gPad->SetLogy(1);
	gStyle->SetOptStat(0);
    h_t_MC->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)]^{2}");
	h_t_MC->GetXaxis()->SetTitle("|t| [GeV/c]^{2}");
    h_t_MC->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    double simu_coherent_lumi = nEvents_coherent/sigma_coherent; // nb^-1 
    double ratio_coherent_phi_preTDR = preTDR_lumi/simu_coherent_lumi;
    for (int i = 1; i <= h_t_MC->GetNbinsX(); ++i) 
    {
        double binError = h_t_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_MC->SetBinError(i, newError);
    }
    h_t_MC->SetLineWidth(1);
    h_t_MC->SetLineColor(kBlack);
	h_t_MC->Draw("HISTsame");

    h_t_REC_EEMC_cut->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    for (int i = 1; i <= h_t_REC_EEMC_cut->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_EEMC_cut->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_EEMC_cut->SetBinError(i, newError);
    }
	h_t_REC_EEMC_cut->SetMarkerStyle(20); // method L RECO
    h_t_REC_EEMC_cut->SetMarkerSize(1.3);
	h_t_REC_EEMC_cut->SetMarkerColor(kP8Blue);
	//h_t_REC_EEMC_cut->Draw("PEsame");

    h_t_REC_wRES_cut_pi12->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_wRES_cut_pi12->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_t_REC_wRES_cut_pi12->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_wRES_cut_pi12->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_wRES_cut_pi12->SetBinError(i, newError);
    }
	h_t_REC_wRES_cut_pi12->SetMarkerStyle(30);
	h_t_REC_wRES_cut_pi12->SetMarkerColor(kP8Pink);
    h_t_REC_wRES_cut_pi12->SetLineColor(kP8Pink);
	h_t_REC_wRES_cut_pi12->Draw("PEsame");

    h_incoherent_phi_noVetoes_proj->Rebin(rebin_width);
    h_incoherent_phi_noVetoes_proj->Scale(sigma_incoherent*(1/nEvents_incoherent)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_incoherent_lumi = nEvents_incoherent/sigma_incoherent; // nb^-1 
    double ratio_incoherent_preTDR = preTDR_lumi/simu_incoherent_lumi;
    for (int i = 1; i <= h_incoherent_phi_noVetoes_proj->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_noVetoes_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_preTDR);
        h_incoherent_phi_noVetoes_proj->SetBinError(i, newError);
    }
    h_incoherent_phi_noVetoes_proj->SetMarkerStyle(40);
    h_incoherent_phi_noVetoes_proj->SetMarkerColor(kP10Violet);
    h_incoherent_phi_noVetoes_proj->SetLineColor(kP10Violet);
    h_incoherent_phi_noVetoes_proj->Draw("PEsame");

    h_incoherent_phi_noVetoes_proj_before->Rebin(rebin_width);
    h_incoherent_phi_noVetoes_proj_before->Scale(sigma_incoherent_before*(1/nEvents_incoherent_before)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_incoherent_lumi_before = nEvents_incoherent_before/sigma_incoherent_before; // nb^-1 
    double ratio_incoherent_preTDR_before = preTDR_lumi/simu_incoherent_lumi_before;
    for (int i = 1; i <= h_incoherent_phi_noVetoes_proj_before->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_noVetoes_proj_before->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_preTDR_before);
        h_incoherent_phi_noVetoes_proj_before->SetBinError(i, newError);
    }
    h_incoherent_phi_noVetoes_proj_before->SetMarkerStyle(24);
    h_incoherent_phi_noVetoes_proj_before->SetMarkerColor(kP8Gray);
    h_incoherent_phi_noVetoes_proj_before->SetLineColor(kP8Gray);
    //h_incoherent_phi_noVetoes_proj_before->Draw("PEsame");

    h_incoherent_phi_onlyDetectors_proj->Rebin(rebin_width);
    h_incoherent_phi_onlyDetectors_proj->Scale(sigma_incoherent*(1/nEvents_incoherent2)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_incoherent_lumi2 = nEvents_incoherent2/sigma_incoherent; // nb^-1 
    double ratio_incoherent_preTDR2 = preTDR_lumi/simu_incoherent_lumi2;
    for (int i = 1; i <= h_incoherent_phi_onlyDetectors_proj->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_onlyDetectors_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_preTDR2);
        h_incoherent_phi_onlyDetectors_proj->SetBinError(i, newError);
    }
    h_incoherent_phi_onlyDetectors_proj->SetMarkerStyle(24);
    h_incoherent_phi_onlyDetectors_proj->SetMarkerColor(kP8Gray);
    h_incoherent_phi_onlyDetectors_proj->SetLineColor(kP8Gray);
    h_incoherent_phi_onlyDetectors_proj->Draw("PEsame");

    h_incoherent_phi_RPonly_proj->Rebin(rebin_width);
    h_incoherent_phi_RPonly_proj->Scale(sigma_incoherent*(1/nEvents_incoherent3)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_incoherent_lumi3 = nEvents_incoherent3/sigma_incoherent; // nb^-1 
    double ratio_incoherent_preTDR3 = preTDR_lumi/simu_incoherent_lumi3;
    for (int i = 1; i <= h_incoherent_phi_RPonly_proj->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_RPonly_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_preTDR3);
        h_incoherent_phi_RPonly_proj->SetBinError(i, newError);
    }
    h_incoherent_phi_RPonly_proj->SetMarkerStyle(31);
    h_incoherent_phi_RPonly_proj->SetMarkerColor(kP8Red);
    h_incoherent_phi_RPonly_proj->SetLineColor(kP8Red);
    h_incoherent_phi_RPonly_proj->Draw("PEsame");

    h_incoherent_phi_OMDonly_proj->Rebin(rebin_width);
    h_incoherent_phi_OMDonly_proj->Scale(sigma_incoherent*(1/nEvents_incoherent4)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_incoherent_lumi4 = nEvents_incoherent4/sigma_incoherent; // nb^-1 
    double ratio_incoherent_preTDR4 = preTDR_lumi/simu_incoherent_lumi4;
    for (int i = 1; i <= h_incoherent_phi_OMDonly_proj->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_OMDonly_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_preTDR4);
        h_incoherent_phi_OMDonly_proj->SetBinError(i, newError);
    }
    h_incoherent_phi_OMDonly_proj->SetMarkerStyle(47);
    h_incoherent_phi_OMDonly_proj->SetMarkerColor(kP8Blue);
    h_incoherent_phi_OMDonly_proj->SetLineColor(kP8Blue);
    h_incoherent_phi_OMDonly_proj->Draw("PEsame");

    h_incoherent_phi_ZDConly_proj->Rebin(rebin_width);
    h_incoherent_phi_ZDConly_proj->Scale(sigma_incoherent*(1/nEvents_incoherent5)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_incoherent_lumi5 = nEvents_incoherent5/sigma_incoherent; // nb^-1 
    double ratio_incoherent_preTDR5 = preTDR_lumi/simu_incoherent_lumi5;
    for (int i = 1; i <= h_incoherent_phi_ZDConly_proj->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_ZDConly_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_preTDR5);
        h_incoherent_phi_ZDConly_proj->SetBinError(i, newError);
    }
    h_incoherent_phi_ZDConly_proj->SetMarkerStyle(33);
    h_incoherent_phi_ZDConly_proj->SetMarkerColor(kP8Green);
    h_incoherent_phi_ZDConly_proj->SetLineColor(kP8Green);
    h_incoherent_phi_ZDConly_proj->Draw("PEsame");

    h_incoherent_phi_etaOnly_proj->Rebin(rebin_width);
    h_incoherent_phi_etaOnly_proj->Scale(sigma_incoherent*(1/nEvents_incoherent6)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_incoherent_lumi6 = nEvents_incoherent6/sigma_incoherent; // nb^-1 
    double ratio_incoherent_preTDR6 = preTDR_lumi/simu_incoherent_lumi6;
    for (int i = 1; i <= h_incoherent_phi_etaOnly_proj->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_etaOnly_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_preTDR6);
        h_incoherent_phi_etaOnly_proj->SetBinError(i, newError);
    }
    h_incoherent_phi_etaOnly_proj->SetMarkerStyle(23);
    h_incoherent_phi_etaOnly_proj->SetMarkerColor(kBlack);
    h_incoherent_phi_etaOnly_proj->SetLineColor(kBlack);
    h_incoherent_phi_etaOnly_proj->Draw("PEsame");

    h_incoherent_phi_allVetoes_proj->Rebin(rebin_width);
    h_incoherent_phi_allVetoes_proj->Scale(sigma_incoherent*(1/nEvents_incoherent7)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_incoherent_lumi7 = nEvents_incoherent7/sigma_incoherent; // nb^-1 
    double ratio_incoherent_preTDR7 = preTDR_lumi/simu_incoherent_lumi7;
    for (int i = 1; i <= h_incoherent_phi_allVetoes_proj->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_allVetoes_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_preTDR7);
        h_incoherent_phi_allVetoes_proj->SetBinError(i, newError);
    }
    h_incoherent_phi_allVetoes_proj->SetMarkerStyle(21);
    h_incoherent_phi_allVetoes_proj->SetMarkerColor(kP8Orange);
    h_incoherent_phi_allVetoes_proj->SetLineColor(kP8Orange);
    h_incoherent_phi_allVetoes_proj->Draw("PEsame");

    h_incoherent_phi_allVetoes_proj_before->Rebin(rebin_width);
    h_incoherent_phi_allVetoes_proj_before->Scale(sigma_incoherent_before*(1/nEvents_incoherent7_before)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_incoherent_lumi7_before = nEvents_incoherent7_before/sigma_incoherent_before; // nb^-1 
    double ratio_incoherent_preTDR7_before = preTDR_lumi/simu_incoherent_lumi7_before;
    for (int i = 1; i <= h_incoherent_phi_allVetoes_proj_before->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_allVetoes_proj_before->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_preTDR7_before);
        h_incoherent_phi_allVetoes_proj_before->SetBinError(i, newError);
    }
    h_incoherent_phi_allVetoes_proj_before->SetMarkerStyle(23);
    h_incoherent_phi_allVetoes_proj_before->SetMarkerColor(kBlack);
    h_incoherent_phi_allVetoes_proj_before->SetLineColor(kBlack);
    //h_incoherent_phi_allVetoes_proj_before->Draw("PEsame");

    TLegend *w14_213 = new TLegend(0.6,0.55,0.72,0.9);
    w14_213->AddEntry(h_t_MC," Coherent #phi: MC", "L");
    w14_213->AddEntry(h_t_REC_wRES_cut_pi12,"Coherent #phi: #theta_{max}= #pi/12 ", "P");    
    //w14_213->AddEntry(h_t_REC_EEMC_cut," Coherent #phi: L", "P");
    w14_213->AddEntry(h_incoherent_phi_noVetoes_proj,"Incoh. #phi: no vetoes", "P");
    //w14_213->AddEntry(h_incoherent_phi_noVetoes_proj_before,"Incoh. #phi: no vetoes (before)", "P");
    w14_213->AddEntry(h_incoherent_phi_RPonly_proj,"Incoh. #phi: RP", "P");
    w14_213->AddEntry(h_incoherent_phi_OMDonly_proj,"Incoh. #phi: OMD", "P");
    w14_213->AddEntry(h_incoherent_phi_ZDConly_proj,"Incoh. #phi: ZDC", "P");
    w14_213->AddEntry(h_incoherent_phi_onlyDetectors_proj,"Incoh. #phi: all detectors", "P");
    w14_213->AddEntry(h_incoherent_phi_etaOnly_proj,"Incoh. #phi: #eta cuts", "P");
    w14_213->AddEntry(h_incoherent_phi_allVetoes_proj,"Incoh. #phi: all vetoes", "P");
    //w14_213->AddEntry(h_incoherent_phi_allVetoes_proj_before,"Incoh. #phi: all vetoes (before)", "P");
	w14_213->SetBorderSize(0);   
    w14_213->SetFillStyle(0);
	w14_213->SetTextSize(30);
	w14_213->SetTextFont(45);
    w14_213->Draw("same");	

    c14_213->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title14_213 = new TLatex();
    title14_213->SetNDC(); 
    title14_213->SetTextSize(0.05);
    title14_213->SetTextAlign(22);  
    //title14_213->DrawLatex(0.5, 0.97, "|t| Distribution with Incoherent #phi No Vetoes");  

    TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

    // Add labels
    TLatex* ep = new TLatex(0.13, 0.85, "ePIC");
    ep->SetNDC();
    ep->SetTextFont(62);
    ep->SetTextSize(.038);
    ep->SetTextColor(kBlack);
    ep->Draw("same");

    TLatex* r421 = new TLatex(0.18, 0.85, " Simulation 25.10.2, 10x100 GeV");
	r421->SetNDC();
	r421->SetTextSize(30);
	r421->SetTextFont(43);
	r421->SetTextColor(kBlack);
	r421->Draw("same");

	TLatex* label_5 = new TLatex(0.17, 0.75, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->SetTextSize(30);
	label_5->SetTextFont(43);
	label_5->SetTextColor(kBlack);
	label_5->Draw("same");

    TLatex* r4211 = new TLatex(0.15, 0.8, "eAu #rightarrow e'Au'#phi");
	r4211->SetNDC();
	r4211->SetTextSize(30);
	r4211->SetTextFont(43);
	r4211->SetTextColor(kBlack);
	r4211->Draw("same"); 

    TLatex* r421111 = new TLatex(0.18, 0.15, "L_{int} = 10 fb^{-1}/A"); // preTDR lumi
	r421111->SetNDC();
	r421111->SetTextSize(25);
	r421111->SetTextFont(43);
	r421111->SetTextColor(kBlack);
	r421111->Draw("same");

    c14_213->Print("./figures/plot_t_incoherent_compareAll_preTDR_proj.pdf");
}

void plot_incoherent_together_ES()
{
    TFile* phi_t_file = TFile::Open("phi_sartre_25_10_2_eAu_10x100v5.root","READ");
	TH1D* h_t_MC = (TH1D*) phi_t_file->Get("h_t_MC");
	TH1D* h_t_REC_EEMC_cut = (TH1D*) phi_t_file->Get("h_t_REC_EEMC_cut");
	TH1D* h_t_REC_wRES_cut_pi12 = (TH1D*) phi_t_file->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_phi_sartre_events = (TH1D*) phi_t_file->Get("h_Nevents");

	TFile* incoherent_phi_noVetoes = TFile::Open("phi_beagle_25_10_3_eAu_noVetoes_10x100v4.root","READ");
    TH1D* h_incoherent_phi_noVetoes_MC = (TH1D*) incoherent_phi_noVetoes->Get("h_t_MC");
    TH1D* h_incoherent_phi_noVetoes_L = (TH1D*) incoherent_phi_noVetoes->Get("h_t_REC_EEMC_cut");
    TH1D* h_incoherent_phi_noVetoes_proj = (TH1D*) incoherent_phi_noVetoes->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_incoherent_events = (TH1D*) incoherent_phi_noVetoes->Get("h_Nevents");

    TFile* incoherent_phi_noVetoes_before = TFile::Open("phi_beagle_25_10_3_eAu_noVetoes_10x100v3.root","READ");
    TH1D* h_incoherent_phi_noVetoes_MC_before = (TH1D*) incoherent_phi_noVetoes_before->Get("h_t_MC");
    TH1D* h_incoherent_phi_noVetoes_L_before = (TH1D*) incoherent_phi_noVetoes_before->Get("h_t_REC_EEMC_cut");
    TH1D* h_incoherent_phi_noVetoes_proj_before = (TH1D*) incoherent_phi_noVetoes_before->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_incoherent_events_before = (TH1D*) incoherent_phi_noVetoes_before->Get("h_Nevents");

    TFile* incoherent_phi_allDetectorVetoes = TFile::Open("phi_beagle_25_10_3_eAu_allDetectorVetoes_10x100v4.root");
    TH1D* h_incoherent_phi_onlyDetectors_MC = (TH1D*) incoherent_phi_allDetectorVetoes->Get("h_t_MC");
    TH1D* h_incoherent_phi_onlyDetectors_L = (TH1D*) incoherent_phi_allDetectorVetoes->Get("h_t_REC_EEMC_cut");
    TH1D* h_incoherent_phi_onlyDetectors_proj = (TH1D*) incoherent_phi_allDetectorVetoes->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_incoherent_events2 = (TH1D*) incoherent_phi_allDetectorVetoes->Get("h_Nevents");

    TFile* incoherent_phi_rpVetoes = TFile::Open("phi_beagle_25_10_3_eAu_rpVetoes_10x100v4.root");
    TH1D* h_incoherent_phi_RPonly_MC = (TH1D*) incoherent_phi_rpVetoes->Get("h_t_MC");
    TH1D* h_incoherent_phi_RPonly_L = (TH1D*) incoherent_phi_rpVetoes->Get("h_t_REC_EEMC_cut");
    TH1D* h_incoherent_phi_RPonly_proj = (TH1D*) incoherent_phi_rpVetoes->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_incoherent_events3 = (TH1D*) incoherent_phi_rpVetoes->Get("h_Nevents");

    TFile* incoherent_phi_omdVetoes = TFile::Open("phi_beagle_25_10_3_eAu_omdVetoes_10x100v4.root");
    TH1D* h_incoherent_phi_OMDonly_MC = (TH1D*) incoherent_phi_omdVetoes->Get("h_t_MC");
    TH1D* h_incoherent_phi_OMDonly_L = (TH1D*) incoherent_phi_omdVetoes->Get("h_t_REC_EEMC_cut");
    TH1D* h_incoherent_phi_OMDonly_proj = (TH1D*) incoherent_phi_omdVetoes->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_incoherent_events4 = (TH1D*) incoherent_phi_omdVetoes->Get("h_Nevents");

    TFile* incoherent_phi_zdcVetoes = TFile::Open("phi_beagle_25_10_3_eAu_zdcVetoes_10x100v4.root");
    TH1D* h_incoherent_phi_ZDConly_MC = (TH1D*) incoherent_phi_zdcVetoes->Get("h_t_MC");
    TH1D* h_incoherent_phi_ZDConly_L = (TH1D*) incoherent_phi_zdcVetoes->Get("h_t_REC_EEMC_cut");
    TH1D* h_incoherent_phi_ZDConly_proj = (TH1D*) incoherent_phi_zdcVetoes->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_incoherent_events5 = (TH1D*) incoherent_phi_zdcVetoes->Get("h_Nevents");

    TFile* incoherent_phi_etaVetoes = TFile::Open("phi_beagle_25_10_3_eAu_etaVetoes_10x100v4.root");
    TH1D* h_incoherent_phi_etaOnly_MC = (TH1D*) incoherent_phi_etaVetoes->Get("h_t_MC");
    TH1D* h_incoherent_phi_etaOnly_L = (TH1D*) incoherent_phi_etaVetoes->Get("h_t_REC_EEMC_cut");
    TH1D* h_incoherent_phi_etaOnly_proj = (TH1D*) incoherent_phi_etaVetoes->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_incoherent_events6 = (TH1D*) incoherent_phi_etaVetoes->Get("h_Nevents");

    TFile* incoherent_phi_allVetoes = TFile::Open("phi_beagle_25_10_3_eAu_allVetoes_10x100v4.root","READ");
    TH1D* h_incoherent_phi_allVetoes_MC = (TH1D*) incoherent_phi_allVetoes->Get("h_t_MC");
    TH1D* h_incoherent_phi_allVetoes_L = (TH1D*) incoherent_phi_allVetoes->Get("h_t_REC_EEMC_cut");
    TH1D* h_incoherent_phi_allVetoes_proj = (TH1D*) incoherent_phi_allVetoes->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_incoherent_events7 = (TH1D*) incoherent_phi_allVetoes->Get("h_Nevents");

    TFile* incoherent_phi_allVetoes_before = TFile::Open("phi_beagle_25_10_3_eAu_allVetoes_10x100v3.root","READ");
    TH1D* h_incoherent_phi_allVetoes_MC_before = (TH1D*) incoherent_phi_allVetoes_before->Get("h_t_MC");
    TH1D* h_incoherent_phi_allVetoes_L_before = (TH1D*) incoherent_phi_allVetoes_before->Get("h_t_REC_EEMC_cut");
    TH1D* h_incoherent_phi_allVetoes_proj_before = (TH1D*) incoherent_phi_allVetoes_before->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_incoherent_events7_before = (TH1D*) incoherent_phi_allVetoes_before->Get("h_Nevents");

    double dN_dt = 0.002;
    double sigma_coherent = 459.05; // nb
    double sigma_incoherent = 261.95; // nb
    double sigma_incoherent_before = 274; // nb
    double nEvents_coherent = h_phi_sartre_events->GetEntries();
    double nEvents_incoherent = h_incoherent_events->GetEntries();
    double nEvents_incoherent_before = h_incoherent_events_before->GetEntries();
    double nEvents_incoherent2 = h_incoherent_events2->GetEntries();
    double nEvents_incoherent3 = h_incoherent_events3->GetEntries();
    double nEvents_incoherent4 = h_incoherent_events4->GetEntries();
    double nEvents_incoherent5 = h_incoherent_events5->GetEntries();
    double nEvents_incoherent6 = h_incoherent_events6->GetEntries();
    double nEvents_incoherent7 = h_incoherent_events7->GetEntries();
    double nEvents_incoherent7_before = h_incoherent_events7_before->GetEntries();
    double early_science_lumi = 8121.83; // nb^-1
    int rebin_width = 2;   

	TCanvas* c14_213 = new TCanvas("c14_213","c14_213",1,1,1200,800);
    c14_213->Divide(1,1,0.01,0.01);
    c14_213->cd(1);
    gPad->SetLogy(1);
	gStyle->SetOptStat(0);
    h_t_MC->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)]^{2}");
	h_t_MC->GetXaxis()->SetTitle("|t| [GeV/c]^{2}");
    h_t_MC->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    double simu_coherent_lumi = nEvents_coherent/sigma_coherent; // nb^-1 
    double ratio_coherent_phi_ES = early_science_lumi/simu_coherent_lumi;
    for (int i = 1; i <= h_t_MC->GetNbinsX(); ++i) 
    {
        double binError = h_t_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_MC->SetBinError(i, newError);
    }
    h_t_MC->SetLineWidth(2);
    h_t_MC->SetLineColor(kBlack);
	h_t_MC->Draw("HISTsame");

    h_t_REC_EEMC_cut->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    for (int i = 1; i <= h_t_REC_EEMC_cut->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_EEMC_cut->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_REC_EEMC_cut->SetBinError(i, newError);
    }
	h_t_REC_EEMC_cut->SetMarkerStyle(20); // method L RECO
    h_t_REC_EEMC_cut->SetMarkerSize(1.3);
	h_t_REC_EEMC_cut->SetMarkerColor(kP8Blue);
	//h_t_REC_EEMC_cut->Draw("PEsame");

    h_t_REC_wRES_cut_pi12->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_wRES_cut_pi12->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_t_REC_wRES_cut_pi12->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_wRES_cut_pi12->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_REC_wRES_cut_pi12->SetBinError(i, newError);
    }
	h_t_REC_wRES_cut_pi12->SetMarkerStyle(30);
	h_t_REC_wRES_cut_pi12->SetMarkerColor(kP8Pink);
    h_t_REC_wRES_cut_pi12->SetLineColor(kP8Pink);
	h_t_REC_wRES_cut_pi12->Draw("PEsame");

    h_incoherent_phi_noVetoes_proj->Rebin(rebin_width);
    h_incoherent_phi_noVetoes_proj->Scale(sigma_incoherent*(1/nEvents_incoherent)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_incoherent_lumi = nEvents_incoherent/sigma_incoherent; // nb^-1 
    double ratio_incoherent_ES = early_science_lumi/simu_incoherent_lumi;
    for (int i = 1; i <= h_incoherent_phi_noVetoes_proj->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_noVetoes_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_ES);
        h_incoherent_phi_noVetoes_proj->SetBinError(i, newError);
    }
    h_incoherent_phi_noVetoes_proj->SetMarkerStyle(40);
    h_incoherent_phi_noVetoes_proj->SetMarkerColor(kP10Violet);
    h_incoherent_phi_noVetoes_proj->SetLineColor(kP10Violet);
    h_incoherent_phi_noVetoes_proj->Draw("PEsame");

    h_incoherent_phi_noVetoes_proj_before->Rebin(rebin_width);
    h_incoherent_phi_noVetoes_proj_before->Scale(sigma_incoherent_before*(1/nEvents_incoherent_before)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_incoherent_lumi_before = nEvents_incoherent_before/sigma_incoherent_before; // nb^-1 
    double ratio_incoherent_ES_before = early_science_lumi/simu_incoherent_lumi_before;
    for (int i = 1; i <= h_incoherent_phi_noVetoes_proj_before->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_noVetoes_proj_before->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_ES_before);
        h_incoherent_phi_noVetoes_proj_before->SetBinError(i, newError);
    }
    h_incoherent_phi_noVetoes_proj_before->SetMarkerStyle(24);
    h_incoherent_phi_noVetoes_proj_before->SetMarkerColor(kP8Gray);
    h_incoherent_phi_noVetoes_proj_before->SetLineColor(kP8Gray);
    h_incoherent_phi_noVetoes_proj_before->Draw("PEsame");

    h_incoherent_phi_onlyDetectors_proj->Rebin(rebin_width);
    h_incoherent_phi_onlyDetectors_proj->Scale(sigma_incoherent*(1/nEvents_incoherent2)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_incoherent_lumi2 = nEvents_incoherent2/sigma_incoherent; // nb^-1 
    double ratio_incoherent_ES2 = early_science_lumi/simu_incoherent_lumi2;
    for (int i = 1; i <= h_incoherent_phi_onlyDetectors_proj->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_onlyDetectors_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_ES2);
        h_incoherent_phi_onlyDetectors_proj->SetBinError(i, newError);
    }
    h_incoherent_phi_onlyDetectors_proj->SetMarkerStyle(24);
    h_incoherent_phi_onlyDetectors_proj->SetMarkerColor(kP8Gray);
    h_incoherent_phi_onlyDetectors_proj->SetLineColor(kP8Gray);
    //h_incoherent_phi_onlyDetectors_proj->Draw("PEsame");

    h_incoherent_phi_RPonly_proj->Rebin(rebin_width);
    h_incoherent_phi_RPonly_proj->Scale(sigma_incoherent*(1/nEvents_incoherent3)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_incoherent_lumi3 = nEvents_incoherent3/sigma_incoherent; // nb^-1 
    double ratio_incoherent_ES3 = early_science_lumi/simu_incoherent_lumi3;
    for (int i = 1; i <= h_incoherent_phi_RPonly_proj->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_RPonly_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_ES3);
        h_incoherent_phi_RPonly_proj->SetBinError(i, newError);
    }
    h_incoherent_phi_RPonly_proj->SetMarkerStyle(31);
    h_incoherent_phi_RPonly_proj->SetMarkerColor(kP8Red);
    h_incoherent_phi_RPonly_proj->SetLineColor(kP8Red);
    //h_incoherent_phi_RPonly_proj->Draw("PEsame");

    h_incoherent_phi_OMDonly_proj->Rebin(rebin_width);
    h_incoherent_phi_OMDonly_proj->Scale(sigma_incoherent*(1/nEvents_incoherent4)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_incoherent_lumi4 = nEvents_incoherent4/sigma_incoherent; // nb^-1 
    double ratio_incoherent_ES4 = early_science_lumi/simu_incoherent_lumi4;
    for (int i = 1; i <= h_incoherent_phi_OMDonly_proj->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_OMDonly_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_ES4);
        h_incoherent_phi_OMDonly_proj->SetBinError(i, newError);
    }
    h_incoherent_phi_OMDonly_proj->SetMarkerStyle(47);
    h_incoherent_phi_OMDonly_proj->SetMarkerColor(kP8Blue);
    h_incoherent_phi_OMDonly_proj->SetLineColor(kP8Blue);
    //h_incoherent_phi_OMDonly_proj->Draw("PEsame");

    h_incoherent_phi_ZDConly_proj->Rebin(rebin_width);
    h_incoherent_phi_ZDConly_proj->Scale(sigma_incoherent*(1/nEvents_incoherent5)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_incoherent_lumi5 = nEvents_incoherent5/sigma_incoherent; // nb^-1 
    double ratio_incoherent_ES5 = early_science_lumi/simu_incoherent_lumi5;
    for (int i = 1; i <= h_incoherent_phi_ZDConly_proj->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_ZDConly_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_ES5);
        h_incoherent_phi_ZDConly_proj->SetBinError(i, newError);
    }
    h_incoherent_phi_ZDConly_proj->SetMarkerStyle(33);
    h_incoherent_phi_ZDConly_proj->SetMarkerColor(kP8Green);
    h_incoherent_phi_ZDConly_proj->SetLineColor(kP8Green);
    //h_incoherent_phi_ZDConly_proj->Draw("PEsame");

    h_incoherent_phi_etaOnly_proj->Rebin(rebin_width);
    h_incoherent_phi_etaOnly_proj->Scale(sigma_incoherent*(1/nEvents_incoherent6)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_incoherent_lumi6 = nEvents_incoherent6/sigma_incoherent; // nb^-1 
    double ratio_incoherent_ES6 = early_science_lumi/simu_incoherent_lumi6;
    for (int i = 1; i <= h_incoherent_phi_etaOnly_proj->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_etaOnly_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_ES6);
        h_incoherent_phi_etaOnly_proj->SetBinError(i, newError);
    }
    h_incoherent_phi_etaOnly_proj->SetMarkerStyle(23);
    h_incoherent_phi_etaOnly_proj->SetMarkerColor(kBlack);
    h_incoherent_phi_etaOnly_proj->SetLineColor(kBlack);
    //h_incoherent_phi_etaOnly_proj->Draw("PEsame");

    h_incoherent_phi_allVetoes_proj->Rebin(rebin_width);
    h_incoherent_phi_allVetoes_proj->Scale(sigma_incoherent*(1/nEvents_incoherent7)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_incoherent_lumi7 = nEvents_incoherent7/sigma_incoherent; // nb^-1 
    double ratio_incoherent_ES7 = early_science_lumi/simu_incoherent_lumi7;
    for (int i = 1; i <= h_incoherent_phi_allVetoes_proj->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_allVetoes_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_ES7);
        h_incoherent_phi_allVetoes_proj->SetBinError(i, newError);
    }
    h_incoherent_phi_allVetoes_proj->SetMarkerStyle(21);
    h_incoherent_phi_allVetoes_proj->SetMarkerColor(kP8Orange);
    h_incoherent_phi_allVetoes_proj->SetLineColor(kP8Orange);
    h_incoherent_phi_allVetoes_proj->Draw("PEsame");

    h_incoherent_phi_allVetoes_proj_before->Rebin(rebin_width);
    h_incoherent_phi_allVetoes_proj_before->Scale(sigma_incoherent_before*(1/nEvents_incoherent7_before)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_incoherent_lumi7_before = nEvents_incoherent7_before/sigma_incoherent_before; // nb^-1 
    double ratio_incoherent_ES7_before = early_science_lumi/simu_incoherent_lumi7_before;
    for (int i = 1; i <= h_incoherent_phi_allVetoes_proj_before->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_allVetoes_proj_before->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_ES7_before);
        h_incoherent_phi_allVetoes_proj_before->SetBinError(i, newError);
    }
    h_incoherent_phi_allVetoes_proj_before->SetMarkerStyle(23);
    h_incoherent_phi_allVetoes_proj_before->SetMarkerColor(kBlack);
    h_incoherent_phi_allVetoes_proj_before->SetLineColor(kBlack);
    h_incoherent_phi_allVetoes_proj_before->Draw("PEsame");

    TLegend *w14_213 = new TLegend(0.6,0.7,0.72,0.9);
    w14_213->AddEntry(h_t_MC," Coherent #phi: MC", "L");
    w14_213->AddEntry(h_t_REC_wRES_cut_pi12,"Coherent #phi: #theta_{max}= #pi/12 ", "P");    
    //w14_213->AddEntry(h_t_REC_EEMC_cut," Coherent #phi: L", "P");
    w14_213->AddEntry(h_incoherent_phi_noVetoes_proj,"Incoh. #phi: no vetoes", "P");
    w14_213->AddEntry(h_incoherent_phi_noVetoes_proj_before,"Incoh. #phi: no vetoes (before)", "P");
    //w14_213->AddEntry(h_incoherent_phi_RPonly_proj,"Incoh. #phi: RP", "P");
    //w14_213->AddEntry(h_incoherent_phi_OMDonly_proj,"Incoh. #phi: OMD", "P");
    //w14_213->AddEntry(h_incoherent_phi_ZDConly_proj,"Incoh. #phi: ZDC", "P");
    //w14_213->AddEntry(h_incoherent_phi_onlyDetectors_proj,"Incoh. #phi: all detectors", "P");
    //w14_213->AddEntry(h_incoherent_phi_etaOnly_proj,"Incoh. #phi: #eta cuts", "P");
    w14_213->AddEntry(h_incoherent_phi_allVetoes_proj,"Incoh. #phi: all vetoes", "P");
    w14_213->AddEntry(h_incoherent_phi_allVetoes_proj_before,"Incoh. #phi: all vetoes (before)", "P");
	w14_213->SetBorderSize(0);   
    w14_213->SetFillStyle(0);
	w14_213->SetTextSize(30);
	w14_213->SetTextFont(45);
    w14_213->Draw("same");	

    c14_213->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title14_213 = new TLatex();
    title14_213->SetNDC(); 
    title14_213->SetTextSize(0.05);
    title14_213->SetTextAlign(22);  
    //title14_213->DrawLatex(0.5, 0.97, "|t| Distribution with Incoherent #phi No Vetoes");  

    TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

    // Add labels
    TLatex* ep = new TLatex(0.13, 0.85, "ePIC");
    ep->SetNDC();
    ep->SetTextFont(62);
    ep->SetTextSize(.038);
    ep->SetTextColor(kBlack);
    ep->Draw("same");

    TLatex* r421 = new TLatex(0.18, 0.85, " Simulation 25.10.2, 10x100 GeV");
	r421->SetNDC();
	r421->SetTextSize(30);
	r421->SetTextFont(43);
	r421->SetTextColor(kBlack);
	r421->Draw("same");

	TLatex* label_5 = new TLatex(0.17, 0.75, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->SetTextSize(30);
	label_5->SetTextFont(43);
	label_5->SetTextColor(kBlack);
	label_5->Draw("same");

    TLatex* r4211 = new TLatex(0.15, 0.8, "eAu #rightarrow e'Au'#phi");
	r4211->SetNDC();
	r4211->SetTextSize(30);
	r4211->SetTextFont(43);
	r4211->SetTextColor(kBlack);
	r4211->Draw("same"); 

    TLatex* r42111 = new TLatex(0.18, 0.15, "L_{int} = 1.6 fb^{-1}/A"); // ES lumi
	r42111->SetNDC();
	r42111->SetTextSize(30);
	r42111->SetTextFont(43);
	r42111->SetTextColor(kBlack);
	r42111->Draw("same");

    c14_213->Print("./figures/plot_t_incoherent_compare_ES_incoherentBEFOREandNow.pdf");
}

void plot_t_rho_preTDR()
{
    TFile* phi_t_file = TFile::Open("phi_sartre_25_10_2_eAu_10x100v5.root","READ");
	TH1D* h_t_MC = (TH1D*) phi_t_file->Get("h_t_MC");
	TH1D* h_t_REC_EEMC_cut = (TH1D*) phi_t_file->Get("h_t_REC_EEMC_cut");
	TH1D* h_t_REC_wRES_cut_pi12 = (TH1D*) phi_t_file->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_t_REC_wRES_cut_pi6 = (TH1D*) phi_t_file->Get("h_t_REC_wRES_cut_pi6");
    TH1D* h_phi_sartre_events = (TH1D*) phi_t_file->Get("h_Nevents");

	TFile* sartre_rho_PID = TFile::Open("rhoPID_sartre_25_10_3_eAu_10x100v4.root","READ");
    TH1D* h_t_rho_MC_PID = (TH1D*) sartre_rho_PID->Get("h_t_MC");
    TH1D* h_t_rho_L_PID = (TH1D*) sartre_rho_PID->Get("h_t_REC_EEMC_cut");
    TH1D* h_t_rho_proj_PID = (TH1D*) sartre_rho_PID->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_rho_events_PID = (TH1D*) sartre_rho_PID->Get("h_Nevents");

    TFile* sartre_rho = TFile::Open("rho_sartre_25_10_3_eAu_10x100v4.root","READ");
    TH1D* h_t_rho_MC = (TH1D*) sartre_rho->Get("h_t_MC");
    TH1D* h_t_rho_L = (TH1D*) sartre_rho->Get("h_t_REC_EEMC_cut");
    TH1D* h_t_rho_proj = (TH1D*) sartre_rho->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_rho_events = (TH1D*) sartre_rho->Get("h_Nevents");

    double dN_dt = 0.002;
    double sigma_coherent = 459.05; // nb
    double sigma_rho = 5418.68; // nb
    double sigma_DIS = 50184.8; // nb
    double nEvents_coherent = h_phi_sartre_events->GetEntries();//6.2e6;
    double nEvents_rho_PID = h_rho_events_PID->GetEntries();//1.1e6;
    double nEvents_rho = h_rho_events->GetEntries();//1.1e6;
    double preTDR_lumi = 50761.4; // nb^-1
    int rebin_width = 2;   

    TCanvas* c14_213 = new TCanvas("c14_213","c14_213",1,1,1200,800);
    gPad->SetLogy(1);
	gStyle->SetOptStat(0);
    h_t_MC->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
	h_t_MC->GetXaxis()->SetTitle("|t| [GeV/c]^{2}");    
    h_t_MC->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    double simu_coherent_lumi = nEvents_coherent/sigma_coherent; // nb^-1 
    double ratio_coherent_phi_preTDR = preTDR_lumi/simu_coherent_lumi;
    for (int i = 1; i <= h_t_MC->GetNbinsX(); ++i) 
    {
        double binError = h_t_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_MC->SetBinError(i, newError);
    }
    h_t_MC->SetLineColor(kBlack);
	h_t_MC->Draw("HISTsame");

    h_t_REC_EEMC_cut->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    for (int i = 1; i <= h_t_REC_EEMC_cut->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_EEMC_cut->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_EEMC_cut->SetBinError(i, newError);
    }
	h_t_REC_EEMC_cut->SetMarkerStyle(20); // method L RECO
	h_t_REC_EEMC_cut->SetMarkerColor(kP8Blue);
	//h_t_REC_EEMC_cut->Draw("PEsame");

    h_t_REC_wRES_cut_pi12->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_wRES_cut_pi12->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_t_REC_wRES_cut_pi12->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_wRES_cut_pi12->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_wRES_cut_pi12->SetBinError(i, newError);
    }
	h_t_REC_wRES_cut_pi12->SetMarkerStyle(30);
	h_t_REC_wRES_cut_pi12->SetMarkerColor(kP8Pink);
    h_t_REC_wRES_cut_pi12->SetLineColor(kP8Pink);
	h_t_REC_wRES_cut_pi12->Draw("PEsame");

    h_t_REC_wRES_cut_pi6->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_wRES_cut_pi6->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/6)));
    for (int i = 1; i <= h_t_REC_wRES_cut_pi6->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_wRES_cut_pi6->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_wRES_cut_pi6->SetBinError(i, newError);
    }
	h_t_REC_wRES_cut_pi6->SetMarkerStyle(4);
	h_t_REC_wRES_cut_pi6->SetMarkerColor(kP8Blue);
    h_t_REC_wRES_cut_pi6->SetLineColor(kP8Blue);
	//h_t_REC_wRES_cut_pi6->Draw("PEsame");

    h_t_rho_proj_PID->Rebin(rebin_width);
    h_t_rho_proj_PID->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_rho_proj_PID->Scale(sigma_rho*(1/nEvents_rho_PID)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_rho_lumi_PID = nEvents_rho_PID/sigma_rho; // nb^-1 
    double ratio_rho_preTDR_PID = preTDR_lumi/simu_rho_lumi_PID;
    for (int i = 1; i <= h_t_rho_proj_PID->GetNbinsX(); ++i) 
    {
        double binError = h_t_rho_proj_PID->GetBinError(i);
        double newError = binError / sqrt(ratio_rho_preTDR_PID);
        h_t_rho_proj_PID->SetBinError(i, newError);
    }
    h_t_rho_proj_PID->SetMarkerStyle(42);
    h_t_rho_proj_PID->SetMarkerColor(kP8Blue);
    h_t_rho_proj_PID->SetLineColor(kP8Blue);
    h_t_rho_proj_PID->Draw("PEsame");

    h_t_rho_MC_PID->Rebin(rebin_width);
    h_t_rho_MC_PID->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_rho_MC_PID->Scale(sigma_rho*(1/nEvents_rho_PID)*(1/dN_dt)/rebin_width);
    for (int i = 1; i <= h_t_rho_MC_PID->GetNbinsX(); ++i) 
    {
        double binError = h_t_rho_MC_PID->GetBinError(i);
        double newError = binError / sqrt(ratio_rho_preTDR_PID);
        h_t_rho_MC_PID->SetBinError(i, newError);
    }
    h_t_rho_MC_PID->SetMarkerStyle(28);
    h_t_rho_MC_PID->SetMarkerColor(kP8Gray);
    h_t_rho_MC_PID->SetLineColor(kP8Gray);
    //h_t_rho_MC_PID->Draw("PEsame");

    h_t_rho_L_PID->Rebin(rebin_width); 
    h_t_rho_L_PID->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_rho_L_PID->Scale(sigma_rho*(1/nEvents_rho_PID)*(1/dN_dt)/rebin_width);
    for (int i = 1; i <= h_t_rho_L_PID->GetNbinsX(); ++i) 
    {
        double binError = h_t_rho_L_PID->GetBinError(i);
        double newError = binError / sqrt(ratio_rho_preTDR_PID);
        h_t_rho_L_PID->SetBinError(i, newError);
    }
    h_t_rho_L_PID->SetMarkerStyle(28);
    h_t_rho_L_PID->SetMarkerColor(kP8Gray);
    h_t_rho_L_PID->SetLineColor(kP8Gray);
    //h_t_rho_L_PID->Draw("PEsame");

    h_t_rho_proj->Rebin(rebin_width);
    h_t_rho_proj->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_rho_proj->Scale(sigma_rho*(1/nEvents_rho)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_rho_lumi = nEvents_rho/sigma_rho; // nb^-1 
    double ratio_rho_preTDR = preTDR_lumi/simu_rho_lumi;
    for (int i = 1; i <= h_t_rho_proj->GetNbinsX(); ++i) 
    {
        double binError = h_t_rho_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_rho_preTDR);
        h_t_rho_proj->SetBinError(i, newError);
    }
    h_t_rho_proj->SetMarkerStyle(28);
    h_t_rho_proj->SetMarkerColor(kP8Gray);
    h_t_rho_proj->SetLineColor(kP8Gray);
    h_t_rho_proj->Draw("PEsame");

    h_t_rho_MC->Rebin(rebin_width);
    h_t_rho_MC->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_rho_MC->Scale(sigma_rho*(1/nEvents_rho)*(1/dN_dt)/rebin_width);
    for (int i = 1; i <= h_t_rho_MC->GetNbinsX(); ++i) 
    {
        double binError = h_t_rho_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_rho_preTDR);
        h_t_rho_MC->SetBinError(i, newError);
    }
    h_t_rho_MC->SetMarkerStyle(26);
    h_t_rho_MC->SetMarkerColor(kP8Orange);
    h_t_rho_MC->SetLineColor(kP8Orange);
    //h_t_rho_MC->Draw("PEsame");

    h_t_rho_L->Rebin(rebin_width); 
    h_t_rho_L->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_rho_L->Scale(sigma_rho*(1/nEvents_rho)*(1/dN_dt)/rebin_width);
    for (int i = 1; i <= h_t_rho_L->GetNbinsX(); ++i) 
    {
        double binError = h_t_rho_L->GetBinError(i);
        double newError = binError / sqrt(ratio_rho_preTDR);
        h_t_rho_L->SetBinError(i, newError);
    }
    h_t_rho_L->SetMarkerStyle(26);
    h_t_rho_L->SetMarkerColor(kP8Orange);
    h_t_rho_L->SetLineColor(kP8Orange);
    //h_t_rho_L->Draw("PEsame");
    
	TLegend *w14_213 = new TLegend(0.63,0.65,0.75,0.85);
    w14_213->AddEntry(h_t_MC,"#phi: MC", "L");
    w14_213->AddEntry(h_t_REC_wRES_cut_pi12,"#phi: #theta_{max}= #pi/12", "P");    
    //w14_213->AddEntry(h_t_REC_EEMC_cut,"#phi: L", "P");
    //w14_213->AddEntry(h_t_rho_MC_PID,"#rho: MC", "P");
    //w14_213->AddEntry(h_t_rho_L_PID,"#rho (w. PID): L", "P");
    w14_213->AddEntry(h_t_rho_proj,"#rho (no PID): #theta_{max}= #pi/12", "P");
    w14_213->AddEntry(h_t_rho_proj_PID,"#rho (w. PID): #theta_{max}= #pi/12", "P");
    //w14_213->AddEntry(h_t_rho_proj_pi6_PID,"#rho: #theta_{max}= #pi/6", "P");
    //w14_213->AddEntry(h_t_rho_proj_pi2_PID,"#rho: #theta_{max}= #pi/2", "P");
    //w14_213->AddEntry(h_t_rho_MC,"#rho: MC", "P");
    //w14_213->AddEntry(h_t_rho_L,"#rho (no PID): L", "P");
	w14_213->SetBorderSize(0);   
    w14_213->SetFillStyle(0);
    w14_213->SetTextSize(30);
	w14_213->SetTextFont(45);
    w14_213->Draw("same");	

    TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

    // Add labels
    TLatex* ep = new TLatex(0.16, 0.85, "ePIC");
    ep->SetNDC();
    ep->SetTextFont(62);
    ep->SetTextSize(.038);
    ep->SetTextColor(kBlack);
    ep->Draw("same");

    TLatex* r421 = new TLatex(0.21, 0.85, " Simulation 25.10.2/3, 10x100 GeV");
	r421->SetNDC();
	r421->SetTextSize(30);
	r421->SetTextFont(43);
	r421->SetTextColor(kBlack);
	r421->Draw("same");

	TLatex* label_5 = new TLatex(0.16, 0.75, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->SetTextSize(30);
	label_5->SetTextFont(43);
	label_5->SetTextColor(kBlack);
	label_5->Draw("same");

    TLatex* r4211 = new TLatex(0.16, 0.8, "eAu #rightarrow e'Au'#phi");
	r4211->SetNDC();
	r4211->SetTextSize(30);
	r4211->SetTextFont(43);
	r4211->SetTextColor(kBlack);
	r4211->Draw("same"); 

    TLatex* r421111 = new TLatex(0.15, 0.15, "L_{int} = 10 fb^{-1}/A"); // preTDR lumi
	r421111->SetNDC();
	r421111->SetTextSize(25);
	r421111->SetTextFont(43);
	r421111->SetTextColor(kBlack);
	r421111->Draw("same");

    c14_213->Print("./figures/plot_t_rho_preTDR_projv2.pdf");
}

void plot_t_rho_ES()
{
    TFile* phi_t_file = TFile::Open("phi_sartre_25_10_2_eAu_10x100v5.root","READ");
	TH1D* h_t_MC = (TH1D*) phi_t_file->Get("h_t_MC");
	TH1D* h_t_REC_EEMC_cut = (TH1D*) phi_t_file->Get("h_t_REC_EEMC_cut");
	TH1D* h_t_REC_wRES_cut_pi12 = (TH1D*) phi_t_file->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_t_REC_wRES_cut_pi6 = (TH1D*) phi_t_file->Get("h_t_REC_wRES_cut_pi6");
    TH1D* h_phi_sartre_events = (TH1D*) phi_t_file->Get("h_Nevents");

	TFile* sartre_rho_PID = TFile::Open("rhoPID_sartre_25_10_3_eAu_10x100v4.root","READ");
    TH1D* h_t_rho_MC_PID = (TH1D*) sartre_rho_PID->Get("h_t_rho_MC_before");
    TH1D* h_t_rho_L_PID = (TH1D*) sartre_rho_PID->Get("h_rho_tREC_L");
    TH1D* h_t_rho_proj_PID = (TH1D*) sartre_rho_PID->Get("h_rho_tREC_cut");
    TH1D* h_t_rho_proj_pi2_PID = (TH1D*) sartre_rho_PID->Get("h_rho_tREC_cut_pi2");
    TH1D* h_t_rho_proj_pi6_PID = (TH1D*) sartre_rho_PID->Get("h_rho_tREC_cut_pi6");
    TH1D* h_rho_events_PID = (TH1D*) sartre_rho_PID->Get("h_Nevents");

    TFile* sartre_rho = TFile::Open("rho_sartre_25_10_3_eAu_10x100v4.root","READ");
    TH1D* h_t_rho_MC = (TH1D*) sartre_rho->Get("h_t_MC");
    TH1D* h_t_rho_L = (TH1D*) sartre_rho->Get("h_t_REC_EEMC_cut");
    TH1D* h_t_rho_proj = (TH1D*) sartre_rho->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_rho_events = (TH1D*) sartre_rho->Get("h_Nevents");

    double dN_dt = 0.002;
    double sigma_coherent = 459.05; // nb
    double sigma_rho = 5418.68; // nb
    double sigma_DIS = 50184.8; // nb
    double nEvents_coherent = h_phi_sartre_events->GetEntries();//6.2e6;
    double nEvents_rho_PID = h_rho_events_PID->GetEntries();//1.1e6;
    double nEvents_rho = h_rho_events->GetEntries();//1.1e6;
    double early_science_lumi = 8121.83; // nb^-1
    int rebin_width = 2;   

    TCanvas* c14_213 = new TCanvas("c14_213","c14_213",1,1,1200,800);
    gPad->SetLogy(1);
	gStyle->SetOptStat(0);
    h_t_MC->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
	h_t_MC->GetXaxis()->SetTitle("|t| [GeV/c]^{2}");    
    h_t_MC->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    double simu_coherent_lumi = nEvents_coherent/sigma_coherent; // nb^-1 
    double ratio_coherent_phi_ES = early_science_lumi/simu_coherent_lumi;
    for (int i = 1; i <= h_t_MC->GetNbinsX(); ++i) 
    {
        double binError = h_t_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_MC->SetBinError(i, newError);
    }
    h_t_MC->SetLineColor(kBlack);
	h_t_MC->Draw("HISTsame");

    h_t_REC_EEMC_cut->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    for (int i = 1; i <= h_t_REC_EEMC_cut->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_EEMC_cut->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_REC_EEMC_cut->SetBinError(i, newError);
    }
	h_t_REC_EEMC_cut->SetMarkerStyle(20); // method L RECO
	h_t_REC_EEMC_cut->SetMarkerColor(kP8Blue);
	//h_t_REC_EEMC_cut->Draw("PEsame");

    h_t_REC_wRES_cut_pi12->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_wRES_cut_pi12->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_t_REC_wRES_cut_pi12->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_wRES_cut_pi12->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_REC_wRES_cut_pi12->SetBinError(i, newError);
    }
	h_t_REC_wRES_cut_pi12->SetMarkerStyle(30);
	h_t_REC_wRES_cut_pi12->SetMarkerColor(kP8Pink);
    h_t_REC_wRES_cut_pi12->SetLineColor(kP8Pink);
	h_t_REC_wRES_cut_pi12->Draw("PEsame");

    h_t_REC_wRES_cut_pi6->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_wRES_cut_pi6->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/6)));
    for (int i = 1; i <= h_t_REC_wRES_cut_pi6->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_wRES_cut_pi6->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_REC_wRES_cut_pi6->SetBinError(i, newError);
    }
	h_t_REC_wRES_cut_pi6->SetMarkerStyle(4);
	h_t_REC_wRES_cut_pi6->SetMarkerColor(kP8Blue);
    h_t_REC_wRES_cut_pi6->SetLineColor(kP8Blue);
	//h_t_REC_wRES_cut_pi6->Draw("PEsame");

    h_t_rho_proj_PID->Rebin(rebin_width);
    h_t_rho_proj_PID->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_rho_proj_PID->Scale(sigma_rho*(1/nEvents_rho_PID)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_rho_lumi_PID = nEvents_rho_PID/sigma_rho; // nb^-1 
    double ratio_rho_ES_PID = early_science_lumi/simu_rho_lumi_PID;
    for (int i = 1; i <= h_t_rho_proj_PID->GetNbinsX(); ++i) 
    {
        double binError = h_t_rho_proj_PID->GetBinError(i);
        double newError = binError / sqrt(ratio_rho_ES_PID);
        h_t_rho_proj_PID->SetBinError(i, newError);
    }
    h_t_rho_proj_PID->SetMarkerStyle(42);
    h_t_rho_proj_PID->SetMarkerColor(kP8Blue);
    h_t_rho_proj_PID->SetLineColor(kP8Blue);
    h_t_rho_proj_PID->Draw("PEsame");

    h_t_rho_MC_PID->Rebin(rebin_width);
    h_t_rho_MC_PID->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_rho_MC_PID->Scale(sigma_rho*(1/nEvents_rho_PID)*(1/dN_dt)/rebin_width);
    for (int i = 1; i <= h_t_rho_MC_PID->GetNbinsX(); ++i) 
    {
        double binError = h_t_rho_MC_PID->GetBinError(i);
        double newError = binError / sqrt(ratio_rho_ES_PID);
        h_t_rho_MC_PID->SetBinError(i, newError);
    }
    h_t_rho_MC_PID->SetMarkerStyle(28);
    h_t_rho_MC_PID->SetMarkerColor(kP8Gray);
    h_t_rho_MC_PID->SetLineColor(kP8Gray);
    //h_t_rho_MC_PID->Draw("PEsame");

    h_t_rho_L_PID->Rebin(rebin_width); 
    h_t_rho_L_PID->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_rho_L_PID->Scale(sigma_rho*(1/nEvents_rho_PID)*(1/dN_dt)/rebin_width);
    for (int i = 1; i <= h_t_rho_L_PID->GetNbinsX(); ++i) 
    {
        double binError = h_t_rho_L_PID->GetBinError(i);
        double newError = binError / sqrt(ratio_rho_ES_PID);
        h_t_rho_L_PID->SetBinError(i, newError);
    }
    h_t_rho_L_PID->SetMarkerStyle(28);
    h_t_rho_L_PID->SetMarkerColor(kP8Gray);
    h_t_rho_L_PID->SetLineColor(kP8Gray);
    //h_t_rho_L_PID->Draw("PEsame");

    h_t_rho_proj_pi2_PID->Rebin(rebin_width); 
    h_t_rho_proj_pi2_PID->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_rho_proj_pi2_PID->Scale(sigma_rho*(1/nEvents_rho_PID)*(1/dN_dt)/rebin_width);
    for (int i = 1; i <= h_t_rho_proj_pi2_PID->GetNbinsX(); ++i) 
    {
        double binError = h_t_rho_proj_pi2_PID->GetBinError(i);
        double newError = binError / sqrt(ratio_rho_ES_PID);
        h_t_rho_proj_pi2_PID->SetBinError(i, newError);
    }
    h_t_rho_proj_pi2_PID->SetMarkerStyle(27);
    h_t_rho_proj_pi2_PID->SetMarkerColor(kP8Cyan);
    h_t_rho_proj_pi2_PID->SetLineColor(kP8Cyan);
    //h_t_rho_proj_pi2_PID->Draw("PEsame");

    h_t_rho_proj_pi6_PID->Rebin(rebin_width); 
    h_t_rho_proj_pi6_PID->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_rho_proj_pi6_PID->Scale(sigma_rho*(1/nEvents_rho_PID)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/6)));
    for (int i = 1; i <= h_t_rho_proj_pi6_PID->GetNbinsX(); ++i) 
    {
        double binError = h_t_rho_proj_pi6_PID->GetBinError(i);
        double newError = binError / sqrt(ratio_rho_ES_PID);
        h_t_rho_proj_pi6_PID->SetBinError(i, newError);
    }
    h_t_rho_proj_pi6_PID->SetMarkerStyle(32);
    h_t_rho_proj_pi6_PID->SetMarkerColor(kP8Orange);
    h_t_rho_proj_pi6_PID->SetLineColor(kP8Orange);
    //h_t_rho_proj_pi6_PID->Draw("PEsame");

    h_t_rho_proj->Rebin(rebin_width);
    h_t_rho_proj->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_rho_proj->Scale(sigma_rho*(1/nEvents_rho)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_rho_lumi = nEvents_rho/sigma_rho; // nb^-1 
    double ratio_rho_ES = early_science_lumi/simu_rho_lumi;
    for (int i = 1; i <= h_t_rho_proj->GetNbinsX(); ++i) 
    {
        double binError = h_t_rho_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_rho_ES);
        h_t_rho_proj->SetBinError(i, newError);
    }
    h_t_rho_proj->SetMarkerStyle(28);
    h_t_rho_proj->SetMarkerColor(kP8Gray);
    h_t_rho_proj->SetLineColor(kP8Gray);
    h_t_rho_proj->Draw("PEsame");

    h_t_rho_MC->Rebin(rebin_width);
    h_t_rho_MC->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_rho_MC->Scale(sigma_rho*(1/nEvents_rho)*(1/dN_dt)/rebin_width);
    for (int i = 1; i <= h_t_rho_MC->GetNbinsX(); ++i) 
    {
        double binError = h_t_rho_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_rho_ES);
        h_t_rho_MC->SetBinError(i, newError);
    }
    h_t_rho_MC->SetMarkerStyle(26);
    h_t_rho_MC->SetMarkerColor(kP8Orange);
    h_t_rho_MC->SetLineColor(kP8Orange);
    //h_t_rho_MC->Draw("PEsame");

    h_t_rho_L->Rebin(rebin_width); 
    h_t_rho_L->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_rho_L->Scale(sigma_rho*(1/nEvents_rho)*(1/dN_dt)/rebin_width);
    for (int i = 1; i <= h_t_rho_L->GetNbinsX(); ++i) 
    {
        double binError = h_t_rho_L->GetBinError(i);
        double newError = binError / sqrt(ratio_rho_ES);
        h_t_rho_L->SetBinError(i, newError);
    }
    h_t_rho_L->SetMarkerStyle(26);
    h_t_rho_L->SetMarkerColor(kP8Orange);
    h_t_rho_L->SetLineColor(kP8Orange);
    //h_t_rho_L->Draw("PEsame");

	TLegend *w14_213 = new TLegend(0.63,0.65,0.75,0.85);
    w14_213->AddEntry(h_t_MC,"#phi: MC", "L");
    w14_213->AddEntry(h_t_REC_wRES_cut_pi12,"#phi: #theta_{max}= #pi/12", "P");    
    //w14_213->AddEntry(h_t_REC_EEMC_cut,"#phi: L", "P");
    //w14_213->AddEntry(h_t_rho_MC_PID,"#rho: MC", "P");
    //w14_213->AddEntry(h_t_rho_L_PID,"#rho (w. PID): L", "P");
    w14_213->AddEntry(h_t_rho_proj,"#rho (no PID): #theta_{max}= #pi/12", "P");
    w14_213->AddEntry(h_t_rho_proj_PID,"#rho (w. PID): #theta_{max}= #pi/12", "P");
    //w14_213->AddEntry(h_t_rho_proj_pi6_PID,"#rho: #theta_{max}= #pi/6", "P");
    //w14_213->AddEntry(h_t_rho_proj_pi2_PID,"#rho: #theta_{max}= #pi/2", "P");
    //w14_213->AddEntry(h_t_rho_MC,"#rho: MC", "P");
    //w14_213->AddEntry(h_t_rho_L,"#rho (no PID): L", "P");
	w14_213->SetBorderSize(0);   
    w14_213->SetFillStyle(0);
    w14_213->SetTextSize(30);
	w14_213->SetTextFont(45);
    w14_213->Draw("same");	

    TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

    // Add labels
    TLatex* ep = new TLatex(0.16, 0.85, "ePIC");
    ep->SetNDC();
    ep->SetTextFont(62);
    ep->SetTextSize(.038);
    ep->SetTextColor(kBlack);
    ep->Draw("same");

    TLatex* r421 = new TLatex(0.21, 0.85, " Simulation 25.10.2/3, 10x100 GeV");
	r421->SetNDC();
	r421->SetTextSize(30);
	r421->SetTextFont(43);
	r421->SetTextColor(kBlack);
	r421->Draw("same");

	TLatex* label_5 = new TLatex(0.16, 0.75, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->SetTextSize(30);
	label_5->SetTextFont(43);
	label_5->SetTextColor(kBlack);
	label_5->Draw("same");

    TLatex* r4211 = new TLatex(0.16, 0.8, "eAu #rightarrow e'Au'#phi");
	r4211->SetNDC();
	r4211->SetTextSize(30);
	r4211->SetTextFont(43);
	r4211->SetTextColor(kBlack);
	r4211->Draw("same"); 

    TLatex* r42111 = new TLatex(0.15, 0.15, "L_{int} = 1.6 fb^{-1}/A"); // ES lumi
	r42111->SetNDC();
	r42111->SetTextSize(25);
	r42111->SetTextFont(43);
	r42111->SetTextColor(kBlack);
	r42111->Draw("same");

    c14_213->Print("./figures/plot_t_rho_preTDR_projv2.pdf");
}

void plot_background_preTDR()
{
	TFile* phi_t_file = TFile::Open("phi_sartre_25_10_2_eAu_10x100v5.root","READ");
	TH1D* h_t_MC = (TH1D*) phi_t_file->Get("h_t_MC");
	TH1D* h_t_REC_EEMC_cut = (TH1D*) phi_t_file->Get("h_t_REC_EEMC_cut");
	TH1D* h_t_REC_wRES_cut_pi12 = (TH1D*) phi_t_file->Get("h_t_REC_wRES_cut_pi12");
    TH1F* h_phi_sartre_events = (TH1F*) phi_t_file->Get("h_Nevents");

	TFile* sartre_rho = TFile::Open("rho_sartre_25_10_3_eAu_10x100v4.root","READ");
    TH1D* h_t_rho_MC = (TH1D*) sartre_rho->Get("h_t_MC");
    TH1D* h_t_rho_L = (TH1D*) sartre_rho->Get("h_t_REC_EEMC_cut");
    TH1D* h_t_rho_proj = (TH1D*) sartre_rho->Get("h_t_REC_wRES_cut_pi12");
    TH1F* h_rho_events = (TH1F*) sartre_rho->Get("h_Nevents");

	TFile* incoherent_phi_noVetoes = TFile::Open("phi_beagle_25_10_3_eAu_noVetoes_10x100v4.root","READ");
    TH1D* h_incoherent_phi_noVetoes_MC = (TH1D*) incoherent_phi_noVetoes->Get("h_t_MC");
    TH1D* h_incoherent_phi_noVetoes_L = (TH1D*) incoherent_phi_noVetoes->Get("h_t_REC_EEMC_cut");
    TH1D* h_incoherent_phi_noVetoes_proj = (TH1D*) incoherent_phi_noVetoes->Get("h_t_REC_wRES_cut_pi12");
    TH1F* h_incoherent_events = (TH1F*) incoherent_phi_noVetoes->Get("h_Nevents");

	TFile* DIS_noVetoes = TFile::Open("DIS_beagle_25_10_2_eAu_noVetoes_10x100v4.root","READ");
    TH1D* h_DIS_noVetoes_MC = (TH1D*) DIS_noVetoes->Get("h_t_MC");
    TH1D* h_DIS_noVetoes_L = (TH1D*) DIS_noVetoes->Get("h_t_REC_EEMC_cut");
    TH1D* h_DIS_noVetoes_proj = (TH1D*) DIS_noVetoes->Get("h_t_REC_wRES_cut_pi12");
    TH1F* h_DIS_events = (TH1F*) DIS_noVetoes->Get("h_Nevents");

    double dN_dt = 0.002;
    double sigma_coherent = 459.05; // nb
    double sigma_incoherent = 261.95; // nb
    double sigma_rho = 5418.68; // nb
    double sigma_DIS = 50184.8; // nb
    double nEvents_coherent = h_phi_sartre_events->GetEntries();
    double nEvents_rho = h_rho_events->GetEntries();
	double nEvents_DIS = h_DIS_events->GetEntries();
    double nEvents_incoherent = h_incoherent_events->GetEntries();
    double preTDR_lumi = 50761.4; // nb^-1
    int rebin_width = 2;   

    TCanvas* c14_213 = new TCanvas("c14_213","c14_213",1,1,1200,800);
    gPad->SetLogy(1);
	gStyle->SetOptStat(0);
 	h_t_MC->GetYaxis()->SetTitle("d#sigma/d|t| [nb(GeV/c)^{2}]");
	h_t_MC->GetXaxis()->SetTitle("|t| [GeV/c]^{2}");
    h_t_MC->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    double simu_coherent_lumi = nEvents_coherent/sigma_coherent; // nb^-1 
    double ratio_coherent_phi_preTDR = preTDR_lumi/simu_coherent_lumi;
    for (int i = 1; i <= h_t_MC->GetNbinsX(); ++i) 
    {
        double binError = h_t_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_MC->SetBinError(i, newError);
    }
    h_t_MC->SetLineColor(kBlack);
	h_t_MC->Draw("HISTsame");
    
    h_t_REC_EEMC_cut->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    for (int i = 1; i <= h_t_REC_EEMC_cut->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_EEMC_cut->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_EEMC_cut->SetBinError(i, newError);
    }
	h_t_REC_EEMC_cut->SetMarkerStyle(20); // method L RECO
	h_t_REC_EEMC_cut->SetMarkerColor(kP8Blue);
	//h_t_REC_EEMC_cut->Draw("Psame");

    h_t_REC_wRES_cut_pi12->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_wRES_cut_pi12->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_t_REC_wRES_cut_pi12->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_wRES_cut_pi12->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_wRES_cut_pi12->SetBinError(i, newError);
    }
	h_t_REC_wRES_cut_pi12->SetMarkerStyle(30);
	h_t_REC_wRES_cut_pi12->SetMarkerColor(kP8Pink);
    h_t_REC_wRES_cut_pi12->SetLineColor(kP8Pink);
	h_t_REC_wRES_cut_pi12->Draw("Psame");

    h_t_rho_proj->Rebin(rebin_width);
    h_t_rho_proj->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_rho_proj->Scale(sigma_rho*(1/nEvents_rho)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_rho_lumi = nEvents_rho/sigma_rho; // nb^-1 
    double ratio_rho_preTDR = preTDR_lumi/simu_rho_lumi;
    for (int i = 1; i <= h_t_rho_proj->GetNbinsX(); ++i) 
    {
        double binError = h_t_rho_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_rho_preTDR);
        h_t_rho_proj->SetBinError(i, newError);
    }
    h_t_rho_proj->SetMarkerStyle(28);
    h_t_rho_proj->SetMarkerColor(kP8Gray);
    h_t_rho_proj->SetLineColor(kP8Gray);
    h_t_rho_proj->Draw("PEsame");

	h_DIS_noVetoes_proj->Rebin(rebin_width);
    h_DIS_noVetoes_proj->Scale(sigma_DIS*(1/nEvents_DIS)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_DIS_lumi = nEvents_DIS/sigma_DIS; // nb^-1 
    double ratio_DIS_preTDR = preTDR_lumi/simu_DIS_lumi;
    for (int i = 1; i <= h_DIS_noVetoes_proj->GetNbinsX(); ++i) 
    {
        double binError = h_DIS_noVetoes_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_DIS_preTDR);
        h_DIS_noVetoes_proj->SetBinError(i, newError);
    }
    h_DIS_noVetoes_proj->SetMarkerStyle(2);
    h_DIS_noVetoes_proj->SetMarkerColor(kP8Cyan);
    h_DIS_noVetoes_proj->SetLineColor(kP8Cyan);
    h_DIS_noVetoes_proj->Draw("PEsame");

    h_incoherent_phi_noVetoes_proj->Rebin(rebin_width);
    h_incoherent_phi_noVetoes_proj->Scale(sigma_incoherent*(1/nEvents_incoherent)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_incoherent_lumi = nEvents_incoherent/sigma_incoherent; // nb^-1 
    double ratio_incoherent_preTDR = preTDR_lumi/simu_incoherent_lumi;
    for (int i = 1; i <= h_incoherent_phi_noVetoes_proj->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_noVetoes_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_preTDR);
        h_incoherent_phi_noVetoes_proj->SetBinError(i, newError);
    }
    h_incoherent_phi_noVetoes_proj->SetMarkerStyle(21);
    h_incoherent_phi_noVetoes_proj->SetMarkerColor(kP8Orange);
    h_incoherent_phi_noVetoes_proj->SetLineColor(kP8Orange);
    h_incoherent_phi_noVetoes_proj->Draw("PEsame");
	
    TLegend *w14_213 = new TLegend(0.6,0.6,0.72,0.88);
    w14_213->AddEntry(h_t_MC,"#phi: MC", "L");
    w14_213->AddEntry(h_t_REC_wRES_cut_pi12,"#phi: #theta_{max}= #pi/12", "P");
	w14_213->AddEntry(h_DIS_noVetoes_proj,"DIS", "P");
    w14_213->AddEntry(h_incoherent_phi_noVetoes_proj,"Incoherent #phi: #theta_{max}= #pi/12", "P");    
    //w14_213->AddEntry(h_t_REC_EEMC_cut,"#phi: L","P");
    w14_213->AddEntry(h_t_rho_proj,"#rho: #theta_{max}= #pi/12", "P");
	w14_213->SetBorderSize(0);   
    w14_213->SetFillStyle(0);
    w14_213->SetTextSize(30);
	w14_213->SetTextFont(45);
    w14_213->Draw("same");	
	
    TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

    // Add labels
    TLatex* ep = new TLatex(0.16, 0.85, "ePIC");
    ep->SetNDC();
    ep->SetTextFont(62);
    ep->SetTextSize(.038);
    ep->SetTextColor(kBlack);
    ep->Draw("same");

    TLatex* r421 = new TLatex(0.21, 0.85, " Simulation 25.10.2/3, 10x100 GeV");
	r421->SetNDC();
	r421->SetTextSize(30);
	r421->SetTextFont(43);
	r421->SetTextColor(kBlack);
	r421->Draw("same");

	TLatex* label_5 = new TLatex(0.16, 0.75, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->SetTextSize(30);
	label_5->SetTextFont(43);
	label_5->SetTextColor(kBlack);
	label_5->Draw("same");

    TLatex* r4211 = new TLatex(0.16, 0.8, "eAu #rightarrow e'Au'#phi");
	r4211->SetNDC();
	r4211->SetTextSize(30);
	r4211->SetTextFont(43);
	r4211->SetTextColor(kBlack);
	r4211->Draw("same");

    TLatex* r421111 = new TLatex(0.15, 0.15, "L_{int} = 10 fb^{-1}/A"); // preTDR lumi
	r421111->SetNDC();
	r421111->SetTextSize(25);
	r421111->SetTextFont(43);
	r421111->SetTextColor(kBlack);
	r421111->Draw("same");

    c14_213->Print("./figures/plot_t_backgroundv5_preTDR.pdf");

	//TFile* fout_preTDR = new TFile("plot_t_background_preTDR.root","RECREATE");
    //c14_213->Write();
    //h_t_MC->Write();
    //h_t_REC_EEMC_cut->Write();
    //h_t_REC_wRES_cut_pi12->Write();
    //h_DIS_noVetoes_proj->Write();
    //h_incoherent_phi_noVetoes_proj->Write();
    //h_t_rho_proj->Write();
    //fout_preTDR->Close();
}

void plot_background_ES()
{
	TFile* phi_t_file = TFile::Open("phi_sartre_25_10_2_eAu_10x100v5.root","READ");
	TH1D* h_t_MC = (TH1D*) phi_t_file->Get("h_t_MC");
	TH1D* h_t_REC_EEMC_cut = (TH1D*) phi_t_file->Get("h_t_REC_EEMC_cut");
	TH1D* h_t_REC_wRES_cut_pi12 = (TH1D*) phi_t_file->Get("h_t_REC_wRES_cut_pi12");
    TH1F* h_phi_sartre_events = (TH1F*) phi_t_file->Get("h_Nevents");

	TFile* sartre_rho = TFile::Open("rho_sartre_25_10_3_eAu_10x100v4.root","READ");
    TH1D* h_t_rho_MC = (TH1D*) sartre_rho->Get("h_t_MC");
    TH1D* h_t_rho_L = (TH1D*) sartre_rho->Get("h_t_REC_EEMC_cut");
    TH1D* h_t_rho_proj = (TH1D*) sartre_rho->Get("h_t_REC_wRES_cut_pi12");
    TH1F* h_rho_events = (TH1F*) sartre_rho->Get("h_Nevents");

	TFile* incoherent_phi_noVetoes = TFile::Open("phi_beagle_25_10_3_eAu_noVetoes_10x100v4.root","READ");
    TH1D* h_incoherent_phi_noVetoes_MC = (TH1D*) incoherent_phi_noVetoes->Get("h_t_MC");
    TH1D* h_incoherent_phi_noVetoes_L = (TH1D*) incoherent_phi_noVetoes->Get("h_t_REC_EEMC_cut");
    TH1D* h_incoherent_phi_noVetoes_proj = (TH1D*) incoherent_phi_noVetoes->Get("h_t_REC_wRES_cut_pi12");
    TH1F* h_incoherent_events = (TH1F*) incoherent_phi_noVetoes->Get("h_Nevents");

	TFile* DIS_noVetoes = TFile::Open("DIS_beagle_25_10_2_eAu_noVetoes_10x100v4.root","READ");
    TH1D* h_DIS_noVetoes_MC = (TH1D*) DIS_noVetoes->Get("h_t_MC");
    TH1D* h_DIS_noVetoes_L = (TH1D*) DIS_noVetoes->Get("h_t_REC_EEMC_cut");
    TH1D* h_DIS_noVetoes_proj = (TH1D*) DIS_noVetoes->Get("h_t_REC_wRES_cut_pi12");
    TH1F* h_DIS_events = (TH1F*) DIS_noVetoes->Get("h_Nevents");

    double dN_dt = 0.002;
    double sigma_coherent = 459.05; // nb
    double sigma_incoherent = 261.95; // nb
    double sigma_rho = 5418.68; // nb
    double sigma_DIS = 50184.8; // nb
    double nEvents_coherent = h_phi_sartre_events->GetEntries();
    double nEvents_rho = h_rho_events->GetEntries();
	double nEvents_DIS = h_DIS_events->GetEntries();
    double nEvents_incoherent = h_incoherent_events->GetEntries();
    double early_science_lumi = 8121.83; // nb^-1
    int rebin_width = 2;   

    TCanvas* c14_213 = new TCanvas("c14_213","c14_213",1,1,1200,800);
    gPad->SetLogy(1);
	gStyle->SetOptStat(0);
 	h_t_MC->GetYaxis()->SetTitle("d#sigma/d|t| [nb(GeV/c)^{2}]");
	h_t_MC->GetXaxis()->SetTitle("|t| [GeV/c]^{2}");
    h_t_MC->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    double simu_coherent_lumi = nEvents_coherent/sigma_coherent; // nb^-1 
    double ratio_coherent_phi_ES = early_science_lumi/simu_coherent_lumi;
    for (int i = 1; i <= h_t_MC->GetNbinsX(); ++i) 
    {
        double binError = h_t_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_MC->SetBinError(i, newError);
    }
    h_t_MC->SetLineColor(kBlack);
	h_t_MC->Draw("HISTsame");
    
    h_t_REC_EEMC_cut->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    for (int i = 1; i <= h_t_REC_EEMC_cut->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_EEMC_cut->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_REC_EEMC_cut->SetBinError(i, newError);
    }
	h_t_REC_EEMC_cut->SetMarkerStyle(20); // method L RECO
	h_t_REC_EEMC_cut->SetMarkerColor(kP8Blue);
	//h_t_REC_EEMC_cut->Draw("Psame");

    h_t_REC_wRES_cut_pi12->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_wRES_cut_pi12->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_t_REC_wRES_cut_pi12->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_wRES_cut_pi12->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_ES);
        h_t_REC_wRES_cut_pi12->SetBinError(i, newError);
    }
	h_t_REC_wRES_cut_pi12->SetMarkerStyle(30);
	h_t_REC_wRES_cut_pi12->SetMarkerColor(kP8Pink);
    h_t_REC_wRES_cut_pi12->SetLineColor(kP8Pink);
	h_t_REC_wRES_cut_pi12->Draw("Psame");

    h_t_rho_proj->Rebin(rebin_width);
    h_t_rho_proj->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_rho_proj->Scale(sigma_rho*(1/nEvents_rho)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_rho_lumi = nEvents_rho/sigma_rho; // nb^-1 
    double ratio_rho_ES = early_science_lumi/simu_rho_lumi;
    for (int i = 1; i <= h_t_rho_proj->GetNbinsX(); ++i) 
    {
        double binError = h_t_rho_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_rho_ES);
        h_t_rho_proj->SetBinError(i, newError);
    }
    h_t_rho_proj->SetMarkerStyle(28);
    h_t_rho_proj->SetMarkerColor(kP8Gray);
    h_t_rho_proj->SetLineColor(kP8Gray);
    h_t_rho_proj->Draw("PEsame");

	h_DIS_noVetoes_proj->Rebin(rebin_width);
    h_DIS_noVetoes_proj->Scale(sigma_DIS*(1/nEvents_DIS)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_DIS_lumi = nEvents_DIS/sigma_DIS; // nb^-1 
    double ratio_DIS_ES = early_science_lumi/simu_DIS_lumi;
    for (int i = 1; i <= h_DIS_noVetoes_proj->GetNbinsX(); ++i) 
    {
        double binError = h_DIS_noVetoes_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_DIS_ES);
        h_DIS_noVetoes_proj->SetBinError(i, newError);
    }
    h_DIS_noVetoes_proj->SetMarkerStyle(2);
    h_DIS_noVetoes_proj->SetMarkerColor(kP8Cyan);
    h_DIS_noVetoes_proj->SetLineColor(kP8Cyan);
    h_DIS_noVetoes_proj->Draw("PEsame");

    h_incoherent_phi_noVetoes_proj->Rebin(rebin_width);
    h_incoherent_phi_noVetoes_proj->Scale(sigma_incoherent*(1/nEvents_incoherent)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_incoherent_lumi = nEvents_incoherent/sigma_incoherent; // nb^-1 
    double ratio_incoherent_ES = early_science_lumi/simu_incoherent_lumi;
    for (int i = 1; i <= h_incoherent_phi_noVetoes_proj->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_noVetoes_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_ES);
        h_incoherent_phi_noVetoes_proj->SetBinError(i, newError);
    }
    h_incoherent_phi_noVetoes_proj->SetMarkerStyle(21);
    h_incoherent_phi_noVetoes_proj->SetMarkerColor(kP8Orange);
    h_incoherent_phi_noVetoes_proj->SetLineColor(kP8Orange);
    h_incoherent_phi_noVetoes_proj->Draw("PEsame");
	
    TLegend *w14_213 = new TLegend(0.6,0.6,0.72,0.85);
    w14_213->AddEntry(h_t_MC,"#phi: MC", "L");
    w14_213->AddEntry(h_t_REC_wRES_cut_pi12,"#phi: #theta_{max}= #pi/12", "P");
	w14_213->AddEntry(h_DIS_noVetoes_proj,"DIS", "P");
    w14_213->AddEntry(h_incoherent_phi_noVetoes_proj,"Incoherent #phi: #theta_{max}= #pi/12", "P");    
    //w14_213->AddEntry(h_t_REC_EEMC_cut,"#phi: L","P");
    w14_213->AddEntry(h_t_rho_proj,"#rho: #theta_{max}= #pi/12", "P");
	w14_213->SetBorderSize(0);   
    w14_213->SetFillStyle(0);
    w14_213->SetTextSize(30);
	w14_213->SetTextFont(45);
    w14_213->Draw("same");	
	
    TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

    // Add labels
    TLatex* ep = new TLatex(0.16, 0.85, "ePIC");
    ep->SetNDC();
    ep->SetTextFont(62);
    ep->SetTextSize(.038);
    ep->SetTextColor(kBlack);
    ep->Draw("same");

    TLatex* r421 = new TLatex(0.21, 0.85, " Simulation 25.10.2/3, 10x100 GeV");
	r421->SetNDC();
	r421->SetTextSize(30);
	r421->SetTextFont(43);
	r421->SetTextColor(kBlack);
	r421->Draw("same");

	TLatex* label_5 = new TLatex(0.16, 0.75, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->SetTextSize(30);
	label_5->SetTextFont(43);
	label_5->SetTextColor(kBlack);
	label_5->Draw("same");

    TLatex* r4211 = new TLatex(0.16, 0.8, "eAu #rightarrow e'Au'#phi");
	r4211->SetNDC();
	r4211->SetTextSize(30);
	r4211->SetTextFont(43);
	r4211->SetTextColor(kBlack);
	r4211->Draw("same"); 

    TLatex* r42111 = new TLatex(0.15, 0.15, "L_{int} = 1.6 fb^{-1}/A"); // ES lumi
	r42111->SetNDC();
	r42111->SetTextSize(25);
	r42111->SetTextFont(43);
	r42111->SetTextColor(kBlack);
	r42111->Draw("same");

    c14_213->Print("./figures/plot_t_background_ES.pdf");

	//TFile* fout_ES = new TFile("plot_t_background_ES.root","RECREATE");
    //c14_213->Write();
    //h_t_MC->Write();
    //h_t_REC_EEMC_cut->Write();
    //h_t_REC_wRES_cut_pi12->Write();
    //h_DIS_noVetoes_proj->Write();
    //h_incoherent_phi_noVetoes_proj->Write();
    //h_t_rho_proj->Write();
    //fout_ES->Close();
}

void plot_transform()
{
	TFile* phi_t_file = TFile::Open("phi_sartre_25_10_2_eAu_10x100v5.root","READ");
	TH1D* h_t_MC = (TH1D*) phi_t_file->Get("h_t_MC");
	TH1D* h_t_REC_EEMC_cut = (TH1D*) phi_t_file->Get("h_t_REC_EEMC_cut");
	TH1D* h_t_REC_wRES_cut_pi12 = (TH1D*) phi_t_file->Get("h_t_REC_wRES_cut_pi12");


	const int nTrials = 1000; 
	const double hbarc = 0.197;
	const double t_cut = 0.2;
	const double bmin = -12;
	const double bmax = 12;
	const int noOfBins = 300;

	const char* file = "phi_sartre_25_10_2_eAu_10x100v5.root";
	TFile* input = new TFile(file);

	TH1D* hdsigmadt_MC = (TH1D*)input->Get("h_t_MC");
	TH1D* hdsigmadt_REC = (TH1D*)input->Get("h_t_REC_EEMC_cut");
	TH1D* hdsigmadt_REC_new = (TH1D*)input->Get("h_t_REC_wRES_cut_pi12");

	int nbins = hdsigmadt_MC->GetNbinsX();
	TRandom3 randGen(0); // Random seed

	TH1D* hF_b_MC_2d = new TH1D("hF_b_MC_2d", "", noOfBins, bmin, bmax);
	TH1D* hF_b_REC_2d = new TH1D("hF_b_REC_2d", "", noOfBins, bmin, bmax);
	TH1D* hF_b_REC_new_2d = new TH1D("hF_b_REC_new_2d", "", noOfBins, bmin, bmax);

	vector<vector<double>> trials_MC(noOfBins, vector<double>(nTrials));
	vector<vector<double>> trials_REC(noOfBins, vector<double>(nTrials));
	vector<vector<double>> trials_REC_new(noOfBins, vector<double>(nTrials));

	for (int trial = 0; trial < nTrials; ++trial) 
	{
    	for (int j = 1; j <= noOfBins; ++j) 
    	{
        	double b_2d = hF_b_MC_2d->GetBinCenter(j);
        	double prefactor = 1.0 / (2*TMath::Pi());

        	double F_b_MC_2d = 0, F_b_REC_2d = 0, F_b_REC_new_2d = 0;

        	for (int i = 1; i <= nbins; ++i) 
        	{
            	double tBinWidth = hdsigmadt_MC->GetBinWidth(i);
            	double t = hdsigmadt_MC->GetBinCenter(i);
            	double delta = sqrt(fabs(t));

            	// Gaussian sampling
            	double dsigmadt_MC = randGen.Gaus(hdsigmadt_MC->GetBinContent(i), sqrt(60)*hdsigmadt_MC->GetBinError(i)) / 1e7;
            	double dsigmadt_REC = randGen.Gaus(hdsigmadt_REC->GetBinContent(i), sqrt(60)*hdsigmadt_REC->GetBinError(i)) / 1e7;
            	double dsigmadt_REC_new = randGen.Gaus(hdsigmadt_REC_new->GetBinContent(i), sqrt(60)*hdsigmadt_REC_new->GetBinError(i)) / 1e7;

            	double bessel = TMath::BesselJ0(b_2d * delta / hbarc);

            	if (t > t_cut) continue;

            	double amp_MC = dsigmadt_MC > 0 ? sqrt(dsigmadt_MC) : 0;
            	double amp_REC = dsigmadt_REC > 0 ? sqrt(dsigmadt_REC) : 0;
            	double amp_REC_new = dsigmadt_REC_new > 0 ? sqrt(dsigmadt_REC_new) : 0;

            	if (t > 0.014) { amp_MC *= -1; amp_REC *= -1; amp_REC_new *= -1; }
            	if (t > 0.048) { amp_MC *= -1; amp_REC *= -1; amp_REC_new *= -1; }
            	if (t > 0.098) { amp_MC *= -1; amp_REC *= -1; amp_REC_new *= -1; }
                if (t > 0.17) { amp_MC *= -1; amp_REC *= -1; amp_REC_new *= -1; }

            	F_b_MC_2d += amp_MC * bessel * tBinWidth / 2;
            	F_b_REC_2d += amp_REC * bessel * tBinWidth / 2;
            	F_b_REC_new_2d += amp_REC_new * bessel * tBinWidth / 2;
        	}

        	F_b_MC_2d *= prefactor / hbarc;
        	F_b_REC_2d *= prefactor / hbarc;
        	F_b_REC_new_2d *= prefactor / hbarc;

        	trials_MC[j - 1][trial] = F_b_MC_2d;
        	trials_REC[j - 1][trial] = F_b_REC_2d;
        	trials_REC_new[j - 1][trial] = F_b_REC_new_2d;
    	}
	}

	for (int j = 1; j <= noOfBins; ++j) 
	{
    	auto& vec_MC = trials_MC[j - 1];
    	auto& vec_REC = trials_REC[j - 1];
    	auto& vec_REC_new = trials_REC_new[j - 1];

    	double mean_MC = TMath::Mean(nTrials, vec_MC.data());
    	double std_MC = TMath::RMS(nTrials, vec_MC.data());
    	double mean_REC = TMath::Mean(nTrials, vec_REC.data());
    	double std_REC = TMath::RMS(nTrials, vec_REC.data());
    	double mean_REC_new = TMath::Mean(nTrials, vec_REC_new.data());
    	double std_REC_new = TMath::RMS(nTrials, vec_REC_new.data());

    	hF_b_MC_2d->SetBinContent(j, mean_MC);
    	hF_b_MC_2d->SetBinError(j, std_MC);

    	hF_b_REC_2d->SetBinContent(j, mean_REC);
    	hF_b_REC_2d->SetBinError(j, std_REC);

    	hF_b_REC_new_2d->SetBinContent(j, mean_REC_new);
    	hF_b_REC_new_2d->SetBinError(j, std_REC_new);

    	double mc = hF_b_MC_2d->GetBinContent(j);
    	double rec = hF_b_REC_2d->GetBinContent(j);
    	double rec_new = hF_b_REC_new_2d->GetBinContent(j);

    	double emc = hF_b_MC_2d->GetBinError(j);
    	double erec = hF_b_REC_2d->GetBinError(j);
    	double erec_new = hF_b_REC_new_2d->GetBinError(j);

    	double ratio_rec_mc = (mc != 0) ? rec / mc : 0;
    	double ratio_rec_new_mc = (mc != 0) ? rec_new / mc : 0;
    	double error_ratio = 0;

		if (mc != 0) 
		{
    		error_ratio = ratio_rec_new_mc * sqrt(pow(erec_new / rec_new, 2) + pow(emc / mc, 2));
		}
	}

    TCanvas* c1 = new TCanvas("c1","c1",800,800);
    gPad->SetTicks();
    gPad->SetLeftMargin(0.18);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.08);
    gPad->SetTopMargin(0.12);  
    gPad->SetLogx(0);
	gStyle->SetOptStat(0);
        
    hF_b_MC_2d->Scale(1.0 / hF_b_MC_2d->Integral("width"));
    hF_b_MC_2d->GetXaxis()->SetTitle("b [fm]");
    hF_b_MC_2d->GetYaxis()->SetTitle("F(b)/#scale[0.6]{#int} F(b) db");
    hF_b_MC_2d->GetXaxis()->SetTitleOffset(1.2);
    hF_b_MC_2d->GetYaxis()->SetTitleOffset(1.5);
    hF_b_MC_2d->GetYaxis()->SetRangeUser(-0.07, 0.16);
    //hF_b_MC_2d->GetYaxis()->SetRangeUser(-5, 5);
    hF_b_MC_2d->SetLineColor(kBlack);
    hF_b_MC_2d->SetLineWidth(4);
    hF_b_MC_2d->GetXaxis()->SetLabelSize(0.04);  
    hF_b_MC_2d->GetYaxis()->SetLabelSize(0.04);
    //hF_b_MC_2d->Draw("HIST");

	hF_b_REC_2d->GetXaxis()->SetTitle("b [fm]");
    hF_b_REC_2d->GetYaxis()->SetTitle("F(b)/#scale[0.6]{#int} F(b) db");
    hF_b_REC_2d->Scale(1.0 / hF_b_REC_2d->Integral("width"));
    hF_b_REC_2d->GetYaxis()->SetRangeUser(-0.07, 0.16);
    hF_b_REC_2d->SetMarkerStyle(20); 
    hF_b_REC_2d->SetMarkerColor(kP8Blue);
    hF_b_REC_2d->SetLineColor(kP8Blue);
    //hF_b_REC_2d->Draw("PEsame"); 

	hF_b_REC_new_2d->GetXaxis()->SetTitle("b [fm]");
    hF_b_REC_new_2d->GetYaxis()->SetTitle("F(b)/#scale[0.6]{#int} F(b) db");
    hF_b_REC_new_2d->Scale(1.0 / hF_b_REC_new_2d->Integral("width"));
    hF_b_REC_new_2d->SetMarkerStyle(30); 
    hF_b_REC_new_2d->SetMarkerColor(kP8Pink);
    hF_b_REC_new_2d->SetLineColor(kP8Pink);
    //hF_b_REC_new_2d->Draw("PEsame"); 

    hF_b_REC_2d->Draw(); 
    hF_b_REC_new_2d->Draw("PEsame"); 
    hF_b_MC_2d->SetMarkerStyle(0);
    hF_b_MC_2d->Draw("Lsame");

    TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

    // Add labels
    TLatex* ep = new TLatex(0.22, 0.3, "ePIC");
    ep->SetNDC();
    ep->SetTextFont(62);
    ep->SetTextSize(0.030);
    ep->SetTextColor(kBlack);
    ep->Draw("same");

    TLatex* r421 = new TLatex(0.28, 0.3, " Simulation 25.10.2, 10x100 GeV");
	r421->SetNDC();
	r421->SetTextSize(25);
	r421->SetTextFont(43);
	r421->SetTextColor(kBlack);
	r421->Draw("same");

	TLatex* label_5 = new TLatex(0.22, 0.2, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->SetTextSize(25);
	label_5->SetTextFont(43);
	label_5->SetTextColor(kBlack);
	label_5->Draw("same");

    TLatex* r4211 = new TLatex(0.22, 0.25, "eAu #rightarrow e'Au'#phi");
	r4211->SetNDC();
	r4211->SetTextSize(25);
	r4211->SetTextFont(43);
	r4211->SetTextColor(kBlack);
	r4211->Draw("same");

    TLegend *leg = new TLegend(0.68,0.7,0.72,0.85);
    leg->AddEntry(hF_b_MC_2d, " MC", "L");
    leg->AddEntry(hF_b_REC_2d, " Method L", "PE"); 
    leg->AddEntry(hF_b_REC_new_2d, " Projection method", "PE"); 
    leg->SetBorderSize(0);      
    leg->SetFillStyle(0);       
    leg->SetTextFont(45);       
    leg->SetTextSize(20);    
    leg->Draw("same");

    c1->Print("./figures/plot_Fb_2d.pdf");

	TFile* fout14 = new TFile("plot_Fb_2d.root","RECREATE");
    c1->Write();
    hF_b_REC_2d->Write();
    hF_b_MC_2d->Write();
    hF_b_REC_new_2d->Write();
    fout14->Close();

    cout << "All done. Bye." << endl;
}

void plot_rho_mass()
{
    TFile* sartre_rho = TFile::Open("rho_sartre_25_10_3_eAu_10x100v4.root","READ");
    TH1D* h_t_rho_MC_mass = (TH1D*) sartre_rho->Get("h_VM_mass_MC");
    TH1D* h_t_rho_REC_mass = (TH1D*) sartre_rho->Get("h_VM_mass_REC");
    TH1D* h_t_rho_REC_mass_beforeCut = (TH1D*) sartre_rho->Get("h_VM_mass_REC_before");

    TFile* phi_t_file = TFile::Open("phi_sartre_25_10_2_eAu_10x100v5.root","READ");
	TH1D* h_t_phi_MC_mass = (TH1D*) phi_t_file->Get("h_VM_mass_MC");
	TH1D* h_t_phi_REC_mass = (TH1D*) phi_t_file->Get("h_VM_mass_REC");

    TCanvas* c14_213 = new TCanvas("c14_213","c14_213",1,1,1200,800);
    c14_213->Divide(1,1,0.01,0.01);
    c14_213->cd(1);
    gPad->SetLogy(1);
	gStyle->SetOptStat(0);
    //h_t_phi_MC_mass->GetYaxis()->SetRangeUser(0,1e6);
    h_t_phi_MC_mass->GetYaxis()->SetTitle("counts");
    h_t_phi_MC_mass->GetXaxis()->SetTitle("VM mass [GeV/c^{2}]");
    h_t_phi_MC_mass->SetTitle("Mass Distribution (no PID)");
    h_t_phi_MC_mass->SetLineColor(kBlack);
    h_t_phi_MC_mass->Draw();
    h_t_phi_REC_mass->SetLineColor(kP8Pink);
    h_t_phi_REC_mass->SetLineStyle(2);
    h_t_phi_REC_mass->Draw("same");
    h_t_rho_MC_mass->SetLineColor(kP8Orange);
    h_t_rho_MC_mass->SetLineStyle(3);
    h_t_rho_MC_mass->Draw("same");
    h_t_rho_REC_mass->SetLineColor(kP8Blue);
    h_t_rho_REC_mass->SetLineStyle(7);
    h_t_rho_REC_mass->Draw("same");
    h_t_rho_REC_mass_beforeCut->SetLineColor(kP8Green);
    h_t_rho_REC_mass_beforeCut->SetLineStyle(9);
    h_t_rho_REC_mass_beforeCut->Draw("same");

    TLegend *w14_213 = new TLegend(0.55,0.68,0.85,0.88);
    w14_213->AddEntry(h_t_rho_MC_mass,"#rho MC", "L");
    w14_213->AddEntry(h_t_rho_REC_mass,"#rho REC", "L");
    w14_213->AddEntry(h_t_rho_REC_mass_beforeCut,"#rho REC before mass cut", "L");    
    w14_213->AddEntry(h_t_phi_MC_mass,"#phi MC", "L");    
    w14_213->AddEntry(h_t_phi_REC_mass,"#phi REC", "L"); 
    w14_213->SetBorderSize(0);   
    w14_213->SetFillStyle(0);
    w14_213->SetTextSize(30);
	w14_213->SetTextFont(45);   
    w14_213->Draw("same");


    c14_213->Print("./figures/plot_VM_mass_noPID.pdf");

}

void plot_rho_massv2()
{
    TFile* sartre_rho = TFile::Open("rhoPID_sartre_25_10_3_eAu_10x100v4.root","READ");
    TH1D* h_rho_mass = (TH1D*) sartre_rho->Get("h_VM_mass_MC");
    TH1D* h_rho_rec_after_mass_cut = (TH1D*) sartre_rho->Get("h_VM_mass_REC");
    TH1D* h_rho_rec_before_mass_cut = (TH1D*) sartre_rho->Get("h_VM_mass_REC_before");

    TFile* sartre_rhoNoPID = TFile::Open("rho_sartre_25_10_3_eAu_10x100v4.root","READ");
    TH1D* h_rho_massNoPID = (TH1D*) sartre_rhoNoPID->Get("h_VM_mass_MC");
    TH1D* h_rho_rec_after_mass_cutNoPID = (TH1D*) sartre_rhoNoPID->Get("h_VM_mass_REC");
    TH1D* h_rho_rec_before_mass_cutNoPID = (TH1D*) sartre_rhoNoPID->Get("h_VM_mass_REC_before");

    TFile* phi_t_file = TFile::Open("phi_sartre_25_10_2_eAu_10x100v5.root","READ");
	TH1D* h_t_phi_MC_mass = (TH1D*) phi_t_file->Get("h_VM_mass_MC");
	TH1D* h_t_phi_REC_mass = (TH1D*) phi_t_file->Get("h_VM_mass_REC");

    TCanvas* c14_213 = new TCanvas("c14_213","c14_213",1,1,1200,800);
    c14_213->Divide(1,1,0.01,0.01);
    c14_213->cd(1);
    //gPad->SetLogy(1);
	gStyle->SetOptStat(0);
    //h_t_phi_MC_mass->GetYaxis()->SetRangeUser(1e-2,1e7);
    h_t_phi_MC_mass->GetYaxis()->SetTitle("counts");
    h_t_phi_MC_mass->GetXaxis()->SetTitle("VM mass [GeV/c^{2}]");
    h_t_phi_MC_mass->SetTitle("Mass Distribution (with PID)");
    h_t_phi_MC_mass->SetLineColor(kBlack);
    h_t_phi_MC_mass->Draw();
    h_t_phi_REC_mass->SetLineColor(kP8Pink);
    h_t_phi_REC_mass->SetLineStyle(2);
    h_t_phi_REC_mass->Draw("same");
    h_rho_mass->SetLineColor(kP8Orange);
    h_rho_mass->SetLineStyle(3);
    h_rho_mass->Draw("same");
    h_rho_rec_after_mass_cut->SetLineColor(kP8Blue);
    h_rho_rec_after_mass_cut->SetLineStyle(7);
    h_rho_rec_after_mass_cut->Draw("same");
    h_rho_rec_before_mass_cut->SetLineColor(kP8Green);
    h_rho_rec_before_mass_cut->SetLineStyle(9);
    h_rho_rec_before_mass_cut->Draw("same");

    TLegend *w14_213 = new TLegend(0.55,0.68,0.85,0.88);
    w14_213->AddEntry(h_rho_mass,"#rho MC", "L");
    w14_213->AddEntry(h_rho_rec_after_mass_cut,"#rho REC", "L");
    w14_213->AddEntry(h_rho_rec_before_mass_cut,"#rho REC before mass cut", "L");     
    w14_213->AddEntry(h_t_phi_MC_mass,"#phi MC", "L");    
    w14_213->AddEntry(h_t_phi_REC_mass,"#phi REC", "L"); 
    w14_213->SetBorderSize(0);   
    w14_213->SetFillStyle(0);
    w14_213->SetTextSize(30);
	w14_213->SetTextFont(45);   
    w14_213->Draw("same");


    c14_213->Print("./figures/plot_VM_mass_wPID.pdf");
}

void plot_ZDC_2D_clusters()
{
    TFile* incoherent_file_noVeto = TFile::Open("phi_beagle_25_10_3_eAu_noVetoes_10x100v4.root","READ");
	TH2D* hZDC_xy = (TH2D*) incoherent_file_noVeto->Get("hZDC_xy");

    TCanvas* c14_213 = new TCanvas("c14_213","c14_213",1,1,1200,800);
    c14_213->Divide(1,1,0.01,0.01);
    c14_213->cd(1);
    gPad->SetLogz(1);
	gStyle->SetOptStat(0);
    hZDC_xy->SetTitle("ZDC Cluster Position");
    hZDC_xy->Draw("colzsame");

    c14_213->Print("./figures/plot_ZDC_xyClusters.pdf");
}

void plot_ZDC_2D_hits()
{
    TFile* incoherent_file_noVeto = TFile::Open("phi_beagle_25_10_3_eAu_noVetoes_10x100v4.root","READ");
	TH2D* hZDC_xy = (TH2D*) incoherent_file_noVeto->Get("hZDC_hits_xy");

    TCanvas* c14_213 = new TCanvas("c14_213","c14_213",1,1,1200,800);
    c14_213->Divide(1,1,0.01,0.01);
    c14_213->cd(1);
    gPad->SetLogz(1);
	gStyle->SetOptStat(0);
    hZDC_xy->GetXaxis()->SetRangeUser(-1200,-600);
    hZDC_xy->SetTitle("ZDC Hits Position");
    hZDC_xy->Draw("colzsame");

    c14_213->Print("./figures/plot_ZDC_2d_hits.pdf");
}

void plot_ZDC_energy()
{
    TFile* incoherent_file_noVeto = TFile::Open("phi_beagle_25_10_3_eAu_noVetoes_10x100v4.root","READ");
	TH1D* hZDC_energy_noVeto = (TH1D*) incoherent_file_noVeto->Get("hZDC_energy");

    TFile* incoherent_file_noVeto_before = TFile::Open("phi_beagle_25_10_3_eAu_noVetoes_10x100v3.root","READ");
	TH1D* hZDC_energy_noVeto_before = (TH1D*) incoherent_file_noVeto_before->Get("hZDC_energy");

    TCanvas* c14_213 = new TCanvas("c14_213","c14_213",1,1,1200,800);
    c14_213->Divide(1,1,0.01,0.01);
    c14_213->cd(1);
    gPad->SetLogy(1);
	gStyle->SetOptStat(0);
    hZDC_energy_noVeto->GetYaxis()->SetRangeUser(1e2,1e5);
    hZDC_energy_noVeto->SetLineColor(kBlack);
    hZDC_energy_noVeto->Draw();
    hZDC_energy_noVeto_before->SetLineColor(kP8Pink);
    hZDC_energy_noVeto_before->SetMarkerStyle(2);
    hZDC_energy_noVeto_before->Draw("same");

    TLegend *w14_213 = new TLegend(0.5,0.78,0.75,0.89);
    w14_213->AddEntry(hZDC_energy_noVeto,"new", "L");
    w14_213->AddEntry(hZDC_energy_noVeto_before,"old", "L");
    w14_213->SetBorderSize(0);   
    w14_213->SetFillStyle(0);
    w14_213->SetTextSize(30);
	w14_213->SetTextFont(45);   
    w14_213->Draw("same");

    c14_213->Print("./figures/plot_ZDC_energy.pdf");
}

void plot_ZDC_hits()
{
    TFile* incoherent_file_noVeto = TFile::Open("phi_beagle_25_10_3_eAu_noVetoes_10x100v4.root","READ");
	TH1D* hZDC_hits_noVeto = (TH1D*) incoherent_file_noVeto->Get("hZDC_hits");

    TFile* incoherent_file_noVeto_before = TFile::Open("phi_beagle_25_10_3_eAu_noVetoes_10x100v3.root","READ");
	TH1D* hZDC_hits_noVeto_before = (TH1D*) incoherent_file_noVeto_before->Get("hZDC_hits");

    TCanvas* c14_213 = new TCanvas("c14_213","c14_213",1,1,1200,800);
    gPad->SetLogy(1);
	gStyle->SetOptStat(0);
    hZDC_hits_noVeto->SetLineColor(kBlack);
    hZDC_hits_noVeto->Draw();
    hZDC_hits_noVeto_before->SetLineColor(kP8Pink);
    hZDC_hits_noVeto_before->SetMarkerStyle(2);
    hZDC_hits_noVeto_before->Draw("same");

    TLegend *w14_213 = new TLegend(0.5,0.78,0.75,0.89);
    w14_213->AddEntry(hZDC_hits_noVeto,"new", "L");
    w14_213->AddEntry(hZDC_hits_noVeto_before,"old", "L");
    w14_213->SetBorderSize(0);   
    w14_213->SetFillStyle(0);
    w14_213->SetTextSize(30);
	w14_213->SetTextFont(45);   
    w14_213->Draw("same");

    c14_213->Print("./figures/plot_ZDC_hits.pdf");
}

void plot_EEMC_energy()
{
    TFile* incoherent_file_noVeto = TFile::Open("phi_beagle_25_10_3_eAu_noVetoes_10x100.root","READ");
	TH1D* hEEMC_energy_noVeto = (TH1D*) incoherent_file_noVeto->Get("hEEMC_energy");

	TFile* incoherent_file_ZDCveto = TFile::Open("phi_beagle_25_10_3_eAu_zdcVetoes_10x100v2.root","READ");
    TH1D* hEEMC_energy_veto = (TH1D*) incoherent_file_ZDCveto->Get("hEEMC_energy");

    TCanvas* c14_213 = new TCanvas("c14_213","c14_213",1,1,1200,800);
    c14_213->Divide(1,1,0.01,0.01);
    c14_213->cd(1);
    gPad->SetLogy(1);
	gStyle->SetOptStat(0);
    hEEMC_energy_noVeto->SetLineColor(kBlack);
    hEEMC_energy_noVeto->Draw();
    hEEMC_energy_veto->SetLineColor(kP8Pink);
    hEEMC_energy_veto->SetLineStyle(2);
    hEEMC_energy_veto->Draw("same");

    TLegend *w14_213 = new TLegend(0.55,0.68,0.85,0.88);
    w14_213->AddEntry(hEEMC_energy_noVeto,"no ZDC veto", "L");
    w14_213->AddEntry(hEEMC_energy_veto,"with ZDC veto", "L");
    w14_213->SetBorderSize(0);   
    w14_213->SetFillStyle(0);
    w14_213->SetTextSize(30);
	w14_213->SetTextFont(45);   
    w14_213->Draw("same");

    c14_213->Print("./figures/plot_EEMC_energy.pdf");
}

void plot_EEMC_2D_noVeto()
{
    TFile* incoherent_file_noVeto = TFile::Open("phi_beagle_25_10_3_eAu_noVetoes_10x100v2.root","READ");
	TH2D* hEEMC_xy = (TH2D*) incoherent_file_noVeto->Get("hEEMC_xy");

    TCanvas* c14_213 = new TCanvas("c14_213","c14_213",1,1,1200,800);
    c14_213->Divide(1,1,0.01,0.01);
    c14_213->cd(1);
    gPad->SetLogz(1);
	gStyle->SetOptStat(0);
    hEEMC_xy->SetTitle("EEMC No ZDC Veto");
    hEEMC_xy->Draw("colzsame");

    c14_213->Print("./figures/plot_EEMC_xy_noZDCveto.pdf");

}

void plot_EEMC_2D_ZDCveto()
{
	TFile* incoherent_file_ZDCveto = TFile::Open("phi_beagle_25_10_3_eAu_zdcVetoes_10x100v2.root","READ");
    TH1D* hEEMC_xy = (TH1D*) incoherent_file_ZDCveto->Get("hEEMC_xy");

    TCanvas* c14_213 = new TCanvas("c14_213","c14_213",1,1,1200,800);
    c14_213->Divide(1,1,0.01,0.01);
    c14_213->cd(1);
    gPad->SetLogz(1);
    hEEMC_xy->SetTitle("EEMC With ZDC Veto");
    hEEMC_xy->Draw("colzsame");

    c14_213->Print("./figures/plot_EEMC_xy_ZDCveto.pdf");
}

void plot_EEMC_hits()
{
    TFile* incoherent_file_noVeto = TFile::Open("phi_beagle_25_10_3_eAu_noVetoes_10x100v2.root","READ");
	TH1D* hEEMC_hits_noVeto = (TH1D*) incoherent_file_noVeto->Get("hEEMC_hits");

	TFile* incoherent_file_ZDCveto = TFile::Open("phi_beagle_25_10_3_eAu_zdcVetoes_10x100v2.root","READ");
    TH1D* hEEMC_hits_veto = (TH1D*) incoherent_file_ZDCveto->Get("hEEMC_hits");

    TCanvas* c14_213 = new TCanvas("c14_213","c14_213",1,1,1200,800);
    c14_213->Divide(1,1,0.01,0.01);
    c14_213->cd(1);
    gPad->SetLogy(1);
	gStyle->SetOptStat(0);
    hEEMC_hits_noVeto->SetLineColor(kBlack);
    hEEMC_hits_noVeto->Draw();
    hEEMC_hits_veto->SetLineColor(kP8Pink);
    hEEMC_hits_veto->SetLineStyle(2);
    hEEMC_hits_veto->Draw("same");

    TLegend *w14_213 = new TLegend(0.55,0.68,0.85,0.88);
    w14_213->AddEntry(hEEMC_hits_noVeto,"EEMC no ZDC veto", "L");
    w14_213->AddEntry(hEEMC_hits_veto,"EEMC with ZDC veto", "L");
    w14_213->SetBorderSize(0);   
    w14_213->SetFillStyle(0);
    w14_213->SetTextSize(30);
	w14_213->SetTextFont(45);   
    w14_213->Draw("same");

    c14_213->Print("./figures/plot_EEMC_hits.pdf");
}

void plot_RP_hits()
{
    TFile* incoherent_file_noVeto = TFile::Open("phi_beagle_25_10_3_eAu_noVetoes_10x100v3.root","READ");
	TH1D* hRP_hits_noVeto = (TH1D*) incoherent_file_noVeto->Get("hRP_hits");

	TFile* incoherent_file_ZDCveto = TFile::Open("phi_beagle_25_10_3_eAu_rpVetoes_10x100v3.root","READ");
    TH1D* hRP_hits_veto = (TH1D*) incoherent_file_ZDCveto->Get("hRP_hits");

    TCanvas* c14_213 = new TCanvas("c14_213","c14_213",1,1,1200,800);
    c14_213->Divide(1,1,0.01,0.01);
    c14_213->cd(1);
    gPad->SetLogy(1);
	gStyle->SetOptStat(0);
    hRP_hits_noVeto->SetTitle("RP cluster multiplicity");
    hRP_hits_noVeto->GetXaxis()->SetTitle("clusters/event");
    hRP_hits_noVeto->GetYaxis()->SetTitle("counts");
    hRP_hits_noVeto->SetLineColor(kBlack);
    hRP_hits_noVeto->Draw();
    hRP_hits_veto->SetLineColor(kP8Pink);
    hRP_hits_veto->SetLineStyle(2);
    //hRP_hits_veto->Draw("same");
   
    TLegend *w14_213 = new TLegend(0.55,0.68,0.85,0.88);
    w14_213->AddEntry(hRP_hits_noVeto,"RP hits", "L");
    w14_213->SetBorderSize(0);   
    w14_213->SetFillStyle(0);
    w14_213->SetTextSize(30);
	w14_213->SetTextFont(45);   
    //w14_213->Draw("same");

    c14_213->Print("./figures/plot_RP_hits.pdf");
}

void plot_OMD_hits()
{
    TFile* incoherent_file_noVeto = TFile::Open("phi_beagle_25_10_3_eAu_noVetoes_10x100v3.root","READ");
	TH1D* hOMD_hits_noVeto = (TH1D*) incoherent_file_noVeto->Get("hOMD_hits");

	TFile* incoherent_file_ZDCveto = TFile::Open("phi_beagle_25_10_3_eAu_zdcVetoes_10x100v3.root","READ");
    TH1D* hOMD_hits_veto = (TH1D*) incoherent_file_ZDCveto->Get("hOMD_hits");

    TCanvas* c14_213 = new TCanvas("c14_213","c14_213",1,1,1200,800);
    c14_213->Divide(1,1,0.01,0.01);
    c14_213->cd(1);
    gPad->SetLogy(1);
	gStyle->SetOptStat(0);
    hOMD_hits_noVeto->SetTitle("OMD cluster multiplicity");
    hOMD_hits_noVeto->GetXaxis()->SetTitle("clusters/event");
    hOMD_hits_noVeto->GetYaxis()->SetTitle("counts");
    hOMD_hits_noVeto->SetLineColor(kBlack);
    hOMD_hits_noVeto->Draw();
    hOMD_hits_veto->SetLineColor(kP8Pink);
    hOMD_hits_veto->SetLineStyle(2);
    //hOMD_hits_veto->Draw("same");

    TLegend *w14_213 = new TLegend(0.55,0.68,0.85,0.88);
    w14_213->AddEntry(hOMD_hits_noVeto,"OMD no ZDC veto", "L");
    w14_213->AddEntry(hOMD_hits_veto,"OMD with ZDC veto", "L");
    w14_213->SetBorderSize(0);   
    w14_213->SetFillStyle(0);
    w14_213->SetTextSize(30);
	w14_213->SetTextFont(45);   
    //w14_213->Draw("same");

    c14_213->Print("./figures/plot_OMD_hits.pdf");
}

void plot_t_2d_MC()
{
    TFile* phi_t_file = TFile::Open("phi_sartre_25_10_2_eAu_10x100v3.root","READ");
	TH1D* h_t_REC_2d = (TH1D*) phi_t_file->Get("h_t_REC_2d");

    TCanvas* c16 = new TCanvas("c16","c16",1,1,800,800);
    c16->Divide(1,1,0.01,0.01);
    c16->cd(1);
    gPad->SetLogz(1);
    h_t_REC_2d->GetYaxis()->SetLabelSize(0.03); 
    h_t_REC_2d->GetYaxis()->SetTitleOffset(1.4);
    h_t_REC_2d->GetXaxis()->SetLabelSize(0.03);
    h_t_REC_2d->GetXaxis()->SetTitle("#sqrt{|t|_{x}} [GeV/c]");
    h_t_REC_2d->GetYaxis()->SetTitle("#sqrt{|t|_{y}} [GeV/c]");
    h_t_REC_2d->Draw("colzsame");

    gStyle->SetOptStat(0);
    gPad->SetLogz(1);

    /*double length = 0.14;
    double theta  = M_PI/12; 
    double x_end  = length * sin(theta);
    double y_end  = length * cos(theta);

    TArrow *arrow = new TArrow(0, 0, x_end, y_end, 0.02, "|>");
    arrow->SetLineColor(kBlack);
    arrow->SetLineWidth(2);
    arrow->Draw();

    double radius = 0.1;
    TEllipse *arc = new TEllipse(0, 0, radius, radius, 75, 90); 
    arc->SetFillStyle(0);
    arc->SetLineColor(kBlack);
    arc->SetLineWidth(2);
    arc->Draw();

    TLatex* r421111 = new TLatex(0.14, 0.56, "#theta_{max}"); 
	r421111->SetNDC();
	r421111->SetTextSize(30);
	r421111->SetTextFont(43);
	r421111->SetTextColor(kBlack);
	r421111->Draw("same");
    gPad->Update();*/

    gStyle->SetOptStat(0);
    gPad->SetLogz(1);

    c16->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title16 = new TLatex();
    title16->SetNDC(); 
    title16->SetTextSize(0.05);
    title16->SetTextAlign(22);  
    title16->DrawLatex(0.5, 0.97, "2D |t| Distribution");  

    c16->Print("./figures/plot_t_2Dmc_noarrow.pdf");
}

void plot_t_2d_wRES()
{
    TFile* phi_t_file = TFile::Open("phi_sartre_25_10_2_eAu_10x100v3.root","READ");
	TH1D* h_t_REC_2d_wRES_cut = (TH1D*) phi_t_file->Get("h_t_REC_2d_wRES_cut");

    TCanvas* c16 = new TCanvas("c16","c16",1,1,800,800);
    c16->Divide(1,1,0.01,0.01);
    c16->cd(1);
    gPad->SetLogz(1);
    h_t_REC_2d_wRES_cut->GetYaxis()->SetLabelSize(0.03); 
    h_t_REC_2d_wRES_cut->GetYaxis()->SetTitleOffset(1.4);
    h_t_REC_2d_wRES_cut->GetXaxis()->SetLabelSize(0.03);
    h_t_REC_2d_wRES_cut->GetXaxis()->SetTitle("#sqrt{|t|_{x}} [GeV/c]");
    h_t_REC_2d_wRES_cut->GetYaxis()->SetTitle("#sqrt{|t|_{y}} [GeV/c]");
    h_t_REC_2d_wRES_cut->Draw("colzsame");

    gStyle->SetOptStat(0);
    gPad->SetLogz(1);

    /*double length = 0.14;
    double theta  = M_PI/12; 
    double x_end  = length * sin(theta);
    double y_end  = length * cos(theta);

    TArrow *arrow = new TArrow(0, 0, x_end, y_end, 0.02, "|>");
    arrow->SetLineColor(kBlack);
    arrow->SetLineWidth(2);
    arrow->Draw();

    double radius = 0.1;
    TEllipse *arc = new TEllipse(0, 0, radius, radius, 75, 90); 
    arc->SetFillStyle(0);
    arc->SetLineColor(kBlack);
    arc->SetLineWidth(2);
    arc->Draw();

    TLatex* r421111 = new TLatex(0.14, 0.56, "#theta_{max}"); 
	r421111->SetNDC();
	r421111->SetTextSize(30);
	r421111->SetTextFont(43);
	r421111->SetTextColor(kBlack);
	r421111->Draw("same");
    gPad->Update();*/

    gStyle->SetOptStat(0);
    gPad->SetLogz(1);

    c16->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title16 = new TLatex();
    title16->SetNDC(); 
    title16->SetTextSize(0.05);
    title16->SetTextAlign(22);  
    title16->DrawLatex(0.5, 0.97, "2D |t| Distribution with Detector Resolution");  

    c16->Print("./figures/plot_t_2D_withRes_noArrow.pdf");
}