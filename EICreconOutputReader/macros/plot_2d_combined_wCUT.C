#include "RiceStyle.h"
#include <cmath>

using namespace std;

void plot_2d_combined_wCUT()
{
	TString filename = "/home/macink/miniconda3/envs/bnl_research/macros/eic/EICreconOutputReader/output/combined_2dhistograms_wCUT.root";
	TFile* file = new TFile(filename);
	TString vm_label="#phi";
	TString angle = "#pi/12";
	TString daug_label="K^{+}K^{-}";
	if(filename=="jpsi") {vm_label="J/#psi";daug_label="e^{+}e^{-}";}
	//t distribution
	TH2D* h_t_REC_2d_wCUT = dynamic_cast<TH2D*>(file->Get("h_t_REC_2d_wCUT_combined"));
	
	TCanvas* c1 = new TCanvas("c1","c1",1,1,600,600);
	gPad->SetLogz(1);
	
	h_t_REC_2d_wCUT->Draw();

	TLatex* r46 = new TLatex(0.45, 0.92, "weight: #pi/2/#theta_{Max}, #theta_{max}= "+angle);
	r46->SetNDC();
	r46->SetTextSize(15);
	r46->SetTextFont(43);
	r46->SetTextColor(kBlack);
	r46->Draw("same");

	c1->Print("/home/macink/miniconda3/envs/bnl_research/macros/eic/EICreconOutputReader/figures/combined/combined_2d_wCUT.pdf");
}
	
