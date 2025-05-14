#include "RiceStyle.h"
#include <cmath>

using namespace std;

void plot_2d_combined()
{
	TString filename = "/home/macink/miniconda3/envs/bnl_research/macros/eic/EICreconOutputReader/output/combined_2dhistograms.root";
	TFile* file = new TFile(filename);
	TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";
	if(filename=="jpsi") {vm_label="J/#psi";daug_label="e^{+}e^{-}";}
	//t distribution
	TH2D* h_t_REC_2d = dynamic_cast<TH2D*>(file->Get("h_t_REC_2d_combined"));

	TCanvas* c1 = new TCanvas("c1","c1",1,1,600,600);
	gPad->SetLogz(1);
	
	h_t_REC_2d->Draw();

	c1->Print("/home/macink/miniconda3/envs/bnl_research/macros/eic/EICreconOutputReader/figures/combined/combined_2d.pdf");
}
	
