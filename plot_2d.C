#include "RiceStyle.h"

using namespace std;

void plot_2d(TString filename="./output/eicrecon-sartre_coherent_phi_output.root")
{
	TFile* file = new TFile(filename);
	TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";
	if(filename=="jpsi") {vm_label="J/#psi";daug_label="e^{+}e^{-}";}
	//t distribution
	TH2D* h_t_REC_2d = (TH2D*) file->Get("h_t_REC_2d");

	TCanvas* c1 = new TCanvas("c1","c1",1,1,600,600);
	gPad->SetLogz(1);

	h_t_REC_2d->Draw();
	
	c1->Print("./figures/tests/2d.pdf");
}