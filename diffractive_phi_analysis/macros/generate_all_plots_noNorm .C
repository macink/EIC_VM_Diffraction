#include "ePICStyle.C"
#include <string>
#include <sstream>
#include <locale>
gROOT->ProcessLine("SetePICStyle()");

/*----------------------------------------------------------------
	Generates the |t| distributions for different angle
	cuts that haven't been normalized

	Input: root file
	Output: saves plots as pdfs
-----------------------------------------------------------------*/

// add commas between the thousands
string formatWithCommas(double value) 
{
    stringstream ss;
    ss.imbue(locale(""));
    ss << fixed << setprecision(0) << value;
    return ss.str();
}

using namespace std;

void generate_all_plots_noNorm(TString filename)
{
	// test
    cout << "in the macro..." << filename << endl;
    if (filename == "") {
        cerr << "Error: No filename provided!" << endl;
        return;
    }
    TFile* file = TFile::Open(filename.Data(), "READ");
    if (!file || file->IsZombie()) {
        cerr << "Error: Unable to open file: " << filename << endl;
        return;
    }
	cout << "Marker style: " << gStyle->GetMarkerStyle() << endl;

	// plot labels
	TString vm_label="#phi";
	TString angle_pi2="#pi/2";
	TString angle_pi3="#pi/3";
	TString angle_pi4="#pi/4";
	TString angle_pi6="#pi/6";
	TString angle_pi9="#pi/9";
	TString angle_pi12="#pi/12";
	TString angle_pi16="#pi/16";
	TString angle_pi20="#pi/20";
	TString angle_pi24="#pi/24";
	TString daug_label="K^{+}K^{-}";
	if(filename=="jpsi") {vm_label="J/#psi";daug_label="e^{+}e^{-}";}
	
	// t distribution
	TH1D* h_t_MC = (TH1D*) file->Get("h_t_MC"); // method E
	TH1D* h_t_REC_EEMC = (TH1D*) file->Get("h_t_REC_EEMC"); // method L
	TH1D* h_t_REC_new_method = (TH1D*) file->Get("h_t_REC_new_method"); // projection method with no resolution

	// t reco with resolution and angle cut
		// these distributions do not have a position threshold
	TH1D* h_t_REC_wRES_cut_pi2 = (TH1D*) file->Get("h_t_REC_wRES_cut_pi2");
	TH1D* h_t_REC_wRES_cut_pi3 = (TH1D*) file->Get("h_t_REC_wRES_cut_pi3");
	TH1D* h_t_REC_wRES_cut_pi4 = (TH1D*) file->Get("h_t_REC_wRES_cut_pi4");
	TH1D* h_t_REC_wRES_cut_pi6 = (TH1D*) file->Get("h_t_REC_wRES_cut_pi6");
	TH1D* h_t_REC_wRES_cut_pi9 = (TH1D*) file->Get("h_t_REC_wRES_cut_pi9");
	TH1D* h_t_REC_wRES_cut_pi12 = (TH1D*) file->Get("h_t_REC_wRES_cut_pi12");
	TH1D* h_t_REC_wRES_cut_pi16 = (TH1D*) file->Get("h_t_REC_wRES_cut_pi16");
	TH1D* h_t_REC_wRES_cut_pi20 = (TH1D*) file->Get("h_t_REC_wRES_cut_pi20");
	TH1D* h_t_REC_wRES_cut_pi24 = (TH1D*) file->Get("h_t_REC_wRES_cut_pi24");
		// these distributions are accepted within the 20 mm position threshold
	TH1D* h_t_REC_wRES_cut_pi2_cut = (TH1D*) file->Get("h_t_REC_wRES_cut_pi2_cut");
	TH1D* h_t_REC_wRES_cut_pi3_cut = (TH1D*) file->Get("h_t_REC_wRES_cut_pi3_cut");
	TH1D* h_t_REC_wRES_cut_pi4_cut = (TH1D*) file->Get("h_t_REC_wRES_cut_pi4_cut");
	TH1D* h_t_REC_wRES_cut_pi6_cut = (TH1D*) file->Get("h_t_REC_wRES_cut_pi6_cut");
	TH1D* h_t_REC_wRES_cut_pi9_cut = (TH1D*) file->Get("h_t_REC_wRES_cut_pi9_cut");
	TH1D* h_t_REC_wRES_cut_pi12_cut = (TH1D*) file->Get("h_t_REC_wRES_cut_pi12_cut");
	TH1D* h_t_REC_wRES_cut_pi16_cut = (TH1D*) file->Get("h_t_REC_wRES_cut_pi16_cut");
	TH1D* h_t_REC_wRES_cut_pi20_cut = (TH1D*) file->Get("h_t_REC_wRES_cut_pi20_cut");
	TH1D* h_t_REC_wRES_cut_pi24_cut = (TH1D*) file->Get("h_t_REC_wRES_cut_pi24_cut");

	// create canvas for each angle
	const int numCanvas = 9; 

	// angle labels
	TString angleLabels[numCanvas] = 
	{
    	angle_pi2, angle_pi3, angle_pi4, angle_pi6, angle_pi9, angle_pi12, angle_pi16, angle_pi20, angle_pi24
	};
	double angles[numCanvas] = 
	{
		M_PI/2, M_PI/3, M_PI/4, M_PI/6, M_PI/9, M_PI/12, M_PI/16, M_PI/20, M_PI/24
	};
	
	// histograms for distributions without position cut
	TH1D* h_t_REC_wRES_cut[numCanvas] = 
	{
		h_t_REC_wRES_cut_pi2,
		h_t_REC_wRES_cut_pi3,
		h_t_REC_wRES_cut_pi4,
		h_t_REC_wRES_cut_pi6,
		h_t_REC_wRES_cut_pi9,
		h_t_REC_wRES_cut_pi12,
		h_t_REC_wRES_cut_pi16,
		h_t_REC_wRES_cut_pi20,
		h_t_REC_wRES_cut_pi24
	};
	// histograms for distributions within the 20 mm position threshold
	TH1D* h_t_REC_wRES_cut_cut[numCanvas] = 
	{
		h_t_REC_wRES_cut_pi2_cut,
		h_t_REC_wRES_cut_pi3_cut,
		h_t_REC_wRES_cut_pi4_cut,
		h_t_REC_wRES_cut_pi6_cut,
		h_t_REC_wRES_cut_pi9_cut,
		h_t_REC_wRES_cut_pi12_cut,
		h_t_REC_wRES_cut_pi16_cut,
		h_t_REC_wRES_cut_pi20_cut,
		h_t_REC_wRES_cut_pi24_cut
	};
	

	// reproduce truth plot with our decomposition
	TCanvas* newMethod_canvases[numCanvas];
	for (int i=0; i<1; i++) 
	{
		TString canvasName = Form("newMethod_c%d", i+1);
		newMethod_canvases[i] = new TCanvas(canvasName, canvasName);
		gPad->SetLogy(1);
		gPad->SetTopMargin(0.08); //This is 0.05 in ePICStyle
		
		// get counts
		double integral_MC = h_t_MC->Integral();
		double integral_REC = h_t_REC_EEMC->Integral();
		double integral_REC_new_method = h_t_REC_new_method->Integral();

		// Draw histograms 
		h_t_MC->SetLineColor(kBlue);
		h_t_MC->Draw("same");
		h_t_REC_EEMC->SetMarkerStyle(20); // method L RECO
		h_t_REC_EEMC->SetMarkerColor(kBlack);
		h_t_REC_EEMC->Draw("PEsame");
		h_t_REC_new_method->SetMarkerStyle(30);
		h_t_REC_new_method->SetMarkerColor(kRed);
		h_t_REC_new_method->Draw("P same");

		// Add labels
		TLatex* r42 = new TLatex(0.15, 0.94, "eAu 10x100 GeV");
		r42->SetNDC();
		r42->Draw("same");

		TLatex* r43 = new TLatex(0.45, 0.94, "Sartre, EPIC (no res. for #theta_{max} method)");
		r43->SetNDC();
		r43->Draw("same");

		TLatex* r44 = new TLatex(0.2, 0.86, "1<Q^{2}<10 GeV^{2}, 0.01 < y < 0.85");
		r44->SetNDC();
		r44->Draw("same");

		TLatex* r44_0 = new TLatex(0.2, 0.79, Form("|y_{%s}|<3.5, |M_{inv} #minus M_{%s}| < 0.02 GeV", vm_label.Data(), vm_label.Data()));
		r44_0->SetNDC();
		r44_0->Draw("same");

		TLatex* r44_2 = new TLatex(0.18, 0.2, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
		r44_2->SetNDC();
		r44_2->Draw("same");

		// Add legend
		TLegend* w7 = new TLegend(0.31, 0.59, 0.65, 0.76);
		w7->AddEntry(h_t_MC, Form("%s_{MC} %s evts", vm_label.Data(),formatWithCommas(integral_MC).c_str()), "L");
		w7->AddEntry(h_t_REC_EEMC, Form("%s_{L: RECO} %s evts", vm_label.Data(),formatWithCommas(integral_REC).c_str()), "P");
		w7->AddEntry(h_t_REC_new_method, Form("%s_{#theta_{max}= #pi/2: RECO} %s evts", vm_label.Data(),formatWithCommas(integral_REC_new_method).c_str()), "P");
		w7->Draw("same");

		// Save figure
		newMethod_canvases[i]->Print("./figures/noNorm_new_method.pdf");
	}

	// plots with resolution and cut
	TCanvas* noNorm_canvases_resCut[numCanvas];
	for (int i=0; i<numCanvas; i++) 
	{
    	TString canvasName = Form("noNorm_resCut_c%d", i+1);
    	noNorm_canvases_resCut[i] = new TCanvas(canvasName, canvasName);
    	gPad->SetLogy(1);
		gPad->SetTopMargin(0.08); //This is 0.05 in ePICStyle
    	
		// get counts
		double integral_MC = h_t_MC->Integral();
		double integral_REC = h_t_REC_EEMC->Integral();
		double integral_REC_wRES_cut[numCanvas];
		integral_REC_wRES_cut[i] = h_t_REC_wRES_cut[i]->Integral();
		double integral_REC_wRES_cut_cut[numCanvas];
		integral_REC_wRES_cut_cut[i] = h_t_REC_wRES_cut_cut[i]->Integral();
		
		// Draw histograms with different angle cuts
		h_t_MC->SetLineColor(kBlack);
		h_t_MC->GetXaxis()->SetTitle("|t| [GeV/c]^{2}");
		h_t_MC->Draw("same");
		h_t_REC_EEMC->SetMarkerStyle(24); // method L (RECO)
		h_t_REC_EEMC->SetMarkerColor(kBlue);
		//h_t_REC_EEMC->Draw("PEsame");
		h_t_REC_wRES_cut[i]->SetMarkerStyle(30);
		h_t_REC_wRES_cut[i]->SetMarkerColor(kRed);
		//h_t_REC_wRES_cut[i]->Draw("P same");
		h_t_REC_wRES_cut_cut[i]->SetMarkerStyle(30);
		h_t_REC_wRES_cut_cut[i]->SetMarkerColor(kRed);
		h_t_REC_wRES_cut_cut[i]->Draw("PEsame");

		// Add labels
		TLatex* r42 = new TLatex(0.18, 0.94, "eAu 10x100 GeV");
		r42->SetNDC();
		r42->Draw("same");

		TLatex* r43 = new TLatex(0.6, 0.94, "Sartre, EPIC");
		r43->SetNDC();
		r43->Draw("same");

		TLatex* r44 = new TLatex(0.3, 0.86, "1<Q^{2}<10 GeV^{2}, 0.01 < y < 0.85");
		r44->SetNDC();
		r44->Draw("same");

		TLatex* r44_0 = new TLatex(0.3, 0.79, Form("|y_{%s}|<3.5, |M_{inv} #minus M_{%s}| < 0.02 GeV", vm_label.Data(), vm_label.Data()));
		r44_0->SetNDC();
		r44_0->Draw("same");

		TLatex* r44_2 = new TLatex(0.18, 0.2, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
		r44_2->SetNDC();
		r44_2->Draw("same");

		// Add legend
		TLegend* w7 = new TLegend(0.41, 0.55, 0.65, 0.72);
		w7->AddEntry(h_t_MC, Form("%s_{MC} %s evts", vm_label.Data(),formatWithCommas(integral_MC).c_str()), "L");
		//w7->AddEntry(h_t_REC_EEMC, Form("%s_{L}  %s evts", vm_label.Data(),formatWithCommas(integral_REC).c_str()), "P");
		//w7->AddEntry(h_t_REC_wRES_cut[i], Form("%s_{#theta_{max}= %s: RECO} %s evts", vm_label.Data(),angleLabels[i].Data(),formatWithCommas(integral_REC_wRES_cut[i]).c_str()), "P");
		w7->AddEntry(h_t_REC_wRES_cut_cut[i], Form("%s_{#theta_{max}= %s} %s evts", vm_label.Data(),angleLabels[i].Data(),formatWithCommas(integral_REC_wRES_cut_cut[i]).c_str()), "P");
		w7->Draw("same");

		// Save figure
		TString cleanLabel = angleLabels[i];
		cleanLabel.ReplaceAll("#", "");
		cleanLabel.ReplaceAll("/", "");
		noNorm_canvases_resCut[i]->Print(Form("./figures/noNorm_angle%s_wRES_cut.pdf",cleanLabel.Data()));
	}
}