#include "RiceStyle.h"

using namespace std;

void generate_all_plots(TString filename)
{
	//test
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

	//plot labels
	TString vm_label="#phi";
	TString angle_pi2="#pi/2";
	TString angle_pi3="#pi/3";
	TString angle_pi4="#pi/4";
	TString angle_pi6="#pi/6";
	TString angle_pi9="#pi/9";
	TString angle_pi12="#pi/12";
	TString daug_label="K^{+}K^{-}";
	if(filename=="jpsi") {vm_label="J/#psi";daug_label="e^{+}e^{-}";}

	//t distribution
	TH1D* h_t_MC = (TH1D*) file->Get("h_t_MC");
	TH1D* h_t_REC = (TH1D*) file->Get("h_t_REC");
	TH1D* h_t_REC_new_method = (TH1D*) file->Get("h_t_REC_new_method");
	TH1D* h_t_REC_wRES_cut_pi2 = (TH1D*) file->Get("h_t_REC_wRES_cut_pi2");
	TH1D* h_t_REC_wRES_cut_pi3 = (TH1D*) file->Get("h_t_REC_wRES_cut_pi3");
	TH1D* h_t_REC_wRES_cut_pi4 = (TH1D*) file->Get("h_t_REC_wRES_cut_pi4");
	TH1D* h_t_REC_wRES_cut_pi6 = (TH1D*) file->Get("h_t_REC_wRES_cut_pi6");
	TH1D* h_t_REC_wRES_cut_pi9 = (TH1D*) file->Get("h_t_REC_wRES_cut_pi9");
	TH1D* h_t_REC_wRES_cut_pi12 = (TH1D*) file->Get("h_t_REC_wRES_cut_pi12");
	TH1D* h_t_REC_wRES = (TH1D*) file->Get("h_t_REC_wRES");
	TH1D* h_t_REC_wCUT_pi2 = (TH1D*) file->Get("h_t_REC_wCUT_pi2");
	TH1D* h_t_REC_wCUT_pi3 = (TH1D*) file->Get("h_t_REC_wCUT_pi3");
	TH1D* h_t_REC_wCUT_pi4 = (TH1D*) file->Get("h_t_REC_wCUT_pi4");
	TH1D* h_t_REC_wCUT_pi6 = (TH1D*) file->Get("h_t_REC_wCUT_pi6");
	TH1D* h_t_REC_wCUT_pi9 = (TH1D*) file->Get("h_t_REC_wCUT_pi9");
	TH1D* h_t_REC_wCUT_pi12 = (TH1D*) file->Get("h_t_REC_wCUT_pi12");
	TH2D* h_t_REC_2d = (TH2D*) file->Get("h_t_REC_2d");
	TH2D* h_t_REC_2d_wRES_cut_pi2 = (TH2D*) file->Get("h_t_REC_2d_wRES_cut_pi2");
	TH2D* h_t_REC_2d_wRES_cut_pi3 = (TH2D*) file->Get("h_t_REC_2d_wRES_cut_pi3");
	TH2D* h_t_REC_2d_wRES_cut_pi4 = (TH2D*) file->Get("h_t_REC_2d_wRES_cut_pi4");
	TH2D* h_t_REC_2d_wRES_cut_pi6 = (TH2D*) file->Get("h_t_REC_2d_wRES_cut_pi6");
	TH2D* h_t_REC_2d_wRES_cut_pi9 = (TH2D*) file->Get("h_t_REC_2d_wRES_cut_pi9");
	TH2D* h_t_REC_2d_wRES_cut_pi12 = (TH2D*) file->Get("h_t_REC_2d_wRES_cut_pi12");
	TH2D* h_t_REC_2d_wRES = (TH2D*) file->Get("h_t_REC_2d_wRES");
	TH2D* h_t_REC_2d_wCUT_pi2 = (TH2D*) file->Get("h_t_REC_2d_wCUT_pi2");
	TH2D* h_t_REC_2d_wCUT_pi3 = (TH2D*) file->Get("h_t_REC_2d_wCUT_pi3");
	TH2D* h_t_REC_2d_wCUT_pi4 = (TH2D*) file->Get("h_t_REC_2d_wCUT_pi4");
	TH2D* h_t_REC_2d_wCUT_pi6 = (TH2D*) file->Get("h_t_REC_2d_wCUT_pi6");
	TH2D* h_t_REC_2d_wCUT_pi9 = (TH2D*) file->Get("h_t_REC_2d_wCUT_pi9");
	TH2D* h_t_REC_2d_wCUT_pi12 = (TH2D*) file->Get("h_t_REC_2d_wCUT_pi12");

	//create canvas for each angle
	const int numCanvas = 6; 

	//angle labels
	TString angleLabels[numCanvas] = 
	{
    	angle_pi2, angle_pi3, angle_pi4, angle_pi6, angle_pi9, angle_pi12
	};
	double angles[numCanvas] = 
	{
		M_PI/2, M_PI/3, M_PI/4, M_PI/6, M_PI/9, M_PI/12
	};
	
	//histograms
	TH1D* h_t_REC_wRES_cut[numCanvas] = 
	{
		h_t_REC_wRES_cut_pi2,
		h_t_REC_wRES_cut_pi3,
		h_t_REC_wRES_cut_pi4,
		h_t_REC_wRES_cut_pi6,
		h_t_REC_wRES_cut_pi9,
		h_t_REC_wRES_cut_pi12
	};
	TH1D* h_t_REC_wCUT[numCanvas] = 
	{
		h_t_REC_wCUT_pi2,
		h_t_REC_wCUT_pi3,
		h_t_REC_wCUT_pi4,
		h_t_REC_wCUT_pi6,
		h_t_REC_wCUT_pi9,
		h_t_REC_wCUT_pi12
	};
	TH2D* h_t_REC_2d_wRES_cut[numCanvas] = 
	{
		h_t_REC_2d_wRES_cut_pi2,
		h_t_REC_2d_wRES_cut_pi3,
		h_t_REC_2d_wRES_cut_pi4,
		h_t_REC_2d_wRES_cut_pi6,
		h_t_REC_2d_wRES_cut_pi9,
		h_t_REC_2d_wRES_cut_pi12
	};
	TH2D* h_t_REC_2d_wCUT[numCanvas] = 
	{
		h_t_REC_2d_wCUT_pi2,
		h_t_REC_2d_wCUT_pi3,
		h_t_REC_2d_wCUT_pi4,
		h_t_REC_2d_wCUT_pi6,
		h_t_REC_2d_wCUT_pi9,
		h_t_REC_2d_wCUT_pi12
	};
	
	//reproduce truth plot with our decomposition
	TCanvas* newMethod_canvases[numCanvas];
	for (int i=0; i<1; i++) 
	{
	TString canvasName = Form("newMethod_c%d", i+1);
	newMethod_canvases[i] = new TCanvas(canvasName, canvasName, 1, 1, 600, 600);
	gPad->SetLogy(1);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	gPad->SetRightMargin(0.01);
	// Generate histogram for each canvas
	TString newMethod_histName = Form("newMethod_base%d", i+1);
	TH1D* newMethod_baseHist = makeHist(newMethod_histName, "", "|#it{t} | (GeV^{2})", "dN/d|#it{t} | (GeV^{-2}) ", 100, 0, 0.18, kBlack);
	newMethod_baseHist->GetYaxis()->SetRangeUser(8e-2, 8e8);
	newMethod_baseHist->GetXaxis()->SetTitleColor(kBlack);
	fixedFontHist1D(newMethod_baseHist, 1., 1.2);
	newMethod_baseHist->GetYaxis()->SetTitleSize(newMethod_baseHist->GetYaxis()->GetTitleSize()*1.5);
	newMethod_baseHist->GetXaxis()->SetTitleSize(newMethod_baseHist->GetXaxis()->GetTitleSize()*1.5);
	newMethod_baseHist->GetYaxis()->SetLabelSize(newMethod_baseHist->GetYaxis()->GetLabelSize()*1.5);
	newMethod_baseHist->GetXaxis()->SetLabelSize(newMethod_baseHist->GetXaxis()->GetLabelSize()*1.5);
	newMethod_baseHist->GetXaxis()->SetNdivisions(4,4,0);
	newMethod_baseHist->GetYaxis()->SetNdivisions(5,5,0);
	newMethod_baseHist->Draw();
	// get counts
	double integral_MC = h_t_MC->Integral();
	double integral_REC = h_t_REC->Integral();
	double integral_REC_new_method = h_t_REC_new_method->Integral();
	// Draw histograms 
	h_t_MC->Draw("same");

	h_t_REC->SetMarkerStyle(20);
	h_t_REC->Draw("PEsame");

	h_t_REC_new_method->SetMarkerStyle(30);
	h_t_REC_new_method->SetMarkerColor(kRed);
	h_t_REC_new_method->Draw("P same");

	// Add labels
	TLatex* r42 = new TLatex(0.18, 0.91, "eAu 18x110 GeV");
	r42->SetNDC();
	r42->SetTextSize(22);
	r42->SetTextFont(43);
	r42->SetTextColor(kBlack);
	r42->Draw("same");

	TLatex* r43 = new TLatex(0.9, 0.91, "EPIC");
	r43->SetNDC();
	r43->SetTextSize(0.04);
	r43->Draw("same");

	TLatex* r44 = new TLatex(0.18, 0.84, "1<Q^{2}<10 GeV^{2}, 0.01 < y < 0.85");
	r44->SetNDC();
	r44->SetTextSize(20);
	r44->SetTextFont(43);
	r44->SetTextColor(kBlack);
	r44->Draw("same");

	TLatex* r44_0 = new TLatex(0.18, 0.79, Form("|y_{%s}|<3.5, |M_{inv} #minus M_{%s}| < 0.02 GeV", vm_label.Data(), vm_label.Data()));
	r44_0->SetNDC();
	r44_0->SetTextSize(20);
	r44_0->SetTextFont(43);
	r44_0->SetTextColor(kBlack);
	r44_0->Draw("same");

	TLatex* r44_2 = new TLatex(0.18, 0.18, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	r44_2->SetNDC();
	r44_2->SetTextSize(30);
	r44_2->SetTextFont(43);
	r44_2->SetTextColor(kBlack);
	r44_2->Draw("same");
	// Add legend
	TLegend* w7 = new TLegend(0.18, 0.68, 0.93, 0.76);
	w7->SetLineColor(kWhite);
	w7->SetFillColor(0);
	w7->SetTextSize(17);
	w7->SetTextFont(45);
	w7->AddEntry(h_t_MC, Form("Sartre %s MC: %.f events", vm_label.Data(),integral_MC), "L");
	w7->AddEntry(h_t_REC, Form("Sartre %s RECO w. EEMC: %.f events", vm_label.Data(),integral_REC), "P");
	w7->AddEntry(h_t_REC_new_method, Form("Sartre %s RECO' new method: %.f events", vm_label.Data(),integral_REC_new_method), "P");
	w7->Draw("same");
	// Save figure
	newMethod_canvases[i]->Print("./figures/new_method.pdf");
	}
	
	//plots with resolution and cut, normalized to scale with method L reco
	TCanvas* totalNorm_canvases_resCut[numCanvas];
	for (int i=0; i<numCanvas; i++) 
	{
    	TString canvasName = Form("totalNorm_resCut_c%d", i+1);
    	totalNorm_canvases_resCut[i] = new TCanvas(canvasName, canvasName, 1, 1, 600, 600);
    	gPad->SetLogy(1);
    	gPad->SetTicks();
    	gPad->SetLeftMargin(0.15);
    	gPad->SetBottomMargin(0.15);
    	gPad->SetRightMargin(0.01);
		// Generate histogram for each canvas
    	TString totalNorm_histName_resCut = Form("totalNorm_resCut_base%d", i+1);
    	TH1D* totalNorm_baseHist_resCut = makeHist(totalNorm_histName_resCut, "", "|#it{t} | (GeV^{2})", "dN/d|#it{t} | (GeV^{-2}) ", 100, 0, 0.18, kBlack);
    	totalNorm_baseHist_resCut->GetYaxis()->SetRangeUser(8e-2, 8e8);
    	totalNorm_baseHist_resCut->GetXaxis()->SetTitleColor(kBlack);
    	fixedFontHist1D(totalNorm_baseHist_resCut, 1., 1.2);
    	totalNorm_baseHist_resCut->GetYaxis()->SetTitleSize(totalNorm_baseHist_resCut->GetYaxis()->GetTitleSize()*1.5);
    	totalNorm_baseHist_resCut->GetXaxis()->SetTitleSize(totalNorm_baseHist_resCut->GetXaxis()->GetTitleSize()*1.5);
    	totalNorm_baseHist_resCut->GetYaxis()->SetLabelSize(totalNorm_baseHist_resCut->GetYaxis()->GetLabelSize()*1.5);
    	totalNorm_baseHist_resCut->GetXaxis()->SetLabelSize(totalNorm_baseHist_resCut->GetXaxis()->GetLabelSize()*1.5);
    	totalNorm_baseHist_resCut->GetXaxis()->SetNdivisions(4,4,0);
    	totalNorm_baseHist_resCut->GetYaxis()->SetNdivisions(5,5,0);
    	totalNorm_baseHist_resCut->Draw();
		// Normalize histograms
		double integral_MC = h_t_MC->Integral();
		double integral_REC = h_t_REC->Integral();
		double integral_REC_wRES_cut[numCanvas];
		integral_REC_wRES_cut[i] = h_t_REC_wRES_cut[i]->Integral();
		h_t_REC_wRES_cut[i]->Scale(M_PI/2/angles[i]);
		// Draw histograms with different angle cuts
		h_t_MC->Draw("same");

		h_t_REC->SetMarkerStyle(20);
		h_t_REC->Draw("PEsame");

		h_t_REC_wRES_cut[i]->SetMarkerStyle(30);
		h_t_REC_wRES_cut[i]->SetMarkerColor(kRed);
		h_t_REC_wRES_cut[i]->Draw("P same");

		// Add labels
		TLatex* r42 = new TLatex(0.18, 0.91, "eAu 18x110 GeV");
		r42->SetNDC();
		r42->SetTextSize(22);
		r42->SetTextFont(43);
		r42->SetTextColor(kBlack);
		r42->Draw("same");
	
		TLatex* r43 = new TLatex(0.9, 0.91, "EPIC");
		r43->SetNDC();
		r43->SetTextSize(0.04);
		r43->Draw("same");
	
		TLatex* r44 = new TLatex(0.18, 0.84, "1<Q^{2}<10 GeV^{2}, 0.01 < y < 0.85");
		r44->SetNDC();
		r44->SetTextSize(20);
		r44->SetTextFont(43);
		r44->SetTextColor(kBlack);
		r44->Draw("same");
	
		TLatex* r44_0 = new TLatex(0.18, 0.79, Form("|y_{%s}|<3.5, |M_{inv} #minus M_{%s}| < 0.02 GeV", vm_label.Data(), vm_label.Data()));
		r44_0->SetNDC();
		r44_0->SetTextSize(20);
		r44_0->SetTextFont(43);
		r44_0->SetTextColor(kBlack);
		r44_0->Draw("same");
	
		TLatex* r44_2 = new TLatex(0.18, 0.18, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
		r44_2->SetNDC();
		r44_2->SetTextSize(30);
		r44_2->SetTextFont(43);
		r44_2->SetTextColor(kBlack);
		r44_2->Draw("same");

		// Add legend
		TLegend* w7 = new TLegend(0.18, 0.68, 0.93, 0.76);
		w7->SetLineColor(kWhite);
		w7->SetFillColor(0);
		w7->SetTextSize(17);
		w7->SetTextFont(45);
		w7->AddEntry(h_t_MC, Form("Sartre %s MC: %.f events", vm_label.Data(),integral_MC), "L");
		w7->AddEntry(h_t_REC, Form("Sartre %s RECO w. EEMC: %.f events", vm_label.Data(),integral_REC), "P");
		w7->AddEntry(h_t_REC_wRES_cut[i], Form("Sartre %s RECO' #theta_{Max}= %s: %.f events", vm_label.Data(), angleLabels[i].Data(),integral_REC_wRES_cut[i]), "P");
		w7->Draw("same");

		TLegend *w8 = new TLegend(0.48,0.64,0.93,0.56);
		w8->SetLineColor(kWhite);
		w8->SetFillColor(0);
		w8->SetTextSize(17);
		w8->SetTextFont(45);
		w8->AddEntry(h_t_REC_wRES_cut[i],"normalization: #frac{#pi}{2} / #theta_{max}","P");
		w8->Draw("same");
	
		// Save figure
		TString cleanLabel = angleLabels[i];
		cleanLabel.ReplaceAll("#", "");
		cleanLabel.ReplaceAll("/", "");
		totalNorm_canvases_resCut[i]->Print(Form("./figures/normL_wRES_cut_angle%s.pdf",cleanLabel.Data()));
	}
	
	//plots with resolution and cut, normalized to scale with true t reco
	TCanvas* Lnorm_canvases_resCut[numCanvas];
	for (int i=0; i<numCanvas; i++) 
	{
    	TString canvasName = Form("Lnorm_resCut_c%d", i+1);
    	Lnorm_canvases_resCut[i] = new TCanvas(canvasName, canvasName, 1, 1, 600, 600);
    	gPad->SetLogy(1);
    	gPad->SetTicks();
    	gPad->SetLeftMargin(0.15);
    	gPad->SetBottomMargin(0.15);
    	gPad->SetRightMargin(0.01);
		// Generate histogram for each canvas
    	TString Lnorm_histName_resCut = Form("Lnorm_resCut_base%d", i+1);
    	TH1D* Lnorm_baseHist_resCut = makeHist(Lnorm_histName_resCut, "", "|#it{t} | (GeV^{2})", "dN/d|#it{t} | (GeV^{-2}) ", 100, 0, 0.18, kBlack);
    	Lnorm_baseHist_resCut->GetYaxis()->SetRangeUser(8e-2, 8e8);
    	Lnorm_baseHist_resCut->GetXaxis()->SetTitleColor(kBlack);
    	fixedFontHist1D(Lnorm_baseHist_resCut, 1., 1.2);
    	Lnorm_baseHist_resCut->GetYaxis()->SetTitleSize(Lnorm_baseHist_resCut->GetYaxis()->GetTitleSize()*1.5);
    	Lnorm_baseHist_resCut->GetXaxis()->SetTitleSize(Lnorm_baseHist_resCut->GetXaxis()->GetTitleSize()*1.5);
    	Lnorm_baseHist_resCut->GetYaxis()->SetLabelSize(Lnorm_baseHist_resCut->GetYaxis()->GetLabelSize()*1.5);
    	Lnorm_baseHist_resCut->GetXaxis()->SetLabelSize(Lnorm_baseHist_resCut->GetXaxis()->GetLabelSize()*1.5);
    	Lnorm_baseHist_resCut->GetXaxis()->SetNdivisions(4,4,0);
    	Lnorm_baseHist_resCut->GetYaxis()->SetNdivisions(5,5,0);
    	Lnorm_baseHist_resCut->Draw();
		// Normalize histograms
		double integral_MC = h_t_MC->Integral();
		double integral_REC = h_t_REC->Integral();
		double integral_REC_wRES_cut[numCanvas];
		integral_REC_wRES_cut[i] = h_t_REC_wRES_cut[i]->Integral();
		if(integral_MC>0 && integral_REC_wRES_cut[i]>0) 
		{
    		h_t_REC_wRES_cut[i]->Scale(integral_MC/integral_REC_wRES_cut[i]);
		}
		// Draw histograms with different angle cuts
		h_t_MC->Draw("same");

		h_t_REC->SetMarkerStyle(20);
		h_t_REC->Draw("PEsame");

		h_t_REC_wRES_cut[i]->SetMarkerStyle(30);
		h_t_REC_wRES_cut[i]->SetMarkerColor(kRed);
		h_t_REC_wRES_cut[i]->Draw("P same");

		// Add labels
		TLatex* r42 = new TLatex(0.18, 0.91, "eAu 18x110 GeV");
		r42->SetNDC();
		r42->SetTextSize(22);
		r42->SetTextFont(43);
		r42->SetTextColor(kBlack);
		r42->Draw("same");
	
		TLatex* r43 = new TLatex(0.9, 0.91, "EPIC");
		r43->SetNDC();
		r43->SetTextSize(0.04);
		r43->Draw("same");
	
		TLatex* r44 = new TLatex(0.18, 0.84, "1<Q^{2}<10 GeV^{2}, 0.01 < y < 0.85");
		r44->SetNDC();
		r44->SetTextSize(20);
		r44->SetTextFont(43);
		r44->SetTextColor(kBlack);
		r44->Draw("same");
	
		TLatex* r44_0 = new TLatex(0.18, 0.79, Form("|y_{%s}|<3.5, |M_{inv} #minus M_{%s}| < 0.02 GeV", vm_label.Data(), vm_label.Data()));
		r44_0->SetNDC();
		r44_0->SetTextSize(20);
		r44_0->SetTextFont(43);
		r44_0->SetTextColor(kBlack);
		r44_0->Draw("same");
	
		TLatex* r44_2 = new TLatex(0.18, 0.18, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
		r44_2->SetNDC();
		r44_2->SetTextSize(30);
		r44_2->SetTextFont(43);
		r44_2->SetTextColor(kBlack);
		r44_2->Draw("same");

		// Add legend
		TLegend* w7 = new TLegend(0.18, 0.68, 0.93, 0.76);
		w7->SetLineColor(kWhite);
		w7->SetFillColor(0);
		w7->SetTextSize(17);
		w7->SetTextFont(45);
		w7->AddEntry(h_t_MC, Form("Sartre %s MC: %.f events", vm_label.Data(),integral_MC), "L");
		w7->AddEntry(h_t_REC, Form("Sartre %s RECO w. EEMC: %.f events", vm_label.Data(),integral_REC), "P");
		w7->AddEntry(h_t_REC_wRES_cut[i], Form("Sartre %s RECO' #theta_{Max}= %s: %.f events", vm_label.Data(), angleLabels[i].Data(),integral_REC_wRES_cut[i]), "P");
		w7->Draw("same");

		TLegend *w8 = new TLegend(0.48,0.64,0.93,0.56);
		w8->SetLineColor(kWhite);
		w8->SetFillColor(0);
		w8->SetTextSize(17);
		w8->SetTextFont(45);
		w8->AddEntry(h_t_REC_wRES_cut[i],"normalization: #int|#it{t}|_{MC}/#int|#it{t}|_{RECO'}","P");
		w8->Draw("same");
	
		// Save figure
		TString cleanLabel = angleLabels[i];
		cleanLabel.ReplaceAll("#", "");
		cleanLabel.ReplaceAll("/", "");
		Lnorm_canvases_resCut[i]->Print(Form("./figures/normTotal_wRES_cut_angle%s.pdf",cleanLabel.Data()));
	}
	
	//plot with resolution, normalized to scale with method L reco
	TCanvas* totalNorm_canvases_res[numCanvas];
	for (int i=0; i<1; i++) 
	{
    	TString canvasName = Form("totalNorm_res_c%d", i+1);
    	totalNorm_canvases_res[i] = new TCanvas(canvasName, canvasName, 1, 1, 600, 600);
    	gPad->SetLogy(1);
    	gPad->SetTicks();
    	gPad->SetLeftMargin(0.15);
    	gPad->SetBottomMargin(0.15);
    	gPad->SetRightMargin(0.01);
		// Generate histogram for each canvas
    	TString totalNorm_histName_res = Form("totalNorm_res_base%d", i+1);
    	TH1D* totalNorm_baseHist_res = makeHist(totalNorm_histName_res, "", "|#it{t} | (GeV^{2})", "dN/d|#it{t} | (GeV^{-2}) ", 100, 0, 0.18, kBlack);
    	totalNorm_baseHist_res->GetYaxis()->SetRangeUser(8e-2, 8e8);
    	totalNorm_baseHist_res->GetXaxis()->SetTitleColor(kBlack);
    	fixedFontHist1D(totalNorm_baseHist_res, 1., 1.2);
    	totalNorm_baseHist_res->GetYaxis()->SetTitleSize(totalNorm_baseHist_res->GetYaxis()->GetTitleSize()*1.5);
    	totalNorm_baseHist_res->GetXaxis()->SetTitleSize(totalNorm_baseHist_res->GetXaxis()->GetTitleSize()*1.5);
    	totalNorm_baseHist_res->GetYaxis()->SetLabelSize(totalNorm_baseHist_res->GetYaxis()->GetLabelSize()*1.5);
    	totalNorm_baseHist_res->GetXaxis()->SetLabelSize(totalNorm_baseHist_res->GetXaxis()->GetLabelSize()*1.5);
    	totalNorm_baseHist_res->GetXaxis()->SetNdivisions(4,4,0);
    	totalNorm_baseHist_res->GetYaxis()->SetNdivisions(5,5,0);
    	totalNorm_baseHist_res->Draw();
		// Normalize histograms
		double integral_MC = h_t_MC->Integral();
		double integral_REC = h_t_REC->Integral();
		double integral_REC_wRES = h_t_REC_wRES->Integral();
		// Draw histograms
		h_t_MC->Draw("same");

		h_t_REC->SetMarkerStyle(20);
		h_t_REC->Draw("PEsame");

		h_t_REC_wRES->SetMarkerStyle(30);
		h_t_REC_wRES->SetMarkerColor(kRed);
		h_t_REC_wRES->Draw("P same");

		// Add labels
		TLatex* r42 = new TLatex(0.18, 0.91, "eAu 18x110 GeV");
		r42->SetNDC();
		r42->SetTextSize(22);
		r42->SetTextFont(43);
		r42->SetTextColor(kBlack);
		r42->Draw("same");
	
		TLatex* r43 = new TLatex(0.9, 0.91, "EPIC");
		r43->SetNDC();
		r43->SetTextSize(0.04);
		r43->Draw("same");
	
		TLatex* r44 = new TLatex(0.18, 0.84, "1<Q^{2}<10 GeV^{2}, 0.01 < y < 0.85");
		r44->SetNDC();
		r44->SetTextSize(20);
		r44->SetTextFont(43);
		r44->SetTextColor(kBlack);
		r44->Draw("same");
	
		TLatex* r44_0 = new TLatex(0.18, 0.79, Form("|y_{%s}|<3.5, |M_{inv} #minus M_{%s}| < 0.02 GeV", vm_label.Data(), vm_label.Data()));
		r44_0->SetNDC();
		r44_0->SetTextSize(20);
		r44_0->SetTextFont(43);
		r44_0->SetTextColor(kBlack);
		r44_0->Draw("same");
	
		TLatex* r44_2 = new TLatex(0.18, 0.18, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
		r44_2->SetNDC();
		r44_2->SetTextSize(30);
		r44_2->SetTextFont(43);
		r44_2->SetTextColor(kBlack);
		r44_2->Draw("same");

		// Add legend
		TLegend* w7 = new TLegend(0.18, 0.68, 0.93, 0.76);
		w7->SetLineColor(kWhite);
		w7->SetFillColor(0);
		w7->SetTextSize(17);
		w7->SetTextFont(45);
		w7->AddEntry(h_t_MC, Form("Sartre %s MC: %.f events", vm_label.Data(),integral_MC), "L");
		w7->AddEntry(h_t_REC, Form("Sartre %s RECO w. EEMC: %.f events", vm_label.Data(),integral_REC), "P");
		w7->AddEntry(h_t_REC_wRES, Form("Sartre %s RECO' new method: %.f events", vm_label.Data(),integral_REC_wRES), "P");
		w7->Draw("same");

		TLegend *w8 = new TLegend(0.48,0.64,0.93,0.56);
		w8->SetLineColor(kWhite);
		w8->SetFillColor(0);
		w8->SetTextSize(17);
		w8->SetTextFont(45);
		w8->AddEntry(h_t_REC_wRES,"normalization: #frac{#pi}{2} / #theta_{max}","P");
		w8->Draw("same");
		//save figure
		totalNorm_canvases_res[i]->Print("./figures/normL_wRES.pdf");
	}
	
	//plot with resolution, normalized to scale with true t reco
	TCanvas* Lnorm_canvases_res[numCanvas];
	for (int i=0; i<1; i++) 
	{
    	TString canvasName = Form("Lnorm_res_c%d", i+1);
    	Lnorm_canvases_res[i] = new TCanvas(canvasName, canvasName, 1, 1, 600, 600);
    	gPad->SetLogy(1);
    	gPad->SetTicks();
    	gPad->SetLeftMargin(0.15);
    	gPad->SetBottomMargin(0.15);
    	gPad->SetRightMargin(0.01);
		// Generate histogram for each canvas
    	TString Lnorm_histName_res = Form("Lnorm_res_base%d", i+1);
    	TH1D* Lnorm_baseHist_res = makeHist(Lnorm_histName_res, "", "|#it{t} | (GeV^{2})", "dN/d|#it{t} | (GeV^{-2}) ", 100, 0, 0.18, kBlack);
    	Lnorm_baseHist_res->GetYaxis()->SetRangeUser(8e-2, 8e8);
    	Lnorm_baseHist_res->GetXaxis()->SetTitleColor(kBlack);
    	fixedFontHist1D(Lnorm_baseHist_res, 1., 1.2);
    	Lnorm_baseHist_res->GetYaxis()->SetTitleSize(Lnorm_baseHist_res->GetYaxis()->GetTitleSize()*1.5);
    	Lnorm_baseHist_res->GetXaxis()->SetTitleSize(Lnorm_baseHist_res->GetXaxis()->GetTitleSize()*1.5);
    	Lnorm_baseHist_res->GetYaxis()->SetLabelSize(Lnorm_baseHist_res->GetYaxis()->GetLabelSize()*1.5);
    	Lnorm_baseHist_res->GetXaxis()->SetLabelSize(Lnorm_baseHist_res->GetXaxis()->GetLabelSize()*1.5);
    	Lnorm_baseHist_res->GetXaxis()->SetNdivisions(4,4,0);
    	Lnorm_baseHist_res->GetYaxis()->SetNdivisions(5,5,0);
    	Lnorm_baseHist_res->Draw();
		// Normalize histograms
		double integral_MC = h_t_MC->Integral();
		double integral_REC = h_t_REC->Integral();
		double integral_REC_wRES = h_t_REC_wRES->Integral();
		// Draw histograms 
		h_t_MC->Draw("same");

		h_t_REC->SetMarkerStyle(20);
		h_t_REC->Draw("PEsame");

		h_t_REC_wRES->SetMarkerStyle(30);
		h_t_REC_wRES->SetMarkerColor(kRed);
		h_t_REC_wRES->Draw("P same");

		// Add labels
		TLatex* r42 = new TLatex(0.18, 0.91, "eAu 18x110 GeV");
		r42->SetNDC();
		r42->SetTextSize(22);
		r42->SetTextFont(43);
		r42->SetTextColor(kBlack);
		r42->Draw("same");
	
		TLatex* r43 = new TLatex(0.9, 0.91, "EPIC");
		r43->SetNDC();
		r43->SetTextSize(0.04);
		r43->Draw("same");
	
		TLatex* r44 = new TLatex(0.18, 0.84, "1<Q^{2}<10 GeV^{2}, 0.01 < y < 0.85");
		r44->SetNDC();
		r44->SetTextSize(20);
		r44->SetTextFont(43);
		r44->SetTextColor(kBlack);
		r44->Draw("same");
	
		TLatex* r44_0 = new TLatex(0.18, 0.79, Form("|y_{%s}|<3.5, |M_{inv} #minus M_{%s}| < 0.02 GeV", vm_label.Data(), vm_label.Data()));
		r44_0->SetNDC();
		r44_0->SetTextSize(20);
		r44_0->SetTextFont(43);
		r44_0->SetTextColor(kBlack);
		r44_0->Draw("same");
	
		TLatex* r44_2 = new TLatex(0.18, 0.18, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
		r44_2->SetNDC();
		r44_2->SetTextSize(30);
		r44_2->SetTextFont(43);
		r44_2->SetTextColor(kBlack);
		r44_2->Draw("same");

		// Add legend
		TLegend* w7 = new TLegend(0.18, 0.68, 0.93, 0.76);
		w7->SetLineColor(kWhite);
		w7->SetFillColor(0);
		w7->SetTextSize(17);
		w7->SetTextFont(45);
		w7->AddEntry(h_t_MC, Form("Sartre %s MC: %.f events", vm_label.Data(),integral_MC), "L");
		w7->AddEntry(h_t_REC, Form("Sartre %s RECO w. EEMC: %.f events", vm_label.Data(),integral_REC), "P");
		w7->AddEntry(h_t_REC_wRES, Form("Sartre %s RECO' new method: %.f events", vm_label.Data(),integral_REC_wRES), "P");
		w7->Draw("same");

		TLegend *w8 = new TLegend(0.48,0.64,0.93,0.56);
		w8->SetLineColor(kWhite);
		w8->SetFillColor(0);
		w8->SetTextSize(17);
		w8->SetTextFont(45);
		w8->AddEntry(h_t_REC_wRES,"normalization: #int|#it{t}|_{MC}/#int|#it{t}|_{RECO'}","P");
		w8->Draw("same");
	
		// Save figure
		Lnorm_canvases_res[i]->Print("./figures/normTotal_wRES.pdf");
	}
	
	//plots with cut, normalized to scale with method L reco
	TCanvas* totalNorm_canvases_Cut[numCanvas];
	for (int i=0; i<numCanvas; i++) 
	{
	TString canvasName = Form("totalNorm_Cut_c%d", i+1);
	totalNorm_canvases_Cut[i] = new TCanvas(canvasName, canvasName, 1, 1, 600, 600);
	gPad->SetLogy(1);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	gPad->SetRightMargin(0.01);
	// Generate histogram for each canvas
	TString totalNorm_histName_Cut = Form("totalNorm_Cut_base%d", i+1);
	TH1D* totalNorm_baseHist_Cut = makeHist(totalNorm_histName_Cut, "", "|#it{t} | (GeV^{2})", "dN/d|#it{t} | (GeV^{-2}) ", 100, 0, 0.18, kBlack);
	totalNorm_baseHist_Cut->GetYaxis()->SetRangeUser(8e-2, 8e8);
	totalNorm_baseHist_Cut->GetXaxis()->SetTitleColor(kBlack);
	fixedFontHist1D(totalNorm_baseHist_Cut, 1., 1.2);
	totalNorm_baseHist_Cut->GetYaxis()->SetTitleSize(totalNorm_baseHist_Cut->GetYaxis()->GetTitleSize()*1.5);
	totalNorm_baseHist_Cut->GetXaxis()->SetTitleSize(totalNorm_baseHist_Cut->GetXaxis()->GetTitleSize()*1.5);
	totalNorm_baseHist_Cut->GetYaxis()->SetLabelSize(totalNorm_baseHist_Cut->GetYaxis()->GetLabelSize()*1.5);
	totalNorm_baseHist_Cut->GetXaxis()->SetLabelSize(totalNorm_baseHist_Cut->GetXaxis()->GetLabelSize()*1.5);
	totalNorm_baseHist_Cut->GetXaxis()->SetNdivisions(4,4,0);
	totalNorm_baseHist_Cut->GetYaxis()->SetNdivisions(5,5,0);
	totalNorm_baseHist_Cut->Draw();
	// Normalize histograms
	double integral_MC = h_t_MC->Integral();
	double integral_REC = h_t_REC->Integral();
	double integral_REC_wCut[numCanvas];
	integral_REC_wCut[i] = h_t_REC_wCUT[i]->Integral();
	h_t_REC_wCUT[i]->Scale(M_PI/2/angles[i]);
	// Draw histograms with different angle cuts
	h_t_MC->Draw("same");

	h_t_REC->SetMarkerStyle(20);
	h_t_REC->Draw("PEsame");

	h_t_REC_wCUT[i]->SetMarkerStyle(30);
	h_t_REC_wCUT[i]->SetMarkerColor(kRed);
	h_t_REC_wCUT[i]->Draw("P same");

	// Add labels
	TLatex* r42 = new TLatex(0.18, 0.91, "eAu 18x110 GeV");
	r42->SetNDC();
	r42->SetTextSize(22);
	r42->SetTextFont(43);
	r42->SetTextColor(kBlack);
	r42->Draw("same");

	TLatex* r43 = new TLatex(0.9, 0.91, "EPIC");
	r43->SetNDC();
	r43->SetTextSize(0.04);
	r43->Draw("same");

	TLatex* r44 = new TLatex(0.18, 0.84, "1<Q^{2}<10 GeV^{2}, 0.01 < y < 0.85");
	r44->SetNDC();
	r44->SetTextSize(20);
	r44->SetTextFont(43);
	r44->SetTextColor(kBlack);
	r44->Draw("same");

	TLatex* r44_0 = new TLatex(0.18, 0.79, Form("|y_{%s}|<3.5, |M_{inv} #minus M_{%s}| < 0.02 GeV", vm_label.Data(), vm_label.Data()));
	r44_0->SetNDC();
	r44_0->SetTextSize(20);
	r44_0->SetTextFont(43);
	r44_0->SetTextColor(kBlack);
	r44_0->Draw("same");

	TLatex* r44_2 = new TLatex(0.18, 0.18, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	r44_2->SetNDC();
	r44_2->SetTextSize(30);
	r44_2->SetTextFont(43);
	r44_2->SetTextColor(kBlack);
	r44_2->Draw("same");

	// Add legend
	TLegend* w7 = new TLegend(0.18, 0.68, 0.93, 0.76);
	w7->SetLineColor(kWhite);
	w7->SetFillColor(0);
	w7->SetTextSize(17);
	w7->SetTextFont(45);
	w7->AddEntry(h_t_MC, Form("Sartre %s MC: %.f events", vm_label.Data(),integral_MC), "L");
	w7->AddEntry(h_t_REC, Form("Sartre %s RECO w. EEMC: %.f events", vm_label.Data(),integral_REC), "P");
	w7->AddEntry(h_t_REC_wCUT[i], Form("Sartre %s RECO' #theta_{Max}= %s: %.f events", vm_label.Data(), angleLabels[i].Data(),integral_REC_wCut[i]), "P");
	w7->Draw("same");

	TLegend *w8 = new TLegend(0.48,0.64,0.93,0.56);
	w8->SetLineColor(kWhite);
	w8->SetFillColor(0);
	w8->SetTextSize(17);
	w8->SetTextFont(45);
	w8->AddEntry(h_t_REC_wCUT[i],"normalization: #frac{#pi}{2} / #theta_{max}","P");
	w8->Draw("same");

	// Save figure
	TString cleanLabel = angleLabels[i];
	cleanLabel.ReplaceAll("#", "");
	cleanLabel.ReplaceAll("/", "");
	totalNorm_canvases_Cut[i]->Print(Form("./figures/normL_wCUT_angle%s.pdf",cleanLabel.Data()));
}
	
	//plots with cut, normalized to scale with true t reco
	TCanvas* Lnorm_canvases_Cut[numCanvas];
	for (int i=0; i<numCanvas; i++) 
{
	TString canvasName = Form("Lnorm_Cut_c%d", i+1);
	Lnorm_canvases_Cut[i] = new TCanvas(canvasName, canvasName, 1, 1, 600, 600);
	gPad->SetLogy(1);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	gPad->SetRightMargin(0.01);
	// Generate histogram for each canvas
	TString Lnorm_histName_Cut = Form("Lnorm_Cut_base%d", i+1);
	TH1D* Lnorm_baseHist_Cut = makeHist(Lnorm_histName_Cut, "", "|#it{t} | (GeV^{2})", "dN/d|#it{t} | (GeV^{-2}) ", 100, 0, 0.18, kBlack);
	Lnorm_baseHist_Cut->GetYaxis()->SetRangeUser(8e-2, 8e8);
	Lnorm_baseHist_Cut->GetXaxis()->SetTitleColor(kBlack);
	fixedFontHist1D(Lnorm_baseHist_Cut, 1., 1.2);
	Lnorm_baseHist_Cut->GetYaxis()->SetTitleSize(Lnorm_baseHist_Cut->GetYaxis()->GetTitleSize()*1.5);
	Lnorm_baseHist_Cut->GetXaxis()->SetTitleSize(Lnorm_baseHist_Cut->GetXaxis()->GetTitleSize()*1.5);
	Lnorm_baseHist_Cut->GetYaxis()->SetLabelSize(Lnorm_baseHist_Cut->GetYaxis()->GetLabelSize()*1.5);
	Lnorm_baseHist_Cut->GetXaxis()->SetLabelSize(Lnorm_baseHist_Cut->GetXaxis()->GetLabelSize()*1.5);
	Lnorm_baseHist_Cut->GetXaxis()->SetNdivisions(4,4,0);
	Lnorm_baseHist_Cut->GetYaxis()->SetNdivisions(5,5,0);
	Lnorm_baseHist_Cut->Draw();
	// Normalize histograms
	double integral_MC = h_t_MC->Integral();
	double integral_REC = h_t_REC->Integral();
	double integral_REC_wCut[numCanvas];
	integral_REC_wCut[i] = h_t_REC_wCUT[i]->Integral();
	if(integral_MC>0 && integral_REC_wCut[i]>0) 
	{
		h_t_REC_wCUT[i]->Scale(integral_MC/integral_REC_wCut[i]);
	}
	// Draw histograms with different angle cuts
	h_t_MC->Draw("same");

	h_t_REC->SetMarkerStyle(20);
	h_t_REC->Draw("PEsame");

	h_t_REC_wCUT[i]->SetMarkerStyle(30);
	h_t_REC_wCUT[i]->SetMarkerColor(kRed);
	h_t_REC_wCUT[i]->Draw("P same");

	// Add labels
	TLatex* r42 = new TLatex(0.18, 0.91, "eAu 18x110 GeV");
	r42->SetNDC();
	r42->SetTextSize(22);
	r42->SetTextFont(43);
	r42->SetTextColor(kBlack);
	r42->Draw("same");

	TLatex* r43 = new TLatex(0.9, 0.91, "EPIC");
	r43->SetNDC();
	r43->SetTextSize(0.04);
	r43->Draw("same");

	TLatex* r44 = new TLatex(0.18, 0.84, "1<Q^{2}<10 GeV^{2}, 0.01 < y < 0.85");
	r44->SetNDC();
	r44->SetTextSize(20);
	r44->SetTextFont(43);
	r44->SetTextColor(kBlack);
	r44->Draw("same");

	TLatex* r44_0 = new TLatex(0.18, 0.79, Form("|y_{%s}|<3.5, |M_{inv} #minus M_{%s}| < 0.02 GeV", vm_label.Data(), vm_label.Data()));
	r44_0->SetNDC();
	r44_0->SetTextSize(20);
	r44_0->SetTextFont(43);
	r44_0->SetTextColor(kBlack);
	r44_0->Draw("same");

	TLatex* r44_2 = new TLatex(0.18, 0.18, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	r44_2->SetNDC();
	r44_2->SetTextSize(30);
	r44_2->SetTextFont(43);
	r44_2->SetTextColor(kBlack);
	r44_2->Draw("same");

	// Add legend
	TLegend* w7 = new TLegend(0.18, 0.68, 0.93, 0.76);
	w7->SetLineColor(kWhite);
	w7->SetFillColor(0);
	w7->SetTextSize(17);
	w7->SetTextFont(45);
	w7->AddEntry(h_t_MC, Form("Sartre %s MC: %.f events", vm_label.Data(),integral_MC), "L");
	w7->AddEntry(h_t_REC, Form("Sartre %s RECO w. EEMC: %.f events", vm_label.Data(),integral_REC), "P");
	w7->AddEntry(h_t_REC_wCUT[i], Form("Sartre %s RECO' #theta_{Max}= %s: %.f events", vm_label.Data(), angleLabels[i].Data(),integral_REC_wCut[i]), "P");
	w7->Draw("same");

	TLegend *w8 = new TLegend(0.48,0.64,0.93,0.56);
	w8->SetLineColor(kWhite);
	w8->SetFillColor(0);
	w8->SetTextSize(17);
	w8->SetTextFont(45);
	w8->AddEntry(h_t_REC_wCUT[i],"normalization: #int|#it{t}|_{MC}/#int|#it{t}|_{RECO'}","P");
	w8->Draw("same");

	// Save figure
	TString cleanLabel = angleLabels[i];
	cleanLabel.ReplaceAll("#", "");
	cleanLabel.ReplaceAll("/", "");
	Lnorm_canvases_Cut[i]->Print(Form("./figures/normTotal_wCUT_angle%s.pdf",cleanLabel.Data()));
}
	
	//2d plot truth
	TCanvas* totalNorm_canvases_2d[numCanvas];
	for (int i=0; i<1; i++) 
	{
    	TString canvasName = Form("totalNorm_2d_c%d", i+1);
    	totalNorm_canvases_2d[i] = new TCanvas(canvasName, canvasName, 1, 1, 600, 600);
    	gPad->SetLogz(1);
		gPad->SetLeftMargin(0.15);
		// Draw histograms with resolution
		h_t_REC_2d->Draw();
		// Save figure
		totalNorm_canvases_2d[i]->Print("./figures/2d.pdf");
	}

	//2d plots with resolution and cut
	TCanvas* totalNorm_canvases_2dresCut[numCanvas];
	for (int i=0; i<numCanvas; i++) 
	{
    	TString canvasName = Form("totalNorm_2dresCut_c%d", i+1);
    	totalNorm_canvases_2dresCut[i] = new TCanvas(canvasName, canvasName, 1, 1, 600, 600);
    	gPad->SetLogz(1);
		gPad->SetLeftMargin(0.15);
		// Draw histograms with different angle cuts
		h_t_REC_2d_wRES_cut[i]->Draw();
		// Add labels
		TLatex* r46 = new TLatex(0.45, 0.92, Form("#theta_{max}= %s",angleLabels[i].Data()));
		r46->SetNDC();
		r46->SetTextSize(15);
		r46->SetTextFont(43);
		r46->SetTextColor(kBlack);
		r46->Draw("same");
		// Save figure
		TString cleanLabel = angleLabels[i];
		cleanLabel.ReplaceAll("#", "");
		cleanLabel.ReplaceAll("/", "");
		totalNorm_canvases_2dresCut[i]->Print(Form("./figures/2d_wRES_cut_angle%s.pdf",cleanLabel.Data()));
	}
	
	//2d plot with resolution 
	TCanvas* totalNorm_canvases_2dres[numCanvas];
	for (int i=0; i<1; i++) 
	{
    	TString canvasName = Form("totalNorm_2dres_c%d", i+1);
    	totalNorm_canvases_2dres[i] = new TCanvas(canvasName, canvasName, 1, 1, 600, 600);
    	gPad->SetLogz(1);
		gPad->SetLeftMargin(0.15);
		// Draw histograms with resolution
		h_t_REC_2d_wRES->Draw();
		// Save figure
		totalNorm_canvases_2dres[i]->Print("./figures/2d_wRES.pdf");
	}
	
	//2d plots with cut
	TCanvas* Lnorm_canvases_2dCut[numCanvas];
	for (int i=0; i<numCanvas; i++) 
	{
    	TString canvasName = Form("Lnorm_2dCut_c%d", i+1);
    	Lnorm_canvases_2dCut[i] = new TCanvas(canvasName, canvasName, 1, 1, 600, 600);
    	gPad->SetLogz(1);
		gPad->SetLeftMargin(0.15);
		// Draw histograms with different angle cuts
		h_t_REC_2d_wCUT[i]->Draw();

		// Add labels
		TLatex* r46 = new TLatex(0.45, 0.92, Form("#theta_{max}= %s",angleLabels[i].Data()));
		r46->SetNDC();
		r46->SetTextSize(15);
		r46->SetTextFont(43);
		r46->SetTextColor(kBlack);
		r46->Draw("same");
		// Save figure
		TString cleanLabel = angleLabels[i];
		cleanLabel.ReplaceAll("#", "");
		cleanLabel.ReplaceAll("/", "");
		Lnorm_canvases_2dCut[i]->Print(Form("./figures/2d_wCUT_angle%s.pdf",cleanLabel.Data()));
	}

}