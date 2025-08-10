#include "RiceStyle.h"
using namespace std;
void plot_diffractive_event_kinematics_part2(TString filename)
{
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

	TH1D* h_MC_eta = (TH1D*) file->Get("h_MC_eta");//new TH1D("h_MC_eta",";#eta",100,-5,5);
	TH1D* h_eta = (TH1D*) file->Get("h_eta");//new TH1D("h_eta",";#eta",100,-5,5);
	TH1D* h_MC_phi = (TH1D*) file->Get("h_MC_phi");
	TH1D* h_Ecal_phi = (TH1D*) file->Get("h_Ecal_phi");

	TH1D* h_VM_mass = (TH1D*) file->Get("h_VM_mass");//new TH1D("h_VM_mass",";mass (GeV)",200,0,4);
	TH1D* h_VM_pt = (TH1D*) file->Get("h_VM_pt");//new TH1D("h_VM_pt",";p_{T} (GeV/c)",200,0,2);
	TH1D* h_VM_mass_REC = (TH1D*) file->Get("h_VM_mass_REC");//new TH1D("h_VM_mass_REC",";mass (GeV)",200,0,4);
	TH1D* h_VM_pt_REC = (TH1D*) file->Get("h_VM_pt_REC");//new TH1D("h_VM_pt_REC",";p_{T} (GeV/c)",200,0,2);

	TH1D* h_scat_phi_diff_cal = (TH1D*) file->Get("h_scat_phi_diff_cal");//new TH1D("h_scat_phi_diff_cal",";E #phi_{rec} - #phi_{MC}",100,0,3.14);
	TH2D* h_scat_phi_diff = (TH2D*) file->Get("h_scat_phi_diff");//new TH2D("h_scat_phi_diff",";E #phi_{rec}; #phi_{MC}",100,0,3.14,100,0,3.14);

    TCanvas* c1 = new TCanvas("c1","c1",1,1,1600,800);
    c1->Divide(4,2,0.01,0.01);
    c1->cd(1);
    gPad->SetLogy(1);
	gPad->SetLeftMargin(0.13);
	gPad->SetBottomMargin(0.15);
	h_MC_eta->GetXaxis()->SetTitleSize(0.8*h_MC_eta->GetXaxis()->GetTitleSize());
	h_MC_eta->GetXaxis()->SetLabelSize(0.8*h_MC_eta->GetXaxis()->GetLabelSize());
	h_MC_eta->GetYaxis()->SetTitleSize(0.8*h_MC_eta->GetYaxis()->GetTitleSize());
	h_MC_eta->GetYaxis()->SetLabelSize(0.8*h_MC_eta->GetYaxis()->GetLabelSize());
	h_MC_eta->GetXaxis()->SetTitleOffset(1.6*h_MC_eta->GetXaxis()->GetTitleOffset());
	h_MC_eta->GetYaxis()->SetTitleOffset(2.5*h_MC_eta->GetYaxis()->GetTitleOffset());
	h_MC_eta->GetYaxis()->SetTitle("counts");
	h_MC_eta->GetXaxis()->SetTitle("#eta");
    h_MC_eta->Draw();
    h_eta->SetMarkerStyle(24);
    h_eta->Draw("PEsame");
    TLegend *w7 = new TLegend(0.28,0.7,0.53,0.86);
	w7->SetLineColor(kWhite);
	w7->SetFillColor(0);
	w7->SetTextSize(17);
	w7->SetTextFont(45);
	w7->AddEntry(h_MC_eta, "MC ", "L");
	w7->AddEntry(h_eta, "RECO", "P");
	w7->Draw("same");



   /* c1->cd(3);
    //gPad->SetLogy(doLog_);
    gPad->SetLeftMargin(0.13);
	gPad->SetBottomMargin(0.15);
	h_MC_eta->GetXaxis()->SetTitleSize(0.8*h_MC_eta->GetXaxis()->GetTitleSize());
	h_MC_eta->GetXaxis()->SetLabelSize(0.8*h_MC_eta->GetXaxis()->GetLabelSize());
	h_MC_eta->GetYaxis()->SetTitleSize(0.8*h_MC_eta->GetYaxis()->GetTitleSize());
	h_MC_eta->GetYaxis()->SetLabelSize(0.8*h_MC_eta->GetYaxis()->GetLabelSize());
	h_MC_eta->GetXaxis()->SetTitleOffset(1.6*h_MC_eta->GetXaxis()->GetTitleOffset());
	h_MC_eta->GetYaxis()->SetTitleOffset(2.5*h_MC_eta->GetYaxis()->GetTitleOffset());
	h_MC_eta->GetYaxis()->SetTitle("counts");
	h_MC_eta->GetXaxis()->SetTitle("#eta");
    h_MC_eta->Draw();
    h_eta->SetMarkerStyle(24);
    h_eta->Draw("PEsame");
    w7->Draw("same");*/

    c1->cd(4);
    gPad->SetLogy(1);
    gPad->SetLeftMargin(0.13);
	gPad->SetBottomMargin(0.15);
	h_VM_mass->GetXaxis()->SetTitleSize(0.8*h_VM_mass->GetXaxis()->GetTitleSize());
	h_VM_mass->GetXaxis()->SetLabelSize(0.8*h_VM_mass->GetXaxis()->GetLabelSize());
	h_VM_mass->GetYaxis()->SetTitleSize(0.8*h_VM_mass->GetYaxis()->GetTitleSize());
	h_VM_mass->GetYaxis()->SetLabelSize(0.8*h_VM_mass->GetYaxis()->GetLabelSize());
	h_VM_mass->GetXaxis()->SetTitleOffset(1.6*h_VM_mass->GetXaxis()->GetTitleOffset());
	h_VM_mass->GetYaxis()->SetTitleOffset(2.5*h_VM_mass->GetYaxis()->GetTitleOffset());
	h_VM_mass->GetYaxis()->SetTitle("counts");
	h_VM_mass->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
	h_VM_mass->Draw();
    h_VM_mass_REC->SetMarkerStyle(24);
    h_VM_mass_REC->Draw("PEsame");
	w7->Draw("same");

    c1->cd(5);
    gPad->SetLogy(1);
    gPad->SetLeftMargin(0.13);
	gPad->SetBottomMargin(0.15);
	h_VM_pt->GetXaxis()->SetTitleSize(0.8*h_VM_pt->GetXaxis()->GetTitleSize());
	h_VM_pt->GetXaxis()->SetLabelSize(0.8*h_VM_pt->GetXaxis()->GetLabelSize());
	h_VM_pt->GetYaxis()->SetTitleSize(0.8*h_VM_pt->GetYaxis()->GetTitleSize());
	h_VM_pt->GetYaxis()->SetLabelSize(0.8*h_VM_pt->GetYaxis()->GetLabelSize());
	h_VM_pt->GetXaxis()->SetTitleOffset(1.6*h_VM_pt->GetXaxis()->GetTitleOffset());
	h_VM_pt->GetYaxis()->SetTitleOffset(2.5*h_VM_pt->GetYaxis()->GetTitleOffset());
	h_VM_pt->GetYaxis()->SetTitle("counts");
	h_VM_pt->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	h_VM_pt->Draw();
    h_VM_pt_REC->SetMarkerStyle(24);
    h_VM_pt_REC->Draw("PEsame");
	TLegend *w8 = new TLegend(0.68,0.2,0.83,0.46);
	w8->SetLineColor(kWhite);
	w8->SetFillColor(0);
	w8->SetTextSize(17);
	w8->SetTextFont(45);
	w8->AddEntry(h_VM_pt, "MC ", "L");
	w8->AddEntry(h_VM_pt_REC, "RECO", "P");
	w8->Draw("same");

    c1->cd(6);
    gPad->SetLogy(1);
    gPad->SetLeftMargin(0.13);
	gPad->SetBottomMargin(0.15);
	h_MC_phi->GetXaxis()->SetTitleSize(0.8*h_MC_phi->GetXaxis()->GetTitleSize());
	h_MC_phi->GetXaxis()->SetLabelSize(0.8*h_MC_phi->GetXaxis()->GetLabelSize());
	h_MC_phi->GetYaxis()->SetTitleSize(0.8*h_MC_phi->GetYaxis()->GetTitleSize());
	h_MC_phi->GetYaxis()->SetLabelSize(0.8*h_MC_phi->GetYaxis()->GetLabelSize());
	h_MC_phi->GetXaxis()->SetTitleOffset(1.6*h_MC_phi->GetXaxis()->GetTitleOffset());
	h_MC_phi->GetYaxis()->SetTitleOffset(2.5*h_MC_phi->GetYaxis()->GetTitleOffset());
	h_MC_phi->GetYaxis()->SetTitle("counts");
	h_MC_phi->GetXaxis()->SetTitle("#phi");
    h_MC_phi->Draw();
    h_Ecal_phi->SetMarkerStyle(24);
    h_Ecal_phi->Draw("PEsame");
    w7->Draw("same");

    c1->cd(7);
    gPad->SetLogy(1);
    gPad->SetLeftMargin(0.13);
	gPad->SetBottomMargin(0.15);
	h_scat_phi_diff_cal->GetXaxis()->SetTitleSize(0.8*h_scat_phi_diff_cal->GetXaxis()->GetTitleSize());
	h_scat_phi_diff_cal->GetXaxis()->SetLabelSize(0.8*h_scat_phi_diff_cal->GetXaxis()->GetLabelSize());
	h_scat_phi_diff_cal->GetYaxis()->SetTitleSize(0.8*h_scat_phi_diff_cal->GetYaxis()->GetTitleSize());
	h_scat_phi_diff_cal->GetYaxis()->SetLabelSize(0.8*h_scat_phi_diff_cal->GetYaxis()->GetLabelSize());
	h_scat_phi_diff_cal->GetXaxis()->SetTitleOffset(1.6*h_scat_phi_diff_cal->GetXaxis()->GetTitleOffset());
	h_scat_phi_diff_cal->GetYaxis()->SetTitleOffset(2.5*h_scat_phi_diff_cal->GetYaxis()->GetTitleOffset());
	h_scat_phi_diff_cal->GetYaxis()->SetTitle("counts");
	h_scat_phi_diff_cal->GetXaxis()->SetTitle("#phi_{RECO} - #phi_{MC}");
	h_scat_phi_diff_cal->Draw("PEsame");

    c1->cd(8);
    gPad->SetLogz(1);
    gPad->SetLeftMargin(0.13);
	gPad->SetBottomMargin(0.15);
	h_scat_phi_diff->GetXaxis()->SetTitleSize(0.8*h_scat_phi_diff->GetXaxis()->GetTitleSize());
	h_scat_phi_diff->GetXaxis()->SetLabelSize(0.8*h_scat_phi_diff->GetXaxis()->GetLabelSize());
	h_scat_phi_diff->GetYaxis()->SetTitleSize(0.8*h_scat_phi_diff->GetYaxis()->GetTitleSize());
	h_scat_phi_diff->GetYaxis()->SetLabelSize(0.8*h_scat_phi_diff->GetYaxis()->GetLabelSize());
	h_scat_phi_diff->GetXaxis()->SetTitleOffset(1.6*h_scat_phi_diff->GetXaxis()->GetTitleOffset());
	h_scat_phi_diff->GetYaxis()->SetTitleOffset(2.*h_scat_phi_diff->GetYaxis()->GetTitleOffset());
	h_scat_phi_diff->GetYaxis()->CenterTitle();
	h_scat_phi_diff->GetYaxis()->SetTitle("#phi_{MC}");
	h_scat_phi_diff->GetXaxis()->CenterTitle();
	h_scat_phi_diff->GetXaxis()->SetTitle("#phi_{RECO}");
    h_scat_phi_diff->Draw("colzsame");

	c1->Print("./figures/benchmark-phi-kinematics_part2.pdf");


}