#include "RiceStyle.h"
using namespace std;
void plot_diffractive_event_kinematics5(TString filename)
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

	TH2D* h_Ecal_phi2d = (TH2D*) file->Get("h_Ecal_phi2d");
	TProfile* p_Ecal_phi_vs_phi = h_Ecal_phi2d->ProfileX("p_Ecal_phi_vs_phi");
	TH1D* h_VM_Epz = (TH1D*) file->Get("h_VM_Epz");
	TH1D* h_VM_Epz_REC = (TH1D*) file->Get("h_VM_Epz_REC");
	TH1D* h_VM_pt = (TH1D*) file->Get("h_VM_pt");
	TH1D* h_VM_pt_REC = (TH1D*) file->Get("h_VM_pt_REC");
	TH1D* h_energy_calibration_REC = (TH1D*) file->Get("h_energy_calibration_REC");

	TH1D* h_t_MC = (TH1D*) file->Get("h_t_MC");
	TH1D* h_t_REC = (TH1D*) file->Get("h_t_REC");
	TH1D* h_t_REC_wRES_cut_pi12 = (TH1D*) file->Get("h_t_REC_wRES_cut_pi12");
	TH1D* h_t_ratio = (TH1D*) h_t_REC->Clone("h_t_ratio");
	h_t_ratio->Divide(h_t_MC);
	TH1D* h_t_ratio2 = (TH1D*) h_t_REC_wRES_cut_pi12->Clone("h_t_ratio2");
	h_t_ratio2->Divide(h_t_MC);



    TCanvas* c1 = new TCanvas("c1","c1",1,1,1600,800);
    c1->Divide(4,2,0.01,0.01);
    c1->cd(1);
    //gPad->SetLogy(doLog_);
    gPad->SetLeftMargin(0.13);
	gPad->SetBottomMargin(0.15);
	h_VM_Epz->GetXaxis()->SetTitleSize(0.8*h_VM_Epz->GetXaxis()->GetTitleSize());
	h_VM_Epz->GetXaxis()->SetLabelSize(0.8*h_VM_Epz->GetXaxis()->GetLabelSize());
	h_VM_Epz->GetYaxis()->SetTitleSize(0.8*h_VM_Epz->GetYaxis()->GetTitleSize());
	h_VM_Epz->GetYaxis()->SetLabelSize(0.8*h_VM_Epz->GetYaxis()->GetLabelSize());
	h_VM_Epz->GetXaxis()->SetTitleOffset(1.6*h_VM_Epz->GetXaxis()->GetTitleOffset());
	h_VM_Epz->GetYaxis()->SetTitleOffset(2.0*h_VM_Epz->GetYaxis()->GetTitleOffset());
	h_VM_Epz->GetYaxis()->SetTitle("counts");
	h_VM_Epz->GetXaxis()->SetTitle("E_{VM} - p_{z,VM} [GeV]");
    h_VM_Epz->Draw();
    h_VM_Epz_REC->SetMarkerStyle(24);
    h_VM_Epz_REC->Draw("PEsame");
    TLegend *w7 = new TLegend(0.28,0.7,0.53,0.86);
	w7->SetLineColor(kWhite);
	w7->SetFillColor(0);
	w7->SetTextSize(17);
	w7->SetTextFont(45);
	w7->AddEntry(h_VM_Epz, "MC ", "L");
	w7->AddEntry(h_VM_Epz_REC, "RECO", "P");
	w7->Draw("same");

    c1->cd(2);
    //gPad->SetLogy(doLog_);
    gPad->SetLeftMargin(0.13);
	gPad->SetBottomMargin(0.15);
	h_VM_pt->GetXaxis()->SetTitleSize(0.8*h_VM_pt->GetXaxis()->GetTitleSize());
	h_VM_pt->GetXaxis()->SetLabelSize(0.8*h_VM_pt->GetXaxis()->GetLabelSize());
	h_VM_pt->GetYaxis()->SetTitleSize(0.8*h_VM_pt->GetYaxis()->GetTitleSize());
	h_VM_pt->GetYaxis()->SetLabelSize(0.8*h_VM_pt->GetYaxis()->GetLabelSize());
	h_VM_pt->GetXaxis()->SetTitleOffset(1.6*h_VM_pt->GetXaxis()->GetTitleOffset());
	h_VM_pt->GetYaxis()->SetTitleOffset(2.0*h_VM_pt->GetYaxis()->GetTitleOffset());
	h_VM_pt->GetYaxis()->SetTitle("counts");
	h_VM_pt->GetXaxis()->SetTitle("p_{t} [GeV/c]");
    h_VM_pt->Draw();
    h_VM_pt_REC->SetMarkerStyle(24);
    h_VM_pt_REC->Draw("PEsame");
    w7->Draw("same");

    c1->cd(3);
    //gPad->SetLogy(doLog_);
    gPad->SetLeftMargin(0.13);
	gPad->SetBottomMargin(0.15);
	h_t_ratio->GetXaxis()->SetTitleSize(0.8*h_t_ratio->GetXaxis()->GetTitleSize());
	h_t_ratio->GetXaxis()->SetLabelSize(0.8*h_t_ratio->GetXaxis()->GetLabelSize());
	h_t_ratio->GetYaxis()->SetTitleSize(0.8*h_t_ratio->GetYaxis()->GetTitleSize());
	h_t_ratio->GetYaxis()->SetLabelSize(0.8*h_t_ratio->GetYaxis()->GetLabelSize());
	h_t_ratio->GetXaxis()->SetTitleOffset(1.6*h_t_ratio->GetXaxis()->GetTitleOffset());
	h_t_ratio->GetYaxis()->SetTitleOffset(2.5*h_t_ratio->GetYaxis()->GetTitleOffset());
	h_t_ratio->GetYaxis()->SetTitle("t_{L RECO}/t_{MC}");
	h_t_ratio->GetXaxis()->SetTitle("t_{MC} [GeV/c]^{2}");
    h_t_ratio->Draw();

    c1->cd(4);
    gPad->SetLogy(1);
    gPad->SetLeftMargin(0.13);
	gPad->SetBottomMargin(0.15);
	h_t_ratio2->GetXaxis()->SetTitleSize(0.8*h_t_ratio2->GetXaxis()->GetTitleSize());
	h_t_ratio2->GetXaxis()->SetLabelSize(0.8*h_t_ratio2->GetXaxis()->GetLabelSize());
	h_t_ratio2->GetYaxis()->SetTitleSize(0.8*h_t_ratio2->GetYaxis()->GetTitleSize());
	h_t_ratio2->GetYaxis()->SetLabelSize(0.8*h_t_ratio2->GetYaxis()->GetLabelSize());
	h_t_ratio2->GetXaxis()->SetTitleOffset(1.6*h_t_ratio2->GetXaxis()->GetTitleOffset());
	h_t_ratio2->GetYaxis()->SetTitleOffset(2.5*h_t_ratio2->GetYaxis()->GetTitleOffset());
	h_t_ratio2->GetYaxis()->SetTitle("t_{#theta RECO}/t_{MC}");
	h_t_ratio2->GetXaxis()->SetTitle("t_{MC} [GeV/c]^{2}");
    h_t_ratio2->SetMarkerStyle(24);
    h_t_ratio2->Draw("PEsame");

    c1->cd(5);
    gPad->SetLogy(1);
    gPad->SetLeftMargin(0.13);
	gPad->SetBottomMargin(0.15);
	h_energy_calibration_REC->GetXaxis()->SetTitleSize(0.8*h_energy_calibration_REC->GetXaxis()->GetTitleSize());
	h_energy_calibration_REC->GetXaxis()->SetLabelSize(0.8*h_energy_calibration_REC->GetXaxis()->GetLabelSize());
	h_energy_calibration_REC->GetYaxis()->SetTitleSize(0.8*h_energy_calibration_REC->GetYaxis()->GetTitleSize());
	h_energy_calibration_REC->GetYaxis()->SetLabelSize(0.8*h_energy_calibration_REC->GetYaxis()->GetLabelSize());
	h_energy_calibration_REC->GetXaxis()->SetTitleOffset(1.6*h_energy_calibration_REC->GetXaxis()->GetTitleOffset());
	h_energy_calibration_REC->GetYaxis()->SetTitleOffset(2.5*h_energy_calibration_REC->GetYaxis()->GetTitleOffset());
	h_energy_calibration_REC->GetYaxis()->SetTitle("counts");
	h_energy_calibration_REC->GetXaxis()->SetTitle("E_{RECO} / E_{MC}");
    h_energy_calibration_REC->SetMarkerStyle(24);
    h_energy_calibration_REC->Draw("PEsame");

    c1->cd(6);
    /*gPad->SetLogy(1);
    gPad->SetLeftMargin(0.13);
	gPad->SetBottomMargin(0.15);
	h_EoverP_REC->GetXaxis()->SetTitleSize(0.8*h_EoverP_REC->GetXaxis()->GetTitleSize());
	h_EoverP_REC->GetXaxis()->SetLabelSize(0.8*h_EoverP_REC->GetXaxis()->GetLabelSize());
	h_EoverP_REC->GetYaxis()->SetTitleSize(0.8*h_EoverP_REC->GetYaxis()->GetTitleSize());
	h_EoverP_REC->GetYaxis()->SetLabelSize(0.8*h_EoverP_REC->GetYaxis()->GetLabelSize());
	h_EoverP_REC->GetXaxis()->SetTitleOffset(1.6*h_EoverP_REC->GetXaxis()->GetTitleOffset());
	h_EoverP_REC->GetYaxis()->SetTitleOffset(2.5*h_EoverP_REC->GetYaxis()->GetTitleOffset());
	h_EoverP_REC->GetYaxis()->SetTitle("counts");
    h_EoverP_REC->SetMarkerStyle(24);
    h_EoverP_REC->Draw("PEsame");
*/
    c1->cd(7);
    gPad->SetLogy(1);
    gPad->SetLeftMargin(0.13);
	gPad->SetBottomMargin(0.15);
	p_Ecal_phi_vs_phi->GetXaxis()->SetTitleSize(0.8*p_Ecal_phi_vs_phi->GetXaxis()->GetTitleSize());
	p_Ecal_phi_vs_phi->GetXaxis()->SetLabelSize(0.8*p_Ecal_phi_vs_phi->GetXaxis()->GetLabelSize());
	p_Ecal_phi_vs_phi->GetYaxis()->SetTitleSize(0.8*p_Ecal_phi_vs_phi->GetYaxis()->GetTitleSize());
	p_Ecal_phi_vs_phi->GetYaxis()->SetLabelSize(0.8*p_Ecal_phi_vs_phi->GetYaxis()->GetLabelSize());
	p_Ecal_phi_vs_phi->GetXaxis()->SetTitleOffset(1.6*p_Ecal_phi_vs_phi->GetXaxis()->GetTitleOffset());
	p_Ecal_phi_vs_phi->GetYaxis()->SetTitleOffset(2.5*p_Ecal_phi_vs_phi->GetYaxis()->GetTitleOffset());
	p_Ecal_phi_vs_phi->GetYaxis()->SetTitle("#phi_{e,RECO}-#phi_{e,MC}");
	p_Ecal_phi_vs_phi->GetXaxis()->SetTitle("#phi_{e,RECO}");
    p_Ecal_phi_vs_phi->Draw("E1");

    c1->cd(8);
    gPad->SetLogz(1);
    gPad->SetLeftMargin(0.13);
	gPad->SetBottomMargin(0.15);
	h_Ecal_phi2d->GetXaxis()->SetTitleSize(0.8*h_Ecal_phi2d->GetXaxis()->GetTitleSize());
	h_Ecal_phi2d->GetXaxis()->SetLabelSize(0.8*h_Ecal_phi2d->GetXaxis()->GetLabelSize());
	h_Ecal_phi2d->GetYaxis()->SetTitleSize(0.8*h_Ecal_phi2d->GetYaxis()->GetTitleSize());
	h_Ecal_phi2d->GetYaxis()->SetLabelSize(0.8*h_Ecal_phi2d->GetYaxis()->GetLabelSize());
	h_Ecal_phi2d->GetXaxis()->SetTitleOffset(1.6*h_Ecal_phi2d->GetXaxis()->GetTitleOffset());
	h_Ecal_phi2d->GetYaxis()->SetTitleOffset(2.*h_Ecal_phi2d->GetYaxis()->GetTitleOffset());
	h_Ecal_phi2d->GetYaxis()->CenterTitle();
	h_Ecal_phi2d->GetYaxis()->SetTitle("#phi_{e,RECO}-#phi_{e,MC}");
	h_Ecal_phi2d->GetXaxis()->SetTitle("#phi_{e,RECO}");
    h_Ecal_phi2d->Draw("colzsame");

	c1->Print("./figures/benchmark-phi-kinematics5.pdf");


}