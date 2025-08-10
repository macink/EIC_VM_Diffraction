#include "RiceStyle.h"
using namespace std;
void plot_diffractive_event_kinematics_part4(TString filename)
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

	TH2D* h_Q2_res = (TH2D*) file->Get("h_Q2_res");
	TProfile* p_Q2_res_vs_Q2 = h_Q2_res->ProfileX("p_Q2_res_vs_Q2");
	TH1D* h_Q2_e = (TH1D*) file->Get("h_Q2_e");
	TH1D* h_Q2REC_e = (TH1D*) file->Get("h_Q2REC_e");
	TH1D* h_Q2_ratio = (TH1D*) h_Q2REC_e->Clone("h_Q2_ratio");
	h_Q2_ratio->Divide(h_Q2_e);
	TH2D* h_y_res = (TH2D*) file->Get("h_y_res");
	TProfile* p_y_res_vs_y = h_y_res->ProfileX("p_y_res_vs_y");
	TH1D* h_y_e = (TH1D*) file->Get("h_y_e");  
	TH1D* h_yREC_e = (TH1D*) file->Get("h_yREC_e");
	TH1D* h_y_ratio = (TH1D*) h_yREC_e->Clone("h_y_ratio");
	h_y_ratio->Divide(h_y_e);
	TH2D* h_Ecal_theta = (TH2D*) file->Get("h_Ecal_theta");
	TProfile* p_Ecal_theta_vs_theta = h_Ecal_theta->ProfileX("p_Ecal_theta_vs_theta");


    TCanvas* c1 = new TCanvas("c1","c1",1,1,1600,800);
    c1->Divide(4,2,0.01,0.01);
    c1->cd(1);
    gPad->SetLogy(1);
    gPad->SetLeftMargin(0.13);
	gPad->SetBottomMargin(0.15);
	h_Q2_ratio->GetXaxis()->SetTitleSize(0.8*h_Q2_ratio->GetXaxis()->GetTitleSize());
	h_Q2_ratio->GetXaxis()->SetLabelSize(0.8*h_Q2_ratio->GetXaxis()->GetLabelSize());
	h_Q2_ratio->GetYaxis()->SetTitleSize(0.8*h_Q2_ratio->GetYaxis()->GetTitleSize());
	h_Q2_ratio->GetYaxis()->SetLabelSize(0.8*h_Q2_ratio->GetYaxis()->GetLabelSize());
	h_Q2_ratio->GetXaxis()->SetTitleOffset(1.6*h_Q2_ratio->GetXaxis()->GetTitleOffset());
	h_Q2_ratio->GetYaxis()->SetTitleOffset(2.0*h_Q2_ratio->GetYaxis()->GetTitleOffset());
	h_Q2_ratio->GetYaxis()->SetTitle("Q^{2}_{e,RECO}/Q^{2}_{e,MC}");
	h_Q2_ratio->GetXaxis()->SetTitle("Q^{2}_{e,MC} [GeV/c]^{2}");
    h_Q2_ratio->Draw();

    c1->cd(2);
    gPad->SetLogy(1);
    gPad->SetLeftMargin(0.13);
	gPad->SetBottomMargin(0.15);
	h_y_ratio->GetXaxis()->SetTitleSize(0.8*h_y_ratio->GetXaxis()->GetTitleSize());
	h_y_ratio->GetXaxis()->SetLabelSize(0.8*h_y_ratio->GetXaxis()->GetLabelSize());
	h_y_ratio->GetYaxis()->SetTitleSize(0.8*h_y_ratio->GetYaxis()->GetTitleSize());
	h_y_ratio->GetYaxis()->SetLabelSize(0.8*h_y_ratio->GetYaxis()->GetLabelSize());
	h_y_ratio->GetXaxis()->SetTitleOffset(1.6*h_y_ratio->GetXaxis()->GetTitleOffset());
	h_y_ratio->GetYaxis()->SetTitleOffset(2.0*h_y_ratio->GetYaxis()->GetTitleOffset());
	h_y_ratio->GetYaxis()->SetTitle("y_{e,RECO}/y_{e,MC}");
	h_y_ratio->GetXaxis()->SetTitle("y_{e,MC}");
    h_y_ratio->Draw();
    
    c1->cd(3);
    gPad->SetLogy(1);
    gPad->SetLeftMargin(0.13);
	gPad->SetBottomMargin(0.15);
	p_Q2_res_vs_Q2->GetXaxis()->SetTitleSize(0.8*p_Q2_res_vs_Q2->GetXaxis()->GetTitleSize());
	p_Q2_res_vs_Q2->GetXaxis()->SetLabelSize(0.8*p_Q2_res_vs_Q2->GetXaxis()->GetLabelSize());
	p_Q2_res_vs_Q2->GetYaxis()->SetTitleSize(0.8*p_Q2_res_vs_Q2->GetYaxis()->GetTitleSize());
	p_Q2_res_vs_Q2->GetYaxis()->SetLabelSize(0.8*p_Q2_res_vs_Q2->GetYaxis()->GetLabelSize());
	p_Q2_res_vs_Q2->GetXaxis()->SetTitleOffset(1.6*p_Q2_res_vs_Q2->GetXaxis()->GetTitleOffset());
	p_Q2_res_vs_Q2->GetYaxis()->SetTitleOffset(2.5*p_Q2_res_vs_Q2->GetYaxis()->GetTitleOffset());
	p_Q2_res_vs_Q2->GetYaxis()->SetTitle("(Q^{2}_{e,MC}-Q^{2}_{REC})/Q^{2}_{e,MC}");
	p_Q2_res_vs_Q2->GetXaxis()->SetTitle("Q^{2}_{e,MC} [GeV/c]^{2}");
    p_Q2_res_vs_Q2->Draw("E1");
    


    c1->cd(4);
    gPad->SetLogy(1);
    gPad->SetLeftMargin(0.13);
	gPad->SetBottomMargin(0.15);
	p_y_res_vs_y->GetXaxis()->SetTitleSize(0.8*p_y_res_vs_y->GetXaxis()->GetTitleSize());
	p_y_res_vs_y->GetXaxis()->SetLabelSize(0.8*p_y_res_vs_y->GetXaxis()->GetLabelSize());
	p_y_res_vs_y->GetYaxis()->SetTitleSize(0.8*p_y_res_vs_y->GetYaxis()->GetTitleSize());
	p_y_res_vs_y->GetYaxis()->SetLabelSize(0.8*p_y_res_vs_y->GetYaxis()->GetLabelSize());
	p_y_res_vs_y->GetXaxis()->SetTitleOffset(1.6*p_y_res_vs_y->GetXaxis()->GetTitleOffset());
	p_y_res_vs_y->GetYaxis()->SetTitleOffset(2.5*p_y_res_vs_y->GetYaxis()->GetTitleOffset());
	p_y_res_vs_y->GetYaxis()->SetTitle("y_{e,MC}-y_{REC}/y_{e,MC}");
	p_y_res_vs_y->GetXaxis()->SetTitle("y_{e,MC}");
	p_y_res_vs_y->Draw("E1");

    c1->cd(5);
    gPad->SetLogy(1);
    gPad->SetLeftMargin(0.13);
	gPad->SetBottomMargin(0.15);
	p_Ecal_theta_vs_theta->GetXaxis()->SetTitleSize(0.8*p_Ecal_theta_vs_theta->GetXaxis()->GetTitleSize());
	p_Ecal_theta_vs_theta->GetXaxis()->SetLabelSize(0.8*p_Ecal_theta_vs_theta->GetXaxis()->GetLabelSize());
	p_Ecal_theta_vs_theta->GetYaxis()->SetTitleSize(0.8*p_Ecal_theta_vs_theta->GetYaxis()->GetTitleSize());
	p_Ecal_theta_vs_theta->GetYaxis()->SetLabelSize(0.8*p_Ecal_theta_vs_theta->GetYaxis()->GetLabelSize());
	p_Ecal_theta_vs_theta->GetXaxis()->SetTitleOffset(1.6*p_Ecal_theta_vs_theta->GetXaxis()->GetTitleOffset());
	p_Ecal_theta_vs_theta->GetYaxis()->SetTitleOffset(2.5*p_Ecal_theta_vs_theta->GetYaxis()->GetTitleOffset());
	p_Ecal_theta_vs_theta->GetYaxis()->SetTitle("#theta_{e,REC}-#theta_{e,MC}");
	p_Ecal_theta_vs_theta->GetXaxis()->SetTitle("#theta_{e,MC}");
	p_Ecal_theta_vs_theta->Draw("E1");

    c1->cd(6);
    gPad->SetLogz(1);
    gPad->SetLeftMargin(0.13);
	gPad->SetBottomMargin(0.15);
	h_Ecal_theta->GetXaxis()->SetTitleSize(0.8*h_Ecal_theta->GetXaxis()->GetTitleSize());
	h_Ecal_theta->GetXaxis()->SetLabelSize(0.8*h_Ecal_theta->GetXaxis()->GetLabelSize());
	h_Ecal_theta->GetYaxis()->SetTitleSize(0.8*h_Ecal_theta->GetYaxis()->GetTitleSize());
	h_Ecal_theta->GetYaxis()->SetLabelSize(0.8*h_Ecal_theta->GetYaxis()->GetLabelSize());
	h_Ecal_theta->GetXaxis()->SetTitleOffset(1.6*h_Ecal_theta->GetXaxis()->GetTitleOffset());
	h_Ecal_theta->GetYaxis()->SetTitleOffset(2.*h_Ecal_theta->GetYaxis()->GetTitleOffset());
	h_Ecal_theta->GetYaxis()->CenterTitle();
	h_Ecal_theta->GetYaxis()->SetTitle("#theta_{e,REC}-#theta_{e,MC}");
	h_Ecal_theta->GetXaxis()->CenterTitle();
	h_Ecal_theta->GetXaxis()->SetTitle("#theta_{e,MC}");
    h_Ecal_theta->Draw("colzsame");
    
    c1->cd(7);
    gPad->SetLogz(1);
    gPad->SetLeftMargin(0.13);
	gPad->SetBottomMargin(0.15);
	h_y_res->GetXaxis()->SetTitleSize(0.8*h_y_res->GetXaxis()->GetTitleSize());
	h_y_res->GetXaxis()->SetLabelSize(0.8*h_y_res->GetXaxis()->GetLabelSize());
	h_y_res->GetYaxis()->SetTitleSize(0.8*h_y_res->GetYaxis()->GetTitleSize());
	h_y_res->GetYaxis()->SetLabelSize(0.8*h_y_res->GetYaxis()->GetLabelSize());
	h_y_res->GetXaxis()->SetTitleOffset(1.6*h_y_res->GetXaxis()->GetTitleOffset());
	h_y_res->GetYaxis()->SetTitleOffset(2.*h_y_res->GetYaxis()->GetTitleOffset());
	h_y_res->GetYaxis()->CenterTitle();
	h_y_res->GetYaxis()->SetTitle("y_{e,MC}-y_{REC}/y_{e,MC}");
	h_y_res->GetXaxis()->CenterTitle();
	h_y_res->GetXaxis()->SetTitle("y_{e,MC}");
    h_y_res->Draw("colzsame");


    c1->cd(8);
    gPad->SetLogz(1);
    gPad->SetLeftMargin(0.13);
	gPad->SetBottomMargin(0.15);
	h_Q2_res->GetXaxis()->SetTitleSize(0.8*h_Q2_res->GetXaxis()->GetTitleSize());
	h_Q2_res->GetXaxis()->SetLabelSize(0.8*h_Q2_res->GetXaxis()->GetLabelSize());
	h_Q2_res->GetYaxis()->SetTitleSize(0.8*h_Q2_res->GetYaxis()->GetTitleSize());
	h_Q2_res->GetYaxis()->SetLabelSize(0.8*h_Q2_res->GetYaxis()->GetLabelSize());
	h_Q2_res->GetXaxis()->SetTitleOffset(1.6*h_Q2_res->GetXaxis()->GetTitleOffset());
	h_Q2_res->GetYaxis()->SetTitleOffset(2.*h_Q2_res->GetYaxis()->GetTitleOffset());
	h_Q2_res->GetYaxis()->CenterTitle();
	h_Q2_res->GetYaxis()->SetTitle("(Q^{2}_{e,MC}-Q^{2}_{REC})/Q^{2}_{e,MC}");
	h_Q2_res->GetXaxis()->CenterTitle();
	h_Q2_res->GetXaxis()->SetTitle("Q^{2}_{e,MC} [GeV/c}^{2}]");
    h_Q2_res->Draw("colzsame");


	c1->Print("./figures/benchmark-phi-kinematics_part4.pdf");


}