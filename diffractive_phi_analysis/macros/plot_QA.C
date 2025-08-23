#include "ePICStyle.C"
gROOT->ProcessLine("SetePICStyle()");

using namespace std;

/*----------------------------------------------------------------
	Generates analysis plots from the generated diffractive
    VM analysis root file

	Input: root file
	Output: saves plots as pdfs
-----------------------------------------------------------------*/

void plot_QA(TString filename)
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
		
// Q2
	TH1D* h_Q2_e = (TH1D*) file->Get("h_Q2_e"); // after MC cut
    TH1D* h_Q2_beforeCut = (TH1D*) file->Get("h_Q2_beforeCut"); // before MC cut
    TH1D* h_Q2REC_e_EEMC = (TH1D*) file->Get("h_Q2REC_e_EEMC"); // before reco cut
    TH1D* h_Q2_afterCut = (TH1D*) file->Get("h_Q2_afterCut"); // after reco cut
    TH1D* h_Q2REC_e_EEMC_cut = (TH1D*) file->Get("h_Q2REC_e_EEMC_cut"); // cut inside position threshold
    TH2D* h_Q2_res = (TH2D*) file->Get("h_Q2_res"); // reco (after cut)
    TH2D* h_Q2_res_cut = (TH2D*) file->Get("h_Q2_res_cut"); // cut inside position threshold
    TH2D* h_Q2_response = (TH2D*) file->Get("h_Q2_response"); // reco (after cut)
    TH2D* h_Q2_response_cut = (TH2D*) file->Get("h_Q2_response_cut"); // cut inside position threshold
    TH1D* h_dQ2overQ2_REC = (TH1D*) file->Get("h_dQ2overQ2_REC"); // reco (after cut)
	TH1D* h_dQ2overQ2_REC_cut = (TH1D*) file->Get("h_dQ2overQ2_REC_cut"); // cut inside position threshold

// y
	TH1D* h_y_e = (TH1D*) file->Get("h_y_e"); // after MC cut
    TH1D* h_y_beforeCut = (TH1D*) file->Get("h_y_beforeCut"); // before MC cut
	TH1D* h_yREC_e_EEMC = (TH1D*) file->Get("h_yREC_e_EEMC"); // before reco cut
    TH1D* h_y_afterCut = (TH1D*) file->Get("h_y_afterCut"); // after reco cut
    TH1D* h_yREC_e_EEMC_cut = (TH1D*) file->Get("h_yREC_e_EEMC_cut"); // cut inside position threshold
    TH2D* h_y_res = (TH2D*) file->Get("h_y_res"); // reco (after cut)
    TH2D* h_y_res_cut = (TH2D*) file->Get("h_y_res_cut"); // cut inside position threshold
    TH2D* h_y_response = (TH2D*) file->Get("h_y_response"); // reco (after cut)
    TH2D* h_y_response_cut = (TH2D*) file->Get("h_y_response_cut"); // cut inside position threshold
    TH1D* h_dyOvery_REC = (TH1D*) file->Get("h_dyOvery_REC"); // reco (after cut)
	TH1D* h_dyOvery_REC_cut = (TH1D*) file->Get("h_dyOvery_REC_cut"); // cut inside position threshold
       
// energy
	TH1D* h_energy_MC = (TH1D*) file->Get("h_energy_MC"); // before MC cut
    TH1D* h_energy_MC_after = (TH1D*) file->Get("h_energy_MC_after"); // after MC cut
	TH1D* h_energy_REC_EEMC = (TH1D*) file->Get("h_energy_REC_EEMC"); // before reco cut
	TH1D* h_energy_REC_EEMC_cut = (TH1D*) file->Get("h_energy_REC_EEMC_cut"); // cut inside position threshold
    TH2D* h_energy_res_EEMC = (TH2D*) file->Get("h_energy_res_EEMC"); // before reco cut
    TH2D* h_energy_res_EEMC_cut = (TH2D*) file->Get("h_energy_res_EEMC_cut"); // cut inside position threshold
	TH2D* h_energy_response_EEMC_cut = (TH2D*) file->Get("h_energy_response_EEMC_cut"); // cut inside position threshold
	TH2D* h_energy_response_EEMC = (TH2D*) file->Get("h_energy_response_EEMC"); // before reco cut

// e' p
    TH1D* h_e_pt_MC = (TH1D*) file->Get("h_e_pt_MC"); // before MC cut
    TH1D* h_e_pt_MC_after = (TH1D*) file->Get("h_e_pt_MC_after"); // after MC cut
    TH1D* h_e_pt_REC_EEMC = (TH1D*) file->Get("h_e_pt_REC_EEMC"); // reco (after cut)
    TH1D* h_e_pt_REC_EEMC_cut = (TH1D*) file->Get("h_e_pt_REC_EEMC_cut"); // cut inside position threshold
    TH1D* h_e_pt_REC_trk = (TH1D*) file->Get("h_e_pt_REC_trk"); // reco (after cut)
    TH1D* h_e_pt_REC_trk_cut = (TH1D*) file->Get("h_e_pt_REC_trk_cut"); // cut inside position threshold
    TH2D* h_e_pt_res = (TH2D*) file->Get("h_e_pt_res"); // reco (after cut)
    TH2D* h_e_pt_res_cut = (TH2D*) file->Get("h_e_pt_res_cut"); // cut inside position threshold
	TH1D* h_e_pz_MC = (TH1D*) file->Get("h_e_pz_MC"); // before MC cut
    TH1D* h_e_pz_MC_after = (TH1D*) file->Get("h_e_pz_MC_after"); // after MC cut
    TH1D* h_e_pz_REC_EEMC = (TH1D*) file->Get("h_e_pz_REC_EEMC"); // reco (after cut)
    TH1D* h_e_pz_REC_EEMC_cut = (TH1D*) file->Get("h_e_pz_REC_EEMC_cut"); // cut inside position threshold
    TH1D* h_e_pz_REC_trk = (TH1D*) file->Get("h_e_pz_REC_trk"); // reco (after cut)
    TH1D* h_e_pz_REC_trk_cut = (TH1D*) file->Get("h_e_pz_REC_trk_cut"); // cut inside position threshold
    TH2D* h_e_pz_res = (TH2D*) file->Get("h_e_pz_res"); // reco (after cut)
    TH2D* h_e_pz_res_cut = (TH2D*) file->Get("h_e_pz_res_cut"); // cut inside position threshold
   	TH2D* h_e_pz_response = (TH2D*) file->Get("h_e_pz_response"); // reco (after cut)
    TH2D* h_e_pz_response_cut = (TH2D*) file->Get("h_e_pz_response_cut"); // cut inside position threshold
	TH1D* h_e_p_MC = (TH1D*) file->Get("h_e_p_MC"); // before MC cut
    TH1D* h_e_p_MC_after = (TH1D*) file->Get("h_e_p_MC_after"); // after MC cut
	TH1D* h_e_p_REC_EEMC = (TH1D*) file->Get("h_e_p_REC_EEMC"); // reco (after cut)
    TH1D* h_e_p_REC_EEMC_cut = (TH1D*) file->Get("h_e_p_REC_EEMC_cut"); // cut inside position threshold
    TH1D* h_e_p_REC_trk = (TH1D*) file->Get("h_e_p_REC_trk"); // reco (after cut)
    TH1D* h_e_p_REC_trk_cut = (TH1D*) file->Get("h_e_p_REC_trk_cut"); // cut inside position threshold
    TH2D* h_e_p_res = (TH2D*) file->Get("h_e_p_res"); // reco (after cut)
    TH2D* h_e_p_res_cut = (TH2D*) file->Get("h_e_p_res_cut"); // cut inside position threshold

// theta
    TH1D* h_theta_MC = (TH1D*) file->Get("h_theta_MC"); // before MC cut
    TH1D* h_theta_MC_after = (TH1D*) file->Get("h_theta_MC_after"); // after MC cut
	TH1D* h_theta_REC_EEMC_cut = (TH1D*) file->Get("h_theta_REC_EEMC_cut"); // cut inside position threshold
    TH1D* h_theta_REC_EEMC_after = (TH1D*) file->Get("h_theta_REC_EEMC_after"); // reco (after cut)
    TH1D* h_theta_diff_cut = (TH1D*) file->Get("h_theta_diff_cut"); // cut inside position threshold
    TH1D* h_theta_diff_after = (TH1D*) file->Get("h_theta_diff_after"); // reco (after cut)
	TH2D* h_theta_response_EEMC_cut = (TH2D*) file->Get("h_theta_response_EEMC_cut"); // cut inside position threshold
    TH2D* h_theta_response_EEMC_after = (TH2D*) file->Get("h_theta_response_EEMC_after"); // reco (after cut)
    TH2D* h_theta_response_EEMC_cut_clone = (TH2D*)h_theta_response_EEMC_cut->Clone("h_theta_response_EEMC_cut_clone");

// eta
    TH1D* h_eta_MC = (TH1D*) file->Get("h_eta_MC"); // before MC cut
    TH1D* h_eta_MC_after = (TH1D*) file->Get("h_eta_MC_after"); // after MC cut
	TH1D* h_eta_REC_EEMC_cut = (TH1D*) file->Get("h_eta_REC_EEMC_cut"); // cut inside position threshold
    TH1D* h_eta_REC_EEMC = (TH1D*) file->Get("h_eta_REC_EEMC"); // before reco cut
    TH1D* h_eta_REC_EEMC_after = (TH1D*) file->Get("h_eta_REC_EEMC_after"); // reco (after cut)
    TH1D* h_eta_diff_cut = (TH1D*) file->Get("h_eta_diff_cut"); // cut inside position threshold
    TH1D* h_eta_diff = (TH1D*) file->Get("h_eta_diff"); // before reco cut
    TH1D* h_eta_diff_after = (TH1D*) file->Get("h_eta_diff_after"); // reco (after cut)

// phi
    TH1D* h_phi_MC = (TH1D*) file->Get("h_phi_MC"); // before MC cut
    TH1D* h_phi_MC_after = (TH1D*) file->Get("h_phi_MC_after"); // after MC cut
	TH1D* h_phi_REC_EEMC_cut = (TH1D*) file->Get("h_phi_REC_EEMC_cut"); // cut inside position threshold
    TH1D* h_phi_REC_EEMC = (TH1D*) file->Get("h_phi_REC_EEMC"); // before reco cut
    TH1D* h_phi_REC_EEMC_after = (TH1D*) file->Get("h_phi_REC_EEMC_after"); // reco (after cut)
	TH1D* h_phi_diff_cut = (TH1D*) file->Get("h_phi_diff_cut"); // cut inside position threshold
    TH1D* h_phi_diff = (TH1D*) file->Get("h_phi_diff"); // before reco cut
    TH1D* h_phi_diff_after = (TH1D*) file->Get("h_phi_diff_after"); // reco (after cut)
	
// E and p
	TH2D* h_EvsP_MC = (TH2D*) file->Get("h_EvsP_MC"); 
	TH2D* h_EvsP_REC_cut = (TH2D*) file->Get("h_EvsP_REC_cut"); // cut inside position threshold
	TH2D* h_EvsP_REC = (TH2D*) file->Get("h_EvsP_REC"); // reco (after cut)
    
// E-pz
    TH1D* h_Epz_MC = (TH1D*) file->Get("h_Epz_MC");
	TH1D* h_Epz_REC_cut = (TH1D*) file->Get("h_Epz_REC_cut"); // cut inside position threshold
    TH1D* h_Epz_REC = (TH1D*) file->Get("h_Epz_REC"); // reco (after cut)
    TH1D* h_Epz_afterCut = (TH1D*) file->Get("h_Epz_afterCut"); // cut inside position threshold
    TH2D* h_Epz_response = (TH2D*) file->Get("h_Epz_response"); // reco (after cut)
    TH2D* h_Epz_response_cut = (TH2D*) file->Get("h_Epz_response_cut"); // cut inside position threshold

// E/p
    TH1D* h_EoverP_MC = (TH1D*) file->Get("h_EoverP_MC"); // before MC cut
    TH1D* h_EoverP_MC_after = (TH1D*) file->Get("h_EoverP_MC_after"); // after MC cut
	TH1D* h_EoverP_afterCut = (TH1D*) file->Get("h_EoverP_afterCut"); // after reco cut
	TH1D* h_EcalOverPtrk_cut = (TH1D*) file->Get("h_EcalOverPtrk_cut"); // cut inside position threshold
	TH1D* h_EcalOverPtrk = (TH1D*) file->Get("h_EcalOverPtrk"); // reco (after cut)
	TH2D* h_EoverP_response_cut = (TH2D*) file->Get("h_EoverP_response_cut"); // cut inside position threshold
    TH2D* h_EoverP_response = (TH2D*) file->Get("h_EoverP_response"); // reco (after cut)

// VM
	TH1D* h_VM_mass_MC = (TH1D*) file->Get("h_VM_mass_MC");
    TH1D* h_VM_mass_REC = (TH1D*) file->Get("h_VM_mass_REC");
    TH1D* h_VM_mass_REC_cut = (TH1D*) file->Get("h_VM_mass_REC_cut"); // cut inside position threshold
	TH1D* h_VM_pt_MC = (TH1D*) file->Get("h_VM_pt_MC");
    TH1D* h_VM_pt_REC = (TH1D*) file->Get("h_VM_pt_REC");
	TH1D* h_VM_pt_REC_cut = (TH1D*) file->Get("h_VM_pt_REC_cut"); // cut inside position threshold
    TH2D* h_VM_pt_response = (TH2D*) file->Get("h_VM_pt_response");
    TH2D* h_VM_pt_response_cut = (TH2D*) file->Get("h_VM_pt_response_cut"); // cut inside position threshold
	TH1D* h_VM_pz_MC = (TH1D*) file->Get("h_VM_pz_MC");
    TH1D* h_VM_pz_REC = (TH1D*) file->Get("h_VM_pz_REC");
    TH1D* h_VM_pz_REC_cut = (TH1D*) file->Get("h_VM_pz_REC_cut"); // cut inside position threshold
	TH1D* h_VM_p_MC = (TH1D*) file->Get("h_VM_p_MC");
    TH1D* h_VM_p_REC = (TH1D*) file->Get("h_VM_p_REC");
    TH1D* h_VM_p_REC_cut = (TH1D*) file->Get("h_VM_p_REC_cut"); // cut inside position threshold
	TH1D* h_VM_Epz_REC_cut = (TH1D*) file->Get("h_VM_Epz_REC_cut"); // cut inside position threshold
	TH1D* h_VM_Epz_MC = (TH1D*) file->Get("h_VM_Epz_MC");
    TH1D* h_VM_Epz_REC = (TH1D*) file->Get("h_VM_Epz_REC");
    TH2D* h_VM_Epz_response = (TH2D*) file->Get("h_VM_Epz_response");
	TH2D* h_VM_Epz_response_cut = (TH2D*) file->Get("h_VM_Epz_response_cut"); // cut inside position threshold

// position
    TH1D* h_Xclus_minus_Xtrk = (TH1D*) file->Get("h_Xclus_minus_Xtrk"); // before reco cut
    TH1D* h_Xclus_minus_Xtrk_cut = (TH1D*) file->Get("h_Xclus_minus_Xtrk_cut"); // cut inside position threshold
    TH1D* h_Xclus_minus_Xtrk_cut2 = (TH1D*) file->Get("h_Xclus_minus_Xtrk_cut2"); // cut outside position threshold
	TH1D* h_Yclus_minus_Ytrk = (TH1D*) file->Get("h_Yclus_minus_Ytrk"); // before reco cut
	TH1D* h_Yclus_minus_Ytrk_cut = (TH1D*) file->Get("h_Yclus_minus_Ytrk_cut"); // cut inside position threshold
    TH1D* h_Yclus_minus_Ytrk_cut2 = (TH1D*) file->Get("h_Yclus_minus_Ytrk_cut2"); // cut outside position threshold
    TH2D* h_emClus_position_REC = (TH2D*)file->Get("h_emClus_position_REC");
    TH2D* h_emClus_position_REC_cut = (TH2D*)file->Get("h_emClus_position_REC_cut");


// t
	TH1D* h_t_MC = (TH1D*) file->Get("h_t_MC");
    TH1D* h_t_REC_EEMC = (TH1D*) file->Get("h_t_REC_EEMC"); // after MC cuts
	TH1D* h_t_REC_EEMC_cut = (TH1D*) file->Get("h_t_REC_EEMC_cut"); // cut inside position threshold
    TH1D* h_t_REC_EEMC_cut2 = (TH1D*) file->Get("h_t_REC_EEMC_cut2"); // cut outside position threshold
	TH2D* h_t_res_EEMC_cut = (TH2D*) file->Get("h_t_res_EEMC_cut"); // cut inside position threshold
    TH2D* h_t_res_EEMC_cut_percent = (TH2D*) file->Get("h_t_res_EEMC_cut_percent"); // cut inside position threshold
    TH2D* h_t_res_EEMC = (TH2D*) file->Get("h_t_res_EEMC"); // after MC cuts
	TH2D* h_t_response_EEMC_cut = (TH2D*) file->Get("h_t_response_EEMC_cut"); // cut inside position threshold
    TH2D* h_t_response_EEMC = (TH2D*) file->Get("h_t_response_EEMC"); // after MC cuts

// t distribution with resolution and angle cut
    TH1D* h_t_REC_wRES_cut_pi9 = (TH1D*) file->Get("h_t_REC_wRES_cut_pi9");
	TH1D* h_t_REC_wRES_cut_pi9_cut = (TH1D*) file->Get("h_t_REC_wRES_cut_pi9_cut"); // cut inside position threshold
    TH1D* h_t_REC_wRES_cut_pi9_cut2 = (TH1D*) file->Get("h_t_REC_wRES_cut_pi9_cut2"); // cut outside position threshold
    TH1D* h_t_REC_wRES_cut_pi12 = (TH1D*) file->Get("h_t_REC_wRES_cut_pi12");
	TH1D* h_t_REC_wRES_cut_pi12_cut = (TH1D*) file->Get("h_t_REC_wRES_cut_pi12_cut"); // cut inside position threshold
    TH1D* h_t_REC_wRES_cut_pi12_cut2 = (TH1D*) file->Get("h_t_REC_wRES_cut_pi12_cut2"); // cut outside position threshold
    TH1D* h_t_REC_wRES_cut_pi16 = (TH1D*) file->Get("h_t_REC_wRES_cut_pi16"); 
	TH1D* h_t_REC_wRES_cut_pi16_cut = (TH1D*) file->Get("h_t_REC_wRES_cut_pi16_cut"); // cut inside position threshold
    TH1D* h_t_REC_wRES_cut_pi16_cut2 = (TH1D*) file->Get("h_t_REC_wRES_cut_pi16_cut2"); // cut outside position threshold

// t distribution 2d
    TH2D* h_t_REC_2d = (TH2D*) file->Get("h_t_REC_2d"); 
	TH2D* h_t_REC_2d_cut = (TH2D*) file->Get("h_t_REC_2d_cut"); // cut inside position 
    TH2D* h_t_REC_2d_cut2 = (TH2D*) file->Get("h_t_REC_2d_cut2"); // cut outside position 
    
// t distribution 2d with resolution only 
    TH2D* h_t_REC_2d_wRES = (TH2D*) file->Get("h_t_REC_2d_wRES");
	TH2D* h_t_REC_2d_wRES_cut = (TH2D*) file->Get("h_t_REC_2d_wRES_cut"); // cut inside position threshold
    TH2D* h_t_REC_2d_wRES_cut2 = (TH2D*) file->Get("h_t_REC_2d_wRES_cut2"); // cut outside position threshold


// Q2
	TCanvas* c1 = new TCanvas("c1","c1",1,1,1600,800);
    c1->Divide(3,1,0.01,0.01);
    c1->cd(1);
    gPad->SetLogy(1);
	h_Q2_e->GetYaxis()->SetTitle("counts");
	h_Q2_e->GetXaxis()->SetTitle("Q^{2} [GeV/c]^{2}");
	h_Q2_e->SetLineColor(kBlack);
    h_Q2_e->Draw();
    //h_Q2_beforeCut->SetMarkerStyle(11);
    //h_Q2_beforeCut->SetMarkerColor(kOrange);
    //h_Q2_beforeCut->Draw("PEsame");
    //h_Q2REC_e_EEMC->SetMarkerColor(kBlue);
    //h_Q2REC_e_EEMC->SetMarkerStyle(24);
    //h_Q2REC_e_EEMC->Draw("PEsame");
    //h_Q2_afterCut->SetMarkerStyle(29);
    //h_Q2_afterCut->SetMarkerColor(kRed);
    //h_Q2_afterCut->Draw("PEsame");
    h_Q2REC_e_EEMC_cut->SetMarkerStyle(24);
    h_Q2REC_e_EEMC_cut->SetMarkerColor(kBlue);
    h_Q2REC_e_EEMC_cut->Draw("PEsame");
    TLegend *w1 = new TLegend(0.42,0.7,0.6,0.85);
	w1->AddEntry(h_Q2_e, "MC", "L");
    //w1->AddEntry(h_Q2_beforeCut, "MC (before cut)", "P");
	//w1->AddEntry(h_Q2REC_e_EEMC, "reco (before cut)", "P");
    //w1->AddEntry(h_Q2_afterCut, "reco", "P");
    w1->AddEntry(h_Q2REC_e_EEMC_cut, "RECO", "P");
	w1->Draw("same");

    /*c1->cd(2);
    gPad->SetLogz(1);
    h_Q2_res->GetYaxis()->SetTitle("(Q^{2}_{MC}-Q^{2}_{RECO})/Q^{2}_{MC}");
    h_Q2_res->GetXaxis()->SetTitle("Q^{2}_{MC} [GeV/c]^{2}");
    h_Q2_res->Draw("colzsame");*/

    c1->cd(2);
    gPad->SetLogz(1);
    h_Q2_res_cut->GetYaxis()->SetTitleOffset(1.5);
    h_Q2_res_cut->GetYaxis()->SetTitle("(Q^{2}_{MC}-Q^{2}_{RECO})/Q^{2}_{MC}");
    h_Q2_res_cut->GetXaxis()->SetTitle("Q^{2}_{MC} [GeV/c]^{2}");
    h_Q2_res_cut->Draw("colzsame");

    /*c1->cd(4);
    gPad->SetLogz(1);
    h_Q2_response->GetYaxis()->SetTitle("Q^{2}_{MC} [GeV/c]^{2}");
    h_Q2_response->GetXaxis()->SetTitle("Q^{2}_{RECO} [GeV/c]^{2}");
    h_Q2_response->Draw("colzsame");*/

    c1->cd(3);
    gPad->SetLogz(1);
    h_Q2_response_cut->GetYaxis()->SetTitleOffset(1.5);
    h_Q2_response_cut->GetYaxis()->SetTitle("Q^{2}_{MC} [GeV/c]^{2}");
    h_Q2_response_cut->GetXaxis()->SetTitle("Q^{2}_{RECO} [GeV/c]^{2}");
    h_Q2_response_cut->Draw("colzsame");

    /*c1->cd(6);
    gPad->SetLogy(1);
    h_dQ2overQ2_REC->GetYaxis()->SetTitle("counts");
    h_dQ2overQ2_REC->GetXaxis()->SetTitle("dQ^{2}/Q^{2}_{MC}");
    h_dQ2overQ2_REC->Draw();

    c1->cd(7);
    gPad->SetLogy(1);
    h_dQ2overQ2_REC_cut->GetYaxis()->SetTitle("counts");
    h_dQ2overQ2_REC_cut->GetXaxis()->SetTitle("dQ^{2}/Q^{2}_{MC}");
    h_dQ2overQ2_REC_cut->Draw();*/

    c1->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title1_1 = new TLatex();
    title1_1->SetNDC(); 
    title1_1->SetTextSize(0.05);
    title1_1->SetTextAlign(22);  
    title1_1->DrawLatex(0.5, 0.97, "Q^{2} Truth vs Reco");  

    c1->cd(2);  
    gPad->SetTopMargin(0.08);  
    TLatex* title1_2 = new TLatex();
    title1_2->SetNDC(); 
    title1_2->SetTextSize(0.05);
    title1_2->SetTextAlign(22);  
    title1_2->DrawLatex(0.5, 0.97, "Q^{2} Resolution");

    c1->cd(3);  
    gPad->SetTopMargin(0.08);  
    TLatex* title1_3 = new TLatex();
    title1_3->SetNDC(); 
    title1_3->SetTextSize(0.05);
    title1_3->SetTextAlign(22);  
    title1_3->DrawLatex(0.5, 0.97, "Q^{2} Response");

    c1->Print("./figures/plot_Q2.pdf");

// y
    TCanvas* c2 = new TCanvas("c2","c2",1,1,1200,800);
    c2->Divide(3,1,0.01,0.01);
    c2->cd(1);
    gPad->SetLogy(1);
	h_y_e->GetYaxis()->SetTitle("counts");
	h_y_e->GetXaxis()->SetTitle("y");
	h_y_e->SetLineColor(kBlack);
    h_y_e->Draw();
    //h_y_beforeCut->SetMarkerStyle(11);
    //h_y_beforeCut->SetMarkerColor(kOrange);
    //h_y_beforeCut->Draw("Psame");
    //h_y_afterCut->SetMarkerStyle(29);
    //h_y_afterCut->SetMarkerColor(kRed);
    //h_y_afterCut->Draw("Psame");
    h_yREC_e_EEMC_cut->SetMarkerColor(kBlue);
    h_yREC_e_EEMC_cut->SetMarkerStyle(24);
    h_yREC_e_EEMC_cut->Draw("PEsame");
    //h_yREC_e_EEMC->SetMarkerStyle(21);
    //h_yREC_e_EEMC->SetMarkerColor(kCyan);
    //h_yREC_e_EEMC->Draw("PEsame");
    TLegend *w2 = new TLegend(0.45,0.7,0.6,0.85);
	w2->AddEntry(h_y_e, "MC", "L");
    //w2->AddEntry(h_y_beforeCut, "MC (before cut)", "P");
	w2->AddEntry(h_yREC_e_EEMC_cut, "RECO", "P");
    //w2->AddEntry(h_y_afterCut, "reco (after cut)", "P");
    //w2->AddEntry(h_yREC_e_EEMC, "reco (before cut)", "P");
	w2->Draw("same");

    /*c2->cd(2);
    gPad->SetLogz(1);
	h_y_res->GetYaxis()->SetTitle("(y_{MC}-y_{RECO})/y_{MC}");
	h_y_res->GetXaxis()->SetTitle("y_{MC}");
    h_y_res->Draw("colzsame");*/

    c2->cd(2);
    gPad->SetLogz(1);
    h_y_res_cut->GetYaxis()->SetTitleOffset(1.5);
	h_y_res_cut->GetYaxis()->SetTitle("(y_{MC}-y_{RECO})/y_{MC}");
	h_y_res_cut->GetXaxis()->SetTitle("y_{MC}");
    h_y_res_cut->Draw("colzsame");

    /*c2->cd(4);
    gPad->SetLogz(1);
    h_y_response->GetYaxis()->SetTitle("y_{MC}");
    h_y_response->GetXaxis()->SetTitle("y_{RECO}");
    h_y_response->Draw("colzsame)");*/

    c2->cd(3);
    gPad->SetLogz(1);
    h_y_response_cut->GetYaxis()->SetTitle("y_{MC}");
    h_y_response_cut->GetXaxis()->SetTitle("y_{RECO}");
    h_y_response_cut->Draw("colzsame)");

    /*c2->cd(6);
    gPad->SetLogy(1);
    h_dyOvery_REC->GetYaxis()->SetTitle("counts");
    h_dyOvery_REC->GetXaxis()->SetTitle("dy/y_{MC}");
    h_dyOvery_REC->Draw("same");*/

    /*c2->cd(7);
    gPad->SetLogy(1);
    h_dyOvery_REC_cut->GetYaxis()->SetTitle("counts");
    h_dyOvery_REC_cut->GetXaxis()->SetTitle("dy/y_{MC}");
    h_dyOvery_REC_cut->Draw("same");*/

    c2->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title2 = new TLatex();
    title2->SetNDC(); 
    title2->SetTextSize(0.05);
    title2->SetTextAlign(22);  
    title2->DrawLatex(0.5, 0.97, "y Truth vs Reco");  

    c2->cd(2);  
    gPad->SetTopMargin(0.08);  
    TLatex* title2_2 = new TLatex();
    title2_2->SetNDC(); 
    title2_2->SetTextSize(0.05);
    title2_2->SetTextAlign(22);  
    title2_2->DrawLatex(0.5, 0.97, "y Resolution");

    c2->cd(3);  
    gPad->SetTopMargin(0.08);  
    TLatex* title2_3 = new TLatex();
    title2_3->SetNDC(); 
    title2_3->SetTextSize(0.05);
    title2_3->SetTextAlign(22);  
    title2_3->DrawLatex(0.5, 0.97, "y Response");
	
    c2->Print("./figures/plot_y.pdf");


// Energy
    TCanvas* c3 = new TCanvas("c3","c3",1,1,1600,800);
    c3->Divide(3,1,0.01,0.01);
    c3->cd(1);
    gPad->SetLogy(1);
	h_energy_MC->GetYaxis()->SetTitle("counts");
	h_energy_MC->GetXaxis()->SetTitle("E [GeV]");
	h_energy_MC->SetLineColor(kBlack);
    h_energy_MC->Draw();
    h_energy_REC_EEMC_cut->SetMarkerColor(kBlue);
    h_energy_REC_EEMC_cut->SetMarkerStyle(24);
    h_energy_REC_EEMC_cut->Draw("PEsame");
    //h_energy_REC_EEMC->SetMarkerColor(kRed);
    //h_energy_REC_EEMC->SetMarkerStyle(29);
    //h_energy_REC_EEMC->Draw("PEsame");
    //h_energy_MC_after->SetMarkerStyle(11);
    //h_energy_MC_after->SetMarkerColor(kCyan);
    //h_energy_MC_after->Draw("PEsame");
    TLegend *w3 = new TLegend(0.22,0.2,0.4,0.3);
	w3->AddEntry(h_energy_MC, "MC (before cut)", "L");
	//w3->AddEntry(h_energy_REC_EEMC, "reco (before cut)", "P");
    w3->AddEntry(h_energy_REC_EEMC_cut, "reco (inside threshold)", "P");
    //w3->AddEntry(h_energy_MC_after, "MC (after cut)", "P");
	w3->Draw("same");

   /* c3->cd(2);
    gPad->SetLogz(1);
    h_energy_res_EEMC->GetXaxis()->SetRangeUser(0,10);
    h_energy_res_EEMC->GetYaxis()->SetTitle("(E_{MC}-E_{RECO})/E_{MC}");
    h_energy_res_EEMC->GetXaxis()->SetTitle("E_{MC} [GeV]");
    h_energy_res_EEMC->Draw("same");*/

    c3->cd(2);
    gPad->SetLogz(1);
    h_energy_res_EEMC_cut->GetXaxis()->SetRangeUser(0,10);
    h_energy_res_EEMC_cut->GetYaxis()->SetTitle("(E_{MC}-E_{RECO})/E_{MC}");
    h_energy_res_EEMC_cut->GetXaxis()->SetTitle("E_{MC} [GeV]");
    h_energy_res_EEMC_cut->Draw("same");

    /*c3->cd(4);
    gPad->SetLogz(1);
    h_energy_response_EEMC->GetXaxis()->SetRangeUser(0,10);
    h_energy_response_EEMC->GetYaxis()->SetRangeUser(0,10);
    h_energy_response_EEMC->GetYaxis()->SetTitle("E_{MC} [GeV]");
    h_energy_response_EEMC->GetXaxis()->SetTitle("E_{RECO} [GeV]");
    h_energy_response_EEMC->Draw("same");*/

    c3->cd(3);
    gPad->SetLogz(1);
    h_energy_response_EEMC_cut->GetXaxis()->SetRangeUser(0,10);
    h_energy_response_EEMC_cut->GetYaxis()->SetRangeUser(0,10);
    h_energy_response_EEMC_cut->GetYaxis()->SetTitle("E_{MC} [GeV]");
    h_energy_response_EEMC_cut->GetXaxis()->SetTitle("E_{RECO} [GeV]");
    h_energy_response_EEMC_cut->Draw("same");

    c3->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title3 = new TLatex();
    title3->SetNDC(); 
    title3->SetTextSize(0.05);
    title3->SetTextAlign(22);  
    title3->DrawLatex(0.5, 0.97, "e' Energy Truth vs. Reco");  

    c3->cd(2);  
    gPad->SetTopMargin(0.08);  
    TLatex* title3_2 = new TLatex();
    title3_2->SetNDC(); 
    title3_2->SetTextSize(0.05);
    title3_2->SetTextAlign(22);  
    title3_2->DrawLatex(0.5, 0.97, "e' Energy Resolution");

    c3->cd(3);  
    gPad->SetTopMargin(0.08);  
    TLatex* title3_3 = new TLatex();
    title3_3->SetNDC(); 
    title3_3->SetTextSize(0.05);
    title3_3->SetTextAlign(22);  
    title3_3->DrawLatex(0.5, 0.97, "e' Energy Response");  

    c3->Print("./figures/plot_energy.pdf");


// e' p
    TCanvas* c32 = new TCanvas("c32","c32",1,1,1200,800);
    c32->Divide(6,2,0.01,0.01);
    c32->cd(1);
    gPad->SetLogy(1);
	h_e_pt_MC->GetYaxis()->SetTitle("counts");
	h_e_pt_MC->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	h_e_pt_MC->SetLineColor(kBlack);
    h_e_pt_MC->Draw();
    h_e_pt_MC_after->SetMarkerColor(kBlue);
    h_e_pt_MC_after->SetMarkerStyle(24);
    h_e_pt_MC_after->Draw("PEsame");
    h_e_pt_REC_EEMC_cut->SetMarkerColor(kRed);
    h_e_pt_REC_EEMC_cut->SetMarkerStyle(29);
    h_e_pt_REC_EEMC_cut->Draw("PEsame");
    //h_e_pt_REC_EEMC->SetMarkerStyle(11);
    //h_e_pt_REC_EEMC->SetMarkerColor(kCyan);
    //h_e_pt_REC_EEMC->Draw("PEsame");
    //h_e_pt_REC_trk->SetMarkerStyle(21);
    //h_e_pt_REC_trk->SetMarkerColor(kOrange);
    //h_e_pt_REC_trk->Draw("PEsame");
    h_e_pt_REC_trk_cut->SetMarkerStyle(13);
    h_e_pt_REC_trk_cut->SetMarkerColor(kOrange);
    h_e_pt_REC_trk_cut->Draw("PEsame");
    TLegend *w32 = new TLegend(0.22,0.2,0.4,0.3);
	w32->AddEntry(h_e_pt_MC, "MC", "L");
	w32->AddEntry(h_e_pt_MC_after, "MC (after cut)", "P");
    w32->AddEntry(h_e_pt_REC_EEMC_cut, "reco (inside threshold)", "P");
    //w32->AddEntry(h_e_pt_REC_EEMC, "reco (after cut)", "P");
    //w32->AddEntry(h_e_pt_REC_trk, "reco (after cut)", "P");
    w32->AddEntry(h_e_pt_REC_trk_cut, "RECO", "P");
	w32->Draw("same");

    c32->cd(2);
    gPad->SetLogz(1);
    h_e_pt_res->Draw("colzsame");

    c32->cd(3);
    gPad->SetLogz(1);
    h_e_pt_res_cut->Draw("colzsame");

    c32->cd(4);
    gPad->SetLogy(1);
	h_e_pz_MC->GetYaxis()->SetTitle("counts");
	h_e_pz_MC->GetXaxis()->SetTitle("p_{z} [GeV/c]");
	h_e_pz_MC->SetLineColor(kBlack);
    h_e_pz_MC->Draw();
    //h_e_pz_MC_after->SetMarkerColor(kBlue);
    //h_e_pz_MC_after->SetMarkerStyle(24);
    //h_e_pz_MC_after->Draw("PEsame");
    //h_e_pz_REC_EEMC_cut->SetMarkerColor(kRed);
    //h_e_pz_REC_EEMC_cut->SetMarkerStyle(29);
    //h_e_pz_REC_EEMC_cut->Draw("PEsame");
    //h_e_pz_REC_EEMC->SetMarkerStyle(11);
    //h_e_pz_REC_EEMC->SetMarkerColor(kCyan);
    //h_e_pz_REC_EEMC->Draw("PEsame");
    //h_e_pz_REC_trk->SetMarkerStyle(21);
    //h_e_pz_REC_trk->SetMarkerColor(kOrange);
    //h_e_pz_REC_trk->Draw("PEsame");
    h_e_pz_REC_trk_cut->SetMarkerStyle(13);
    h_e_pz_REC_trk_cut->SetMarkerColor(kBlack);
    h_e_pz_REC_trk_cut->Draw("same");
    w32->Draw("same");

    c32->cd(5);
    gPad->SetLogz(1);
    h_e_pz_res->Draw("colzsame");

    c32->cd(6);
    gPad->SetLogz(1);
    h_e_pz_res_cut->Draw("colzsame");

    c32->cd(7);
    gPad->SetLogz(1);
    h_e_pz_response->Draw("colzsame");

    c32->cd(8);
    gPad->SetLogz(1);
    h_e_pz_response_cut->Draw("colzsame");

    c32->cd(9);
    gPad->SetLogy(1);
	h_e_p_MC->GetYaxis()->SetTitle("counts");
	h_e_p_MC->GetXaxis()->SetTitle("|p| [GeV/c]");
	h_e_p_MC->SetLineColor(kBlack);
    h_e_p_MC->Draw();
    //h_e_p_MC_after->SetMarkerColor(kBlue);
    //h_e_p_MC_after->SetMarkerStyle(24);
    //h_e_p_MC_after->Draw("PEsame");
    //h_e_p_REC_EEMC_cut->SetMarkerColor(kRed);
    //h_e_p_REC_EEMC_cut->SetMarkerStyle(29);
    //h_e_p_REC_EEMC_cut->Draw("PEsame");
    //h_e_p_REC_EEMC->SetMarkerStyle(11);
    //h_e_p_REC_EEMC->SetMarkerColor(kCyan);
    //h_e_p_REC_EEMC->Draw("PEsame");
    //h_e_p_REC_trk->SetMarkerStyle(21);
    //h_e_p_REC_trk->SetMarkerColor(kOrange);
    //h_e_p_REC_trk->Draw("PEsame");
    h_e_p_REC_trk_cut->SetMarkerStyle(13);
    h_e_p_REC_trk_cut->SetMarkerColor(kBlack);
    h_e_p_REC_trk_cut->Draw("same");
    w32->Draw("same");

    c32->cd(10);
    gPad->SetLogz(1);
    h_e_p_res->Draw("colzsame");

    c32->cd(11);
    gPad->SetLogz(1);
    h_e_p_res_cut->Draw("colzsame");
   
    c32->Print("./figures/plot_eScat_p.pdf");


// E and p
    TCanvas* c31 = new TCanvas("c31","c31",1,1,1200,800);
    c31->Divide(2,1,0.01,0.01);
    c31->cd(1);
    gPad->SetLogz(1);
    h_EvsP_MC->GetYaxis()->SetRangeUser(0,11);
    h_EvsP_MC->GetXaxis()->SetRangeUser(0,11);
    h_EvsP_MC->GetXaxis()->SetTitle("|p|_{MC} [GeV]");
    h_EvsP_MC->GetYaxis()->SetTitle("E_{MC} [GeV]");
    h_EvsP_MC->Draw("same");

    /*c31->cd(2);
    gPad->SetLogz(1);
    h_EvsP_REC->GetXaxis()->SetTitle("|p|_{trk} [GeV]");
    h_EvsP_REC->GetYaxis()->SetTitle("E_{EMCal} [GeV]");
    h_EvsP_REC->Draw("same");*/

    c31->cd(2);
    gPad->SetLogz(1);
    h_EvsP_REC_cut->GetXaxis()->SetTitle("|p|_{trk} [GeV]");
    h_EvsP_REC_cut->GetYaxis()->SetTitle("E_{EMCal} [GeV]");
    h_EvsP_REC_cut->Draw("same");

    c31->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title31 = new TLatex();
    title31->SetNDC(); 
    title31->SetTextSize(0.05);
    title31->SetTextAlign(22);  
    title31->DrawLatex(0.5, 0.97, "E vs. P (MC)");  

    c31->cd(2);  
    gPad->SetTopMargin(0.08);  
    TLatex* title31_2 = new TLatex();
    title31_2->SetNDC(); 
    title31_2->SetTextSize(0.05);
    title31_2->SetTextAlign(22);  
    title31_2->DrawLatex(0.5, 0.97, "E vs. P (RECO)");

    c31->Print("./figures/plot_EvsP.pdf");


// eta
    TCanvas* c4 = new TCanvas("c4","c4",1,1,1200,800);
    c4->Divide(4,1,0.01,0.01);
    c4->cd(1);
    gPad->SetLogy(1);
	h_eta_MC->GetYaxis()->SetTitle("counts");
	h_eta_MC->GetXaxis()->SetTitle("#eta");
	h_eta_MC->SetLineColor(kBlack);
    h_eta_MC->Draw();
    h_eta_REC_EEMC->SetMarkerColor(kBlue);
    h_eta_REC_EEMC->SetMarkerStyle(24);
    h_eta_REC_EEMC->Draw("PEsame");
    //h_eta_MC_after->SetMarkerStyle(11);
    //h_eta_MC_after->SetMarkerColor(kOrange);
    //h_eta_MC_after->Draw("PEsame");
    h_eta_REC_EEMC_cut->SetMarkerStyle(29);
    h_eta_REC_EEMC_cut->SetMarkerColor(kRed);
    h_eta_REC_EEMC_cut->Draw("PEsame");
    h_eta_REC_EEMC_after->SetMarkerStyle(31);
    h_eta_REC_EEMC_after->SetMarkerColor(kCyan);
    h_eta_REC_EEMC_after->Draw("PEsame");
    TLegend *w4 = new TLegend(0.22,0.2,0.4,0.35);
	w4->AddEntry(h_eta_MC, "MC (before cut)", "L");
	w4->AddEntry(h_eta_REC_EEMC, "reco (before cut)", "P");
    //w4->AddEntry(h_eta_MC_after, "MC (after cut)", "P");
    w4->AddEntry(h_eta_REC_EEMC_cut, "reco (inside threshold)", "P");
    w4->AddEntry(h_eta_REC_EEMC_after, "reco (after cut)", "P");
	w4->Draw("same");

    c4->cd(2);
    gPad->SetLogy(1);
	h_eta_diff->GetYaxis()->SetTitle("counts");
	h_eta_diff->GetXaxis()->SetTitle("#eta_{RECO}-#eta_{MC}");
	h_eta_diff->SetLineColor(kBlack);
    h_eta_diff->Draw();

    c4->cd(3);
    gPad->SetLogy(1);
	h_eta_diff_cut->GetYaxis()->SetTitle("counts");
	h_eta_diff_cut->GetXaxis()->SetTitle("#eta_{RECO}-#eta_{MC}");
	h_eta_diff_cut->SetLineColor(kBlack);
    h_eta_diff_cut->Draw();

    c4->cd(4);
    gPad->SetLogy(1);
	h_eta_diff_after->GetYaxis()->SetTitle("counts");
	h_eta_diff_after->GetXaxis()->SetTitle("#eta_{RECO}-#eta_{MC}");
	h_eta_diff_after->SetLineColor(kBlack);
    h_eta_diff_after->Draw();
       
    c4->Print("./figures/plot_eta.pdf");


// phi
    TCanvas* c5 = new TCanvas("c5","c5",1,1,1200,800);
    c5->Divide(2,2,0.01,0.01);
    c5->cd(1);
    gPad->SetLogy(1);
	h_phi_MC->GetYaxis()->SetTitle("counts");
	h_phi_MC->GetXaxis()->SetTitle("#phi");
	h_phi_MC->SetLineColor(kBlack);
    h_phi_MC->Draw();
    h_phi_REC_EEMC->SetMarkerColor(kBlue);
    h_phi_REC_EEMC->SetMarkerStyle(24);
    h_phi_REC_EEMC->Draw("PEsame");
    h_phi_MC_after->SetMarkerStyle(11);
    h_phi_MC_after->SetMarkerColor(kCyan);
    h_phi_MC_after->Draw("PEsame");
    h_phi_REC_EEMC_cut->SetMarkerStyle(29);
    h_phi_REC_EEMC_cut->SetMarkerColor(kRed);
    h_phi_REC_EEMC_cut->Draw("PEsame");
    TLegend *w5 = new TLegend(0.22,0.7,0.4,0.85);
    h_phi_REC_EEMC_after->SetMarkerColor(kOrange);
    h_phi_REC_EEMC_after->SetMarkerStyle(31);
    h_phi_REC_EEMC_after->Draw("PEsame");
	w5->AddEntry(h_phi_MC, "MC (before cut)", "L");
	w5->AddEntry(h_phi_REC_EEMC, "reco (before cut)", "P");
    w5->AddEntry(h_phi_MC_after, "MC (after cut)", "P");
    w5->AddEntry(h_phi_REC_EEMC_cut, "reco (inside threshold)", "P");
    w5->AddEntry(h_phi_REC_EEMC_after, "reco (after cut)", "P");
	w5->Draw("same");

    c5->cd(2);
    gPad->SetLogy(1);
	h_phi_diff->GetYaxis()->SetTitle("counts");
	h_phi_diff->GetXaxis()->SetTitle("#phi_{RECO}-#phi_{MC}");
	h_phi_diff->SetLineColor(kBlack);
    h_phi_diff->Draw();

    c5->cd(3);
    gPad->SetLogy(1);
	h_phi_diff_cut->GetYaxis()->SetTitle("counts");
	h_phi_diff_cut->GetXaxis()->SetTitle("#phi_{RECO}-#phi_{MC}");
	h_phi_diff_cut->SetLineColor(kBlack);
    h_phi_diff_cut->Draw();

    c5->cd(4);
    gPad->SetLogy(1);
	h_phi_diff_after->GetYaxis()->SetTitle("counts");
	h_phi_diff_after->GetXaxis()->SetTitle("#phi_{RECO}-#phi_{MC}");
	h_phi_diff_after->SetLineColor(kBlack);
    h_phi_diff_after->Draw();
   
    c5->Print("./figures/plot_phi.pdf");


// theta
    TCanvas* c6 = new TCanvas("c6","c6",1,1,1200,800);
    c6->Divide(2,2,0.01,0.01);
    c6->cd(1);
    gPad->SetLogy(1);
	//h_theta_MC->GetYaxis()->SetTitle("counts");
	//h_theta_MC->GetXaxis()->SetTitle("#phi");
	//h_theta_MC->SetLineColor(kBlack);
    //h_theta_MC->Draw();
    h_theta_MC_after->GetXaxis()->SetRangeUser(-1,3.5);
    h_theta_MC_after->GetYaxis()->SetTitle("counts");
	h_theta_MC_after->GetXaxis()->SetTitle("#theta");
    h_theta_MC_after->SetLineColor(kBlack);
    //h_theta_MC_after->SetMarkerStyle(11);
    h_theta_MC_after->Draw();
    h_theta_REC_EEMC_cut->SetMarkerColor(kBlue);
    h_theta_REC_EEMC_cut->SetMarkerStyle(24);
    h_theta_REC_EEMC_cut->Draw("PEsame");
    //h_theta_REC_EEMC_after->SetMarkerColor(kOrange);
    //h_theta_REC_EEMC_after->SetMarkerStyle(29);
    //h_theta_REC_EEMC_after->Draw("PEsame");
    TLegend *w6 = new TLegend(0.22,0.7,0.4,0.85);
	//w6->AddEntry(h_theta_MC, "MC (before cut)", "L");
    w6->AddEntry(h_theta_MC_after, "MC", "L");
	w6->AddEntry(h_theta_REC_EEMC_cut, "RECO", "P");
    //w6->AddEntry(h_theta_REC_EEMC_after, "reco (after cut)", "P");
	w6->Draw("same");

    /*c6->cd(2);
    gPad->SetLogy(1);
	h_theta_diff_after->GetYaxis()->SetTitle("counts");
	h_theta_diff_after->GetXaxis()->SetTitle("#theta_{RECO}-#theta_{MC}");
	h_theta_diff_after->SetLineColor(kBlack);
    h_theta_diff_after->Draw();*/

    c6->cd(2);
    //gPad->SetLogy(1);
	h_theta_diff_cut->GetYaxis()->SetTitle("counts");
	h_theta_diff_cut->GetXaxis()->SetTitle("#theta_{RECO}-#theta_{MC}");
	h_theta_diff_cut->SetLineColor(kBlack);
    h_theta_diff_cut->Draw();

    /*c6->cd(4);
    gPad->SetLogz(1);
    h_theta_response_EEMC_after->GetYaxis()->SetTitle("#theta_{MC} [GeV]");
    h_theta_response_EEMC_after->GetXaxis()->SetTitle("#theta_{RECO} [GeV]");
    h_theta_response_EEMC_after->Draw("same");*/

    c6->cd(3);
    gPad->SetLogz(1);
    h_theta_response_EEMC_cut->GetXaxis()->SetRangeUser(-0.5,3.1);
    h_theta_response_EEMC_cut->GetYaxis()->SetRangeUser(2.3,3.14);
    h_theta_response_EEMC_cut->GetYaxis()->SetTitle("#theta_{MC}");
    h_theta_response_EEMC_cut->GetXaxis()->SetTitle("#theta_{RECO}");
    h_theta_response_EEMC_cut->Draw("same");

    c6->cd(4);
    gPad->SetLogz(1);
    h_theta_response_EEMC_cut_clone->GetXaxis()->SetRangeUser(2.5,3.1);
    h_theta_response_EEMC_cut_clone->GetYaxis()->SetRangeUser(2.3,3.14);
    h_theta_response_EEMC_cut_clone->GetYaxis()->SetTitle("#theta_{MC}");
    h_theta_response_EEMC_cut_clone->GetXaxis()->SetTitle("#theta_{RECO}");
    h_theta_response_EEMC_cut_clone->Draw("same");

    c6->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title6 = new TLatex();
    title6->SetNDC(); 
    title6->SetTextSize(0.05);
    title6->SetTextAlign(22);  
    title6->DrawLatex(0.5, 0.97, "#theta Truth vs Reco");  

    c6->cd(2);  
    gPad->SetTopMargin(0.08);  
    TLatex* title6_2 = new TLatex();
    title6_2->SetNDC(); 
    title6_2->SetTextSize(0.05);
    title6_2->SetTextAlign(22);  
    title6_2->DrawLatex(0.5, 0.97, "#theta Difference");

    c6->cd(3);  
    gPad->SetTopMargin(0.08);  
    TLatex* title6_3 = new TLatex();
    title6_3->SetNDC(); 
    title6_3->SetTextSize(0.05);
    title6_3->SetTextAlign(22);  
    title6_3->DrawLatex(0.5, 0.97, "#theta Response");  
	
    c6->Print("./figures/plot_theta.pdf");


// VM (E-pz)
    TCanvas* c7 = new TCanvas("c7","c7",1,1,1200,800);
    c7->Divide(2,1,0.01,0.01);
    c7->cd(1);
    gPad->SetLogy(1);
	h_VM_Epz_MC->GetYaxis()->SetTitle("counts");
	h_VM_Epz_MC->GetXaxis()->SetTitle("(E_{VM}-p_{z,VM}) [GeV]");
	h_VM_Epz_MC->SetLineColor(kBlack);
    h_VM_Epz_MC->Draw();
    //h_VM_Epz_REC->SetMarkerStyle(11);
    //h_VM_Epz_REC->SetMarkerColor(kRed);
    //h_VM_Epz_REC->Draw("PEsame");
    h_VM_Epz_REC_cut->SetMarkerColor(kBlue);
    h_VM_Epz_REC_cut->SetMarkerStyle(24);
    h_VM_Epz_REC_cut->Draw("PEsame");
    TLegend *w7 = new TLegend(0.42,0.75,0.6,0.85);
	w7->AddEntry(h_VM_Epz_MC, "MC", "L");
	w7->AddEntry(h_VM_Epz_REC_cut, "RECO", "P");
    //w7->AddEntry(h_VM_Epz_REC, "reco ", "P");
	w7->Draw("same");

    /*c7->cd(2);
    gPad->SetLogz(1);
	h_VM_Epz_response->GetYaxis()->SetTitle("(E_{VM,MC}-p_{z,VM,MC}) [GeV]");
	h_VM_Epz_response->GetXaxis()->SetTitle("(E_{VM,RECO}-p_{z,VM,RECO}) [GeV]");
    h_VM_Epz_response->Draw("colzsame");*/

    c7->cd(2);
    gPad->SetLogz(1);
	h_VM_Epz_response_cut->GetYaxis()->SetTitle("(E_{VM,MC}-p_{z,VM,MC}) [GeV]");
	h_VM_Epz_response_cut->GetXaxis()->SetTitle("(E_{VM,RECO}-p_{z,VM,RECO}) [GeV]");
    h_VM_Epz_response_cut->Draw("colzsame");

    c7->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title7 = new TLatex();
    title7->SetNDC(); 
    title7->SetTextSize(0.05);
    title7->SetTextAlign(22);  
    title7->DrawLatex(0.5, 0.97, "E-p_{z} Truth vs Reco (VM)");  

    c7->cd(2);  
    gPad->SetTopMargin(0.08);  
    TLatex* title7_2 = new TLatex();
    title7_2->SetNDC(); 
    title7_2->SetTextSize(0.05);
    title7_2->SetTextAlign(22);  
    title7_2->DrawLatex(0.5, 0.97, "E-p_{z} Response (VM)");

    c7->Print("./figures/plot_VM_Epz.pdf");


// VM pT
    TCanvas* c8 = new TCanvas("c8","c8",1,1,1200,800);
    c8->Divide(2,1,0.01,0.01);
    c8->cd(1);
    gPad->SetLogy(1);
	h_VM_pt_MC->GetYaxis()->SetTitle("counts");
	h_VM_pt_MC->GetXaxis()->SetTitle("p_{T,VM} [GeV/c]");
	h_VM_pt_MC->SetLineColor(kBlack);
    h_VM_pt_MC->Draw();
    //h_VM_pt_REC->SetMarkerColor(kRed);
    //h_VM_pt_REC->SetMarkerStyle(11);
    //h_VM_pt_REC->Draw("PEsame");
    h_VM_pt_REC_cut->SetMarkerColor(kBlue);
    h_VM_pt_REC_cut->SetMarkerStyle(24);
    h_VM_pt_REC_cut->Draw("PEsame");
    TLegend *w8 = new TLegend(0.45,0.2,0.6,0.35);
	w8->AddEntry(h_VM_pt_MC, "MC", "L");
	w8->AddEntry(h_VM_pt_REC_cut, "RECO", "P");
    //w8->AddEntry(h_VM_pt_REC, "reco", "P");
	w8->Draw("same");

    /*c8->cd(2);
    gPad->SetLogz(1);
	h_VM_pt_response->GetYaxis()->SetTitle("p_{T,VM,MC} [GeV/c]");
	h_VM_pt_response->GetXaxis()->SetTitle("p_{T,VM,RECO} [GeV/c]");
    h_VM_pt_response->Draw("colzsame");*/

    c8->cd(2);
    gPad->SetLogz(1);
	h_VM_pt_response_cut->GetYaxis()->SetTitle("p_{T,VM,MC} [GeV/c]");
	h_VM_pt_response_cut->GetXaxis()->SetTitle("p_{T,VM,RECO} [GeV/c]");
    h_VM_pt_response_cut->Draw("colzsame");

    c8->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title8 = new TLatex();
    title8->SetNDC(); 
    title8->SetTextSize(0.05);
    title8->SetTextAlign(22);  
    title8->DrawLatex(0.5, 0.97, "p_{T} Truth vs Reco (VM)");  

    c8->cd(2);  
    gPad->SetTopMargin(0.08);  
    TLatex* title8_2 = new TLatex();
    title8_2->SetNDC(); 
    title8_2->SetTextSize(0.05);
    title8_2->SetTextAlign(22);  
    title8_2->DrawLatex(0.5, 0.97, "p_{T} Response (VM)");

    c8->Print("./figures/plot_VM_pt.pdf");


// VM 
    TCanvas* c9 = new TCanvas("c9","c9",1,1,1200,800);
    c9->Divide(3,1,0.01,0.01);
    c9->cd(1);
    gPad->SetLogy(1);
	h_VM_p_MC->GetYaxis()->SetTitle("counts");
	h_VM_p_MC->GetXaxis()->SetTitle("|p|_{VM} [GeV/c]");
	h_VM_p_MC->SetLineColor(kBlack);
    h_VM_p_MC->Draw();
    h_VM_p_REC->SetMarkerColor(kRed);
    h_VM_p_REC->SetMarkerStyle(11);
    h_VM_p_REC->Draw("PEsame");
    h_VM_p_REC_cut->SetMarkerColor(kBlue);
    h_VM_p_REC_cut->SetMarkerStyle(24);
    h_VM_p_REC_cut->Draw("PEsame");
    TLegend *w9 = new TLegend(0.22,0.7,0.4,0.85);
	w9->AddEntry(h_VM_p_MC, "MC", "L");
    w9->AddEntry(h_VM_p_REC, "reco", "P");
	w9->AddEntry(h_VM_p_REC_cut, "reco (inside threshold)", "P");
	w9->Draw("same");

    c9->cd(2);
    gPad->SetLogy(1);
	h_VM_pz_MC->GetYaxis()->SetTitle("counts");
	h_VM_pz_MC->GetXaxis()->SetTitle("p_{z,VM} [GeV/c]");
	h_VM_pz_MC->SetLineColor(kBlack);
    h_VM_pz_MC->Draw();
    h_VM_pz_REC->SetMarkerColor(kRed);
    h_VM_pz_REC->SetMarkerStyle(11);
    h_VM_pz_REC->Draw("PEsame");
    h_VM_pz_REC_cut->SetMarkerColor(kBlue);
    h_VM_pz_REC_cut->SetMarkerStyle(24);
    h_VM_pz_REC_cut->Draw("PEsame");
    w9->Draw("same");

    c9->cd(3);
    gPad->SetLogy(1);
	h_VM_mass_MC->GetYaxis()->SetTitle("counts");
	h_VM_mass_MC->GetXaxis()->SetTitle("VM mass [GeV/c^{2}]");
	h_VM_mass_MC->SetLineColor(kBlack);
    h_VM_mass_MC->Draw();
    h_VM_mass_REC->SetLineColor(kRed);
    h_VM_mass_REC->SetLineStyle(9);
    h_VM_mass_REC->Draw("same");
    h_VM_mass_REC->SetMarkerStyle(11);
    h_VM_mass_REC->Draw("PEsame");
    h_VM_mass_REC_cut->SetLineColor(kBlue);
    h_VM_mass_REC_cut->SetLineStyle(3);
    h_VM_mass_REC_cut->Draw("same");
    h_VM_mass_REC_cut->SetMarkerStyle(24);
    h_VM_mass_REC_cut->Draw("PEsame");
    w9->Draw("same");

    c9->Print("./figures/plot_VM.pdf");


// position
    TCanvas* c10 = new TCanvas("c10","c10",1,1,1600,800);
    c10->Divide(3,1,0.01,0.01);
    c10->cd(1);
    gPad->SetLogy(1);
    h_Xclus_minus_Xtrk->GetYaxis()->SetTitle("counts");
	h_Xclus_minus_Xtrk->GetXaxis()->SetTitle("x_{clus}-x_{trk} [mm]");
    h_Xclus_minus_Xtrk->SetLineColor(kBlack);
    h_Xclus_minus_Xtrk->Draw();
    h_Xclus_minus_Xtrk_cut->SetMarkerColor(kBlue);
    h_Xclus_minus_Xtrk_cut->SetMarkerStyle(24);
    h_Xclus_minus_Xtrk_cut->Draw("PEsame");
    h_Xclus_minus_Xtrk_cut2->SetMarkerColor(kRed);
    h_Xclus_minus_Xtrk_cut2->SetMarkerStyle(29);
    h_Xclus_minus_Xtrk_cut2->Draw("PEsame");
    TLegend *w99 = new TLegend(0.22,0.7,0.4,0.85);
    w99->AddEntry(h_Xclus_minus_Xtrk, "RECO", "L");
	w99->AddEntry(h_Xclus_minus_Xtrk_cut, "Accept", "P");
    w99->AddEntry(h_Xclus_minus_Xtrk_cut2, "Reject", "P");
	w99->Draw("same");
   
    c10->cd(2);
    gPad->SetLogy(1);
    h_Yclus_minus_Ytrk->GetYaxis()->SetTitle("counts");
	h_Yclus_minus_Ytrk->GetXaxis()->SetTitle("y_{clus}-y_{trk} [mm]");
    h_Yclus_minus_Ytrk->SetLineColor(kBlack);
    h_Yclus_minus_Ytrk->Draw();
	h_Yclus_minus_Ytrk_cut->SetMarkerColor(kBlue);
    h_Yclus_minus_Ytrk_cut->SetMarkerStyle(24);
    h_Yclus_minus_Ytrk_cut->Draw("PEsame");
    h_Yclus_minus_Ytrk_cut2->SetMarkerColor(kRed);
    h_Yclus_minus_Ytrk_cut2->SetMarkerStyle(29);
    h_Yclus_minus_Ytrk_cut2->Draw("PEsame");
    w99->Draw("same");

    c10->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title10 = new TLatex();
    title10->SetNDC(); 
    title10->SetTextSize(0.05);
    title10->SetTextAlign(22);  
    title10->DrawLatex(0.5, 0.97, "x Position Difference");  

    c10->cd(2);  
    gPad->SetTopMargin(0.08); 
    TLatex* title10_2 = new TLatex();
    title10_2->SetNDC(); 
    title10_2->SetTextSize(0.05);
    title10_2->SetTextAlign(22);  
    title10_2->DrawLatex(0.5, 0.97, "y Position Difference"); 

    c10->Print("./figures/plot_position.pdf");


    TCanvas* c10_1 = new TCanvas("c10_1","c10_1",1,1,1200,800);
    c10_1->Divide(2,1,0.01,0.01);
    c10_1->cd(1);
    gPad->SetLogz(1);
    h_emClus_position_REC->GetYaxis()->SetTitleOffset(1.5);
    h_emClus_position_REC->GetXaxis()->SetTitle("x [mm]");
    h_emClus_position_REC->GetYaxis()->SetTitle("y [mm]");
    h_emClus_position_REC->Draw("colzsame");

    c10_1->cd(2);
    gPad->SetLogz(1);
    h_emClus_position_REC_cut->GetYaxis()->SetTitleOffset(1.5);
    h_emClus_position_REC_cut->GetXaxis()->SetTitle("x [mm]");
    h_emClus_position_REC_cut->GetYaxis()->SetTitle("y [mm]");
    h_emClus_position_REC_cut->Draw("colzsame"); 

    c10_1->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title10_3 = new TLatex();
    title10_3->SetNDC(); 
    title10_3->SetTextSize(0.05);
    title10_3->SetTextAlign(22);  
    title10_3->DrawLatex(0.5, 0.97, "Cluster Positions");  

    c10_1->cd(2);  
    gPad->SetTopMargin(0.08); 
    TLatex* title10_4 = new TLatex();
    title10_4->SetNDC(); 
    title10_4->SetTextSize(0.05);
    title10_4->SetTextAlign(22);  
    title10_4->DrawLatex(0.5, 0.97, "Cluster Positions After Threshold Cut");  

    c10_1->Print("./figures/plot_position_2d.pdf");


// E-pz
    TCanvas* c11 = new TCanvas("c11","c11",1,1,1200,800);
    c11->Divide(2,2,0.01,0.01);
    c11->cd(1);
    gPad->SetLogy(1);
    h_Epz_MC->GetXaxis()->SetRangeUser(0,40);
	h_Epz_MC->GetYaxis()->SetTitle("counts");
	h_Epz_MC->GetXaxis()->SetTitle("(E - p_{z}) [GeV]");
	h_Epz_MC->SetLineColor(kBlack);
    h_Epz_MC->Draw();
    //h_Epz_REC->SetMarkerColor(kBlue);
    //h_Epz_REC->SetMarkerStyle(24);
    //h_Epz_REC->Draw("PEsame");
    h_Epz_REC_cut->SetMarkerStyle(29);
    h_Epz_REC_cut->SetMarkerColor(kRed);
    h_Epz_REC_cut->Draw("PEsame");
    //h_Epz_afterCut->SetMarkerStyle(24);
    //h_Epz_afterCut->SetMarkerColor(kBlue);
    //h_Epz_afterCut->Draw("PEsame");
    TLegend *w10 = new TLegend(0.22,0.7,0.4,0.85);
	w10->AddEntry(h_Epz_MC, "MC", "L");
	//w10->AddEntry(h_Epz_REC, "reco (after cut)", "P");
	w10->AddEntry(h_Epz_REC_cut, "RECO", "P");
    //w10->AddEntry(h_Epz_afterCut, "RECO", "P");
	w10->Draw("same");

  /*  c11->cd(2);
    gPad->SetLogz(1);
	h_Epz_response->GetYaxis()->SetTitle("(E_{MC}-p_{z,MC}) [GeV]");
	h_Epz_response->GetXaxis()->SetTitle("(E_{EMCal} - p_{z,trk}) [GeV]");
    h_Epz_response->Draw();*/

    c11->cd(2);
    gPad->SetLogz(1);
    h_Epz_response_cut->GetYaxis()->SetRangeUser(2,21);
	h_Epz_response_cut->GetYaxis()->SetTitle("(E_{MC}-p_{z,MC}) [GeV]");
	h_Epz_response_cut->GetXaxis()->SetTitle("(E_{EMCal} - p_{z,trk}) [GeV]");
    h_Epz_response_cut->Draw();

    c11->cd(3);
    gPad->SetLogy(1);
    h_EoverP_MC->GetXaxis()->SetRangeUser(0,2);
    h_EoverP_MC->SetTitle("E/|p| Truth vs. Reco");
	h_EoverP_MC->GetYaxis()->SetTitle("counts");
    h_EoverP_MC->GetXaxis()->SetTitle("E/|p|");
	h_EoverP_MC->GetXaxis()->SetRangeUser(0,6);
	h_EoverP_MC->SetLineColor(kBlack);
    h_EoverP_MC->Draw();
    //h_EcalOverPtrk->SetMarkerColor(kBlue);
    //h_EcalOverPtrk->SetMarkerStyle(24);
    //h_EcalOverPtrk->Draw("PEsame");
    h_EcalOverPtrk_cut->SetMarkerStyle(24);
    h_EcalOverPtrk_cut->SetMarkerColor(kBlue);
    h_EcalOverPtrk_cut->Draw("PEsame");
    //h_EoverP_MC_after->SetMarkerStyle(11);
    //h_EoverP_MC_after->SetMarkerColor(kOrange);
    //h_EoverP_MC_after->Draw("PEsame");
    TLegend *w111 = new TLegend(0.52,0.7,0.7,0.85);
	w111->AddEntry(h_EoverP_MC, "MC", "L");
	//w111->AddEntry(h_EcalOverPtrk, "reco (after cut)", "P");
    w111->AddEntry(h_EcalOverPtrk_cut, "RECO", "P");
    //w111->AddEntry(h_EoverP_MC_after, "MC (after cut)", "P");
	w111->Draw("same");

    c11->cd(4);
    gPad->SetLogz(1);
    h_EoverP_response_cut->SetTitle("E/|p| Response (After Threshold Cut)");   
    h_EoverP_response_cut->GetYaxis()->SetTitleOffset(1.5);
    h_EoverP_response_cut->GetYaxis()->SetRangeUser(1,1.02);
	h_EoverP_response_cut->GetYaxis()->SetTitle("E_{MC}/|p|_{MC}");
	h_EoverP_response_cut->GetXaxis()->SetTitle("E_{EMCal}/|p|_{trk}");
	h_EoverP_response_cut->SetLineColor(kBlack);
    h_EoverP_response_cut->Draw();

    c11->cd(4);  
    gPad->SetTopMargin(0.08);  
    TLatex* title12_22 = new TLatex();
    title12_22->SetNDC(); 
    title12_22->SetTextSize(0.05);
    title12_22->SetTextAlign(22);  
    title12_22->DrawLatex(0.5, 0.97, "E/|p| Response");  

    c11->cd(3);  
    gPad->SetTopMargin(0.08); 
    TLatex* title11_3 = new TLatex();
    title11_3->SetNDC(); 
    title11_3->SetTextSize(0.05);
    title11_3->SetTextAlign(22);  
    title11_3->DrawLatex(0.5, 0.97, "E/|p| Truth vs. Reco");  


    c11->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title11 = new TLatex();
    title11->SetNDC(); 
    title11->SetTextSize(0.05);
    title11->SetTextAlign(22);  
    title11->DrawLatex(0.5, 0.97, "E-p_{z} Truth vs. Reco (e'+HFS)");  

    c11->cd(2);  
    gPad->SetTopMargin(0.08);  
    TLatex* title11_2 = new TLatex();
    title11_2->SetNDC(); 
    title11_2->SetTextSize(0.05);
    title11_2->SetTextAlign(22);  
    title11_2->DrawLatex(0.5, 0.97, "E-p_{z} Response (e'+HFS)");
   
    c11->Print("./figures/plot_Epz.pdf");


// E/p
    TCanvas* c12 = new TCanvas("c12","c12",1,1,1600,800);
    c12->Divide(2,1,0.01,0.01);
    c12->cd(1);
    gPad->SetLogy(1);
    h_EoverP_MC->SetTitle("E/|p| Truth vs. Reco");
	h_EoverP_MC->GetYaxis()->SetTitle("counts");
	h_EoverP_MC->GetXaxis()->SetTitle("E/|p|");
	h_EoverP_MC->SetLineColor(kBlack);
    h_EoverP_MC->Draw();
    //h_EcalOverPtrk->SetMarkerColor(kCyan);
    //h_EcalOverPtrk->SetMarkerStyle(20);
    //h_EcalOverPtrk->Draw("PEsame");
    h_EcalOverPtrk_cut->SetMarkerStyle(24);
    h_EcalOverPtrk_cut->SetMarkerColor(kBlue);
    h_EcalOverPtrk_cut->Draw("PEsame");
    //h_EoverP_MC_after->SetMarkerStyle(11);
    //h_EoverP_MC_after->SetMarkerColor(kOrange);
    //h_EoverP_MC_after->Draw("PEsame");
    TLegend *w11 = new TLegend(0.22,0.7,0.4,0.85);
	w11->AddEntry(h_EoverP_MC, "MC", "L");
	//w11->AddEntry(h_EcalOverPtrk, "reco (after cut)", "P");
    w11->AddEntry(h_EcalOverPtrk_cut, "RECO", "P");
    //w11->AddEntry(h_EoverP_MC_after, "MC (after cut)", "P");
	w11->Draw("same");

    /*c12->cd(2);
    gPad->SetLogz(1);    
    h_EoverP_response->SetTitle("E/|p| Response");   
    h_EoverP_response->GetYaxis()->SetTitleOffset(1.5);
    h_EoverP_response->GetYaxis()->SetRangeUser(1, 1.2);
	h_EoverP_response->GetYaxis()->SetTitle("E_{MC}/|p|_{MC}");
	h_EoverP_response->GetXaxis()->SetTitle("E_{EMCal}/|p|_{trk}");
	h_EoverP_response->SetLineColor(kBlack);
    h_EoverP_response->Draw();*/

    c12->cd(2);
    gPad->SetLogz(1);
    h_EoverP_response_cut->SetTitle("E/|p| Response (After Threshold Cut)");   
    h_EoverP_response_cut->GetYaxis()->SetTitleOffset(1.5);
    h_EoverP_response_cut->GetYaxis()->SetRangeUser(1,1.02);
	h_EoverP_response_cut->GetYaxis()->SetTitle("E_{MC}/|p|_{MC}");
	h_EoverP_response_cut->GetXaxis()->SetTitle("E_{EMCal}/|p|_{trk}");
	h_EoverP_response_cut->SetLineColor(kBlack);
    h_EoverP_response_cut->Draw();

    c12->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title1 = new TLatex();
    title1->SetNDC(); 
    title1->SetTextSize(0.05);
    title1->SetTextAlign(22);  
    title1->DrawLatex(0.5, 0.97, "E/|p| Truth vs. Reco");  

    c12->cd(2);  
    gPad->SetTopMargin(0.08);  
    TLatex* title12_2 = new TLatex();
    title12_2->SetNDC(); 
    title12_2->SetTextSize(0.05);
    title12_2->SetTextAlign(22);  
    title12_2->DrawLatex(0.5, 0.97, "E/|p| Response");  

    c12->Print("./figures/plot_EoverP.pdf");



// t distribution
    TCanvas* c14 = new TCanvas("c14","c14",1,1,1800,800);
    c14->Divide(2,1,0.01,0.01);
    c14->cd(1);
    gPad->SetLogy(1);
	h_t_MC->GetYaxis()->SetTitle("counts");
    h_t_MC->GetXaxis()->SetLabelSize(0.03);
	h_t_MC->GetXaxis()->SetTitle("|t| [GeV/c]^{2}");
	h_t_MC->SetLineColor(kBlack);
    h_t_MC->Draw();
    h_t_REC_wRES_cut_pi12->SetLineStyle(2);
    h_t_REC_wRES_cut_pi12->SetLineColor(kCyan);
    h_t_REC_wRES_cut_pi12->Draw("same");
    h_t_REC_wRES_cut_pi12_cut->SetMarkerStyle(24);
    h_t_REC_wRES_cut_pi12_cut->SetMarkerColor(kOrange);
    h_t_REC_wRES_cut_pi12_cut->Draw("PEsame");
    h_t_REC_wRES_cut_pi12_cut2->SetMarkerStyle(29);
    h_t_REC_wRES_cut_pi12_cut2->SetMarkerColor(kBlue);
    h_t_REC_wRES_cut_pi12_cut2->Draw("PEsame");
    TLegend *w13 = new TLegend(0.38,0.7,0.65,0.85);
	w13->AddEntry(h_t_MC, "MC", "L");
	w13->AddEntry(h_t_REC_wRES_cut_pi12, "#pi/12", "L");
    w13->AddEntry(h_t_REC_wRES_cut_pi12_cut, "#pi/12 (inside)", "P");
    w13->AddEntry(h_t_REC_wRES_cut_pi12_cut2, "#pi/12 (outside)", "P");
	w13->Draw("same");

    c14->cd(2);
    gPad->SetLogy(1);
    h_t_MC->GetXaxis()->SetLabelSize(0.03);
	h_t_MC->GetYaxis()->SetTitle("counts");
	h_t_MC->GetXaxis()->SetTitle("|t| [GeV/c]^{2}");
	h_t_MC->SetLineColor(kBlack);
    h_t_MC->Draw();
    //h_t_REC_wRES_cut_pi9->SetLineColor(kBlue);
    //h_t_REC_wRES_cut_pi9->SetLineStyle(2);
    //h_t_REC_wRES_cut_pi9->Draw("same");
    h_t_REC_wRES_cut_pi9_cut->SetMarkerColor(kBlue);
    h_t_REC_wRES_cut_pi9_cut->SetMarkerStyle(24);
    h_t_REC_wRES_cut_pi9_cut->Draw("PEsame");
    //h_t_REC_wRES_cut_pi12->SetLineStyle(2);
    //h_t_REC_wRES_cut_pi12->SetLineColor(kRed);
    //h_t_REC_wRES_cut_pi12->Draw("same");
    h_t_REC_wRES_cut_pi12_cut->SetMarkerStyle(29);
    h_t_REC_wRES_cut_pi12_cut->SetMarkerColor(kRed);
    h_t_REC_wRES_cut_pi12_cut->Draw("PEsame");
    //h_t_REC_wRES_cut_pi16->SetLineStyle(2);
    //h_t_REC_wRES_cut_pi16->SetLineColor(kYellow);
    //h_t_REC_wRES_cut_pi16->Draw("same");
    h_t_REC_wRES_cut_pi16_cut->SetMarkerStyle(21);
    h_t_REC_wRES_cut_pi16_cut->SetMarkerColor(kOrange);
    h_t_REC_wRES_cut_pi16_cut->Draw("PEsame");
    TLegend *w14 = new TLegend(0.52,0.7,0.7,0.85);
	w14->AddEntry(h_t_MC, "MC", "L");
	//w14->AddEntry(h_t_REC_wRES_cut_pi9_cut, "#pi/9", "L");
    w14->AddEntry(h_t_REC_wRES_cut_pi9_cut, "#theta_{max} = #pi/9", "P");
	//w14->AddEntry(h_t_REC_wRES_cut_pi12, "#pi/12", "L");
    w14->AddEntry(h_t_REC_wRES_cut_pi12_cut, "#theta_{max} = #pi/12", "P");
	//w14->AddEntry(h_t_REC_wRES_cut_pi16, "#pi/16", "L");
    w14->AddEntry(h_t_REC_wRES_cut_pi16_cut, "#theta_{max} = #pi/16", "P");
	w14->Draw("same");

    c14->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title14 = new TLatex();
    title14->SetNDC(); 
    title14->SetTextSize(0.05);
    title14->SetTextAlign(22);  
    title14->DrawLatex(0.5, 0.97, "|t| Distribution Cut on Postion");  

    c14->cd(2);  
    gPad->SetTopMargin(0.08);  
    TLatex* title14_2 = new TLatex();
    title14_2->SetNDC(); 
    title14_2->SetTextSize(0.05);
    title14_2->SetTextAlign(22);  
    title14_2->DrawLatex(0.5, 0.97, "Cut |t| Distribution, Different Angles");  

    c14->Print("./figures/plot_t_distribution.pdf");



    TCanvas* c15 = new TCanvas("c15","c15",1,1,1800,800);
    c15->Divide(3,1,0.01,0.01);

    /*c14->cd(3);
    gPad->SetLogz(1);
    h_t_res_EEMC->GetXaxis()->SetTitle("|t|_{MC} [GeV/c]^{2}");
    h_t_res_EEMC->GetYaxis()->SetTitle("(|t|_{MC}-|t|_{EEMC})/|t|_{MC}");
    h_t_res_EEMC->Draw("colzsame");*/

    c15->cd(1);
    gPad->SetLogy(1);
    h_t_MC->GetXaxis()->SetLabelSize(0.03);
	h_t_MC->GetYaxis()->SetTitle("counts");
	h_t_MC->GetXaxis()->SetTitle("|t| [GeV/c]^{2}");
	h_t_MC->SetLineColor(kBlack);
    h_t_MC->Draw();
    h_t_REC_EEMC_cut->SetMarkerStyle(24);
    h_t_REC_EEMC_cut->SetMarkerColor(kBlue);
    h_t_REC_EEMC_cut->Draw("PEsame");
    h_t_REC_wRES_cut_pi12_cut->SetMarkerStyle(29);
    h_t_REC_wRES_cut_pi12_cut->SetMarkerColor(kRed);
    h_t_REC_wRES_cut_pi12_cut->Draw("PEsame");
    TLegend *w15 = new TLegend(0.52,0.7,0.7,0.85);
	w15->AddEntry(h_t_MC, "MC", "L");
    w15->AddEntry(h_t_REC_EEMC_cut, "method L", "P");
	w15->AddEntry(h_t_REC_wRES_cut_pi12_cut, "#pi/12", "P");
	w15->Draw("same");

    c15->cd(2);
    gPad->SetLogz(1);
    h_t_res_EEMC_cut_percent->GetYaxis()->SetRangeUser(-1,1);
    h_t_res_EEMC_cut_percent->GetXaxis()->SetLabelSize(0.03);
    h_t_res_EEMC_cut_percent->GetXaxis()->SetTitle("|t|_{MC} [GeV/c]^{2}");
    h_t_res_EEMC_cut_percent->GetYaxis()->SetTitle("(|t|_{MC}-|t|_{EEMC})/|t|_{MC}");
    h_t_res_EEMC_cut_percent->Draw("colzsame");

    /*
    h_t_response_EEMC_cut_0: mean of each slice
    h_t_response_EEMC_cut_1: sigma (width) of each slice
    h_t_response_EEMC_cut_2: chi2 of each fit
    */

    c15->cd(2);
    h_t_res_EEMC_cut_percent->FitSlicesX(0,0,1000);
    TH1D* h_sigma15 = (TH1D*)gDirectory->Get("h_t_res_EEMC_cut_percent_1");
    h_sigma15->Draw("same");

    /*c14->cd(5);
    gPad->SetLogz(1);
    h_t_response_EEMC->GetXaxis()->SetTitle("|t|_{EEMC} [GeV/c]^{2}");
    h_t_response_EEMC->GetYaxis()->SetTitle("|t|_{MC} [GeV/c]^{2}");
    h_t_response_EEMC->Draw("colzsame");*/

    c15->cd(3);
    gPad->SetLogz(1);
    h_t_response_EEMC_cut->GetYaxis()->SetLabelSize(0.03); 
    h_t_response_EEMC_cut->GetYaxis()->SetTitleOffset(1.4);
    h_t_response_EEMC_cut->GetXaxis()->SetLabelSize(0.03);
    h_t_response_EEMC_cut->GetXaxis()->SetTitle("|t|_{EEMC} [GeV/c]^{2}");
    h_t_response_EEMC_cut->GetYaxis()->SetTitle("|t|_{MC} [GeV/c]^{2}");
    //h_t_response_EEMC_cut->Draw("colzsame");
    h_t_response_EEMC_cut->Draw();

    c15->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title15 = new TLatex();
    title15->SetNDC(); 
    title15->SetTextSize(0.05);
    title15->SetTextAlign(22);  
    title15->DrawLatex(0.5, 0.97, "|t| Distribution");  

    c15->cd(2);  
    gPad->SetTopMargin(0.08);  
    TLatex* title15_2 = new TLatex();
    title15_2->SetNDC(); 
    title15_2->SetTextSize(0.05);
    title15_2->SetTextAlign(22);  
    title15_2->DrawLatex(0.5, 0.97, "|t| Resolution"); 
    
    c15->cd(3);  
    gPad->SetTopMargin(0.08);  
    TLatex* title15_3 = new TLatex();
    title15_3->SetNDC(); 
    title15_3->SetTextSize(0.05);
    title15_3->SetTextAlign(22);  
    title15_3->DrawLatex(0.5, 0.97, "|t| Response"); 

    c15->Print("./figures/plot_t_QAdistribution.pdf");




    TCanvas* c16 = new TCanvas("c16","c16",1,1,1200,800);
    c16->Divide(2,1,0.01,0.01);
    c16->cd(1);
    gPad->SetLogz(1);
    h_t_REC_2d->GetYaxis()->SetLabelSize(0.03); 
    h_t_REC_2d->GetYaxis()->SetTitleOffset(1.4);
    h_t_REC_2d->GetXaxis()->SetLabelSize(0.03);
    h_t_REC_2d->GetXaxis()->SetTitle("#sqrt{|t|_{x}} [GeV/c]");
    h_t_REC_2d->GetYaxis()->SetTitle("#sqrt{|t|_{y}} [GeV/c]");
    h_t_REC_2d->Draw("colzsame");

    c16->cd(2);
    gPad->SetLogz(1);
    h_t_REC_2d_wRES->GetYaxis()->SetLabelSize(0.03); 
    h_t_REC_2d_wRES->GetYaxis()->SetTitleOffset(1.4);
    h_t_REC_2d_wRES->GetXaxis()->SetLabelSize(0.03);
    h_t_REC_2d_wRES->GetXaxis()->SetTitle("#sqrt{|t|_{x}} [GeV/c]");
    h_t_REC_2d_wRES->GetYaxis()->SetTitle("#sqrt{|t|_{y}} [GeV/c]");
    h_t_REC_2d_wRES->Draw("colzsame");

    c16->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title16 = new TLatex();
    title16->SetNDC(); 
    title16->SetTextSize(0.05);
    title16->SetTextAlign(22);  
    title16->DrawLatex(0.5, 0.97, "2D |t| Distribution");  

    c16->cd(2);  
    gPad->SetTopMargin(0.08);  
    TLatex* title16_2 = new TLatex();
    title16_2->SetNDC(); 
    title16_2->SetTextSize(0.05);
    title16_2->SetTextAlign(22);  
    title16_2->DrawLatex(0.5, 0.97, "2D |t| Distribution with Resolution");  

    c16->Print("./figures/plot_t_2D.pdf");
}
