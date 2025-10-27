#include "ePICStyle.C"
gROOT->ProcessLine("SetePICStyle()");

using namespace std;

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

// cross section
    //TH1D* h_sigma = (TH1D*) file->Get("h_sigma");

// Bjorken-x
    TH2D* h_x_res_cut = (TH2D*) file->Get("h_x_res_cut");
    TH2D* h_x_response_cut = (TH2D*) file->Get("h_x_response_cut");
    TH1D* h_dxOverx_cut = (TH1D*) file->Get("h_dxOverx_cut");
    TH1D* h_x_purity = (TH1D*) file->Get("h_x_purity");
    TH1D* h_x_corrected = (TH1D*) file->Get("h_x_corrected");
    TH1D* h_x_acceptance = (TH1D*) file->Get("h_x_acceptance");
    TH1D* h_x_cut = (TH1D*) file->Get("h_x_cut");
    TH2D* h_x_migration = (TH2D*) file->Get("h_x_migration");
    TH1D* h_x_beforeCut = (TH1D*) file->Get("h_x_beforeCut");
    TH1D* h_x = (TH1D*) file->Get("h_x");
    TH1D* h_x_REC = (TH1D*) file->Get("h_x_REC");
    TH1D* h_x_afterCut = (TH1D*) file->Get("h_x_afterCut");
		
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
    TH2D* h_Q2_migration = (TH2D*) file->Get("h_Q2_migration");//new TH2D("h_Q2_migration",";Q^{2}_{RECO} [GeV/c]^{2}; Q^{2}_{MC} [GeV/c]^{2}",100,1,10,100,1,10);
    TH1D* h_Q2_purity = (TH1D*) file->Get("h_Q2_purity");//new TH1D("h_Q2_purity",";Q^{2} [GeV/c]^{2}",100,1,10);
    TH1D* h_Q2_acceptance = (TH1D*) file->Get("h_Q2_acceptance");//new TH1D("h_Q2_acceptance",";Q^{2} [GeV/c]^{2}",100,1,10);
    TH1D* h_Q2_corrected = (TH1D*) file->Get("h_Q2_corrected");//new TH1D("h_Q2_corrected",";Q^{2} [GeV/c]^{2}",100,1,10);

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
    TH2D* h_y_migration = (TH2D*) file->Get("h_y_migration"); //new TH2D("h_y_migration",";y_{RECO};y_{MC}",100,0.01,0.85,100,0.01,0.85);
    TH1D* h_y_purity = (TH1D*) file->Get("h_y_purity"); //new TH1D("h_y_purity",";y",100,0.01,0.85);
    TH1D* h_y_acceptance = (TH1D*) file->Get("h_y_acceptance"); //new TH1D("h_y_acceptance",";y",100,0.01,0.85);
    TH1D* h_y_corrected = (TH1D*) file->Get("h_y_corrected"); //new TH1D("h_y_corrected",";y",100,0.01,0.85);

// energy
	TH1D* h_energy_MC = (TH1D*) file->Get("h_energy_MC"); // before MC cut
    TH1D* h_energy_MC_after = (TH1D*) file->Get("h_energy_MC_after"); // after MC cut
	TH1D* h_energy_REC_EEMC = (TH1D*) file->Get("h_energy_REC_EEMC"); // before reco cut
	TH1D* h_energy_REC_EEMC_cut = (TH1D*) file->Get("h_energy_REC_EEMC_cut"); // cut inside position threshold
    TH2D* h_energy_res_EEMC = (TH2D*) file->Get("h_energy_res_EEMC"); // before reco cut
    TH2D* h_energy_res_EEMC_cut = (TH2D*) file->Get("h_energy_res_EEMC_cut"); // cut inside position threshold
	TH2D* h_energy_response_EEMC_cut = (TH2D*) file->Get("h_energy_response_EEMC_cut"); // cut inside position threshold
	TH2D* h_energy_response_EEMC = (TH2D*) file->Get("h_energy_response_EEMC"); // before reco cut
    TH2D* h_energy_migration = (TH2D*) file->Get("h_energy_migration"); //new TH2D("h_energy_migration",";E_{RECO} [GeV];E_{MC} [GeV]",100,0,20,100,0,20);
    TH1D* h_energy_purity = (TH1D*) file->Get("h_energy_purity"); //new TH1D("h_energy_purity",";E [GeV]",100,0,20);
    TH1D* h_energy_acceptance = (TH1D*) file->Get("h_energy_acceptance"); //new TH1D("h_energy_acceptance",";E [GeV]",100,0,20);
    TH1D* h_energy_corrected = (TH1D*) file->Get("h_energy_corrected"); //new TH1D("h_energy_corrected",";E [GeV]",100,0,20);
    TH2D* h_energy_res_EEMC_after = (TH2D*) file->Get("h_energy_res_EEMC_after");
    TH2D* h_energy_response_EEMC_after = (TH2D*) file->Get("h_energy_response_EEMC_after");
    TH1D* h_energy_REC_EEMC_after = (TH1D*) file->Get("h_energy_REC_EEMC_after");
    TH1D* h_energy_MC_mig = (TH1D*) file->Get("h_energy_MC_mig");

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
    TH2D* h_e_p_trk_res = (TH2D*) file->Get("h_e_p_trk_res");
    TH2D* h_e_pz_trk_res = (TH2D*) file->Get("h_e_pz_trk_res");

// theta
    TH1D* h_theta_MC = (TH1D*) file->Get("h_theta_MC"); // before MC cut
    TH1D* h_theta_MC_after = (TH1D*) file->Get("h_theta_MC_after"); // after MC cut
    TH1D* h_theta_REC_EEMC = (TH1D*) file->Get("h_theta_REC_EEMC");
	TH1D* h_theta_REC_EEMC_cut = (TH1D*) file->Get("h_theta_REC_EEMC_cut"); // cut inside position threshold
    TH1D* h_theta_REC_EEMC_after = (TH1D*) file->Get("h_theta_REC_EEMC_after"); // reco (after cut)
    TH1D* h_theta_diff_cut = (TH1D*) file->Get("h_theta_diff_cut"); // cut inside position threshold
    TH1D* h_theta_diff_after = (TH1D*) file->Get("h_theta_diff_after"); // reco (after cut)
	TH2D* h_theta_response_EEMC_cut = (TH2D*) file->Get("h_theta_response_EEMC_cut"); // cut inside position threshold
    TH2D* h_theta_response_EEMC_after = (TH2D*) file->Get("h_theta_response_EEMC_after"); // reco (after cut)
    TH2D* h_theta_response_EEMC_cut_clone = (TH2D*)h_theta_response_EEMC_cut->Clone("h_theta_response_EEMC_cut_clone");
    TH2D* h_theta_migration = (TH2D*) file->Get("h_theta_migration"); //new TH2D("h_theta_migration",";#theta_{RECO}; #theta_{MC}",100,0,3.14,100,0,3.14);
    TH1D* h_theta_purity = (TH1D*) file->Get("h_theta_purity"); //new TH1D("h_theta_purity",";#theta",100,0,3.14);
    TH1D* h_theta_acceptance = (TH1D*) file->Get("h_theta_acceptance"); //new TH1D("h_theta_acceptance",";#theta",100,0,3.14);
    TH1D* h_theta_corrected = (TH1D*) file->Get("h_theta_corrected"); //new TH1D("h_theta_corrected",";#theta",100,0,3.14);
    TH2D* h_theta_resolution = (TH2D*) file->Get("h_theta_resolution");
    TH1D* h_theta_resolution_Y = h_theta_resolution->ProjectionY();
    TH1D* h_theta_resolution_X = h_theta_resolution->ProjectionX();
    TH2D* h_theta_res_diff2d = (TH2D*) file->Get("h_theta_res_diff2d");
    TH1D* h_theta_res_diff = (TH1D*) file->Get("h_theta_res_diff");



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
    TH2D* h_phi_resolution = (TH2D*) file->Get("h_phi_resolution");
    TH1D* h_phi_resolution_Y = h_phi_resolution->ProjectionY();
    TH1D* h_phi_resolution_X = h_phi_resolution->ProjectionX();
    TH2D* h_phi_res_diff2d = (TH2D*) file->Get("h_phi_res_diff2d");
    TH1D* h_phi_res_diff = (TH1D*) file->Get("h_phi_res_diff");

	
// E and p
	TH2D* h_EvsP_MC = (TH2D*) file->Get("h_EvsP_MC"); 
	TH2D* h_EvsP_REC_cut = (TH2D*) file->Get("h_EvsP_REC_cut"); // cut inside position threshold
	TH2D* h_EvsP_REC = (TH2D*) file->Get("h_EvsP_REC"); // reco (after cut)
    
// E-pz
    TH1D* h_Epz_MC = (TH1D*) file->Get("h_Epz_MC"); 
    TH1D* h_Epz_MC_after = (TH1D*) file->Get("h_Epz_MC_after");
	TH1D* h_Epz_REC_cut = (TH1D*) file->Get("h_Epz_REC_cut"); // cut inside position threshold
    TH1D* h_Epz_REC = (TH1D*) file->Get("h_Epz_REC"); // reco (after cut)
    TH1D* h_Epz_afterCut = (TH1D*) file->Get("h_Epz_afterCut"); // cut inside position threshold
    TH2D* h_Epz_response = (TH2D*) file->Get("h_Epz_response"); // reco (after cut)
    TH2D* h_Epz_res = (TH2D*) file->Get("h_Epz_res"); // reco (after cut)
    TH2D* h_Epz_response_cut = (TH2D*) file->Get("h_Epz_response_cut"); // cut inside position threshold
    TH2D* h_Epz_migration = (TH2D*) file->Get("h_Epz_migration"); //new TH2D("h_Epz_migration",";(E_{RECO} - p_{z,RECO}) [GeV];(E_{MC}-p_{z,MC}) [GeV]",100,0,45,100,0,45);
    TH1D* h_Epz_purity = (TH1D*) file->Get("h_Epz_purity"); //new TH1D("h_Epz_purity",";(E - p_{z}) [GeV]",100,0,45);
    TH1D* h_Epz_acceptance = (TH1D*) file->Get("h_Epz_acceptance"); //new TH1D("h_Epz_acceptance",";(E - p_{z}) [GeV]",100,0,45);
    TH1D* h_Epz_corrected = (TH1D*) file->Get("h_Epz_corrected"); //new TH1D("h_Epz_corrected",";(E - p_{z}) [GeV]",100,0,45);
    TH1D* h_Epz_beforeCut = (TH1D*) file->Get("h_Epz_beforeCut");

// E/p
    TH1D* h_EoverP_MC = (TH1D*) file->Get("h_EoverP_MC"); // before MC cut
    TH1D* h_EoverP_MC_after = (TH1D*) file->Get("h_EoverP_MC_after"); // after MC cut
	TH1D* h_EoverP_afterCut = (TH1D*) file->Get("h_EoverP_afterCut");
	TH1D* h_EcalOverPtrk_cut = (TH1D*) file->Get("h_EcalOverPtrk_cut"); // cut inside position threshold
	TH1D* h_EcalOverPtrk = (TH1D*) file->Get("h_EcalOverPtrk"); // reco (after cut)
	TH2D* h_EoverP_response_cut = (TH2D*) file->Get("h_EoverP_response_cut"); // cut inside position threshold
    TH2D* h_EoverP_response = (TH2D*) file->Get("h_EoverP_response"); // reco (after cut)
    TH2D* h_EoverP_res = (TH2D*) file->Get("h_EoverP_res"); // reco (after cut)
    TH2D* h_EoverP_migration = (TH2D*) file->Get("h_EoverP_migration"); //new TH2D("h_EoverP_migration",";E_{EEMC}/|p|_{trk};E_{MC}/|p|_{MC}",100,0,2,100,0,2);
    TH1D* h_EoverP_purity = (TH1D*) file->Get("h_EoverP_purity"); //new TH1D("h_EoverP_purity",";E/|p|",100,0,2);
    TH1D* h_EoverP_acceptance = (TH1D*) file->Get("h_EoverP_acceptance"); //new TH1D("h_EoverP_acceptance",";E/|p|",100,0,2);
    TH1D* h_EoverP_corrected = (TH1D*) file->Get("h_EoverP_corrected"); //new TH1D("h_EoverP_corrected",";E/|p|",100,0,2);
    TH1D* h_EoverP_beforeCut = (TH1D*) file->Get("h_Epz_beforeCut");

// VM
	TH1D* h_VM_mass_MC = (TH1D*) file->Get("h_VM_mass_MC");
    TH1D* h_VM_mass_REC = (TH1D*) file->Get("h_VM_mass_REC");
    TH1D* h_VM_mass_REC_cut = (TH1D*) file->Get("h_VM_mass_REC_cut"); // cut inside position threshold
	TH1D* h_VM_pt_MC = (TH1D*) file->Get("h_VM_pt_MC");
    TH1D* h_VM_pt_REC = (TH1D*) file->Get("h_VM_pt_REC");
	TH1D* h_VM_pt_REC_cut = (TH1D*) file->Get("h_VM_pt_REC_cut"); // cut inside position threshold
    TH2D* h_VM_pt_response = (TH2D*) file->Get("h_VM_pt_response");
    TH2D* h_VM_pt_response_cut = (TH2D*) file->Get("h_VM_pt_response_cut"); // cut inside position threshold
	TH2D* h_VM_pt_migration = (TH2D*) file->Get("h_VM_pt_migration"); //new TH2D("h_VM_pt_migration",";p_{T,VM,RECO} [GeV/c];p_{T,VM,MC} [GeV/c]",200,0,2,200,0,2);
    TH1D* h_VM_pt_purity = (TH1D*) file->Get("h_VM_pt_purity"); //new TH1D("h_VM_pt_purity",";p_{T,VM} [GeV/c]",200,0,2);
    TH1D* h_VM_pt_acceptance = (TH1D*) file->Get("h_VM_pt_acceptance"); //new TH1D("h_VM_pt_acceptance",";p_{T,VM} [GeV/c]",200,0,2);
    TH1D* h_VM_pt_corrected = (TH1D*) file->Get("h_VM_pt_corrected"); //new TH1D("h_VM_pt_corrected",";p_{T,VM} [GeV/c]",200,0,2);
    TH1D* h_VM_pz_MC = (TH1D*) file->Get("h_VM_pz_MC");
    TH1D* h_VM_pz_REC = (TH1D*) file->Get("h_VM_pz_REC");
    TH1D* h_VM_pz_REC_cut = (TH1D*) file->Get("h_VM_pz_REC_cut"); // cut inside position threshold
	TH1D* h_VM_p_MC = (TH1D*) file->Get("h_VM_p_MC");
    TH1D* h_VM_p_REC = (TH1D*) file->Get("h_VM_p_REC");
    TH1D* h_VM_p_REC_cut = (TH1D*) file->Get("h_VM_p_REC_cut"); // cut inside position threshold
	TH1D* h_VM_Epz_REC_cut = (TH1D*) file->Get("h_VM_Epz_REC_cut"); // cut inside position threshold
	TH1D* h_VM_Epz_MC = (TH1D*) file->Get("h_VM_Epz_MC");
    TH1D* h_VM_Epz_REC = (TH1D*) file->Get("h_VM_Epz_REC");
    TH2D* h_VM_Epz_res = (TH2D*) file->Get("h_VM_Epz_res");
    TH2D* h_VM_Epz_response = (TH2D*) file->Get("h_VM_Epz_response");
	TH2D* h_VM_Epz_response_cut = (TH2D*) file->Get("h_VM_Epz_response_cut"); // cut inside position threshold
    TH2D* h_VM_Epz_migration = (TH2D*) file->Get("h_VM_Epz_migration");//new TH2D("h_VM_Epz_migration",";(E_{VM,RECO}-p_{z,VM,RECO}) [GeV];(E_{VM,MC}-p_{z,VM,MC} [GeV])",100,0,20,100,0,20);
    TH1D* h_VM_Epz_purity = (TH1D*) file->Get("h_VM_Epz_purity"); //new TH1D("h_VM_Epz_purity",";(E_{VM}-p_{z,VM}) [GeV]",100,0,20);
    TH1D* h_VM_Epz_acceptance = (TH1D*) file->Get("h_VM_Epz_acceptance"); //new TH1D("h_VM_Epz_acceptance",";(E_{VM}-p_{z,VM}) [GeV]",100,0,20);
    TH1D* h_VM_Epz_corrected = (TH1D*) file->Get("h_VM_Epz_corrected"); //new TH1D("h_VM_Epz_corrected",";(E_{VM}-p_{z,VM}) [GeV]",100,0,20);
    TH1D* h_VM_Epz_MC_before = (TH1D*) file->Get("h_VM_Epz_MC_before");
    TH1D* h_VM_pt_MC_before = (TH1D*) file->Get("h_VM_pt_MC_before");
    TH1D* h_VM_mass_MC_before = (TH1D*) file->Get("h_VM_mass_MC_before");
    TH1D* h_VM_pz_MC_before = (TH1D*) file->Get("h_VM_pz_MC_before");
    TH1D* h_VM_p_MC_before = (TH1D*) file->Get("h_VM_p_MC_before");
    TH1D* h_VM_mass_REC_before = (TH1D*) file->Get("h_VM_mass_REC_before"); //new TH1D("h_VM_mass_REC_before",";m_{VM,MC} [GeV/c]",200,0,2);
    TH1D* h_VM_pt_REC_before = (TH1D*) file->Get("h_VM_pt_REC_before");//new TH1D("h_VM_pt_REC_before",";p_{T,VM,MC} [GeV/c]",200,0,20);
    TH1D* h_VM_pz_REC_before = (TH1D*) file->Get("h_VM_pz_REC_before");//new TH1D("h_VM_pz_REC_before",";p_{z,VM,MC} [GeV/c]",200,0,20);
    TH1D* h_VM_p_REC_before = (TH1D*) file->Get("h_VM_p_REC_before");//new TH1D("h_VM_p_REC_before",";|p|_{VM,MC} [GeV/c]",200,0,20);
    TH1D* h_VM_Epz_REC_before = (TH1D*) file->Get("h_VM_Epz_REC_before");//new TH1D("h_VM_Epz_MC_before",";(E_{VM,MC}-p_{z,VM,MC}) [GeV]",200,0,20);
    TH2D* h_VM_pt_resolution = (TH2D*) file->Get("h_VM_pt_resolution");
    TH1D* h_VM_pt_resolution_Y = h_VM_pt_resolution->ProjectionY();
    TH1D* h_VM_pt_resolution_X = h_VM_pt_resolution->ProjectionX();
    TH2D* h_VM_pt_res_diff2d = (TH2D*) file->Get("h_VM_pt_res_diff2d");
    TH1D* h_VM_pt_res_diff = (TH1D*) file->Get("h_VM_pt_res_diff");
    TH2D* h_VM_pz_resolution = (TH2D*) file->Get("h_VM_pz_resolution");
    TH1D* h_VM_pz_resolution_Y = h_VM_pz_resolution->ProjectionY();
    TH1D* h_VM_pz_resolution_X = h_VM_pz_resolution->ProjectionX();
    TH2D* h_VM_pz_res_diff2d = (TH2D*) file->Get("h_VM_pz_res_diff2d");
    TH1D* h_VM_pz_res_diff = (TH1D*) file->Get("h_VM_pz_res_diff");
    TH2D* h_VM_p_resolution = (TH2D*) file->Get("h_VM_p_resolution");
    TH1D* h_VM_p_resolution_Y = h_VM_p_resolution->ProjectionY();
    TH1D* h_VM_p_resolution_X = h_VM_p_resolution->ProjectionX();
    TH2D* h_VM_p_res_diff2d = (TH2D*) file->Get("h_VM_p_res_diff2d");
    TH1D* h_VM_p_res_diff = (TH1D*) file->Get("h_VM_p_res_diff");
    TH1D* h_VM_theta_MC = (TH1D*) file->Get("h_VM_theta_MC");
    TH1D* h_VM_theta_REC = (TH1D*) file->Get("h_VM_theta_REC");



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
    TH1D* h_t_MC_before = (TH1D*) file->Get("h_t_MC_before");
    TH1D* h_t_REC_EEMC = (TH1D*) file->Get("h_t_REC_EEMC"); // after MC cuts
	TH1D* h_t_REC_EEMC_cut = (TH1D*) file->Get("h_t_REC_EEMC_cut"); // cut inside position threshold
    TH1D* h_t_REC_EEMC_cut2 = (TH1D*) file->Get("h_t_REC_EEMC_cut2"); // cut outside position threshold
	TH2D* h_t_res_EEMC_cut = (TH2D*) file->Get("h_t_res_EEMC_cut"); // cut inside position threshold
    TH2D* h_t_res_EEMC_cut_percent = (TH2D*) file->Get("h_t_res_EEMC_cut_percent"); // cut inside position threshold
    TH2D* h_t_res_EEMC = (TH2D*) file->Get("h_t_res_EEMC"); // after MC cuts
	TH2D* h_t_response_EEMC_cut = (TH2D*) file->Get("h_t_response_EEMC_cut"); // cut inside position threshold
    TH2D* h_t_response_EEMC = (TH2D*) file->Get("h_t_response_EEMC"); // after MC cuts
    TH2D* h_t_migration = (TH2D*) file->Get("h_t_migration"); //new TH2D("h_t_migration",";|t|_{RECO} [GeV/c]^{2};|t|_{MC} [GeV/c]^{2}",100,0,0.2,100,0,0.2);
    TH1D* h_t_purity = (TH1D*) file->Get("h_t_purity"); //new TH1D("h_t_purity",";|t| [GeV/c]^{2}",100,0,0.2);
    TH1D* h_t_acceptance = (TH1D*) file->Get("h_t_acceptance"); //new TH1D("h_t_acceptance",";|t| [GeV/c]^{2}",100,0,0.2);
    TH1D* h_t_corrected = (TH1D*) file->Get("h_t_corrected"); //new TH1D("h_t_corrected",";|t| [GeV/c]^{2}",100,0,0.2);
    TH2D* h_t_response_cut = (TH2D*) file->Get("h_t_response_cut");
    TH2D* h_t_res = (TH2D*) file->Get("h_t_res");

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

// t incoherent
    TH1D* h_t_incoherent = (TH1D*) file->Get("h_t_incoherent_REC");
    TH1D* h_t_incoherent_REC = (TH1D*) file->Get("h_t_incoherent_REC");
    TH1D* h_t_incoherent_REC_cut = (TH1D*) file->Get("h_t_incoherent_REC_cut");

// extra analysis
    TH1D* h_angle_resolution = (TH1D*) file->Get("h_angle_resolution");
    TH1D* h_solid_angle_REC = (TH1D*) file->Get("h_solid_angle_REC");
    TH2D* h_phi_response = (TH2D*) file->Get("h_phi_response");
    TH2D* h_theta_response = (TH2D*) file->Get("h_theta_response");
    TH1D* h_omega_MC = (TH1D*) file->Get("h_omega_MC");
    TH1D* h_omega_REC = (TH1D*) file->Get("h_omega_REC");
    TH2D* h_Q2_vs_x_MC = (TH2D*) file->Get("h_Q2_vs_x_MC");
    TH2D* h_Q2_vs_x_REC = (TH2D*) file->Get("h_Q2_vs_x_REC");


// check gap in y
    TCanvas* c1234 = new TCanvas("c1234","c1234",1,1,800,800);
    c1234->Divide(2,3,0.01,0.01);
    c1234->cd(1);
    //gPad->SetLogy(1);
	h_y_e->GetYaxis()->SetTitle("counts");
	h_y_e->GetXaxis()->SetTitle("y");
	h_y_e->SetLineColor(kBlack);
    h_y_e->Draw();
    h_y_beforeCut->SetMarkerStyle(11);
    h_y_beforeCut->SetMarkerColor(kOrange);
    //h_y_beforeCut->Draw("Psame");
    h_y_afterCut->SetMarkerStyle(24);
    h_y_afterCut->SetMarkerColor(kBlue);
    h_y_afterCut->Draw("PEsame");
    h_yREC_e_EEMC_cut->SetMarkerColor(kBlue);
    h_yREC_e_EEMC_cut->SetMarkerStyle(24);
    //h_yREC_e_EEMC_cut->Draw("PEsame");
    h_yREC_e_EEMC->SetMarkerStyle(7);
    h_yREC_e_EEMC->SetMarkerColor(kRed);
    h_yREC_e_EEMC->Draw("PEsame");
    TLegend *w1c1234 = new TLegend(0.32,0.7,0.6,0.85);
	w1c1234->AddEntry(h_y_e, " MC", "L");
    //w1->AddEntry(h_Q2_beforeCut, "MC (before cut)", "P");
	w1c1234->AddEntry(h_yREC_e_EEMC, " RECO without cut", "P");
    w1c1234->AddEntry(h_y_afterCut, " RECO", "P");
    //w1->AddEntry(h_Q2REC_e_EEMC_cut, "RECO", "P");
	w1c1234->Draw("same");

    c1234->cd(2);
    gPad->SetLogy(1);
	h_y_e->GetYaxis()->SetTitle("counts");
	h_y_e->GetXaxis()->SetTitle("y");
	h_y_e->SetLineColor(kBlack);
    h_y_e->Draw();
    h_y_beforeCut->SetMarkerStyle(11);
    h_y_beforeCut->SetMarkerColor(kOrange);
    //h_y_beforeCut->Draw("Psame");
    h_y_afterCut->SetMarkerStyle(24);
    h_y_afterCut->SetMarkerColor(kBlue);
    h_y_afterCut->Draw("PEsame");
    h_yREC_e_EEMC_cut->SetMarkerColor(kBlue);
    h_yREC_e_EEMC_cut->SetMarkerStyle(24);
    //h_yREC_e_EEMC_cut->Draw("PEsame");
    h_yREC_e_EEMC->SetMarkerStyle(7);
    h_yREC_e_EEMC->SetMarkerColor(kRed);
    h_yREC_e_EEMC->Draw("PEsame");
	w1c1234->Draw("same");


    c1234->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title1_1_123c1234 = new TLatex();
    title1_1_123c1234->SetNDC(); 
    title1_1_123c1234->SetTextSize(0.05);
    title1_1_123c1234->SetTextAlign(22);  
    title1_1_123c1234->DrawLatex(0.5, 0.97, "y Truth vs Reco"); 

    c1234->cd(2);  
    gPad->SetTopMargin(0.08); 
    TLatex* title21_1_123c1234 = new TLatex();
    title21_1_123c1234->SetNDC(); 
    title21_1_123c1234->SetTextSize(0.05);
    title21_1_123c1234->SetTextAlign(22);  
    title21_1_123c1234->DrawLatex(0.5, 0.97, "y Truth vs Reco (Log Scale)"); 
    
    c1234->cd(3);
    //gPad->SetLogy(1);
    h_Q2_e->GetXaxis()->SetRangeUser(1,10);
	h_Q2_e->GetYaxis()->SetTitle("counts");
	h_Q2_e->GetXaxis()->SetTitle("Q^{2} [GeV/c]^{2}");
	h_Q2_e->SetLineColor(kBlack);
    h_Q2_e->Draw();
    h_Q2_beforeCut->SetMarkerStyle(11);
    h_Q2_beforeCut->SetMarkerColor(kOrange);
    //h_Q2_beforeCut->Draw("PEsame");
    h_Q2REC_e_EEMC->SetMarkerColor(kRed);
    h_Q2REC_e_EEMC->SetMarkerStyle(7);
    h_Q2REC_e_EEMC->Draw("PEsame");
    h_Q2_afterCut->SetMarkerStyle(24);
    h_Q2_afterCut->SetMarkerColor(kBlue);
    h_Q2_afterCut->Draw("PEsame");
   // h_Q2REC_e_EEMC_cut->SetMarkerStyle(21);
   // h_Q2REC_e_EEMC_cut->SetMarkerColor(kMagenta);
   // h_Q2REC_e_EEMC_cut->Draw("PEsame");
    w1c1234->Draw("same");

    c1234->cd(4);
    gPad->SetLogy(1);
    h_Q2_e->GetXaxis()->SetRangeUser(1,10);
	h_Q2_e->GetYaxis()->SetTitle("counts");
	h_Q2_e->GetXaxis()->SetTitle("Q^{2} [GeV/c]^{2}");
	h_Q2_e->SetLineColor(kBlack);
    h_Q2_e->Draw();
    h_Q2_beforeCut->SetMarkerStyle(11);
    h_Q2_beforeCut->SetMarkerColor(kOrange);
    //h_Q2_beforeCut->Draw("PEsame");
    h_Q2REC_e_EEMC->SetMarkerColor(kRed);
    h_Q2REC_e_EEMC->SetMarkerStyle(7);
    h_Q2REC_e_EEMC->Draw("PEsame");
    h_Q2_afterCut->SetMarkerStyle(24);
    h_Q2_afterCut->SetMarkerColor(kBlue);
    h_Q2_afterCut->Draw("PEsame");
   // h_Q2REC_e_EEMC_cut->SetMarkerStyle(21);
   // h_Q2REC_e_EEMC_cut->SetMarkerColor(kMagenta);
   // h_Q2REC_e_EEMC_cut->Draw("PEsame");
    w1c1234->Draw("same");

    c1234->cd(3);  
    gPad->SetTopMargin(0.08); 
    TLatex* title1_1_123c123 = new TLatex();
    title1_1_123c123->SetNDC(); 
    title1_1_123c123->SetTextSize(0.05);
    title1_1_123c123->SetTextAlign(22);  
    title1_1_123c123->DrawLatex(0.5, 0.97, "Q^{2} Truth vs Reco"); 

    c1234->cd(4);  
    gPad->SetTopMargin(0.08); 
    TLatex* title1_14_123c123 = new TLatex();
    title1_14_123c123->SetNDC(); 
    title1_14_123c123->SetTextSize(0.05);
    title1_14_123c123->SetTextAlign(22);  
    title1_14_123c123->DrawLatex(0.5, 0.97, "Q^{2} Truth vs Reco (Log Scale)"); 

    c1234->cd(5);
    //gPad->SetLogy(1);
    h_x->GetXaxis()->SetNdivisions(505);
	h_x->GetYaxis()->SetTitle("counts");
	h_x->GetXaxis()->SetTitle("x");
	h_x->SetLineColor(kBlack);
    h_x->Draw();
    h_x_beforeCut->SetMarkerStyle(11);
    h_x_beforeCut->SetMarkerColor(kOrange);
    //h_x_beforeCut->Draw("PEsame");
    h_x_afterCut->SetMarkerStyle(24);
    h_x_afterCut->SetMarkerColor(kBlue);
    h_x_afterCut->Draw("PEsame");
    h_x_cut->SetMarkerStyle(24);
    h_x_cut->SetMarkerColor(kBlue);
    //h_x_cut->Draw("PEsame");
    h_x_REC->SetMarkerStyle(7);
    h_x_REC->SetMarkerColor(kRed);
    h_x_REC->Draw("PEsame");
    w1c1234->Draw("same");

    c1234->cd(6);
    gPad->SetLogy(1);
    h_x->GetXaxis()->SetNdivisions(505);
	h_x->GetYaxis()->SetTitle("counts");
	h_x->GetXaxis()->SetTitle("x");
	h_x->SetLineColor(kBlack);
    h_x->Draw();
    h_x_beforeCut->SetMarkerStyle(11);
    h_x_beforeCut->SetMarkerColor(kOrange);
    //h_x_beforeCut->Draw("PEsame");
    h_x_afterCut->SetMarkerStyle(24);
    h_x_afterCut->SetMarkerColor(kBlue);
    h_x_afterCut->Draw("PEsame");
    h_x_cut->SetMarkerStyle(24);
    h_x_cut->SetMarkerColor(kBlue);
    //h_x_cut->Draw("PEsame");
    h_x_REC->SetMarkerStyle(7);
    h_x_REC->SetMarkerColor(kRed);
    h_x_REC->Draw("PEsame");
    w1c1234->Draw("same");

    c1234->cd(5);  
    gPad->SetTopMargin(0.08); 
    TLatex* title1_1_123c12 = new TLatex();
    title1_1_123c12->SetNDC(); 
    title1_1_123c12->SetTextSize(0.05);
    title1_1_123c12->SetTextAlign(22);  
    title1_1_123c12->DrawLatex(0.5, 0.97, "x Truth vs Reco");
    
    c1234->cd(6);  
    gPad->SetTopMargin(0.08); 
    TLatex* title1_156_123c12 = new TLatex();
    title1_156_123c12->SetNDC(); 
    title1_156_123c12->SetTextSize(0.05);
    title1_156_123c12->SetTextAlign(22);  
    title1_156_123c12->DrawLatex(0.5, 0.97, "x Truth vs Reco (Log Scale)");

    //c1234->cd(4);
    //gPad->SetLogy(1);
    h_Epz_MC->GetXaxis()->SetRangeUser(15,25);
	h_Epz_MC->GetYaxis()->SetTitle("counts");
	h_Epz_MC->GetXaxis()->SetTitle("(E - p_{z}) [GeV]");
	h_Epz_MC->SetLineColor(kBlack);
    //h_Epz_MC->Draw();
    //h_Epz_MC_after->SetMarkerStyle(6);
    h_Epz_MC_after->SetLineColor(kBlack);
    //h_Epz_MC_after->Draw();
    h_Epz_REC->SetMarkerColor(kBlue);
    h_Epz_REC->SetMarkerStyle(24);
    //h_Epz_REC->Draw("PEsame");
    TLegend *w1c12344 = new TLegend(0.55,0.7,0.65,0.85);
	w1c12344->AddEntry(h_Epz_MC_after, " MC", "L");
    //w1->AddEntry(h_Q2_beforeCut, "MC (before cut)", "P");
	//w1->AddEntry(h_Q2REC_e_EEMC, " RECO before cut", "P");
    w1c12344->AddEntry(h_Epz_REC, " RECO", "P");
    //w1->AddEntry(h_Q2REC_e_EEMC_cut, "RECO", "P");
	//w1c12344->Draw("same");

    //c1234->cd(4);  
    gPad->SetTopMargin(0.08); 
    TLatex* title1_1_123c1 = new TLatex();
    title1_1_123c1->SetNDC(); 
    title1_1_123c1->SetTextSize(0.05);
    title1_1_123c1->SetTextAlign(22);  
    //title1_1_123c1->DrawLatex(0.5, 0.97, "(E-p_{z}) Truth vs Reco"); 

    //c1234->cd(5);
    //gPad->SetLogy(1);
    h_EoverP_MC->SetTitle("E/|p| Truth vs. Reco");
	h_EoverP_MC->GetYaxis()->SetTitle("counts");
	h_EoverP_MC->GetXaxis()->SetTitle("E/|p|");
	h_EoverP_MC->SetLineColor(kBlack);
    //h_EoverP_MC->Draw();
    h_EoverP_MC_after->SetLineColor(kBlack);
    //h_EoverP_MC_after->Draw();
    h_EcalOverPtrk->SetMarkerColor(kBlue);
    h_EcalOverPtrk->SetMarkerStyle(24);
    //h_EcalOverPtrk->Draw("PEsame");
    //h_EcalOverPtrk_cut->SetMarkerStyle(24);
    //h_EcalOverPtrk_cut->SetMarkerColor(kBlue);
    //h_EcalOverPtrk_cut->Draw("PEsame");
    //h_EoverP_beforeCut->SetMarkerStyle(7);
    //h_EoverP_beforeCut->SetMarkerColor(kRed);
    //h_EoverP_beforeCut->Draw("PEsame");
    h_EoverP_afterCut->SetMarkerStyle(24);
    h_EoverP_afterCut->SetMarkerColor(kBlue);
    //h_EoverP_afterCut->Draw("PEsame");
    //w1c12344->Draw("same");

    //c1234->cd(5);  
    gPad->SetTopMargin(0.08); 
    TLatex* title1_1_123c = new TLatex();
    title1_1_123c->SetNDC(); 
    title1_1_123c->SetTextSize(0.05);
    title1_1_123c->SetTextAlign(22);  
    //title1_1_123c->DrawLatex(0.5, 0.97, "E/p Truth vs Reco"); 

    //c1234->cd(4);
    gPad->SetLogy(1);
	h_eta_MC->GetYaxis()->SetTitle("counts");
	h_eta_MC->GetXaxis()->SetTitle("#eta");
	h_eta_MC->SetLineColor(kBlack);
    //h_eta_MC->Draw();
    //h_eta_MC_after->SetMarkerStyle(11);
    h_eta_MC_after->SetLineColor(kBlack);
    //h_eta_MC_after->Draw();
    h_eta_REC_EEMC->SetMarkerColor(kRed);
    h_eta_REC_EEMC->SetMarkerStyle(7);
    //h_eta_REC_EEMC->Draw("PEsame");
    //h_eta_REC_EEMC_cut->SetMarkerStyle(29);
    //h_eta_REC_EEMC_cut->SetMarkerColor(kRed);
    //h_eta_REC_EEMC_cut->Draw("PEsame");
    h_eta_REC_EEMC_after->SetMarkerStyle(24);
    h_eta_REC_EEMC_after->SetMarkerColor(kBlue);
    //h_eta_REC_EEMC_after->Draw("PEsame");
    //w1c1234->Draw("same");

    //c1234->cd(4);  
    gPad->SetTopMargin(0.08); 
    TLatex* title1_1_123c12341 = new TLatex();
    title1_1_123c12341->SetNDC(); 
    title1_1_123c12341->SetTextSize(0.05);
    title1_1_123c12341->SetTextAlign(22);  
    //title1_1_123c12341->DrawLatex(0.5, 0.97, "#eta Truth vs Reco"); 

    //c1234->cd(5);
    gPad->SetLogy(1);
	h_phi_MC->GetYaxis()->SetTitle("counts");
	h_phi_MC->GetXaxis()->SetTitle("#phi");
	h_phi_MC->SetLineColor(kBlack);
    //h_phi_MC->Draw();
    h_phi_REC_EEMC->SetMarkerColor(kRed);
    h_phi_REC_EEMC->SetMarkerStyle(7);
    //h_phi_REC_EEMC->Draw("PEsame");
    //h_phi_MC_after->SetMarkerStyle(11);
    h_phi_MC_after->SetLineColor(kBlack);
    //h_phi_MC_after->Draw();
    //h_phi_REC_EEMC_cut->SetMarkerStyle(29);
    //h_phi_REC_EEMC_cut->SetMarkerColor(kRed);
    //h_phi_REC_EEMC_cut->Draw("PEsame");
    h_phi_REC_EEMC_after->SetMarkerColor(kBlue);
    h_phi_REC_EEMC_after->SetMarkerStyle(24);
    //h_phi_REC_EEMC_after->Draw("PEsame");
    //w1c1234->Draw("same");

    //c1234->cd(5);  
    gPad->SetTopMargin(0.08); 
    TLatex* title1_1_123c12342 = new TLatex();
    title1_1_123c12342->SetNDC(); 
    title1_1_123c12342->SetTextSize(0.05);
    title1_1_123c12342->SetTextAlign(22);  
    //title1_1_123c12342->DrawLatex(0.5, 0.97, "#phi Truth vs Reco"); 

    //c1234->cd(6);
    gPad->SetLogy(1);
	h_theta_MC_after->GetYaxis()->SetTitle("counts");
	h_theta_MC->GetXaxis()->SetTitle("#phi");
	h_theta_MC->SetLineColor(kBlack);
    //h_theta_MC->Draw();
    //h_theta_MC_after->GetXaxis()->SetRangeUser(-1,3.5);
    //h_theta_MC_after->GetYaxis()->SetTitle("counts");
	//h_theta_MC_after->GetXaxis()->SetTitle("#theta");
    h_theta_MC_after->SetLineColor(kBlack);
    //h_theta_MC_after->SetMarkerStyle(11);
    //h_theta_MC_after->Draw();
    //h_theta_REC_EEMC_cut->SetMarkerColor(kBlue);
    //h_theta_REC_EEMC_cut->SetMarkerStyle(24);
    //h_theta_REC_EEMC_cut->Draw("PEsame");
    h_theta_REC_EEMC->SetMarkerColor(kRed);
    h_theta_REC_EEMC->SetMarkerStyle(7);
    //h_theta_REC_EEMC->Draw("PEsame");
    h_theta_REC_EEMC_after->SetMarkerColor(kBlue);
    h_theta_REC_EEMC_after->SetMarkerStyle(24);
    //h_theta_REC_EEMC_after->Draw("PEsame");
    //w1c1234->Draw("same");

    //c1234->cd(6);  
    gPad->SetTopMargin(0.08); 
    TLatex* title1_1_123c12343 = new TLatex();
    title1_1_123c12343->SetNDC(); 
    title1_1_123c12343->SetTextSize(0.05);
    title1_1_123c12343->SetTextAlign(22);  
    //title1_1_123c12343->DrawLatex(0.5, 0.97, "#theta Truth vs Reco"); 

    c1234->Print("./figures/plot_compare_with_y.pdf");


// cross section
    //TCanvas* canva = new TCanvas("canva","canva",1,1,1600,800);
    //gPad->SetLogy(1);
    //h_sigma->Draw();
    //canva->Print("./figures/plot_sigma.pdf");

// Bjorken-x
	TCanvas* c1_123 = new TCanvas("c1_123","c1_123",1,1,1600,600);
    c1_123->Divide(3,1,0.01,0.01);
    c1_123->cd(1);
    gPad->SetLogy(1);
	h_x->GetYaxis()->SetTitle("counts");
	h_x->GetXaxis()->SetTitle("x");
	h_x->SetLineColor(kBlack);
    h_x->Draw();
    h_x_beforeCut->SetMarkerStyle(11);
    h_x_beforeCut->SetMarkerColor(kOrange);
    //h_x_beforeCut->Draw("PEsame");
    h_x_afterCut->SetMarkerStyle(24);
    h_x_afterCut->SetMarkerColor(kBlue);
    h_x_afterCut->Draw("PEsame");
    h_x_cut->SetMarkerStyle(24);
    h_x_cut->SetMarkerColor(kBlue);
    //h_x_cut->Draw("PEsame");
    h_x_REC->SetMarkerStyle(6);
    h_x_REC->SetMarkerColor(kRed);
    //h_x_REC->Draw("PEsame");
    TLegend *w1_123 = new TLegend(0.42,0.7,0.6,0.85);
	w1_123->AddEntry(h_x, " MC", "L");
    //w1_123->AddEntry(h_x_beforeCut, "MC (before cut)", "P");
    //w1_123->AddEntry(h_x_cut, "RECO", "P");
    w1_123->AddEntry(h_x_afterCut, " RECO", "P");
    //w1_123->AddEntry(h_x_REC, "RECO before cut", "P");
	w1_123->Draw("same");

    c1_123->cd(2);
    gPad->SetLogz(1);
    h_x_res_cut->GetXaxis()->SetRangeUser(0,0.15);
    h_x_res_cut->GetYaxis()->SetTitleOffset(1.5);
    h_x_res_cut->GetYaxis()->SetTitle("(x_{MC}-x_{RECO})/x_{MC}");
    h_x_res_cut->GetXaxis()->SetTitle("x_{MC}");
    h_x_res_cut->Draw("colzsame");

    c1_123->cd(3);
    gPad->SetLogz(1);
    h_x_response_cut->GetXaxis()->SetRangeUser(0,0.25);
    h_x_response_cut->GetYaxis()->SetRangeUser(0,0.15);
    h_x_response_cut->GetYaxis()->SetTitle("x_{MC}");
    h_x_response_cut->GetXaxis()->SetTitle("x_{RECO}");
    h_x_response_cut->Draw("colzsame");

    //c1_123->cd(4);
    //gPad->SetLogy(1);
    //h_dxOverx_cut->GetYaxis()->SetTitle("counts");
    //h_dxOverx_cut->GetXaxis()->SetTitle("dx/x_{MC}");
    //h_dxOverx_cut->Draw();

    //c1_123->cd(5);
    //gPad->SetLogy(1);
    //h_x_purity->GetYaxis()->SetTitle("bins");
    //h_x_purity->Draw("same");

    //c1_123->cd(6);
    //gPad->SetLogz(1);
    //h_x_migration->Draw("colzsame");

    //c1_123->cd(7);
    //gPad->SetLogy(1);
    //h_x->GetYaxis()->SetTitle("counts");
	//h_x->GetXaxis()->SetTitle("x");
	//h_x->SetLineColor(kBlack);
    //h_x->Draw();
    //h_x_corrected->SetMarkerStyle(24);
    //h_x_corrected->SetMarkerColor(kBlue);
    //h_x_corrected->Draw("PEsame");
    //TLegend *w1_1_123 = new TLegend(0.42,0.7,0.6,0.85);
	//w1_1_123->AddEntry(h_x, "MC", "L");
    //w1_1_123->AddEntry(h_x_corrected, "Corrected", "P");
	//w1_1_123->Draw("same");

    //c1_123->cd(8);
    //gPad->SetLogy(1);
    //h_x_acceptance->Draw("same");

    c1_123->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title1_1_123 = new TLatex();
    title1_1_123->SetNDC(); 
    title1_1_123->SetTextSize(0.05);
    title1_1_123->SetTextAlign(22);  
    title1_1_123->DrawLatex(0.5, 0.97, "x Truth vs Reco");  

    c1_123->cd(2);  
    gPad->SetTopMargin(0.08);  
    TLatex* title1_2_123 = new TLatex();
    title1_2_123->SetNDC(); 
    title1_2_123->SetTextSize(0.05);
    title1_2_123->SetTextAlign(22);  
    title1_2_123->DrawLatex(0.5, 0.97, "x Resolution");

    c1_123->cd(3);  
    gPad->SetTopMargin(0.08);  
    TLatex* title1_3_123 = new TLatex();
    title1_3_123->SetNDC(); 
    title1_3_123->SetTextSize(0.05);
    title1_3_123->SetTextAlign(22);  
    title1_3_123->DrawLatex(0.5, 0.97, "x Response");

    c1_123->cd(5);  
    gPad->SetTopMargin(0.08);  
    TLatex* title1_4_123 = new TLatex();
    title1_4_123->SetNDC(); 
    title1_4_123->SetTextSize(0.05);
    title1_4_123->SetTextAlign(22);  
    //title1_4_123->DrawLatex(0.5, 0.97, "Purity");

    c1_123->cd(6);  
    gPad->SetTopMargin(0.08);  
    TLatex* title1_5_123 = new TLatex();
    title1_5_123->SetNDC(); 
    title1_5_123->SetTextSize(0.05);
    title1_5_123->SetTextAlign(22);  
    //title1_5_123->DrawLatex(0.5, 0.97, "Bin Migration");

    c1_123->cd(7);  
    gPad->SetTopMargin(0.08);  
    TLatex* title1_6_123 = new TLatex();
    title1_6_123->SetNDC(); 
    title1_6_123->SetTextSize(0.05);
    title1_6_123->SetTextAlign(22);  
    //title1_6_123->DrawLatex(0.5, 0.97, "Corrected Acceptance");

    c1_123->cd(8);  
    gPad->SetTopMargin(0.08);  
    TLatex* title1_7_123 = new TLatex();
    title1_7_123->SetNDC(); 
    title1_7_123->SetTextSize(0.05);
    title1_7_123->SetTextAlign(22);  
    //title1_7_123->DrawLatex(0.5, 0.97, "Acceptance");

    c1_123->Print("./figures/plot_x.pdf");


// Q2
	TCanvas* c1 = new TCanvas("c1","c1",1,1,1800,600);
    c1->Divide(3,1,0.01,0.01);
    c1->cd(1);
    gPad->SetLogy(1);
    h_Q2_e->GetXaxis()->SetRangeUser(1,10);
	h_Q2_e->GetYaxis()->SetTitle("counts");
	h_Q2_e->GetXaxis()->SetTitle("Q^{2} [GeV/c]^{2}");
	h_Q2_e->SetLineColor(kBlack);
    h_Q2_e->Draw();
    h_Q2_beforeCut->SetMarkerStyle(11);
    h_Q2_beforeCut->SetMarkerColor(kOrange);
    //h_Q2_beforeCut->Draw("PEsame");
    h_Q2REC_e_EEMC->SetMarkerColor(kRed);
    h_Q2REC_e_EEMC->SetMarkerStyle(6);
    //h_Q2REC_e_EEMC->Draw("PEsame");
    h_Q2_afterCut->SetMarkerStyle(24);
    h_Q2_afterCut->SetMarkerColor(kBlue);
    h_Q2_afterCut->Draw("PEsame");
   // h_Q2REC_e_EEMC_cut->SetMarkerStyle(21);
   // h_Q2REC_e_EEMC_cut->SetMarkerColor(kMagenta);
   // h_Q2REC_e_EEMC_cut->Draw("PEsame");
    TLegend *w1 = new TLegend(0.42,0.7,0.6,0.85);
	w1->AddEntry(h_Q2_e, " MC", "L");
    //w1->AddEntry(h_Q2_beforeCut, "MC (before cut)", "P");
	//w1->AddEntry(h_Q2REC_e_EEMC, " RECO before cut", "P");
    w1->AddEntry(h_Q2_afterCut, " RECO", "P");
    //w1->AddEntry(h_Q2REC_e_EEMC_cut, "RECO", "P");
	w1->Draw("same");

    c1->cd(2);
    gPad->SetLogz(1);
    h_Q2_res->GetXaxis()->SetRangeUser(1,10);
    //h_Q2_res->GetYaxis()->SetRangeUser(-0.4,0.4);
    h_Q2_res->GetYaxis()->SetTitle("(Q^{2}_{MC}-Q^{2}_{RECO})/Q^{2}_{MC}");
    h_Q2_res->GetXaxis()->SetTitle("Q^{2}_{MC} [GeV/c]^{2}");
    h_Q2_res->Draw("colzsame");

    //c1->cd(2);
    //gPad->SetLogz(1);
    //h_Q2_res_cut->GetYaxis()->SetTitleOffset(1.5);
    //h_Q2_res_cut->GetYaxis()->SetTitle("(Q^{2}_{MC}-Q^{2}_{RECO})/Q^{2}_{MC}");
    //h_Q2_res_cut->GetXaxis()->SetTitle("Q^{2}_{MC} [GeV/c]^{2}");
    //h_Q2_res_cut->Draw("colzsame");

    c1->cd(3);
    gPad->SetLogz(1);
    h_Q2_response->GetXaxis()->SetRangeUser(1,10);
    h_Q2_response->GetYaxis()->SetRangeUser(1,10);
    h_Q2_response->GetYaxis()->SetTitleOffset(1.5);
    h_Q2_response->GetYaxis()->SetTitle("Q^{2}_{MC} [GeV/c]^{2}");
    h_Q2_response->GetXaxis()->SetTitle("Q^{2}_{RECO} [GeV/c]^{2}");
    h_Q2_response->Draw("colzsame");

    //c1->cd(4);
    //gPad->SetLogz(1);
    //h_Q2_vs_x_MC->Draw("colzsame");

    //c1->cd(5);
    //gPad->SetLogz(1);
    //h_Q2_vs_x_REC->Draw("colzsame");

    //c1->cd(3);
    //gPad->SetLogz(1);
    //h_Q2_response_cut->GetYaxis()->SetTitleOffset(1.5);
    //h_Q2_response_cut->GetYaxis()->SetTitle("Q^{2}_{MC} [GeV/c]^{2}");
    //h_Q2_response_cut->GetXaxis()->SetTitle("Q^{2}_{RECO} [GeV/c]^{2}");
    //h_Q2_response_cut->Draw("colzsame");

    //c1->cd(8);
    //gPad->SetLogy(1);
    //h_dQ2overQ2_REC->GetYaxis()->SetTitle("counts");
    //h_dQ2overQ2_REC->GetXaxis()->SetTitle("dQ^{2}/Q^{2}_{MC}");
    //h_dQ2overQ2_REC->Draw();

    //c1->cd(8);
    //gPad->SetLogy(1);
    //h_dQ2overQ2_REC_cut->GetYaxis()->SetTitle("counts");
    //h_dQ2overQ2_REC_cut->GetXaxis()->SetTitle("dQ^{2}/Q^{2}_{MC}");
    //h_dQ2overQ2_REC_cut->Draw();

    //c1->cd(4);
    //gPad->SetLogy(1);
    //h_Q2_purity->SetMinimum(0);
    //h_Q2_purity->SetMaximum(1);
    //h_Q2_purity->GetYaxis()->SetTitle("Purity");
    //h_Q2_purity->Draw();

    //c1->cd(5);
    //gPad->SetLogz(1);
    //h_Q2_migration->Draw("colzsame");

    //c1->cd(6);
    //gPad->SetLogy(1);
    //h_Q2_e->GetYaxis()->SetTitle("Counts");
	//h_Q2_e->GetXaxis()->SetTitle("Q^{2} [GeV/c]^{2}");
	//h_Q2_e->SetLineColor(kBlack);
    //h_Q2_e->Draw();
    //h_Q2_corrected->SetMarkerStyle(24);
    //h_Q2_corrected->SetMarkerColor(kBlue);
    //h_Q2_corrected->Draw("PEsame");
    //TLegend *w1_1 = new TLegend(0.42,0.7,0.6,0.85);
	//w1_1->AddEntry(h_Q2_e, "MC", "L");
    //w1_1->AddEntry(h_Q2_corrected, "Corrected", "P");
	//w1_1->Draw("same");

    //c1->cd(7);
    //gPad->SetLogy(1);
    //h_Q2_acceptance->SetMinimum(0);
    //h_Q2_acceptance->SetMaximum(1);
    //h_Q2_acceptance->GetYaxis()->SetTitle("Acceptance");
    //h_Q2_acceptance->Draw("same");

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

    c1->cd(4);  
    gPad->SetTopMargin(0.08);  
    TLatex* title1_4 = new TLatex();
    title1_4->SetNDC(); 
    title1_4->SetTextSize(0.05);
    title1_4->SetTextAlign(22);  
    //title1_4->DrawLatex(0.5, 0.97, "Purity");

    c1->cd(5);  
    gPad->SetTopMargin(0.08);  
    TLatex* title1_5 = new TLatex();
    title1_5->SetNDC(); 
    title1_5->SetTextSize(0.05);
    title1_5->SetTextAlign(22);  
    //title1_5->DrawLatex(0.5, 0.97, "Bin Migration");

    c1->cd(6);  
    gPad->SetTopMargin(0.08);  
    TLatex* title1_6 = new TLatex();
    title1_6->SetNDC(); 
    title1_6->SetTextSize(0.05);
    title1_6->SetTextAlign(22);  
    //title1_6->DrawLatex(0.5, 0.97, "Corrected Acceptance");

    c1->cd(7);  
    gPad->SetTopMargin(0.08);  
    TLatex* title1_7 = new TLatex();
    title1_7->SetNDC(); 
    title1_7->SetTextSize(0.05);
    title1_7->SetTextAlign(22);  
    //title1_7->DrawLatex(0.5, 0.97, "Acceptance");

    c1->Print("./figures/plot_Q2.pdf");

// Q2 and y
    TCanvas* c111 = new TCanvas("c111","c111",1,1,1600,800);
    c111->Divide(3,2,0.01,0.01);
    c111->cd(1);
    gPad->SetLogy(1);
    //h_Q2_e->GetXaxis()->SetRangeUser(1,10);
	h_Q2_e->GetYaxis()->SetTitle("counts");
	h_Q2_e->GetXaxis()->SetTitle("Q^{2} [GeV/c]^{2}");
	h_Q2_e->SetLineColor(kBlack);
    h_Q2_e->Draw();
    //h_Q2_beforeCut->SetMarkerStyle(11);
    //h_Q2_beforeCut->SetMarkerColor(kOrange);
    //h_Q2_beforeCut->Draw("PEsame");
    h_Q2REC_e_EEMC->SetMarkerColor(kRed);
    h_Q2REC_e_EEMC->SetMarkerStyle(6);
    //h_Q2REC_e_EEMC->Draw("PEsame");
    h_Q2_afterCut->SetMarkerStyle(24);
    h_Q2_afterCut->SetMarkerColor(kBlue);
    h_Q2_afterCut->Draw("PEsame");
    //h_Q2REC_e_EEMC_cut->SetMarkerStyle(21);
    //h_Q2REC_e_EEMC_cut->SetMarkerColor(kMagenta);
    //h_Q2REC_e_EEMC_cut->Draw("PEsame");
    TLegend *w111 = new TLegend(0.42,0.7,0.6,0.85);
	w111->AddEntry(h_Q2_e, "MC", "L");
    //w1->AddEntry(h_Q2_beforeCut, "MC (before cut)", "P");
	//w111->AddEntry(h_Q2REC_e_EEMC, "RECO before cut", "P");
    w111->AddEntry(h_Q2_afterCut, "RECO", "P");
    //w111->AddEntry(h_Q2REC_e_EEMC_cut, "RECO", "P");
	w111->Draw("same");

    c111->cd(2);
    gPad->SetLogz(1);
    h_Q2_res->GetYaxis()->SetTitleOffset(1.5);
    h_Q2_res->GetYaxis()->SetTitle("(Q^{2}_{MC}-Q^{2}_{RECO})/Q^{2}_{MC}");
    h_Q2_res->GetXaxis()->SetTitle("Q^{2}_{MC} [GeV/c]^{2}");
    h_Q2_res->Draw("colzsame");

    c111->cd(3);
    gPad->SetLogz(1);
    h_Q2_response->GetYaxis()->SetTitleOffset(1.5);
    h_Q2_response->GetYaxis()->SetTitle("Q^{2}_{MC} [GeV/c]^{2}");
    h_Q2_response->GetXaxis()->SetTitle("Q^{2}_{RECO} [GeV/c]^{2}");
    h_Q2_response->Draw("colzsame");

    c111->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title1_111 = new TLatex();
    title1_111->SetNDC(); 
    title1_111->SetTextSize(0.05);
    title1_111->SetTextAlign(22);  
    title1_111->DrawLatex(0.5, 0.97, "Q^{2} Truth vs Reco");  

    c111->cd(2);  
    gPad->SetTopMargin(0.08);  
    TLatex* title1_211 = new TLatex();
    title1_211->SetNDC(); 
    title1_211->SetTextSize(0.05);
    title1_211->SetTextAlign(22);  
    title1_211->DrawLatex(0.5, 0.97, "Q^{2} Resolution");

    c111->cd(3);  
    gPad->SetTopMargin(0.08);  
    TLatex* title1_311 = new TLatex();
    title1_311->SetNDC(); 
    title1_311->SetTextSize(0.05);
    title1_311->SetTextAlign(22);  
    title1_311->DrawLatex(0.5, 0.97, "Q^{2} Response");

    c111->cd(4);
    gPad->SetLogy(1);
	h_y_e->GetYaxis()->SetTitle("counts");
	h_y_e->GetXaxis()->SetTitle("y");
	h_y_e->SetLineColor(kBlack);
    h_y_e->Draw();
    //h_y_beforeCut->SetMarkerStyle(11);
    //h_y_beforeCut->SetMarkerColor(kOrange);
    //h_y_beforeCut->Draw("Psame");
    h_y_afterCut->SetMarkerStyle(24);
    h_y_afterCut->SetMarkerColor(kBlue);
    h_y_afterCut->Draw("Psame");
    //h_yREC_e_EEMC_cut->SetMarkerColor(kBlue);
    //h_yREC_e_EEMC_cut->SetMarkerStyle(24);
    //h_yREC_e_EEMC_cut->Draw("PEsame");
    h_yREC_e_EEMC->SetMarkerStyle(6);
    h_yREC_e_EEMC->SetMarkerColor(kRed);
    //h_yREC_e_EEMC->Draw("PEsame");
    TLegend *w211 = new TLegend(0.45,0.7,0.6,0.85);
	w211->AddEntry(h_y_e, "MC", "L");
    //w2->AddEntry(h_y_beforeCut, "MC (before cut)", "P");
	//w211->AddEntry(h_yREC_e_EEMC, "RECO before cut", "P");
    w211->AddEntry(h_y_afterCut, "RECO", "P");
    //w211->AddEntry(h_yREC_e_EEMC, "RECO", "P");
	w211->Draw("same");

    c111->cd(5);
    gPad->SetLogz(1);
    h_y_res->GetXaxis()->SetRangeUser(0,0.4);
    h_y_res->GetYaxis()->SetTitleOffset(1.5);
	h_y_res->GetYaxis()->SetTitle("(y_{MC}-y_{RECO})/y_{MC}");
	h_y_res->GetXaxis()->SetTitle("y_{MC}");
    h_y_res->Draw("colzsame");

    c111->cd(6);
    gPad->SetLogz(1);
    h_y_response->GetXaxis()->SetRangeUser(0,0.4);
    h_y_response->GetYaxis()->SetRangeUser(0,0.4);
    h_y_response->GetYaxis()->SetTitle("y_{MC}");
    h_y_response->GetXaxis()->SetTitle("y_{RECO}");
    h_y_response->Draw("colzsame");

    c111->cd(4);  
    gPad->SetTopMargin(0.08); 
    TLatex* title1_411 = new TLatex();
    title1_411->SetNDC(); 
    title1_411->SetTextSize(0.05);
    title1_411->SetTextAlign(22);  
    title1_411->DrawLatex(0.5, 0.97, "y Truth vs Reco");  

    c111->cd(5);  
    gPad->SetTopMargin(0.08);  
    TLatex* title1_511 = new TLatex();
    title1_511->SetNDC(); 
    title1_511->SetTextSize(0.05);
    title1_511->SetTextAlign(22);  
    title1_511->DrawLatex(0.5, 0.97, "y Resolution");

    c111->cd(6);  
    gPad->SetTopMargin(0.08);  
    TLatex* title1_611 = new TLatex();
    title1_611->SetNDC(); 
    title1_611->SetTextSize(0.05);
    title1_611->SetTextAlign(22);  
    title1_611->DrawLatex(0.5, 0.97, "y Response");

    c111->Print("./figures/plot_Q2andy.pdf");

// y
    TCanvas* c2 = new TCanvas("c2","c2",1,1,1800,600);
    c2->Divide(3,1,0.01,0.01);
    c2->cd(1);
    gPad->SetLogy(1);
	h_y_e->GetYaxis()->SetTitle("counts");
	h_y_e->GetXaxis()->SetTitle("y");
	h_y_e->SetLineColor(kBlack);
    h_y_e->Draw();
    h_y_beforeCut->SetMarkerStyle(11);
    h_y_beforeCut->SetMarkerColor(kOrange);
    //h_y_beforeCut->Draw("Psame");
    h_y_afterCut->SetMarkerStyle(24);
    h_y_afterCut->SetMarkerColor(kBlue);
    h_y_afterCut->Draw("PEsame");
    h_yREC_e_EEMC_cut->SetMarkerColor(kBlue);
    h_yREC_e_EEMC_cut->SetMarkerStyle(24);
    //h_yREC_e_EEMC_cut->Draw("PEsame");
    h_yREC_e_EEMC->SetMarkerStyle(7);
    h_yREC_e_EEMC->SetMarkerColor(kRed);
    //h_yREC_e_EEMC->Draw("PEsame");
    TLegend *w2 = new TLegend(0.3,0.7,0.6,0.85);
	w2->AddEntry(h_y_e, " MC", "L");
    //w2->AddEntry(h_y_beforeCut, "MC (before cut)", "P");
	//w2->AddEntry(h_yREC_e_EEMC_cut, "RECO", "P");
    w2->AddEntry(h_y_afterCut, " RECO", "P");
    //w2->AddEntry(h_yREC_e_EEMC, "RECO before cut", "P");
	w2->Draw("same");

    c2->cd(2);
    gPad->SetLogz(1);
    h_y_res->GetXaxis()->SetRangeUser(0,0.8);
	h_y_res->GetYaxis()->SetTitle("(y_{MC}-y_{RECO})/y_{MC}");
	h_y_res->GetXaxis()->SetTitle("y_{MC}");
    h_y_res->Draw("colzsame");

    //c2->cd(2);
    //gPad->SetLogz(1);
    //h_y_res_cut->GetYaxis()->SetTitleOffset(1.5);
	//h_y_res_cut->GetYaxis()->SetTitle("(y_{MC}-y_{RECO})/y_{MC}");
	//h_y_res_cut->GetXaxis()->SetTitle("y_{MC}");
    //h_y_res_cut->Draw("colzsame");

    c2->cd(3);
    gPad->SetLogz(1);
    h_y_response->GetXaxis()->SetRangeUser(0,0.8);
    h_y_response->GetYaxis()->SetRangeUser(0,0.8);
    h_y_response->GetYaxis()->SetTitle("y_{MC}");
    h_y_response->GetXaxis()->SetTitle("y_{RECO}");
    h_y_response->Draw("colzsame");

    //c2->cd(3);
    //gPad->SetLogz(1);
    //h_y_response_cut->GetYaxis()->SetTitle("y_{MC}");
    //h_y_response_cut->GetXaxis()->SetTitle("y_{RECO}");
    //h_y_response_cut->Draw("colzsame)");

    //c2->cd(6);
    //gPad->SetLogy(1);
    //h_dyOvery_REC->GetYaxis()->SetTitle("counts");
    //h_dyOvery_REC->GetXaxis()->SetTitle("dy/y_{MC}");
    //h_dyOvery_REC->Draw("same");

    //c2->cd(8);
    //gPad->SetLogy(1);
    //h_dyOvery_REC_cut->GetYaxis()->SetTitle("counts");
    //h_dyOvery_REC_cut->GetXaxis()->SetTitle("dy/y_{MC}");
    //h_dyOvery_REC_cut->Draw("same");

    //c2->cd(4);
    //gPad->SetLogy(1);
    //h_y_purity->GetYaxis()->SetTitle("bins");
    //h_y_purity->Draw();

    //c2->cd(5);
    //gPad->SetLogz(1);
    //h_y_migration->Draw("colzsame");

    //c2->cd(6);
    //gPad->SetLogy(1);
    //h_y_e->GetYaxis()->SetTitle("counts");
	//h_y_e->GetXaxis()->SetTitle("y");
	//h_y_e->SetLineColor(kBlack);
    //h_y_e->Draw();
    //h_y_corrected->SetMarkerStyle(24);
    //h_y_corrected->SetMarkerColor(kBlue);
    //h_y_corrected->Draw("PEsame");
    //TLegend *w2_1 = new TLegend(0.42,0.7,0.6,0.85);
	//w2_1->AddEntry(h_y_e, "MC", "L");
    //w2_1->AddEntry(h_y_corrected, "Corrected", "P");
	//w2_1->Draw("same");

    //c2->cd(7);
    //gPad->SetLogy(1);
    //h_y_acceptance->Draw("same");

    c2->cd(4);  
    gPad->SetTopMargin(0.08);  
    TLatex* title2_4 = new TLatex();
    title2_4->SetNDC(); 
    title2_4->SetTextSize(0.05);
    title2_4->SetTextAlign(22);  
    //title2_4->DrawLatex(0.5, 0.97, "Purity");

    c2->cd(5);  
    gPad->SetTopMargin(0.08);  
    TLatex* title2_5 = new TLatex();
    title2_5->SetNDC(); 
    title2_5->SetTextSize(0.05);
    title2_5->SetTextAlign(22);  
    //title2_5->DrawLatex(0.5, 0.97, "Bin Migration");

    c2->cd(6);  
    gPad->SetTopMargin(0.08);  
    TLatex* title2_6 = new TLatex();
    title2_6->SetNDC(); 
    title2_6->SetTextSize(0.05);
    title2_6->SetTextAlign(22);  
    //title2_6->DrawLatex(0.5, 0.97, "Corrected Acceptance");

    c2->cd(7);  
    gPad->SetTopMargin(0.08);  
    TLatex* title2_7 = new TLatex();
    title2_7->SetNDC(); 
    title2_7->SetTextSize(0.05);
    title2_7->SetTextAlign(22);  
    //title2_7->DrawLatex(0.5, 0.97, "Acceptance");

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
    TCanvas* c3 = new TCanvas("c3","c3",1,1,1600,600);
    c3->Divide(3,1,0.01,0.01);
    c3->cd(1);
    gPad->SetLogy(1);
	h_energy_MC->GetYaxis()->SetTitle("counts");
	h_energy_MC->GetXaxis()->SetTitle("E [GeV]");
	h_energy_MC->SetLineColor(kBlack);
    //h_energy_MC->Draw();
    //h_energy_MC_after->SetMarkerStyle(11);
    h_energy_MC_after->SetLineColor(kBlack);
    h_energy_MC_after->Draw();
    //h_energy_REC_EEMC_cut->SetMarkerColor(kBlue);
    //h_energy_REC_EEMC_cut->SetMarkerStyle(24);
    //h_energy_REC_EEMC_cut->Draw("PEsame");
    h_energy_REC_EEMC_after->SetMarkerStyle(24);
    h_energy_REC_EEMC_after->SetMarkerColor(kBlue);
    h_energy_REC_EEMC_after->Draw("PEsame");
    h_energy_REC_EEMC->SetMarkerColor(kRed);
    h_energy_REC_EEMC->SetMarkerStyle(30);
    //h_energy_REC_EEMC->Draw("PEsame");
    TLegend *w3 = new TLegend(0.56,0.2,0.65,0.4);
	//w3->AddEntry(h_energy_MC, "MC", "L");
    w3->AddEntry(h_energy_MC_after, " MC", "L");
	//w3->AddEntry(h_energy_REC_EEMC, "RECO bc", "P");
    w3->AddEntry(h_energy_REC_EEMC_after, " RECO", "P");
    //w3->AddEntry(h_energy_REC_EEMC_cut, "RECO", "P");
	w3->Draw("same");

    c3->cd(2);
    gPad->SetLogz(1);
    h_energy_res_EEMC_after->GetXaxis()->SetRangeUser(1,10);
    h_energy_res_EEMC_after->GetYaxis()->SetRangeUser(-0.2,0.2);
    h_energy_res_EEMC_after->GetYaxis()->SetTitle("(E_{MC}-E_{RECO})/E_{MC}");
    h_energy_res_EEMC_after->GetXaxis()->SetTitle("E_{MC} [GeV]");
    h_energy_res_EEMC_after->Draw("colzsame");

    //c3->cd(2);
    //gPad->SetLogz(1);
    //h_energy_res_EEMC_cut->GetXaxis()->SetRangeUser(0,10);
    //h_energy_res_EEMC_cut->GetYaxis()->SetTitle("(E_{MC}-E_{RECO})/E_{MC}");
    //h_energy_res_EEMC_cut->GetXaxis()->SetTitle("E_{MC} [GeV]");
    //h_energy_res_EEMC_cut->Draw("same");

    c3->cd(3);
    gPad->SetLogz(1);
    h_energy_response_EEMC_after->GetXaxis()->SetRangeUser(0,12);
    h_energy_response_EEMC_after->GetYaxis()->SetRangeUser(0,12);
    h_energy_response_EEMC_after->GetYaxis()->SetTitle("E_{MC} [GeV]");
    h_energy_response_EEMC_after->GetXaxis()->SetTitle("E_{RECO} [GeV]");
    h_energy_response_EEMC_after->Draw("colzsame");

    //c3->cd(3);
    //gPad->SetLogz(1);
    //h_energy_response_EEMC_cut->GetXaxis()->SetRangeUser(0,10);
    //h_energy_response_EEMC_cut->GetYaxis()->SetRangeUser(0,10);
    //h_energy_response_EEMC_cut->GetYaxis()->SetTitle("E_{MC} [GeV]");
    //h_energy_response_EEMC_cut->GetXaxis()->SetTitle("E_{RECO} [GeV]");
    //h_energy_response_EEMC_cut->Draw("same");

    //c3->cd(4);
    //gPad->SetLogy(1);
    //h_energy_purity->GetYaxis()->SetTitle("bins");
    //h_energy_purity->Draw();

    //c3->cd(5);
    //gPad->SetLogz(1);
    //h_energy_migration->Draw("colzsame");

    //c3->cd(6);
    //gPad->SetLogy(1);
    h_energy_MC->GetYaxis()->SetTitle("counts");
	h_energy_MC->GetXaxis()->SetTitle("E [GeV]");
	h_energy_MC->SetLineColor(kBlack);
    //h_energy_MC->Draw();
    //h_energy_corrected->SetMarkerStyle(24);
    //h_energy_corrected->SetMarkerColor(kBlue);
    ////h_energy_corrected->Draw("PEsame");
    //TLegend *w3_1 = new TLegend(0.57,0.7,0.75,0.85);
	//w3_1->AddEntry(h_energy_MC, "MC", "L");
    //w3_1->AddEntry(h_energy_corrected, "Corrected", "P");
	//w3_1->Draw("same");

    //c3->cd(7);
    //gPad->SetLogy(1);
    //h_energy_acceptance->Draw("same");

    //c3->cd(8);
    //gPad->SetLogz(1);
    //h_energy_res_EEMC_after->Draw("colzsame");

    //c3->cd(9);
    //gPad->SetLogz(1);
    //h_energy_response_EEMC_after->Draw("colzsame");

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

    c3->cd(4);  
    gPad->SetTopMargin(0.08);  
    TLatex* title3_4 = new TLatex();
    title3_4->SetNDC(); 
    title3_4->SetTextSize(0.05);
    title3_4->SetTextAlign(22);  
    //title3_4->DrawLatex(0.5, 0.97, "Purity");  

    c3->cd(5);  
    gPad->SetTopMargin(0.08);  
    TLatex* title3_5 = new TLatex();
    title3_5->SetNDC(); 
    title3_5->SetTextSize(0.05);
    title3_5->SetTextAlign(22);  
    //title3_5->DrawLatex(0.5, 0.97, "Bin Migration");  

    c3->cd(6);  
    gPad->SetTopMargin(0.08);  
    TLatex* title3_6 = new TLatex();
    title3_6->SetNDC(); 
    title3_6->SetTextSize(0.05);
    title3_6->SetTextAlign(22);  
    //title3_6->DrawLatex(0.5, 0.97, "Corrected Acceptance");  

    c3->cd(7);  
    gPad->SetTopMargin(0.08);  
    TLatex* title3_7 = new TLatex();
    title3_7->SetNDC(); 
    title3_7->SetTextSize(0.05);
    title3_7->SetTextAlign(22);  
    //title3_7->DrawLatex(0.5, 0.97, "Acceptance");  

    c3->Print("./figures/plot_energy.pdf");

// e' p
    TCanvas* c32 = new TCanvas("c32","c32",1,1,1200,800);
    c32->Divide(3,3,0.01,0.01);
    c32->cd(1);
    gPad->SetLogy(1);
	h_e_pt_MC->GetYaxis()->SetTitle("counts");
	h_e_pt_MC->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	h_e_pt_MC->SetLineColor(kBlack);
    //h_e_pt_MC->Draw();
    //h_e_pt_MC_after->SetMarkerColor(kBlue);
    //h_e_pt_MC_after->SetMarkerStyle(24);
    h_e_pt_MC_after->SetLineColor(kBlack);
    h_e_pt_MC_after->Draw();
    //h_e_pt_REC_EEMC_cut->SetMarkerColor(kRed);
    //h_e_pt_REC_EEMC_cut->SetMarkerStyle(29);
    //h_e_pt_REC_EEMC_cut->Draw("PEsame");
    h_e_pt_REC_EEMC->SetMarkerStyle(24);
    h_e_pt_REC_EEMC->SetMarkerColor(kBlue);
    h_e_pt_REC_EEMC->Draw("PEsame");
    h_e_pt_REC_trk->SetMarkerStyle(7);
    h_e_pt_REC_trk->SetMarkerColor(kOrange);
    h_e_pt_REC_trk->Draw("PEsame");
    //h_e_pt_REC_trk_cut->SetMarkerStyle(13);
    //h_e_pt_REC_trk_cut->SetMarkerColor(kOrange);
    //h_e_pt_REC_trk_cut->Draw("PEsame");
    TLegend *w32 = new TLegend(0.6,0.7,0.75,0.85);
	//w32->AddEntry(h_e_pt_MC, "MC", "L");
	w32->AddEntry(h_e_pt_MC_after, " MC", "P");
    //w32->AddEntry(h_e_pt_REC_EEMC_cut, "reco (inside threshold)", "P");
    w32->AddEntry(h_e_pt_REC_EEMC, " RECO cal", "P");
    w32->AddEntry(h_e_pt_REC_trk, "RECO trk", "P");
    //w32->AddEntry(h_e_pt_REC_trk_cut, "RECO", "P");
	w32->Draw("same");

    c32->cd(1);  
    gPad->SetTopMargin(0.08);  
    TLatex* title32_1 = new TLatex();
    title32_1->SetNDC(); 
    title32_1->SetTextSize(0.05);
    title32_1->SetTextAlign(22);  
    title32_1->DrawLatex(0.5, 0.97, "p_{e,T} Compare Reco with MC"); 

    c32->cd(2);
    gPad->SetLogz(1);
    h_e_pt_res->GetYaxis()->SetTitleOffset(1.5);
    h_e_pt_res->GetXaxis()->SetRangeUser(0,4);
    h_e_pt_res->GetYaxis()->SetRangeUser(-0.2,0.2);
    h_e_pt_res->Draw("colzsame");

    c32->cd(2);  
    gPad->SetTopMargin(0.08);  
    TLatex* title32_2 = new TLatex();
    title32_2->SetNDC(); 
    title32_2->SetTextSize(0.05);
    title32_2->SetTextAlign(22);  
    title32_2->DrawLatex(0.5, 0.97, "p_{e,T} Resolution"); 

    c32->cd(3);
    gPad->SetLogy(0);  
    TProfile* prof_res_e_pt = h_e_pt_res->ProfileX();  
    prof_res_e_pt->SetLineColor(kBlue);
    prof_res_e_pt->GetXaxis()->SetRangeUser(0,4);
    prof_res_e_pt->SetMarkerStyle(24);
    prof_res_e_pt->SetMarkerColor(kBlue);
    prof_res_e_pt->Draw();

    c32->cd(3);  
    gPad->SetTopMargin(0.08);  
    TLatex* title32_3 = new TLatex();
    title32_3->SetNDC(); 
    title32_3->SetTextSize(0.05);
    title32_3->SetTextAlign(22);  
    title32_3->DrawLatex(0.5, 0.97, "p_{e,T} Resolution ProfileX");

    //c32->cd(3);
    //gPad->SetLogz(1);
    //h_e_pt_res_cut->Draw("colzsame");

    c32->cd(4);
    gPad->SetLogy(1);
	h_e_pz_MC->GetYaxis()->SetTitle("counts");
	h_e_pz_MC->GetXaxis()->SetTitle("p_{z} [GeV/c]");
	h_e_pz_MC->SetLineColor(kBlack);
    //h_e_pz_MC->Draw();
    //h_e_pz_MC_after->SetMarkerColor(kBlue);
    //h_e_pz_MC_after->SetMarkerStyle(24);
    h_e_pz_MC_after->SetLineColor(kBlack);
    h_e_pz_MC_after->Draw();
    //h_e_pz_REC_EEMC_cut->SetMarkerColor(kRed);
    //h_e_pz_REC_EEMC_cut->SetMarkerStyle(29);
    //h_e_pz_REC_EEMC_cut->Draw("PEsame");
    h_e_pz_REC_EEMC->SetMarkerStyle(24);
    h_e_pz_REC_EEMC->SetMarkerColor(kBlue);
    h_e_pz_REC_EEMC->Draw("PEsame");
    h_e_pz_REC_trk->SetMarkerStyle(7);
    h_e_pz_REC_trk->SetMarkerColor(kOrange);
    h_e_pz_REC_trk->Draw("PEsame");
    //h_e_pz_REC_trk_cut->SetMarkerStyle(13);
    //h_e_pz_REC_trk_cut->SetMarkerColor(kBlack);
    //h_e_pz_REC_trk_cut->Draw("same");
    w32->Draw("same");

    c32->cd(5);
    gPad->SetLogz(1);
    h_e_pz_res->Draw("colzsame");

    c32->cd(6);
    TProfile* prof_res_e_pz = h_e_pz_res->ProfileX();  
    prof_res_e_pz->SetLineColor(kBlue);
    prof_res_e_pz->SetMarkerStyle(24);
    prof_res_e_pz->SetMarkerColor(kBlue);
    prof_res_e_pz->Draw();

    //c32->cd(4);
    //gPad->SetLogz(1);
    //h_e_pz_trk_res->Draw("colzsame");

    //c32->cd(6);
    //gPad->SetLogz(1);
    //h_e_pz_response->Draw("colzsame");

    //c32->cd(8);
    //gPad->SetLogz(1);
    //h_e_pz_response_cut->Draw("colzsame");

    c32->cd(7);
    gPad->SetLogy(1);
	//h_e_p_MC->GetYaxis()->SetTitle("counts");
	//h_e_p_MC->GetXaxis()->SetTitle("|p| [GeV/c]");
	//h_e_p_MC->SetLineColor(kBlack);
    //h_e_p_MC->Draw();
    //h_e_p_MC_after->SetMarkerColor(kBlue);
    //h_e_p_MC_after->SetMarkerStyle(24);
    h_e_p_MC_after->SetMarkerColor(kBlack);
    h_e_p_MC_after->Draw();
    //h_e_p_REC_EEMC_cut->SetMarkerColor(kRed);
    //h_e_p_REC_EEMC_cut->SetMarkerStyle(29);
    //h_e_p_REC_EEMC_cut->Draw("PEsame");
    h_e_p_REC_EEMC->SetMarkerStyle(24);
    h_e_p_REC_EEMC->SetMarkerColor(kBlue);
    h_e_p_REC_EEMC->Draw("PEsame");
    h_e_p_REC_trk->SetMarkerStyle(7);
    h_e_p_REC_trk->SetMarkerColor(kOrange);
    h_e_p_REC_trk->Draw("PEsame");
    //h_e_p_REC_trk_cut->SetMarkerStyle(13);
    //h_e_p_REC_trk_cut->SetMarkerColor(kBlack);
    //h_e_p_REC_trk_cut->Draw("same");
    w32->Draw("same"); 

    c32->cd(4);  
    gPad->SetTopMargin(0.08);  
    TLatex* title32_4 = new TLatex();
    title32_4->SetNDC(); 
    title32_4->SetTextSize(0.05);
    title32_4->SetTextAlign(22);  
    //title32_4->DrawLatex(0.5, 0.97, "|p|_{e} Compare Reco with MC"); 

    c32->cd(8);
    gPad->SetLogz(1);
    h_e_p_res->Draw("colzsame");

    c32->cd(9);
    gPad->SetLogy(0);  
    TProfile* prof_res_e_p = h_e_p_res->ProfileX();  
    prof_res_e_p->SetLineColor(kBlue);
    prof_res_e_p->GetXaxis()->SetRangeUser(0,4);
    prof_res_e_p->SetMarkerStyle(24);
    prof_res_e_p->SetMarkerColor(kBlue);
    prof_res_e_p->Draw();

    c32->cd(5);  
    gPad->SetTopMargin(0.08);  
    TLatex* title32_5 = new TLatex();
    title32_5->SetNDC(); 
    title32_5->SetTextSize(0.05);
    title32_5->SetTextAlign(22);  
    //title32_5->DrawLatex(0.5, 0.97, "|p|_{e} Resolution");

    c32->cd(6);  
    gPad->SetTopMargin(0.08);  
    TLatex* title32_6 = new TLatex();
    title32_6->SetNDC(); 
    title32_6->SetTextSize(0.05);
    title32_6->SetTextAlign(22);  
    //title32_6->DrawLatex(0.5, 0.97, "|p|_{e} Resolution ProfileX"); 

    cout << "|p|_{e} resolution: " << prof_res_e_p->GetMean() << endl;

    //c32->cd(11);
    //gPad->SetLogz(1);
    //h_e_p_res_cut->Draw("colzsame");
   
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

    c31->cd(2);
    //gPad->SetLogz(1);
    //h_EvsP_REC->GetXaxis()->SetTitle("|p|_{trk} [GeV]");
    //h_EvsP_REC->GetYaxis()->SetTitle("E_{EMCal} [GeV]");
    //h_EvsP_REC->Draw("same");

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
    c4->Divide(2,1,0.01,0.01);
    c4->cd(1);
    gPad->SetLogy(1);
	h_eta_MC->GetYaxis()->SetTitle("counts");
	h_eta_MC->GetXaxis()->SetTitle("#eta");
	h_eta_MC->SetLineColor(kBlack);
    //h_eta_MC->Draw();
    //h_eta_MC_after->SetMarkerStyle(11);
    h_eta_MC_after->SetLineColor(kBlack);
    h_eta_MC_after->Draw();
    h_eta_REC_EEMC->SetMarkerColor(kBlue);
    h_eta_REC_EEMC->SetMarkerStyle(24);
    //h_eta_REC_EEMC->Draw("PEsame");
    //h_eta_REC_EEMC_cut->SetMarkerStyle(29);
    //h_eta_REC_EEMC_cut->SetMarkerColor(kRed);
    //h_eta_REC_EEMC_cut->Draw("PEsame");
    h_eta_REC_EEMC_after->SetMarkerStyle(24);
    h_eta_REC_EEMC_after->SetMarkerColor(kBlue);
    h_eta_REC_EEMC_after->Draw("PEsame");
    TLegend *w4 = new TLegend(0.62,0.6,0.7,0.7);
	//w4->AddEntry(h_eta_MC, "MC (before cut)", "L");
    w4->AddEntry(h_eta_MC_after, " MC", "P");
    //w4->AddEntry(h_eta_REC_EEMC, "RECO", "P");
    //w4->AddEntry(h_eta_REC_EEMC_cut, "reco (inside threshold)", "P");
    w4->AddEntry(h_eta_REC_EEMC_after, " RECO", "P");
	w4->Draw("same");

    //c4->cd(2);
    //gPad->SetLogy(1);
	//h_eta_diff->GetYaxis()->SetTitle("counts");
	//h_eta_diff->GetXaxis()->SetTitle("#eta_{RECO}-#eta_{MC}");
	//h_eta_diff->SetLineColor(kBlack);
    //h_eta_diff->Draw();

    //c4->cd(3);
    //gPad->SetLogy(1);
	//h_eta_diff_cut->GetYaxis()->SetTitle("counts");
	//h_eta_diff_cut->GetXaxis()->SetTitle("#eta_{RECO}-#eta_{MC}");
	//h_eta_diff_cut->SetLineColor(kBlack);
    //h_eta_diff_cut->Draw();

    c4->cd(2);
    //gPad->SetLogy(1);
	h_eta_diff_after->GetYaxis()->SetTitle("counts");
	h_eta_diff_after->GetXaxis()->SetTitle("#eta_{RECO}-#eta_{MC}");
	h_eta_diff_after->SetLineColor(kBlack);
    h_eta_diff_after->Draw();
       
    c4->Print("./figures/plot_eta.pdf");


// phi
    TCanvas* c5 = new TCanvas("c5","c5",1,1,1600,600);
    c5->Divide(4,1,0.01,0.01);
    c5->cd(1);
    gPad->SetLogy(1);
	h_phi_MC->GetYaxis()->SetTitle("counts");
	h_phi_MC->GetXaxis()->SetTitle("#phi");
	h_phi_MC->SetLineColor(kBlack);
    //h_phi_MC->Draw();
    h_phi_REC_EEMC->SetMarkerColor(kBlue);
    h_phi_REC_EEMC->SetMarkerStyle(24);
    //h_phi_REC_EEMC->Draw("PEsame");
    //h_phi_MC_after->SetMarkerStyle(11);
    h_phi_MC_after->SetLineColor(kBlack);
    h_phi_MC_after->Draw();
    //h_phi_REC_EEMC_cut->SetMarkerStyle(29);
    //h_phi_REC_EEMC_cut->SetMarkerColor(kRed);
    //h_phi_REC_EEMC_cut->Draw("PEsame");
    h_phi_REC_EEMC_after->SetMarkerColor(kBlue);
    h_phi_REC_EEMC_after->SetMarkerStyle(24);
    h_phi_REC_EEMC_after->Draw("PEsame");
    TLegend *w5 = new TLegend(0.22,0.7,0.4,0.85);
	//w5->AddEntry(h_phi_MC, "MC (before cut)", "L");
	//w5->AddEntry(h_phi_REC_EEMC, "reco (before cut)", "P");
    w5->AddEntry(h_phi_MC_after, " MC", "P");
    //w5->AddEntry(h_phi_REC_EEMC_cut, "reco (inside threshold)", "P");
    w5->AddEntry(h_phi_REC_EEMC_after, " RECO", "P");
	w5->Draw("same");

    c5->cd(2);
    //gPad->SetLogy(1);
	h_phi_diff_after->GetYaxis()->SetTitle("counts");
	h_phi_diff_after->GetXaxis()->SetTitle("#phi_{RECO}-#phi_{MC}");
	h_phi_diff_after->SetLineColor(kBlack);
    h_phi_diff_after->Draw();

    c5->cd(3);
    gPad->SetLogz(1);
    h_phi_resolution->Draw("colzsame");

    //c5->cd(3);
    //gPad->SetLogz(1);
    //h_phi_res_diff2d->Draw("colzsame");

    c5->cd(4);
    gPad->SetLogz(1);
    h_phi_response->Draw("PEsame");
    
    c5->cd(2);  
    gPad->SetTopMargin(0.08);  
    TLatex* title5_3 = new TLatex();
    title5_3->SetNDC(); 
    title5_3->SetTextSize(0.05);
    title5_3->SetTextAlign(22);  
    //title5_3->DrawLatex(0.5, 0.97, "#phi Resolution");  

    c5->cd(3);  
    gPad->SetTopMargin(0.08);  
    TLatex* title5_4 = new TLatex();
    title5_4->SetNDC(); 
    title5_4->SetTextSize(0.05);
    title5_4->SetTextAlign(22);  
    //title5_4->DrawLatex(0.5, 0.97, "#phi Resolution");  

    c5->cd(5);  
    gPad->SetTopMargin(0.08);  
    TLatex* title5_5 = new TLatex();
    title5_5->SetNDC(); 
    title5_5->SetTextSize(0.05);
    title5_5->SetTextAlign(22);  
    //title5_5->DrawLatex(0.5, 0.97, "#phi Resolution");  

    c5->Print("./figures/plot_phi.pdf");


// theta
    TCanvas* c6 = new TCanvas("c6","c6",1,1,1200,800);
    c6->Divide(3,1,0.01,0.01);
    c6->cd(1);
    gPad->SetLogy(1);
	h_theta_MC_after->GetYaxis()->SetTitle("counts");
	h_theta_MC->GetXaxis()->SetTitle("#phi");
	h_theta_MC->SetLineColor(kBlack);
    //h_theta_MC->Draw();
    //h_theta_MC_after->GetXaxis()->SetRangeUser(-1,3.5);
    //h_theta_MC_after->GetYaxis()->SetTitle("counts");
	//h_theta_MC_after->GetXaxis()->SetTitle("#theta");
    h_theta_MC_after->SetLineColor(kBlack);
    //h_theta_MC_after->SetMarkerStyle(11);
    h_theta_MC_after->Draw();
    //h_theta_REC_EEMC_cut->SetMarkerColor(kBlue);
    //h_theta_REC_EEMC_cut->SetMarkerStyle(24);
    //h_theta_REC_EEMC_cut->Draw("PEsame");
    //h_theta_REC_EEMC->SetMarkerColor(kRed);
    //h_theta_REC_EEMC->SetMarkerStyle(6);
    //h_theta_REC_EEMC->Draw("PEsame");
    h_theta_REC_EEMC_after->SetMarkerColor(kBlue);
    h_theta_REC_EEMC_after->SetMarkerStyle(24);
    h_theta_REC_EEMC_after->Draw("PEsame");
    TLegend *w6 = new TLegend(0.22,0.7,0.4,0.85);
	//w6->AddEntry(h_theta_MC, "MC (before cut)", "L");
    w6->AddEntry(h_theta_MC_after, "MC", "L");
	//w6->AddEntry(h_theta_REC_EEMC_cut, "RECO", "P");
    //w6->AddEntry(h_theta_REC_EEMC, "RECO before cut", "P");
    w6->AddEntry(h_theta_REC_EEMC_after, "RECO after cut", "P");
	w6->Draw("same");

    cout << "theta MC = " << h_theta_MC_after->GetMean() << " theta RECO = " << h_theta_REC_EEMC_after->GetMean() << endl;

    c6->cd(2);
    gPad->SetLogz(1);
    h_theta_response_EEMC_after->GetYaxis()->SetTitleOffset(1.5);
    h_theta_response_EEMC_after->GetXaxis()->SetRangeUser(2.7,3.1);
    h_theta_response_EEMC_after->GetYaxis()->SetRangeUser(2.7,3.1);
    h_theta_response_EEMC_after->GetYaxis()->SetTitle("#theta_{MC} [GeV]");
    h_theta_response_EEMC_after->GetXaxis()->SetTitle("#theta_{RECO} [GeV]");
    h_theta_response_EEMC_after->Draw("same");

    c6->cd(3);
    gPad->SetLogz(1);
    h_theta_resolution->GetXaxis()->SetRangeUser(2.7,3.1);
    h_theta_resolution->GetYaxis()->SetRangeUser(-0.02,0.02);
    h_theta_resolution->Draw("colzsame");


    //c6->cd(4);
    //gPad->SetLogx(1);
	//h_theta_diff_after->GetYaxis()->SetTitle("counts");
	//h_theta_diff_after->GetXaxis()->SetTitle("#theta_{RECO}-#theta_{MC}");
	//h_theta_diff_after->SetLineColor(kBlack);
    //h_theta_diff_after->Draw();

    //c6->cd(5);
    //gPad->SetLogz(1);
    //h_theta_res_diff2d->Draw("colzsame");

    //c6->cd(6);
    //gPad->SetLogy(1);
    //h_theta_res_diff->Draw("PEsame");

    c6->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title6 = new TLatex();
    title6->SetNDC(); 
    title6->SetTextSize(0.05);
    title6->SetTextAlign(22);  
    title6->DrawLatex(0.5, 0.97, "#theta Truth vs Reco");  

    c6->cd(4);  
    gPad->SetTopMargin(0.08);  
    TLatex* title6_2 = new TLatex();
    title6_2->SetNDC(); 
    title6_2->SetTextSize(0.05);
    title6_2->SetTextAlign(22);  
    //title6_2->DrawLatex(0.5, 0.97, "#theta Difference");

    c6->cd(2);  
    gPad->SetTopMargin(0.08);  
    TLatex* title6_3 = new TLatex();
    title6_3->SetNDC(); 
    title6_3->SetTextSize(0.05);
    title6_3->SetTextAlign(22);  
    title6_3->DrawLatex(0.5, 0.97, "#theta Response");  

    c6->cd(3);  
    gPad->SetTopMargin(0.08);  
    TLatex* title6_4 = new TLatex();
    title6_4->SetNDC(); 
    title6_4->SetTextSize(0.05);
    title6_4->SetTextAlign(22);  
    title6_4->DrawLatex(0.5, 0.97, "#theta Resolution");  

    c6->cd(5);  
    gPad->SetTopMargin(0.08);  
    TLatex* title6_5 = new TLatex();
    title6_5->SetNDC(); 
    title6_5->SetTextSize(0.05);
    title6_5->SetTextAlign(22);  
    //title6_5->DrawLatex(0.5, 0.97, "#theta Resolution");  

    c6->cd(6);  
    gPad->SetTopMargin(0.08);  
    TLatex* title6_6 = new TLatex();
    title6_6->SetNDC(); 
    title6_6->SetTextSize(0.05);
    title6_6->SetTextAlign(22);  
    //title6_6->DrawLatex(0.5, 0.97, "#theta Resolution");  
	
    c6->Print("./figures/plot_theta.pdf");


// VM (E-pz)
    TCanvas* c7 = new TCanvas("c7","c7",1,1,1600,600);
    c7->Divide(3,1,0.01,0.01);
    c7->cd(1);
    gPad->SetLogy(1);
	h_VM_Epz_MC->GetYaxis()->SetTitle("counts");
	h_VM_Epz_MC->GetXaxis()->SetTitle("(E_{VM}-p_{z,VM}) [GeV]");
	h_VM_Epz_MC->SetLineColor(kBlack);
    h_VM_Epz_MC->Draw();
    h_VM_Epz_REC->SetMarkerStyle(24);
    h_VM_Epz_REC->SetMarkerColor(kBlue);
    h_VM_Epz_REC->Draw("PEsame");
    //h_VM_Epz_MC_before->SetMarkerStyle(21);
    //h_VM_Epz_MC_before->SetMarkerColor(kRed);
    //h_VM_Epz_MC_before->Draw("PEsame");
    //h_VM_Epz_REC_before->SetMarkerStyle(6);
    //h_VM_Epz_REC_before->SetMarkerColor(kRed);
    //h_VM_Epz_REC_before->Draw("PEsame");
    //h_VM_Epz_REC_cut->SetMarkerColor(kBlue);
    //h_VM_Epz_REC_cut->SetMarkerStyle(24);
    //h_VM_Epz_REC_cut->Draw("PEsame");
    TLegend *w7 = new TLegend(0.42,0.75,0.6,0.85);
	w7->AddEntry(h_VM_Epz_MC, "MC", "L");
	//w7->AddEntry(h_VM_Epz_REC_cut, "RECO", "P");
    w7->AddEntry(h_VM_Epz_REC, "RECO after cut", "P");
    //w7->AddEntry(h_VM_Epz_REC_before, "RECO before cut", "P");
	w7->Draw("same");

    c7->cd(2);
    gPad->SetLogz(1);
    h_VM_Epz_response->GetXaxis()->SetRangeUser(0,20);
    h_VM_Epz_response->GetYaxis()->SetRangeUser(0,20);
	h_VM_Epz_response->GetYaxis()->SetTitle("(E_{VM,MC}-p_{z,VM,MC}) [GeV]");
	h_VM_Epz_response->GetXaxis()->SetTitle("(E_{VM,RECO}-p_{z,VM,RECO}) [GeV]");
    h_VM_Epz_response->Draw("colzsame");

    c7->cd(3);
    h_VM_Epz_res->Draw("colzsame");

    //c7->cd(2);
    //gPad->SetLogz(1);
	//h_VM_Epz_response_cut->GetYaxis()->SetTitle("(E_{VM,MC}-p_{z,VM,MC}) [GeV]");
	//h_VM_Epz_response_cut->GetXaxis()->SetTitle("(E_{VM,RECO}-p_{z,VM,RECO}) [GeV]");
    //h_VM_Epz_response_cut->Draw("colzsame");

    //c7->cd(3);
    //gPad->SetLogy(1);
    //h_VM_Epz_purity->GetYaxis()->SetTitle("bins");
    //h_VM_Epz_purity->Draw("same");

    //c7->cd(4);
    //gPad->SetLogz(1);
    //h_VM_Epz_migration->Draw("colzsame");

    //c7->cd(5);
    //gPad->SetLogy(1);
    //h_VM_Epz_MC->GetYaxis()->SetTitle("counts");
	//h_VM_Epz_MC->GetXaxis()->SetTitle("(E_{VM}-p_{z,VM}) [GeV]");
	//h_VM_Epz_MC->SetLineColor(kBlack);
    //h_VM_Epz_MC->Draw();
    //h_VM_Epz_corrected->SetMarkerStyle(24);
    //h_VM_Epz_corrected->SetMarkerColor(kBlue);
    //h_VM_Epz_corrected->Draw("PEsame");
    //TLegend *w7_1 = new TLegend(0.57,0.7,0.75,0.85);
	//w7_1->AddEntry(h_VM_Epz_MC, "MC", "L");
    //w7_1->AddEntry(h_VM_Epz_corrected, "Corrected", "P");
	//w7_1->Draw("same");

    //c7->cd(6);
    //gPad->SetLogy(1);
    //h_VM_Epz_acceptance->Draw("same");

    //c7->cd(3);  
    //gPad->SetTopMargin(0.08);  
    //TLatex* title7_3 = new TLatex();
    //title7_3->SetNDC(); 
    //title7_3->SetTextSize(0.05);
    //title7_3->SetTextAlign(22);  
    //title7_3->DrawLatex(0.5, 0.97, "Purity");  

    c7->cd(4);  
    gPad->SetTopMargin(0.08);  
    TLatex* title7_4 = new TLatex();
    title7_4->SetNDC(); 
    title7_4->SetTextSize(0.05);
    title7_4->SetTextAlign(22);  
    //title7_4->DrawLatex(0.5, 0.97, "Bin Migration"); 

    c7->cd(5);  
    gPad->SetTopMargin(0.08);  
    TLatex* title7_5 = new TLatex();
    title7_5->SetNDC(); 
    title7_5->SetTextSize(0.05);
    title7_5->SetTextAlign(22);  
    //title7_5->DrawLatex(0.5, 0.97, "Corrected Acceptance"); 

    c7->cd(6);  
    gPad->SetTopMargin(0.08);  
    TLatex* title7_6 = new TLatex();
    title7_6->SetNDC(); 
    title7_6->SetTextSize(0.05);
    title7_6->SetTextAlign(22);  
    //title7_6->DrawLatex(0.5, 0.97, "Acceptance"); 

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

    c7->cd(3);  
    gPad->SetTopMargin(0.08);  
    TLatex* title7_3 = new TLatex();
    title7_3->SetNDC(); 
    title7_3->SetTextSize(0.05);
    title7_3->SetTextAlign(22);  
    title7_3->DrawLatex(0.5, 0.97, "E-p_{z} Resolution (VM)");  


    c7->Print("./figures/plot_VM_Epz.pdf");

// VM all plots
    TCanvas* c777 = new TCanvas("c777","c777",1,1,1800,600);
    c777->Divide(3,1,0.01,0.01);
    c777->cd(1);
    gPad->SetLogy(1);
    //h_VM_Epz_MC->GetXaxis()->SetRangeUser(0,20);
    //h_VM_Epz_MC->GetYaxis()->SetRangeUser(1e1,1e7);
    h_VM_Epz_MC->GetYaxis()->SetTitle("counts");
	h_VM_Epz_MC->GetXaxis()->SetTitle("(E_{VM}-p_{z,VM}) [GeV]");
	h_VM_Epz_MC->SetLineColor(kBlack);
    h_VM_Epz_MC->Draw();
    //h_VM_Epz_MC_before->SetMarkerStyle(11);
    //h_VM_Epz_MC_before->SetMarkerColor(kMagenta);
    //h_VM_Epz_MC_before->Draw("PEsame");
    h_VM_Epz_REC->SetMarkerStyle(24);
    h_VM_Epz_REC->SetMarkerColor(kBlue);
    h_VM_Epz_REC->Draw("PEsame");
    h_Epz_afterCut->SetMarkerStyle(30);
    h_Epz_afterCut->SetMarkerColor(kOrange);
    //h_Epz_afterCut->Draw("PEsame");
    h_VM_Epz_REC_cut->SetMarkerColor(kBlue);
    h_VM_Epz_REC_cut->SetMarkerStyle(24);
    //h_VM_Epz_REC_cut->Draw("PEsame");
    TLegend *w777 = new TLegend(0.42,0.75,0.6,0.85);
	w777->AddEntry(h_VM_Epz_MC, " MC", "L");
    //w777->AddEntry(h_VM_Epz_MC_before, "MC before cut", "L");
	//w777->AddEntry(h_VM_Epz_REC_cut, " RECO", "P");
    w777->AddEntry(h_VM_Epz_REC, " RECO", "P");
    //w777->AddEntry(h_Epz_afterCut, "RECO after cut", "P");
	w777->Draw("same");

    //c777->cd(2);
    //gPad->SetLogz(1);
    //h_VM_Epz_response_cut->GetXaxis()->SetRangeUser(0,20);
    //h_VM_Epz_response_cut->GetYaxis()->SetRangeUser(0,20);
	//h_VM_Epz_response_cut->GetYaxis()->SetTitle("(E_{VM,MC}-p_{z,VM,MC}) [GeV]");
	//h_VM_Epz_response_cut->GetXaxis()->SetTitle("(E_{VM,RECO}-p_{z,VM,RECO}) [GeV]");
    //h_VM_Epz_response_cut->Draw("colzsame");

    c777->cd(2);
    gPad->SetLogz(1);
    //h_VM_Epz_response->GetXaxis()->SetRangeUser(0,20);
    //h_VM_Epz_response->GetYaxis()->SetRangeUser(0,20);
	h_VM_Epz_response->GetYaxis()->SetTitle("(E_{VM,MC}-p_{z,VM,MC}) [GeV]");
	h_VM_Epz_response->GetXaxis()->SetTitle("(E_{VM,RECO}-p_{z,VM,RECO}) [GeV]");
    h_VM_Epz_response->Draw("colzsame");

    c777->cd(3);
    gPad->SetLogy(1);
	//h_VM_p_MC->GetXaxis()->SetRangeUser(0,2);
	h_VM_p_MC->GetYaxis()->SetTitle("counts");
	h_VM_p_MC->GetXaxis()->SetTitle("|p|_{VM} [GeV/c]");
	h_VM_p_MC->SetLineColor(kBlack);
    h_VM_p_MC->Draw();
    h_VM_p_REC->SetMarkerColor(kBlue);
    h_VM_p_REC->SetMarkerStyle(24);
    h_VM_p_REC->Draw("PEsame");
    //h_VM_p_REC_cut->SetMarkerColor(kBlue);
    //h_VM_p_REC_cut->SetMarkerStyle(24);
    //h_VM_p_REC_cut->Draw("PEsame");
    //TLegend *w977111 = new TLegend(0.52,0.2,0.6,0.35);
	//w977111->AddEntry(h_VM_p_MC, "MC", "L");
    //w977111->AddEntry(h_VM_p_REC, "RECO", "P");
	//w977111->AddEntry(h_VM_p_REC_cut, "RECO", "P");
	//w977111->Draw("same");
    w777->Draw("same");

    c777->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title777 = new TLatex();
    title777->SetNDC(); 
    title777->SetTextSize(0.05);
    title777->SetTextAlign(22);  
    title777->DrawLatex(0.5, 0.97, "E-p_{z} Truth vs Reco (VM)");  

    c777->cd(2);  
    gPad->SetTopMargin(0.08);  
    TLatex* title7_277 = new TLatex();
    title7_277->SetNDC(); 
    title7_277->SetTextSize(0.05);
    title7_277->SetTextAlign(22);  
    title7_277->DrawLatex(0.5, 0.97, "E-p_{z} Response (VM)");

    c777->cd(3);  
    gPad->SetTopMargin(0.08); 
    TLatex* title7_3773 = new TLatex();
    title7_3773->SetNDC(); 
    title7_3773->SetTextSize(0.05);
    title7_3773->SetTextAlign(22);  
    title7_3773->DrawLatex(0.5, 0.97, "|p|_{VM} Truth vs Reco");  

    c777->Print("./figures/plot_VM_all.pdf");


    TCanvas* c7771 = new TCanvas("c7771","c7771",1,1,1800,600);
    c7771->Divide(3,1,0.01,0.01);

    c7771->cd(1);
    gPad->SetLogy(1);
    //h_VM_pt_MC->GetXaxis()->SetRangeUser(0,2);
	h_VM_pt_MC->GetYaxis()->SetTitle("counts");
	h_VM_pt_MC->GetXaxis()->SetTitle("p_{T,VM} [GeV/c]");
	h_VM_pt_MC->SetLineColor(kBlack);
    h_VM_pt_MC->Draw();
    h_VM_pt_REC->SetMarkerColor(kBlue);
    h_VM_pt_REC->SetMarkerStyle(24);
    h_VM_pt_REC->Draw("PEsame");
    //h_VM_pt_REC_cut->SetMarkerColor(kBlue);
    //h_VM_pt_REC_cut->SetMarkerStyle(24);
    //h_VM_pt_REC_cut->Draw("PEsame");
    TLegend *w877111 = new TLegend(0.45,0.2,0.6,0.35);
	w877111->AddEntry(h_VM_pt_MC, "MC", "L");
	//w877111->AddEntry(h_VM_pt_REC_cut, "RECO", "P");
    w877111->AddEntry(h_VM_pt_REC, "RECO", "P");
	w877111->Draw("same");

    //c7771->cd(2);
    //gPad->SetLogz(1);
	//h_VM_pt_response_cut->GetYaxis()->SetTitle("p_{T,VM,MC} [GeV/c]");
	//h_VM_pt_response_cut->GetXaxis()->SetTitle("p_{T,VM,RECO} [GeV/c]");
    //h_VM_pt_response_cut->Draw("colzsame");

    c7771->cd(2);
    gPad->SetLogz(1);
    //h_VM_pt_response->GetXaxis()->SetRangeUser(0,2);
    //h_VM_pt_response->GetYaxis()->SetRangeUser(0,2);
	h_VM_pt_response->GetYaxis()->SetTitle("p_{T,VM,MC} [GeV/c]");
	h_VM_pt_response->GetXaxis()->SetTitle("p_{T,VM,RECO} [GeV/c]");
    h_VM_pt_response->Draw("colzsame");

    c7771->cd(3);
    gPad->SetLogy(1);
    //h_VM_mass_MC->GetXaxis()->SetRangeUser(0,4);
	h_VM_mass_MC->GetYaxis()->SetTitle("counts");
	h_VM_mass_MC->GetXaxis()->SetTitle("VM mass [GeV/c^{2}]");
	h_VM_mass_MC->SetLineColor(kBlack);
    h_VM_mass_MC->Draw();
    h_VM_mass_REC->SetMarkerColor(kBlue);
    //h_VM_mass_REC->SetLineStyle(9);
    //h_VM_mass_REC->Draw("same");
    h_VM_mass_REC->SetMarkerStyle(24);
    h_VM_mass_REC->Draw("PEsame");
    //h_VM_mass_REC_cut->SetMarkerColor(kBlue);
    //h_VM_mass_REC_cut->SetLineStyle(3);
    //h_VM_mass_REC_cut->Draw("same");
    //h_VM_mass_REC_cut->SetMarkerStyle(24);
    //h_VM_mass_REC_cut->Draw("PEsame");
    TLegend *w9771 = new TLegend(0.55,0.7,0.7,0.85);
	w9771->AddEntry(h_VM_mass_MC, " MC", "L");
    w9771->AddEntry(h_VM_mass_REC, " RECO", "P");
	//w9771->AddEntry(h_VM_mass_REC_cut, "RECO", "P");
	w9771->Draw("same");

    c7771->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title7_4771 = new TLatex();
    title7_4771->SetNDC(); 
    title7_4771->SetTextSize(0.05);
    title7_4771->SetTextAlign(22);  
    title7_4771->DrawLatex(0.5, 0.97, "p_{T} Truth vs Reco (VM)");  

    c7771->cd(2);  
    gPad->SetTopMargin(0.08);  
    TLatex* title7_5772 = new TLatex();
    title7_5772->SetNDC(); 
    title7_5772->SetTextSize(0.05);
    title7_5772->SetTextAlign(22);  
    title7_5772->DrawLatex(0.5, 0.97, "p_{T} Response (VM)");

    c7771->cd(3);  
    gPad->SetTopMargin(0.08);  
    TLatex* title7_677 = new TLatex();
    title7_677->SetNDC(); 
    title7_677->SetTextSize(0.05);
    title7_677->SetTextAlign(22);  
    title7_677->DrawLatex(0.5, 0.97, "m_{VM} Truth vs Reco");

    c7771->Print("./figures/plot_VM_all2.pdf");


// VM pT
    TCanvas* c8 = new TCanvas("c8","c8",1,1,1200,800);
    c8->Divide(3,1,0.01,0.01);
    c8->cd(1);
    gPad->SetLogy(1);
	h_VM_pt_MC->GetYaxis()->SetTitle("counts");
	h_VM_pt_MC->GetXaxis()->SetTitle("p_{T,VM} [GeV/c]");
	h_VM_pt_MC->SetLineColor(kBlack);
    h_VM_pt_MC->Draw();
    h_VM_pt_REC->SetMarkerColor(kBlue);
    h_VM_pt_REC->SetMarkerStyle(24);
    h_VM_pt_REC->Draw("PEsame");
    //h_VM_pt_REC_before->SetMarkerStyle(6);
    //h_VM_pt_REC_before->SetMarkerColor(kRed);
    //h_VM_pt_REC_before->Draw("PEsame");
    //h_VM_pt_MC_before->SetMarkerStyle(21);
    //h_VM_pt_MC_before->SetMarkerColor(kRed);
    //h_VM_pt_MC_before->Draw("PEsame");
    h_VM_pt_REC_cut->SetMarkerColor(kCyan);
    h_VM_pt_REC_cut->SetMarkerStyle(24);
    //h_VM_pt_REC_cut->Draw("PEsame");
    TLegend *w8 = new TLegend(0.55,0.3,0.65,0.45);
	w8->AddEntry(h_VM_pt_MC, " MC", "L");
	//w8->AddEntry(h_VM_pt_REC_cut, "RECO", "P");
    w8->AddEntry(h_VM_pt_REC, " RECO", "P");
    //w8->AddEntry(h_VM_pt_REC_before, "RECO before cut", "P");
    //w8->AddEntry(h_VM_pt_MC, Form("MC = %.0f", h_VM_pt_MC->GetEntries()), "L");
    //w8->AddEntry(h_VM_pt_REC, Form("RECO = %.0f", h_VM_pt_REC->GetEntries()), "P");
	w8->Draw("same");

    c8->cd(2);
    gPad->SetLogz(1);
	h_VM_pt_response->GetYaxis()->SetTitle("p_{T,VM,MC} [GeV/c]");
	h_VM_pt_response->GetXaxis()->SetTitle("p_{T,VM,RECO} [GeV/c]");
    h_VM_pt_response->Draw("colzsame");

    c8->cd(3);
    h_VM_pt_resolution->Draw("colzsame");

    //c8->cd(2);
    //gPad->SetLogz(1);
	//h_VM_pt_response_cut->GetYaxis()->SetTitle("p_{T,VM,MC} [GeV/c]");
	//h_VM_pt_response_cut->GetXaxis()->SetTitle("p_{T,VM,RECO} [GeV/c]");
    //h_VM_pt_response_cut->Draw("colzsame");

    //c8->cd(3);
    //gPad->SetLogy(1);
    //h_VM_pt_purity->GetYaxis()->SetTitle("bins");
    //h_VM_pt_purity->Draw("same");

    //c8->cd(4);
    //gPad->SetLogz(1);
    //h_VM_pt_migration->Draw("colzsame");

    //c8->cd(5);
    //gPad->SetLogy(1);
    //h_VM_pt_MC->GetYaxis()->SetTitle("counts");
	//h_VM_pt_MC->GetXaxis()->SetTitle("p_{T,VM} [GeV/c]");
	//h_VM_pt_MC->SetLineColor(kBlack);
    //h_VM_pt_MC->Draw();
    //h_VM_pt_corrected->SetMarkerStyle(24);
    //h_VM_pt_corrected->SetMarkerColor(kBlue);
    //h_VM_pt_corrected->Draw("PEsame");
    //TLegend *w8_1 = new TLegend(0.57,0.2,0.75,0.35);
	//w8_1->AddEntry(h_VM_pt_MC, "MC", "L");
    //w8_1->AddEntry(h_VM_pt_corrected, "Corrected", "P");
	//w8_1->Draw("same");

    //c8->cd(6);
    //gPad->SetLogy(1);
    //h_VM_pt_acceptance->Draw("same");

    //c8->cd(3);  
    //gPad->SetTopMargin(0.08);  
    //TLatex* title8_3 = new TLatex();
    //title8_3->SetNDC(); 
    //title8_3->SetTextSize(0.05);
    //title8_3->SetTextAlign(22);  
    //title8_3->DrawLatex(0.5, 0.97, "Purity");  

    c8->cd(4);  
    gPad->SetTopMargin(0.08);  
    TLatex* title8_4 = new TLatex();
    title8_4->SetNDC(); 
    title8_4->SetTextSize(0.05);
    title8_4->SetTextAlign(22);  
    //title8_4->DrawLatex(0.5, 0.97, "Bin Migration"); 

    c8->cd(5);  
    gPad->SetTopMargin(0.08);  
    TLatex* title8_5 = new TLatex();
    title8_5->SetNDC(); 
    title8_5->SetTextSize(0.05);
    title8_5->SetTextAlign(22);  
    //title8_5->DrawLatex(0.5, 0.97, "Corrected Acceptance"); 

    c8->cd(6);  
    gPad->SetTopMargin(0.08);  
    TLatex* title8_6 = new TLatex();
    title8_6->SetNDC(); 
    title8_6->SetTextSize(0.05);
    title8_6->SetTextAlign(22);  
    //title8_6->DrawLatex(0.5, 0.97, "Acceptance"); 

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

    c8->cd(3);  
    gPad->SetTopMargin(0.08);  
    TLatex* title8_3 = new TLatex();
    title8_3->SetNDC(); 
    title8_3->SetTextSize(0.05);
    title8_3->SetTextAlign(22);  
    title8_3->DrawLatex(0.5, 0.97, "p_{T} Resolution (VM)");  

    c8->Print("./figures/plot_VM_pt.pdf");


// VM 
    TCanvas* c9 = new TCanvas("c9","c9",1,1,1600,600);
    c9->Divide(3,1,0.01,0.01);
    c9->cd(1);
    gPad->SetLogy(1);
	h_VM_p_MC->GetYaxis()->SetTitle("counts");
	h_VM_p_MC->GetXaxis()->SetTitle("|p|_{VM} [GeV/c]");
	h_VM_p_MC->SetLineColor(kBlack);
    h_VM_p_MC->Draw();
    h_VM_p_REC->SetMarkerColor(kBlue);
    h_VM_p_REC->SetMarkerStyle(24);
    h_VM_p_REC->Draw("PEsame");
    //h_VM_p_REC_before->SetMarkerStyle(6);
    //h_VM_p_REC_before->SetMarkerColor(kRed);
    //h_VM_p_REC_before->Draw("PEsame");
    //h_VM_p_MC_before->SetMarkerStyle(21);
    //h_VM_p_MC_before->SetMarkerColor(kRed);
    //h_VM_p_MC_before->Draw("PEsame");
    //h_VM_p_REC_cut->SetMarkerColor(kBlue);
    //h_VM_p_REC_cut->SetMarkerStyle(24);
    //h_VM_p_REC_cut->Draw("PEsame");
    TLegend *w9 = new TLegend(0.52,0.7,0.7,0.85);
	w9->AddEntry(h_VM_p_MC, " MC", "L");
    w9->AddEntry(h_VM_p_REC, " RECO", "P");
    //w9->AddEntry(h_VM_p_REC_before, "RECO bc", "P");
   // w9->AddEntry(h_VM_p_MC, Form("MC = %.0f", h_VM_p_MC->GetEntries()), "");
   // w9->AddEntry(h_VM_p_REC, Form("RECO = %.0f", h_VM_p_REC->GetEntries()), "");
	//w9->AddEntry(h_VM_p_REC_cut, "reco (inside threshold)", "P");
	w9->Draw("same");

    cout << "VM MC |p|: " << h_VM_p_MC->GetMean() << " VM rec |p|: " << h_VM_p_REC->GetMean() << endl;
    cout << "VM MC theta = " << h_VM_theta_MC->GetMean() << " VM rec theta = " << h_VM_theta_REC->GetMean() << endl;

    c9->cd(2);
    gPad->SetLogy(1);
	h_VM_pz_MC->GetYaxis()->SetTitle("counts");
	h_VM_pz_MC->GetXaxis()->SetTitle("p_{z,VM} [GeV/c]");
	h_VM_pz_MC->SetLineColor(kBlack);
    h_VM_pz_MC->Draw();
    h_VM_pz_REC->SetMarkerColor(kBlue);
    h_VM_pz_REC->SetMarkerStyle(24);
    h_VM_pz_REC->Draw("PEsame");
    //h_VM_pz_REC_before->SetMarkerStyle(6);
    //h_VM_pz_REC_before->SetMarkerColor(kRed);
    //h_VM_pz_REC_before->Draw("PEsame");
    //h_VM_pz_MC_before->SetMarkerStyle(21);
    //h_VM_pz_MC_before->SetMarkerColor(kRed);
    //h_VM_pz_MC_before->Draw("PEsame");
    //h_VM_pz_REC_cut->SetMarkerColor(kBlue);
    //h_VM_pz_REC_cut->SetMarkerStyle(24);
    //h_VM_pz_REC_cut->Draw("PEsame");
    //TLegend *w919 = new TLegend(0.32,0.7,0.45,0.85);
    //w919->AddEntry(h_VM_pz_MC, Form("MC = %.0f", h_VM_pz_MC->GetEntries()), "L");
    //w919->AddEntry(h_VM_pz_REC, Form("RECO = %.0f", h_VM_pz_REC->GetEntries()), "P");
	//w919->AddEntry(h_VM_pz_MC, " MC", "P");
    //w919->AddEntry(h_VM_pz_REC, " RECO", "P");
    //w919->AddEntry(h_VM_pz_REC_before, "RECO before cut", "P");
	w9->Draw("same");

    c9->cd(3);
    gPad->SetLogy(1);
	h_VM_mass_MC->GetYaxis()->SetTitle("counts");
	h_VM_mass_MC->GetXaxis()->SetTitle("VM mass [GeV/c^{2}]");
	h_VM_mass_MC->SetLineColor(kBlack);
    h_VM_mass_MC->Draw();
    h_VM_mass_REC->SetMarkerColor(kBlue);
    //h_VM_mass_REC->SetLineStyle(9);
    //h_VM_mass_REC->Draw("same");
    h_VM_mass_REC->SetMarkerStyle(24);
    h_VM_mass_REC->Draw("PEsame");
    //h_VM_mass_REC_before->SetMarkerStyle(6);
    //h_VM_mass_REC_before->SetMarkerColor(kRed);
    //h_VM_mass_REC_before->Draw("PEsame");
    //h_VM_mass_MC_before->SetMarkerStyle(21);
    //h_VM_mass_MC_before->SetMarkerColor(kRed);
    //h_VM_mass_MC_before->Draw("PEsame");
    //h_VM_mass_REC_cut->SetLineColor(kBlue);
    //h_VM_mass_REC_cut->SetLineStyle(3);
    //h_VM_mass_REC_cut->Draw("same");
    //h_VM_mass_REC_cut->SetMarkerStyle(24);
    //h_VM_mass_REC_cut->Draw("PEsame");
    w9->Draw("same");

    c9->Print("./figures/plot_VM.pdf");

// VM compare
    //TFile* truthSeeded = TFile::Open("merged_sartre_TruthSeeded_10x100eAu_v1.root","READ");
    //TH1D* h_truthSeeded_VM_pt_REC = (TH1D*) truthSeeded->Get("h_VM_pt_REC");
    //TH1D* h_truthSeeded_VM_pt_MC = (TH1D*) truthSeeded->Get("h_VM_pt_MC");
    //TH2D* h_truthpt_resolution = (TH2D*) truthSeeded->Get("h_VM_pt_resolution");
    //TH1D* h_truthpt_res_diff = (TH1D*) truthSeeded->Get("h_VM_pt_res_diff");
    //TH2D* h_truthpt_res_diff2d = (TH2D*) truthSeeded->Get("h_VM_pt_res_diff2d");
    //TH1D* h_truthpt_resolution_Y = h_truthpt_resolution->ProjectionY();
    //TH1D* h_truthpt_resolution_X = h_truthpt_resolution->ProjectionX();
    //TH1D* h_truthSeeded_VM_pz_REC = (TH1D*) truthSeeded->Get("h_VM_pz_REC");
    //TH1D* h_truthSeeded_VM_pz_MC = (TH1D*) truthSeeded->Get("h_VM_pz_MC");
    //TH2D* h_truthpz_resolution = (TH2D*) truthSeeded->Get("h_VM_pz_resolution");
    //TH1D* h_truthpz_res_diff = (TH1D*) truthSeeded->Get("h_VM_pz_res_diff");
    //TH2D* h_truthpz_res_diff2d = (TH2D*) truthSeeded->Get("h_VM_pz_res_diff2d");
    //TH1D* h_truthpz_resolution_Y = h_truthpz_resolution->ProjectionY();
    //TH1D* h_truthpz_resolution_X = h_truthpz_resolution->ProjectionX();
    //TH1D* h_truthSeeded_VM_p_MC = (TH1D*) truthSeeded->Get("h_VM_p_MC");
    //TH1D* h_truthSeeded_VM_p_REC = (TH1D*) truthSeeded->Get("h_VM_p_REC");
    //TH2D* h_truthp_resolution = (TH2D*) truthSeeded->Get("h_VM_p_resolution");
    //TH1D* h_truthp_res_diff = (TH1D*) truthSeeded->Get("h_VM_p_res_diff");
    //TH2D* h_truthp_res_diff2d = (TH2D*) truthSeeded->Get("h_VM_p_res_diff2d");
    //TH1D* h_truthp_resolution_Y = h_truthp_resolution->ProjectionY();
    //TH1D* h_truthp_resolution_X = h_truthp_resolution->ProjectionX();

    TCanvas* c9_213 = new TCanvas("c9_213","c9_213",1,1,1200,600);
    c9_213->Divide(3,1,0.01,0.01);
    c9_213->cd(1);
    gPad->SetLogy(1);
    h_VM_pt_MC->GetXaxis()->SetRangeUser(0,2);
	h_VM_pt_MC->GetYaxis()->SetTitle("counts");
	h_VM_pt_MC->GetXaxis()->SetTitle("p_{T,VM} [GeV/c]");
	h_VM_pt_MC->SetLineColor(kBlack);
    h_VM_pt_MC->Draw();
    h_VM_pt_REC->SetMarkerColor(kBlue);
    h_VM_pt_REC->SetMarkerStyle(24);
    h_VM_pt_REC->Draw("PEsame");
    h_VM_pt_REC_cut->SetMarkerColor(kBlue);
    h_VM_pt_REC_cut->SetMarkerStyle(24);
    //h_VM_pt_REC_cut->Draw("PEsame");
    //h_truthSeeded_VM_pt_MC->SetMarkerColor(kOrange);
    //h_truthSeeded_VM_pt_MC->SetMarkerStyle(21);
    //h_truthSeeded_VM_pt_MC->Draw("PEsame");
    //h_truthSeeded_VM_pt_REC->SetMarkerColor(kRed);
    //h_truthSeeded_VM_pt_REC->SetMarkerStyle(30);
    //h_truthSeeded_VM_pt_REC->Draw("PEsame");
    TLegend *w8771 = new TLegend(0.45,0.2,0.6,0.35);
    w8771->AddEntry(h_VM_pt_MC, " MC","L");//Form("MC real= %.0f", h_VM_pt_MC->GetEntries()), "L");
    w8771->AddEntry(h_VM_pt_REC, " RECO","L");//Form("real= %.0f", h_VM_pt_REC->GetEntries()), "P");
    //w8771->AddEntry(h_truthSeeded_VM_pt_REC, " truth seeded","P");//Form("truth= %.0f", h_truthSeeded_VM_pt_REC->GetEntries()), "P");
    //w8771->AddEntry(h_truthSeeded_VM_pt_MC,Form("MC truth= %.0f", h_truthSeeded_VM_pt_MC->GetEntries()), "P");
    w8771->Draw("same");

    c9_213->cd(2);
    gPad->SetLogy(1);
    h_VM_pz_MC->GetXaxis()->SetRangeUser(0,2);
	h_VM_pz_MC->GetYaxis()->SetTitle("counts");
	h_VM_pz_MC->GetXaxis()->SetTitle("p_{z,VM} [GeV/c]");
	h_VM_pz_MC->SetLineColor(kBlack);
    h_VM_pz_MC->Draw();
    h_VM_pz_REC->SetMarkerColor(kBlue);
    h_VM_pz_REC->SetMarkerStyle(24);
    h_VM_pz_REC->Draw("PEsame");
    //h_VM_pz_MC_before->SetMarkerStyle(21);
    //h_VM_pz_MC_before->SetMarkerColor(kRed);
    //h_VM_pz_MC_before->Draw("PEsame");
    h_VM_pz_REC_cut->SetMarkerColor(kBlue);
    h_VM_pz_REC_cut->SetMarkerStyle(24);
    //h_VM_pz_REC_cut->Draw("PEsame");
    //h_truthSeeded_VM_pz_REC->SetMarkerStyle(30);
    //h_truthSeeded_VM_pz_REC->SetMarkerColor(kRed);
    //h_truthSeeded_VM_pz_REC->Draw("PEsame");
    TLegend *w91999 = new TLegend(0.42,0.7,0.55,0.85);
    w91999->AddEntry(h_VM_pz_MC, " MC", "L");
    w91999->AddEntry(h_VM_pz_REC, " RECO", "P");
    //w91999->AddEntry(h_truthSeeded_VM_pz_REC," truth seeded", "P");
	w91999->Draw("same");

    c9_213->cd(3);
    gPad->SetLogy(1);
    h_VM_p_MC->GetXaxis()->SetRangeUser(0,2);
	h_VM_p_MC->GetYaxis()->SetTitle("counts");
	h_VM_p_MC->GetXaxis()->SetTitle("|p|_{VM} [GeV/c]");
	h_VM_p_MC->SetLineColor(kBlack);
    h_VM_p_MC->Draw();
    h_VM_p_REC->SetMarkerColor(kBlue);
    h_VM_p_REC->SetMarkerStyle(24);
    h_VM_p_REC->Draw("PEsame");
    h_VM_p_REC_cut->SetMarkerColor(kBlue);
    h_VM_p_REC_cut->SetMarkerStyle(24);
    //h_VM_p_REC_cut->Draw("PEsame");
    //h_truthSeeded_VM_p_REC->SetMarkerStyle(30);
    //h_truthSeeded_VM_p_REC->SetMarkerColor(kRed);
    //h_truthSeeded_VM_p_REC->Draw("PEsame");
    //h_truthSeeded_VM_p_MC->SetMarkerColor(kOrange);
    //h_truthSeeded_VM_p_MC->SetMarkerStyle(11);
    //h_truthSeeded_VM_p_MC->Draw("PEsame");
    TLegend *w1977 = new TLegend(0.45,0.2,0.6,0.3);
    w1977->AddEntry(h_VM_p_MC, " MC","L");//Form("MC real = %.0f", h_VM_p_MC->GetEntries()), "L");
    w1977->AddEntry(h_VM_p_REC, " RECO","P");//Form("real = %.0f", h_VM_p_REC->GetEntries()), "P");
    //w1977->AddEntry(h_truthSeeded_VM_p_REC, " truth seeded","P");//Form("truth = %.0f", h_truthSeeded_VM_p_REC->GetEntries()), "P");
    //w1977->AddEntry(h_truthSeeded_VM_p_MC, Form("MC truth = %.0f", h_truthSeeded_VM_p_MC->GetEntries()), "P");
	w1977->Draw("same");
	//w8771->Draw("same");

    c9_213->Print("./figures/plot_VM_compare.pdf");



    TCanvas* c9_213_1 = new TCanvas("c9_213_1","c9_213_1",1,1,1200,600);
    c9_213_1->Divide(3,1,0.01,0.01);
    c9_213_1->cd(1);
    //gPad->SetLogy(1);
    h_VM_pt_MC->GetXaxis()->SetRangeUser(0,10);
    h_VM_pt_MC->GetYaxis()->SetTitle("counts");
	h_VM_pt_MC->GetXaxis()->SetTitle("p_{T,VM} [GeV/c]");
	h_VM_pt_MC->SetLineColor(kBlack);
    h_VM_pt_MC->Draw();
    h_VM_pt_REC->SetMarkerColor(kBlue);
    h_VM_pt_REC->SetMarkerStyle(24);
    h_VM_pt_REC->Draw("PEsame");
    //h_truthSeeded_VM_pt_REC->SetMarkerColor(kRed);
    //h_truthSeeded_VM_pt_REC->SetMarkerStyle(30);
    //h_truthSeeded_VM_pt_REC->Draw("PEsame");
	w8771->Draw("same");

    c9_213_1->cd(2);
    //gPad->SetLogy(1);
    h_VM_pz_MC->GetXaxis()->SetRangeUser(0,10);
	h_VM_pz_MC->GetYaxis()->SetTitle("counts");
	h_VM_pz_MC->GetXaxis()->SetTitle("p_{z,VM} [GeV/c]");
	h_VM_pz_MC->SetLineColor(kBlack);
    h_VM_pz_MC->Draw();
    h_VM_pz_REC->SetMarkerColor(kBlue);
    h_VM_pz_REC->SetMarkerStyle(24);
    h_VM_pz_REC->Draw("PEsame");
    //h_truthSeeded_VM_pz_REC->SetMarkerStyle(30);
    //h_truthSeeded_VM_pz_REC->SetMarkerColor(kRed);
    //h_truthSeeded_VM_pz_REC->Draw("PEsame");
    //h_VM_pz_MC_before->SetMarkerStyle(21);
    //h_VM_pz_MC_before->SetMarkerColor(kRed);
    //h_VM_pz_MC_before->Draw("PEsame");
    //h_VM_pz_REC_cut->SetMarkerColor(kBlue);
    //h_VM_pz_REC_cut->SetMarkerStyle(24);
    //h_VM_pz_REC_cut->Draw("PEsame");
    //TLegend *w91999 = new TLegend(0.32,0.75,0.45,0.9);
    //w91999->AddEntry(h_VM_pz_MC, Form("MC = %.0f", h_VM_pz_MC->GetEntries()), "L");
    //w91999->AddEntry(h_VM_pz_REC, Form("RECO = %.0f", h_VM_pz_REC->GetEntries()), "P");
	//w91999->Draw("same");
	w91999->Draw("same");

    c9_213_1->cd(3);
    //gPad->SetLogy(1);
    h_VM_p_MC->GetXaxis()->SetRangeUser(0,10);
	h_VM_p_MC->GetYaxis()->SetTitle("counts");
	h_VM_p_MC->GetXaxis()->SetTitle("|p|_{VM} [GeV/c]");
	h_VM_p_MC->SetLineColor(kBlack);
    h_VM_p_MC->Draw();
    h_VM_p_REC->SetMarkerColor(kBlue);
    h_VM_p_REC->SetMarkerStyle(24);
    h_VM_p_REC->Draw("PEsame");
    //h_truthSeeded_VM_p_REC->SetMarkerStyle(30);
    //h_truthSeeded_VM_p_REC->SetMarkerColor(kRed);
    //h_truthSeeded_VM_p_REC->Draw("PEsame");
    //h_VM_p_REC_cut->SetMarkerColor(kBlue);
    //h_VM_p_REC_cut->SetMarkerStyle(24);
    //h_VM_p_REC_cut->Draw("PEsame");
    TLegend *w19777 = new TLegend(0.45,0.8,0.55,0.9);
    w19777->AddEntry(h_VM_p_MC, " MC", "L");
    w19777->AddEntry(h_VM_p_REC, " real seeded", "P");
    //w19777->AddEntry(h_truthSeeded_VM_p_REC, " truth seeded", "P");
	w19777->Draw("same");

    c9_213_1->Print("./figures/plot_VM_compare_1.pdf");

// VM resolution
    TCanvas* c9_213_2 = new TCanvas("c9_213_2","c9_213_2",1,1,1200,600);
    c9_213_2->Divide(3,2,0.01,0.01);
    c9_213_2->cd(1);
    gPad->SetLogz(1);
    //h_VM_pt_resolution->GetYaxis()->SetRangeUser(-0.1,0.1);
    //h_VM_pt_resolution->GetXaxis()->SetRangeUser(0,8);
    //h_VM_pt_resolution->Draw("colzsame");
   
    c9_213_2->cd(1);  
    gPad->SetTopMargin(0.08);  
    TLatex* title9_213_2 = new TLatex();
    title9_213_2->SetNDC(); 
    title9_213_2->SetTextSize(0.05);
    title9_213_2->SetTextAlign(22);  
    title9_213_2->DrawLatex(0.5, 0.97, "VM p_{T} Resolution (real seeded)");  
        
    
    c9_213_2->cd(4);
    gPad->SetLogz(1);
    //h_truthpt_resolution->Draw("colzsame");

    c9_213_2->cd(4);  
    gPad->SetTopMargin(0.08);  
    TLatex* title9_213_22 = new TLatex();
    title9_213_22->SetNDC(); 
    title9_213_22->SetTextSize(0.05);
    title9_213_22->SetTextAlign(22);  
    title9_213_22->DrawLatex(0.5, 0.97, "VM p_{T} Resolution (truth seeded)");  
    
    c9_213_2->cd(2);
    gPad->SetLogz(1);
    //h_VM_pz_resolution->Draw("colzsame");

    c9_213_2->cd(2);  
    gPad->SetTopMargin(0.08);  
    TLatex* title9_213_23 = new TLatex();
    title9_213_23->SetNDC(); 
    title9_213_23->SetTextSize(0.05);
    title9_213_23->SetTextAlign(22);  
    title9_213_23->DrawLatex(0.5, 0.97, "VM p_{z} Resolution (real seeded)"); 
    
    
    c9_213_2->cd(5);
    gPad->SetLogz(1);  
    //h_truthpz_resolution->Draw("colzsame");

    c9_213_2->cd(5);  
    gPad->SetTopMargin(0.08);  
    TLatex* title9_213_24 = new TLatex();
    title9_213_24->SetNDC(); 
    title9_213_24->SetTextSize(0.05);
    title9_213_24->SetTextAlign(22);  
    title9_213_24->DrawLatex(0.5, 0.97, "VM p_{z} Resolution (truth seeded)"); 

    c9_213_2->cd(3);
    gPad->SetLogz(1);
    //h_VM_p_resolution->Draw("colzsame");

    c9_213_2->cd(3);  
    gPad->SetTopMargin(0.08);  
    TLatex* title9_213_25 = new TLatex();
    title9_213_25->SetNDC(); 
    title9_213_25->SetTextSize(0.05);
    title9_213_25->SetTextAlign(22);  
    title9_213_25->DrawLatex(0.5, 0.97, "VM |p| Resolution (real seeded)"); 
        
    
    c9_213_2->cd(6);
    gPad->SetLogz(1);
    //h_truthp_resolution->Draw("colzsame");

    c9_213_2->cd(6);  
    gPad->SetTopMargin(0.08);  
    TLatex* title9_213_26 = new TLatex();
    title9_213_26->SetNDC(); 
    title9_213_26->SetTextSize(0.05);
    title9_213_26->SetTextAlign(22);  
    title9_213_26->DrawLatex(0.5, 0.97, "VM |p| Resolution (truth seeded)"); 

    c9_213_2->Print("./figures/plot_VM_compare_resolution.pdf");




// VM resolution profiles
    TCanvas* c9_213_22_2 = new TCanvas("c9_213_22_2","c9_213_22_2",1,1,1600,600);
    c9_213_22_2->Divide(3,1,0.01,0.01);
    c9_213_22_2->cd(1);
    gPad->SetLogy(0);
    gPad->Clear();
    TProfile* prof_res_VM_pt = h_VM_pt_resolution->ProfileY();  
    //prof_res_VM_pt->GetYaxis()->SetRangeUser(-0.004,0.004);
    prof_res_VM_pt->SetLineColor(kBlue);
    prof_res_VM_pt->SetMarkerStyle(24);
    prof_res_VM_pt->SetMarkerColor(kBlue);
    prof_res_VM_pt->Draw();

    c9_213_22_2->cd(1);  
    gPad->SetTopMargin(0.08);  
    TLatex* title9_213_1_2 = new TLatex();
    title9_213_1_2->SetNDC(); 
    title9_213_1_2->SetTextSize(0.05);
    title9_213_1_2->SetTextAlign(22);  
    title9_213_1_2->DrawLatex(0.5, 0.97, "VM p_{T} Resolution Profile Real Seeded"); 
    
    c9_213_22_2->cd(2);
    gPad->SetLogy(0);
    TProfile* prof_res_VM_pz = h_VM_pz_resolution->ProfileY();  
    //prof_res_VM_pz->GetYaxis()->SetRangeUser(-0.01,0.01);
    prof_res_VM_pz->SetLineColor(kBlue);
    prof_res_VM_pz->SetMarkerStyle(24);
    prof_res_VM_pz->SetMarkerColor(kBlue);
    prof_res_VM_pz->Draw();


    c9_213_22_2->cd(2);  
    gPad->SetTopMargin(0.08);  
    TLatex* title9_213_12_22 = new TLatex();
    title9_213_12_22->SetNDC(); 
    title9_213_12_22->SetTextSize(0.05);
    title9_213_12_22->SetTextAlign(22);  
    title9_213_12_22->DrawLatex(0.5, 0.97, "VM p_{z} Resolution Profile Real Seeded"); 

    c9_213_22_2->cd(3);
    gPad->SetLogy(0);
    TProfile* prof_res_VM_p = h_VM_p_resolution->ProfileY();  
    //prof_res_VM_p->GetYaxis()->SetRangeUser(-0.01,0.01);
    prof_res_VM_p->SetLineColor(kBlue);
    prof_res_VM_p->SetMarkerStyle(24);
    prof_res_VM_p->SetMarkerColor(kBlue);
    prof_res_VM_p->Draw();
    
    c9_213_22_2->cd(3);  
    gPad->SetTopMargin(0.08);  
    TLatex* title9_213_13_23 = new TLatex();
    title9_213_13_23->SetNDC(); 
    title9_213_13_23->SetTextSize(0.05);
    title9_213_13_23->SetTextAlign(22);  
    title9_213_13_23->DrawLatex(0.5, 0.97, "VM |p| Resolution Profile Real Seeded"); 

    c9_213_22_2->Print("./figures/plot_VM_compare_resolution_profiles_realSeeded.pdf");

    TCanvas* c9_213_22_211 = new TCanvas("c9_213_22_211","c9_213_22_211",1,1,1600,600);
    c9_213_22_211->Divide(3,1,0.01,0.01);
    c9_213_22_211->cd(1);
    //TProfile* prof_res_VM_pt_truth = h_truthpt_resolution->ProfileY();  
    gPad->SetLogy(0);
    //prof_res_VM_pt_truth->GetYaxis()->SetRangeUser(-0.004,0.004);
    //prof_res_VM_pt_truth->SetLineColor(kRed);
    //prof_res_VM_pt_truth->SetMarkerColor(kRed);
    //prof_res_VM_pt_truth->SetMarkerStyle(24);
    //prof_res_VM_pt_truth->Draw();

    c9_213_22_211->cd(1);  
    gPad->SetTopMargin(0.08);  
    TLatex* title9_213_1_2411 = new TLatex();
    title9_213_1_2411->SetNDC(); 
    title9_213_1_2411->SetTextSize(0.05);
    title9_213_1_2411->SetTextAlign(22);  
    title9_213_1_2411->DrawLatex(0.5, 0.97, "VM p_{T} Resolution Profile Truth Seeded"); 
    
    c9_213_22_211->cd(2);
    gPad->SetLogy(0);
    //TProfile* prof_res_VM_pz_truth = h_truthpz_resolution->ProfileY();  
    //prof_res_VM_pz_truth->GetYaxis()->SetRangeUser(-0.01,0.01);
    //prof_res_VM_pz_truth->SetLineColor(kRed);
    //prof_res_VM_pz_truth->SetMarkerColor(kRed);
    //prof_res_VM_pz_truth->SetMarkerStyle(24);
    //prof_res_VM_pz_truth->Draw();

    c9_213_22_211->cd(2);  
    gPad->SetTopMargin(0.08);  
    TLatex* title9_213_12_2511 = new TLatex();
    title9_213_12_2511->SetNDC(); 
    title9_213_12_2511->SetTextSize(0.05);
    title9_213_12_2511->SetTextAlign(22);  
    title9_213_12_2511->DrawLatex(0.5, 0.97, "VM p_{z} Resolution Profile Truth Seeded"); 

    c9_213_22_211->cd(3);
    gPad->SetLogy(0);
    //TProfile* prof_res_VM_p_truth = h_truthp_resolution->ProfileY();  
    //prof_res_VM_p_truth->GetYaxis()->SetRangeUser(-0.01,0.01);
    //prof_res_VM_p_truth->SetLineColor(kRed);
    //prof_res_VM_p_truth->SetMarkerColor(kRed);
    //prof_res_VM_p_truth->SetMarkerStyle(24);
    //prof_res_VM_p_truth->Draw();
    
    c9_213_22_211->cd(3);  
    gPad->SetTopMargin(0.08);  
    TLatex* title9_213_13_23611 = new TLatex();
    title9_213_13_23611->SetNDC(); 
    title9_213_13_23611->SetTextSize(0.05);
    title9_213_13_23611->SetTextAlign(22);  
    title9_213_13_23611->DrawLatex(0.5, 0.97, "VM |p| Resolution Profile Truth Seeded"); 

    c9_213_22_211->Print("./figures/plot_VM_compare_resolution_profiles_TruthSeeded.pdf");



// position
    TCanvas* c10 = new TCanvas("c10","c10",1,1,1600,800);
    c10->Divide(2,1,0.01,0.01);
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

// geometric position and radius
    TCanvas* c1011 = new TCanvas("c1011","c1011",1,1,1600,800);
    c1011->Divide(3,1,0.01,0.01);
    c1011->cd(1);
    gPad->SetLogy(1);
    h_Xclus_minus_Xtrk->GetYaxis()->SetTitle("counts");
	h_Xclus_minus_Xtrk->GetXaxis()->SetTitle("x_{clus}-x_{trk} [mm]");
    h_Xclus_minus_Xtrk->SetLineColor(kBlack);
    //h_Xclus_minus_Xtrk->SetMarkerColor(kBlack);
    h_Xclus_minus_Xtrk->Draw();
    h_Xclus_minus_Xtrk_cut->SetMarkerColor(kBlue);
    h_Xclus_minus_Xtrk_cut->SetMarkerStyle(24);
    h_Xclus_minus_Xtrk_cut->Draw("PEsame");
    h_Xclus_minus_Xtrk_cut2->SetMarkerColor(kRed);
    h_Xclus_minus_Xtrk_cut2->SetMarkerStyle(29);
    h_Xclus_minus_Xtrk_cut2->Draw("PEsame");
    TLegend *w9911 = new TLegend(0.22,0.7,0.4,0.85);
    w9911->AddEntry(h_Xclus_minus_Xtrk, "RECO", "L");
	w9911->AddEntry(h_Xclus_minus_Xtrk_cut, "Accept", "P");
    w9911->AddEntry(h_Xclus_minus_Xtrk_cut2, "Reject", "P");
	w9911->Draw("same");
   
    c1011->cd(2);
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
    w9911->Draw("same");

    c1011->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title1011 = new TLatex();
    title1011->SetNDC(); 
    title1011->SetTextSize(0.05);
    title1011->SetTextAlign(22);  
    title1011->DrawLatex(0.5, 0.97, "x Position Difference");  

    c1011->cd(2);  
    gPad->SetTopMargin(0.08); 
    TLatex* title10_211 = new TLatex();
    title10_211->SetNDC(); 
    title10_211->SetTextSize(0.05);
    title10_211->SetTextAlign(22);  
    title10_211->DrawLatex(0.5, 0.97, "y Position Difference"); 

    c1011->cd(3);
    gPad->SetLogz(1);
    h_emClus_position_REC_cut->GetYaxis()->SetTitleOffset(1.5);
    h_emClus_position_REC_cut->GetXaxis()->SetTitle("x [mm]");
    h_emClus_position_REC_cut->GetYaxis()->SetTitle("y [mm]");
    h_emClus_position_REC_cut->Draw("colzsame"); 

    c1011->cd(3);  
    gPad->SetTopMargin(0.08); 
    TLatex* title10_411 = new TLatex();
    title10_411->SetNDC(); 
    title10_411->SetTextSize(0.05);
    title10_411->SetTextAlign(22);  
    title10_411->DrawLatex(0.5, 0.97, "Cluster Positions After Threshold Cut");  

    c1011->Print("./figures/plot_geometric.pdf");


// cluster positions
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
    TCanvas* c11 = new TCanvas("c11","c11",1,1,1600,600);
    c11->Divide(3,1,0.01,0.01);
    c11->cd(1);
    //gPad->SetLogy(1);
    h_Epz_MC->GetXaxis()->SetRangeUser(15,25);
	h_Epz_MC->GetYaxis()->SetTitle("counts");
	h_Epz_MC->GetXaxis()->SetTitle("(E - p_{z}) [GeV]");
	h_Epz_MC->SetLineColor(kBlack);
    //h_Epz_MC->Draw();
    //h_Epz_MC_after->SetMarkerStyle(6);
    h_Epz_MC_after->SetLineColor(kBlack);
    h_Epz_MC_after->Draw();
    h_Epz_REC->SetMarkerColor(kBlue);
    h_Epz_REC->SetMarkerStyle(24);
    h_Epz_REC->Draw("PEsame");
    //h_Epz_beforeCut->SetMarkerStyle(6);
    //h_Epz_beforeCut->SetMarkerColor(kRed);
    //h_Epz_beforeCut->Draw("PEsame");
    //h_Epz_REC_cut->SetMarkerStyle(29);
    //h_Epz_REC_cut->SetMarkerColor(kRed);
    //h_Epz_REC_cut->Draw("PEsame");
    //h_Epz_afterCut->SetMarkerStyle(24);
    //h_Epz_afterCut->SetMarkerColor(kBlue);
    //h_Epz_afterCut->Draw("PEsame");
    TLegend *w10 = new TLegend(0.65,0.5,0.8,0.6);
	w10->AddEntry(h_Epz_MC_after, " MC", "L");
    w10->AddEntry(h_Epz_REC, " RECO", "P");
	//w10->AddEntry(h_Epz_beforeCut, "RECO bc", "P");
    //w10->AddEntry(h_Epz_afterCut, "RECO", "P");
	w10->Draw("same");

    c11->cd(2);
    gPad->SetLogz(1);
    //h_Epz_response->GetYaxis()->SetRangeUser(18,22);
    //h_Epz_response->GetXaxis()->SetRangeUser(15,25);
	h_Epz_response->GetYaxis()->SetTitle("(E_{MC}-p_{z,MC}) [GeV]");
	h_Epz_response->GetXaxis()->SetTitle("(E_{EMCal} - p_{z,trk}) [GeV]");
    h_Epz_response->Draw("colzsame");

    c11->cd(3);
    h_Epz_res->Draw("colzsame");

    //c11->cd(2);
    //gPad->SetLogz(1);
    //h_Epz_response_cut->GetYaxis()->SetRangeUser(2,21);
	//h_Epz_response_cut->GetYaxis()->SetTitle("(E_{MC}-p_{z,MC}) [GeV]");
	//h_Epz_response_cut->GetXaxis()->SetTitle("(E_{EMCal} - p_{z,trk}) [GeV]");
    //h_Epz_response_cut->Draw();

    //c11->cd(3);
    //gPad->SetLogy(1);
    //h_Epz_purity->GetYaxis()->SetTitle("bins");
    //h_Epz_purity->Draw("same");

    //c11->cd(4);
    //gPad->SetLogz(1);
    //h_Epz_migration->Draw("colzsame");

    //c11->cd(5);
    //gPad->SetLogy(1);
    //h_Epz_MC->GetYaxis()->SetTitle("counts");
	//h_Epz_MC->GetXaxis()->SetTitle("(E-p_{z}) [GeV]");
	//h_Epz_MC->SetLineColor(kBlack);
    //h_Epz_MC->Draw();
    //h_Epz_corrected->SetMarkerStyle(24);
    //h_Epz_corrected->SetMarkerColor(kBlue);
    //h_Epz_corrected->Draw("PEsame");
    //TLegend *w11_1 = new TLegend(0.57,0.2,0.75,0.35);
	//w11_1->AddEntry(h_Epz_MC, "MC", "L");
    //w11_1->AddEntry(h_Epz_corrected, "Corrected", "P");
	//w11_1->Draw("same");

    //c11->cd(6);
    //gPad->SetLogy(1);
    //h_Epz_acceptance->Draw("same");

    //c11->cd(3);  
    //gPad->SetTopMargin(0.08);  
    //TLatex* title11_3 = new TLatex();
    //title11_3->SetNDC(); 
    //title11_3->SetTextSize(0.05);
    //title11_3->SetTextAlign(22);  
    //title11_3->DrawLatex(0.5, 0.97, "Purity");  

    c11->cd(4);  
    gPad->SetTopMargin(0.08);  
    TLatex* title11_4 = new TLatex();
    title11_4->SetNDC(); 
    title11_4->SetTextSize(0.05);
    title11_4->SetTextAlign(22);  
    //title11_4->DrawLatex(0.5, 0.97, "Bin Migration"); 

    c11->cd(5);  
    gPad->SetTopMargin(0.08);  
    TLatex* title11_5 = new TLatex();
    title11_5->SetNDC(); 
    title11_5->SetTextSize(0.05);
    title11_5->SetTextAlign(22);  
    //title11_5->DrawLatex(0.5, 0.97, "Corrected Acceptance"); 

    c11->cd(6);  
    gPad->SetTopMargin(0.08);  
    TLatex* title11_6 = new TLatex();
    title11_6->SetNDC(); 
    title11_6->SetTextSize(0.05);
    title11_6->SetTextAlign(22);  
    //title11_6->DrawLatex(0.5, 0.97, "Acceptance"); 

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

    c11->cd(3);  
    gPad->SetTopMargin(0.08);  
    TLatex* title11_3 = new TLatex();
    title11_3->SetNDC(); 
    title11_3->SetTextSize(0.05);
    title11_3->SetTextAlign(22);  
    title11_3->DrawLatex(0.5, 0.97, "E-p_{z} Resolution (e'+HFS)");  
   
    c11->Print("./figures/plot_Epz.pdf");

// E-pz and E/p
    TCanvas* c1122 = new TCanvas("c1122","c1122",1,1,1200,800);
    c1122->Divide(2,2,0.01,0.01);
    c1122->cd(1);
    gPad->SetLogy(1);
    //h_Epz_MC->GetXaxis()->SetRangeUser(15,25);
	h_Epz_MC->GetYaxis()->SetTitle("counts");
	h_Epz_MC->GetXaxis()->SetTitle("(E - p_{z}) [GeV]");
	h_Epz_MC->SetLineColor(kBlack);
    h_Epz_MC->Draw();
    h_Epz_REC->SetMarkerColor(kBlue);
    h_Epz_REC->SetMarkerStyle(24);
    h_Epz_REC->Draw("PEsame");
    //h_Epz_REC_cut->SetMarkerStyle(24);
    //h_Epz_REC_cut->SetMarkerColor(kBlue);
    //h_Epz_REC_cut->Draw("PEsame");
    //h_Epz_afterCut->SetMarkerStyle(24);
    //h_Epz_afterCut->SetMarkerColor(kBlue);
    //h_Epz_afterCut->Draw("PEsame");
    TLegend *w1022 = new TLegend(0.65,0.7,0.75,0.85);
	w1022->AddEntry(h_Epz_MC, "MC", "L");
	w1022->AddEntry(h_Epz_REC, "RECO", "P");
	//w1022->AddEntry(h_Epz_REC_cut, "RECO", "P");
    //w10->AddEntry(h_Epz_afterCut, "RECO", "P");
	w1022->Draw("same");

    //c1122->cd(2);
    //gPad->SetLogz(1);
    //h_Epz_response_cut->GetXaxis()->SetRangeUser(15,25);
    //h_Epz_response_cut->GetYaxis()->SetRangeUser(2,21);
	//h_Epz_response_cut->GetYaxis()->SetTitle("(E_{MC}-p_{z,MC}) [GeV]");
	//h_Epz_response_cut->GetXaxis()->SetTitle("(E_{EMCal} - p_{z,trk}) [GeV]");
    //h_Epz_response_cut->Draw();

     c1122->cd(2);
    gPad->SetLogz(1);
    h_Epz_response->GetYaxis()->SetRangeUser(10,25);
    h_Epz_response->GetXaxis()->SetRangeUser(10,30);
	h_Epz_response->GetYaxis()->SetTitle("(E_{MC}-p_{z,MC}) [GeV]");
	h_Epz_response->GetXaxis()->SetTitle("(E_{EMCal} - p_{z,trk}) [GeV]");
    h_Epz_response->Draw();

    c1122->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title1122 = new TLatex();
    title1122->SetNDC(); 
    title1122->SetTextSize(0.05);
    title1122->SetTextAlign(22);  
    title1122->DrawLatex(0.5, 0.97, "E-p_{z} Truth vs. Reco");  

    c1122->cd(2);  
    gPad->SetTopMargin(0.08);  
    TLatex* title11_222 = new TLatex();
    title11_222->SetNDC(); 
    title11_222->SetTextSize(0.05);
    title11_222->SetTextAlign(22);  
    title11_222->DrawLatex(0.5, 0.97, "E-p_{z} Response");

    c1122->cd(3);
    gPad->SetLogy(1);
    h_EoverP_MC->SetTitle("E/|p| Truth vs. Reco");
	h_EoverP_MC->GetYaxis()->SetTitle("counts");
	h_EoverP_MC->GetXaxis()->SetTitle("E/|p|");
	h_EoverP_MC->SetLineColor(kBlack);
    h_EoverP_MC->Draw();
    h_EcalOverPtrk->SetMarkerColor(kBlue);
    h_EcalOverPtrk->SetMarkerStyle(24);
    h_EcalOverPtrk->Draw("PEsame");
    //h_EcalOverPtrk_cut->SetMarkerStyle(24);
    //h_EcalOverPtrk_cut->SetMarkerColor(kBlue);
    //h_EcalOverPtrk_cut->Draw("PEsame");
    //h_EoverP_MC_after->SetMarkerStyle(11);
    //h_EoverP_MC_after->SetMarkerColor(kOrange);
    //h_EoverP_MC_after->Draw("PEsame");
    TLegend *w1122 = new TLegend(0.22,0.7,0.4,0.85);
	w1122->AddEntry(h_EoverP_MC, "MC", "L");
	w1122->AddEntry(h_EcalOverPtrk, "RECO", "P");
    //w1122->AddEntry(h_EcalOverPtrk_cut, "RECO", "P");
    //w11->AddEntry(h_EoverP_MC_after, "MC (after cut)", "P");
	w1122->Draw("same");

    //c1122->cd(4);
    //gPad->SetLogz(1);
    //h_EoverP_response_cut->SetTitle("E/|p| Response (After Threshold Cut)");   
    //h_EoverP_response_cut->GetYaxis()->SetTitleOffset(1.5);
    //h_EoverP_response_cut->GetYaxis()->SetRangeUser(1,1.02);
	//h_EoverP_response_cut->GetYaxis()->SetTitle("E_{MC}/|p|_{MC}");
	//h_EoverP_response_cut->GetXaxis()->SetTitle("E_{EMCal}/|p|_{trk}");
	//h_EoverP_response_cut->SetLineColor(kBlack);
    //h_EoverP_response_cut->Draw();

    c1122->cd(4);
    gPad->SetLogz(1);
    h_EoverP_response->SetTitle("E/|p| Response (After Threshold Cut)");   
    h_EoverP_response->GetYaxis()->SetTitleOffset(1.5);
    h_EoverP_response->GetYaxis()->SetRangeUser(1,1.02);
	h_EoverP_response->GetYaxis()->SetTitle("E_{MC}/|p|_{MC}");
	h_EoverP_response->GetXaxis()->SetTitle("E_{EMCal}/|p|_{trk}");
	h_EoverP_response->SetLineColor(kBlack);
    h_EoverP_response->Draw();

    c1122->cd(3);  
    gPad->SetTopMargin(0.08); 
    TLatex* title112 = new TLatex();
    title112->SetNDC(); 
    title112->SetTextSize(0.05);
    title112->SetTextAlign(22);  
    title112->DrawLatex(0.5, 0.97, "E/|p| Truth vs. Reco");  

    c1122->cd(4);  
    gPad->SetTopMargin(0.08);  
    TLatex* title12_212 = new TLatex();
    title12_212->SetNDC(); 
    title12_212->SetTextSize(0.05);
    title12_212->SetTextAlign(22);  
    title12_212->DrawLatex(0.5, 0.97, "E/|p| Response");  
   
    c1122->Print("./figures/plot_Epz_EoverP.pdf");

// E/p
    TCanvas* c12 = new TCanvas("c12","c12",1,1,1600,600);
    c12->Divide(3,1,0.01,0.01);
    c12->cd(1);
    gPad->SetLogy(1);
    h_EoverP_MC->SetTitle("E/|p| Truth vs. Reco");
	h_EoverP_MC->GetYaxis()->SetTitle("counts");
	h_EoverP_MC->GetXaxis()->SetTitle("E/|p|");
	h_EoverP_MC->SetLineColor(kBlack);
    //h_EoverP_MC->Draw();
    h_EoverP_MC_after->SetLineColor(kBlack);
    h_EoverP_MC_after->Draw();
    h_EcalOverPtrk->SetMarkerColor(kBlue);
    h_EcalOverPtrk->SetMarkerStyle(24);
    //h_EcalOverPtrk->Draw("PEsame");
    //h_EcalOverPtrk_cut->SetMarkerStyle(24);
    //h_EcalOverPtrk_cut->SetMarkerColor(kBlue);
    //h_EcalOverPtrk_cut->Draw("PEsame");
    //h_EoverP_beforeCut->SetMarkerStyle(7);
    //h_EoverP_beforeCut->SetMarkerColor(kRed);
    //h_EoverP_beforeCut->Draw("PEsame");
    h_EoverP_afterCut->SetMarkerStyle(24);
    h_EoverP_afterCut->SetMarkerColor(kBlue);
    h_EoverP_afterCut->Draw("PEsame");
    TLegend *w11 = new TLegend(0.22,0.7,0.4,0.85);
	w11->AddEntry(h_EoverP_MC_after, " MC", "L");
	w11->AddEntry(h_EoverP_afterCut, " RECO", "P");
    //w11->AddEntry(h_EcalOverPtrk_cut, "RECO", "P");
    //w11->AddEntry(h_EoverP_beforeCut, "RECO before cut", "P");
	w11->Draw("same");

    //c12->cd(2);
    //gPad->SetLogy(0);
    //h_EoverP_beforeCut->SetTitle("E/|p|");
	//h_EoverP_beforeCut->GetYaxis()->SetTitle("counts");
	//h_EoverP_beforeCut->GetXaxis()->SetTitle("E/|p|");
    //h_EoverP_beforeCut->SetMarkerStyle(7);
    //h_EoverP_beforeCut->SetMarkerColor(kRed);
    //h_EoverP_beforeCut->Draw("PEsame");
    h_EoverP_MC_after->SetMarkerStyle(11);
    h_EoverP_MC_after->SetMarkerColor(kOrange);
    //h_EoverP_MC_after->Draw("PEsame");
    TLegend *w223311 = new TLegend(0.22,0.2,0.4,0.3);
	//w11->AddEntry(h_EoverP_MC, "MC", "L");
	//w223311->AddEntry(h_EcalOverPtrk, "RECO ac", "P");
    //w11->AddEntry(h_EcalOverPtrk_cut, "RECO", "P");
    w223311->AddEntry(h_EoverP_beforeCut, "RECO before cut", "P");
	//w223311->Draw("same");


    c12->cd(2);
    gPad->SetLogz(1);    
    h_EoverP_response->GetYaxis()->SetTitleOffset(1.5);
    h_EoverP_response->SetTitle("E/|p| Response");   
    //h_EoverP_response->GetYaxis()->SetTitleOffset(1.5);
    //h_EoverP_response->GetYaxis()->SetRangeUser(1, 1.2);
	h_EoverP_response->GetYaxis()->SetTitle("E_{MC}/|p|_{MC}");
	h_EoverP_response->GetXaxis()->SetTitle("E_{EMCal}/|p|_{trk}");
	h_EoverP_response->SetLineColor(kBlack);
    h_EoverP_response->Draw();

    c12->cd(3);
    h_EoverP_res->Draw("colzsame");

    //c12->cd(2);
    //gPad->SetLogz(1);
    //h_EoverP_response_cut->SetTitle("E/|p| Response (After Threshold Cut)");   
    //h_EoverP_response_cut->GetYaxis()->SetTitleOffset(1.5);
    //h_EoverP_response_cut->GetYaxis()->SetRangeUser(1,1.02);
	//h_EoverP_response_cut->GetYaxis()->SetTitle("E_{MC}/|p|_{MC}");
	//h_EoverP_response_cut->GetXaxis()->SetTitle("E_{EMCal}/|p|_{trk}");
	//h_EoverP_response_cut->SetLineColor(kBlack);
    //h_EoverP_response_cut->Draw();

    //c12->cd(3);
    //gPad->SetLogy(1);
    //h_EoverP_purity->GetYaxis()->SetTitle("bins");
    //h_EoverP_purity->Draw("same");

    //c12->cd(3);
    //gPad->SetLogz(1);
    //h_EoverP_migration->Draw("colzsame");

    //c12->cd(5);
    //gPad->SetLogy(1);
    //h_EoverP_MC->GetYaxis()->SetTitle("counts");
	//h_EoverP_MC->GetXaxis()->SetTitle("E/|p|");
	//h_EoverP_MC->SetLineColor(kBlack);
    //h_EoverP_MC->Draw();
    //h_EoverP_corrected->SetMarkerStyle(24);
    //h_EoverP_corrected->SetMarkerColor(kBlue);
    //h_EoverP_corrected->Draw("PEsame");
    //TLegend *w12_1 = new TLegend(0.57,0.2,0.75,0.35);
	//w12_1->AddEntry(h_EoverP_MC, "MC", "L");
    //w12_1->AddEntry(h_EoverP_corrected, "Corrected", "P");
	//w12_1->Draw("same");

    //c12->cd(6);
    //gPad->SetLogy(1);
    //h_EoverP_acceptance->Draw("same");

    //c12->cd(3);  
    //gPad->SetTopMargin(0.08);  
    //TLatex* title12_3 = new TLatex();
    //title12_3->SetNDC(); 
    //title12_3->SetTextSize(0.05);
    //title12_3->SetTextAlign(22);  
    //title12_3->DrawLatex(0.5, 0.97, "Purity");  

    c12->cd(4);  
    gPad->SetTopMargin(0.08);  
    TLatex* title12_4 = new TLatex();
    title12_4->SetNDC(); 
    title12_4->SetTextSize(0.05);
    title12_4->SetTextAlign(22);  
    //title12_4->DrawLatex(0.5, 0.97, "Bin Migration"); 

    c12->cd(5);  
    gPad->SetTopMargin(0.08);  
    TLatex* title12_5 = new TLatex();
    title12_5->SetNDC(); 
    title12_5->SetTextSize(0.05);
    title12_5->SetTextAlign(22);  
    //title12_5->DrawLatex(0.5, 0.97, "Corrected Acceptance"); 

    c12->cd(6);  
    gPad->SetTopMargin(0.08);  
    TLatex* title12_6 = new TLatex();
    title12_6->SetNDC(); 
    title12_6->SetTextSize(0.05);
    title12_6->SetTextAlign(22);  
    //title12_6->DrawLatex(0.5, 0.97, "Acceptance"); 

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

    c12->cd(3);  
    gPad->SetTopMargin(0.08);  
    TLatex* title12_3 = new TLatex();
    title12_3->SetNDC(); 
    title12_3->SetTextSize(0.05);
    title12_3->SetTextAlign(22);  
    title12_3->DrawLatex(0.5, 0.97, "E/|p| Resolution");  

    c12->Print("./figures/plot_EoverP.pdf");

// t distribution
    TCanvas* c14 = new TCanvas("c14","c14",1,1,1200,800);
    c14->Divide(1,1,0.01,0.01);
    c14->cd(1);
    gPad->SetLogy(1);
	h_t_MC->GetYaxis()->SetTitle("counts");
    h_t_MC->GetXaxis()->SetLabelSize(0.03);
	h_t_MC->GetXaxis()->SetTitle("|t| [GeV/c]^{2}");
	h_t_MC->SetLineColor(kBlack);
    h_t_MC->Draw();
    h_t_REC_wRES_cut_pi12->Scale((M_PI/2)/(M_PI/12));
    h_t_REC_wRES_cut_pi12->SetMarkerStyle(30);
    h_t_REC_wRES_cut_pi12->SetMarkerSize(1.5);
    h_t_REC_wRES_cut_pi12->SetMarkerColor(kRed);
    h_t_REC_wRES_cut_pi12->SetLineColor(kRed);
    h_t_REC_wRES_cut_pi12->Draw("PEsame");
    h_t_REC_EEMC_cut->SetMarkerStyle(24); // method L RECO
    h_t_REC_EEMC_cut->SetMarkerSize(1.5);
	h_t_REC_EEMC_cut->SetMarkerColor(kBlue);
    h_t_REC_EEMC_cut->Draw("PEsame");
    TLegend *w13 = new TLegend(0.3,0.7,0.65,0.85);
	w13->AddEntry(h_t_MC, Form(" MC %.1f evts",h_t_MC->GetEntries()), "L");
	w13->AddEntry(h_t_REC_EEMC_cut, Form(" Method L reco %.1f evts",h_t_REC_EEMC_cut->GetEntries()), "P");    
    w13->AddEntry(h_t_REC_wRES_cut_pi12, Form(" #theta_{max}= #pi/12 %.1f evts",h_t_REC_wRES_cut_pi12->GetEntries()), "P");    
	w13->Draw("same");

    c14->cd(1);  
    gPad->SetTopMargin(0.08);  
    TLatex* title14_3 = new TLatex();
    title14_3->SetNDC(); 
    title14_3->SetTextSize(0.05);
    title14_3->SetTextAlign(22);  
    title14_3->DrawLatex(0.5, 0.97, "|t| Reconstruction Methods");

    c14->Print("./figures/plot_t_distribution.pdf");


// t DIS distribution 
  /*  TCanvas* c14_213 = new TCanvas("c14_213","c14_213",1,1,1200,800);
    TFile* DIS = TFile::Open("BeAGLE_DIS_noVetoes_10x100_eAu.root","READ");
    TH1D* h_DIS = (TH1D*) DIS->Get("h_t_MC");
    TFile* DIS_vetoes = TFile::Open("BeAGLE_DIS_vetoes_10x100_eAu.root","READ");
    TH1D* h_DIS_vetoes = (TH1D*) DIS_vetoes->Get("h_t_REC_wRES_cut_pi12");
    //TFile* DIS_vetoes_noDetectors = TFile::Open("merged_Vetoes_DIS_noDetectors.root");
    //TH1D* h_DIS_vetoes_noDetector = (TH1D*) DIS_vetoes_noDetectors->Get("h_t_REC_wRES_cut_pi12_cut");
    //TFile* DIS_vetoes_noDetectors_noEta = TFile::Open("merged_Vetoes_DIS_noDetectors_noEta.root");
    //TH1D* h_DIS_vetoes_noDetector_noEta = (TH1D*) DIS_vetoes_noDetectors_noEta->Get("h_t_REC_wRES_cut_pi12_cut");
    //TFile* DIS_vetoes_noDetectors_noExtraHFS = TFile::Open("merged_Vetoes_DIS_noDetectors_noExtraHFS.root");
    //TH1D* h_DIS_vetoes_noDetector_noExtraHFS = (TH1D*) DIS_vetoes_noDetectors_noExtraHFS->Get("h_t_REC_wRES_cut_pi12_cut");
    TFile* DIS_vetoes_onlyDetectors = TFile::Open("BeAGLE_DIS_onlyDetectorVetoes_10x100_eAu.root");
    TH1D* h_DIS_vetoes_onlyDetectors = (TH1D*) DIS_vetoes_onlyDetectors->Get("h_t_REC_wRES_cut_pi12");

    c14_213->Divide(1,1,0.01,0.01);

    // get counts
		double integral_MC = h_t_MC->Integral();
		double integral_REC = h_t_REC_EEMC->Integral();
		double integral_REC_wRES_cut_pi12 = h_t_REC_wRES_cut_pi12_cut->Integral();
		double integral_DIS = h_DIS->Integral();
		double integral_DIS_vetoes = h_DIS_vetoes->Integral();
        double dN_dt = 0.002;
        double sigma_coherent = 447; // nb
        double sigma_incoherent = 247; // nb
        double sigma_DIS = 67236; // nb
        double nEvents_coherent = h_t_MC->GetEntries();//6.2e6;
        double nEvents_incoherent = 10000000-nEvents_coherent; //3.8e6; total events 10M
        double nEvents_DIS = h_DIS->GetEntries();//1.1e6;
        double nEvents_DIS_vetoes = h_DIS_vetoes->GetEntries();//14837;
        int rebin_width = 2;        
		// Draw histograms with different angle cuts
        c14_213->cd(1);
        gPad->SetLogy(1);
        h_t_MC->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
        h_t_MC->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
        h_t_MC->SetLineColor(kBlack);
		h_t_MC->Draw("HISTsame");
		h_t_REC_EEMC->SetMarkerStyle(24); // method L RECO
		h_t_REC_EEMC->SetMarkerColor(kBlue);
		//h_t_REC_EEMC->Draw("PEsame");
        h_t_REC_wRES_cut_pi12->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));		
		h_t_REC_wRES_cut_pi12->SetMarkerStyle(30);
		h_t_REC_wRES_cut_pi12->SetMarkerColor(kRed);
        h_t_REC_wRES_cut_pi12->SetLineColor(kRed);
		h_t_REC_wRES_cut_pi12->Draw("Psame");
        h_DIS->Rebin(rebin_width);
        h_DIS->Scale(sigma_DIS*(1/nEvents_DIS)*(1/dN_dt)/rebin_width);
        //h_DIS->Scale((sigma_DIS)*(1/nEvents_DIS)*(dN_dt));
        //h_DIS->Scale(1/integral_DIS);
        //h_DIS->Scale(sigma_DIS/sigma_coherent);
        h_DIS->SetMarkerStyle(29);
        h_DIS->SetMarkerColor(kOrange);
        h_DIS->SetLineColor(kOrange);
        h_DIS->Draw("Psame");
        h_DIS_vetoes->Rebin(rebin_width); 
        h_DIS_vetoes->Scale((sigma_DIS)*(1/nEvents_DIS)*(1/dN_dt)/rebin_width);
        h_DIS_vetoes->SetMarkerStyle(21);
        h_DIS_vetoes->SetMarkerColor(kMagenta);
        h_DIS_vetoes->SetLineColor(kMagenta);
        //h_DIS_vetoes->SetFillColorAlpha(kMagenta, 0.3);
        //h_DIS_vetoes->Draw("E2same");
        h_DIS_vetoes->Draw("Psame");
        //h_DIS_vetoes_noDetector_noEta->Rebin(rebin_width);
        //h_DIS_vetoes_noDetector_noEta->Scale((sigma_DIS)*(1/nEvents_DIS)*(1/dN_dt)/rebin_width);
        //h_DIS_vetoes_noDetector_noEta->Sumw2();
        //h_DIS_vetoes_noDetector_noEta->SetMarkerStyle(21);
        //h_DIS_vetoes_noDetector_noEta->SetMarkerColor(kMagenta);
        //h_DIS_vetoes_noDetector_noEta->SetLineColor(kMagenta);
        //h_DIS_vetoes_noDetector_noEta->Draw("Psame");
        //h_DIS_vetoes_noDetector->Rebin(rebin_width);
        //h_DIS_vetoes_noDetector->Scale((sigma_DIS)*(1/nEvents_DIS)*(1/dN_dt)/rebin_width);
        //h_DIS_vetoes_noDetector->Sumw2();
        //h_DIS_vetoes_noDetector->SetMarkerStyle(21);
        //h_DIS_vetoes_noDetector->SetMarkerColor(kMagenta);
        //h_DIS_vetoes_noDetector->SetLineColor(kMagenta);
        //h_DIS_vetoes_noDetector->Draw("Psame");
        //h_DIS_vetoes_noDetector_noExtraHFS->Rebin(rebin_width);
        //h_DIS_vetoes_noDetector_noExtraHFS->Scale((sigma_DIS)*(1/nEvents_DIS)*(1/dN_dt)/rebin_width);
        //h_DIS_vetoes_noDetector_noExtraHFS->Sumw2();
        //h_DIS_vetoes_noDetector_noExtraHFS->SetMarkerStyle(21);
        //h_DIS_vetoes_noDetector_noExtraHFS->SetMarkerColor(kGreen+2);
        //h_DIS_vetoes_noDetector_noExtraHFS->SetLineColor(kGreen+2);
        //h_DIS_vetoes_noDetector_noExtraHFS->Draw("Psame");
        h_DIS_vetoes_onlyDetectors->Rebin(rebin_width);
        h_DIS_vetoes_onlyDetectors->Scale((sigma_DIS)*(1/nEvents_DIS)*(1/dN_dt)/rebin_width);
        h_DIS_vetoes_onlyDetectors->SetMarkerStyle(26);
        h_DIS_vetoes_onlyDetectors->SetMarkerColor(kBlue);
        h_DIS_vetoes_onlyDetectors->SetLineColor(kBlue);
        h_DIS_vetoes_onlyDetectors->Draw("PEsame");
        TLegend *w14_213 = new TLegend(0.25,0.67,0.7,0.9);
	w14_213->AddEntry(h_t_MC,Form("#phi: MC %.1f evts",h_t_MC->GetEntries()), "L");
    w14_213->AddEntry(h_t_REC_wRES_cut_pi12,Form("#phi #theta_{max}= #pi/12 %.1f evts",h_t_REC_wRES_cut_pi12->GetEntries()), "P");    
    w14_213->AddEntry(h_DIS,Form("DIS: no vetoes %.1f evts",h_DIS->GetEntries()), "P");
    //w14_213->AddEntry(h_DIS_vetoes_noDetector_noEta,Form("DIS: HFS cuts %.1f evts",h_DIS_vetoes_noDetector_noEta->GetEntries()), "P");
    //w14_213->AddEntry(h_DIS_vetoes_noDetector_noExtraHFS,Form("DIS: #eta cuts %.1f evts",h_DIS_vetoes_noDetector_noExtraHFS->GetEntries()), "P");
    //w14_213->AddEntry(h_DIS_vetoes_onlyDetectors,Form("DIS: Only detectors %.1f evts",h_DIS_vetoes_onlyDetectors->GetEntries()), "P");
    w14_213->AddEntry(h_DIS_vetoes,Form("DIS: after vetoes %.1f evts",h_DIS_vetoes->GetEntries()), "P");
    //w14_213->AddEntry(h_DIS_vetoes_noDetector,Form("DIS: #eta + HFS cuts %.1f evts",h_DIS_vetoes_noDetector->GetEntries()), "P");


	w14_213->Draw("same");	

    c14_213->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title14_213 = new TLatex();
    title14_213->SetNDC(); 
    title14_213->SetTextSize(0.05);
    title14_213->SetTextAlign(22);  
    title14_213->DrawLatex(0.5, 0.97, "|t| Distribution");  

    c14_213->Print("./figures/plot_t_DIS.pdf");
*/
// t with background incoherent and DIS distribution
   /* TCanvas* c142_213 = new TCanvas("c142_213","c142_213",1,1,800,600);
    c142_213->Divide(1,1,0.01,0.01);

    const int incoherent_bins = 91;
    double incoherent_min = 0.001, incoherent_max = 0.18; 
    double incoherent_xvals[incoherent_bins] = {0.0010133167966835457,0.0030211660950511393,0.0030211660950511393,0.003221951024887905,0.0052298003232554985,0.005631370182929016,
        0.008241574270806892,0.008241574270806892,0.008241574270806892,0.010450208499011251,0.011855703007868562,0.013863552306236163,0.016473756394114032,0.018481605692481633,
        0.020690239920685992,0.022898874148890344,0.025308293306931462,0.02912320697382989,0.031733411061707774,0.034945969939095926,0.037355389097137044,0.03936323839550463,
        0.04117030276403547,0.042977367132566315,0.04538678629060742,0.04759542051881178,0.04960326981717938,0.051410334185710206,0.05361896841391458,0.055626817712282166,
        0.057835451940486525, 0.05964251630901737,0.06185115053722171,0.06405978476542606, 0.0658668491339569,0.06827626829199804,0.0712880422395494,0.0712880422395494,
        0.07309510660808025,0.07450060111693757,0.07691002027497867,0.07911865450318303,0.08092571887171388,0.08293356817008146,0.08735083662649018,0.08735083662649018,
        0.08554377225795934,0.0893586859248578,0.09196889001273567,0.09477987903045029,0.09638615846914436,0.0981932228376752,0.10000028720620605,0.1022089214344104,
        0.10421677073277802,0.10602383510130883,0.10863403918918671,0.11044110355771755,0.11305130764559543,0.11485837201412627,0.11706700624233063,0.11927564047053499,
        0.1210827048390658,0.12349212399710696,0.12529918836563775,0.12811017738335242,0.12991724175188324,0.1321258759800876,0.1341337252784552,0.13634235950665954,0.13835020880502713,
        0.14035805810339474,0.14236590740176236,0.14417297177029317,0.14618082106866076,0.14818867036702837,0.15039730459523273,0.15200358403392678,0.15461378812180468,
        0.1566216374201723,0.1584287017887031,0.1604365510870707,0.16264518531527508,0.16465303461364264,0.16706245377168377,0.1690703030700514,0.17067658250874546,
        0.1728852167369498,0.17509385096515417,0.17710170026352176,0.17910954956188935};
    double incoherent_yvals[incoherent_bins] = {26.4259758563496,26.4259758563496,32.1883134302757,39.9882250899722,39.9882250899722,50.6678823310276,50.6678823310276,
        50.6678823310276,62.9457858323086,61.7163084604715,69.4704726439587,75.1739081002721,75.1739081002721,82.9661089645002,82.9661089645002,82.9661089645002,82.9661089645002,
        82.9661089645002,82.9661089645002,84.6189128366981,91.5660155319823,88.0239551889403,95.2506074330171,91.5660155319823,91.5660155319823,91.5660155319823,
        95.2506074330171,93.3901418761231,89.7775186113819,91.5660155319823,91.5660155319823,89.7775186113819,89.7775186113819,89.7775186113819,89.7775186113819,91.5660155319823,
        88.0239551889403,88.0239551889403,88.0239551889403,89.7775186113819,88.0239551889403,84.6189128366981,84.6189128366981,84.6189128366981,84.6189128366981,84.6189128366981,
        84.6189128366981,84.6189128366981,84.6189128366981,84.6189128366981,84.6189128366981,84.6189128366981,82.9661089645002,79.7567199496749,78.1988859956862,
        81.3455881901156,79.7567199496749,78.1988859956862,78.1988859956862,78.1988859956862,76.671480156466,78.1988859956862,76.671480156466,76.671480156466,
        75.1739081002721,73.7055871040409,76.671480156466,75.1739081002721,75.1739081002721,75.1739081002721,72.2659458266437,72.2659458266437,69.4704726439587,
        68.1135529868728,69.4704726439587,69.4704726439587,69.4704726439587,68.1135529868728,68.1135529868728,64.1997561565852,66.783137121768,66.783137121768,66.783137121768,
        64.1997561565852,64.1997561565852,64.1997561565852,61.7163084604715,61.7163084604715,65.4787073680391,62.9457858323086,60.5108456368344};

    TGraph *incoherent_hist = new TGraph(incoherent_bins, incoherent_xvals, incoherent_yvals);
	
		// Draw histograms with different angle cuts
        c142_213->cd(1);
        gPad->SetLogy(1);
        h_t_MC->GetYaxis()->SetTitle("counts");
        h_t_MC->GetXaxis()->SetRangeUser(0,0.18);
        h_t_MC->Scale(1/incoherent_hist->Integral());
		h_t_MC->SetLineColor(kBlack);
        //h_t_MC->Scale((sigma_coherent)*(1/nEvents_coherent)*(1/dN_dt));
		h_t_MC->Draw("HISTsame");
        //h_t_REC_EEMC->Scale((sigma_coherent)*(1/nEvents_coherent)*(1/dN_dt));
        h_t_REC_EEMC->Scale(h_t_MC->Integral()/h_t_REC_EEMC->Integral());
		h_t_REC_EEMC->SetMarkerStyle(24); // method L RECO
		h_t_REC_EEMC->SetMarkerColor(kBlue);
		h_t_REC_EEMC->Draw("Psame");
        //incoherent_hist->Scale(sigma_incoherent*(1/nEvents_incoherent)*(1/dN_dt));
        //incoherent_hist->Scale(h_t_MC->Integral()/incoherent_hist->Integral(),"width");
        incoherent_hist->SetMarkerStyle(21);
        incoherent_hist->SetMarkerColor(kMagenta);
        incoherent_hist->Draw("Psame");
        h_DIS->Scale(1/incoherent_hist->Integral());
        //h_DIS->Scale(1/rebin_width);
         h_DIS->SetMarkerStyle(29);
        h_DIS->SetMarkerColor(kOrange);
        h_DIS->Draw("Psame");
        TLegend *w142_213 = new TLegend(0.55,0.65,0.75,0.9);
	w142_213->AddEntry(h_t_MC, " #phi_{MC}", "L");
    w142_213->AddEntry(h_t_REC_EEMC, " #phi_{L}", "P");
    w142_213->AddEntry(incoherent_hist, " Incoherent", "P");
        w142_213->AddEntry(h_DIS,"DIS","P");//Form("DIS %f evts",h_DIS->GetEntries()), "P");
	w142_213->Draw("same");

    TLegend *w142_2131 = new TLegend(0.25,0.75,0.3,0.9);
    w142_2131->AddEntry((TObject*)0,"eAu 10x100 GeV", "");
    w142_2131->AddEntry((TObject*)0, " #phi #rightarrow K^{+}K^{-}", "");
	w142_2131->Draw("same");

    c142_213->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title142_213 = new TLatex();
    title142_213->SetNDC(); 
    title142_213->SetTextSize(0.05);
    title142_213->SetTextAlign(22);  
    title142_213->DrawLatex(0.5, 0.97, "|t| Distribution");  

    c142_213->Print("./figures/plot_t_incoherent.pdf");
*/
// t QA
    /*TCanvas* c15 = new TCanvas("c15","c15",1,1,1800,800);
    c15->Divide(3,2,0.01,0.01);

    c15->cd(2);
    gPad->SetLogz(1);
    h_t_res_EEMC->GetYaxis()->SetRangeUser(-1,1);
    h_t_res_EEMC->GetXaxis()->SetTitle("|t|_{MC} [GeV/c]^{2}");
    h_t_res_EEMC->GetYaxis()->SetTitle("(|t|_{MC}-|t|_{EEMC})/|t|_{MC}");
    h_t_res_EEMC->Draw("colzsame");

    c15->cd(1);
    gPad->SetLogy(1);
    h_t_MC->GetXaxis()->SetLabelSize(0.03);
	h_t_MC->GetYaxis()->SetTitle("counts");
	h_t_MC->GetXaxis()->SetTitle("|t| [GeV/c]^{2}");
	h_t_MC->SetLineColor(kBlack);
    h_t_MC->Draw();
    //h_t_MC_before->SetMarkerStyle(21);
    //h_t_MC_before->SetMarkerColor(kCyan);
    //h_t_MC_before->Draw("PEsame");
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

    //c15->cd(2);
    //gPad->SetLogz(1);
    //h_t_res_EEMC_cut_percent->GetYaxis()->SetRangeUser(-1,1);
    //h_t_res_EEMC_cut_percent->GetXaxis()->SetLabelSize(0.03);
    //h_t_res_EEMC_cut_percent->GetXaxis()->SetTitle("|t|_{MC} [GeV/c]^{2}");
    //h_t_res_EEMC_cut_percent->GetYaxis()->SetTitle("(|t|_{MC}-|t|_{EEMC})/|t|_{MC}");
    //h_t_res_EEMC_cut_percent->Draw("colzsame");

    
    //h_t_response_EEMC_cut_0: mean of each slice
    //h_t_response_EEMC_cut_1: sigma (width) of each slice
    //h_t_response_EEMC_cut_2: chi2 of each fit
    

    //c15->cd(2);
    //h_t_res_EEMC_cut_percent->FitSlicesX(0,0,1000);
    //TH1D* h_sigma15 = (TH1D*)gDirectory->Get("h_t_res_EEMC_cut_percent_1");
    //h_sigma15->Draw("same");

    //c14->cd(5);
    //gPad->SetLogz(1);
    //h_t_response_EEMC->GetXaxis()->SetTitle("|t|_{EEMC} [GeV/c]^{2}");
    //h_t_response_EEMC->GetYaxis()->SetTitle("|t|_{MC} [GeV/c]^{2}");
    //h_t_response_EEMC->Draw("colzsame");

    //c15->cd(3);
    //gPad->SetLogz(1);
    //h_t_response_EEMC_cut->GetYaxis()->SetLabelSize(0.03); 
    //h_t_response_EEMC_cut->GetYaxis()->SetTitleOffset(1.4);
    //h_t_response_EEMC_cut->GetXaxis()->SetLabelSize(0.03);
    //h_t_response_EEMC_cut->GetXaxis()->SetTitle("|t|_{EEMC} [GeV/c]^{2}");
    //h_t_response_EEMC_cut->GetYaxis()->SetTitle("|t|_{MC} [GeV/c]^{2}");
    //h_t_response_EEMC_cut->Draw("colzsame");
    //h_t_response_EEMC_cut->Draw();

    c15->cd(3);
    gPad->SetLogz(1);
    h_t_response_EEMC->GetYaxis()->SetLabelSize(0.03); 
    h_t_response_EEMC->GetYaxis()->SetTitleOffset(1.4);
    h_t_response_EEMC->GetXaxis()->SetLabelSize(0.03);
    h_t_response_EEMC->GetXaxis()->SetTitle("|t|_{EEMC} [GeV/c]^{2}");
    h_t_response_EEMC->GetYaxis()->SetTitle("|t|_{MC} [GeV/c]^{2}");
    h_t_response_EEMC->Draw("colzsame");

    c15->cd(4);
    gPad->SetLogz(1);
    //h_t_response_cut->Draw("colzsame");

    c15->cd(5);
    gPad->SetLogz(1);
    //h_t_res->Draw("colzsame");

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
    h_t_REC_2d_wRES_cut->GetYaxis()->SetLabelSize(0.03); 
    h_t_REC_2d_wRES_cut->GetYaxis()->SetTitleOffset(1.4);
    h_t_REC_2d_wRES_cut->GetXaxis()->SetLabelSize(0.03);
    h_t_REC_2d_wRES_cut->GetXaxis()->SetTitle("#sqrt{|t|_{x}} [GeV/c]");
    h_t_REC_2d_wRES_cut->GetYaxis()->SetTitle("#sqrt{|t|_{y}} [GeV/c]");
    h_t_REC_2d_wRES_cut->Draw("colzsame");

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
*/
// migration, acceptance, purity
// for t
    TCanvas* c17 = new TCanvas("c17","c17",1,1,1400,800);
    c17->Divide(3,2,0.01,0.01);
    c17->cd(1);
    h_t_migration->GetXaxis()->SetNdivisions(505);
    h_t_migration->Draw("colzsame");

    TH2D* h_t_mig_norm = (TH2D*)h_t_migration->Clone("h_t_mig_norm");
    int nx = h_t_mig_norm->GetNbinsX();
    int ny = h_t_mig_norm->GetNbinsY();
    for (int i=1; i<=nx; ++i) 
    {
        double rowSum = 0;
        for (int j=1; j<=ny; ++j) rowSum += h_t_mig_norm->GetBinContent(i,j);
        if (rowSum > 0) 
        {
            for (int j=1; j<=ny; ++j) 
            {
                double v = h_t_mig_norm->GetBinContent(i,j);
                h_t_mig_norm->SetBinContent(i,j, v / rowSum);
            }
        }
    }
    c17->cd(2);
    h_t_mig_norm->GetXaxis()->SetNdivisions(505);
    h_t_mig_norm->Draw("colzsame");

    for (int j=1; j<=ny; ++j) 
    {
        double recoColSum = 0;
        for (int i=1; i<=nx; ++i) recoColSum += h_t_migration->GetBinContent(i,j);
        double diag = h_t_migration->GetBinContent(j,j);
        if (recoColSum>0) h_t_purity->SetBinContent(j, diag / recoColSum);
    }
    c17->cd(3);
    h_t_purity->GetXaxis()->SetNdivisions(505);
    h_t_purity->GetXaxis()->SetTitle("|t|_{reco bin} [GeV/c]^{2}");
    h_t_purity->GetYaxis()->SetTitle("purity");
    h_t_purity->Draw("same");

    TH1D* h_t_stability = new TH1D("h_t_stability","Stability;|t|_{true bin};stability",nx,0,0.2);
    for (int i=1; i<=nx; ++i) 
    {
        double trueRowSum = 0;
        for (int j=1; j<=ny; ++j) trueRowSum += h_t_migration->GetBinContent(i,j);
        double diag = h_t_migration->GetBinContent(i,i);
        if (trueRowSum>0) h_t_stability->SetBinContent(i, diag / trueRowSum);
    }
    c17->cd(4);
    h_t_stability->GetXaxis()->SetNdivisions(505);
    h_t_stability->Draw("same");

    int nb = h_t_MC->GetNbinsX();
    for (int i=1; i<=nb; ++i) 
    {
        double Ntrue = h_t_MC_before->GetBinContent(i);   
        double Naccepted = 0;
        for (int j=1; j<=nb; ++j) Naccepted += h_t_migration->GetBinContent(i,j);
        double acceptance = (Ntrue>0) ? Naccepted / Ntrue : 0;
        h_t_acceptance->SetBinContent(i, acceptance);
        double Nreco_in_bin_i = h_t_REC_wRES_cut_pi12->GetBinContent(i);
        double corrected = (acceptance>0) ? Nreco_in_bin_i / acceptance : 0;
        h_t_corrected->SetBinContent(i, corrected);
    }
    c17->cd(5);
    h_t_acceptance->GetXaxis()->SetNdivisions(505);
    h_t_acceptance->Draw("same");

    c17->cd(6);
    gPad->SetLogy(1);
    h_t_corrected->GetXaxis()->SetNdivisions(505);
    h_t_corrected->Draw("same");

    c17->cd(1);  
    gPad->SetTopMargin(0.08);  
    TLatex* title17_1 = new TLatex();
    title17_1->SetNDC(); 
    title17_1->SetTextSize(0.05);
    title17_1->SetTextAlign(22);  
    title17_1->DrawLatex(0.5, 0.97, "Bin Migration");  

    c17->cd(2);  
    gPad->SetTopMargin(0.08);  
    TLatex* title17_2 = new TLatex();
    title17_2->SetNDC(); 
    title17_2->SetTextSize(0.05);
    title17_2->SetTextAlign(22);  
    title17_2->DrawLatex(0.5, 0.97, "Normalized Migration (Probability)");  

    c17->cd(3);  
    gPad->SetTopMargin(0.08);  
    TLatex* title17_3 = new TLatex();
    title17_3->SetNDC(); 
    title17_3->SetTextSize(0.05);
    title17_3->SetTextAlign(22);  
    title17_3->DrawLatex(0.5, 0.97, "Purity");  

    c17->cd(4);  
    gPad->SetTopMargin(0.08);  
    TLatex* title17_4 = new TLatex();
    title17_4->SetNDC(); 
    title17_4->SetTextSize(0.05);
    title17_4->SetTextAlign(22);  
    title17_4->DrawLatex(0.5, 0.97, "Stability");

    c17->cd(5);  
    gPad->SetTopMargin(0.08);  
    TLatex* title17_5 = new TLatex();
    title17_5->SetNDC(); 
    title17_5->SetTextSize(0.05);
    title17_5->SetTextAlign(22);  
    title17_5->DrawLatex(0.5, 0.97, "Acceptance");

    c17->cd(6);  
    gPad->SetTopMargin(0.08);  
    TLatex* title17_6 = new TLatex();
    title17_6->SetNDC(); 
    title17_6->SetTextSize(0.05);
    title17_6->SetTextAlign(22);  
    title17_6->DrawLatex(0.5, 0.97, "Bin by Bin Corrected");

    c17->Print("./figures/plot_t_migrationCorrectedPurityAcceptance.pdf");
    
// for energy
    TCanvas* c18 = new TCanvas("c18","c18",1,1,1400,800);
    c18->Divide(3,2,0.01,0.01);
    c18->cd(1);
    h_energy_migration->GetXaxis()->SetRangeUser(1,11);
    h_energy_migration->GetYaxis()->SetRangeUser(1,11);
    h_energy_migration->Draw("colzsame");

    // normalize give migration probability per true bin
    TH2D* h_energy_mig_norm = (TH2D*)h_energy_migration->Clone("h_energy_mig_norm");
    int nx_energy = h_energy_mig_norm->GetNbinsX();
    int ny_energy = h_energy_mig_norm->GetNbinsY();
    for (int i=1; i<=nx_energy; ++i) 
    {
        double rowSum = 0;
        for (int j=1; j<=ny_energy; ++j) rowSum += h_energy_mig_norm->GetBinContent(i,j);
        if (rowSum > 0) 
        {
            for (int j=1; j<=ny_energy; ++j) 
            {
                double v = h_energy_mig_norm->GetBinContent(i,j);
                h_energy_mig_norm->SetBinContent(i,j, v / rowSum);
            }
        }
    }
    c18->cd(2);
    h_energy_mig_norm->Draw("colzsame");

    // purity for reco bin j: fraction of events in reco bin j that came from same true bin j
    for (int j=1; j<=ny_energy; ++j) 
    {
        double recoColSum = 0;
        for (int i=1; i<=nx_energy; ++i) recoColSum += h_energy_migration->GetBinContent(i,j);
        double diag = h_energy_migration->GetBinContent(j,j);
        if (recoColSum>0) h_energy_purity->SetBinContent(j, diag / recoColSum);
    }
    c18->cd(3);
    h_energy_purity->GetXaxis()->SetTitle("E_{reco bin} [GeV]");
    h_energy_purity->GetYaxis()->SetTitle("purity");
    h_energy_purity->Draw("same");

    // stability for true bin i: fraction of true events reconstructed in same bin
        // values close to 1 in well reconstructed bins
    TH1D* h_energy_stability = new TH1D("h_energy_stability","Stability;E_{true bin};stability",nx_energy,0,20);
    for (int i=1; i<=nx_energy; ++i) 
    {
        double trueRowSum = 0;
        for (int j=1; j<=ny_energy; ++j) trueRowSum += h_energy_migration->GetBinContent(i,j);
        double diag = h_energy_migration->GetBinContent(i,i);
        if (trueRowSum>0) h_energy_stability->SetBinContent(i, diag / trueRowSum);
    }
    c18->cd(4);
    h_energy_stability->Draw("same");

    // acceptance: fraction of true events in a true bin that survives selection and reco cuts
        // y axis: acceptance fraction (0-1) per true bin
            // want values near 1 in efficient regions
    // bin by bin correction: estimate of true bin content
        // y axis: estimated true counts (yield) per bin after dividing reco/acceptance
    int nb_energy = h_energy_MC->GetNbinsX();
    for (int i=1; i<=nb_energy; ++i) 
    {
        double Ntrue = h_energy_MC->GetBinContent(i);   
        double Naccepted = 0;
        for (int j=1; j<=nb_energy; ++j) Naccepted += h_energy_migration->GetBinContent(i,j);
        double acceptance = (Ntrue>0) ? Naccepted / Ntrue : 0;
        h_energy_acceptance->SetBinContent(i, acceptance);
        double Nreco_in_bin_i = h_energy_REC_EEMC_after->GetBinContent(i);
        double corrected = (acceptance>0) ? Nreco_in_bin_i / acceptance : 0;
        h_energy_corrected->SetBinContent(i, corrected);
    }
    c18->cd(5);
    h_energy_acceptance->Draw("same");

    c18->cd(6);
    gPad->SetLogy(0);
    h_energy_corrected->Draw("same");

    c18->cd(1);  
    gPad->SetTopMargin(0.08);  
    TLatex* title18_1 = new TLatex();
    title18_1->SetNDC(); 
    title18_1->SetTextSize(0.05);
    title18_1->SetTextAlign(22);  
    title18_1->DrawLatex(0.5, 0.97, "Bin Migration");  

    c18->cd(2);  
    gPad->SetTopMargin(0.08);  
    TLatex* title18_2 = new TLatex();
    title18_2->SetNDC(); 
    title18_2->SetTextSize(0.05);
    title18_2->SetTextAlign(22);  
    title18_2->DrawLatex(0.5, 0.97, "Normalized Migration (Probability)");  

    c18->cd(3);  
    gPad->SetTopMargin(0.08);  
    TLatex* title18_3 = new TLatex();
    title18_3->SetNDC(); 
    title18_3->SetTextSize(0.05);
    title18_3->SetTextAlign(22);  
    title18_3->DrawLatex(0.5, 0.97, "Purity");  

    c18->cd(4);  
    gPad->SetTopMargin(0.08);  
    TLatex* title18_4 = new TLatex();
    title18_4->SetNDC(); 
    title18_4->SetTextSize(0.05);
    title18_4->SetTextAlign(22);  
    title18_4->DrawLatex(0.5, 0.97, "Stability");

    c18->cd(5);  
    gPad->SetTopMargin(0.08);  
    TLatex* title18_5 = new TLatex();
    title18_5->SetNDC(); 
    title18_5->SetTextSize(0.05);
    title18_5->SetTextAlign(22);  
    title18_5->DrawLatex(0.5, 0.97, "Acceptance");

    c18->cd(6);  
    gPad->SetTopMargin(0.08);  
    TLatex* title18_6 = new TLatex();
    title18_6->SetNDC(); 
    title18_6->SetTextSize(0.05);
    title18_6->SetTextAlign(22);  
    title18_6->DrawLatex(0.5, 0.97, "Bin by Bin Corrected");
        
    c18->Print("./figures/plot_energy_migrationCorrectedPurityAcceptance.pdf");

// for theta
    TCanvas* c19 = new TCanvas("c19","c19",1,1,1400,800);
    c19->Divide(3,2,0.01,0.01);
    c19->cd(1);
    h_theta_migration->GetXaxis()->SetRangeUser(2.7,3.14);
    h_theta_migration->GetYaxis()->SetRangeUser(2.7,3.14);
    h_theta_migration->Draw("colzsame");

    TH2D* h_theta_mig_norm = (TH2D*)h_theta_migration->Clone("h_theta_mig_norm");
    int nx_theta = h_theta_mig_norm->GetNbinsX();
    int ny_theta = h_theta_mig_norm->GetNbinsY();
    for (int i=1; i<=nx_theta; ++i) 
    {
        double rowSum = 0;
        for (int j=1; j<=ny_theta; ++j) rowSum += h_theta_mig_norm->GetBinContent(i,j);
        if (rowSum > 0) 
        {
            for (int j=1; j<=ny_theta; ++j) 
            {
                double v = h_theta_mig_norm->GetBinContent(i,j);
                h_theta_mig_norm->SetBinContent(i,j, v / rowSum);
            }
        }
    }
    c19->cd(2);
    h_theta_mig_norm->Draw("colzsame");

    for (int j=1; j<=ny_theta; ++j) 
    {
        double recoColSum = 0;
        for (int i=1; i<=nx_theta; ++i) recoColSum += h_theta_migration->GetBinContent(i,j);
        double diag = h_theta_migration->GetBinContent(j,j);
        if (recoColSum>0) h_theta_purity->SetBinContent(j, diag / recoColSum);
    }
    c19->cd(3);
    h_theta_purity->GetXaxis()->SetTitle("#theta_{reco bin}");
    h_theta_purity->GetYaxis()->SetTitle("purity");
    h_theta_purity->Draw("same");

    TH1D* h_theta_stability = new TH1D("h_theta_stability","Stability;#theta_{true bin};stability",nx_theta,0,3.14);
    for (int i=1; i<=nx_theta; ++i) 
    {
        double trueRowSum = 0;
        for (int j=1; j<=ny_theta; ++j) trueRowSum += h_theta_migration->GetBinContent(i,j);
        double diag = h_theta_migration->GetBinContent(i,i);
        if (trueRowSum>0) h_theta_stability->SetBinContent(i, diag / trueRowSum);
    }
    c19->cd(4);
    h_theta_stability->Draw("same");
    
    int nb_theta = h_theta_MC->GetNbinsX();
    for (int i=1; i<=nb_theta; ++i) 
    {
        double Ntrue = h_theta_MC->GetBinContent(i);   
        double Naccepted = 0;
        for (int j=1; j<=nb_theta; ++j) Naccepted += h_theta_migration->GetBinContent(i,j);
        double acceptance = (Ntrue>0) ? Naccepted / Ntrue : 0;
        h_theta_acceptance->SetBinContent(i, acceptance);
        double Nreco_in_bin_i = h_theta_REC_EEMC_after->GetBinContent(i);
        double corrected = (acceptance>0) ? Nreco_in_bin_i / acceptance : 0;
        h_theta_corrected->SetBinContent(i, corrected);
    }
    c19->cd(5);
    h_theta_acceptance->Draw("same");

    c19->cd(6);
    gPad->SetLogy(0);
    h_theta_corrected->Draw("same");

    c19->cd(1);  
    gPad->SetTopMargin(0.08);  
    TLatex* title19_1 = new TLatex();
    title19_1->SetNDC(); 
    title19_1->SetTextSize(0.05);
    title19_1->SetTextAlign(22);  
    title19_1->DrawLatex(0.5, 0.97, "Bin Migration");  

    c19->cd(2);  
    gPad->SetTopMargin(0.08);  
    TLatex* title19_2 = new TLatex();
    title19_2->SetNDC(); 
    title19_2->SetTextSize(0.05);
    title19_2->SetTextAlign(22);  
    title19_2->DrawLatex(0.5, 0.97, "Normalized Migration (Probability)");  

    c19->cd(3);  
    gPad->SetTopMargin(0.08);  
    TLatex* title19_3 = new TLatex();
    title19_3->SetNDC(); 
    title19_3->SetTextSize(0.05);
    title19_3->SetTextAlign(22);  
    title19_3->DrawLatex(0.5, 0.97, "Purity");  

    c19->cd(4);  
    gPad->SetTopMargin(0.08);  
    TLatex* title19_4 = new TLatex();
    title19_4->SetNDC(); 
    title19_4->SetTextSize(0.05);
    title19_4->SetTextAlign(22);  
    title19_4->DrawLatex(0.5, 0.97, "Stability");

    c19->cd(5);  
    gPad->SetTopMargin(0.08);  
    TLatex* title19_5 = new TLatex();
    title19_5->SetNDC(); 
    title19_5->SetTextSize(0.05);
    title19_5->SetTextAlign(22);  
    title19_5->DrawLatex(0.5, 0.97, "Acceptance");

    c19->cd(6);  
    gPad->SetTopMargin(0.08);  
    TLatex* title19_6 = new TLatex();
    title19_6->SetNDC(); 
    title19_6->SetTextSize(0.05);
    title19_6->SetTextAlign(22);  
    title19_6->DrawLatex(0.5, 0.97, "Bin by Bin Corrected");

    c19->Print("./figures/plot_theta_migrationCorrectedPurityAcceptance.pdf");

// for VM pT        
    TCanvas* c20 = new TCanvas("c20","c20",1,1,1400,800);
    c20->Divide(3,2,0.01,0.01);
    c20->cd(1);    
    h_VM_pt_migration->GetXaxis()->SetRangeUser(0,4);
    h_VM_pt_migration->GetYaxis()->SetRangeUser(0,4);
    h_VM_pt_migration->Draw("colzsame");

    TH2D* h_VM_pt_mig_norm = (TH2D*)h_VM_pt_migration->Clone("h_VM_pt_mig_norm");
    int nx_VMpt = h_VM_pt_mig_norm->GetNbinsX();
    int ny_VMpt = h_VM_pt_mig_norm->GetNbinsY();
    for (int i=1; i<=nx_VMpt; ++i) 
    {
        double rowSum = 0;
        for (int j=1; j<=ny_VMpt; ++j) rowSum += h_VM_pt_mig_norm->GetBinContent(i,j);
        if (rowSum > 0) 
        {
            for (int j=1; j<=ny_VMpt; ++j) 
            {
                double v = h_VM_pt_mig_norm->GetBinContent(i,j);
                h_VM_pt_mig_norm->SetBinContent(i,j, v / rowSum);
            }
        }
    }
    c20->cd(2);
    h_VM_pt_mig_norm->Draw("colzsame");

    for (int j=1; j<=ny_VMpt; ++j) 
    {
        double recoColSum = 0;
        for (int i=1; i<=nx_VMpt; ++i) recoColSum += h_VM_pt_migration->GetBinContent(i,j);
        double diag = h_VM_pt_migration->GetBinContent(j,j);
        if (recoColSum>0) h_VM_pt_purity->SetBinContent(j, diag / recoColSum);
    }
    c20->cd(3);
    h_VM_pt_purity->GetXaxis()->SetTitle("VM p_{T, reco bin}");
    h_VM_pt_purity->GetYaxis()->SetTitle("purity");
    h_VM_pt_purity->Draw("same");

    TH1D* h_VM_pt_stability = new TH1D("h_VM_pt_stability","Stability;VM p_{T, true bin};stability",nx_VMpt,0,10);
    for (int i=1; i<=nx_VMpt; ++i) 
    {
        double trueRowSum = 0;
        for (int j=1; j<=ny_VMpt; ++j) trueRowSum += h_VM_pt_migration->GetBinContent(i,j);
        double diag = h_VM_pt_migration->GetBinContent(i,i);
        if (trueRowSum>0) h_VM_pt_stability->SetBinContent(i, diag / trueRowSum);
    }
    c20->cd(4);
    h_VM_pt_stability->Draw("same");
    
    int nb_VMpt = h_VM_pt_MC_before->GetNbinsX();
    for (int i=1; i<=nb_VMpt; ++i) 
    {
        double Ntrue = h_VM_pt_MC_before->GetBinContent(i);   
        double Naccepted = 0;
        for (int j=1; j<=nb_VMpt; ++j) Naccepted += h_VM_pt_migration->GetBinContent(i,j);
        double acceptance = (Ntrue>0) ? Naccepted / Ntrue : 0;
        h_VM_pt_acceptance->SetBinContent(i, acceptance);
        double Nreco_in_bin_i = h_VM_pt_REC->GetBinContent(i);
        double corrected = (acceptance>0) ? Nreco_in_bin_i / acceptance : 0;
        h_VM_pt_corrected->SetBinContent(i, corrected);
    }
    c20->cd(5);
    h_VM_pt_acceptance->Draw("same");
    
    c20->cd(6);
    gPad->SetLogy(0);
    h_VM_pt_corrected->Draw("same");

    c20->cd(1);  
    gPad->SetTopMargin(0.08);  
    TLatex* title20_1 = new TLatex();
    title20_1->SetNDC(); 
    title20_1->SetTextSize(0.05);
    title20_1->SetTextAlign(22);  
    title20_1->DrawLatex(0.5, 0.97, "Bin Migration");  

    c20->cd(2);  
    gPad->SetTopMargin(0.08);  
    TLatex* title20_2 = new TLatex();
    title20_2->SetNDC(); 
    title20_2->SetTextSize(0.05);
    title20_2->SetTextAlign(22);  
    title20_2->DrawLatex(0.5, 0.97, "Normalized Migration (Probability)");  

    c20->cd(3);  
    gPad->SetTopMargin(0.08);  
    TLatex* title20_3 = new TLatex();
    title20_3->SetNDC(); 
    title20_3->SetTextSize(0.05);
    title20_3->SetTextAlign(22);  
    title20_3->DrawLatex(0.5, 0.97, "Purity");  

    c20->cd(4);  
    gPad->SetTopMargin(0.08);  
    TLatex* title20_4 = new TLatex();
    title20_4->SetNDC(); 
    title20_4->SetTextSize(0.05);
    title20_4->SetTextAlign(22);  
    title20_4->DrawLatex(0.5, 0.97, "Stability");

    c20->cd(5);  
    gPad->SetTopMargin(0.08);  
    TLatex* title20_5 = new TLatex();
    title20_5->SetNDC(); 
    title20_5->SetTextSize(0.05);
    title20_5->SetTextAlign(22);  
    title20_5->DrawLatex(0.5, 0.97, "Acceptance");

    c20->cd(6);  
    gPad->SetTopMargin(0.08);  
    TLatex* title20_6 = new TLatex();
    title20_6->SetNDC(); 
    title20_6->SetTextSize(0.05);
    title20_6->SetTextAlign(22);  
    title20_6->DrawLatex(0.5, 0.97, "Bin by Bin Corrected");
    
    c20->Print("./figures/plot_VM_migrationCorrectedPurityAcceptance.pdf");
    
// for VM e-pz
    TCanvas* c21 = new TCanvas("c21","c21",1,1,1400,800);
    c21->Divide(3,2,0.01,0.01);
    c21->cd(1);
    h_VM_Epz_migration->GetXaxis()->SetRangeUser(0,10);
    h_VM_Epz_migration->GetYaxis()->SetRangeUser(0,10);
    h_VM_Epz_migration->Draw("colzsame");

    TH2D* h_VM_Epz_mig_norm = (TH2D*)h_VM_Epz_migration->Clone("h_VM_Epz_mig_norm");
    int nx_VMEpz = h_VM_Epz_mig_norm->GetNbinsX();
    int ny_VMEpz = h_VM_Epz_mig_norm->GetNbinsY();
    for (int i=1; i<=nx_VMEpz; ++i) 
    {
        double rowSum = 0;
        for (int j=1; j<=ny_VMEpz; ++j) rowSum += h_VM_Epz_mig_norm->GetBinContent(i,j);
        if (rowSum > 0) 
        {
            for (int j=1; j<=ny_VMEpz; ++j) 
            {
                double v = h_VM_Epz_mig_norm->GetBinContent(i,j);
                h_VM_Epz_mig_norm->SetBinContent(i,j, v / rowSum);
            }
        }
    }
    c21->cd(2);
    h_VM_Epz_mig_norm->Draw("colzsame");

    for (int j=1; j<=ny_VMEpz; ++j) 
    {
        double recoColSum = 0;
        for (int i=1; i<=nx_VMEpz; ++i) recoColSum += h_VM_Epz_migration->GetBinContent(i,j);
        double diag = h_VM_Epz_migration->GetBinContent(j,j);
        if (recoColSum>0) h_VM_Epz_purity->SetBinContent(j, diag / recoColSum);
    }
    c21->cd(3);
    h_VM_Epz_purity->GetXaxis()->SetTitle("VM (E-p_{z})_{reco bin}");
    h_VM_Epz_purity->GetYaxis()->SetTitle("purity");
    h_VM_Epz_purity->Draw("same");

    TH1D* h_VM_Epz_stability = new TH1D("h_VM_Epz_stability","Stability;VM (E-p_{z})_{true bin};stability",nx_VMEpz,0,20);
    for (int i=1; i<=nx_VMEpz; ++i) 
    {
        double trueRowSum = 0;
        for (int j=1; j<=ny_VMEpz; ++j) trueRowSum += h_VM_Epz_migration->GetBinContent(i,j);
        double diag = h_VM_Epz_migration->GetBinContent(i,i);
        if (trueRowSum>0) h_VM_Epz_stability->SetBinContent(i, diag / trueRowSum);
    }
    c21->cd(4);
    h_VM_Epz_stability->Draw("same");
    
    int nb_VMEpz = h_VM_Epz_MC_before->GetNbinsX();
    for (int i=1; i<=nb_VMEpz; ++i) 
    {
        double Ntrue = h_VM_Epz_MC_before->GetBinContent(i);   
        double Naccepted = 0;
        for (int j=1; j<=nb_VMEpz; ++j) Naccepted += h_VM_Epz_migration->GetBinContent(i,j);
        double acceptance = (Ntrue>0) ? Naccepted / Ntrue : 0;
        h_VM_Epz_acceptance->SetBinContent(i, acceptance);
        double Nreco_in_bin_i = h_VM_Epz_REC->GetBinContent(i);
        double corrected = (acceptance>0) ? Nreco_in_bin_i / acceptance : 0;
        h_VM_Epz_corrected->SetBinContent(i, corrected);
    }
    c21->cd(5);
    h_VM_Epz_acceptance->Draw("same");
    
    c21->cd(6);
    gPad->SetLogy(0);
    h_VM_Epz_corrected->Draw("same");

    c21->cd(1);  
    gPad->SetTopMargin(0.08);  
    TLatex* title21_1 = new TLatex();
    title21_1->SetNDC(); 
    title21_1->SetTextSize(0.05);
    title21_1->SetTextAlign(22);  
    title21_1->DrawLatex(0.5, 0.97, "Bin Migration");  

    c21->cd(2);  
    gPad->SetTopMargin(0.08);  
    TLatex* title21_2 = new TLatex();
    title21_2->SetNDC(); 
    title21_2->SetTextSize(0.05);
    title21_2->SetTextAlign(22);  
    title21_2->DrawLatex(0.5, 0.97, "Normalized Migration (Probability)");  

    c21->cd(3);  
    gPad->SetTopMargin(0.08);  
    TLatex* title21_3 = new TLatex();
    title21_3->SetNDC(); 
    title21_3->SetTextSize(0.05);
    title21_3->SetTextAlign(22);  
    title21_3->DrawLatex(0.5, 0.97, "Purity");  

    c21->cd(4);  
    gPad->SetTopMargin(0.08);  
    TLatex* title21_4 = new TLatex();
    title21_4->SetNDC(); 
    title21_4->SetTextSize(0.05);
    title21_4->SetTextAlign(22);  
    title21_4->DrawLatex(0.5, 0.97, "Stability");

    c21->cd(5);  
    gPad->SetTopMargin(0.08);  
    TLatex* title21_5 = new TLatex();
    title21_5->SetNDC(); 
    title21_5->SetTextSize(0.05);
    title21_5->SetTextAlign(22);  
    title21_5->DrawLatex(0.5, 0.97, "Acceptance");

    c21->cd(6);  
    gPad->SetTopMargin(0.08);  
    TLatex* title21_6 = new TLatex();
    title21_6->SetNDC(); 
    title21_6->SetTextSize(0.05);
    title21_6->SetTextAlign(22);  
    title21_6->DrawLatex(0.5, 0.97, "Bin by Bin Corrected");
    
    c21->Print("./figures/plot_VM_Epz_migrationCorrectedPurityAcceptance.pdf");

// for x
    TCanvas* c22 = new TCanvas("c22","c22",1,1,1400,800);
    c22->Divide(3,2,0.01,0.01);
    c22->cd(1);
    h_x_migration->Draw("colzsame");   

    TH2D* h_x_mig_norm = (TH2D*)h_x_migration->Clone("h_x_mig_norm");
    int nx_x = h_x_mig_norm->GetNbinsX();
    int ny_x = h_x_mig_norm->GetNbinsY();
    for (int i=1; i<=nx_x; ++i) 
    {
        double rowSum = 0;
        for (int j=1; j<=ny_x; ++j) rowSum += h_x_mig_norm->GetBinContent(i,j);
        if (rowSum > 0) 
        {
            for (int j=1; j<=ny_x; ++j) 
            {
                double v = h_x_mig_norm->GetBinContent(i,j);
                h_x_mig_norm->SetBinContent(i,j, v / rowSum);
            }
        }
    }
    c22->cd(2);
    h_x_mig_norm->Draw("colzsame");

    for (int j=1; j<=ny_x; ++j) 
    {
        double recoColSum = 0;
        for (int i=1; i<=nx_x; ++i) recoColSum += h_x_migration->GetBinContent(i,j);
        double diag = h_x_migration->GetBinContent(j,j);
        if (recoColSum>0) h_x_purity->SetBinContent(j, diag / recoColSum);
    }
    c22->cd(3);
    h_x_purity->GetXaxis()->SetTitle("x_{reco bin}");
    h_x_purity->GetYaxis()->SetTitle("purity");
    h_x_purity->Draw("same");

    TH1D* h_x_stability = new TH1D("h_x_stability","Stability;x_{true bin};stability",nx_x,0,0.5);
    for (int i=1; i<=nx_x; ++i) 
    {
        double trueRowSum = 0;
        for (int j=1; j<=ny_x; ++j) trueRowSum += h_x_migration->GetBinContent(i,j);
        double diag = h_x_migration->GetBinContent(i,i);
        if (trueRowSum>0) h_x_stability->SetBinContent(i, diag / trueRowSum);
    }
    c22->cd(4);
    h_x_stability->Draw("same");
    
    int nb_x = h_x_beforeCut->GetNbinsX();
    for (int i=1; i<=nb_x; ++i) 
    {
        double Ntrue = h_x_beforeCut->GetBinContent(i);   
        double Naccepted = 0;
        for (int j=1; j<=nb_x; ++j) Naccepted += h_x_migration->GetBinContent(i,j);
        double acceptance = (Ntrue>0) ? Naccepted / Ntrue : 0;
        h_x_acceptance->SetBinContent(i, acceptance);
        double Nreco_in_bin_i = h_x_afterCut->GetBinContent(i);
        double corrected = (acceptance>0) ? Nreco_in_bin_i / acceptance : 0;
        h_x_corrected->SetBinContent(i, corrected);
    }
    c22->cd(5);
    h_x_acceptance->Draw("same");
    
    c22->cd(6);
    gPad->SetLogy(0);
    h_x_corrected->Draw("same");

    c22->cd(1);  
    gPad->SetTopMargin(0.08);  
    TLatex* title22_1 = new TLatex();
    title22_1->SetNDC(); 
    title22_1->SetTextSize(0.05);
    title22_1->SetTextAlign(22);  
    title22_1->DrawLatex(0.5, 0.97, "Bin Migration");  

    c22->cd(2);  
    gPad->SetTopMargin(0.08);  
    TLatex* title22_2 = new TLatex();
    title22_2->SetNDC(); 
    title22_2->SetTextSize(0.05);
    title22_2->SetTextAlign(22);  
    title22_2->DrawLatex(0.5, 0.97, "Normalized Migration (Probability)");  

    c22->cd(3);  
    gPad->SetTopMargin(0.08);  
    TLatex* title22_3 = new TLatex();
    title22_3->SetNDC(); 
    title22_3->SetTextSize(0.05);
    title22_3->SetTextAlign(22);  
    title22_3->DrawLatex(0.5, 0.97, "Purity");  

    c22->cd(4);  
    gPad->SetTopMargin(0.08);  
    TLatex* title22_4 = new TLatex();
    title22_4->SetNDC(); 
    title22_4->SetTextSize(0.05);
    title22_4->SetTextAlign(22);  
    title22_4->DrawLatex(0.5, 0.97, "Stability");

    c22->cd(5);  
    gPad->SetTopMargin(0.08);  
    TLatex* title22_5 = new TLatex();
    title22_5->SetNDC(); 
    title22_5->SetTextSize(0.05);
    title22_5->SetTextAlign(22);  
    title22_5->DrawLatex(0.5, 0.97, "Acceptance");

    c22->cd(6);  
    gPad->SetTopMargin(0.08);  
    TLatex* title22_6 = new TLatex();
    title22_6->SetNDC(); 
    title22_6->SetTextSize(0.05);
    title22_6->SetTextAlign(22);  
    title22_6->DrawLatex(0.5, 0.97, "Bin by Bin Corrected");

    c22->Print("./figures/plot_x_migrationCorrectedPurityAcceptance.pdf");
    
// for Q2
    TCanvas* c23 = new TCanvas("c23","c23",1,1,1400,800);
    c23->Divide(3,2,0.01,0.01);
    c23->cd(1);
    h_Q2_migration->Draw("colzsame");
    
    TH2D* h_Q2_mig_norm = (TH2D*)h_Q2_migration->Clone("h_Q2_mig_norm");
    int nx_Q2 = h_Q2_mig_norm->GetNbinsX();
    int ny_Q2 = h_Q2_mig_norm->GetNbinsY();
    for (int i=1; i<=nx_Q2; ++i) 
    {
        double rowSum = 0;
        for (int j=1; j<=ny_Q2; ++j) rowSum += h_Q2_mig_norm->GetBinContent(i,j);
        if (rowSum > 0) 
        {
            for (int j=1; j<=ny_Q2; ++j) 
            {
                double v = h_Q2_mig_norm->GetBinContent(i,j);
                h_Q2_mig_norm->SetBinContent(i,j, v / rowSum);
            }
        }
    }
    c23->cd(2);
    h_Q2_mig_norm->Draw("colzsame");

    for (int j=1; j<=ny_Q2; ++j) 
    {
        double recoColSum = 0;
        for (int i=1; i<=nx_Q2; ++i) recoColSum += h_Q2_migration->GetBinContent(i,j);
        double diag = h_Q2_migration->GetBinContent(j,j);
        if (recoColSum>0) h_Q2_purity->SetBinContent(j, diag / recoColSum);
    }
    c23->cd(3);
    h_Q2_purity->GetXaxis()->SetTitle("Q^{2}_{reco bin}");
    h_Q2_purity->GetYaxis()->SetTitle("purity");
    h_Q2_purity->Draw("same");

    TH1D* h_Q2_stability = new TH1D("h_Q2_stability","Stability;Q^{2}_{true bin};stability",nx_Q2,0,10);
    for (int i=1; i<=nx_Q2; ++i) 
    {
        double trueRowSum = 0;
        for (int j=1; j<=ny_Q2; ++j) trueRowSum += h_Q2_migration->GetBinContent(i,j);
        double diag = h_Q2_migration->GetBinContent(i,i);
        if (trueRowSum>0) h_Q2_stability->SetBinContent(i, diag / trueRowSum);
    }
    c23->cd(4);
    h_Q2_stability->Draw("same");
    
    int nb_Q2 = h_Q2_beforeCut->GetNbinsX();
    for (int i=1; i<=nb_Q2; ++i) 
    {
        double Ntrue = h_Q2_beforeCut->GetBinContent(i);   
        double Naccepted = 0;
        for (int j=1; j<=nb_Q2; ++j) Naccepted += h_Q2_migration->GetBinContent(i,j);
        double acceptance = (Ntrue>0) ? Naccepted / Ntrue : 0;
        h_Q2_acceptance->SetBinContent(i, acceptance);
        double Nreco_in_bin_i = h_Q2_afterCut->GetBinContent(i);
        double corrected = (acceptance>0) ? Nreco_in_bin_i / acceptance : 0;
        h_Q2_corrected->SetBinContent(i, corrected);
    }
    c23->cd(5);
    h_Q2_acceptance->Draw("same");
    
    c23->cd(6);
    gPad->SetLogy(0);
    h_Q2_corrected->Draw("same");

    c23->cd(1);  
    gPad->SetTopMargin(0.08);  
    TLatex* title23_1 = new TLatex();
    title23_1->SetNDC(); 
    title23_1->SetTextSize(0.05);
    title23_1->SetTextAlign(22);  
    title23_1->DrawLatex(0.5, 0.97, "Bin Migration");  

    c23->cd(2);  
    gPad->SetTopMargin(0.08);  
    TLatex* title23_2 = new TLatex();
    title23_2->SetNDC(); 
    title23_2->SetTextSize(0.05);
    title23_2->SetTextAlign(22);  
    title23_2->DrawLatex(0.5, 0.97, "Normalized Migration (Probability)");  

    c23->cd(3);  
    gPad->SetTopMargin(0.08);  
    TLatex* title23_3 = new TLatex();
    title23_3->SetNDC(); 
    title23_3->SetTextSize(0.05);
    title23_3->SetTextAlign(22);  
    title23_3->DrawLatex(0.5, 0.97, "Purity");  

    c23->cd(4);  
    gPad->SetTopMargin(0.08);  
    TLatex* title23_4 = new TLatex();
    title23_4->SetNDC(); 
    title23_4->SetTextSize(0.05);
    title23_4->SetTextAlign(22);  
    title23_4->DrawLatex(0.5, 0.97, "Stability");

    c23->cd(5);  
    gPad->SetTopMargin(0.08);  
    TLatex* title23_5 = new TLatex();
    title23_5->SetNDC(); 
    title23_5->SetTextSize(0.05);
    title23_5->SetTextAlign(22);  
    title23_5->DrawLatex(0.5, 0.97, "Acceptance"); 

    c23->cd(6);  
    gPad->SetTopMargin(0.08);  
    TLatex* title23_6 = new TLatex();
    title23_6->SetNDC(); 
    title23_6->SetTextSize(0.05);
    title23_6->SetTextAlign(22);  
    title23_6->DrawLatex(0.5, 0.97, "Bin by Bin Corrected"); 
    
    c23->Print("./figures/plot_Q2_migrationCorrectedPurityAcceptance.pdf");
    
// for y
    TCanvas* c24 = new TCanvas("c24","c24",1,1,1400,800);
    c24->Divide(3,2,0.01,0.01);
    c24->cd(1);
    h_y_migration->Draw("colzsame");

    TH2D* h_y_mig_norm = (TH2D*)h_y_migration->Clone("h_y_mig_norm");
    int nx_y = h_y_mig_norm->GetNbinsX();
    int ny_y = h_y_mig_norm->GetNbinsY();
    for (int i=1; i<=nx_y; ++i) 
    {
        double rowSum = 0;
        for (int j=1; j<=ny_y; ++j) rowSum += h_y_mig_norm->GetBinContent(i,j);
        if (rowSum > 0) 
        {
            for (int j=1; j<=ny_y; ++j) 
            {
                double v = h_y_mig_norm->GetBinContent(i,j);
                h_y_mig_norm->SetBinContent(i,j, v / rowSum);
            }
        }
    }
    c24->cd(2);
    h_y_mig_norm->Draw("colzsame");

    for (int j=1; j<=ny_y; ++j) 
    {
        double recoColSum = 0;
        for (int i=1; i<=nx_y; ++i) recoColSum += h_y_migration->GetBinContent(i,j);
        double diag = h_y_migration->GetBinContent(j,j);
        if (recoColSum>0) h_y_purity->SetBinContent(j, diag / recoColSum);
    }
    c24->cd(3);
    h_y_purity->GetXaxis()->SetTitle("y_{reco bin}");
    h_y_purity->GetYaxis()->SetTitle("purity");
    h_y_purity->Draw("same");

    TH1D* h_y_stability = new TH1D("h_y_stability","Stability;y_{true bin};stability",nx_y,0.1,0.85);
    for (int i=1; i<=nx_y; ++i) 
    {
        double trueRowSum = 0;
        for (int j=1; j<=ny_y; ++j) trueRowSum += h_y_migration->GetBinContent(i,j);
        double diag = h_y_migration->GetBinContent(i,i);
        if (trueRowSum>0) h_y_stability->SetBinContent(i, diag / trueRowSum);
    }
    c24->cd(4);
    h_y_stability->Draw("same");
    
    int nb_y = h_y_beforeCut->GetNbinsX();
    for (int i=1; i<=nb_y; ++i) 
    {
        double Ntrue = h_y_beforeCut->GetBinContent(i);   
        double Naccepted = 0;
        for (int j=1; j<=nb_y; ++j) Naccepted += h_y_migration->GetBinContent(i,j);
        double acceptance = (Ntrue>0) ? Naccepted / Ntrue : 0;
        h_y_acceptance->SetBinContent(i, acceptance);
        double Nreco_in_bin_i = h_y_afterCut->GetBinContent(i);
        double corrected = (acceptance>0) ? Nreco_in_bin_i / acceptance : 0;
        h_y_corrected->SetBinContent(i, corrected);
    }
    c24->cd(5);
    h_y_acceptance->Draw("same");
    
    c24->cd(6);
    gPad->SetLogy(0);
    h_y_corrected->Draw("same");

    c24->cd(1);  
    gPad->SetTopMargin(0.08);  
    TLatex* title24_1 = new TLatex();
    title24_1->SetNDC(); 
    title24_1->SetTextSize(0.05);
    title24_1->SetTextAlign(22);  
    title24_1->DrawLatex(0.5, 0.97, "Bin Migration");  

    c24->cd(2);  
    gPad->SetTopMargin(0.08);  
    TLatex* title24_2 = new TLatex();
    title24_2->SetNDC(); 
    title24_2->SetTextSize(0.05);
    title24_2->SetTextAlign(22);  
    title24_2->DrawLatex(0.5, 0.97, "Normalized Migration (Probability)");  

    c24->cd(3);  
    gPad->SetTopMargin(0.08);  
    TLatex* title24_3 = new TLatex();
    title24_3->SetNDC(); 
    title24_3->SetTextSize(0.05);
    title24_3->SetTextAlign(22);  
    title24_3->DrawLatex(0.5, 0.97, "Purity");  

    c24->cd(4);  
    gPad->SetTopMargin(0.08);  
    TLatex* title24_4 = new TLatex();
    title24_4->SetNDC(); 
    title24_4->SetTextSize(0.05);
    title24_4->SetTextAlign(22);  
    title24_4->DrawLatex(0.5, 0.97, "Stability");

    c24->cd(5);  
    gPad->SetTopMargin(0.08);  
    TLatex* title24_5 = new TLatex();
    title24_5->SetNDC(); 
    title24_5->SetTextSize(0.05);
    title24_5->SetTextAlign(22);  
    title24_5->DrawLatex(0.5, 0.97, "Acceptance");

    c24->cd(6);  
    gPad->SetTopMargin(0.08);  
    TLatex* title24_6 = new TLatex();
    title24_6->SetNDC(); 
    title24_6->SetTextSize(0.05);
    title24_6->SetTextAlign(22);  
    title24_6->DrawLatex(0.5, 0.97, "Bin by Bin Corrected");
    
    c24->Print("./figures/plot_y_migrationCorrectedPurityAcceptance.pdf");

// for E-pz
    TCanvas* c25 = new TCanvas("c25","c25",1,1,1400,800);
    c25->Divide(3,2,0.01,0.01);
    c25->cd(1);
    h_Epz_migration->GetXaxis()->SetRangeUser(12,27);
    h_Epz_migration->GetYaxis()->SetRangeUser(12,27);
    h_Epz_migration->Draw("colzsame");

    TH2D* h_Epz_mig_norm = (TH2D*)h_Epz_migration->Clone("h_Epz_mig_norm");
    int nx_Epz = h_Epz_mig_norm->GetNbinsX();
    int ny_Epz = h_Epz_mig_norm->GetNbinsY();
    for (int i=1; i<=nx_Epz; ++i) 
    {
        double rowSum = 0;
        for (int j=1; j<=ny_Epz; ++j) rowSum += h_Epz_mig_norm->GetBinContent(i,j);
        if (rowSum > 0) 
        {
            for (int j=1; j<=ny_Epz; ++j) 
            {
                double v = h_Epz_mig_norm->GetBinContent(i,j);
                h_Epz_mig_norm->SetBinContent(i,j, v / rowSum);
            }
        }
    }
    c25->cd(2);
    h_Epz_mig_norm->Draw("colzsame");

    for (int j=1; j<=ny_Epz; ++j) 
    {
        double recoColSum = 0;
        for (int i=1; i<=nx_Epz; ++i) recoColSum += h_Epz_migration->GetBinContent(i,j);
        double diag = h_Epz_migration->GetBinContent(j,j);
        if (recoColSum>0) h_Epz_purity->SetBinContent(j, diag / recoColSum);
    }
    c25->cd(3);
    h_Epz_purity->GetXaxis()->SetTitle("(E-p_{z})_{reco bin}");
    h_Epz_purity->GetYaxis()->SetTitle("purity");
    h_Epz_purity->Draw("same");

    TH1D* h_Epz_stability = new TH1D("h_Epz_stability","Stability;(E-p_{z})_{true bin};stability",nx_Epz,0,45);
    for (int i=1; i<=nx_Epz; ++i) 
    {
        double trueRowSum = 0;
        for (int j=1; j<=ny_Epz; ++j) trueRowSum += h_Epz_migration->GetBinContent(i,j);
        double diag = h_Epz_migration->GetBinContent(i,i);
        if (trueRowSum>0) h_Epz_stability->SetBinContent(i, diag / trueRowSum);
    }
    c25->cd(4);
    h_Epz_stability->Draw("same");
    
    int nb_Epz = h_Epz_MC->GetNbinsX();
    for (int i=1; i<=nb_Epz; ++i) 
    {
        double Ntrue = h_Epz_MC->GetBinContent(i);   
        double Naccepted = 0;
        for (int j=1; j<=nb_Epz; ++j) Naccepted += h_Epz_migration->GetBinContent(i,j);
        double acceptance = (Ntrue>0) ? Naccepted / Ntrue : 0;
        h_Epz_acceptance->SetBinContent(i, acceptance);
        double Nreco_in_bin_i = h_Epz_afterCut->GetBinContent(i);
        double corrected = (acceptance>0) ? Nreco_in_bin_i / acceptance : 0;
        h_Epz_corrected->SetBinContent(i, corrected);
    }
    c25->cd(5);
    h_Epz_acceptance->Draw("same");
    
    c25->cd(6);
    gPad->SetLogy(0);
    h_Epz_corrected->Draw("same");    

    c25->cd(1);  
    gPad->SetTopMargin(0.08);  
    TLatex* title25_1 = new TLatex();
    title25_1->SetNDC(); 
    title25_1->SetTextSize(0.05);
    title25_1->SetTextAlign(22);  
    title25_1->DrawLatex(0.5, 0.97, "Bin Migration");  

    c25->cd(2);  
    gPad->SetTopMargin(0.08);  
    TLatex* title25_2 = new TLatex();
    title25_2->SetNDC(); 
    title25_2->SetTextSize(0.05);
    title25_2->SetTextAlign(22);  
    title25_2->DrawLatex(0.5, 0.97, "Normalized Migration (Probability)");  

    c25->cd(3);  
    gPad->SetTopMargin(0.08);  
    TLatex* title25_3 = new TLatex();
    title25_3->SetNDC(); 
    title25_3->SetTextSize(0.05);
    title25_3->SetTextAlign(22);  
    title25_3->DrawLatex(0.5, 0.97, "Purity");  

    c25->cd(4);  
    gPad->SetTopMargin(0.08);  
    TLatex* title25_4 = new TLatex();
    title25_4->SetNDC(); 
    title25_4->SetTextSize(0.05);
    title25_4->SetTextAlign(22);  
    title25_4->DrawLatex(0.5, 0.97, "Stability"); 

    c25->cd(5);  
    gPad->SetTopMargin(0.08);  
    TLatex* title25_5 = new TLatex();
    title25_5->SetNDC(); 
    title25_5->SetTextSize(0.05);
    title25_5->SetTextAlign(22);  
    title25_5->DrawLatex(0.5, 0.97, "Acceptance"); 

    c25->cd(6);  
    gPad->SetTopMargin(0.08);  
    TLatex* title25_6 = new TLatex();
    title25_6->SetNDC(); 
    title25_6->SetTextSize(0.05);
    title25_6->SetTextAlign(22);  
    title25_6->DrawLatex(0.5, 0.97, "Bin by Bin Corrected"); 
    
    c25->Print("./figures/plot_Epz_migrationCorrectedPurityAcceptance.pdf");

// for E/p
    TCanvas* c26 = new TCanvas("c26","c26",1,1,1400,800);
    c26->Divide(3,2,0.01,0.01);
    c26->cd(1);
    h_EoverP_migration->GetXaxis()->SetRangeUser(0.95,1.05);
    h_EoverP_migration->GetYaxis()->SetRangeUser(0.95,1.05);
    h_EoverP_migration->Draw("colzsame");

    TH2D* h_EoverP_mig_norm = (TH2D*)h_EoverP_migration->Clone("h_EoverP_mig_norm");
    int nx_EoverP = h_EoverP_mig_norm->GetNbinsX();
    int ny_EoverP = h_EoverP_mig_norm->GetNbinsY();
    for (int i=1; i<=nx_EoverP; ++i) 
    {
        double rowSum = 0;
        for (int j=1; j<=ny_EoverP; ++j) rowSum += h_EoverP_mig_norm->GetBinContent(i,j);
        if (rowSum > 0) 
        {
            for (int j=1; j<=ny_EoverP; ++j) 
            {
                double v = h_EoverP_mig_norm->GetBinContent(i,j);
                h_EoverP_mig_norm->SetBinContent(i,j, v / rowSum);
            }
        }
    }
    c26->cd(2);
    h_EoverP_mig_norm->Draw("colzsame");

    for (int j=1; j<=ny_EoverP; ++j) 
    {
        double recoColSum = 0;
        for (int i=1; i<=nx_EoverP; ++i) recoColSum += h_EoverP_migration->GetBinContent(i,j);
        double diag = h_EoverP_migration->GetBinContent(j,j);
        if (recoColSum>0) h_EoverP_purity->SetBinContent(j, diag / recoColSum);
    }
    c26->cd(3);
    h_EoverP_purity->GetXaxis()->SetTitle("(E/|p|)_{reco bin}");
    h_EoverP_purity->GetYaxis()->SetTitle("purity");
    h_EoverP_purity->Draw("same");

    TH1D* h_EoverP_stability = new TH1D("h_EoverP_stability","Stability;(E/|p|)_{true bin};stability",nx_EoverP,0,2);
    for (int i=1; i<=nx_EoverP; ++i) 
    {
        double trueRowSum = 0;
        for (int j=1; j<=ny_EoverP; ++j) trueRowSum += h_EoverP_migration->GetBinContent(i,j);
        double diag = h_EoverP_migration->GetBinContent(i,i);
        if (trueRowSum>0) h_EoverP_stability->SetBinContent(i, diag / trueRowSum);
    }
    c26->cd(4);
    h_EoverP_stability->Draw("same");
    
    int nb_EoverP = h_EoverP_MC->GetNbinsX();
    for (int i=1; i<=nb_EoverP; ++i) 
    {
        double Ntrue = h_EoverP_MC->GetBinContent(i);   
        double Naccepted = 0;
        for (int j=1; j<=nb_EoverP; ++j) Naccepted += h_EoverP_migration->GetBinContent(i,j);
        double acceptance = (Ntrue>0) ? Naccepted / Ntrue : 0;
        h_EoverP_acceptance->SetBinContent(i, acceptance);
        double Nreco_in_bin_i = h_EoverP_afterCut->GetBinContent(i);
        double corrected = (acceptance>0) ? Nreco_in_bin_i / acceptance : 0;
        h_EoverP_corrected->SetBinContent(i, corrected);
    }
    c26->cd(5);
    h_EoverP_acceptance->Draw("same");
    
    c26->cd(6);
    gPad->SetLogy(0);
    h_EoverP_corrected->Draw("same");    

    c26->cd(1);  
    gPad->SetTopMargin(0.08);  
    TLatex* title26_1 = new TLatex();
    title26_1->SetNDC(); 
    title26_1->SetTextSize(0.05);
    title26_1->SetTextAlign(22);  
    title26_1->DrawLatex(0.5, 0.97, "Bin Migration");  

    c26->cd(2);  
    gPad->SetTopMargin(0.08);  
    TLatex* title26_2 = new TLatex();
    title26_2->SetNDC(); 
    title26_2->SetTextSize(0.05);
    title26_2->SetTextAlign(22);  
    title26_2->DrawLatex(0.5, 0.97, "Normalized Migration (Probability)");  

    c26->cd(3);  
    gPad->SetTopMargin(0.08);  
    TLatex* title26_3 = new TLatex();
    title26_3->SetNDC(); 
    title26_3->SetTextSize(0.05);
    title26_3->SetTextAlign(22);  
    title26_3->DrawLatex(0.5, 0.97, "Purity");  

    c26->cd(4);  
    gPad->SetTopMargin(0.08);  
    TLatex* title26_4 = new TLatex();
    title26_4->SetNDC(); 
    title26_4->SetTextSize(0.05);
    title26_4->SetTextAlign(22);  
    title26_4->DrawLatex(0.5, 0.97, "Stability");

    c26->cd(5);  
    gPad->SetTopMargin(0.08);  
    TLatex* title26_5 = new TLatex();
    title26_5->SetNDC(); 
    title26_5->SetTextSize(0.05);
    title26_5->SetTextAlign(22);  
    title26_5->DrawLatex(0.5, 0.97, "Acceptance"); 

    c26->cd(6);  
    gPad->SetTopMargin(0.08);  
    TLatex* title26_6 = new TLatex();
    title26_6->SetNDC(); 
    title26_6->SetTextSize(0.05);
    title26_6->SetTextAlign(22);  
    title26_6->DrawLatex(0.5, 0.97, "Bin by Bin Corrected");
    
    c26->Print("./figures/plot_EoverP_migrationCorrectedPurityAcceptance.pdf");
}

