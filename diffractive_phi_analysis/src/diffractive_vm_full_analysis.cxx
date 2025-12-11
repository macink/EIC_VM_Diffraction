#include "pleaseIncludeMe.h"
#include <TVector3.h>
#include <TVector.h>
#include <cmath>
#include "TF1.h"
#include "TF2.h"
#include "TSystem.h"

/*-----------------------------------------------------
    Analysis for Exclusive Diffractive VM Production
- Purpose: reconstuct kinematic variables and the t-distribution
    using methods E, L, and projection with wedge cut
- Input: root file/s
- Output: root file
------------------------------------------------------*/

using namespace std;


/*-----------------------------
    For incoherent analysis
------------------------------*/
struct Cluster_EEMC
{
    float x, y, energy;
};

struct Cluster_ZDC 
{
    float x, y, z, energy;
};

struct Hit_RP 
{
    float x, y, z, energy;
};

struct Hit_OMD 
{
    float x, y, z, energy;
};

struct Event 
{
    vector<Cluster_EEMC> clusters_eemc;
    vector<Cluster_ZDC> clusters_zdc;
    vector<Hit_RP> hit_rp;
    vector<Hit_OMD> hit_omd;
};

auto giveme_t_method_L(TLorentzVector eIn, TLorentzVector eOut, TLorentzVector pIn, TLorentzVector vmOut)
{
    /*---------------------------------------------------------
        Function to calculate method L for t reconstruction
    Input: (four momenta) P_e, P_e', P_A, P_VM
    Returns: double, |t| = -(P_A'corr - P_A)^2 
    ----------------------------------------------------------*/
	TLorentzVector aInVec(pIn.Px()*197,pIn.Py()*197,pIn.Pz()*197,sqrt(pIn.Px()*197*pIn.Px()*197 + pIn.Py()*197*pIn.Py()*197 + pIn.Pz()*197*pIn.Pz()*197 + MASS_AU197*MASS_AU197) );
	double method_L = 0;
	TLorentzVector a_beam_scattered = aInVec-(vmOut+eOut-eIn);
	double p_Aplus = a_beam_scattered.E()+a_beam_scattered.Pz();
	double p_TAsquared = TMath::Power(a_beam_scattered.Pt(),2);
	double p_Aminus = (MASS_AU197*MASS_AU197 + p_TAsquared) / p_Aplus;
	TLorentzVector a_beam_scattered_corr; 
	a_beam_scattered_corr.SetPxPyPzE(a_beam_scattered.Px(),a_beam_scattered.Py(),(p_Aplus-p_Aminus)/2., (p_Aplus+p_Aminus)/2. );
	method_L = -(a_beam_scattered_corr-aInVec).Mag2();
	return method_L; 
}

auto giveme_t_new_method(TLorentzVector eIn, TLorentzVector eOut, TLorentzVector pIn, TLorentzVector vmOut)
{
    /*---------------------------------------------------------
        Function to calculate method L for t reconstruction
    Input: four momenta, P_e, P_e', P_A, P_VM
    Returns: four vector, P = sqrt(|t|) = -(P_A'corr - P_A)
    ----------------------------------------------------------*/
	TLorentzVector aInVec(pIn.Px()*197,pIn.Py()*197,pIn.Pz()*197,sqrt(pIn.Px()*197*pIn.Px()*197 + pIn.Py()*197*pIn.Py()*197 + pIn.Pz()*197*pIn.Pz()*197 + MASS_AU197*MASS_AU197) );
	double method_L = 0;
	TLorentzVector a_beam_scattered = aInVec-(vmOut+eOut-eIn);
	double p_Aplus = a_beam_scattered.E()+a_beam_scattered.Pz();
	double p_TAsquared = TMath::Power(a_beam_scattered.Pt(),2);
	double p_Aminus = (MASS_AU197*MASS_AU197 + p_TAsquared) / p_Aplus;
	TLorentzVector a_beam_scattered_corr; 
	a_beam_scattered_corr.SetPxPyPzE(a_beam_scattered.Px(),a_beam_scattered.Py(),(p_Aplus-p_Aminus)/2., (p_Aplus+p_Aminus)/2. );
	method_L = -(a_beam_scattered_corr-aInVec).Mag2();
	TLorentzVector method_L_4vect = -(a_beam_scattered_corr-aInVec); 
	return method_L_4vect;
}

int diffractive_vm_full_analysis(TString rec_file, TString outputfile)
{	
    // reads in lists of root files	
    TString listname = rec_file;
    TString outname = outputfile;
    TChain *chain = new TChain("events");
    int nfiles = 0;
    char filename[512];
    ifstream *inputstream = new ifstream;
    inputstream->open(listname.Data());
    if(!inputstream)
    {
      printf("[e] Cannot open file list: %s\n", listname.Data());
    }
    while(inputstream->good())
    {
        inputstream->getline(filename, 512);
        if(inputstream->good())
        {
            TFile *ftmp = TFile::Open(filename, "read");
            if(!ftmp||!(ftmp->IsOpen())||!(ftmp->GetNkeys())) 
            {
                printf("[e] Could you open file: %s\n", filename);
            } 
            else
            {
                cout<<"[i] Add "<<nfiles<<"th file: "<<filename<<endl;
                chain->Add(filename);
                nfiles++;
            }
        }
    }
    inputstream->close();
    printf("[i] Read in %d files with %lld events in total\n", nfiles, chain->GetEntries());
    TTreeReader tree_reader(chain);

    // MC particle arrays for each MC particle
    TTreeReaderArray<int> mc_genStatus_array = {tree_reader, "MCParticles.generatorStatus"};
    TTreeReaderArray<double> mc_px_array = {tree_reader, "MCParticles.momentum.x"};
    TTreeReaderArray<double> mc_py_array = {tree_reader, "MCParticles.momentum.y"};
    TTreeReaderArray<double> mc_pz_array = {tree_reader, "MCParticles.momentum.z"};
    TTreeReaderArray<double> mc_mass_array = {tree_reader, "MCParticles.mass"};
    TTreeReaderArray<int> mc_pdg_array = {tree_reader, "MCParticles.PDG"};

    // Reconstructed EcalEndcapNClusters arrays
    TTreeReaderArray<float> em_energy_array = {tree_reader, "EcalEndcapNClusters.energy"};
    TTreeReaderArray<float> em_x_array = {tree_reader, "EcalEndcapNClusters.position.x"};
    TTreeReaderArray<float> em_y_array = {tree_reader, "EcalEndcapNClusters.position.y"};
    TTreeReaderArray<float> emhits_x_array = {tree_reader, "EcalEndcapNRecHits.position.x"};
    TTreeReaderArray<float> emhits_y_array = {tree_reader, "EcalEndcapNRecHits.position.y"};
    TTreeReaderArray<float> emhits_energy_array = {tree_reader, "EcalEndcapNRecHits.energy"};
    TTreeReaderArray<unsigned int> em_rec_id_array = {tree_reader, "EcalEndcapNClusterAssociations.recID"};
    TTreeReaderArray<unsigned int> em_sim_id_array = {tree_reader, "EcalEndcapNClusterAssociations.simID"};

    // Reconstructed Calorimeter track arrays
    TTreeReaderArray<float> emCal_trk_x_array = {tree_reader, "_CalorimeterTrackProjections_points.position.x"};
    TTreeReaderArray<float> emCal_trk_y_array = {tree_reader, "_CalorimeterTrackProjections_points.position.y"};
    TTreeReaderArray<float> emCal_trk_px_array = {tree_reader, "_CalorimeterTrackProjections_points.momentum.x"};
    TTreeReaderArray<float> emCal_trk_py_array = {tree_reader, "_CalorimeterTrackProjections_points.momentum.y"};
    TTreeReaderArray<float> emCal_trk_pz_array = {tree_reader, "_CalorimeterTrackProjections_points.momentum.y"};

    // Reconstructed particles pz array for each reconstructed particle
    TTreeReaderArray<float> reco_px_array = {tree_reader, "ReconstructedChargedParticles.momentum.x"};
    TTreeReaderArray<float> reco_py_array = {tree_reader, "ReconstructedChargedParticles.momentum.y"};
    TTreeReaderArray<float> reco_pz_array = {tree_reader, "ReconstructedChargedParticles.momentum.z"};
    TTreeReaderArray<int>   reco_pdg_array = {tree_reader, "ReconstructedChargedParticles.PDG"};
    TTreeReaderArray<unsigned int> reco_clusterE_array = {tree_reader, "ReconstructedChargedParticles.clusters_begin"};
    TTreeReaderArray<float> reco_charge_array = {tree_reader, "ReconstructedChargedParticles.charge"};
    TTreeReaderArray<unsigned int> rec_id = {tree_reader, "ReconstructedChargedParticleAssociations.recID"};
    TTreeReaderArray<unsigned int> sim_id = {tree_reader, "ReconstructedChargedParticleAssociations.simID"};

    // Reconstructed Far Forward ZDC Neutrals 
    TTreeReaderArray<float> zdc_x_array = {tree_reader, "ReconstructedFarForwardZDCNeutrals.referencePoint.x"};
    TTreeReaderArray<float> zdc_y_array = {tree_reader, "ReconstructedFarForwardZDCNeutrals.referencePoint.y"};
    TTreeReaderArray<float> zdc_z_array = {tree_reader, "ReconstructedFarForwardZDCNeutrals.referencePoint.z"};
    TTreeReaderArray<float> zdc_energy_array = {tree_reader, "ReconstructedFarForwardZDCNeutrals.energy"};

    // Reconstructed ECal ZDC Clusters
    TTreeReaderArray<float> Ecal_zdc_x_array = {tree_reader, "EcalFarForwardZDCClusters.position.x"};
    TTreeReaderArray<float> Ecal_zdc_y_array = {tree_reader, "EcalFarForwardZDCClusters.position.y"};
    TTreeReaderArray<float> Ecal_zdc_z_array = {tree_reader, "EcalFarForwardZDCClusters.position.z"};
    TTreeReaderArray<float> Ecal_zdc_energy_array = {tree_reader, "EcalFarForwardZDCClusters.energy"};

    // Reconstructed HCal ZDC Clusters
    TTreeReaderArray<float> Hcal_zdc_x_array = {tree_reader, "HcalFarForwardZDCClusters.position.x"};
    TTreeReaderArray<float> Hcal_zdc_y_array = {tree_reader, "HcalFarForwardZDCClusters.position.y"};
    TTreeReaderArray<float> Hcal_zdc_z_array = {tree_reader, "HcalFarForwardZDCClusters.position.z"};
    TTreeReaderArray<float> Hcal_zdc_energy_array = {tree_reader, "HcalFarForwardZDCClusters.energy"};

    // Forward Roman Pot Rec Hits
    TTreeReaderArray<float> rp_x_array = {tree_reader, "ForwardRomanPotRecHits.position.x"};
    TTreeReaderArray<float> rp_y_array = {tree_reader, "ForwardRomanPotRecHits.position.y"};
    TTreeReaderArray<float> rp_z_array = {tree_reader, "ForwardRomanPotRecHits.position.z"};

    // Forward Off Momentum Tracker Rec Hits
    TTreeReaderArray<float> omd_x_array = {tree_reader, "ForwardOffMTrackerRecHits.position.x"};
    TTreeReaderArray<float> omd_y_array = {tree_reader, "ForwardOffMTrackerRecHits.position.y"};
    TTreeReaderArray<float> omd_z_array = {tree_reader, "ForwardOffMTrackerRecHits.position.z"};

    // Include DIS Kinematics
    TTreeReaderArray<float> DIS_y = {tree_reader, "InclusiveKinematicsElectron.y"};
    TTreeReaderArray<float> DIS_x = {tree_reader, "InclusiveKinematicsElectron.x"};
    TTreeReaderArray<float> DIS_Q2 = {tree_reader, "InclusiveKinematicsElectron.Q2"};

    TString output_name_dir = outputfile+"_output.root";
    cout << "Output file = " << output_name_dir << endl;
    TFile* output = new TFile(output_name_dir,"RECREATE");

    TTree* outputTree = new TTree("miniTree", "Tree with structured event data");
    Event event;
    outputTree->Branch("event", &event);

    /*--------------------------------------------------------------------------------------
        - Below are all of the histograms for analysis, they are labelled accordingly
        - "_after" refers the histogram being generated after the event cut, etc
    ----------------------------------------------------------------------------------------*/

    // MC events
    TH1D* h_Q2_e = new TH1D("h_Q2_e",";Q^{2}_{MC} [GeV/c]^{2}",100,0,10);
	TH1D* h_y_e = new TH1D("h_y_e",";y_{MC}",100,0.01,0.85);
	TH1D* h_energy_MC = new TH1D("h_energy_MC",";E_{MC} [GeV]",100,0,20);
    TH1D* h_energy_MC_after = new TH1D("h_energy_MC_after",";E_{MC} [GeV]",100,0,20);
	TH1D* h_e_pt_MC = new TH1D("h_e_pt_MC",";p_{T,e,MC} [GeV/c]",200,0,20);
    TH1D* h_e_pt_MC_after = new TH1D("h_e_pt_MC_after",";p_{T,e,MC} [GeV/c]",200,0,20);
	TH1D* h_e_pz_MC = new TH1D("h_e_pz_MC",";p_{z,e,MC} [GeV/c]",200,0,20);
    TH1D* h_e_pz_MC_after = new TH1D("h_e_pz_MC_after",";p_{z,e,MC} [GeV/c]",200,0,20);
	TH1D* h_e_p_MC = new TH1D("h_e_p_MC",";p_{e,MC} [GeV/c]",200,0,20);
    TH1D* h_e_p_MC_after = new TH1D("h_e_p_MC_after",";p_{e,MC} [GeV/c]",200,0,20);
    TH1D* h_eta_MC = new TH1D("h_eta_MC",";#eta_{MC}",100,-4,4);
    TH1D* h_eta_MC_after = new TH1D("h_eta_MC_after",";#eta_{MC}",100,-4,4);
	TH1D* h_phi_MC = new TH1D("h_phi_MC",";#phi_{MC}",100,-3.14,3.14);
    TH1D* h_phi_MC_after = new TH1D("h_phi_MC_after",";#phi_{MC}",100,-3.14,3.14);
	TH1D* h_theta_MC = new TH1D("h_theta_MC",";#theta_{MC}",100,0,3.14);
    TH1D* h_theta_MC_after = new TH1D("h_theta_MC_after",";#theta_{MC}",100,0,3.14);
    TH1D* h_VM_mass_MC = new TH1D("h_VM_mass_MC",";VM_{MC} mass [GeV/c^{2}]",200,0,2);
    TH1D* h_VM_pt_MC = new TH1D("h_VM_pt_MC",";p_{T,VM,MC} [GeV/c]",200,0,10);
	TH1D* h_VM_pz_MC = new TH1D("h_VM_pz_MC",";p_{z,VM,MC} [GeV/c]",200,0,20);
	TH1D* h_VM_p_MC = new TH1D("h_VM_p_MC",";p_{VM,MC} [GeV/c]",200,0,20);
    TH1D* h_VM_Epz_MC = new TH1D("h_VM_Epz_MC",";(E_{VM,MC}-p_{z,VM,MC}) [GeV]",200,0,20);
    TH1D* h_Epz_MC = new TH1D("h_Epz_MC", ";(E_{MC} - p_{z,MC}) [GeV]",100,15,25);
    TH1D* h_Epz_MC_after = new TH1D("h_Epz_MC_after", ";(E_{MC} - p_{z,MC}) [GeV]",100,15,25);
    TH2D* h_EvsP_MC = new TH2D("h_EvsP_MC",";|p|_{MC} [GeV]; E_{MC} [GeV]",100,0,20,100,0,20);
    TH1D* h_EoverP_MC = new TH1D("h_EoverP_MC",";E_{MC}/|p|_{MC}",100,0,2);
    TH1D* h_EoverP_MC_after = new TH1D("h_EoverP_MC_after",";E_{MC}/|p|_{MC}",100,0,2);
    TH1D* h_t_MC = new TH1D("h_t_MC",";|t|_{MC} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_MC_before = new TH1D("h_t_MC_before",";|t|_{MC} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_sigma = new TH1D("h_sigma",":#sigma; counts",100,0,20);

    // RECO track
	TH1D* h_energy_REC_trk = new TH1D("h_energy_REC_trk",";E_{trk} [GeV]",100,0,20);
    TH2D* h_energy_res_trk = new TH2D("h_energy_res_trk",";E_{MC} [GeV]; (E_{MC}-E_{trk})/E_{MC}",100,0,20,1000,-1,1);
	TH1D* h_eta_REC_trk = new TH1D("h_eta_REC_trk",";#eta_{trk}",100,-4,4);
	TH1D* h_e_pt_REC_trk = new TH1D("h_e_pt_REC_trk",";p_{T,e,trk} [GeV/c]",200,0,20);
	TH1D* h_e_pz_REC_trk = new TH1D("h_e_pz_REC_trk",";p_{z,e,trk} [GeV/c]",200,0,20);
	TH1D* h_e_p_REC_trk = new TH1D("h_e_p_REC_trk",";p_{e,trk} [GeV/c]",200,0,20);
    TH1D* h_t_REC_trk_cut = new TH1D("h_t_REC_trk_cut",";|t|_{trk} [GeV/c]^{2}; counts",100,0,0.2);
    TH2D* h_t_res_trk = new TH2D("h_t_res_trk",";|t|_{MC} [GeV/c]^{2}; (|t|_{MC}-|t|_{trk})/|t|_{MC}",100,0,0.2,1000,-10,10);
    TH2D* h_t_res_trk_cut = new TH2D("h_t_res_trk_cut",";|t|_{MC} [GeV/c]^{2}; (|t|_{MC}-|t|_{trk})/|t|_{MC}",100,0,0.2,1000,-10,10);
    TH2D* h_t_response_trk = new TH2D("h_t_response_trk","; |t|_{trk} [GeV/c]^{2};|t|_{MC} [GeV/c]^{2}",100,0,0.2,100,0,0.2);
    TH2D* h_t_response_trk_cut = new TH2D("h_t_response_trk_cut","; |t|_{trk} [GeV/c]^{2};|t|_{MC} [GeV/c]^{2}",100,0,0.2,100,0,0.2);
    TH1D* h_trk_position_x_REC = new TH1D("h_trk_position_x_REC",";x_{trk} [mm]",80,-800,800);
    TH1D* h_trk_position_y_REC = new TH1D("h_trk_position_y_REC",";y_{trk} [mm]",80,-800,800);
    TH2D* h_XvsY_trk = new TH2D("h_XvsY_trk",";x [mm]; y [mm]",80,-800,800,80,-800,800);

	// RECO Q2
    TH1D* h_Q2REC_e_EEMC = new TH1D("h_Q2REC_e_EEMC",";Q^{2}_{EEMC} [GeV/c]^{2}",100,0,10);
    TH2D* h_Q2_res = new TH2D("h_Q2_res",";Q^{2}_{MC} [GeV/c]^{2}; (Q^{2}_{MC}-Q^{2}_{EEMC})/Q^{2}_{MC}",100,0,10,1000,-1,1);
    TH1D* h_Q2_res2 = new TH1D("h_Q2_res2",";(Q^{2}_{REC}-Q^{2}_{MC})/Q^{2}_{MC}",100,-0.4,0.4);
    TH2D* h_Q2_response = new TH2D("h_Q2_response","; Q^{2}_{EEMC} [GeV/c]^{2};Q^{2}_{MC} [GeV/c]^{2}",100,0,10,100,0,10);
	TH1D* h_dQ2overQ2_REC = new TH1D("h_dQ2overQ2_REC",";dQ^{2}_{REC}/Q^{2}_{MC}",100,-2,2);
    TH2D* h_Q2_migration = new TH2D("h_Q2_migration",";Q^{2}_{MC} [GeV/c]^{2};Q^{2}_{RECO} [GeV/c]^{2}",100,0,10,100,0,10);
    
    // RECO Bjorken-x
    TH2D* h_x_res_cut = new TH2D("h_x_res_cut",";x_{MC} ; x_{MC}-x_{RECO}/x_{MC}",100,0,0.5,100,0,0.5);
    TH1D* h_x_res_cut2 = new TH1D("h_x_res_cut2",";x_{REC}-x_{MC}/x_{MC}",100,-0.4,0.4);
    TH2D* h_x_response_cut = new TH2D("h_x_response_cut",";x_{RECO} ; x_{MC}",100,0,0.5,100,0,0.5);
    TH2D* h_x_migration = new TH2D("h_x_migration",";x_{MC};x_{RECO} ",100,0,0.5,100,0,0.5);
    TH1D* h_x_beforeCut = new TH1D("h_x_beforeCut",";x",100,0,0.5);
    TH1D* h_x = new TH1D("h_x",";x_{MC}",100,0,0.5);
    TH1D* h_x_REC = new TH1D("h_x_REC",";x_{RECO}",100,0,0.5);
    TH1D* h_x_afterCut = new TH1D("h_x_afterCut",";x_{RECO}",100,0,0.5);

    // RECO y
	TH1D* h_yREC_e_EEMC = new TH1D("h_yREC_e_EEMC",";y_{EEMC}",100,0.01,0.85);
    TH2D* h_y_res = new TH2D("h_y_res",";y_{e,MC} ;(y_{e,MC}-y_{e,EEMC})/y_{e,MC}",100,0.01,0.85,1000,-1,1);
    TH1D* h_y_res2 = new TH1D("h_y_res2",";(y_{e,REC}-y_{e,MC})/y_{e,MC}",100,-0.4,0.4);
    TH2D* h_y_response = new TH2D("h_y_response"," ; y_{EEMC};y_{MC}",100,0.01,0.85,100,0.01,0.85);
    TH1D* h_dyOvery_REC = new TH1D("h_dyOvery_REC",";dy/y",100,-2,2);
    TH2D* h_y_migration = new TH2D("h_y_migration",";y_{MC};y_{RECO}",100,0.01,0.85,100,0.01,0.85);

    // RECO energy
	TH1D* h_energy_REC_EEMC = new TH1D("h_energy_REC_EEMC",";E_{EEMC} [GeV]",100,0,20);
    TH2D* h_energy_res_EEMC = new TH2D("h_energy_res_EEMC",";E_{MC} [GeV]; (E_{MC}-E_{EEMC})/E_{MC}",100,0,20,1000,-1,1);
    TH2D* h_energy_res_EEMC_after = new TH2D("h_energy_res_EEMC_after",";E_{MC} [GeV]; (E_{MC}-E_{EEMC})/E_{MC}",100,0,20,1000,-1,1);
    TH2D* h_energy_response_EEMC = new TH2D("h_energy_response_EEMC","; E_{EEMC} [GeV];E_{MC} [GeV]",100,0,20,100,0,20);
    TH2D* h_energy_response_EEMC_after = new TH2D("h_energy_response_EEMC_after","; E_{EEMC} [GeV];E_{MC} [GeV]",100,0,20,100,0,20);
    TH2D* h_energy_migration = new TH2D("h_energy_migration",";E_{MC} [GeV];E_{RECO} [GeV]",100,0,20,100,0,20);
    TH1D* h_energy_REC_EEMC_after = new TH1D("h_energy_REC_EEMC_after",";E_{EEMC} [GeV]",100,0,20);

    // RECO eta
	TH1D* h_eta_REC_EEMC = new TH1D("h_eta_REC_EEMC",";#eta_{EEMC}",100,-4,4);
    TH1D* h_eta_REC_EEMC_after = new TH1D("h_eta_REC_EEMC_after",";#eta_{EEMC}",100,-4,4);
    TH1D* h_eta_diff = new TH1D("h_eta_diff",";#eta_{EEMC}-#eta_{MC}",100,-4,4);
    TH1D* h_eta_diff_after = new TH1D("h_eta_diff_after",";#eta_{EEMC}-#eta_{MC}",100,-4,4);

    // RECO phi
    TH1D* h_phi_REC_EEMC = new TH1D("h_phi_REC_EEMC",";#phi_{EEMC}",100,-3.14,3.14);
    TH1D* h_phi_REC_EEMC_after = new TH1D("h_phi_REC_EEMC_after",";#phi_{EEMC}",100,-3.14,3.14);
    TH1D* h_phi_diff = new TH1D("h_phi_diff",";#phi_{EEMC}-#phi_{MC}",100,-3.14,3.14);
    TH1D* h_phi_diff_cut = new TH1D("h_phi_diff_cut",";#phi_{EEMC}-#phi_{MC}",100,-3.14,3.14);
    TH1D* h_phi_diff_after = new TH1D("h_phi_diff_after",";#phi_{EEMC}-#phi_{MC}",100,-3.14,3.14);
    TH2D* h_phi_resolution = new TH2D("h_phi_resolution",";#phi_{MC}; (#phi_{RECO}-#phi_{MC})/#phi_{MC}",100,-3.14,3.14,1000,-1,1);
    TH1D* h_phi_res_diff = new TH1D("h_phi_res_diff",";#phi_{RECO}-#phi_{MC}",100,-3.14,3.14);

    // RECO theta
    TH1D* h_theta_REC_EEMC = new TH1D("h_theta_REC_EEMC",";#theta_{EEMC}",100,0,3.14);
    TH1D* h_theta_REC_EEMC_after = new TH1D("h_theta_REC_EEMC_after",";#theta_{EEMC}",100,-0,3.14);
    TH2D* h_theta_response_EEMC = new TH2D("h_theta_response_EEMC","; #theta_{EEMC}; #theta_{MC}",100,0,3.14,100,0,3.14);
    TH2D* h_theta_response_EEMC_after = new TH2D("h_theta_response_EEMC_after","; #theta_{EEMC}; #theta_{MC}",100,0,3.14,100,0,3.14);
    TH1D* h_theta_diff = new TH1D("h_theta_diff",";#theta_{EEMC}-#theta_{MC}",100,0,3.14);
    TH1D* h_theta_diff_after = new TH1D("h_theta_diff_after",";#theta_{EEMC}-#theta_{MC}",100,0,3.14);
    TH2D* h_theta_migration = new TH2D("h_theta_migration","; #theta_{MC};#theta_{RECO}",100,0,3.14,100,0,3.14);
    TH2D* h_theta_resolution = new TH2D("h_theta_resolution",";#theta_{MC}; (#theta_{RECO}-#theta_{MC})/#theta_{MC}",100,0,3.14,1000,-1,1);

    // RECO e' P
	TH1D* h_e_pt_REC_EEMC = new TH1D("h_e_pt_REC_EEMC",";p_{T,e,EEMC} [GeV/c]",200,0,20);
    TH2D* h_e_pt_res = new TH2D("h_e_pt_res",";p_{T,e,MC} [GeV/c]; (p_{T,e,EEMC}-p_{T,e,MC})/p_{T,e,MC}",200,0,20,2000,-1,1);
	TH1D* h_e_pz_REC_EEMC = new TH1D("h_e_pz_REC_EEMC",";p_{z,e,EEMC} [GeV/c]",200,0,20);
    TH2D* h_e_pz_res = new TH2D("h_e_pz_res",";p_{z,e,MC} [GeV/c]; (p_z,e,EEMC}-p_{z,e,MC})/p_{z,e,MC}",200,0,20,2000,-1,1);
    TH2D* h_e_pz_response = new TH2D("h_e_pz_response","; p_{z,e,EEMC} [GeV/c];p_{z,e,MC} [GeV/c]",200,0,20,200,0,2);
	TH1D* h_e_p_REC_EEMC = new TH1D("h_e_p_REC_EEMC",";p_{e,EEMC} [GeV/c]",200,0,20);
	TH2D* h_e_p_res = new TH2D("h_e_p_res",";p_{e,MC} [GeV/c]; (p_{e,EEMC}-p_{e,MC})/p_{e,MC}",200,0,20,2000,-1,1);
    
    // RECO VM 
    TH1D* h_VM_mass_REC = new TH1D("h_VM_mass_REC",";VM_{RECO} mass [GeV/c^{2}]",200,0,2);
    TH1D* h_VM_p_REC = new TH1D("h_VM_p_REC",";p_{VM,RECO} [GeV/c]",200,0,20);
    TH1D* h_VM_pt_REC = new TH1D("h_VM_pt_REC",";p_{T,VM,RECO} [GeV/c]",200,0,10);
    TH2D* h_VM_pt_response = new TH2D("h_VM_pt_response","; p_{T,VM,RECO} [GeV/c];p_{T,VM,MC} [GeV/c]",200,0,10,200,0,10);
	TH1D* h_VM_pz_REC = new TH1D("h_VM_pz_REC",";p_{z,VM,RECO} [GeV/c]",200,0,20);
    TH2D* h_VM_Epz_res = new TH2D("h_VM_Epz_res",";(E-p_{z})_{MC} [GeV]; ((E-p_{z})_{REC}-(E-p_{z})_{MC})/(E-p_{z})_{MC}",200,0,20,200,-1,1);
    TH2D* h_VM_Epz_response = new TH2D("h_VM_Epz_response","; (E_{VM,RECO}-p_{z,VM,RECO}) [GeV];(E_{VM,MC}-p_{z,VM,MC} [GeV])",200,0,20,200,0,20);
    TH1D* h_VM_Epz_REC = new TH1D("h_VM_Epz_REC",";(E_{VM,REC}-p_{z,VM,REC}) [GeV]",200,0,20);
    TH2D* h_VM_Epz_migration = new TH2D("h_VM_Epz_migration",";(E_{VM,MC}-p_{z,VM,MC}) [GeV];(E_{VM,RECO}-p_{z,VM,RECO}) [GeV]",200,0,20,200,0,20);
    TH2D* h_VM_pt_migration = new TH2D("h_VM_pt_migration",";p_{T,VM,MC} [GeV/c];p_{T,VM,RECO} [GeV/c]",200,0,10,200,0,10);
    TH2D* h_VM_pt_resolution = new TH2D("h_VM_pt_resolution",";p_{T,MC} [GeV/c]; (p_{T,RECO}-p_{T,MC})/p_{T,MC}",200,0,10,1000,-1,1);
    TH2D* h_VM_pz_resolution = new TH2D("h_VM_pz_resolution",";p_{z,MC} [GeV/c]; (p_{z,RECO}-p_{z,MC})/p_{z,MC}",200,0,20,1000,-1,1);
    TH2D* h_VM_p_resolution = new TH2D("h_VM_p_resolution",";p_{MC} [GeV/c]; (p_{RECO}-p_{MC})/p_{MC}",200,0,20,1000,-1,1);
    TH1D* h_VM_theta_MC = new TH1D("h_VM_theta_MC",";#theta_{VM MC}",100,0,3.14);
    TH1D* h_VM_theta_REC = new TH1D("h_VM_theta_REC",";#theta_{VM RECO}",100,0,3.14);
    
    // RECO position
	TH1D* h_emHits_position_x_REC = new TH1D("h_emHits_position_x_REC","x [mm]",80,-800,800);
	TH1D* h_emHits_position_y_REC = new TH1D("h_emHits_position_y_REC","y [mm]",80,-800,800);
	TH1D* h_emClus_position_x_REC = new TH1D("h_emClus_position_x_REC","x [mm]",80,-800,800);
	TH1D* h_emClus_position_y_REC = new TH1D("h_emClus_position_y_REC","y [mm]",80,-800,800);
    TH2D* h_emClus_position_REC = new TH2D("h_emClus_position_REC",";x (mm);y (mm)",80,-800,800,80,-800,800);
    TH2D* h_emClus_position_REC_cut = new TH2D("h_emClus_position_REC_cut",";x (mm);y (mm)",80,-800,800,80,-800,800);
    TH2D* h_XvsY_hits = new TH2D("h_XvsY_hits",";x [mm]; y [mm]",80,-800,800,80,-800,800);
    TH2D* h_XvsY_clus = new TH2D("h_XvsY_clus",";x [mm]; y [mm]",80,-800,800,80,-800,800);
    TH1D* h_Xclus_minus_Xtrk = new TH1D("h_Xclus_minus_Xtrk",";x_{clus}-x_{trk} [mm]",80,-800,800);
    TH1D* h_Yclus_minus_Ytrk = new TH1D("h_Yclus_minus_Ytrk",";y_{clus}-y_{trk} [mm]",80,-800,800);
    TH1D* h_Xhits_minus_Xtrk = new TH1D("h_Xhits_minus_Xtrk",";x_{hits}-x_{trk} [mm]",80,-800,800);
    TH1D* h_Yhits_minus_Ytrk = new TH1D("h_Yhits_minus_Ytrk",";y_{hits}-y_{trk} [mm]",80,-800,800);
    
    // RECO E-pz
    TH1D* h_Epz_REC = new TH1D("h_Epz_REC", ";(E_{EEMC} - p_{z,EEMC}) [GeV]",100,15,25);
    TH2D* h_Epz_res = new TH2D("h_Epz_res",";(E-p_{z})_{MC} [GeV]; ((E-p_{z})_{REC}-(E-p_{z})_{MC})/(E-p_{z})_{MC}",100,15,25,100,-1,1);
    TH2D* h_Epz_response = new TH2D("h_Epz_response", ";(E_{EEMC} - p_{z,EEMC}) [GeV];(E_{MC}-p_{z,MC}) [GeV]",100,15,25,100,15,25);
    TH2D* h_Epz_migration = new TH2D("h_Epz_migration",";(E_{MC}-p_{z,MC}) [GeV];(E_{RECO} - p_{z,RECO}) [GeV]",100,15,25,100,15,25);
    
    // RECO E/P
    TH1D* h_EcalOverPtrk = new TH1D("h_EcalOverPtrk",";E_{EEMC}/|p|_{trk}",100,0,2);
    TH1D* h_EtrkOverPcal = new TH1D("h_EtrkOverPcal",";E_{trk}/|p|_{EEMC}",100,0,2);
    TH2D* h_EoverP_res = new TH2D("h_EoverP_res",";(E/p)_{MC}; ((E/p)_{REC}-(E/p)_{MC})/(E/p)_{MC}",100,0,2,100,-1,1);
    TH2D* h_EoverP_response = new TH2D("h_EoverP_response","; E_{EEMC}/|p|_{trk};E_{MC}/|p|_{MC}",100,0,2,100,0,2);
    TH2D* h_EoverP_migration = new TH2D("h_EoverP_migration",";E_{MC}/|p|_{MC};E_{EEMC}/|p|_{trk}",100,0,2,100,0,2);

    // rho
    TH1D* h_rho_rec_before_mass_cut = new TH1D("h_rho_rec_before_mass_cut",";m_{#rho,REC} [GeV/c]",200,0,2);
    TH1D* h_rho_rec_after_mass_cut = new TH1D("h_rho_rec_after_mass_cut",";m_{#rho,REC} [GeV/c]",200,0,2);
    TH1D* h_rho_mass = new TH1D("h_rho_mass",";m_{#rho,MC} [GeV/c]",200,0,2);
    
    // RECO t
    TH1D* h_t_REC_EEMC = new TH1D("h_t_REC_EEMC",";|t|_{EEMC} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_EEMC_cut = new TH1D("h_t_REC_EEMC_cut",";|t|_{EEMC} [GeV/c]^{2}; counts",100,0,0.2);
    TH2D* h_t_res_EEMC_cut_percent = new TH2D("h_t_res_EEMC_cut_percent",";|t|_{MC} [GeV/c]^{2}; (|t|_{MC}-|t|_{EEMC})/|t|_{MC}",100,0,0.2,1000,-10,10);
    TH2D* h_t_res_EEMC = new TH2D("h_t_res_EEMC",";|t|_{MC} [GeV/c]^{2}; (|t|_{MC}-|t|_{EEMC}) [GeV/c]^{2}",100,0,0.2,1000,-10,10);
    TH2D* h_t_res = new TH2D("h_t_res",";|t|_{MC} [GeV/c]^{2}; (|t|_{MC}-|t|_{RECO}) [GeV/c]^{2}",100,0,0.2,1000,-10,10);
    TH2D* h_t_res_EEMC_cut = new TH2D("h_t_res_EEMC_cut",";|t|_{MC} [GeV/c]^{2}; (|t|_{MC}-|t|_{EEMC}) [GeV/c]^{2}",100,0,0.2,1000,-10,10);
    TH2D* h_t_migration = new TH2D("h_t_migration",";|t|_{RECO} [GeV/c]^{2};|t|_{MC} [GeV/c]^{2}",100,0,0.2,100,0,0.2);
    TH2D* h_t_response_cut = new TH2D("h_t_response_cut","; |t|_{RECO} [GeV/c]^{2};|t|_{MC} [GeV/c]^{2}",100,0,0.2,100,0,0.2);
    TH2D* h_t_response_EEMC_cut = new TH2D("h_t_response_EEMC_cut","; |t|_{EEMC} [GeV/c]^{2};|t|_{MC} [GeV/c]^{2}",100,0,0.2,100,0,0.2);
    TH2D* h_t_res_proj_percent_pi12 = new TH2D("h_t_res_proj_percent_pi12",";|t|_{MC} [GeV/c]^{2}; (|t|_{MC}-|t|_{RECO})/|t|_{MC}",100,0,0.2,1000,-10,10);
    TH2D* h_t_res_proj_12 = new TH2D("h_t_res_proj_12",";|t|_{MC} [GeV/c]^{2}; (|t|_{MC}-|t|_{RECO}) [GeV/c]^{2}",100,0,0.2,1000,-10,10);
    TH1D* h_t_REC_new_method = new TH1D("h_t_REC_new_method",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wRES = new TH1D("h_t_REC_wRES",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wRES_cut = new TH1D("h_t_REC_wRES_cut",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wRES_cut_pi2 = new TH1D("h_t_REC_wRES_cut_pi2",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wRES_cut_pi3 = new TH1D("h_t_REC_wRES_cut_pi3",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wRES_cut_pi4 = new TH1D("h_t_REC_wRES_cut_pi4",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wRES_cut_pi6 = new TH1D("h_t_REC_wRES_cut_pi6",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wRES_cut_pi9 = new TH1D("h_t_REC_wRES_cut_pi9",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wRES_cut_pi12 = new TH1D("h_t_REC_wRES_cut_pi12",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wRES_cut_pi16 = new TH1D("h_t_REC_wRES_cut_pi16",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wRES_cut_pi20 = new TH1D("h_t_REC_wRES_cut_pi20",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wRES_cut_pi24 = new TH1D("h_t_REC_wRES_cut_pi24",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH2D* h_t_REC_2d = new TH2D("h_t_REC_2d",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_cut = new TH2D("h_t_REC_2d_cut",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_cut2 = new TH2D("h_t_REC_2d_cut2",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wRES = new TH2D("h_t_REC_2d_wRES",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wRES_cut = new TH2D("h_t_REC_2d_wRES_cut",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wRES_cut_pi2 = new TH2D("h_t_REC_2d_wRES_cut_pi2",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wRES_cut_pi3 = new TH2D("h_t_REC_2d_wRES_cut_pi3",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wRES_cut_pi4 = new TH2D("h_t_REC_2d_wRES_cut_pi4",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wRES_cut_pi6 = new TH2D("h_t_REC_2d_wRES_cut_pi6",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wRES_cut_pi9 = new TH2D("h_t_REC_2d_wRES_cut_pi9",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wRES_cut_pi12 = new TH2D("h_t_REC_2d_wRES_cut_pi12",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wRES_cut_pi16 = new TH2D("h_t_REC_2d_wRES_cut_pi16",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wRES_cut_pi20 = new TH2D("h_t_REC_2d_wRES_cut_pi20",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wRES_cut_pi24 = new TH2D("h_t_REC_2d_wRES_cut_pi24",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wCUT = new TH2D("h_t_REC_2d_wCUT",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wCUT_pi2 = new TH2D("h_t_REC_2d_wCUT_pi2",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wCUT_pi3 = new TH2D("h_t_REC_2d_wCUT_pi3",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wCUT_pi4 = new TH2D("h_t_REC_2d_wCUT_pi4",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wCUT_pi6 = new TH2D("h_t_REC_2d_wCUT_pi6",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wCUT_pi9 = new TH2D("h_t_REC_2d_wCUT_pi9",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wCUT_pi12 = new TH2D("h_t_REC_2d_wCUT_pi12",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wCUT_pi16 = new TH2D("h_t_REC_2d_wCUT_pi16",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wCUT_pi20 = new TH2D("h_t_REC_2d_wCUT_pi20",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wCUT_pi24 = new TH2D("h_t_REC_2d_wCUT_pi24",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    
        /*-----------------------------------------------------------------------------------------
            - t distribution with cut only (for testing)
            - distributions should all look the same since no resolution has been added
        ------------------------------------------------------------------------------------------*/
    TH1D* h_t_REC_wCUT_pi2 = new TH1D("h_t_REC_wCUT_pi2",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wCUT_pi3 = new TH1D("h_t_REC_wCUT_pi3",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wCUT_pi4 = new TH1D("h_t_REC_wCUT_pi4",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wCUT_pi6 = new TH1D("h_t_REC_wCUT_pi6",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wCUT_pi9 = new TH1D("h_t_REC_wCUT_pi9",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wCUT_pi12 = new TH1D("h_t_REC_wCUT_pi12",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wCUT_pi16 = new TH1D("h_t_REC_wCUT_pi16",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wCUT_pi20 = new TH1D("h_t_REC_wCUT_pi20",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wCUT_pi24 = new TH1D("h_t_REC_wCUT_pi24",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    
    // Compare cuts and extra histograms for analysis
    TH1D* h_Epz_beforeCut = new TH1D("h_Epz_beforeCut", ";(E_{EEMC} - p_{z,trk}) [GeV]",100,15,25);
    TH1D* h_EoverP_afterCut = new TH1D("h_EoverP_afterCut",";E_{EEMC}/|p|_{trk}",100,0,2);
    TH1D* h_EoverP_beforeCut = new TH1D("h_EoverP_beforeCut",";E_{EEMC}/|p|_{trk}",100,0,2);
    TH1D* h_Q2_beforeCut = new TH1D("h_Q2_beforeCut",";Q^{2} [GeV/c]^{2}",100,0,10);
    TH1D* h_Q2_afterCut = new TH1D("h_Q2_afterCut",";Q^{2} [GeV/c]^{2}",100,0,10);
    TH1D* h_y_beforeCut = new TH1D("h_y_beforeCut",";y",100,0.01,0.85);
    TH1D* h_y_afterCut = new TH1D("h_y_afterCut",";y",100,0.01,0.85);
    TH1D* h_VM_Epz_MC_before = new TH1D("h_VM_Epz_MC_before",";(E_{VM,MC}-p_{z,VM,MC}) [GeV]",200,0,20);
    TH1D* h_VM_pt_MC_before = new TH1D("h_VM_pt_MC_before",";p_{T,VM,MC} [GeV/c]",200,0,10);
    TH1D* h_VM_mass_MC_before = new TH1D("h_VM_mass_MC_before",";m_{VM,MC} [GeV/c]",200,0,2);
    TH1D* h_VM_pz_MC_before = new TH1D("h_VM_pz_MC_before",";p_{z,VM,MC} [GeV/c]",200,0,20);
    TH1D* h_VM_p_MC_before = new TH1D("h_VM_p_MC_before",";|p|_{VM,MC} [GeV/c]",200,0,20);
    TH1D* h_VM_mass_REC_before = new TH1D("h_VM_mass_REC_before",";m_{VM,MC} [GeV/c]",200,0,2);
    TH1D* h_VM_pt_REC_before = new TH1D("h_VM_pt_REC_before",";p_{T,VM,MC} [GeV/c]",200,0,10);
    TH1D* h_VM_pz_REC_before = new TH1D("h_VM_pz_REC_before",";p_{z,VM,MC} [GeV/c]",200,0,20);
    TH1D* h_VM_p_REC_before = new TH1D("h_VM_p_REC_before",";|p|_{VM,MC} [GeV/c]",200,0,20);
    TH1D* h_VM_Epz_REC_before = new TH1D("h_VM_Epz_REC_before",";(E_{VM,MC}-p_{z,VM,MC}) [GeV]",200,0,20);
    TH1D* h_solid_angle_REC = new TH1D("h_solid_angle_REC",";#Omega between #theta_{e'} and #phi_{e'}",1000,-0.1,0.1);
    TH1D* h_angle_resolution = new TH1D("h_angle_resolution",";Angular Resolution",1000,-0.1,0.1);
    TH2D* h_phi_response = new TH2D("h_phi_response",";#phi_{RECO}; #phi_{MC}",100,-3.14,3.14,100,-3.14,3.14);
    TH2D* h_theta_response = new TH2D("h_theta_response",";#theta_{RECO}; #theta_{MC}",100,0,3.14,100,0,3.14);
    TH1D* h_omega_MC = new TH1D("h_omega_MC",";#Omega_{MC}",100,-3.14,3.14);
    TH1D* h_omega_REC = new TH1D("h_omega_REC",";#Omega_{REC}",100,-3.14,3.14);
    TH1F* h_Nevents = new TH1F("h_Nevents", "Total Number of Events", 1, 0, 1);
	TH2D* h_e_pz_trk_res = new TH2D("h_e_pz_trk_res",";p_{z,e,MC} [GeV/c]; (p_{z,e,trk}-p_{z,e,MC})/p_{z,e,MC}",200,0,20,2000,-1,1);
	TH2D* h_e_p_trk_res = new TH2D("h_e_p_trk_res",";|p|_{e,MC} [GeV/c]; (|p|_{e,trk}-|p|_{e,MC})/|p|_{e,MC}",200,0,20,2000,-1,1);
    TH2D* h_Q2_vs_x_REC = new TH2D("h_Q2_vs_x_REC",";x_{RECO}; Q^{2}_{RECO} [GeV/c]^{2}",100,0,0.5,100,0,10);
    TH2D* h_Q2_vs_x_MC = new TH2D("h_Q2_vs_x_MC",";x_{MC}; Q^{2}_{MC} [GeV/c]^{2}",100,0,0.5,100,0,10);
    TH2D* h_Q2_res_beforecut = new TH2D("h_Q2_res_beforecut",";Q^{2}_{MC} [GeV/c]^{2}; (Q^{2}_{MC}-Q^{2}_{EEMC})/Q^{2}_{MC}",100,0,10,1000,-1,1);
    TH2D* h_Q2_response_beforecut = new TH2D("h_Q2_response_beforecut","; Q^{2}_{EEMC} [GeV/c]^{2};Q^{2}_{MC} [GeV/c]^{2}",100,0,10,100,0,10);     
    TH2D* h_EvsP_REC = new TH2D("h_EvsP_REC",";|p|_{trk} [GeV]; E_{EEMC} [GeV]",100,0,20,100,0,20);

    // detectors
    TH1D* hEEMC_hits = new TH1D("hEEMC_hits", "Number of EEMC clusters per event", 50, 0, 50);
    TH1D* hZDC_hits  = new TH1D("hZDC_hits", "ZDC cluster multiplicity;clusters/event;Counts", 5, -0.5, 4.5);
    TH1D* hRP_hits   = new TH1D("hRP_hits",   "Number of RP hits per event", 50, 0, 50);
    TH1D* hOMD_hits  = new TH1D("hOMD_hits",  "Number of OMD hits per event", 50, 0, 50);
    TH1D* hEEMC_energy = new TH1D("hEEMC_energy", "EEMC Cluster Energy", 50,0,30);
    TH1D* hZDC_energy = new TH1D("hZDC_energy", "ZDC Energy;E [GeV];Counts", 200, 0, 200);
    TH2D* hEEMC_xy = new TH2D("hEEMC_xy", "EEMC hit map;x [mm];y [mm]", 100, -200, 200, 100, -200, 200);
    TH2D* hZDC_xy = new TH2D("hZDC_xy", "ZDC XY;x [mm];y [mm]", 200, -1200, -600, 200, -300, 300);
    TH2D* hRP_xy = new TH2D("hRP_xy", "RP XY;x [mm];y [mm]", 200, -1200, -600, 200, -300, 300);
    TH2D* hOMD_xy = new TH2D("hOMD_xy", "OMD XY;x [mm];y [mm]", 200, -1200, -600, 200, -300, 300);
    TH2D* hZDC_hits_xy = new TH2D("hZDC_hits_xy", "ZDC Hit Positions;x [mm];y [mm]",200, -1200, 600, 200, -300, 300);

    chain->GetEntries();
    while (tree_reader.Next()) 
    {    
         /*---------------
            Event Loop
        ----------------*/
        h_Nevents->Fill(0.5);
        cout << "Total events processed" << h_Nevents->GetEntries() << endl;

        event.clusters_eemc.clear();
        event.clusters_zdc.clear();
        event.hit_rp.clear();
        event.hit_omd.clear();

        /*---------------------------
            Get incoherent events
        ----------------------------*/
       
		//EEMC
        for(int iclus=0;iclus<em_energy_array.GetSize();iclus++)
        {
      		Cluster_EEMC cluster;
			cluster.energy=em_energy_array[iclus];
			cluster.x=em_x_array[iclus];
			cluster.y=em_y_array[iclus];
    	    event.clusters_eemc.push_back(cluster);
    	}

    	//ZDC
        for (int ihit = 0; ihit < zdc_x_array.GetSize(); ihit++) 
        {
            hZDC_hits_xy->Fill(zdc_x_array[ihit], zdc_y_array[ihit]);
        }

        double Ecal_total = 0.0;
        double Hcal_total = 0.0;

        for (int i = 0; i < Ecal_zdc_energy_array.GetSize(); ++i) 
        {
            Ecal_total += Ecal_zdc_energy_array[i];
        }
        for (int i = 0; i < Hcal_zdc_energy_array.GetSize(); ++i) 
        {
            Hcal_total += Hcal_zdc_energy_array[i];
        }
        bool has_ECAL_cluster_energy = (Ecal_total > 0.0);
        bool has_HCAL_cluster_energy = (Hcal_total > 0.0);

        double veto_energy = 0.0;
        if (has_ECAL_cluster_energy) 
        {
            veto_energy = Ecal_total;
        }
        else if (has_HCAL_cluster_energy) 
        {
            veto_energy = Hcal_total;  
        }

        double threshold = 21.0; // GeV
        if (veto_energy > threshold) 
        {
            Cluster_ZDC cluster;
            cluster.energy = veto_energy;
            if (has_ECAL_cluster_energy) 
            {
                cluster.x = Ecal_zdc_x_array[0];
                cluster.y = Ecal_zdc_y_array[0];
                cluster.z = Ecal_zdc_z_array[0];
            } 
            else if (has_HCAL_cluster_energy) 
            {
                cluster.x = Hcal_zdc_x_array[0];
                cluster.y = Hcal_zdc_y_array[0];
                cluster.z = Hcal_zdc_z_array[0];
            }
            event.clusters_zdc.push_back(cluster);
            hZDC_energy->Fill(cluster.energy);
            hZDC_xy->Fill(cluster.x, cluster.y);
        }
        hZDC_hits->Fill(event.clusters_zdc.size());
        
    	//RP
    	for(int ihit=0;ihit<rp_x_array.GetSize();ihit++)
        {
      		Hit_RP hit;
			hit.x=rp_x_array[ihit];
			hit.y=rp_y_array[ihit];
			hit.z=rp_z_array[ihit];
    	    event.hit_rp.push_back(hit);
    	}

    	//OMD
    	for(int ihit=0;ihit<omd_x_array.GetSize();ihit++)
        {
      		Hit_OMD hit;
			hit.x=omd_x_array[ihit];
			hit.y=omd_y_array[ihit];
			hit.z=omd_z_array[ihit];
    	    event.hit_omd.push_back(hit);
    	}            

        hEEMC_hits->Fill(event.clusters_eemc.size());
        hRP_hits->Fill(event.hit_rp.size());
        hOMD_hits->Fill(event.hit_omd.size());

        for (auto& cluster : event.clusters_eemc) 
        {
            hEEMC_energy->Fill(cluster.energy);
            hEEMC_xy->Fill(cluster.x, cluster.y, cluster.energy);
        }
        //if(event.clusters_zdc.size()>0) continue;
    	if(event.hit_rp.size()>0) continue;
    	//if(event.hit_omd.size()>0) continue;
        cout << "I survived" << endl;
        
    	// Beam particles
    	TLorentzVector ebeam(0,0,0,0);
    	TLorentzVector pbeam(0,0,0,0);
        TLorentzVector abeam(0,0,0,0);
        TLorentzVector vmMC(0,0,0,0);
        TLorentzVector rho_MC(0,0,0,0);
        TLorentzVector vmMC_before(0,0,0,0);
    	TLorentzVector kplusMC(0,0,0,0);
    	TLorentzVector kminusMC(0,0,0,0);
        TLorentzVector piplusMC(0,0,0,0);
    	TLorentzVector piminusMC(0,0,0,0);
        TLorentzVector hfsMC(0,0,0,0);
        TLorentzVector particleMC(0,0,0,0);

    	//MC level
    	TLorentzVector scatMC(0,0,0,0);
    	int mc_elect_index=-1;
    	double maxPt=-99.;
        int incoherent=0;

        // loop over all MC particles in the event
    	for(int imc=0;imc<mc_px_array.GetSize();imc++)
        {
            TVector3 mctrk(static_cast<double>(mc_px_array[imc]), 
                static_cast<double>(mc_py_array[imc]), 
                static_cast<double>(mc_pz_array[imc]));
            particleMC.SetVectM(mctrk,MASS_PION); // assume pions;
    		if(mc_genStatus_array[imc]==4) // 4 is Sartre.
            {
                /*--------------------------------------
                    -genStatus = 4 is Sartre
                    -genStatus = 1 is stable particle
                    -PDG 11 = electron
                    -PDG 2212 = proton
                    -PDG 321 = kaon
                    -PDG 211 = pion
                    -All other PDG set as ion
                ---------------------------------------*/
                if(mc_pdg_array[imc]==11) ebeam.SetVectM(mctrk, MASS_ELECTRON);
				if(mc_pdg_array[imc]==2212) pbeam.SetVectM(mctrk, MASS_PROTON);
                else
                {
                    double MASS_A = mc_mass_array[imc];
                    abeam.SetVectM(mctrk,MASS_A);
                    //abeam.SetVectM(mctrk,MASS_AU197);
                    cout<<"genStatus=4 candidate: PDG="<<mc_pdg_array[imc]<<" Mass="<< abeam.M()<<endl;
                }
            }
         	if(mc_genStatus_array[imc]!=1) continue; 
    		if(mc_pdg_array[imc]==11 && mctrk.Perp()>maxPt)
            {
                /*---------------------------------------------------
                        electron with highest pT = scattered electron
                    ----------------------------------------------------*/
    			maxPt=mctrk.Perp();
    			mc_elect_index=imc;
    			scatMC.SetVectM(mctrk,mc_mass_array[imc]); 
    		}
            // if particle is stable and kaon PDG == 321
    		if(mc_pdg_array[imc]==321 
                && mc_genStatus_array[imc]==1) kplusMC.SetVectM(mctrk,MASS_KAON);
    		if(mc_pdg_array[imc]==-321 
                && mc_genStatus_array[imc]==1) kminusMC.SetVectM(mctrk,MASS_KAON);

            /*-----------------------------------------------------
                    for rho background
                        - uncomment to set daughters as pions
                        - be sure to comment out the kaons above
                ----------------------------------------------------*/
            //if particle is stable and pion PDG == 211
    		//if(mc_pdg_array[imc]==211 
            //    && mc_genStatus_array[imc]==1) piplusMC.SetVectM(mctrk,MASS_PION);
    		//if(mc_pdg_array[imc]==-211 
            //    && mc_genStatus_array[imc]==1) piminusMC.SetVectM(mctrk,MASS_PION);
            cout << "PDG for MC VM is " << mc_pdg_array[imc] << " particle mass is "<< mc_mass_array[imc] << endl;
            if(mctrk.Eta()>3.5&&(mc_mass_array[imc]>0.9383||mc_mass_array[imc]<0.938))incoherent++;
            if(imc!=mc_elect_index) hfsMC += particleMC;
        }
        cout << "# particles eta>3.5 " << incoherent << endl; 
        //if(incoherent) continue;

        // checks
        cout<<"A energy: "<<abeam.E()<<" p Energy: "<<pbeam.E()<<" e Energy: "<<ebeam.E()<<endl;
        cout<<"Electron Beam: "<<" px: "<<ebeam.Px()<<" py: "<<ebeam.Py()<<" pz: "<<ebeam.Pz()<< " E: "<<ebeam.E()<< "e beam theta: " << ebeam.Theta() << endl;
		cout<<"p Beam: "<<" px: "<<pbeam.Px()<<" py: "<<pbeam.Py()<<" pz: "<<pbeam.Pz()<<" E: "<<pbeam.E()<< endl;
        cout <<"A beam: "<<" px: "<< abeam.Px()<<" py: "<< abeam.Py()<<" pz: "<< abeam.Pz()<< " E: "<< abeam.E()<< endl;
		cout <<"Scattered Electron: "<<" px: "<<scatMC.Px()<<" py: "<<scatMC.Py()<<" pz: "<< scatMC.Pz()<<" E: "<< scatMC.E()<<endl;
    
    	// protection.
    	if(ebeam.E()==abeam.E() && ebeam.E()==0) 
        {
    		cout << "problem with MC incoming beams" << endl;
    		continue;
    	}
        
        double phi_MC = scatMC.Phi();
        double theta_MC = scatMC.Theta();
        double eta_MC = scatMC.Eta();
        double e_pT_MC = scatMC.Pt();
        double e_pz_MC = scatMC.Pz();
        double e_p_MC = scatMC.P();
        double energy_MC = scatMC.E();
        double EoverP_MC = scatMC.E()/scatMC.P();
        double Epz_MC = (scatMC+hfsMC).E()-(scatMC+hfsMC).Pz();
        //double Epz_MC = scatMC.E()-scatMC.Pz();
        
        h_phi_MC->Fill(phi_MC);
        h_theta_MC->Fill(theta_MC);
        h_eta_MC->Fill(eta_MC);
        h_e_pt_MC->Fill(e_pT_MC);
		h_e_pz_MC->Fill(e_pz_MC);
		h_e_p_MC->Fill(e_p_MC);
		h_energy_MC->Fill(energy_MC);
        h_EoverP_MC->Fill(EoverP_MC);
        h_Epz_MC->Fill(Epz_MC);

    	TLorentzVector qbeam = ebeam-scatMC; // p_e - p_e'
        // electron method
    	double Q2 = -(qbeam).Mag2();  
    	double pq = abeam.Dot(qbeam); 
    	double y = pq/abeam.Dot(ebeam);
        double BjorkenX = Q2/(2*abeam.Dot(qbeam));
                
        h_Q2_beforeCut->Fill(Q2);
        h_y_beforeCut->Fill(y);
        h_x_beforeCut->Fill(BjorkenX);

        /*----------------------------------------------------
            - If studying rho background, use vmMC as pions
            - Comment out the vmMC with kaons
        -------------------------------------------------------*/
		vmMC=kplusMC+kminusMC;
        //vmMC=piplusMC+piminusMC;
        cout << "made it this far before VM energy which is " << vmMC.E() << endl;
        
    	if(vmMC.E()==0) continue;
        double vm_Epz_MC = vmMC.E()-vmMC.Pz();
        double vm_pT_MC = vmMC.Pt();
        double vm_mass_MC = vmMC.M();
        double vm_pz_MC = vmMC.Pz();
        double vm_p_MC = vmMC.P();
        
        h_VM_Epz_MC_before->Fill(vm_Epz_MC);
        h_VM_pt_MC_before->Fill(vm_pT_MC);
        h_VM_mass_MC_before->Fill(vm_mass_MC);
        h_VM_pz_MC_before->Fill(vm_pz_MC);
        h_VM_p_MC_before->Fill(vm_p_MC);

        double method_E_before = -(qbeam-vmMC).Mag2();  // t = -(p_e - p_e' - p_v)^2
    	h_t_MC_before->Fill(method_E_before);

        // MC level phase space cut
    	if(Q2<1.||Q2>10.) continue;
    	if(y<0.01||y>0.85) continue;

        cout << "made it past DIS cuts" << endl;
        
		// t dist
        double t_MC = 0;
    	if(vmMC.E()!=0 && fabs(vmMC.Rapidity())<3.5)
    	{
            cout << "made it past VM energy and rapidity cuts" << endl;
            h_phi_MC_after->Fill(phi_MC);
            h_theta_MC_after->Fill(theta_MC);
            h_eta_MC_after->Fill(eta_MC);
            h_e_pt_MC_after->Fill(e_pT_MC);
    		h_e_pz_MC_after->Fill(e_pz_MC);
    		h_e_p_MC_after->Fill(e_p_MC);
    		h_energy_MC_after->Fill(energy_MC);
            h_EoverP_MC_after->Fill(EoverP_MC);
            h_Epz_MC_after->Fill(Epz_MC);
            h_EvsP_MC->Fill(e_p_MC,energy_MC);
            
            h_Q2_e->Fill(Q2);
        	h_y_e->Fill(y);
            h_x->Fill(BjorkenX);
            h_Q2_vs_x_MC->Fill(BjorkenX,Q2);

            double theta_VM_MC = vmMC.Theta();
            h_VM_theta_MC->Fill(theta_VM_MC);
            h_VM_mass_MC->Fill(vm_mass_MC);
        	h_VM_pt_MC->Fill(vm_pT_MC);
    		h_VM_pz_MC->Fill(vm_pz_MC);
    		h_VM_p_MC->Fill(vm_p_MC);
            h_VM_Epz_MC->Fill(vm_Epz_MC);
            
    		double method_E = -(qbeam-vmMC).Mag2();  // t = -(p_e - p_e' - p_v)^2
    		t_MC = method_E;
    		h_t_MC->Fill(method_E);
          
    		/*------------------------------------------------
                projection method
                - Same as method L but takes a theta_max cut
            -------------------------------------------------*/
    		TLorentzVector T = -ebeam + scatMC + vmMC;  // true t (method E)
            
    		// define nHat direction -> direction perpendicular to electron scattering plane
    		TVector3 eScattered_momentum = scatMC.Vect(); 
    		TVector3 e_momentum = ebeam.Vect();
    		TVector3 nHat = e_momentum.Cross(eScattered_momentum); // nHat = p_e x p_e'
    		nHat = nHat.Unit();
        
            // define projected VM (y-component of t)
    		double pv = vmMC.Vect().Dot(nHat);
    		TVector3 py_vector = nHat*pv;  // VM in nhat direction
    		TLorentzVector py(py_vector.X(),py_vector.Y(),0,0);
            double ty = -(py).Mag2();
            double qy = sqrt(ty);
        
            // define x and z-components
            TLorentzVector px = T - py; 
            TLorentzVector pz(0,0,px.Z(),px.E());
    		TLorentzVector px_true = T - py - pz;
            double tx_true = -(px_true).Mag2();
            double qx_true = sqrt(tx_true);
            double tz = -(pz).Mag2();
            double t_total = tx_true+ty+tz; // no resolution added so should still be t_MC 

    		h_t_REC_2d->Fill(qx_true,qy);
            h_t_REC_new_method->Fill(t_total);
       
    		// apply cut with qz and E subtracted out of qx-component (just for testing, distributions should match method E)
    		double theta = atan(fabs(qx_true)/fabs(qy));
    		if(fabs(theta)<PI/2) 
    		{
    			h_t_REC_wCUT_pi2->Fill(t_total);
    			h_t_REC_2d_wCUT_pi2->Fill(qx_true,qy);
    		}
    		if(fabs(theta)<PI/3) 
    		{
    			h_t_REC_wCUT_pi3->Fill(t_total);
    			h_t_REC_2d_wCUT_pi3->Fill(qx_true,qy);
    		}
    		if(fabs(theta)<PI/4) 
    		{
    			h_t_REC_wCUT_pi4->Fill(t_total);
    			h_t_REC_2d_wCUT_pi4->Fill(qx_true,qy);
    		}
    		if(fabs(theta)<PI/6) 
    		{
    			h_t_REC_wCUT_pi6->Fill(t_total);
    			h_t_REC_2d_wCUT_pi6->Fill(qx_true,qy);
    		}
    		if(fabs(theta)<PI/9) 
    		{
    			h_t_REC_wCUT_pi9->Fill(t_total);
    			h_t_REC_2d_wCUT_pi9->Fill(qx_true,qy);
    		}
    		if(fabs(theta)<PI/12) 
    		{
    			h_t_REC_wCUT_pi12->Fill(t_total);
    			h_t_REC_2d_wCUT_pi12->Fill(qx_true,qy);
    		}
            if(fabs(theta)<PI/16) 
    		{
                h_t_REC_wCUT_pi16->Fill(t_total);
    			h_t_REC_2d_wCUT_pi16->Fill(qx_true,qy);
    		}
            if(fabs(theta)<PI/20) 
    		{
                h_t_REC_wCUT_pi20->Fill(t_total);
    			h_t_REC_2d_wCUT_pi20->Fill(qx_true,qy);
    		}
            if(fabs(theta)<PI/24) 
    		{
                h_t_REC_wCUT_pi24->Fill(t_total);
    			h_t_REC_2d_wCUT_pi24->Fill(qx_true,qy);
    		}
    	}


    	// rec level 
        /*--------------------------------------------------------------------------------
            -find highest energy cluster in negative endcap to find scattered electron
        ---------------------------------------------------------------------------------*/
    	double maxEnergy=-99.;
    	double xpos=-999.;
    	double ypos=-999.;
                
        // find highest energy cluster in negative endcap to find scattered electron
    	for(int iclus=0;iclus<em_energy_array.GetSize();iclus++)
        {
        	if(em_energy_array[iclus]>maxEnergy)
            {
                /*---------------------------------------------------------------
                    -check if cluster has higher energy
                    -if it does, define as maxEnergy and define that position
                ----------------------------------------------------------------*/
        		maxEnergy=em_energy_array[iclus];
    			xpos=em_x_array[iclus]; // built from EE clusters
    			ypos=em_y_array[iclus];
    		}
    	}
                
    	/*----------------------------------------------------------------------------------
            find highest energy hits (above threshold) and store that position and index
        -----------------------------------------------------------------------------------*/
    	double maxHitEnergy=0.01;  // threshold 10 MeV
    	double xhitpos=-999.;
    	double yhitpos=-999.;
    	int hit_index=-1;
                
        // individual calorimeter cells that registered energy
    	for(int ihit=0;ihit<emhits_energy_array.GetSize();ihit++)
        {	
    		if(emhits_energy_array[ihit]>maxHitEnergy)
            {
                // find highest energy hits (above threshold) and store that position and index
    			maxHitEnergy=emhits_energy_array[ihit];
			    xhitpos=emhits_x_array[ihit];
			    yhitpos=emhits_y_array[ihit];
			    hit_index=ihit;
    		}
    	}

    	/*----------------------------------------------------------------------------------------------
            -sum over all 3x3 towers around the leading tower
            -agregate nearby hits to form a shower to reconstruct particles energy and impact point
        -----------------------------------------------------------------------------------------------*/
    	double xClus=xhitpos*maxHitEnergy;
    	double yClus=yhitpos*maxHitEnergy;
    	for(int ihit=0;ihit<emhits_energy_array.GetSize();ihit++)
        {
    		double hitenergy=emhits_energy_array[ihit];
    		double x=emhits_x_array[ihit];
    		double y=emhits_y_array[ihit];
            // sum over energies within 70 mm of maxHitEnergy
    		double d=sqrt( (x-xhitpos)*(x-xhitpos) + (y-yhitpos)*(y-yhitpos));
    		if(d<70. && ihit!=hit_index && hitenergy>0.01)  
            {
    			maxHitEnergy+=hitenergy;
    			xClus+=x*hitenergy;
    			yClus+=y*hitenergy;
    		}
    	}
	
    	// weighted average cluster position.
    	xClus = xClus/maxHitEnergy; // built from hits
    	yClus = yClus/maxHitEnergy;
    	double radius=sqrt(xClus*xClus+yClus*yClus);
                
        // remove clusters outside of detector bounds
    	if(radius>550.) continue; 
    	double clusEnergy=1.044*maxHitEnergy; // 4.4% energy calibration.
    	h_energy_REC_EEMC->Fill(clusEnergy);
    	h_emClus_position_REC->Fill(xClus,yClus); 

        double min_distance = 1e6;
        double xtrk = -999.;
        double ytrk = -999.;
                
        // get positions from calorimeter track
        for(int itrk=0; itrk<emCal_trk_x_array.GetSize(); itrk++) 
        {
            double dx = xpos - emCal_trk_x_array[itrk];
            double dy = ypos - emCal_trk_y_array[itrk];
            double dist = sqrt(dx*dx + dy*dy);
            if(dist<min_distance) 
            {
                min_distance = dist;
                xtrk = emCal_trk_x_array[itrk];
                ytrk = emCal_trk_y_array[itrk];
            }
        }
        
        double x_clus_trk_diff = xpos-xtrk; // built from EE cluster
        double x_hit_trk_diff = xClus-xtrk; // built from hits
        double y_clus_trk_diff = ypos-ytrk;
        double y_hit_trk_diff = yClus-ytrk;
        
        // energy resolution
        double res= (energy_MC-clusEnergy)/energy_MC;
        h_energy_res_EEMC->Fill(energy_MC,res);
        h_energy_response_EEMC->Fill(energy_MC,clusEnergy);
        
        /*---------------------------------------------------------------------
                -association of rec level scat' e
                -sim_id = index of simulated MC particle 
                -rec_id = index of reconstructed particle
        		-find the rec_id that matches MC scattered electron (sim_id)
            ----------------------------------------------------------------------*/
        int rec_elect_index=-1;
        for(int i=0;i<sim_id.GetSize();i++)
        {
        	if(sim_id[i]==mc_elect_index) // mc_elect_index found when defining scatMC
            {
        		rec_elect_index = rec_id[i]; 
        	}
        }

        // reconstructed particles
        TLorentzVector scatMCmatchREC(0,0,0,0);
        TLorentzVector scatREC(0,0,0,0);
        TLorentzVector scatClusEREC(0,0,0,0);
        TLorentzVector hfs(0,0,0,0);
        TLorentzVector particle(0,0,0,0);
        TLorentzVector kplusREC(0,0,0,0);
        TLorentzVector kminusREC(0,0,0,0);
        TLorentzVector piplusREC(0,0,0,0);
        TLorentzVector piminusREC(0,0,0,0);
        TLorentzVector vmREC(0,0,0,0);
        TLorentzVector rho_REC(0,0,0,0);
        TLorentzVector vmREC_before(0,0,0,0);

        double maxP=-1.;
        
        // track loop
        for(int itrk=0;itrk<reco_pz_array.GetSize();itrk++)
        {
        	TVector3 trk(reco_px_array[itrk],reco_py_array[itrk],reco_pz_array[itrk]);
            int q = reco_charge_array[itrk];  
            cout << "reco charge " << q << endl;
        	if(rec_elect_index!=-1 && itrk==rec_elect_index && trk.Mag()>maxP && q<0) 
            {
                cout << "made it past e requirements" << endl;
                // if stable and matches electron index
        		scatMCmatchREC.SetVectM(trk,MASS_ELECTRON); 
        		maxP=trk.Mag();
        		scatREC.SetVectM(trk,MASS_ELECTRON);
        		// use emcal energy to define 4 vector
        		double p = sqrt(clusEnergy*clusEnergy- MASS_ELECTRON*MASS_ELECTRON );
                double pt = TMath::Sin(scatREC.Theta())*p;
        		scatClusEREC.SetPtEtaPhiM(pt,scatREC.Eta(),scatREC.Phi(),MASS_ELECTRON);
        	}
        }
    	
        double eta_REC = scatREC.Eta();
        double phi_REC = scatREC.Phi();
        double theta_REC = scatREC.Theta();
        
        h_eta_REC_EEMC->Fill(eta_REC);
        h_phi_REC_EEMC->Fill(phi_REC);
        h_theta_REC_EEMC->Fill(theta_REC);
        h_theta_response_EEMC->Fill(theta_REC,theta_MC);
        h_theta_diff->Fill(theta_REC-theta_MC);
        h_phi_diff->Fill(phi_REC-phi_MC);
        h_eta_diff->Fill(eta_REC-eta_MC);

        double e_pT_REC_cal = scatClusEREC.Pt();
        double e_pz_REC_cal = scatClusEREC.Pz();
        double e_p_REC_cal = scatClusEREC.P();
        double e_pT_REC_trk = scatMCmatchREC.Pt();
        double e_pz_REC_trk = scatMCmatchREC.Pz();
        double e_p_REC_trk = scatMCmatchREC.P();
            
        // default clustering position
    	h_emClus_position_x_REC->Fill(xpos);
    	h_emClus_position_y_REC->Fill(ypos);
        h_XvsY_clus->Fill(xpos,ypos);
        
    	// self clustering position
    	h_emHits_position_x_REC->Fill(xClus);
    	h_emHits_position_y_REC->Fill(yClus);
        h_XvsY_hits->Fill(xClus,yClus);

        // track positions
        h_XvsY_trk->Fill(xtrk,ytrk);
        h_trk_position_x_REC->Fill(xtrk);
        h_trk_position_y_REC->Fill(ytrk);

        h_Xclus_minus_Xtrk->Fill(x_clus_trk_diff);
        h_Yclus_minus_Ytrk->Fill(y_clus_trk_diff);
        h_Xhits_minus_Xtrk->Fill(x_hit_trk_diff);
        h_Yhits_minus_Ytrk->Fill(y_hit_trk_diff);
        
    	// track-base e' energy
        double energy_REC_trk = scatMCmatchREC.E();
    	double res_trk= (energy_MC-energy_REC_trk)/energy_MC;
        h_energy_res_trk->Fill(energy_MC,res_trk);
    	h_energy_REC_trk->Fill(energy_REC_trk);

        int extra_hfs = 0;
        int eta_out = 0;
        int incoherent_rec = 0;
        /*-----------------------------------------------------------------
            for PID
            - If studying rho production, uncomment this block for PID
            - Comment out the other block
                - Other block is reco of kaons with no PID
        ------------------------------------------------------------------*/
        /*for(int itrk=0;itrk<reco_pz_array.GetSize();itrk++)
        {   
            TVector3 trk(reco_px_array[itrk],reco_py_array[itrk],reco_pz_array[itrk]);
        	particle.SetVectM(trk,MASS_PION); // assume pions;
            //------------------------------------------------------
            //        -exclude e' from hadronic final state and kaons
            //        -select phi->kk daughters;
            //-------------------------------------------------------
        	if(itrk!=rec_elect_index) 
            {
            	hfs += particle;
                extra_hfs++;
            	h_eta_REC_trk->Fill(trk.Eta());
            	if(fabs(trk.Eta())<3.0)
                {
                    // if particle is kaon PDG == 321
            		if(reco_pdg_array[itrk]==321) kplusREC.SetVectM(trk,MASS_KAON);
            		if(reco_pdg_array[itrk]==-321) kminusREC.SetVectM(trk,MASS_KAON);
                    cout << "PDG for reco VM is " << reco_pdg_array[itrk] << endl;
            	}
        	}
            //incoherent
        	if(trk.Eta()>3.5 || extra_hfs>2) incoherent_rec++;
            if(trk.Eta()>3.5) eta_out++;
        }*/

        // for no PID
        for(int itrk=0;itrk<reco_pz_array.GetSize();itrk++)
        {
        	TVector3 trk(reco_px_array[itrk],reco_py_array[itrk],reco_pz_array[itrk]);
        	particle.SetVectM(trk,MASS_PION); // assume pions;
            //------------------------------------------------------
            //        -exclude e' from hadronic final state and kaons
            //        -select phi->kk daughters;
            //-------------------------------------------------------
        	if(itrk!=rec_elect_index) 
            {
            	hfs += particle;
                extra_hfs++;
            	h_eta_REC_trk->Fill(trk.Eta());
            	if(fabs(trk.Eta())<3.0)
                {
            		if(reco_charge_array[itrk]>0) kplusREC.SetVectM(trk,MASS_KAON);
            		if(reco_charge_array[itrk]<0) kminusREC.SetVectM(trk,MASS_KAON);
            	}
        	}
            //incoherent
        	if(trk.Eta()>3.5 || extra_hfs>2) incoherent_rec++;
            if(trk.Eta()>3.5) eta_out++;
        }
        cout << "# reco paticles eta > 3.5 " << eta_out << "# extra HFS particles " << extra_hfs << endl;  
        cout << "incoherent: " << incoherent_rec << endl;
        //if(incoherent_rec) continue;
        
        // 4vector of VM;
        if(kplusREC.E()!=0. && kminusREC.E()!=0.)
        {
         	vmREC=kplusREC+kminusREC;
        }
        double vm_mass_REC = vmREC.M();
        double vm_pT_REC = vmREC.Pt();
        double vm_pz_REC = vmREC.Pz();
        double vm_p_REC = vmREC.P();
        double vm_Epz_REC = vmREC.E()-vmREC.Pz();
       
        h_VM_mass_REC_before->Fill(vm_mass_REC);
        h_VM_pt_REC_before->Fill(vm_pT_REC);
        h_VM_pz_REC_before->Fill(vm_pz_REC);
        h_VM_p_REC_before->Fill(vm_p_REC);
        h_VM_Epz_REC_before->Fill(vm_Epz_REC);
        
        // cluster-base DIS kine;
        TLorentzVector qbeamREC=ebeam-scatClusEREC; 
        
        // electron method
        double Q2REC = -(qbeamREC).Mag2();
        double pqREC = abeam.Dot(qbeamREC); 
        double yREC = pqREC/abeam.Dot(ebeam);;
        double BjorkenX_REC = Q2REC/(2*abeam.Dot(qbeamREC));

        h_Q2REC_e_EEMC->Fill(Q2REC);
        h_Q2_res_beforecut->Fill(Q2,(Q2-Q2REC)/Q2);
        h_Q2_response_beforecut->Fill(Q2REC,Q2);
        h_yREC_e_EEMC->Fill(yREC);
        h_x_REC->Fill(BjorkenX_REC);

        double energy_REC = scatClusEREC.E();
        double Epz_REC = (scatClusEREC+hfs).E() - (scatMCmatchREC+hfs).Pz();
        double EoverP_REC = energy_REC/e_p_REC_trk;
        h_Epz_beforeCut->Fill(Epz_REC);
        h_EoverP_beforeCut->Fill(EoverP_REC);

        // Event selection:
        if( Epz_REC>15&&Epz_REC<25 ) 
        {
            if( EoverP_REC>0.9&&EoverP_REC<1.2 )
            {
                if(Q2REC>1.&&Q2REC<10.) 
                {
                    if(yREC>0.01&&yREC<0.85) 
                    {
                        // VM rec
                        if(vmREC.E()!=0) 
                        {
                            // select phi mass and rapidity window 
                            if( fabs(vmREC.M()-1.02)<0.02&& fabs(vmREC.Rapidity())<3.5)
                            {
                                // Reco e'
                                double x_REC = sin(theta_REC)*cos(phi_REC);
                                double y_REC = sin(theta_REC)*sin(phi_REC);
                                double z_REC = cos(theta_REC);

                                // MC e'
                                double x_MC = sin(theta_MC)*cos(phi_MC);
                                double y_MC = sin(theta_MC)*sin(phi_MC);
                                double z_MC = cos(theta_MC);

                                // Get angle between (Omega)
                                double dot_product = x_REC*x_MC + y_REC*y_MC + z_REC*z_MC;
                                // Clamp to avoid domain errors in acos
                                dot_product = std::max(-1.0, std::min(1.0, dot_product));
                                double angle_resolution_rad = acos(dot_product); // in radians
                                double solid_angle_REC = 2 * M_PI * (1 - cos(angle_resolution_rad));
                                double omega_MC = z_MC/dot_product;
                                double omega_REC = z_REC/dot_product;
                                
                                // omega
                                h_omega_MC->Fill(omega_MC);
                                h_omega_REC->Fill(omega_REC);
                                h_angle_resolution->Fill(angle_resolution_rad);
                                h_solid_angle_REC->Fill(solid_angle_REC);

                                // Q2 emcal reco resolution
                                h_Q2_afterCut->Fill(Q2REC);
                            	double resQ2 = (Q2-Q2REC)/Q2;
                            	h_Q2_res->Fill(Q2,resQ2);
                                h_Q2_res2->Fill((Q2REC-Q2)/Q2);
                                h_Q2_response->Fill(Q2REC,Q2);
                                h_dQ2overQ2_REC->Fill(resQ2);
                                h_Q2_migration->Fill(Q2,Q2REC);
                                h_Q2_vs_x_REC->Fill(BjorkenX_REC,Q2REC);
                                
                                // x resolution
                                h_x_afterCut->Fill(BjorkenX_REC);
                                double resx = (BjorkenX-BjorkenX_REC)/BjorkenX;
                                h_x_res_cut->Fill(BjorkenX,resx);
                                h_x_res_cut2->Fill((BjorkenX_REC-BjorkenX)/BjorkenX);
                                h_x_response_cut->Fill(BjorkenX_REC,BjorkenX);       
                                h_x_migration->Fill(BjorkenX,BjorkenX_REC);       
        
                                // y emcal reco resolution
                                double resy = (y-yREC)/y;
                                h_y_res->Fill(y,resy);
                                h_y_res2->Fill((yREC-y)/y);
                                h_y_response->Fill(yREC,y);
                                h_dyOvery_REC->Fill(resy);
                                h_y_migration->Fill(y,yREC);
                                h_y_afterCut->Fill(yREC);

                                // VM theta
                                double theta_VM_REC = vmREC.Theta();
                                h_VM_theta_REC->Fill(theta_VM_REC);
                                h_VM_mass_REC->Fill(vm_mass_REC);
                                
                                // VM E-pz
                                double resEpzVM = (vm_Epz_MC-vm_Epz_REC)/vm_Epz_MC;
                                h_VM_Epz_REC->Fill(vm_Epz_REC);
                                h_VM_Epz_response->Fill(vm_Epz_REC,vm_Epz_MC);
                                h_VM_Epz_migration->Fill(vm_Epz_MC,vm_Epz_REC);
                                h_VM_Epz_res->Fill(vm_Epz_MC,resEpzVM);

                                // VM pT
                                h_VM_pt_REC->Fill(vm_pT_REC);
                                h_VM_pt_response->Fill(vm_pT_REC,vm_pT_MC);    
                                double respTVM = (vm_pT_MC-vm_pT_REC)/vm_pT_MC;
                                h_VM_pt_resolution->Fill(vm_pT_MC,respTVM);
                                h_VM_pt_migration->Fill(vm_pT_MC,vm_pT_REC);

                                // VM pz 
                                h_VM_pz_REC->Fill(vm_pz_REC);
                                double respzVM = (vm_pz_MC-vm_pz_REC)/vm_pz_MC;
                                h_VM_pz_resolution->Fill(vm_pz_MC,respzVM);

                                // VM |p| resolution
                                h_VM_p_REC->Fill(vm_p_REC);
                                double respVM = (vm_p_MC-vm_p_REC)/vm_p_MC;
                                h_VM_p_resolution->Fill(vm_p_MC,respVM);    
                                
                                // e' p from cal
                                h_e_pt_REC_EEMC->Fill(e_pT_REC_cal);
                            	h_e_pz_REC_EEMC->Fill(e_pz_REC_cal);
                            	h_e_p_REC_EEMC->Fill(e_p_REC_cal);
                                h_e_pt_res->Fill(e_pT_MC,(e_pT_REC_cal-e_pT_MC)/e_pT_MC);
                                h_e_pz_res->Fill(e_pz_MC,(e_pz_REC_cal-e_pz_MC)/e_pz_MC);
                                h_e_p_res->Fill(e_p_MC,(e_p_REC_cal-e_p_MC)/e_p_MC);
                                h_e_pz_response->Fill(e_pz_REC_cal,e_pz_MC);

                                // e' p from track
                                h_e_pt_REC_trk->Fill(e_pT_REC_trk);
                            	h_e_pz_REC_trk->Fill(e_pz_REC_trk);
                            	h_e_p_REC_trk->Fill(e_p_REC_trk);
                                h_e_pz_trk_res->Fill(e_pz_MC,(e_pz_REC_trk-e_pz_MC)/e_pz_MC);
                                h_e_p_trk_res->Fill(e_p_MC,(e_p_REC_trk-e_p_MC)/e_p_MC);

                                // E_cal vs p_trk
                                h_EvsP_REC->Fill(e_p_REC_trk,energy_REC);
                                
                            	// E-pz scat' e
                                double resEpz = (Epz_MC-Epz_REC)/Epz_MC;
                                h_Epz_res->Fill(Epz_MC,resEpz);
                                h_Epz_REC->Fill(Epz_REC);
                                h_Epz_response->Fill(Epz_REC,Epz_MC);
                                h_Epz_migration->Fill(Epz_MC,Epz_REC);

                            	// E over p
                                double resEoP = (EoverP_MC-EoverP_REC)/EoverP_MC;
                                h_EoverP_res->Fill(EoverP_MC,resEoP);
                                h_EoverP_afterCut->Fill(EoverP_REC);
                                h_EoverP_response->Fill(EoverP_REC,EoverP_MC);
                                h_EoverP_migration->Fill(EoverP_MC,EoverP_REC);

                                // eta
                                h_eta_REC_EEMC_after->Fill(eta_REC);
                                h_eta_diff_after->Fill(eta_REC-eta_MC);

                                // theta
                                h_theta_REC_EEMC_after->Fill(theta_REC);
                                h_theta_response_EEMC_after->Fill(theta_REC,theta_MC);
                                h_theta_diff_after->Fill(theta_REC-theta_MC);
                                double restheta = (theta_MC-theta_REC)/theta_MC;
                                h_theta_resolution->Fill(theta_MC,restheta);
                                h_theta_migration->Fill(theta_MC,theta_REC);

                                // phi
                                h_phi_REC_EEMC_after->Fill(phi_REC);
                                h_phi_diff_after->Fill(phi_REC-phi_MC);
                                double resphi = (phi_MC-phi_REC)/phi_MC;
                                h_phi_resolution->Fill(phi_MC,resphi);
                                h_phi_response->Fill(phi_REC,phi_MC);

                                // energy
                                double resE= (energy_MC-energy_REC)/energy_MC;
                                h_energy_res_EEMC_after->Fill(energy_MC, resE);
                                h_energy_response_EEMC_after->Fill(energy_REC,energy_MC);
                                h_energy_REC_EEMC_after->Fill(energy_REC);
                                h_energy_migration->Fill(energy_MC,energy_REC);
                                h_emClus_position_REC_cut->Fill(xClus,yClus); 
                            
                                // 2 versions: track and energy cluster:
                                double t_trk_REC = giveme_t_method_L(ebeam,scatMCmatchREC,abeam,vmREC);  // method L (track e')
                                double t_REC = giveme_t_method_L(ebeam,scatClusEREC,abeam,vmREC); // method L (EEMC e')
                                h_t_REC_trk_cut->Fill(t_trk_REC);
                                h_t_REC_EEMC_cut->Fill(t_REC);
                                h_t_response_EEMC_cut->Fill(t_REC,t_MC);
                                h_t_response_trk_cut->Fill(t_trk_REC,t_MC);

                                // t track resolution 
                        		double resttrk = (t_MC-t_trk_REC)/t_MC;
                        		h_t_res_trk_cut->Fill(t_MC, resttrk);
	
                            	// t EEMC resolution;
                        		double res_percent = (t_MC-t_REC)/t_MC;
                        		h_t_res_EEMC_cut_percent->Fill(t_MC,res_percent);
                                double rest = (t_MC-t_REC);
                        		h_t_res_EEMC_cut->Fill(t_MC,rest);
                                
                        		// projection method
                                TLorentzVector T_rec = giveme_t_new_method(ebeam,scatClusEREC,abeam,vmREC); // method L (EEMC e')
          
                        		// define nHat direction
                        		TVector3 eScattered_momentum_rec = scatClusEREC.Vect(); 
                        		TVector3 e_momentum = ebeam.Vect();
                        		TVector3 nHat_rec = e_momentum.Cross(eScattered_momentum_rec); // nHat = p_e x p_e'
                        		nHat_rec = nHat_rec.Unit();

                                // define projected VM (y-component)
                        		double pv_rec = vmREC.Vect().Dot(nHat_rec);
                                TVector3 py_vector_rec = nHat_rec*pv_rec;
                        		TLorentzVector py_rec(py_vector_rec.X(),py_vector_rec.Y(),0,0);
                                double ty_rec = -(py_rec).Mag2();
                                double qy_rec = sqrt(ty_rec);
                                
                                // define x and z-components 
                                TLorentzVector px_rec = T_rec - py_rec; 
                                TLorentzVector pz_rec(0,0,px_rec.Z(),px_rec.E());
                        		TLorentzVector px_rec_true = T_rec - py_rec - pz_rec;
                                double tz_rec = -(pz_rec).Mag2();
                                double tx_rec = -(px_rec_true).Mag2();
                        		double qx_rec = sqrt(tx_rec);

                                // plot t = tx+ty+t//
                        		double t_total_rec = tx_rec+ty_rec+tz_rec;
                        		h_t_REC_2d_wRES_cut->Fill(qx_rec,qy_rec);
                        		h_t_REC_wRES_cut->Fill(t_total_rec);

                        		// apply cut with // subtracted out of t
                        		double theta_rec = atan(fabs(qx_rec)/fabs(qy_rec));
                        		if(fabs(theta_rec)<PI/2)
                        		{
                        			h_t_REC_wRES_cut_pi2->Fill(t_total_rec);
                        			h_t_REC_2d_wRES_cut_pi2->Fill(qx_rec,qy_rec);
                        		}
                        		if(fabs(theta_rec)<PI/3)
                        		{
                        			h_t_REC_wRES_cut_pi3->Fill(t_total_rec);
                        			h_t_REC_2d_wRES_cut_pi3->Fill(qx_rec,qy_rec);
                        		}
                        		if(fabs(theta_rec)<PI/4)
                        		{
                        			h_t_REC_wRES_cut_pi4->Fill(t_total_rec);
                        			h_t_REC_2d_wRES_cut_pi4->Fill(qx_rec,qy_rec);
                        		}
                        		if(fabs(theta_rec)<PI/6)
                        		{
                        			h_t_REC_wRES_cut_pi6->Fill(t_total_rec);
                        			h_t_REC_2d_wRES_cut_pi6->Fill(qx_rec,qy_rec);
                        		}
                        		if(fabs(theta_rec)<PI/9)
                        		{
                        			h_t_REC_wRES_cut_pi9->Fill(t_total_rec);
                        			h_t_REC_2d_wRES_cut_pi9->Fill(qx_rec,qy_rec);
                        		}
                        		if(fabs(theta_rec)<PI/12)
                        		{
                        			h_t_REC_wRES_cut_pi12->Fill(t_total_rec);
                        			h_t_REC_2d_wRES_cut_pi12->Fill(qx_rec,qy_rec);
                                    h_t_response_cut->Fill(t_total_rec,t_MC);
                                    double restcut = (t_MC-t_total_rec)/t_MC;
                                	h_t_res->Fill(t_MC,restcut);
                                    h_t_migration->Fill(t_MC,t_total_rec);
                            		double res_percent_proj_pi12 = (t_MC-t_total_rec)/t_MC;
                            		h_t_res_proj_percent_pi12->Fill(t_MC,res_percent_proj_pi12);
                                    double res_t_proj_pi12 = (t_MC-t_total_rec);
                            		h_t_res_proj_12->Fill(t_MC,res_t_proj_pi12);
                                }
                            	if(fabs(theta_rec)<PI/16)
                            	{
                            		h_t_REC_wRES_cut_pi16->Fill(t_total_rec);
                            		h_t_REC_2d_wRES_cut_pi16->Fill(qx_rec,qy_rec);
                            	}
                                if(fabs(theta_rec)<PI/20)
                            	{
                            		h_t_REC_wRES_cut_pi20->Fill(t_total_rec);
                            		h_t_REC_2d_wRES_cut_pi20->Fill(qx_rec,qy_rec);
                            	}
                                if(fabs(theta_rec)<PI/24)
                            	{
                            		h_t_REC_wRES_cut_pi24->Fill(t_total_rec);
                            		h_t_REC_2d_wRES_cut_pi24->Fill(qx_rec,qy_rec);
                            	}
                            }
                        }
                    }
                }
            }
        }
    }

    output->Write();
    output->Close();
    return 0;
}