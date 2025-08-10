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
    float x, y, z;
};

struct Hit_OMD 
{
    float x, y, z;
};

struct Event 
{
    vector<Cluster_EEMC> clusters_eemc;
    vector<Cluster_ZDC> clusters_zdc;
    vector<Hit_RP> hit_rp;
    vector<Hit_OMD> hit_omd;

    void clear() 
    {
        clusters_eemc.clear();
        clusters_zdc.clear();
        hit_rp.clear();
        hit_omd.clear();
    }
};

int diffractive_vm_full_analysis(TString rec_file, TString outputfile)
{	
    // read in multiple root files at a time	
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

    // Uncomment this block to run single files at a time
    // and comment out everything above
    /*
    TString name_of_input = (TString) rec_file;
    std::cout << "Input file = " << name_of_input << endl;
    auto tree = new TChain("events");
    tree->Add(name_of_input);
    TTreeReader tree_reader(tree);       // !the tree reader
    */

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

    // Reconstructed particle arrays for each reconstructed particle
    TTreeReaderArray<float> reco_px_array = {tree_reader, "ReconstructedChargedParticles.momentum.x"};
    TTreeReaderArray<float> reco_py_array = {tree_reader, "ReconstructedChargedParticles.momentum.y"};
    TTreeReaderArray<float> reco_pz_array = {tree_reader, "ReconstructedChargedParticles.momentum.z"};
    TTreeReaderArray<float> reco_charge_array = {tree_reader, "ReconstructedChargedParticles.charge"};
    TTreeReaderArray<unsigned int> rec_id = {tree_reader, "ReconstructedChargedParticleAssociations.recID"};
    TTreeReaderArray<unsigned int> sim_id = {tree_reader, "ReconstructedChargedParticleAssociations.simID"};

    // Reconstructed Far Forward ZDC Neutrals 
    TTreeReaderArray<float> zdc_x_array = {tree_reader, "ReconstructedFarForwardZDCNeutrals.referencePoint.x"};
    TTreeReaderArray<float> zdc_y_array = {tree_reader, "ReconstructedFarForwardZDCNeutrals.referencePoint.y"};
    TTreeReaderArray<float> zdc_z_array = {tree_reader, "ReconstructedFarForwardZDCNeutrals.referencePoint.z"};
    TTreeReaderArray<float> zdc_energy_array = {tree_reader, "ReconstructedFarForwardZDCNeutrals.energy"};

    // Forward Roman Pot Rec Hits
    TTreeReaderArray<float> rp_x_array = {tree_reader, "ForwardRomanPotRecHits.position.x"};
    TTreeReaderArray<float> rp_y_array = {tree_reader, "ForwardRomanPotRecHits.position.y"};
    TTreeReaderArray<float> rp_z_array = {tree_reader, "ForwardRomanPotRecHits.position.z"};

    // Forward Off Momentum Tracker Rec Hits
    TTreeReaderArray<float> omd_x_array = {tree_reader, "ForwardOffMTrackerRecHits.position.x"};
    TTreeReaderArray<float> omd_y_array = {tree_reader, "ForwardOffMTrackerRecHits.position.y"};
    TTreeReaderArray<float> omd_z_array = {tree_reader, "ForwardOffMTrackerRecHits.position.z"};

    //Reconstructed all central + B0 tracks
    //not sure if i want these or not??
    /*TTreeReaderArray<float> reco_px_array = {tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.x"};
    TTreeReaderArray<float> reco_py_array = {tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.y"};
    TTreeReaderArray<float> reco_pz_array = {tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.z"};
    TTreeReaderArray<float> reco_charge_array = {tree_reader, "ReconstructedTruthSeededChargedParticles.charge"};*/

    TString output_name_dir = outputfile+"_output.root";
    cout << "Output file = " << output_name_dir << endl;
    TFile* output = new TFile(output_name_dir,"RECREATE");

    TTree* outputTree = new TTree("miniTree", "Tree with structured event data");
    Event* event = nullptr;
    outputTree->Branch("event", &event);

    // MC events
    TH1D* h_t_MC = new TH1D("h_t_MC",";|t|_{MC} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_Q2_e = new TH1D("h_Q2_e",";Q^{2}_{MC} [GeV/c]^{2}",100,1,10);
	TH1D* h_y_e = new TH1D("h_y_e",";y_{MC}",100,0.01,0.85);
	TH1D* h_energy_MC = new TH1D("h_energy_MC",";E_{MC} [GeV]",100,0,20);
	TH1D* h_eta_MC = new TH1D("h_eta_MC",";#eta_{MC}",100,-3.14,3.14);
	TH1D* h_e_pt_MC = new TH1D("h_e_pt_MC",";p_{T,e,MC} [GeV/c]",200,0,2);
	TH1D* h_e_pz_MC = new TH1D("h_e_pz_MC",";p_{z,e,MC} [GeV/c]",200,0,2);
	TH1D* h_e_p_MC = new TH1D("h_e_p_MC",";p_{e,MC} [GeV/c]",200,0,2);
	TH1D* h_phi_MC = new TH1D("h_phi_MC",";#phi_{MC}",100,-3.14,3.14);
	TH1D* h_theta_MC = new TH1D("h_theta_MC",";#theta_{MC}",100,-3.14,3.14);
	// RECO EEMC
    TH1D* h_Q2REC_e_EEMC = new TH1D("h_Q2REC_e_EEMC",";Q^{2}_{EEMC} [GeV/c]^{2}",100,1,10);
	TH1D* h_yREC_e_EEMC = new TH1D("h_yREC_e_EEMC",";y_{EEMC}",100,0.01,0.85);
	TH1D* h_energy_REC_EEMC = new TH1D("h_energy_REC_EEMC",";E_{EEMC} [GeV]",100,0,20);
	TH1D* h_eta_REC_EEMC = new TH1D("h_eta_REC_EEMC",";#eta_{EEMC}",100,-3.14,3.14);
	TH1D* h_e_pt_REC_EEMC = new TH1D("h_e_pt_REC_EEMC",";p_{T,e,EEMC} [GeV/c]",200,0,2);
	TH1D* h_e_pz_REC_EEMC = new TH1D("h_e_pz_REC_EEMC",";p_{z,e,EEMC} [GeV/c]",200,0,2);
	TH1D* h_e_p_REC_EEMC = new TH1D("h_e_p_REC_EEMC",";p_{e,EEMC} [GeV/c]",200,0,2);
	TH1D* h_phi_REC_EEMC = new TH1D("h_phi_REC_EEMC",";#phi_{EEMC}",100,-3.14,3.14);
    TH1D* h_theta_REC_EEMC = new TH1D("h_theta_REC_EEMC",";#theta_{EEMC}",100,-3.14,3.14);
	// RECO track
    TH1D* h_Q2REC_e_trk = new TH1D("h_Q2REC_e_trk",";Q^{2}_{trk} [GeV/c]^{2}",100,1,10);
    TH1D* h_yREC_e_trk = new TH1D("h_yREC_e_trk",";y_{trk}",100,0.01,0.85);
	TH1D* h_energy_REC_trk = new TH1D("h_energy_REC_trk",";E_{trk} [GeV]",100,0,20);
	TH1D* h_eta_REC_trk = new TH1D("h_eta_REC_trk",";#eta_{trk}",100,-3.14,3.14);
	TH1D* h_e_pt_REC_trk = new TH1D("h_e_pt_REC_trk",";p_{T,e,trk} [GeV/c]",200,0,2);
	TH1D* h_e_pz_REC_trk = new TH1D("h_e_pz_REC_trk",";p_{z,e,trk} [GeV/c]",200,0,2);
	TH1D* h_e_p_REC_trk = new TH1D("h_e_p_REC_trk",";p_{e,trk} [GeV/c]",200,0,2);
    // VM 
    TH1D* h_VM_mass_MC = new TH1D("h_VM_mass_MC",";VM_{MC} mass [GeV/c^{2}]",200,0,4);
    TH1D* h_VM_pt_MC = new TH1D("h_VM_pt_MC",";p_{T,VM,MC} [GeV/c]",200,0,2);
	TH1D* h_VM_pz_MC = new TH1D("h_VM_pz_MC",";p_{z,VM,MC} [GeV/c]",200,0,2);
	TH1D* h_VM_p_MC = new TH1D("h_VM_p_MC",";p_{VM,MC} [GeV/c]",200,0,2);
    TH1D* h_VM_mass_REC = new TH1D("h_VM_mass_REC",";VM_{RECO} mass [GeV/c^{2}]",200,0,4);
    TH1D* h_VM_pt_REC = new TH1D("h_VM_pt_REC",";p_{T,VM,RECO} [GeV/c]",200,0,2);
    TH2D* h_VM_pt_response = new TH2D("h_VM_pt_response","; p_{T,VM,RECO} [GeV/c];p_{T,VM,MC} [GeV/c]",200,0,2,200,0,2);
	TH1D* h_VM_pz_REC = new TH1D("h_VM_pz_REC",";p_{z,VM,RECO} [GeV/c]",200,0,2);
    TH2D* h_VM_Epz_response = new TH2D("h_VM_Epz_response","; (E_{VM,RECO}-p_{z,VM,RECO}) [GeV];(E_{VM,MC}-p_{z,VM,MC} [GeV])",100,0,20,100,0,20);
    TH1D* h_VM_Epz_MC = new TH1D("h_VM_Epz_MC",";(E_{VM,MC}-p_{z,VM,MC}) [GeV]",100,0,20);
    TH1D* h_VM_Epz_REC = new TH1D("h_VM_Epz_REC",";(E_{VM,REC}-p_{z,VM,REC}) [GeV]",100,0,20);
	TH1D* h_VM_p_REC = new TH1D("h_VM_p_REC",";p_{VM,RECO} [GeV/c]",200,0,2);
    // Position
	TH1D* h_emHits_position_x_REC = new TH1D("h_emHits_position_x_REC","x [mm]",80,-800,800);
	TH1D* h_emHits_position_y_REC = new TH1D("h_emHits_position_y_REC","y [mm]",80,-800,800);
	TH1D* h_emClus_position_x_REC = new TH1D("h_emClus_position_x_REC","x [mm]",80,-800,800);
	TH1D* h_emClus_position_y_REC = new TH1D("h_emClus_position_y_REC","y [mm]",80,-800,800);
	// Analysis hists
    TH2D* h_Q2_res = new TH2D("h_Q2_res",";Q^{2}_{MC} [GeV/c]^{2}; (Q^{2}_{MC}-Q^{2}_{EEMC})/Q^{2}_{MC}",100,1,10,1000,-1,1);
    TH2D* h_Q2_response = new TH2D("h_Q2_response","; Q^{2}_{EEMC} [GeV/c]^{2};Q^{2}_{MC} [GeV/c]^{2}",100,1,10,1000,1,10);
	TH1D* h_dQ2overQ2_REC = new TH1D("h_dQ2overQ2_REC",";dQ^{2}_{REC}/Q^{2}_REC",100,1,10);
    TH2D* h_y_res = new TH2D("h_y_res",";y_{e,MC} ;(y_{e,MC}-y_{e,EEMC})/y_{e,MC}",100,0.01,0.85,1000,-1,1);
    TH2D* h_y_response = new TH2D("h_y_response"," ; y_{EEMC};y_{MC}",100,0.01,0.85,1000,0.01,0.85);
    TH1D* h_dyOvery_REC = new TH1D("h_dyOvery_REC",";dy/y",100,0.01,0.85);
	TH2D* h_energy_res_EEMC = new TH2D("h_energy_res_EEMC",";E_{MC} [GeV]; (E_{MC}-E_{EEMC})/E_{MC}",100,0,20,1000,-1,1);
    TH2D* h_energy_res_trk = new TH2D("h_energy_res_trk",";E_{MC} [GeV]; (E_{MC}-E_{trk})/E_{MC}",100,-3.14,3.14,1000,-3.14,3.14);
    TH2D* h_energy_response_EEMC = new TH2D("h_energy_response_EEMC","; E_{EEMC} [GeV];E_{MC} [GeV]",100,0,20,1000,0,20);
    TH2D* h_e_pt_res = new TH2D("h_e_pt_res",";p_{T,e,MC} [GeV/c]; (p_{T,e,EEMC}-p_{T,e,MC})/p_{T,e,MC}",200,0,2,200,0,2);
	TH2D* h_e_pz_res = new TH2D("h_e_pz_res",";p_{z,e,MC} [GeV/c]; (p_z,e,EEMC}-p_{z,e,MC})/p_{z,e,MC}",200,0,2,200,0,2);
    TH2D* h_e_pz_response = new TH2D("h_e_pz_response","; p_{z,e,EEMC} [GeV/c];p_{z,e,MC} [GeV/c]",200,0,2,200,0,2);
	TH2D* h_e_p_res = new TH2D("h_e_p_res",";p_{e,MC} [GeV/c]; (p_{e,EEMC}-p_{e,MC})/p_{e,MC}",200,0,2,200,0,2);
    TH2D* h_XvsY_hits = new TH2D("h_XvsY_hits",";x [mm]; y [mm]",80,-800,800,80,-800,800);
    TH2D* h_XvsY_clus = new TH2D("h_XvsY_clus",";x [mm]; y [mm]",80,-800,800,80,-800,800);
    TH2D* h_XvsY_trk = new TH2D("h_XvsY_trk",";x [mm]; y [mm]",80,-800,800,80,-800,800);
    TH1D* h_Xclus_minus_Xtrk = new TH1D("h_Xclus_minus_Xtrk",";x_{clus}-x_{trk} [mm]",80,-80,80);
    TH1D* h_Yclus_minus_Ytrk = new TH1D("h_Yclus_minus_Ytrk",";y_{clus}-y_{trk} [mm]",80,-80,80);
    TH1D* h_Xhits_minus_Xtrk = new TH1D("h_Xhits_minus_Xtrk",";x_{hits}-x_{trk} [mm]",80,-80,80);
    TH1D* h_Yhits_minus_Ytrk = new TH1D("h_Yhits_minus_Ytrk",";y_{hits}-y_{trk} [mm]",80,-80,80);
    TH1D* h_trk_position_x_REC = new TH1D("h_trk_position_x_REC",";x_{trk} [mm]",80,-800,800);
    TH1D* h_trk_position_y_REC = new TH1D("h_trk_position_y_REC",";y_{trk} [mm]",80,-800,800);
    TH2D* h_EvsP_REC = new TH2D("h_EvsP_REC",";p_{trk} [GeV]; E_{EEMC} [GeV]",100,0,20,100,0,20);
    TH2D* h_EvsP_MC = new TH2D("h_EvsP_MC",";p_{MC} [GeV]; E_{MC} [GeV]",100,0,20,100,0,20);
    TH1D* h_Epz_REC = new TH1D("h_Epz_REC", ";(E_{EEMC} - p_{z,EEMC}) [GeV]",100,0,20);
    TH1D* h_Epz_REC_trk = new TH1D("h_Epz_REC_trk", ";(E_{trk} - p_{z,trk}) [GeV]",100,0,20);
    TH1D* h_Epz_MC = new TH1D("h_Epz_MC", ";(E_{MC} - p_{z,MC}) [GeV]",100,0,20);
    TH2D* h_Epz_response = new TH2D("h_Epz_response", ";(E_{EEMC} - p_{z,EEMC}) [GeV];(E_{MC}-p_{z,MC}) [GeV]",100,0,20,100,0,20);
    TH1D* h_EoverP_MC = new TH1D("h_EoverP_MC",";E_{MC}/|p|_{MC}",100,0,20);
    TH1D* h_EcalOverPtrk = new TH1D("h_EcalOverPtrk",";E_{EEMC}/|p|_{trk}",100,0,20);
    TH1D* h_EtrkOverPcal = new TH1D("h_EtrkOverPcal",";E_{trk}/|p|_{EEMC}",100,0,20);
    TH2D* h_EoverP_response = new TH2D("h_EoverP_response","; E_{EEMC}/|p|_{trk};E_{MC}/|p|_{MC}",100,0,20,100,0,20);
    TH2D* h_theta_response_EEMC = new TH2D("h_theta_response_EEMC","; theta_{EEMC} [GeV];theta_{MC} [GeV]",100,-3.14,3.14,100,-3.14,3.14);
    TH1D* h_theta_diff = new TH1D("h_theta_diff",";#theta_{EEMC}-#theta_{MC}",100,-3.14,3.14);
    TH1D* h_phi_diff = new TH1D("h_phi_diff",";#phi_{EEMC}-#phi_{MC}",100,-3.14,3.14);
    TH1D* h_eta_diff = new TH1D("h_eta_diff",";#eta_{EEMC}-#eta_{MC}",100,-3.14,3.14);
    // t
    TH1D* h_t_REC_EEMC = new TH1D("h_t_REC_EEMC",";|t|_{EEMC} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_trk = new TH1D("h_t_REC_trk",";|t|_{trk} [GeV/c]^{2}; counts",100,0,0.2);
    TH2D* h_t_res_EEMC = new TH2D("h_t_res_EEMC",";|t|_{MC} [GeV/c]^{2}; (|t|_{MC}-|t|_{EEMC})/|t|_{MC}",100,0,0.2,1000,-10,10);
    TH2D* h_t_res_trk = new TH2D("h_t_res_trk",";|t|_{MC} [GeV/c]^{2}; (|t|_{MC}-|t|_{trk})/|t|_{MC}",100,0,0.2,1000,-10,10);
    TH2D* h_t_response_EEMC = new TH2D("h_t_response_EEMC","; |t|_{EEMC} [GeV/c]^{2};|t|_{MC} [GeV/c]^{2}",100,0,0.2,1000,0,0.2);
	TH2D* h_t_response_trk = new TH2D("h_t_response_trk","; |t|_{trk} [GeV/c]^{2};|t|_{MC} [GeV/c]^{2}",100,0,0.2,1000,0,0.2);
	TH2D* h_t_response_trk_EEMC = new TH2D("h_t_response_trk_EEMC","; |t|_{trk} [GeV/c]^{2};|t|_{EEMC} [GeV/c]^{2}",100,0,0.2,1000,0,0.2);
    // t distribution new method reco
    TH1D* h_t_REC_new_method = new TH1D("h_t_REC_new_method",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    // t distribution with resolution only
    TH1D* h_t_REC_wRES = new TH1D("h_t_REC_wRES",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    // t distribution with cut only (for testing)
        // distributions should all look the same since no resolution has been added
    TH1D* h_t_REC_wCUT_pi2 = new TH1D("h_t_REC_wCUT_pi2",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wCUT_pi3 = new TH1D("h_t_REC_wCUT_pi3",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wCUT_pi4 = new TH1D("h_t_REC_wCUT_pi4",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wCUT_pi6 = new TH1D("h_t_REC_wCUT_pi6",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wCUT_pi9 = new TH1D("h_t_REC_wCUT_pi9",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wCUT_pi12 = new TH1D("h_t_REC_wCUT_pi12",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wCUT_pi16 = new TH1D("h_t_REC_wCUT_pi16",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wCUT_pi20 = new TH1D("h_t_REC_wCUT_pi20",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wCUT_pi24 = new TH1D("h_t_REC_wCUT_pi24",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    // t distribution with resolution and angle cut (EEMC RECO)
    TH1D* h_t_REC_wRES_cut_pi2 = new TH1D("h_t_REC_wRES_cut_pi2",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wRES_cut_pi3 = new TH1D("h_t_REC_wRES_cut_pi3",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wRES_cut_pi4 = new TH1D("h_t_REC_wRES_cut_pi4",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wRES_cut_pi6 = new TH1D("h_t_REC_wRES_cut_pi6",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wRES_cut_pi9 = new TH1D("h_t_REC_wRES_cut_pi9",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wRES_cut_pi12 = new TH1D("h_t_REC_wRES_cut_pi12",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wRES_cut_pi16 = new TH1D("h_t_REC_wRES_cut_pi16",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wRES_cut_pi20 = new TH1D("h_t_REC_wRES_cut_pi20",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wRES_cut_pi24 = new TH1D("h_t_REC_wRES_cut_pi24",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    // t distribution 2d
    TH2D* h_t_REC_2d = new TH2D("h_t_REC_2d",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    // t distribution 2d with resolution only 
    TH2D* h_t_REC_2d_wRES = new TH2D("h_t_REC_2d_wRES",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    // t distribution 2d with resolution and angle cut (EEMC RECO)
    TH2D* h_t_REC_2d_wRES_cut_pi2 = new TH2D("h_t_REC_2d_wRES_cut_pi2",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wRES_cut_pi3 = new TH2D("h_t_REC_2d_wRES_cut_pi3",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wRES_cut_pi4 = new TH2D("h_t_REC_2d_wRES_cut_pi4",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wRES_cut_pi6 = new TH2D("h_t_REC_2d_wRES_cut_pi6",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wRES_cut_pi9 = new TH2D("h_t_REC_2d_wRES_cut_pi9",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wRES_cut_pi12 = new TH2D("h_t_REC_2d_wRES_cut_pi12",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wRES_cut_pi16 = new TH2D("h_t_REC_2d_wRES_cut_pi16",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wRES_cut_pi20 = new TH2D("h_t_REC_2d_wRES_cut_pi20",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wRES_cut_pi24 = new TH2D("h_t_REC_2d_wRES_cut_pi24",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    // t distribution 2d with angle cut only
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
    // Compare cuts
    TH1D* h_Epz_afterCut = new TH1D("h_Epz_afterCut", ";(E_{EEMC} - p_{z,trk}) [GeV]",100,0,20);
    TH1D* h_EoverP_afterCut = new TH1D("h_EoverP_afterCut",";E_{EEMC}/|p|_{trk}",100,0,20);
    TH1D* h_Q2_afterCut = new TH1D("h_Q2_afterCut",";Q^{2}_{RECO} [GeV/c]^{2}",100,1,10);
    TH1D* h_yREC_afterCut = new TH1D("h_yREC_afterCut",";y_{RECO}",100,0.01,0.85);

    // RECO inside of position threshold of 20 mm for comparisons
    // RECO EEMC
    TH1D* h_Q2REC_e_EEMC_cut = new TH1D("h_Q2REC_e_EEMC_cut",";Q^{2}_{EEMC} [GeV/c]^{2}",100,1,10);
	TH1D* h_yREC_e_EEMC_cut = new TH1D("h_yREC_e_EEMC_cut",";y_{EEMC}",100,0.01,0.85);
	TH1D* h_energy_REC_EEMC_cut = new TH1D("h_energy_REC_EEMC_cut",";E_{EEMC} [GeV]",100,0,20);
	TH1D* h_eta_REC_EEMC_cut = new TH1D("h_eta_REC_EEMC_cut",";#eta_{EEMC}",100,-3.14,3.14);
	TH1D* h_e_pt_REC_EEMC_cut = new TH1D("h_e_pt_REC_EEMC_cut",";p_{T,e,EEMC} [GeV/c]",200,0,2);
	TH1D* h_e_pz_REC_EEMC_cut = new TH1D("h_e_pz_REC_EEMC_cut",";p_{z,e,EEMC} [GeV/c]",200,0,2);
	TH1D* h_e_p_REC_EEMC_cut = new TH1D("h_e_p_REC_EEMC_cut",";p_{e,EEMC} [GeV/c]",200,0,2);
	TH1D* h_phi_REC_EEMC_cut = new TH1D("h_phi_REC_EEMC_cut",";#phi_{EEMC}",100,-3.14,3.14);
    TH1D* h_theta_REC_EEMC_cut = new TH1D("h_theta_REC_EEMC_cut",";#theta_{EEMC}",100,-3.14,3.14);
	// RECO track
    TH1D* h_Q2REC_e_trk_cut = new TH1D("h_Q2REC_e_trk_cut",";Q^{2}_{trk} [GeV/c]^{2}",100,1,10);
    TH1D* h_yREC_e_trk_cut = new TH1D("h_yREC_e_trk_cut",";y_{trk}",100,0.01,0.85);
	TH1D* h_energy_REC_trk_cut = new TH1D("h_energy_REC_trk_cut",";E_{trk} [GeV]",100,0,20);
	TH1D* h_eta_REC_trk_cut = new TH1D("h_eta_REC_trk_cut",";#eta_{trk}",100,-3.14,3.14);
	TH1D* h_e_pt_REC_trk_cut = new TH1D("h_e_pt_REC_trk_cut",";p_{T,e,trk} [GeV/c]",200,0,2);
	TH1D* h_e_pz_REC_trk_cut = new TH1D("h_e_pz_REC_trk_cut",";p_{z,e,trk} [GeV/c]",200,0,2);
	TH1D* h_e_p_REC_trk_cut = new TH1D("h_e_p_REC_trk_cut",";p_{e,trk} [GeV/c]",200,0,2);
    // VM
    TH1D* h_VM_mass_REC_cut = new TH1D("h_VM_mass_REC_cut",";VM_{RECO} mass [GeV/c^{2}]",200,0,4);
    TH1D* h_VM_pt_REC_cut = new TH1D("h_VM_pt_REC_cut",";p_{T,VM,RECO} [GeV/c]",200,0,2);
    TH2D* h_VM_pt_response_cut = new TH2D("h_VM_pt_response_cut","; p_{T,VM,RECO} [GeV/c];p_{T,VM,MC} [GeV/c]",200,0,2,200,0,2);
	TH1D* h_VM_pz_REC_cut = new TH1D("h_VM_pz_REC_cut",";p_{z,VM,RECO} [GeV/c]",200,0,2);
    TH2D* h_VM_Epz_response_cut = new TH2D("h_VM_Epz_response_cut","; (E_{VM,RECO}-p_{z,VM,RECO}) [GeV];(E_{VM,MC}-p_{z,VM,MC} [GeV])",100,0,20,100,0,20);
    TH1D* h_VM_Epz_MC_cut = new TH1D("h_VM_Epz_MC_cut",";(E_{VM,MC}-p_{z,VM,MC}) [GeV]",100,0,20);
    TH1D* h_VM_Epz_REC_cut = new TH1D("h_VM_Epz_REC_cut",";(E_{VM,REC}-p_{z,VM,REC}) [GeV]",100,0,20);
	TH1D* h_VM_p_REC_cut = new TH1D("h_VM_p_REC_cut",";p_{VM,RECO} [GeV/c]",200,0,2);
    // energy clus
	TH1D* h_emHits_position_x_REC_cut = new TH1D("h_emHits_position_x_REC_cut","x [mm]",80,-800,800);
	TH1D* h_emHits_position_y_REC_cut = new TH1D("h_emHits_position_y_REC_cut","y [mm]",80,-800,800);
	TH1D* h_emClus_position_x_REC_cut = new TH1D("h_emClus_position_x_REC_cut","x [mm]",80,-800,800);
	TH1D* h_emClus_position_y_REC_cut = new TH1D("h_emClus_position_y_REC_cut","y [mm]",80,-800,800);
	// Analysis hists
    TH2D* h_Q2_res_cut = new TH2D("h_Q2_res_cut",";Q^{2}_{MC} [GeV/c]^{2}; (Q^{2}_{MC}-Q^{2}_{EEMC})/Q^{2}_{MC}",100,1,10,1000,-1,1);
    TH2D* h_Q2_response_cut = new TH2D("h_Q2_response_cut","; Q^{2}_{EEMC} [GeV/c]^{2};Q^{2}_{MC} [GeV/c]^{2}",100,1,10,1000,1,10);
	TH1D* h_dQ2overQ2_REC_cut = new TH1D("h_dQ2overQ2_REC_cut",";dQ^{2}_{REC}/Q^{2}_REC",100,1,10);
    TH2D* h_y_res_cut = new TH2D("h_y_res_cut",";y_{e,MC} ;(y_{e,MC}-y_{e,EEMC})/y_{e,MC}",100,0.01,0.85,1000,-1,1);
    TH2D* h_y_response_cut = new TH2D("h_y_response_cut"," ; y_{EEMC};y_{MC}",100,0.01,0.85,1000,0.01,0.85);
    TH1D* h_dyOvery_REC_cut = new TH1D("h_dyOvery_REC_cut",";dy/y",100,0.01,0.85);
	TH2D* h_energy_res_EEMC_cut = new TH2D("h_energy_res_EEMC_cut",";E_{MC} [GeV]; (E_{MC}-E_{EEMC})/E_{MC}",100,0,20,1000,-1,1);
    TH2D* h_energy_response_EEMC_cut = new TH2D("h_energy_response_EEMC_cut","; E_{EEMC} [GeV];E_{MC} [GeV]",100,0,20,1000,0,20);
    TH2D* h_theta_response_EEMC_cut = new TH2D("h_theta_response_EEMC_cut","; theta_{EEMC} [GeV];theta_{MC} [GeV]",100,-3.14,3.14,100,-3.14,3.14);
    TH2D* h_energy_res_trk_cut = new TH2D("h_energy_res_trk_cut",";E_{MC} [GeV]; (E_{MC}-E_{trk})/E_{MC}",100,-3.14,3.14,1000,-3.14,3.14);
    TH2D* h_e_pt_res_cut = new TH2D("h_e_pt_res_cut",";p_{T,e,MC} [GeV/c]; (p_{T,e,EEMC}-p_{T,e,MC})/p_{T,e,MC}",200,0,2,200,0,2);
	TH2D* h_e_pz_res_cut = new TH2D("h_e_pz_res_cut",";p_{z,e,MC} [GeV/c]; (p_z,e,EEMC}-p_{z,e,MC})/p_{z,e,MC}",200,0,2,200,0,2);
    TH2D* h_e_pz_response_cut = new TH2D("h_e_pz_response_cut","; p_{z,e,EEMC} [GeV/c];p_{z,e,MC} [GeV/c]",200,0,2,200,0,2);
	TH2D* h_e_p_res_cut = new TH2D("h_e_p_res_cut",";p_{e,MC} [GeV/c]; (p_{e,EEMC}-p_{e,MC})/p_{e,MC}",200,0,2,200,0,2);
    TH2D* h_XvsY_hits_cut = new TH2D("h_XvsY_hits_cut",";x [mm]; y [mm]",80,-800,800,80,-800,800);
    TH2D* h_XvsY_clus_cut = new TH2D("h_XvsY_clus_cut",";x [mm]; y [mm]",80,-800,800,80,-800,800);
    TH2D* h_XvsY_trk_cut = new TH2D("h_XvsY_trk_cut",";x [mm]; y [mm]",80,-800,800,80,-800,800);
    TH1D* h_Xclus_minus_Xtrk_cut = new TH1D("h_Xclus_minus_Xtrk_cut",";x_{clus}-x_{trk} [mm]",80,-800,800);
    TH1D* h_Yclus_minus_Ytrk_cut = new TH1D("h_Yclus_minus_Ytrk_cut",";y_{clus}-y_{trk} [mm]",80,-800,800);
    TH1D* h_Xhits_minus_Xtrk_cut = new TH1D("h_Xhits_minus_Xtrk_cut",";x_{hits}-x_{trk} [mm]",80,-800,800);
    TH1D* h_Yhits_minus_Ytrk_cut = new TH1D("h_Yhits_minus_Ytrk_cut",";y_{hits}-y_{trk} [mm]",80,-800,800);
    TH1D* h_trk_position_x_REC_cut = new TH1D("h_trk_position_x_REC_cut",";x_{trk} [mm]",80,-800,800);
    TH1D* h_trk_position_y_REC_cut = new TH1D("h_trk_position_y_REC_cut",";y_{trk} [mm]",80,-800,800);
    TH2D* h_EvsP_REC_cut = new TH2D("h_EvsP_REC_cut",";p_{trk} [GeV]; E_{EEMC} [GeV]",100,0,20,100,0,20);
    TH2D* h_EvsP_MC_cut = new TH2D("h_EvsP_MC_cut",";p_{MC} [GeV]; E_{MC} [GeV]",100,0,20,100,0,20);
    TH1D* h_Epz_REC_cut = new TH1D("h_Epz_REC_cut", ";(E_{EEMC} - p_{z,EEMC}) [GeV]",100,0,20);
    TH1D* h_Epz_REC_trk_cut = new TH1D("h_Epz_REC_trk_cut", ";(E_{trk} - p_{z,trk}) [GeV]",100,0,20);
    TH1D* h_Epz_MC_cut = new TH1D("h_Epz_MC_cut", ";(E_{MC} - p_{z,MC}) [GeV]",100,0,20);
    TH2D* h_Epz_response_cut = new TH2D("h_Epz_response_cut", ";(E_{EEMC} - p_{z,EEMC}) [GeV];(E_{MC}-p_{z,MC}) [GeV]",100,0,20,100,0,20);
    TH1D* h_EcalOverPtrk_cut = new TH1D("h_EcalOverPtrk_cut",";E_{EEMC}/|p|_{trk}",100,0,20);
    TH1D* h_EtrkOverPcal_cut = new TH1D("h_EtrkOverPcal_cut",";E_{trk}/|p|_{EEMC}",100,0,20);
    TH2D* h_EoverP_response_cut = new TH2D("h_EoverP_response_cut","; E_{EEMC}/|p|_{trk};E_{MC}/|p|_{MC}",100,0,20,100,0,20);
    TH1D* h_theta_diff_cut = new TH1D("h_theta_diff_cut",";#theta_{EEMC}-#theta_{MC}",100,-3.14,3.14);
    TH1D* h_phi_diff_cut = new TH1D("h_phi_diff_cut",";#phi_{EEMC}-#phi_{MC}",100,-3.14,3.14);
    TH1D* h_eta_diff_cut = new TH1D("h_eta_diff_cut",";#eta_{EEMC}-#eta_{MC}",100,-3.14,3.14);
    // t distributions
    TH1D* h_t_REC_EEMC_cut = new TH1D("h_t_REC_EEMC_cut",";|t|_{EEMC} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_trk_cut = new TH1D("h_t_REC_trk_cut",";|t|_{trk} [GeV/c]^{2}; counts",100,0,0.2);
    TH2D* h_t_res_EEMC_cut = new TH2D("h_t_res_EEMC_cut",";|t|_{MC} [GeV/c]^{2}; (|t|_{MC}-|t|_{EEMC})/|t|_{MC}",100,0,0.2,1000,-10,10);
    TH2D* h_t_res_trk_cut = new TH2D("h_t_res_trk_cut",";|t|_{MC} [GeV/c]^{2}; (|t|_{MC}-|t|_{trk})/|t|_{MC}",100,0,0.2,1000,-10,10);
    TH2D* h_t_response_EEMC_cut = new TH2D("h_t_response_EEMC_cut","; |t|_{EEMC} [GeV/c]^{2};|t|_{MC} [GeV/c]^{2}",100,0,0.2,1000,0,0.2);
	TH2D* h_t_response_trk_cut = new TH2D("h_t_response_trk_cut","; |t|_{trk} [GeV/c]^{2};|t|_{MC} [GeV/c]^{2}",100,0,0.2,1000,0,0.2);
	TH2D* h_t_response_trk_EEMC_cut = new TH2D("h_t_response_trk_EEMC_cut","; |t|_{trk} [GeV/c]^{2};|t|_{EEMC} [GeV/c]^{2}",100,0,0.2,1000,0,0.2);
    // t distribution with resolution only
    TH1D* h_t_REC_wRES_cut = new TH1D("h_t_REC_wRES_cut",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    // t distribution with resolution and angle cut (EEMC RECO)
    TH1D* h_t_REC_wRES_cut_pi2_cut = new TH1D("h_t_REC_wRES_cut_pi2_cut",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wRES_cut_pi3_cut = new TH1D("h_t_REC_wRES_cut_pi3_cut",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wRES_cut_pi4_cut = new TH1D("h_t_REC_wRES_cut_pi4_cut",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wRES_cut_pi6_cut = new TH1D("h_t_REC_wRES_cut_pi6_cut",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wRES_cut_pi9_cut = new TH1D("h_t_REC_wRES_cut_pi9_cut",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wRES_cut_pi12_cut = new TH1D("h_t_REC_wRES_cut_pi12_cut",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wRES_cut_pi16_cut = new TH1D("h_t_REC_wRES_cut_pi16_cut",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wRES_cut_pi20_cut = new TH1D("h_t_REC_wRES_cut_pi20_cut",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wRES_cut_pi24_cut = new TH1D("h_t_REC_wRES_cut_pi24_cut",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    // t distribution 2d
    TH2D* h_t_REC_2d_cut = new TH2D("h_t_REC_2d_cut",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    // t distribution 2d with resolution only 
    TH2D* h_t_REC_2d_wRES_cut = new TH2D("h_t_REC_2d_wRES_cut",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    // t distribution 2d with resolution and angle cut (EEMC RECO)
    TH2D* h_t_REC_2d_wRES_cut_pi2_cut = new TH2D("h_t_REC_2d_wRES_cut_pi2_cut",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wRES_cut_pi3_cut = new TH2D("h_t_REC_2d_wRES_cut_pi3_cut",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wRES_cut_pi4_cut = new TH2D("h_t_REC_2d_wRES_cut_pi4_cut",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wRES_cut_pi6_cut = new TH2D("h_t_REC_2d_wRES_cut_pi6_cut",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wRES_cut_pi9_cut = new TH2D("h_t_REC_2d_wRES_cut_pi9_cut",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wRES_cut_pi12_cut = new TH2D("h_t_REC_2d_wRES_cut_pi12_cut",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wRES_cut_pi16_cut = new TH2D("h_t_REC_2d_wRES_cut_pi16_cut",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wRES_cut_pi20_cut = new TH2D("h_t_REC_2d_wRES_cut_pi20_cut",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wRES_cut_pi24_cut = new TH2D("h_t_REC_2d_wRES_cut_pi24_cut",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    // RECO inside of position threshold of 20 mm for comparisons
    // RECO EEMC
    TH1D* h_Q2REC_e_EEMC_cut2 = new TH1D("h_Q2REC_e_EEMC_cut2",";Q^{2}_{EEMC} [GeV/c]^{2}",100,1,10);
	TH1D* h_yREC_e_EEMC_cut2 = new TH1D("h_yREC_e_EEMC_cut2",";y_{EEMC}",100,0.01,0.85);
	TH1D* h_energy_REC_EEMC_cut2 = new TH1D("h_energy_REC_EEMC_cut2",";E_{EEMC} [GeV]",100,0,20);
	TH1D* h_eta_REC_EEMC_cut2 = new TH1D("h_eta_REC_EEMC_cut2",";#eta_{EEMC}",100,-3.14,3.14);
	TH1D* h_e_pt_REC_EEMC_cut2 = new TH1D("h_e_pt_REC_EEMC_cut2",";p_{T,e,EEMC} [GeV/c]",200,0,2);
	TH1D* h_e_pz_REC_EEMC_cut2 = new TH1D("h_e_pz_REC_EEMC_cut2",";p_{z,e,EEMC} [GeV/c]",200,0,2);
	TH1D* h_e_p_REC_EEMC_cut2 = new TH1D("h_e_p_REC_EEMC_cut2",";p_{e,EEMC} [GeV/c]",200,0,2);
	TH1D* h_phi_REC_EEMC_cut2 = new TH1D("h_phi_REC_EEMC_cut2",";#phi_{EEMC}",100,-3.14,3.14);
    TH1D* h_theta_REC_EEMC_cut2 = new TH1D("h_theta_REC_EEMC_cut2",";#theta_{EEMC}",100,-3.14,3.14);
	// RECO track
    TH1D* h_Q2REC_e_trk_cut2 = new TH1D("h_Q2REC_e_trk_cut2",";Q^{2}_{trk} [GeV/c]^{2}",100,1,10);
    TH1D* h_yREC_e_trk_cut2 = new TH1D("h_yREC_e_trk_cut2",";y_{trk}",100,0.01,0.85);
	TH1D* h_energy_REC_trk_cut2 = new TH1D("h_energy_REC_trk_cut2",";E_{trk} [GeV]",100,0,20);
	TH1D* h_eta_REC_trk_cut2 = new TH1D("h_eta_REC_trk_cut2",";#eta_{trk}",100,-3.14,3.14);
	TH1D* h_e_pt_REC_trk_cut2 = new TH1D("h_e_pt_REC_trk_cut2",";p_{T,e,trk} [GeV/c]",200,0,2);
	TH1D* h_e_pz_REC_trk_cut2 = new TH1D("h_e_pz_REC_trk_cut2",";p_{z,e,trk} [GeV/c]",200,0,2);
	TH1D* h_e_p_REC_trk_cut2 = new TH1D("h_e_p_REC_trk_cut2",";p_{e,trk} [GeV/c]",200,0,2);
    // VM
    TH1D* h_VM_mass_REC_cut2 = new TH1D("h_VM_mass_REC_cut2",";VM_{RECO} mass [GeV/c^{2}]",200,0,4);
    TH1D* h_VM_pt_REC_cut2 = new TH1D("h_VM_pt_REC_cut2",";p_{T,VM,RECO} [GeV/c]",200,0,2);
    TH2D* h_VM_pt_response_cut2 = new TH2D("h_VM_pt_response_cut2","; p_{T,VM,RECO} [GeV/c];p_{T,VM,MC} [GeV/c]",200,0,2,200,0,2);
	TH1D* h_VM_pz_REC_cut2 = new TH1D("h_VM_pz_REC_cut2",";p_{z,VM,RECO} [GeV/c]",200,0,2);
    TH2D* h_VM_Epz_response_cut2 = new TH2D("h_VM_Epz_response_cut2","; (E_{VM,RECO}-p_{z,VM,RECO}) [GeV];(E_{VM,MC}-p_{z,VM,MC} [GeV])",100,0,20,100,0,20);
    TH1D* h_VM_Epz_MC_cut2 = new TH1D("h_VM_Epz_MC_cut2",";(E_{VM,MC}-p_{z,VM,MC}) [GeV]",100,0,20);
    TH1D* h_VM_Epz_REC_cut2 = new TH1D("h_VM_Epz_REC_cut2",";(E_{VM,REC}-p_{z,VM,REC}) [GeV]",100,0,20);
	TH1D* h_VM_p_REC_cut2 = new TH1D("h_VM_p_REC_cut2",";p_{VM,RECO} [GeV/c]",200,0,2);
    // energy clus
	TH1D* h_emHits_position_x_REC_cut2 = new TH1D("h_emHits_position_x_REC_cut2","x [mm]",800,-800,800);
	TH1D* h_emHits_position_y_REC_cut2 = new TH1D("h_emHits_position_y_REC_cut2","y [mm]",800,-800,800);
	TH1D* h_emClus_position_x_REC_cut2 = new TH1D("h_emClus_position_x_REC_cut2","x [mm]",800,-800,800);
	TH1D* h_emClus_position_y_REC_cut2 = new TH1D("h_emClus_position_y_REC_cut2","y [mm]",800,-800,800);
	// Analysis hists
    TH2D* h_Q2_res_cut2 = new TH2D("h_Q2_res_cut2",";Q^{2}_{MC} [GeV/c]^{2}; (Q^{2}_{MC}-Q^{2}_{EEMC})/Q^{2}_{MC}",100,1,10,1000,-1,1);
    TH2D* h_Q2_response_cut2 = new TH2D("h_Q2_response_cut2","; Q^{2}_{EEMC} [GeV/c]^{2};Q^{2}_{MC} [GeV/c]^{2}",100,1,10,1000,1,10);
	TH1D* h_dQ2overQ2_REC_cut2 = new TH1D("h_dQ2overQ2_REC_cut2",";dQ^{2}_{REC}/Q^{2}_REC",100,1,10);
    TH2D* h_y_res_cut2 = new TH2D("h_y_res_cut2",";y_{e,MC} ;(y_{e,MC}-y_{e,EEMC})/y_{e,MC}",100,0.01,0.85,1000,-1,1);
    TH2D* h_y_response_cut2 = new TH2D("h_y_response_cut2","; y_{EEMC};y_{MC}",100,0.01,0.85,1000,0.01,0.85);
    TH1D* h_dyOvery_REC_cut2 = new TH1D("h_dyOvery_REC_cut2",";dy/y",100,0.01,0.85);
	TH2D* h_energy_res_EEMC_cut2 = new TH2D("h_energy_res_EEMC_cut2",";E_{MC} [GeV]; (E_{MC}-E_{EEMC})/E_{MC}",100,0,20,1000,-1,1);
    TH2D* h_energy_response_EEMC_cut2 = new TH2D("h_energy_response_EEMC_cut2","; E_{EEMC} [GeV];E_{MC} [GeV]",100,0,20,1000,0,20);
    TH2D* h_theta_response_EEMC_cut2 = new TH2D("h_theta_response_EEMC_cut2","; theta_{EEMC} [GeV];theta_{MC} [GeV]",100,-3.14,3.14,100,-3.14,3.14);
    TH2D* h_energy_res_trk_cut2 = new TH2D("h_energy_res_trk_cut2",";E_{MC} [GeV]; (E_{MC}-E_{trk})/E_{MC}",100,-3.14,3.14,1000,-3.14,3.14);
    TH2D* h_e_pt_res_cut2 = new TH2D("h_e_pt_res_cut2",";p_{T,e,MC} [GeV/c]; (p_{T,e,EEMC}-p_{T,e,MC})/p_{T,e,MC}",200,0,2,200,0,2);
	TH2D* h_e_pz_res_cut2 = new TH2D("h_e_pz_res_cut2",";p_{z,e,MC} [GeV/c]; (p_z,e,EEMC}-p_{z,e,MC})/p_{z,e,MC}",200,0,2,200,0,2);
    TH2D* h_e_pz_response_cut2 = new TH2D("h_e_pz_response_cut2",";p_{z,e,EEMC} [GeV/c];p_{z,e,MC} [GeV/c]",200,0,2,200,0,2);
	TH2D* h_e_p_res_cut2 = new TH2D("h_e_p_res_cut2",";p_{e,MC} [GeV/c]; (p_{e,EEMC}-p_{e,MC})/p_{e,MC}",200,0,2,200,0,2);
    TH2D* h_XvsY_hits_cut2 = new TH2D("h_XvsY_hits_cut2",";x [mm]; y [mm]",800,-800,800,80,-800,800);
    TH2D* h_XvsY_clus_cut2 = new TH2D("h_XvsY_clus_cut2",";x [mm]; y [mm]",800,-800,800,80,-800,800);
    TH2D* h_XvsY_trk_cut2 = new TH2D("h_XvsY_trk_cut2",";x [mm]; y [mm]",800,-800,800,80,-800,800);
    TH1D* h_Xclus_minus_Xtrk_cut2 = new TH1D("h_Xclus_minus_Xtrk_cut2",";x_{clus}-x_{trk} [cm]",800,-80,80);
    TH1D* h_Yclus_minus_Ytrk_cut2 = new TH1D("h_Yclus_minus_Ytrk_cut2",";y_{clus}-y_{trk} [cm]",800,-80,80);
    TH1D* h_Xhits_minus_Xtrk_cut2 = new TH1D("h_Xhits_minus_Xtrk_cut2",";x_{hits}-x_{trk} [cm]",800,-80,80);
    TH1D* h_Yhits_minus_Ytrk_cut2 = new TH1D("h_Yhits_minus_Ytrk_cut2",";y_{hits}-y_{trk} [cm]",800,-80,80);
    TH1D* h_trk_position_x_REC_cut2 = new TH1D("h_trk_position_x_REC_cut2",";x_{trk} [mm]",800,-800,800);
    TH1D* h_trk_position_y_REC_cut2 = new TH1D("h_trk_position_y_REC_cut2",";y_{trk} [mm]",800,-800,800);
    TH2D* h_EvsP_REC_cut2 = new TH2D("h_EvsP_REC_cut2",";p_{trk} [GeV]; E_{EEMC} [GeV]",100,0,20,100,0,20);
    TH2D* h_EvsP_MC_cut2 = new TH2D("h_EvsP_MC_cut2",";p_{MC} [GeV]; E_{MC} [GeV]",100,0,20,100,0,20);
    TH1D* h_Epz_REC_cut2 = new TH1D("h_Epz_REC_cut2", ";(E_{EEMC} - p_{z,EEMC}) [GeV]",100,0,20);
    TH1D* h_Epz_REC_trk_cut2 = new TH1D("h_Epz_REC_trk_cut2", ";(E_{trk} - p_{z,trk}) [GeV]",100,0,20);
    TH1D* h_Epz_MC_cut2 = new TH1D("h_Epz_MC_cut2", ";(E_{MC} - p_{z,MC}) [GeV]",100,0,20);
    TH2D* h_Epz_response_cut2 = new TH2D("h_Epz_response_cut2", ";(E_{EEMC} - p_{z,EEMC}) [GeV];(E_{MC}-p_{z,MC}) [GeV]",100,0,20,100,0,20);
    TH1D* h_EcalOverPtrk_cut2 = new TH1D("h_EcalOverPtrk_cut2",";E_{EEMC}/|p|_{trk}",100,0,20);
    TH1D* h_EtrkOverPcal_cut2 = new TH1D("h_EtrkOverPcal_cut2",";E_{trk}/|p|_{EEMC}",100,0,20);
    TH2D* h_EoverP_response_cut2 = new TH2D("h_EoverP_response_cut2",";E_{EEMC}/|p|_{trk};E_{MC}/|p|_{MC}",100,0,20,100,0,20);
    TH1D* h_theta_diff_cut2 = new TH1D("h_theta_diff_cut2",";#theta_{EEMC}-#theta_{MC}",100,-3.14,3.14);
    TH1D* h_phi_diff_cut2 = new TH1D("h_phi_diff_cut2",";#phi_{EEMC}-#phi_{MC}",100,-3.14,3.14);
    TH1D* h_eta_diff_cut2 = new TH1D("h_eta_diff_cut2",";#eta_{EEMC}-#eta_{MC}",100,-3.14,3.14);
    // t
    TH1D* h_t_REC_EEMC_cut2 = new TH1D("h_t_REC_EEMC_cut2",";|t|_{EEMC} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_trk_cut2 = new TH1D("h_t_REC_trk_cut2",";|t|_{trk} [GeV/c]^{2}; counts",100,0,0.2);
    TH2D* h_t_res_EEMC_cut2 = new TH2D("h_t_res_EEMC_cut2",";|t|_{MC} [GeV/c]^{2}; (|t|_{MC}-|t|_{EEMC})/|t|_{MC}",100,0,0.2,1000,-10,10);
    TH2D* h_t_res_trk_cut2 = new TH2D("h_t_res_trk_cut2",";|t|_{MC} [GeV/c]^{2}; (|t|_{MC}-|t|_{trk})/|t|_{MC}",100,0,0.2,1000,-10,10);
    TH2D* h_t_response_EEMC_cut2 = new TH2D("h_t_response_EEMC_cut2","; |t|_{EEMC} [GeV/c]^{2};|t|_{MC} [GeV/c]^{2}",100,0,0.2,1000,0,0.2);
	TH2D* h_t_response_trk_cut2 = new TH2D("h_t_response_trk_cut2","; |t|_{trk} [GeV/c]^{2};|t|_{MC} [GeV/c]^{2}",100,0,0.2,1000,0,0.2);
	TH2D* h_t_response_trk_EEMC_cut2 = new TH2D("h_t_response_trk_EEMC_cut2","; |t|_{trk} [GeV/c]^{2};|t|_{EEMC} [GeV/c]^{2}",100,0,0.2,1000,0,0.2);
    // t distribution with resolution only
    TH1D* h_t_REC_wRES_cut2 = new TH1D("h_t_REC_wRES_cut2",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    // t distribution with resolution and angle cut (EEMC RECO)
    TH1D* h_t_REC_wRES_cut_pi2_cut2 = new TH1D("h_t_REC_wRES_cut_pi2_cut2",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wRES_cut_pi3_cut2 = new TH1D("h_t_REC_wRES_cut_pi3_cut2",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wRES_cut_pi4_cut2 = new TH1D("h_t_REC_wRES_cut_pi4_cut2",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wRES_cut_pi6_cut2 = new TH1D("h_t_REC_wRES_cut_pi6_cut2",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wRES_cut_pi9_cut2 = new TH1D("h_t_REC_wRES_cut_pi9_cut2",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wRES_cut_pi12_cut2 = new TH1D("h_t_REC_wRES_cut_pi12_cut2",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wRES_cut_pi16_cut2 = new TH1D("h_t_REC_wRES_cut_pi16_cut2",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wRES_cut_pi20_cut2 = new TH1D("h_t_REC_wRES_cut_pi20_cut2",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wRES_cut_pi24_cut2 = new TH1D("h_t_REC_wRES_cut_pi24_cut2",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    // t distribution 2d
    TH2D* h_t_REC_2d_cut2 = new TH2D("h_t_REC_2d_cut2",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    // t distribution 2d with resolution only 
    TH2D* h_t_REC_2d_wRES_cut2 = new TH2D("h_t_REC_2d_wRES_cut2",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    // t distribution 2d with resolution and angle cut (EEMC RECO)
    TH2D* h_t_REC_2d_wRES_cut_pi2_cut2 = new TH2D("h_t_REC_2d_wRES_cut_pi2_cut2",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wRES_cut_pi3_cut2 = new TH2D("h_t_REC_2d_wRES_cut_pi3_cut2",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wRES_cut_pi4_cut2 = new TH2D("h_t_REC_2d_wRES_cut_pi4_cut2",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wRES_cut_pi6_cut2 = new TH2D("h_t_REC_2d_wRES_cut_pi6_cut2",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wRES_cut_pi9_cut2 = new TH2D("h_t_REC_2d_wRES_cut_pi9_cut2",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wRES_cut_pi12_cut2 = new TH2D("h_t_REC_2d_wRES_cut_pi12_cut2",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wRES_cut_pi16_cut2 = new TH2D("h_t_REC_2d_wRES_cut_pi16_cut2",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wRES_cut_pi20_cut2 = new TH2D("h_t_REC_2d_wRES_cut_pi20_cut2",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wRES_cut_pi24_cut2 = new TH2D("h_t_REC_2d_wRES_cut_pi24_cut2",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);

    // Uncomment this block to run single files at a time
    // and comment out "chain->GetEntries();""
    //tree_reader.SetEntriesRange(0, tree->GetEntries()); 
    chain->GetEntries();
    while (tree_reader.Next()) 
    {
         /*---------------
            Event Loop
        ----------------*/
        
        event->clear(); 
        
    	// Beam particles
    	TLorentzVector ebeam(0,0,0,0);
    	TLorentzVector pbeam(0,0,0,0);
        TLorentzVector abeam(0,0,0,0);
        TLorentzVector vmMC(0,0,0,0);
    	TLorentzVector kplusMC(0,0,0,0);
    	TLorentzVector kminusMC(0,0,0,0);

        // incoherent check
        /*if(event->clusters_zdc.size()>0) continue;
		if(event->hit_rp.size()>0) continue;
		if(event->hit_omd.size()>0) continue;*/

    	//MC level
    	TLorentzVector scatMC(0,0,0,0);
    	int mc_elect_index=-1;
    	double maxPt=-99.;
        //int incoherent=0;
        
        // loop over all MC particles in the event
    	for(int imc=0;imc<mc_px_array.GetSize();imc++)
        {
            TVector3 mctrk(static_cast<double>(mc_px_array[imc]), 
                static_cast<double>(mc_py_array[imc]), 
                static_cast<double>(mc_pz_array[imc]));
    		if(mc_genStatus_array[imc]==4) // 4 is Sartre.
            {
                /*--------------------------------------
                    -genStatus = 4 is Sartre
                    -genStatus = 1 is stable particle
                    -PDG 11 = electron
                    -PDG 2212 = proton
                    -PDG 321 = kaon
                    -All other PDG set as ion
                ---------------------------------------*/
    			if(mc_pdg_array[imc]==11) ebeam.SetVectM(mctrk, MASS_ELECTRON);
				if(mc_pdg_array[imc]==2212) pbeam.SetVectM(mctrk, MASS_PROTON);
                else
                {
                    double MASS_A = mc_mass_array[imc];
                    abeam.SetVectM(mctrk,MASS_A);
                    cout<<"genStatus=4 candidate: PDG="<<mc_pdg_array[imc]<<" Mass="<< MASS_A<<endl;
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
             //if(mctrk.Eta()>3.5&&(mc_mass_array[imc]>0.9383||mc_mass_array[imc]<0.938))incoherent++;
    	}

         /*---------------------------
            Get incoherent events
        ----------------------------*/
       
		//EEMC
       /* for(int iclus=0;iclus<em_energy_array.GetSize();iclus++)
        {
      		Cluster_EEMC cluster;
			cluster.energy=em_energy_array[iclus];
			cluster.x=em_x_array[iclus];
			cluster.y=em_y_array[iclus];
    	    event->clusters_eemc.push_back(cluster);
    	}

    	//ZDC
        for(int iclus=0;iclus<zdc_energy_array.GetSize();iclus++)
        {
      		Cluster_ZDC cluster;
			cluster.energy=zdc_energy_array[iclus];
			cluster.x=zdc_x_array[iclus];
			cluster.y=zdc_y_array[iclus];
			cluster.z=zdc_z_array[iclus];
    	    event->clusters_zdc.push_back(cluster);
    	}

    	//RP
    	for(int ihit=0;ihit<rp_x_array.GetSize();ihit++)
        {
      		Hit_RP hit;
			hit.x=rp_x_array[ihit];
			hit.y=rp_y_array[ihit];
			hit.z=rp_z_array[ihit];
    	    event->hit_rp.push_back(hit);
    	}

    	//OMD
    	for(int ihit=0;ihit<omd_x_array.GetSize();ihit++)
        {
      		Hit_OMD hit;
			hit.x=omd_x_array[ihit];
			hit.y=omd_y_array[ihit];
			hit.z=omd_z_array[ihit];
    	    event->hit_omd.push_back(hit);
    	}            

        if(incoherent) continue;*/
    
        h_phi_MC->Fill(scatMC.Phi());
        h_theta_MC->Fill(scatMC.Theta());
        h_eta_MC->Fill(scatMC.Eta());
        h_e_pt_MC->Fill(scatMC.Pt());
		h_e_pz_MC->Fill(scatMC.Pz());
		h_e_p_MC->Fill(scatMC.P());
		h_energy_MC->Fill(scatMC.E());
        h_EoverP_MC->Fill(scatMC.E()/scatMC.P());

		// checks
        cout<<"A energy: "<<abeam.E()<<" p Energy: "<<pbeam.E()<<" e Energy: "<<ebeam.E()<<endl;
        cout<<"Electron Beam: "<<" px: "<<ebeam.Px()<<" py: "<<ebeam.Py()<<" pz: "<<ebeam.Pz()<< " E: "<<ebeam.E()<< endl;
		cout<<"p Beam: "<<" px: "<<pbeam.Px()<<" py: "<<pbeam.Py()<<" pz: "<<pbeam.Pz()<<" E: "<<pbeam.E()<< endl;
        cout <<"A beam: "<<" px: "<< abeam.Px()<<" py: "<< abeam.Py()<<" pz: "<< abeam.Pz()<< " E: "<< abeam.E()<< endl;
		cout <<"Scattered Electron: "<<" px: "<<scatMC.Px()<<" py: "<<scatMC.Py()<<" pz: "<< scatMC.Pz()<<" E: "<< scatMC.E()<<endl;
    
    	// protection.
    	if(ebeam.E()==abeam.E() && ebeam.E()==0) 
        {
    		cout << "problem with MC incoming beams" << endl;
    		continue;
    	}
        
    	TLorentzVector qbeam=ebeam-scatMC; // p_e - p_e'
    	double Q2=-(qbeam).Mag2();  
    	double pq=abeam.Dot(qbeam); 
    	double y=pq/abeam.Dot(ebeam);

        // MC level phase space cut
    	if(Q2<1.||Q2>10.) continue;
    	if(y<0.01||y>0.85) continue;
    	h_Q2_e->Fill(Q2);
    	h_y_e->Fill(y);

        // VM 
		vmMC=kplusMC+kminusMC;
    	if(vmMC.E()==0) continue;
        h_VM_mass_MC->Fill(vmMC.M());
    	h_VM_pt_MC->Fill(vmMC.Pt());
		h_VM_pz_MC->Fill(vmMC.Pz());
		h_VM_p_MC->Fill(vmMC.P());
        h_VM_Epz_MC->Fill(vmMC.E()-vmMC.Pz());
        

		// t dist
        double t_MC = 0;
    	if(vmMC.E()!=0 && fabs(vmMC.Rapidity())<3.5)
    	{
    		double method_E = -(qbeam-vmMC).Mag2();  // t = -(p_e - p_e' - p_v)^2
    		t_MC=method_E;
    		h_t_MC->Fill( method_E );

    		// new method
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

            // no resolution added so should still be t_MC 
            double t_total = tx_true+ty+tz;

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
        // find highest energy cluster in negative endcap
        // to find scattered electron
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
    	if( radius>550. ) continue; 
    	double clusEnergy=1.044*maxHitEnergy; // 4.4% energy calibration.
    	h_energy_REC_EEMC->Fill(clusEnergy);

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
        double res= (scatMC.E()-clusEnergy)/scatMC.E();
        h_energy_res_EEMC->Fill(scatMC.E(), res);
        h_energy_response_EEMC->Fill(clusEnergy,scatMC.E());
        
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
        TLorentzVector vmREC(0,0,0,0);

        double maxP=-1.;
        // track loop
        for(int itrk=0;itrk<reco_pz_array.GetSize();itrk++)
        {
        	TVector3 trk(reco_px_array[itrk],reco_py_array[itrk],reco_pz_array[itrk]);
        	if(rec_elect_index!=-1 && itrk==rec_elect_index) 
            {
                // if stable and matches electron index
        		scatMCmatchREC.SetVectM(trk,MASS_ELECTRON); 
        	}
        	if(trk.Mag()>maxP)
            {
        		maxP=trk.Mag();
        		scatREC.SetVectM(trk,MASS_ELECTRON);
        		// use emcal energy to define 4 vector
        		double p = sqrt(clusEnergy*clusEnergy- MASS_ELECTRON*MASS_ELECTRON );
                double pt = TMath::Sin(scatREC.Theta())*p;
        		scatClusEREC.SetPtEtaPhiM(pt,scatREC.Eta(),scatREC.Phi(),MASS_ELECTRON);
        	}
        }
        
        h_eta_REC_EEMC->Fill(scatClusEREC.Eta());
        h_phi_REC_EEMC->Fill(scatClusEREC.Phi());
        h_theta_REC_EEMC->Fill(scatREC.Theta());
        h_theta_response_EEMC->Fill(scatMC.Theta(),scatREC.Theta());
        h_theta_diff->Fill(scatREC.Theta()-scatMC.Theta());
        h_phi_diff->Fill(scatClusEREC.Phi()-scatMC.Phi());
        h_eta_diff->Fill(scatClusEREC.Eta()-scatMC.Eta());
            
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
    	double res_trk= (scatMC.E()-scatMCmatchREC.E())/scatMC.E();
        h_energy_res_trk->Fill(scatMC.E(), res_trk);
    	h_energy_REC_trk->Fill(scatMCmatchREC.E());

        //int incoherent_rec=0;
        // loop over track again;
        for(int itrk=0;itrk<reco_pz_array.GetSize();itrk++)
        {
        	TVector3 trk(reco_px_array[itrk],reco_py_array[itrk],reco_pz_array[itrk]);
        	particle.SetVectM(trk,MASS_PION); // assume pions;
            /*------------------------------------------------------
                    -exclude e' from hadronic final state and kaons
                    -select phi->kk daughters;
            -------------------------------------------------------*/
        	if(itrk!=rec_elect_index) 
            {
            	hfs += particle;
            	h_eta_REC_trk->Fill(trk.Eta());
            	if(fabs(trk.Eta())<3.0)
                {
            		if(reco_charge_array[itrk]>0) kplusREC.SetVectM(trk,MASS_KAON);
            		if(reco_charge_array[itrk]<0) kminusREC.SetVectM(trk,MASS_KAON);
            	}
        	}
            //incoherent
        	//if(trk.Eta()>3.5) incoherent_rec++;
        }

        //selection
    	//if(incoherent_rec) continue;
        
        // 4vector of VM;
        if(kplusREC.E()!=0. && kminusREC.E()!=0.)
        {
        	vmREC=kplusREC+kminusREC;
        }
        h_VM_mass_REC->Fill(vmREC.M());
        h_VM_pt_REC->Fill(vmREC.Pt());
    	h_VM_pz_REC->Fill(vmREC.Pz());
    	h_VM_p_REC->Fill(vmREC.P());
        h_VM_Epz_REC->Fill(vmREC.E()-vmREC.Pz());
        h_VM_Epz_response->Fill(vmREC.E()-vmREC.Pz(),vmMC.E()-vmMC.Pz());
        h_VM_pt_response->Fill(vmREC.Pt(),vmMC.Pt());

        // e' p from emcal
        h_e_pt_REC_EEMC->Fill(scatClusEREC.Pt());
    	h_e_pz_REC_EEMC->Fill(scatClusEREC.Pz());
    	h_e_p_REC_EEMC->Fill(scatClusEREC.P());
        h_e_pt_res->Fill(scatMC.Pt(),(scatClusEREC.Pt()-scatMC.Pt())/scatMC.Pt());
        h_e_pz_res->Fill(scatMC.Pz(),(scatClusEREC.Pz()-scatMC.Pz())/scatMC.Pz());
        h_e_p_res->Fill(scatMC.P(),(scatClusEREC.P()-scatMC.P())/scatMC.P());
        h_e_pz_response->Fill(scatClusEREC.Pz(),scatMC.Pz());
        
    	// e' p from track
        h_e_pt_REC_trk->Fill(scatMCmatchREC.Pt());
    	h_e_pz_REC_trk->Fill(scatMCmatchREC.Pz());
    	h_e_p_REC_trk->Fill(scatMCmatchREC.P());

        // E_cal vs p_trk
        h_EvsP_REC->Fill(scatMCmatchREC.P(),scatREC.E());
        h_EvsP_MC->Fill(scatMC.P(),scatREC.E());

        // cluster-base DIS kine;
        TLorentzVector qbeamREC=ebeam-scatClusEREC; // emcal e' reco
        TLorentzVector qbeamREC_trk=ebeam-scatMCmatchREC; // track e' reco
        double Q2REC=-(qbeamREC).Mag2();  
        double Q2REC_trk = -(qbeamREC_trk).Mag2();
        double pqREC=abeam.Dot(qbeamREC); 
        double yREC=pqREC/abeam.Dot(ebeam); 
        double pqREC_trk=abeam.Dot(qbeamREC_trk); 
        double yREC_trk=pqREC_trk/abeam.Dot(ebeam); 
        h_Q2REC_e_EEMC->Fill(Q2REC);
        h_Q2REC_e_trk->Fill(Q2REC_trk);
        h_yREC_e_EEMC->Fill(yREC);
        h_yREC_e_trk->Fill(yREC_trk);

        // Q2 emcal reco resolution
    	res = (Q2 - Q2REC)/Q2;
    	h_Q2_res->Fill(Q2, res);
        h_Q2_response->Fill(Q2REC,Q2);
        h_dQ2overQ2_REC->Fill(res);
        
        // y emcal reco resolution
        res = (y - yREC)/y;
        h_y_res->Fill(y,res);
        h_y_response->Fill(yREC,y);
        h_dyOvery_REC->Fill(res);

    	// track-base Epz scat' e
        double EpzREC_trk = (scatMCmatchREC+hfs).E() - (scatMCmatchREC+hfs).Pz();
        h_Epz_REC_trk->Fill(EpzREC_trk);
        double EpzREC = (scatClusEREC+hfs).E() - (scatMCmatchREC+hfs).Pz();
        h_Epz_REC->Fill(EpzREC);
        double EpzMC = scatMC.E()-scatMC.Pz();
        h_Epz_MC->Fill(EpzMC);
        h_Epz_response->Fill(EpzREC,EpzMC);
    
    	// E over p  
        double EoverP = scatClusEREC.E() / scatMCmatchREC.P();
        h_EcalOverPtrk->Fill(scatClusEREC.E()/scatMCmatchREC.P());
        h_EtrkOverPcal->Fill(scatClusEREC.E()/scatMCmatchREC.P());
        h_EoverP_response->Fill(scatClusEREC.E()/scatMCmatchREC.P(),scatMC.E()/scatMC.P());

        // Event selection:
        if( EpzREC>27||EpzREC<40 ) 
        {
            h_Epz_afterCut->Fill(EpzREC);
            if( EoverP>0.8||EoverP<1.18 )
            {
                h_EoverP_afterCut->Fill(EoverP);
                if(Q2REC>1.||Q2REC<10.) 
                {
                    h_Q2_afterCut->Fill(Q2REC);
                    if(yREC>0.01||yREC<0.85) 
                    {
                        h_yREC_afterCut->Fill(yREC);
                        // VM rec
                        if(vmREC.E()!=0) 
                            {
                        // select phi mass and rapidity window 
                        if( fabs(vmREC.M()-1.02)<0.02&& fabs(vmREC.Rapidity())<3.5)
                        {
                            // 2 versions: track and energy cluster:
                        	double t_trk_REC = giveme_t_method_L(ebeam,scatMCmatchREC,abeam,vmREC);  // method L (track e')
                            double t_REC = giveme_t_method_L(ebeam,scatClusEREC,abeam,vmREC); // method L (EEMC e')
                            h_t_REC_trk->Fill(t_trk_REC);
                            h_t_REC_EEMC->Fill(t_REC);
                            h_t_response_EEMC->Fill(t_REC,t_MC);
                            h_t_response_trk->Fill(t_trk_REC,t_MC);
                            h_t_response_trk_EEMC->Fill(t_trk_REC,t_REC);

                        	// new method
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
                        	h_t_REC_2d_wRES->Fill(qx_rec,qy_rec);
                        	h_t_REC_wRES->Fill(t_total_rec);

                            //t track resolution 
                        	res = (t_MC-t_trk_REC)/t_MC;
                        	h_t_res_trk->Fill(t_MC, res);
	
                            //t EEMC resolution;
                        	res = (t_MC-t_REC)/t_MC;
                        	h_t_res_EEMC->Fill(t_MC, res);

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

        /*---------------------------------------------------------
            -Reject position that are outside of the threshold
            -Repeat the reco loop same as above
        ----------------------------------------------------------*/
        double position_threshold = 20; // mm
        if(fabs(x_clus_trk_diff)<=position_threshold && fabs(y_clus_trk_diff)<=position_threshold)
        {
            // energy resolution
            double res= (scatMC.E()-clusEnergy)/scatMC.E();
            h_energy_res_EEMC_cut->Fill(scatMC.E(), res);
            h_energy_response_EEMC_cut->Fill(clusEnergy,scatMC.E());
            h_energy_REC_EEMC_cut->Fill(clusEnergy);
        
            // association of rec level scat' e
            int rec_elect_index=-1;
            for(int i=0;i<sim_id.GetSize();i++)
            {
            	if(sim_id[i]==mc_elect_index) 
                {
            		rec_elect_index = rec_id[i]; 
            	}
            }
    
            TLorentzVector scatMCmatchREC(0,0,0,0);
            TLorentzVector scatREC(0,0,0,0);
            TLorentzVector scatClusEREC(0,0,0,0);
            TLorentzVector hfs(0,0,0,0);
            TLorentzVector particle(0,0,0,0);
            TLorentzVector kplusREC(0,0,0,0);
            TLorentzVector kminusREC(0,0,0,0);
            TLorentzVector vmREC(0,0,0,0);

            double maxP=-1.;
            // track loop
            for(int itrk=0;itrk<reco_pz_array.GetSize();itrk++)
            {
            	TVector3 trk(reco_px_array[itrk],reco_py_array[itrk],reco_pz_array[itrk]);
            	if(rec_elect_index!=-1 && itrk==rec_elect_index)
                {
            		scatMCmatchREC.SetVectM(trk,MASS_ELECTRON);  
            	}
            	if(trk.Mag()>maxP)
                {
            		maxP=trk.Mag();
            		scatREC.SetVectM(trk,MASS_ELECTRON);
            		double p = sqrt(clusEnergy*clusEnergy- MASS_ELECTRON*MASS_ELECTRON );
                    double pt = TMath::Sin(scatREC.Theta())*p;
            		scatClusEREC.SetPtEtaPhiM(pt,scatREC.Eta(),scatREC.Phi(),MASS_ELECTRON);
            	}
            }

            h_eta_REC_EEMC_cut->Fill(scatClusEREC.Eta());
            h_phi_REC_EEMC_cut->Fill(scatClusEREC.Phi());
            h_theta_REC_EEMC_cut->Fill(scatClusEREC.Theta());
            h_theta_response_EEMC_cut->Fill(scatClusEREC.Theta(),scatMC.Theta());
            h_theta_diff_cut->Fill(scatClusEREC.Theta()-scatMC.Theta());
            h_phi_diff_cut->Fill(scatClusEREC.Phi()-scatMC.Phi());
            h_eta_diff_cut->Fill(scatClusEREC.Eta()-scatMC.Eta());
            
            // default clustering position
        	h_emClus_position_x_REC_cut->Fill(xpos);
        	h_emClus_position_y_REC_cut->Fill(ypos);
            h_XvsY_clus_cut->Fill(xpos,ypos);
        
        	// self clustering position
        	h_emHits_position_x_REC_cut->Fill(xClus);
        	h_emHits_position_y_REC_cut->Fill(yClus);
            h_XvsY_hits_cut->Fill(xClus,yClus);

            // track positions
            h_XvsY_trk_cut->Fill(xtrk,ytrk);
            h_trk_position_x_REC_cut->Fill(xtrk);
            h_trk_position_y_REC_cut->Fill(ytrk);

            h_Xclus_minus_Xtrk_cut->Fill(x_clus_trk_diff);
            h_Yclus_minus_Ytrk_cut->Fill(y_clus_trk_diff);
            h_Xhits_minus_Xtrk_cut->Fill(x_hit_trk_diff);
            h_Yhits_minus_Ytrk_cut->Fill(y_hit_trk_diff);
        
        	// track-base e' energy
        	double res_trk= (scatMC.E()-scatMCmatchREC.E())/scatMC.E();
            h_energy_res_trk_cut->Fill(scatMC.E(), res_trk);
        	h_energy_REC_trk_cut->Fill(scatMCmatchREC.E());

            // loop over track again;
            for(int itrk=0;itrk<reco_pz_array.GetSize();itrk++)
            {
            	TVector3 trk(reco_px_array[itrk],reco_py_array[itrk],reco_pz_array[itrk]);
            	particle.SetVectM(trk,MASS_PION);
            	if(itrk!=rec_elect_index) 
                {
                  	hfs += particle; 
                	h_eta_REC_trk_cut->Fill(trk.Eta());
                	if(fabs(trk.Eta())<3.0)
                    {
                		if(reco_charge_array[itrk]>0) kplusREC.SetVectM(trk,MASS_KAON);
                		if(reco_charge_array[itrk]<0) kminusREC.SetVectM(trk,MASS_KAON);
                	}
            	}
            }
        
            // 4vector of VM;
            if(kplusREC.E()!=0. && kminusREC.E()!=0.)
             {
            	vmREC=kplusREC+kminusREC;
            }
            h_VM_mass_REC_cut->Fill(vmREC.M());
            h_VM_pt_REC_cut->Fill(vmREC.Pt());
        	h_VM_pz_REC_cut->Fill(vmREC.Pz());
        	h_VM_p_REC_cut->Fill(vmREC.P());
            h_VM_Epz_MC_cut->Fill(vmMC.E()-vmMC.Pz());
            h_VM_Epz_REC_cut->Fill(vmREC.E()-vmREC.Pz());
            h_VM_Epz_response_cut->Fill(vmREC.E()-vmREC.Pz(),vmMC.E()-vmMC.Pz());
            h_VM_pt_response_cut->Fill(vmREC.Pt(),vmMC.Pt());

            // e' p from emcal
            h_e_pt_REC_EEMC_cut->Fill(scatClusEREC.Pt());
        	h_e_pz_REC_EEMC_cut->Fill(scatClusEREC.Pz());
        	h_e_p_REC_EEMC_cut->Fill(scatClusEREC.P());
            h_e_pt_res_cut->Fill(scatMC.Pt(),(scatClusEREC.Pt()-scatMC.Pt())/scatMC.Pt());
            h_e_pz_res_cut->Fill(scatMC.Pz(),(scatClusEREC.Pz()-scatMC.Pz())/scatMC.Pz());
            h_e_p_res_cut->Fill(scatMC.P(),(scatClusEREC.P()-scatMC.P())/scatMC.P());
            h_e_pz_response_cut->Fill(scatClusEREC.Pz(),scatMC.Pz());
        
        	// e' p from track
            h_e_pt_REC_trk_cut->Fill(scatMCmatchREC.Pt());
        	h_e_pz_REC_trk_cut->Fill(scatMCmatchREC.Pz());
        	h_e_p_REC_trk_cut->Fill(scatMCmatchREC.P());

            // E_cal vs p_trk
            h_EvsP_REC_cut->Fill(scatMCmatchREC.P(),scatREC.E());
            h_EvsP_MC_cut->Fill(scatMC.P(),scatREC.E());

            // cluster-base DIS kine;
            TLorentzVector qbeamREC=ebeam-scatClusEREC; // emcal e' reco
            TLorentzVector qbeamREC_trk=ebeam-scatMCmatchREC; // track e' reco
            double Q2REC=-(qbeamREC).Mag2();  
            double Q2REC_trk = -(qbeamREC_trk).Mag2();
            double pqREC=abeam.Dot(qbeamREC); 
            double yREC=pqREC/abeam.Dot(ebeam); 
            double pqREC_trk=abeam.Dot(qbeamREC_trk); 
            double yREC_trk=pqREC_trk/abeam.Dot(ebeam); 
            h_Q2REC_e_EEMC_cut->Fill(Q2REC);
            h_Q2REC_e_trk_cut->Fill(Q2REC_trk);
            h_yREC_e_EEMC_cut->Fill(yREC);
            h_yREC_e_trk_cut->Fill(yREC_trk);

            // Q2 emcal reco resolution
        	res = (Q2 - Q2REC)/Q2;
        	h_Q2_res_cut->Fill(Q2, res);
            h_Q2_response_cut->Fill(Q2REC,Q2);
            h_dQ2overQ2_REC->Fill(res);

            // y emcal reco resolution
            res = (y - yREC)/y;
            h_y_res_cut->Fill(y,res);
            h_y_response_cut->Fill(yREC,y);
            h_dyOvery_REC_cut->Fill(res);
    
        	// track-base Epz scat' e
            double EpzREC_trk = (scatMCmatchREC+hfs).E() - (scatMCmatchREC+hfs).Pz();
            h_Epz_REC_trk_cut->Fill(EpzREC_trk);
            double EpzREC = (scatClusEREC+hfs).E() - (scatMCmatchREC+hfs).Pz();
            h_Epz_REC_cut->Fill(EpzREC);
            h_Epz_response_cut->Fill(EpzREC,scatMC.E()-scatMC.Pz());
            h_Epz_MC_cut->Fill(EpzMC);
    
        	// E over p  
            double EoverP = scatClusEREC.E() / scatMCmatchREC.P();
            h_EcalOverPtrk_cut->Fill(scatClusEREC.E()/scatMCmatchREC.P());
            h_EtrkOverPcal_cut->Fill(scatMCmatchREC.E()/scatClusEREC.P());
            h_EoverP_response_cut->Fill(scatClusEREC.E()/scatMCmatchREC.P(),scatMC.E()/scatMC.P());

            // Event selection:
            if( EpzREC>27||EpzREC<40 )
            {	
                if( EoverP>0.8||EoverP<1.18 )
                {
                if(Q2REC>1.||Q2REC<10.)
                {
                    if(yREC>0.01||yREC<0.85)
                    {
                        // VM rec
                        if(vmREC.E()!=0){
                    	// select phi mass and rapidity window 
                    	if( fabs(vmREC.M()-1.02)<0.02&& fabs(vmREC.Rapidity())<3.5)
                    	{
                        	// 2 versions: track and energy cluster:
                    		double t_trk_REC = giveme_t_method_L(ebeam,scatMCmatchREC,abeam,vmREC);  // method L (track e')
                        	double t_REC = giveme_t_method_L(ebeam,scatClusEREC,abeam,vmREC); // method L (EEMC e')
                        	h_t_REC_trk_cut->Fill(t_trk_REC);
                        	h_t_REC_EEMC_cut->Fill(t_REC);
                            h_t_response_EEMC_cut->Fill(t_REC,t_MC);
                            h_t_response_trk_cut->Fill(t_trk_REC,t_MC);
                            h_t_response_trk_EEMC_cut->Fill(t_trk_REC,t_REC);

                    		// new method
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
        
                        	// t track resolution 
                    		res = (t_MC-t_trk_REC)/t_MC;
                    		h_t_res_trk_cut->Fill(t_MC, res);
	
                        	// t EEMC resolution;
                    		res = (t_MC-t_REC)/t_MC;
                    		h_t_res_EEMC_cut->Fill(t_MC, res);

                    		// apply cut with // subtracted out of t
                    		double theta_rec = atan(fabs(qx_rec)/fabs(qy_rec));
                    		if(fabs(theta_rec)<PI/2)
                    		{
                    			h_t_REC_wRES_cut_pi2_cut->Fill(t_total_rec);
                    			h_t_REC_2d_wRES_cut_pi2_cut->Fill(qx_rec,qy_rec);
                    		}
                    		if(fabs(theta_rec)<PI/3)
                    		{
                    			h_t_REC_wRES_cut_pi3_cut->Fill(t_total_rec);
                    			h_t_REC_2d_wRES_cut_pi3_cut->Fill(qx_rec,qy_rec);
                    		}
                    		if(fabs(theta_rec)<PI/4)
                    		{
                    			h_t_REC_wRES_cut_pi4_cut->Fill(t_total_rec);
                    			h_t_REC_2d_wRES_cut_pi4_cut->Fill(qx_rec,qy_rec);
                    		}
                    		if(fabs(theta_rec)<PI/6)
                    		{
                    			h_t_REC_wRES_cut_pi6_cut->Fill(t_total_rec);
                    			h_t_REC_2d_wRES_cut_pi6_cut->Fill(qx_rec,qy_rec);
                    		}
                    		if(fabs(theta_rec)<PI/9)
                    		{
                    			h_t_REC_wRES_cut_pi9_cut->Fill(t_total_rec);
                    			h_t_REC_2d_wRES_cut_pi9_cut->Fill(qx_rec,qy_rec);
                    		}
                    		if(fabs(theta_rec)<PI/12)
                    		{
                    			h_t_REC_wRES_cut_pi12_cut->Fill(t_total_rec);
                    			h_t_REC_2d_wRES_cut_pi12_cut->Fill(qx_rec,qy_rec);
                    		}
                    		if(fabs(theta_rec)<PI/16)
                    		{
                    			h_t_REC_wRES_cut_pi16_cut->Fill(t_total_rec);
                    			h_t_REC_2d_wRES_cut_pi16_cut->Fill(qx_rec,qy_rec);
                    		}
                            if(fabs(theta_rec)<PI/20)
                    		{
                    			h_t_REC_wRES_cut_pi20_cut->Fill(t_total_rec);
                    			h_t_REC_2d_wRES_cut_pi20_cut->Fill(qx_rec,qy_rec);
                    		}
                            if(fabs(theta_rec)<PI/24)
                    		{
                    			h_t_REC_wRES_cut_pi24_cut->Fill(t_total_rec);
                    			h_t_REC_2d_wRES_cut_pi24_cut->Fill(qx_rec,qy_rec);
                    		}
                        }
                        }
                    }
                }
                }
            }
        }
           
        /*---------------------------------------------------------
            -Reject position that are inside of the threshold
            -Repeat the reco loop same as above
            -Just for completeness and further analysis
        ----------------------------------------------------------*/
        else
        {
            // energy resolution
            double res= (scatMC.E()-clusEnergy)/scatMC.E();
            h_energy_res_EEMC_cut2->Fill(scatMC.E(), res);
            h_energy_response_EEMC_cut2->Fill(clusEnergy,scatMC.E());
            h_energy_REC_EEMC_cut2->Fill(clusEnergy);
        
            // association of rec level scat' e
            int rec_elect_index=-1;
            for(int i=0;i<sim_id.GetSize();i++)
            {
            	if(sim_id[i]==mc_elect_index) 
                {
            		rec_elect_index = rec_id[i]; 
            	}
            }
    
            TLorentzVector scatMCmatchREC(0,0,0,0);
            TLorentzVector scatREC(0,0,0,0);
            TLorentzVector scatClusEREC(0,0,0,0);
            TLorentzVector hfs(0,0,0,0);
            TLorentzVector particle(0,0,0,0);
            TLorentzVector kplusREC(0,0,0,0);
            TLorentzVector kminusREC(0,0,0,0);
            TLorentzVector vmREC(0,0,0,0);

            double maxP=-1.;
            // track loop
            for(int itrk=0;itrk<reco_pz_array.GetSize();itrk++)
            {
            	TVector3 trk(reco_px_array[itrk],reco_py_array[itrk],reco_pz_array[itrk]);
            	if(rec_elect_index!=-1 && itrk==rec_elect_index) 
                {
            		scatMCmatchREC.SetVectM(trk,MASS_ELECTRON);  
            	}
            	if(trk.Mag()>maxP)
                {
            		maxP=trk.Mag();
            		scatREC.SetVectM(trk,MASS_ELECTRON);
            		double p = sqrt(clusEnergy*clusEnergy- MASS_ELECTRON*MASS_ELECTRON );
                    double pt = TMath::Sin(scatREC.Theta())*p;
            		scatClusEREC.SetPtEtaPhiM(pt,scatREC.Eta(),scatREC.Phi(),MASS_ELECTRON);
            	}
            }

            h_eta_REC_EEMC_cut2->Fill(scatClusEREC.Eta());
            h_phi_REC_EEMC_cut2->Fill(scatClusEREC.Phi());
            h_theta_REC_EEMC_cut2->Fill(scatClusEREC.Theta());
            h_theta_response_EEMC_cut2->Fill(scatClusEREC.Theta(),scatMC.Theta());
            h_theta_diff_cut2->Fill(scatClusEREC.Theta()-scatMC.Theta());
            h_phi_diff_cut2->Fill(scatClusEREC.Phi()-scatMC.Phi());
            h_eta_diff_cut2->Fill(scatClusEREC.Eta()-scatMC.Eta());
            
            // default clustering position
        	h_emClus_position_x_REC_cut2->Fill(xpos);
        	h_emClus_position_y_REC_cut2->Fill(ypos);
            h_XvsY_clus_cut2->Fill(xpos,ypos);
        
        	// self clustering position
        	h_emHits_position_x_REC_cut2->Fill(xClus);
        	h_emHits_position_y_REC_cut2->Fill(yClus);
            h_XvsY_hits_cut2->Fill(xClus,yClus);

            // track positions
            h_XvsY_trk_cut2->Fill(xtrk,ytrk);
            h_trk_position_x_REC_cut2->Fill(xtrk);
            h_trk_position_y_REC_cut2->Fill(ytrk);

            h_Xclus_minus_Xtrk_cut2->Fill(x_clus_trk_diff);
            h_Yclus_minus_Ytrk_cut2->Fill(y_clus_trk_diff);
            h_Xhits_minus_Xtrk_cut2->Fill(x_hit_trk_diff);
            h_Yhits_minus_Ytrk_cut2->Fill(y_hit_trk_diff);
        
        	// track-base e' energy
        	double res_trk= (scatMC.E()-scatMCmatchREC.E())/scatMC.E();
            h_energy_res_trk_cut2->Fill(scatMC.E(), res_trk);
        	h_energy_REC_trk_cut2->Fill(scatMCmatchREC.E());

            // loop over track again;
            for(int itrk=0;itrk<reco_pz_array.GetSize();itrk++)
            {
            	TVector3 trk(reco_px_array[itrk],reco_py_array[itrk],reco_pz_array[itrk]);
            	particle.SetVectM(trk,MASS_PION); 
            	if(itrk!=rec_elect_index) 
                {
                	hfs += particle; 
                	h_eta_REC_trk->Fill(trk.Eta());
                	if(fabs(trk.Eta())<3.0)
                    {
                		if(reco_charge_array[itrk]>0) kplusREC.SetVectM(trk,MASS_KAON);
                		if(reco_charge_array[itrk]<0) kminusREC.SetVectM(trk,MASS_KAON);
                	}
            	}
            }
        
            // 4vector of VM;
            if(kplusREC.E()!=0. && kminusREC.E()!=0.)
            {
            	vmREC=kplusREC+kminusREC;
            }
            h_VM_mass_REC_cut2->Fill(vmREC.M());
            h_VM_pt_REC_cut2->Fill(vmREC.Pt());
        	h_VM_pz_REC_cut2->Fill(vmREC.Pz());
        	h_VM_p_REC_cut2->Fill(vmREC.P());
            h_VM_Epz_MC_cut2->Fill(vmMC.E()-vmMC.Pz());
            h_VM_Epz_REC_cut2->Fill(vmREC.E()-vmREC.Pz());
            h_VM_Epz_response_cut2->Fill(vmREC.E()-vmREC.Pz(),vmMC.E()-vmMC.Pz());
            h_VM_pt_response_cut2->Fill(vmREC.Pt(),vmMC.Pt());

            // e' p from emcal
            h_e_pt_REC_EEMC_cut2->Fill(scatClusEREC.Pt());
        	h_e_pz_REC_EEMC_cut2->Fill(scatClusEREC.Pz());
        	h_e_p_REC_EEMC_cut2->Fill(scatClusEREC.P());
            h_e_pt_res_cut2->Fill(scatMC.Pt(),(scatClusEREC.Pt()-scatMC.Pt())/scatMC.Pt());
            h_e_pz_res_cut2->Fill(scatMC.Pz(),(scatClusEREC.Pz()-scatMC.Pz())/scatMC.Pz());
            h_e_p_res_cut2->Fill(scatMC.P(),(scatClusEREC.P()-scatMC.P())/scatMC.P());
            h_e_pz_response_cut2->Fill(scatClusEREC.Pz(),scatMC.Pz());
        
        	// e' p from track
            h_e_pt_REC_trk_cut2->Fill(scatMCmatchREC.Pt());
        	h_e_pz_REC_trk_cut2->Fill(scatMCmatchREC.Pz());
        	h_e_p_REC_trk_cut2->Fill(scatMCmatchREC.P());

            // E_cal vs p_trk
            h_EvsP_REC_cut2->Fill(scatMCmatchREC.P(),scatREC.E());
            h_EvsP_MC_cut2->Fill(scatMC.P(),scatREC.E());

            // cluster-base DIS kine;
            TLorentzVector qbeamREC=ebeam-scatClusEREC; // emcal e' reco
            TLorentzVector qbeamREC_trk=ebeam-scatMCmatchREC; // track e' reco
            double Q2REC=-(qbeamREC).Mag2();  
            double Q2REC_trk = -(qbeamREC_trk).Mag2();
            double pqREC=abeam.Dot(qbeamREC); 
            double yREC=pqREC/abeam.Dot(ebeam); 
            double pqREC_trk=abeam.Dot(qbeamREC_trk); 
            double yREC_trk=pqREC_trk/abeam.Dot(ebeam); 
            h_Q2REC_e_EEMC_cut2->Fill(Q2REC);
            h_Q2REC_e_trk_cut2->Fill(Q2REC_trk);
            h_yREC_e_EEMC_cut2->Fill(yREC);
            h_yREC_e_trk_cut2->Fill(yREC_trk);

            // Q2 emcal reco resolution
        	res = (Q2 - Q2REC)/Q2;
        	h_Q2_res_cut2->Fill(Q2, res);
            h_Q2_response_cut2->Fill(Q2REC,Q2);
            h_dQ2overQ2_REC_cut2->Fill(res);

            // y emcal reco resolution
            res = (y - yREC)/y;
            h_y_res_cut2->Fill(y,res);
            h_y_response_cut2->Fill(yREC,y);
            h_dyOvery_REC_cut2->Fill(res);

        	// track-base Epz scat' e
            double EpzREC_trk = (scatMCmatchREC+hfs).E() - (scatMCmatchREC+hfs).Pz();
            h_Epz_REC_trk_cut2->Fill(EpzREC_trk);
            double EpzREC = (scatClusEREC+hfs).E() - (scatMCmatchREC+hfs).Pz();
            h_Epz_REC_cut2->Fill(EpzREC);
            h_Epz_response_cut2->Fill(EpzREC,scatMC.E()-scatMC.Pz());
            h_Epz_MC_cut2->Fill(scatMC.E()-scatMC.Pz());
    
        	// E over p  
            double EoverP = scatClusEREC.E() / scatMCmatchREC.P();
            h_EcalOverPtrk_cut2->Fill(scatClusEREC.E()/scatMCmatchREC.P());
            h_EtrkOverPcal_cut2->Fill(scatMCmatchREC.E()/scatClusEREC.P());
            h_EoverP_response_cut2->Fill(scatClusEREC.E()/scatMCmatchREC.P(),scatMC.E()/scatMC.P());
    
            if( EpzREC>27||EpzREC<40 )
            {
                if( EoverP>0.8||EoverP<1.18 )
                {	
                    if(Q2REC>1.||Q2REC<10.)
                    {
                        if(yREC>0.01||yREC<0.85)
                        {
                            // VM rec
                            if(vmREC.E()!=0){
                            // select phi mass and rapidity window 
                        	if( fabs(vmREC.M()-1.02)<0.02&& fabs(vmREC.Rapidity())<3.5)
                        	{
                            	// 2 versions: track and energy cluster:
                        		double t_trk_REC = giveme_t_method_L(ebeam,scatMCmatchREC,abeam,vmREC);  // method L (track e')
                            	double t_REC = giveme_t_method_L(ebeam,scatClusEREC,abeam,vmREC); // method L (EEMC e')
                            	h_t_REC_trk_cut2->Fill(t_trk_REC);
                            	h_t_REC_EEMC_cut2->Fill(t_REC);
                                h_t_response_EEMC_cut2->Fill(t_REC,t_MC);
                                h_t_response_trk_cut2->Fill(t_trk_REC,t_MC);
                                h_t_response_trk_EEMC_cut2->Fill(t_trk_REC,t_REC);

                        		// new method
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
                    		h_t_REC_2d_wRES_cut2->Fill(qx_rec,qy_rec);
                    		h_t_REC_wRES_cut2->Fill(t_total_rec);

                        	// t track resolution 
                    		res = (t_MC-t_trk_REC)/t_MC;
                    		h_t_res_trk_cut2->Fill(t_MC, res);
	
                        	// t EEMC resolution;
                    		res = (t_MC-t_REC)/t_MC;
                    		h_t_res_EEMC_cut2->Fill(t_MC, res);

                    		// apply cut with // subtracted out of t
                    		double theta_rec = atan(fabs(qx_rec)/fabs(qy_rec));
                    		if(fabs(theta_rec)<PI/2)
                    		{
                    			h_t_REC_wRES_cut_pi2_cut2->Fill(t_total_rec);
                    			h_t_REC_2d_wRES_cut_pi2_cut2->Fill(qx_rec,qy_rec);
                    		}
                    		if(fabs(theta_rec)<PI/3)
                    		{
                    			h_t_REC_wRES_cut_pi3_cut2->Fill(t_total_rec);
                    			h_t_REC_2d_wRES_cut_pi3_cut2->Fill(qx_rec,qy_rec);
                    		}
                    		if(fabs(theta_rec)<PI/4)
                    		{
                    			h_t_REC_wRES_cut_pi4_cut2->Fill(t_total_rec);
                    			h_t_REC_2d_wRES_cut_pi4_cut2->Fill(qx_rec,qy_rec);
                    		}
                    		if(fabs(theta_rec)<PI/6)
                    		{
                    			h_t_REC_wRES_cut_pi6_cut2->Fill(t_total_rec);
                    			h_t_REC_2d_wRES_cut_pi6_cut2->Fill(qx_rec,qy_rec);
                    		}
                    		if(fabs(theta_rec)<PI/9)
                    		{
                    			h_t_REC_wRES_cut_pi9_cut2->Fill(t_total_rec);
                    			h_t_REC_2d_wRES_cut_pi9_cut2->Fill(qx_rec,qy_rec);
                    		}
                    		if(fabs(theta_rec)<PI/12)
                    		{
                    			h_t_REC_wRES_cut_pi12_cut2->Fill(t_total_rec);
                    			h_t_REC_2d_wRES_cut_pi12_cut2->Fill(qx_rec,qy_rec);
                    		}
                    		if(fabs(theta_rec)<PI/16)
                    		{
                    			h_t_REC_wRES_cut_pi16_cut2->Fill(t_total_rec);
                    			h_t_REC_2d_wRES_cut_pi16_cut2->Fill(qx_rec,qy_rec);
                    		}
                            if(fabs(theta_rec)<PI/20)
                    		{
                    			h_t_REC_wRES_cut_pi20_cut2->Fill(t_total_rec);
                    			h_t_REC_2d_wRES_cut_pi20_cut2->Fill(qx_rec,qy_rec);
                    		}
                            if(fabs(theta_rec)<PI/24)
                    		{
                    			h_t_REC_wRES_cut_pi24_cut2->Fill(t_total_rec);
                    			h_t_REC_2d_wRES_cut_pi24_cut2->Fill(qx_rec,qy_rec);
                    		}
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