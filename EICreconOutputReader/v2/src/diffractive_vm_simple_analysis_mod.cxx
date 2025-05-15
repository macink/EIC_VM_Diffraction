#include "pleaseIncludeMe.h"
#include <TVector3.h>
#include <TVector.h>
#include <cmath>
#include "TF1.h"
#include "TF2.h"
#include "TSystem.h"

auto giveme_t_method_L(TLorentzVector eIn, TLorentzVector eOut, TLorentzVector pIn, TLorentzVector vmOut)
{
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

auto giveme_t_new_method(TLorentzVector eIn, 
	TLorentzVector eOut, 
	TLorentzVector pIn, 
	TLorentzVector vmOut)
{
	TLorentzVector aInVec(pIn.Px()*197,pIn.Py()*197,pIn.Pz()*197,sqrt(pIn.Px()*197*pIn.Px()*197 + pIn.Py()*197*pIn.Py()*197 + pIn.Pz()*197*pIn.Pz()*197 + MASS_AU197*MASS_AU197) );
	double method_L = 0;
	TLorentzVector a_beam_scattered = aInVec-(vmOut+eOut-eIn);
	double p_Aplus = a_beam_scattered.E()+a_beam_scattered.Pz();
	double p_TAsquared = TMath::Power(a_beam_scattered.Pt(),2);
	double p_Aminus = (MASS_AU197*MASS_AU197 + p_TAsquared) / p_Aplus;
	TLorentzVector a_beam_scattered_corr; 
	a_beam_scattered_corr.SetPxPyPzE(a_beam_scattered.Px(),a_beam_scattered.Py(),(p_Aplus-p_Aminus)/2., (p_Aplus+p_Aminus)/2. );
	method_L = -(a_beam_scattered_corr-aInVec).Mag2();
	// Added to return 4 vector 
	TLorentzVector method_L_4vect = -(a_beam_scattered_corr-aInVec);
	return method_L_4vect;
}


int diffractive_vm_simple_analysis_mod(TString rec_file, TString outputfile)
{	
// read our configuration	
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

    /*
TString name_of_input = (TString) rec_file;
std::cout << "Input file = " << name_of_input << endl;
auto tree = new TChain("events");
tree->Add(name_of_input);
TTreeReader tree_reader(tree);       // !the tree reader
*/
    
TTreeReaderArray<int> mc_genStatus_array = {tree_reader, "MCParticles.generatorStatus"};
// MC particle pz array for each MC particle
TTreeReaderArray<float> mc_px_array = {tree_reader, "MCParticles.momentum.x"};
TTreeReaderArray<float> mc_py_array = {tree_reader, "MCParticles.momentum.y"};
TTreeReaderArray<float> mc_pz_array = {tree_reader, "MCParticles.momentum.z"};
TTreeReaderArray<double> mc_mass_array = {tree_reader, "MCParticles.mass"};
TTreeReaderArray<int> mc_pdg_array = {tree_reader, "MCParticles.PDG"};

//Reconstructed EcalEndcapNClusters
TTreeReaderArray<float> em_energy_array = {tree_reader, "EcalEndcapNClusters.energy"};
TTreeReaderArray<float> em_x_array = {tree_reader, "EcalEndcapNClusters.position.x"};
TTreeReaderArray<float> em_y_array = {tree_reader, "EcalEndcapNClusters.position.y"};
TTreeReaderArray<float> emhits_x_array = {tree_reader, "EcalEndcapNRecHits.position.x"};
TTreeReaderArray<float> emhits_y_array = {tree_reader, "EcalEndcapNRecHits.position.y"};
TTreeReaderArray<float> emhits_energy_array = {tree_reader, "EcalEndcapNRecHits.energy"};

TTreeReaderArray<unsigned int> em_rec_id_array = {tree_reader, "EcalEndcapNClusterAssociations.recID"};
TTreeReaderArray<unsigned int> em_sim_id_array = {tree_reader, "EcalEndcapNClusterAssociations.simID"};

// Reconstructed particles pz array for each reconstructed particle
TTreeReaderArray<float> reco_px_array = {tree_reader, "ReconstructedChargedParticles.momentum.x"};
TTreeReaderArray<float> reco_py_array = {tree_reader, "ReconstructedChargedParticles.momentum.y"};
TTreeReaderArray<float> reco_pz_array = {tree_reader, "ReconstructedChargedParticles.momentum.z"};
TTreeReaderArray<float> reco_charge_array = {tree_reader, "ReconstructedChargedParticles.charge"};

TTreeReaderArray<unsigned int> rec_id = {tree_reader, "ReconstructedChargedParticleAssociations.recID"};
TTreeReaderArray<unsigned int> sim_id = {tree_reader, "ReconstructedChargedParticleAssociations.simID"};

TString output_name_dir = outputfile+"_output.root";
cout << "Output file = " << output_name_dir << endl;
TFile* output = new TFile(output_name_dir,"RECREATE");

//events
TH1D* h_Q2_e = new TH1D("h_Q2_e",";Q^{2}_{e,MC}",100,0,20);
TH1D* h_y_e = new TH1D("h_y_e",";y_{e,MC}",100,0,1);
TH1D* h_energy_MC = new TH1D("h_energy_MC",";E_{MC} (GeV)",100,0,20);
TH1D* h_t_MC = new TH1D("h_t_MC",";t_{MC}; counts",100,0,0.2);

TH1D* h_Q2REC_e = new TH1D("h_Q2REC_e",";Q^{2}_{e,REC}",100,0,20);
TH1D* h_yREC_e = new TH1D("h_yREC_e",";y_{e,REC}",100,0,1);
TH1D* h_energy_REC = new TH1D("h_energy_REC",";E_{REC} (GeV)",100,0,20);
TH1D* h_trk_energy_REC = new TH1D("h_trk_energy_REC",";E_{REC} (GeV)",100,0,20);
TH1D* h_trk_Epz_REC = new TH1D("h_trk_Epz_REC",";E - p_{z} (GeV)",200,0,50);

//track
TH1D* h_eta = new TH1D("h_eta",";#eta",100,-5,5);
TH2D* h_trk_energy_res = new TH2D("h_trk_energy_res",";E_{MC} (GeV); E_{MC}-E_{REC}/E_{MC} track-base ",100,0,20,1000,-1,1);
TH2D* h_trk_Pt_res = new TH2D("h_trk_Pt_res",";p_{T,MC} (GeV); P_{T,MC}-P_{T,REC}/P_{T,MC} track-base ",100,0,15,1000,-1,1);
TH1D* h_Epz_REC = new TH1D("h_Epz_REC",";E - p_{z} (GeV)",200,0,50);

//VM & t
TH1D* h_VM_mass_REC = new TH1D("h_VM_mass_REC",";mass (GeV)",200,0,4);
TH1D* h_VM_pt_REC = new TH1D("h_VM_pt_REC",";p_{T} (GeV/c)",200,0,2);
TH2D* h_VM_res = new TH2D("h_VM_res",";p_{T,MC} (GeV); p_{T,MC}-E_{T,REC}/p_{T,MC}",100,0,2,1000,-1,1);
TH1D* h_t_REC = new TH1D("h_t_REC",";t_{REC} (GeV^{2}); counts",100,0,0.2);
TH1D* h_t_trk_REC = new TH1D("h_t_trk_REC",";t_{REC}(GeV^{2}) track-base; counts",100,0,0.2);
TH1D* h_t_combo_REC = new TH1D("h_t_combo_REC",";t_{combo,REC}(GeV^{2}); counts",100,0,0.2);
TH2D* h_t_res = new TH2D("h_t_res",";t_{MC} (GeV^{2}); t_{MC}-t_{REC}/t_{MC}",100,0,0.2,1000,-10,10);
TH2D* h_trk_t_res = new TH2D("h_trk_t_res",";t_{MC} (GeV^{2}); t_{MC}-t_{REC}/t_{MC} track-base",100,0,0.2,1000,-10,10);
TH2D* h_t_2D = new TH2D("h_t_2D",";t_{MC} (GeV^{2}); t_{REC} (GeV^{2}) track-base",100,0,0.2,100,0,0.2);
TH2D* h_t_REC_2D = new TH2D("h_t_REC_2D",";t_{trk,REC} (GeV^{2}); t_{EEMC,REC} (GeV^{2})",100,0,0.2,100,0,0.2);
TH2D* h_t_RECMC_2D = new TH2D("h_t_RECMC_2D",";t_{MC} (GeV^{2}); t_{trk,REC} / t_{EEMC,REC} ",100,0,0.2,200,-10,10);

// combined plots
//t distribution
TH1D* h_t_REC_new_method = new TH1D("h_t_REC_new_method",";t_{REC} (GeV^{2}); counts",100,0,0.2);
TH1D* h_t_REC_new_method_wZ = new TH1D("h_t_REC_new_method_wZ",";t_{REC} (GeV^{2}); counts",100,0,0.2);
//t distribution with resolution only and with angle cut only (for testing)
TH1D* h_t_REC_wRES = new TH1D("h_t_REC_wRES",";t_{REC} (GeV^{2}); counts",100,0,0.2);
TH1D* h_t_REC_wCUT = new TH1D("h_t_REC_wCUT",";t_{REC} (GeV^{2}); counts",100,0,0.2);
//t distribution with cut only (for testing)
TH1D* h_t_REC_wCUT_pi2 = new TH1D("h_t_REC_wCUT_pi2",";t_{REC} (GeV^{2}); counts",100,0,0.2);
TH1D* h_t_REC_wCUT_pi3 = new TH1D("h_t_REC_wCUT_pi3",";t_{REC} (GeV^{2}); counts",100,0,0.2);
TH1D* h_t_REC_wCUT_pi4 = new TH1D("h_t_REC_wCUT_pi4",";t_{REC} (GeV^{2}); counts",100,0,0.2);
TH1D* h_t_REC_wCUT_pi6 = new TH1D("h_t_REC_wCUT_pi6",";t_{REC} (GeV^{2}); counts",100,0,0.2);
TH1D* h_t_REC_wCUT_pi9 = new TH1D("h_t_REC_wCUT_pi9",";t_{REC} (GeV^{2}); counts",100,0,0.2);
TH1D* h_t_REC_wCUT_pi12 = new TH1D("h_t_REC_wCUT_pi12",";t_{REC} (GeV^{2}); counts",100,0,0.2);
//t distribution with resolution and angle cut
TH1D* h_t_REC_wRES_cut = new TH1D("h_t_REC_wRES_cut",";t_{REC} (GeV^{2}); counts",100,0,0.2);
TH1D* h_t_REC_wRES_cut_pi2 = new TH1D("h_t_REC_wRES_cut_pi2",";t_{REC} (GeV^{2}); counts",100,0,0.2);
TH1D* h_t_REC_wRES_cut_pi3 = new TH1D("h_t_REC_wRES_cut_pi3",";t_{REC} (GeV^{2}); counts",100,0,0.2);
TH1D* h_t_REC_wRES_cut_pi4 = new TH1D("h_t_REC_wRES_cut_pi4",";t_{REC} (GeV^{2}); counts",100,0,0.2);
TH1D* h_t_REC_wRES_cut_pi6 = new TH1D("h_t_REC_wRES_cut_pi6",";t_{REC} (GeV^{2}); counts",100,0,0.2);
TH1D* h_t_REC_wRES_cut_pi9 = new TH1D("h_t_REC_wRES_cut_pi9",";t_{REC} (GeV^{2}); counts",100,0,0.2);
TH1D* h_t_REC_wRES_cut_pi12 = new TH1D("h_t_REC_wRES_cut_pi12",";t_{REC} (GeV^{2}); counts",100,0,0.2);
TH1D* h_t_REC_wRES_cut_pi16 = new TH1D("h_t_REC_wRES_cut_pi16",";t_{REC} (GeV^{2}); counts",100,0,0.2);
TH1D* h_t_REC_wRES_cut_pi20 = new TH1D("h_t_REC_wRES_cut_pi20",";t_{REC} (GeV^{2}); counts",100,0,0.2);
TH1D* h_t_REC_wRES_cut_pi24 = new TH1D("h_t_REC_wRES_cut_pi24",";t_{REC} (GeV^{2}); counts",100,0,0.2);
//t distribution 2d
TH2D* h_t_REC_2d = new TH2D("h_t_REC_2d",";#sqrt{t_{x}} (GeV); #sqrt{t_{y}} (GeV)",100,0,0.2,100,0,0.2);
TH2D* h_t_REC_2d_wZ = new TH2D("h_t_REC_2d_wZ",";#sqrt{t_{x}} (GeV); #sqrt{t_{y}} (GeV)",100,0,0.2,100,0,0.2);
//t distribution 2d with resolution only 
TH2D* h_t_REC_2d_wRES = new TH2D("h_t_REC_2d_wRES",";#sqrt{t_{x}} (GeV); #sqrt{t_{y}} (GeV)",100,0,0.2,100,0,0.2);
//t distribution 2d with resolution and angle cut
TH2D* h_t_REC_2d_wRES_cut = new TH2D("h_t_REC_2d_wRES_cut",";#sqrt{t_{x}} (GeV); #sqrt{t_{y}} (GeV)",100,0,0.2,100,0,0.2);
TH2D* h_t_REC_2d_wRES_cut_pi2 = new TH2D("h_t_REC_2d_wRES_cut_pi2",";#sqrt{t_{x}} (GeV); #sqrt{t_{y}} (GeV)",100,0,0.2,100,0,0.2);
TH2D* h_t_REC_2d_wRES_cut_pi3 = new TH2D("h_t_REC_2d_wRES_cut_pi3",";#sqrt{t_{x}} (GeV); #sqrt{t_{y}} (GeV)",100,0,0.2,100,0,0.2);
TH2D* h_t_REC_2d_wRES_cut_pi4 = new TH2D("h_t_REC_2d_wRES_cut_pi4",";#sqrt{t_{x}} (GeV); #sqrt{t_{y}} (GeV)",100,0,0.2,100,0,0.2);
TH2D* h_t_REC_2d_wRES_cut_pi6 = new TH2D("h_t_REC_2d_wRES_cut_pi6",";#sqrt{t_{x}} (GeV); #sqrt{t_{y}} (GeV)",100,0,0.2,100,0,0.2);
TH2D* h_t_REC_2d_wRES_cut_pi9 = new TH2D("h_t_REC_2d_wRES_cut_pi9",";#sqrt{t_{x}} (GeV); #sqrt{t_{y}} (GeV)",100,0,0.2,100,0,0.2);
TH2D* h_t_REC_2d_wRES_cut_pi12 = new TH2D("h_t_REC_2d_wRES_cut_pi12",";#sqrt{t_{x}} (GeV); #sqrt{t_{y}} (GeV)",100,0,0.2,100,0,0.2);
TH2D* h_t_REC_2d_wRES_cut_pi16 = new TH2D("h_t_REC_2d_wRES_cut_pi16",";#sqrt{t_{x}} (GeV); #sqrt{t_{y}} (GeV)",100,0,0.2,100,0,0.2);
TH2D* h_t_REC_2d_wRES_cut_pi20 = new TH2D("h_t_REC_2d_wRES_cut_pi20",";#sqrt{t_{x}} (GeV); #sqrt{t_{y}} (GeV)",100,0,0.2,100,0,0.2);
TH2D* h_t_REC_2d_wRES_cut_pi24 = new TH2D("h_t_REC_2d_wRES_cut_pi24",";#sqrt{t_{x}} (GeV); #sqrt{t_{y}} (GeV)",100,0,0.2,100,0,0.2);
//t distribution with angle cut only
TH2D* h_t_REC_2d_wCUT = new TH2D("h_t_REC_2d_wCUT",";#sqrt{t_{x}} (GeV); #sqrt{t_{y}} (GeV)",100,0,0.2,100,0,0.2);
TH2D* h_t_REC_2d_wCUT_pi2 = new TH2D("h_t_REC_2d_wCUT_pi2",";#sqrt{t_{x}} (GeV); #sqrt{t_{y}} (GeV)",100,0,0.2,100,0,0.2);
TH2D* h_t_REC_2d_wCUT_pi3 = new TH2D("h_t_REC_2d_wCUT_pi3",";#sqrt{t_{x}} (GeV); #sqrt{t_{y}} (GeV)",100,0,0.2,100,0,0.2);
TH2D* h_t_REC_2d_wCUT_pi4 = new TH2D("h_t_REC_2d_wCUT_pi4",";#sqrt{t_{x}} (GeV); #sqrt{t_{y}} (GeV)",100,0,0.2,100,0,0.2);
TH2D* h_t_REC_2d_wCUT_pi6 = new TH2D("h_t_REC_2d_wCUT_pi6",";#sqrt{t_{x}} (GeV); #sqrt{t_{y}} (GeV)",100,0,0.2,100,0,0.2);
TH2D* h_t_REC_2d_wCUT_pi9 = new TH2D("h_t_REC_2d_wCUT_pi9",";#sqrt{t_{x}} (GeV); #sqrt{t_{y}} (GeV)",100,0,0.2,100,0,0.2);
TH2D* h_t_REC_2d_wCUT_pi12 = new TH2D("h_t_REC_2d_wCUT_pi12",";#sqrt{t_{x}} (GeV); #sqrt{t_{y}} (GeV)",100,0,0.2,100,0,0.2);
TH2D* h_t_REC_2d_wCUT_pi16 = new TH2D("h_t_REC_2d_wCUT_pi16",";#sqrt{t_{x}} (GeV); #sqrt{t_{y}} (GeV)",100,0,0.2,100,0,0.2);
TH2D* h_t_REC_2d_wCUT_pi20 = new TH2D("h_t_REC_2d_wCUT_pi20",";#sqrt{t_{x}} (GeV); #sqrt{t_{y}} (GeV)",100,0,0.2,100,0,0.2);
TH2D* h_t_REC_2d_wCUT_pi24 = new TH2D("h_t_REC_2d_wCUT_pi24",";#sqrt{t_{x}} (GeV); #sqrt{t_{y}} (GeV)",100,0,0.2,100,0,0.2);

//energy clus
TH2D* h_emClus_position_REC = new TH2D("h_emClus_position_REC",";x (mm);y (mm)",80,-800,800,80,-800,800);
TH2D* h_emHits_position_REC = new TH2D("h_emHits_position_REC",";x (mm);y (mm)",80,-800,800,80,-800,800);
TH2D* h_energy_res = new TH2D("h_energy_res",";E_{MC} (GeV); E_{MC}-E_{REC}/E_{MC} emcal",100,0,20,1000,-1,1);
TH1D* h_energy_calibration_REC = new TH1D("h_energy_calibration_REC",";E (GeV)",200,0,2);
TH1D* h_EoverP_REC = new TH1D("h_EoverP_REC",";E/p",200,0,2);
TH1D* h_ClusOverHit_REC = new TH1D("h_ClusOverHit_REC",";cluster energy / new cluster energy",200,0,2);

//tree_reader.SetEntriesRange(0, tree->GetEntries());
chain->GetEntries();
while (tree_reader.Next()) {
	/*
	Beam particles
	*/
	TLorentzVector ebeam(0,0,0,0);
	TLorentzVector pbeam(0,0,0,0);

	TLorentzVector vmMC(0,0,0,0);
	TLorentzVector kplusMC(0,0,0,0);
	TLorentzVector kminusMC(0,0,0,0);

	//MC level
	TLorentzVector scatMC(0,0,0,0);
	int mc_elect_index=-1;
	double maxPt=-99.;

	for(int imc=0;imc<mc_px_array.GetSize();imc++){
        TVector3 mctrk(static_cast<double>(mc_px_array[imc]), 
               static_cast<double>(mc_py_array[imc]), 
               static_cast<double>(mc_pz_array[imc]));
		//TVector3 mctrk(mc_px_array[imc],mc_py_array[imc],mc_pz_array[imc]);	
		if(mc_genStatus_array[imc]==4){//4 is Sartre.
			if(mc_pdg_array[imc]==11) ebeam.SetVectM(mctrk, MASS_ELECTRON);
				if(mc_pdg_array[imc]==2212) pbeam.SetVectM(mctrk, MASS_PROTON);
		}
		if(mc_genStatus_array[imc]!=1) continue;
		if(mc_pdg_array[imc]==11 	
			&& mctrk.Perp()>maxPt){
			maxPt=mctrk.Perp();
			mc_elect_index=imc;
			scatMC.SetVectM(mctrk,mc_mass_array[imc]);
		}
		if(mc_pdg_array[imc]==321
			&& mc_genStatus_array[imc]==1) kplusMC.SetVectM(mctrk,MASS_KAON);
        //cout << mc_pdg_array[imc] << endl;
		if(mc_pdg_array[imc]==-321
			&& mc_genStatus_array[imc]==1) kminusMC.SetVectM(mctrk,MASS_KAON);
        //cout << mc_pdg_array[imc] << endl;
	}
	vmMC=kplusMC+kminusMC; 
	cout << "this is working!! " << endl;

	//protection.
	if(ebeam.E()==pbeam.E() && ebeam.E()==0) {
		std::cout << "problem with MC incoming beams" << std::endl;
		continue;
	}
	TLorentzVector qbeam=ebeam-scatMC;
	double Q2=-(qbeam).Mag2();  
	double pq=pbeam.Dot(qbeam);
	double y= pq/pbeam.Dot(ebeam);

    //MC level phase space cut
	if(Q2<1.||Q2>10.) continue;
	if(y<0.01||y>0.85) continue;
	h_Q2_e->Fill(Q2);
	h_y_e->Fill(y);
	h_energy_MC->Fill(scatMC.E());
	
	double keep = 0;
	double t_MC=0.;
	if(vmMC.E()!=0 
		&& fabs(vmMC.Rapidity())<3.5)
	{
		double method_E = -(qbeam-vmMC).Mag2();  // t = -(p_e - p_e' -p_v)^2
		t_MC=method_E;
		h_t_MC->Fill( method_E );

		// new method
		TLorentzVector T = -ebeam + scatMC + vmMC;  // true t (method E)
		
		// define nHat direction
		TVector3 eScattered_momentum = scatMC.Vect(); 
		TVector3 e_momentum = ebeam.Vect();
		TVector3 nHat = e_momentum.Cross(eScattered_momentum); // nHat = p_e x p_e'
		nHat = nHat.Unit();
		double pv = vmMC.Vect().Dot(nHat);

		TVector3 py_vector = nHat*pv;
		TLorentzVector py(py_vector.X(),py_vector.Y(),0,0);
        TLorentzVector px = T - py;
        TLorentzVector pz(0,0,px.Z(),px.E());
		TLorentzVector px_true = T - py - pz;

        double tx = -(px).Mag2();
		double tx_true = -(px_true).Mag2();
		double ty = -(py).Mag2();
        double t_true = tx_true+ty;
		double t = tx+ty;
		double qx = sqrt(tx);
        double qx_true = sqrt(tx_true);
		double qy = sqrt(ty);

		h_t_REC_2d->Fill(qx_true,qy);
        h_t_REC_2d_wZ->Fill(qx,qy);
        h_t_REC_new_method_wZ->Fill(t);
		h_t_REC_new_method->Fill(t_true);
        //cout << "t:  " << t << endl;
		
		// apply cut
		double theta = atan(fabs(qx_true)/fabs(qy));
        //test with cut only -> should be same as truth with out resolution
		if(fabs(theta)<PI/9) 
		{
			h_t_REC_wCUT->Fill(t_true); 
            h_t_REC_2d_wCUT->Fill(qx_true,qy);
		}
		if(fabs(theta)<PI/2) 
		{
			h_t_REC_wCUT_pi2->Fill(t_true);
			h_t_REC_2d_wCUT_pi2->Fill(qx_true,qy);
		}
		if(fabs(theta)<PI/3) 
		{
			h_t_REC_wCUT_pi3->Fill(t_true);
			h_t_REC_2d_wCUT_pi3->Fill(qx_true,qy);
		}
		if(fabs(theta)<PI/4) 
		{
			h_t_REC_wCUT_pi4->Fill(t_true);
			h_t_REC_2d_wCUT_pi4->Fill(qx_true,qy);
		}
		if(fabs(theta)<PI/6) 
		{
			h_t_REC_wCUT_pi6->Fill(t_true);
			h_t_REC_2d_wCUT_pi6->Fill(qx_true,qy);
		}
		if(fabs(theta)<PI/9) 
		{
			h_t_REC_wCUT_pi9->Fill(t_true);
			h_t_REC_2d_wCUT_pi9->Fill(qx_true,qy);
		}
		if(fabs(theta)<PI/12) 
		{
			h_t_REC_wCUT_pi12->Fill(t_true);
			h_t_REC_2d_wCUT_pi12->Fill(qx_true,qy);
		}
        if(fabs(theta)<PI/16) 
		{
			h_t_REC_2d_wCUT_pi16->Fill(qx_true,qy);
		}
        if(fabs(theta)<PI/20) 
		{
			h_t_REC_2d_wCUT_pi20->Fill(qx_true,qy);
		}
        if(fabs(theta)<PI/24) 
		{
			h_t_REC_2d_wCUT_pi24->Fill(qx_true,qy);
		}
	
	}
	

	//rec level
	//leading cluster
	double maxEnergy=-99.;
	double xpos=-999.;
	double ypos=-999.;
	for(int iclus=0;iclus<em_energy_array.GetSize();iclus++){
		if(em_energy_array[iclus]>maxEnergy){
			maxEnergy=em_energy_array[iclus];
			xpos=em_x_array[iclus];
			ypos=em_y_array[iclus];
		}
	}
	//leading hit energy
	double maxHitEnergy=0.01;//threshold 10 MeV
	double xhitpos=-999.;
	double yhitpos=-999.;
	int hit_index=-1;
	for(int ihit=0;ihit<emhits_energy_array.GetSize();ihit++){	
		if(emhits_energy_array[ihit]>maxHitEnergy){
			maxHitEnergy=emhits_energy_array[ihit];
			        xhitpos=emhits_x_array[ihit];
			        yhitpos=emhits_y_array[ihit];
			        hit_index=ihit;
		}
	}
	//sum over all 3x3 towers around the leading tower
	double xClus=xhitpos*maxHitEnergy;
	double yClus=yhitpos*maxHitEnergy;
	for(int ihit=0;ihit<emhits_energy_array.GetSize();ihit++){
		double hitenergy=emhits_energy_array[ihit];
		double x=emhits_x_array[ihit];
		double y=emhits_y_array[ihit];
		double d=sqrt( (x-xhitpos)*(x-xhitpos) + (y-yhitpos)*(y-yhitpos));
		if(d<70. && ihit!=hit_index && hitenergy>0.01)  {
			maxHitEnergy+=hitenergy;//clustering around leading tower 3 crystal = 60mm.
			xClus+=x*hitenergy;
			yClus+=y*hitenergy;
		}
	}
	
	h_ClusOverHit_REC->Fill( maxEnergy / maxHitEnergy );
	//weighted average cluster position.
	xClus = xClus/maxHitEnergy;
	yClus = yClus/maxHitEnergy;
	double radius=sqrt(xClus*xClus+yClus*yClus);
	if( radius>550. ) continue; //geometric acceptance cut
	//4.4% energy calibration.
	double clusEnergy=1.044*maxHitEnergy; 

	h_energy_REC->Fill(clusEnergy);
	//ratio of reco / truth Energy
	h_energy_calibration_REC->Fill( clusEnergy / scatMC.E() );
	//energy resolution
	double res= (scatMC.E()-clusEnergy)/scatMC.E();
	h_energy_res->Fill(scatMC.E(), res);
	h_emClus_position_REC->Fill(xpos,ypos);//default clustering position
	h_emHits_position_REC->Fill(xClus,yClus); //self clustering position
	
	//association of rec level scat' e
	int rec_elect_index=-1;
	for(int i=0;i<sim_id.GetSize();i++){
		if(sim_id[i]==mc_elect_index){
			//find the rec_id
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
    //track loop
	for(int itrk=0;itrk<reco_pz_array.GetSize();itrk++){
		TVector3 trk(reco_px_array[itrk],reco_py_array[itrk],reco_pz_array[itrk]);
		if(rec_elect_index!=-1
			&& itrk==rec_elect_index){
			scatMCmatchREC.SetVectM(trk,MASS_ELECTRON);//Reserved to calculate t.
		}
		if(trk.Mag()>maxP){
			//track-base 4 vector
			maxP=trk.Mag();
			scatREC.SetVectM(trk,MASS_ELECTRON);

			//use emcal energy to define 4 vector
			double p = sqrt(clusEnergy*clusEnergy- MASS_ELECTRON*MASS_ELECTRON );
			double eta=scatREC.Eta();
			double phi=scatREC.Phi();
			double pt = TMath::Sin(scatREC.Theta())*p;
			scatClusEREC.SetPtEtaPhiM(pt,eta,phi,MASS_ELECTRON);
		}
	}
	//loop over track again;
	for(int itrk=0;itrk<reco_pz_array.GetSize();itrk++){
		TVector3 trk(reco_px_array[itrk],reco_py_array[itrk],reco_pz_array[itrk]);
		particle.SetVectM(trk,MASS_PION);//assume pions;
		if(itrk!=rec_elect_index) {
    		hfs += particle; //hfs 4vector sum.
    		//selecting phi->kk daughters;
    		h_eta->Fill(trk.Eta());
    		if(fabs(trk.Eta())<3.0){
    			if(reco_charge_array[itrk]>0) kplusREC.SetVectM(trk,MASS_KAON);
    			if(reco_charge_array[itrk]<0) kminusREC.SetVectM(trk,MASS_KAON);
    		}
		}
	}
	//4vector of VM;
	if(kplusREC.E()!=0. && kminusREC.E()!=0.){
		vmREC=kplusREC+kminusREC;
	}

	//track-base e' energy REC;
	h_trk_energy_REC->Fill(scatMCmatchREC.E());
	
	//track-base e' energy resolution;
	res= (scatMC.E()-scatMCmatchREC.E())/scatMC.E();
	h_trk_energy_res->Fill(scatMC.E(), res);
	
	//track-base e' pt resolution;
	res= (scatMC.Pt()-scatMCmatchREC.Pt())/scatMC.Pt();
	h_trk_Pt_res->Fill(scatMC.Pt(), res);

	//track-base Epz scat' e
	double EpzREC= (scatMCmatchREC+hfs).E() - (scatMCmatchREC+hfs).Pz();
	h_trk_Epz_REC->Fill( EpzREC );

	//EEMC cluster Epz scat' e
	EpzREC= (scatClusEREC+hfs).E() - (scatClusEREC+hfs).Pz();
	h_Epz_REC->Fill( EpzREC );

	//E over p
	double EoverP=scatClusEREC.E() / scatMCmatchREC.P();
	h_EoverP_REC->Fill( EoverP );

	//cluster-base DIS kine;
	TLorentzVector qbeamREC=ebeam-scatClusEREC;
	double Q2REC=-(qbeamREC).Mag2();  
	double pqREC=pbeam.Dot(qbeamREC);
	double yREC= pqREC/pbeam.Dot(ebeam);
	h_Q2REC_e->Fill(Q2REC);
	h_yREC_e->Fill(yREC);

	//Event selection:
	if( EpzREC<27||EpzREC>40 ) continue;
	if( EoverP<0.8||EoverP>1.18 ) continue;		

	//MC level phase space cut
	if(Q2REC<1.||Q2REC>10.) continue;
	if(yREC<0.01||yREC>0.85) continue;

	//VM rec
	if(vmREC.E()==0) continue;
	double phi_mass = vmREC.M();
	h_VM_mass_REC->Fill(phi_mass);
	h_VM_pt_REC->Fill(vmREC.Pt());

	
	//select phi mass and rapidity window 
	if( fabs(phi_mass-1.02)<0.02&& fabs(vmREC.Rapidity())<3.5)
	{
    	//2 versions: track and energy cluster:
		double t_trk_REC = giveme_t_method_L(ebeam,scatMCmatchREC,pbeam,vmREC);
    	double t_REC = giveme_t_method_L(ebeam,scatClusEREC,pbeam,vmREC);
    	h_t_trk_REC->Fill( t_trk_REC );
    	h_t_REC->Fill( t_REC);
		//h_t_REC->Fill( t_REC, 2);
    	h_t_REC_2D->Fill(t_trk_REC,t_REC);
    	if((t_trk_REC/t_REC) > 0.5 && (t_trk_REC/t_REC) < 1.5)
		{
    		h_t_combo_REC->Fill( (t_trk_REC+t_REC)/2. );//w=1./(fabs(1.0-(t_trk_REC/t_REC)))
    	} 
    	h_t_RECMC_2D->Fill(t_MC,t_trk_REC/t_REC);

		// new method
		TLorentzVector T_rec = giveme_t_new_method(ebeam,scatClusEREC,pbeam,vmREC); // method L

		// define nHat direction
		TVector3 eScattered_momentum_rec = scatClusEREC.Vect(); 
		TVector3 e_momentum = ebeam.Vect();
		TVector3 nHat_rec = e_momentum.Cross(eScattered_momentum_rec); // nHat = p_e x p_e'
		nHat_rec = nHat_rec.Unit();
		double pv_rec = vmREC.Vect().Dot(nHat_rec);

		TVector3 py_vector_rec = nHat_rec*pv_rec;
		TLorentzVector py_rec(py_vector_rec.X(),py_vector_rec.Y(),0,0);
        TLorentzVector px_rec = T_rec - py_rec;
        TLorentzVector pz_rec(0,0,px_rec.Z(),px_rec.E());
		TLorentzVector px_rec_true = T_rec - py_rec - pz_rec;
		
        double tx_rec = -(px_rec_true).Mag2();
		double ty_rec = -(py_rec).Mag2();
		double t_rec = tx_rec+ty_rec;
		double qx_rec = sqrt(tx_rec);
		double qy_rec = sqrt(ty_rec);

		h_t_REC_2d_wRES->Fill(qx_rec,qy_rec);
		h_t_REC_wRES->Fill(t_rec);

	
    	//t track resolution 
		res= (t_MC-t_trk_REC)/t_MC;
		h_trk_t_res->Fill(t_MC, res);
	
    	//t EEMC resolution;
		res= (t_MC-t_REC)/t_MC;
		h_t_res->Fill(t_MC, res);
	
		//2D t
		h_t_2D->Fill(t_MC,t_trk_REC);

    	//VM pt resolution;
    	res= (vmMC.Pt()-vmREC.Pt())/vmMC.Pt();
		h_VM_res->Fill(vmMC.Pt(), res);

		// apply cut
		double theta_rec = atan(fabs(qx_rec)/fabs(qy_rec));
		if(fabs(theta_rec)<PI/6)
		{
    		h_t_REC_wRES_cut->Fill(t_rec);
			h_t_REC_2d_wRES_cut->Fill(qx_rec,qy_rec);
		}
		if(fabs(theta_rec)<PI/2)
		{
			h_t_REC_wRES_cut_pi2->Fill(t_rec);
			h_t_REC_2d_wRES_cut_pi2->Fill(qx_rec,qy_rec);
		}
		if(fabs(theta_rec)<PI/3)
		{
			h_t_REC_wRES_cut_pi3->Fill(t_rec);
			h_t_REC_2d_wRES_cut_pi3->Fill(qx_rec,qy_rec);
		}
		if(fabs(theta_rec)<PI/4)
		{
    		h_t_REC_wRES_cut_pi4->Fill(t_rec);
			h_t_REC_2d_wRES_cut_pi4->Fill(qx_rec,qy_rec);
		}
		if(fabs(theta_rec)<PI/6)
		{
    		h_t_REC_wRES_cut_pi6->Fill(t_rec);
			h_t_REC_2d_wRES_cut_pi6->Fill(qx_rec,qy_rec);
		}
		if(fabs(theta_rec)<PI/9)
		{
    		h_t_REC_wRES_cut_pi9->Fill(t_rec);
			h_t_REC_2d_wRES_cut_pi9->Fill(qx_rec,qy_rec);
		}
		if(fabs(theta_rec)<PI/12)
		{
    		h_t_REC_wRES_cut_pi12->Fill(t_rec);
			h_t_REC_2d_wRES_cut_pi12->Fill(qx_rec,qy_rec);
		}
		if(fabs(theta_rec)<PI/16)
		{
    		h_t_REC_wRES_cut_pi16->Fill(t_rec);
			h_t_REC_2d_wRES_cut_pi16->Fill(qx_rec,qy_rec);
		}
        if(fabs(theta_rec)<PI/20)
		{
    		h_t_REC_wRES_cut_pi20->Fill(t_rec);
			h_t_REC_2d_wRES_cut_pi20->Fill(qx_rec,qy_rec);
		}
        if(fabs(theta_rec)<PI/24)
		{
    		h_t_REC_wRES_cut_pi24->Fill(t_rec);
			h_t_REC_2d_wRES_cut_pi24->Fill(qx_rec,qy_rec);
		}
    }

}

output->Write();
output->Close();

return 0;
}
