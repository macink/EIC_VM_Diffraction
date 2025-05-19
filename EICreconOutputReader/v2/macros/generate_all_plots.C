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
	TString angle_pi16="#pi/16";
	TString angle_pi20="#pi/20";
	TString angle_pi24="#pi/24";
	TString daug_label="K^{+}K^{-}";
	if(filename=="jpsi") {vm_label="J/#psi";daug_label="e^{+}e^{-}";}
	
	//t distribution
	TH1D* h_t_MC = (TH1D*) file->Get("h_t_MC");
	TH1D* h_t_REC = (TH1D*) file->Get("h_t_REC");
	TH1D* h_t_REC_new_method = (TH1D*) file->Get("h_t_REC_new_method");
	TH1D* h_t_REC_new_method_wZ = (TH1D*) file->Get("h_t_REC_new_method_wZ");
	TH1D* h_t_REC_new_method_wE = (TH1D*) file->Get("h_t_REC_new_method_wE");
	//t reco with resolution and angle cut
	TH1D* h_t_REC_wRES_cut_pi2 = (TH1D*) file->Get("h_t_REC_wRES_cut_pi2");
	TH1D* h_t_REC_wRES_cut_pi3 = (TH1D*) file->Get("h_t_REC_wRES_cut_pi3");
	TH1D* h_t_REC_wRES_cut_pi4 = (TH1D*) file->Get("h_t_REC_wRES_cut_pi4");
	TH1D* h_t_REC_wRES_cut_pi6 = (TH1D*) file->Get("h_t_REC_wRES_cut_pi6");
	TH1D* h_t_REC_wRES_cut_pi9 = (TH1D*) file->Get("h_t_REC_wRES_cut_pi9");
	TH1D* h_t_REC_wRES_cut_pi12 = (TH1D*) file->Get("h_t_REC_wRES_cut_pi12");
	TH1D* h_t_REC_wRES_cut_pi16 = (TH1D*) file->Get("h_t_REC_wRES_cut_pi16");
	TH1D* h_t_REC_wRES_cut_pi20 = (TH1D*) file->Get("h_t_REC_wRES_cut_pi20");
	TH1D* h_t_REC_wRES_cut_pi24 = (TH1D*) file->Get("h_t_REC_wRES_cut_pi24");
	//t reco with resolution and angle cut with qz absorbed into qx
	TH1D* h_t_REC_wRES_cut_pi2_wZ = (TH1D*) file->Get("h_t_REC_wRES_cut_pi2_wZ");
	TH1D* h_t_REC_wRES_cut_pi3_wZ = (TH1D*) file->Get("h_t_REC_wRES_cut_pi3_wZ");
	TH1D* h_t_REC_wRES_cut_pi4_wZ = (TH1D*) file->Get("h_t_REC_wRES_cut_pi4_wZ");
	TH1D* h_t_REC_wRES_cut_pi6_wZ = (TH1D*) file->Get("h_t_REC_wRES_cut_pi6_wZ");
	TH1D* h_t_REC_wRES_cut_pi9_wZ = (TH1D*) file->Get("h_t_REC_wRES_cut_pi9_wZ");
	TH1D* h_t_REC_wRES_cut_pi12_wZ = (TH1D*) file->Get("h_t_REC_wRES_cut_pi12_wZ");
	TH1D* h_t_REC_wRES_cut_pi16_wZ = (TH1D*) file->Get("h_t_REC_wRES_cut_pi16_wZ");
	TH1D* h_t_REC_wRES_cut_pi20_wZ = (TH1D*) file->Get("h_t_REC_wRES_cut_pi20_wZ");
	TH1D* h_t_REC_wRES_cut_pi24_wZ = (TH1D*) file->Get("h_t_REC_wRES_cut_pi24_wZ");
	//t reco with resolution and angle cut with E absorbed into qx
	TH1D* h_t_REC_wRES_cut_pi2_wE = (TH1D*) file->Get("h_t_REC_wRES_cut_pi2_wE");
	TH1D* h_t_REC_wRES_cut_pi3_wE = (TH1D*) file->Get("h_t_REC_wRES_cut_pi3_wE");
	TH1D* h_t_REC_wRES_cut_pi4_wE = (TH1D*) file->Get("h_t_REC_wRES_cut_pi4_wE");
	TH1D* h_t_REC_wRES_cut_pi6_wE = (TH1D*) file->Get("h_t_REC_wRES_cut_pi6_wE");
	TH1D* h_t_REC_wRES_cut_pi9_wE = (TH1D*) file->Get("h_t_REC_wRES_cut_pi9_wE");
	TH1D* h_t_REC_wRES_cut_pi12_wE = (TH1D*) file->Get("h_t_REC_wRES_cut_pi12_wE");
	TH1D* h_t_REC_wRES_cut_pi16_wE = (TH1D*) file->Get("h_t_REC_wRES_cut_pi16_wE");
	TH1D* h_t_REC_wRES_cut_pi20_wE = (TH1D*) file->Get("h_t_REC_wRES_cut_pi20_wE");
	TH1D* h_t_REC_wRES_cut_pi24_wE = (TH1D*) file->Get("h_t_REC_wRES_cut_pi24_wE");
	//t reco with only resolution no cut and with only cut no resolution (for testing)
	TH1D* h_t_REC_wRES = (TH1D*) file->Get("h_t_REC_wRES");
	TH1D* h_t_REC_wRES_wZ = (TH1D*) file->Get("h_t_REC_wRES_wZ");
	TH1D* h_t_REC_wRES_wE = (TH1D*) file->Get("h_t_REC_wRES_wE");
	TH1D* h_t_REC_wCUT = (TH1D*) file->Get("h_t_REC_wCUT"); // all angle cuts with no resolution will be the same
	TH1D* h_t_REC_wCUT_wZ = (TH1D*) file->Get("h_t_REC_wCUT_wZ");
	TH1D* h_t_REC_wCUT_wE = (TH1D*) file->Get("h_t_REC_wCUT_wE");
	//t reco 2d
	TH2D* h_t_REC_2d = (TH2D*) file->Get("h_t_REC_2d");
	TH2D* h_t_REC_2d_wZ = (TH2D*) file->Get("h_t_REC_2d_wZ");
	TH2D* h_t_REC_2d_wE = (TH2D*) file->Get("h_t_REC_2d_wE");
	// t reco 2d with resolution and angle cut
	TH2D* h_t_REC_2d_wRES_cut_pi2 = (TH2D*) file->Get("h_t_REC_2d_wRES_cut_pi2");
	TH2D* h_t_REC_2d_wRES_cut_pi3 = (TH2D*) file->Get("h_t_REC_2d_wRES_cut_pi3");
	TH2D* h_t_REC_2d_wRES_cut_pi4 = (TH2D*) file->Get("h_t_REC_2d_wRES_cut_pi4");
	TH2D* h_t_REC_2d_wRES_cut_pi6 = (TH2D*) file->Get("h_t_REC_2d_wRES_cut_pi6");
	TH2D* h_t_REC_2d_wRES_cut_pi9 = (TH2D*) file->Get("h_t_REC_2d_wRES_cut_pi9");
	TH2D* h_t_REC_2d_wRES_cut_pi12 = (TH2D*) file->Get("h_t_REC_2d_wRES_cut_pi12");
	TH2D* h_t_REC_2d_wRES_cut_pi16 = (TH2D*) file->Get("h_t_REC_2d_wRES_cut_pi16");
	TH2D* h_t_REC_2d_wRES_cut_pi20 = (TH2D*) file->Get("h_t_REC_2d_wRES_cut_pi20");
	TH2D* h_t_REC_2d_wRES_cut_pi24 = (TH2D*) file->Get("h_t_REC_2d_wRES_cut_pi24");
	// t reco 2d with resolution and angle cut with qz absorbed into qx
	TH2D* h_t_REC_2d_wRES_cut_pi2_wZ = (TH2D*) file->Get("h_t_REC_2d_wRES_cut_pi2_wZ");
	TH2D* h_t_REC_2d_wRES_cut_pi3_wZ = (TH2D*) file->Get("h_t_REC_2d_wRES_cut_pi3_wZ");
	TH2D* h_t_REC_2d_wRES_cut_pi4_wZ = (TH2D*) file->Get("h_t_REC_2d_wRES_cut_pi4_wZ");
	TH2D* h_t_REC_2d_wRES_cut_pi6_wZ = (TH2D*) file->Get("h_t_REC_2d_wRES_cut_pi6_wZ");
	TH2D* h_t_REC_2d_wRES_cut_pi9_wZ = (TH2D*) file->Get("h_t_REC_2d_wRES_cut_pi9_wZ");
	TH2D* h_t_REC_2d_wRES_cut_pi12_wZ = (TH2D*) file->Get("h_t_REC_2d_wRES_cut_pi12_wZ");
	TH2D* h_t_REC_2d_wRES_cut_pi16_wZ = (TH2D*) file->Get("h_t_REC_2d_wRES_cut_pi16_wZ");
	TH2D* h_t_REC_2d_wRES_cut_pi20_wZ = (TH2D*) file->Get("h_t_REC_2d_wRES_cut_pi20_wZ");
	TH2D* h_t_REC_2d_wRES_cut_pi24_wZ = (TH2D*) file->Get("h_t_REC_2d_wRES_cut_pi24_wZ");
	// t reco 2d with resolution and angle cut with E absorbed into qx
	TH2D* h_t_REC_2d_wRES_cut_pi2_wE = (TH2D*) file->Get("h_t_REC_2d_wRES_cut_pi2_wE");
	TH2D* h_t_REC_2d_wRES_cut_pi3_wE = (TH2D*) file->Get("h_t_REC_2d_wRES_cut_pi3_wE");
	TH2D* h_t_REC_2d_wRES_cut_pi4_wE = (TH2D*) file->Get("h_t_REC_2d_wRES_cut_pi4_wE");
	TH2D* h_t_REC_2d_wRES_cut_pi6_wE = (TH2D*) file->Get("h_t_REC_2d_wRES_cut_pi6_wE");
	TH2D* h_t_REC_2d_wRES_cut_pi9_wE = (TH2D*) file->Get("h_t_REC_2d_wRES_cut_pi9_wE");
	TH2D* h_t_REC_2d_wRES_cut_pi12_wE = (TH2D*) file->Get("h_t_REC_2d_wRES_cut_pi12_wE");
	TH2D* h_t_REC_2d_wRES_cut_pi16_wE = (TH2D*) file->Get("h_t_REC_2d_wRES_cut_pi16_wE");
	TH2D* h_t_REC_2d_wRES_cut_pi20_wE = (TH2D*) file->Get("h_t_REC_2d_wRES_cut_pi20_wE");
	TH2D* h_t_REC_2d_wRES_cut_pi24_wE = (TH2D*) file->Get("h_t_REC_2d_wRES_cut_pi24_wE");
	//t reco 2d with resolution only
	TH2D* h_t_REC_2d_wRES = (TH2D*) file->Get("h_t_REC_2d_wRES");
	TH2D* h_t_REC_2d_wRES_wZ = (TH2D*) file->Get("h_t_REC_2d_wRES_wZ");
	TH2D* h_t_REC_2d_wRES_wE = (TH2D*) file->Get("h_t_REC_2d_wRES_wE");
	//t reco 2d with angle cut only
	TH2D* h_t_REC_2d_wCUT_pi2 = (TH2D*) file->Get("h_t_REC_2d_wCUT_pi2");
	TH2D* h_t_REC_2d_wCUT_pi3 = (TH2D*) file->Get("h_t_REC_2d_wCUT_pi3");
	TH2D* h_t_REC_2d_wCUT_pi4 = (TH2D*) file->Get("h_t_REC_2d_wCUT_pi4");
	TH2D* h_t_REC_2d_wCUT_pi6 = (TH2D*) file->Get("h_t_REC_2d_wCUT_pi6");
	TH2D* h_t_REC_2d_wCUT_pi9 = (TH2D*) file->Get("h_t_REC_2d_wCUT_pi9");
	TH2D* h_t_REC_2d_wCUT_pi12 = (TH2D*) file->Get("h_t_REC_2d_wCUT_pi12");
	TH2D* h_t_REC_2d_wCUT_pi16 = (TH2D*) file->Get("h_t_REC_2d_wCUT_pi16");
	TH2D* h_t_REC_2d_wCUT_pi20 = (TH2D*) file->Get("h_t_REC_2d_wCUT_pi20");
	TH2D* h_t_REC_2d_wCUT_pi24 = (TH2D*) file->Get("h_t_REC_2d_wCUT_pi24");
	//t reco 2d with angle cut only with qz absorbed into qx
	TH2D* h_t_REC_2d_wCUT_pi2_wZ = (TH2D*) file->Get("h_t_REC_2d_wCUT_pi2_wZ");
	TH2D* h_t_REC_2d_wCUT_pi3_wZ = (TH2D*) file->Get("h_t_REC_2d_wCUT_pi3_wZ");
	TH2D* h_t_REC_2d_wCUT_pi4_wZ = (TH2D*) file->Get("h_t_REC_2d_wCUT_pi4_wZ");
	TH2D* h_t_REC_2d_wCUT_pi6_wZ = (TH2D*) file->Get("h_t_REC_2d_wCUT_pi6_wZ");
	TH2D* h_t_REC_2d_wCUT_pi9_wZ = (TH2D*) file->Get("h_t_REC_2d_wCUT_pi9_wZ");
	TH2D* h_t_REC_2d_wCUT_pi12_wZ = (TH2D*) file->Get("h_t_REC_2d_wCUT_pi12_wZ");
	TH2D* h_t_REC_2d_wCUT_pi16_wZ = (TH2D*) file->Get("h_t_REC_2d_wCUT_pi16_wZ");
	TH2D* h_t_REC_2d_wCUT_pi20_wZ = (TH2D*) file->Get("h_t_REC_2d_wCUT_pi20_wZ");
	TH2D* h_t_REC_2d_wCUT_pi24_wZ = (TH2D*) file->Get("h_t_REC_2d_wCUT_pi24_wZ");
	//t reco 2d with angle cut only with E absorbed into qx
	TH2D* h_t_REC_2d_wCUT_pi2_wE = (TH2D*) file->Get("h_t_REC_2d_wCUT_pi2_wE");
	TH2D* h_t_REC_2d_wCUT_pi3_wE = (TH2D*) file->Get("h_t_REC_2d_wCUT_pi3_wE");
	TH2D* h_t_REC_2d_wCUT_pi4_wE = (TH2D*) file->Get("h_t_REC_2d_wCUT_pi4_wE");
	TH2D* h_t_REC_2d_wCUT_pi6_wE = (TH2D*) file->Get("h_t_REC_2d_wCUT_pi6_wE");
	TH2D* h_t_REC_2d_wCUT_pi9_wE = (TH2D*) file->Get("h_t_REC_2d_wCUT_pi9_wE");
	TH2D* h_t_REC_2d_wCUT_pi12_wE = (TH2D*) file->Get("h_t_REC_2d_wCUT_pi12_wE");
	TH2D* h_t_REC_2d_wCUT_pi16_wE = (TH2D*) file->Get("h_t_REC_2d_wCUT_pi16_wE");
	TH2D* h_t_REC_2d_wCUT_pi20_wE = (TH2D*) file->Get("h_t_REC_2d_wCUT_pi20_wE");
	TH2D* h_t_REC_2d_wCUT_pi24_wE = (TH2D*) file->Get("h_t_REC_2d_wCUT_pi24_wE");

	//create canvas for each angle
	const int numCanvas = 9; 

	//angle labels
	TString angleLabels[numCanvas] = 
	{
    	angle_pi2, angle_pi3, angle_pi4, angle_pi6, angle_pi9, angle_pi12, angle_pi16, angle_pi20, angle_pi24
	};
	double angles[numCanvas] = 
	{
		M_PI/2, M_PI/3, M_PI/4, M_PI/6, M_PI/9, M_PI/12, M_PI/16, M_PI/20, M_PI/24
	};
	
	//histograms
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
	TH1D* h_t_REC_wRES_cut_wZ[numCanvas] = 
	{
		h_t_REC_wRES_cut_pi2_wZ,
		h_t_REC_wRES_cut_pi3_wZ,
		h_t_REC_wRES_cut_pi4_wZ,
		h_t_REC_wRES_cut_pi6_wZ,
		h_t_REC_wRES_cut_pi9_wZ,
		h_t_REC_wRES_cut_pi12_wZ,
		h_t_REC_wRES_cut_pi16_wZ,
		h_t_REC_wRES_cut_pi20_wZ,
		h_t_REC_wRES_cut_pi24_wZ
	};
	TH1D* h_t_REC_wRES_cut_wE[numCanvas] = 
	{
		h_t_REC_wRES_cut_pi2_wE,
		h_t_REC_wRES_cut_pi3_wE,
		h_t_REC_wRES_cut_pi4_wE,
		h_t_REC_wRES_cut_pi6_wE,
		h_t_REC_wRES_cut_pi9_wE,
		h_t_REC_wRES_cut_pi12_wE,
		h_t_REC_wRES_cut_pi16_wE,
		h_t_REC_wRES_cut_pi20_wE,
		h_t_REC_wRES_cut_pi24_wE
	};
	TH2D* h_t_REC_2d_wRES_cut[numCanvas] = 
	{
		h_t_REC_2d_wRES_cut_pi2,
		h_t_REC_2d_wRES_cut_pi3,
		h_t_REC_2d_wRES_cut_pi4,
		h_t_REC_2d_wRES_cut_pi6,
		h_t_REC_2d_wRES_cut_pi9,
		h_t_REC_2d_wRES_cut_pi12,
		h_t_REC_2d_wRES_cut_pi16,
		h_t_REC_2d_wRES_cut_pi20,
		h_t_REC_2d_wRES_cut_pi24
	};
	TH2D* h_t_REC_2d_wRES_cut_wZ[numCanvas] = 
	{
		h_t_REC_2d_wRES_cut_pi2_wZ,
		h_t_REC_2d_wRES_cut_pi3_wZ,
		h_t_REC_2d_wRES_cut_pi4_wZ,
		h_t_REC_2d_wRES_cut_pi6_wZ,
		h_t_REC_2d_wRES_cut_pi9_wZ,
		h_t_REC_2d_wRES_cut_pi12_wZ,
		h_t_REC_2d_wRES_cut_pi16_wZ,
		h_t_REC_2d_wRES_cut_pi20_wZ,
		h_t_REC_2d_wRES_cut_pi24_wZ
	};
	TH2D* h_t_REC_2d_wRES_cut_wE[numCanvas] = 
	{
		h_t_REC_2d_wRES_cut_pi2_wE,
		h_t_REC_2d_wRES_cut_pi3_wE,
		h_t_REC_2d_wRES_cut_pi4_wE,
		h_t_REC_2d_wRES_cut_pi6_wE,
		h_t_REC_2d_wRES_cut_pi9_wE,
		h_t_REC_2d_wRES_cut_pi12_wE,
		h_t_REC_2d_wRES_cut_pi16_wE,
		h_t_REC_2d_wRES_cut_pi20_wE,
		h_t_REC_2d_wRES_cut_pi24_wE
	};
	TH2D* h_t_REC_2d_wCUT[numCanvas] = 
	{
		h_t_REC_2d_wCUT_pi2,
		h_t_REC_2d_wCUT_pi3,
		h_t_REC_2d_wCUT_pi4,
		h_t_REC_2d_wCUT_pi6,
		h_t_REC_2d_wCUT_pi9,
		h_t_REC_2d_wCUT_pi12,
		h_t_REC_2d_wCUT_pi16,
		h_t_REC_2d_wCUT_pi20,
		h_t_REC_2d_wCUT_pi24
	};
	TH2D* h_t_REC_2d_wCUT_wZ[numCanvas] = 
	{
		h_t_REC_2d_wCUT_pi2_wZ,
		h_t_REC_2d_wCUT_pi3_wZ,
		h_t_REC_2d_wCUT_pi4_wZ,
		h_t_REC_2d_wCUT_pi6_wZ,
		h_t_REC_2d_wCUT_pi9_wZ,
		h_t_REC_2d_wCUT_pi12_wZ,
		h_t_REC_2d_wCUT_pi16_wZ,
		h_t_REC_2d_wCUT_pi20_wZ,
		h_t_REC_2d_wCUT_pi24_wZ
	};
	TH2D* h_t_REC_2d_wCUT_wE[numCanvas] = 
	{
		h_t_REC_2d_wCUT_pi2_wE,
		h_t_REC_2d_wCUT_pi3_wE,
		h_t_REC_2d_wCUT_pi4_wE,
		h_t_REC_2d_wCUT_pi6_wE,
		h_t_REC_2d_wCUT_pi9_wE,
		h_t_REC_2d_wCUT_pi12_wE,
		h_t_REC_2d_wCUT_pi16_wE,
		h_t_REC_2d_wCUT_pi20_wE,
		h_t_REC_2d_wCUT_pi24_wE
	};

/*	//reproduce truth plot with our decomposition, no norm
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
		double integral_REC_new_method_wZ = h_t_REC_new_method_wZ->Integral();
		double integral_REC_new_method_wE = h_t_REC_new_method_wE->Integral();
		// Draw histograms 
		h_t_MC->Draw("same");

		h_t_REC->SetMarkerStyle(20);
		h_t_REC->Draw("PEsame");

		h_t_REC_new_method->SetMarkerStyle(30);
		h_t_REC_new_method->SetMarkerColor(kRed);
		h_t_REC_new_method->Draw("P same");

		h_t_REC_new_method_wZ->SetMarkerStyle(25);
		h_t_REC_new_method_wZ->SetMarkerColor(kGreen+1);
		//h_t_REC_new_method_wZ->Draw("P same");

		h_t_REC_wRES_cut_wE[i]->SetMarkerStyle(5);
		h_t_REC_wRES_cut_wE[i]->SetMarkerColor(kOrange-3);
		//h_t_REC_wRES_cut_wE[i]->Draw("P same");

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
		TLegend* w7 = new TLegend(0.26, 0.64, 0.93, 0.76);
		w7->SetLineColor(kWhite);
		w7->SetFillColor(0);
		w7->SetTextSize(17);
		w7->SetTextFont(45);
		w7->AddEntry(h_t_MC, Form("Sartre %s MC %.f evts", vm_label.Data(),integral_MC), "L");
		w7->AddEntry(h_t_REC, Form("Sartre %s RECO w. EEMC %.f evts", vm_label.Data(),integral_REC), "P");
		w7->AddEntry(h_t_REC_new_method, Form("Sartre %s RECO' new meth. %.f evts", vm_label.Data(),integral_REC_new_method), "P");
		//w7->AddEntry(h_t_REC_new_method_wZ, Form("Sartre %s RECO' new meth. w. z %.f evts", vm_label.Data(),integral_REC_new_method_wZ), "P");
		//w7->AddEntry(h_t_REC_new_method_wE, Form("Sartre %s RECO' new meth. w. E %.f evts", vm_label.Data(),integral_REC_new_method_wE), "P");
		w7->Draw("same");
		// Save figure
		newMethod_canvases[i]->Print("./figures/noNorm_new_method.pdf");
	}
*/

	//reproduce truth plot with our decomposition, total norm
	TCanvas* totalNorm_newMethod_canvases[numCanvas];
	for (int i=0; i<1; i++) 
	{
		TString canvasName = Form("totalNorm_newMethod_c%d", i+1);
		totalNorm_newMethod_canvases[i] = new TCanvas(canvasName, canvasName, 1, 1, 600, 600);
		gPad->SetLogy(1);
		gPad->SetTicks();
		gPad->SetLeftMargin(0.15);
		gPad->SetBottomMargin(0.15);
		gPad->SetRightMargin(0.01);
		// Generate histogram for each canvas
		TString totalNorm_newMethod_histName = Form("totalNorm_newMethod_base%d", i+1);
		TH1D* totalNorm_newMethod_baseHist = makeHist(totalNorm_newMethod_histName, "", "|#it{t} | (GeV^{2})", "dN/d|#it{t} | (GeV^{-2}) ", 100, 0, 0.18, kBlack);
		totalNorm_newMethod_baseHist->GetYaxis()->SetRangeUser(8e-2, 8e8);
		totalNorm_newMethod_baseHist->GetXaxis()->SetTitleColor(kBlack);
		fixedFontHist1D(totalNorm_newMethod_baseHist, 1., 1.2);
		totalNorm_newMethod_baseHist->GetYaxis()->SetTitleSize(totalNorm_newMethod_baseHist->GetYaxis()->GetTitleSize()*1.5);
		totalNorm_newMethod_baseHist->GetXaxis()->SetTitleSize(totalNorm_newMethod_baseHist->GetXaxis()->GetTitleSize()*1.5);
		totalNorm_newMethod_baseHist->GetYaxis()->SetLabelSize(totalNorm_newMethod_baseHist->GetYaxis()->GetLabelSize()*1.5);
		totalNorm_newMethod_baseHist->GetXaxis()->SetLabelSize(totalNorm_newMethod_baseHist->GetXaxis()->GetLabelSize()*1.5);
		totalNorm_newMethod_baseHist->GetXaxis()->SetNdivisions(4,4,0);
		totalNorm_newMethod_baseHist->GetYaxis()->SetNdivisions(5,5,0);
		totalNorm_newMethod_baseHist->Draw();
		// get counts
		double integral_MC = h_t_MC->Integral();
		double integral_REC = h_t_REC->Integral();
		double integral_REC_new_method = h_t_REC_new_method->Integral();
		double integral_REC_new_method_wZ = h_t_REC_new_method_wZ->Integral();
		double integral_REC_new_method_wE = h_t_REC_new_method_wE->Integral();
		// normalize to total events
		if(integral_MC>0 && integral_REC_new_method>0 && integral_REC_new_method_wZ>0 && integral_REC_new_method_wE>0) 
		{
    		h_t_REC_new_method->Scale(integral_MC/integral_REC_new_method);
			h_t_REC_new_method_wZ->Scale(integral_MC/integral_REC_new_method_wZ);
			h_t_REC_new_method_wE->Scale(integral_MC/integral_REC_new_method_wE);
		}

		// Draw histograms 
		h_t_MC->Draw("same");

		h_t_REC->SetMarkerStyle(20);
		h_t_REC->Draw("PEsame");

		h_t_REC_new_method->SetMarkerStyle(30);
		h_t_REC_new_method->SetMarkerColor(kRed);
		h_t_REC_new_method->Draw("P same");

		h_t_REC_new_method_wZ->SetMarkerStyle(25);
		h_t_REC_new_method_wZ->SetMarkerColor(kGreen+1);
		//h_t_REC_new_method_wZ->Draw("P same");

		h_t_REC_wRES_cut_wE[i]->SetMarkerStyle(5);
		h_t_REC_wRES_cut_wE[i]->SetMarkerColor(kOrange-3);
		//h_t_REC_wRES_cut_wE[i]->Draw("P same");

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
		TLegend* w7 = new TLegend(0.26, 0.64, 0.93, 0.76);
		w7->SetLineColor(kWhite);
		w7->SetFillColor(0);
		w7->SetTextSize(17);
		w7->SetTextFont(45);
		w7->AddEntry(h_t_MC, Form("Sartre %s MC %.f evts", vm_label.Data(),integral_MC), "L");
		w7->AddEntry(h_t_REC, Form("Sartre %s RECO w. EEMC %.f evts", vm_label.Data(),integral_REC), "P");
		w7->AddEntry(h_t_REC_new_method, Form("Sartre %s RECO' new meth. %.f evts", vm_label.Data(),integral_REC_new_method), "P");
		//w7->AddEntry(h_t_REC_new_method_wZ, Form("Sartre %s RECO' new meth. w. z %.f evts", vm_label.Data(),integral_REC_new_method_wZ), "P");
		//w7->AddEntry(h_t_REC_new_method_wE, Form("Sartre %s RECO' new meth. w. E %.f evts", vm_label.Data(),integral_REC_new_method_wE), "P");
		w7->Draw("same");

		TLegend *w8 = new TLegend(0.48,0.61,0.93,0.56);
		w8->SetLineColor(kWhite);
		w8->SetFillColor(0);
		w8->SetTextSize(17);
		w8->SetTextFont(45);
		w8->AddEntry(h_t_REC_new_method,"normalization: #int|#it{t}|_{MC}/#int|#it{t}|_{RECO'}","P");
		//w8->AddEntry(h_t_REC_new_method_wZ,"","P");
		//w8->AddEntry(h_t_REC_new_method_wE,"","P");
		w8->Draw("same");
		// Save figure
		totalNorm_newMethod_canvases[i]->Print("./figures/totalNorm_new_method.pdf");
	}


/*  //reproduce truth plot with our decomposition, L norm
	TCanvas* Lnorm_newMethod_canvases[numCanvas];
	for (int i=0; i<1; i++) 
	{
		TString canvasName = Form("Lnorm_newMethod_c%d", i+1);
		Lnorm_newMethod_canvases[i] = new TCanvas(canvasName, canvasName, 1, 1, 600, 600);
		gPad->SetLogy(1);
		gPad->SetTicks();
		gPad->SetLeftMargin(0.15);
		gPad->SetBottomMargin(0.15);
		gPad->SetRightMargin(0.01);
		// Generate histogram for each canvas
		TString Lnorm_newMethod_histName = Form("Lnorm_newMethod_base%d", i+1);
		TH1D* Lnorm_newMethod_baseHist = makeHist(Lnorm_newMethod_histName, "", "|#it{t} | (GeV^{2})", "dN/d|#it{t} | (GeV^{-2}) ", 100, 0, 0.18, kBlack);
		Lnorm_newMethod_baseHist->GetYaxis()->SetRangeUser(8e-2, 8e8);
		Lnorm_newMethod_baseHist->GetXaxis()->SetTitleColor(kBlack);
		fixedFontHist1D(Lnorm_newMethod_baseHist, 1., 1.2);
		Lnorm_newMethod_baseHist->GetYaxis()->SetTitleSize(Lnorm_newMethod_baseHist->GetYaxis()->GetTitleSize()*1.5);
		Lnorm_newMethod_baseHist->GetXaxis()->SetTitleSize(Lnorm_newMethod_baseHist->GetXaxis()->GetTitleSize()*1.5);
		Lnorm_newMethod_baseHist->GetYaxis()->SetLabelSize(Lnorm_newMethod_baseHist->GetYaxis()->GetLabelSize()*1.5);
		Lnorm_newMethod_baseHist->GetXaxis()->SetLabelSize(Lnorm_newMethod_baseHist->GetXaxis()->GetLabelSize()*1.5);
		Lnorm_newMethod_baseHist->GetXaxis()->SetNdivisions(4,4,0);
		Lnorm_newMethod_baseHist->GetYaxis()->SetNdivisions(5,5,0);
		Lnorm_newMethod_baseHist->Draw();
		// get counts
		double integral_MC = h_t_MC->Integral();
		double integral_REC = h_t_REC->Integral();
		double integral_REC_new_method = h_t_REC_new_method->Integral();
		double integral_REC_new_method_wZ = h_t_REC_new_method_wZ->Integral();
		double integral_REC_new_method_wE = h_t_REC_new_method_wE->Integral();
		//normalize
		h_t_REC_new_method->Scale((M_PI/2)/(M_PI/2));
		h_t_REC_new_method_wZ->Scale((M_PI/2)/(M_PI/2));
		h_t_REC_new_method_wE->Scale((M_PI/2)/(M_PI/2));

		// Draw histograms 
		h_t_MC->Draw("same");

		h_t_REC->SetMarkerStyle(20);
		h_t_REC->Draw("PEsame");

		h_t_REC_new_method->SetMarkerStyle(30);
		h_t_REC_new_method->SetMarkerColor(kRed);
		h_t_REC_new_method->Draw("P same");

		h_t_REC_new_method_wZ->SetMarkerStyle(25);
		h_t_REC_new_method_wZ->SetMarkerColor(kGreen+1);
		//h_t_REC_new_method_wZ->Draw("P same");

		h_t_REC_wRES_cut_wE[i]->SetMarkerStyle(5);
		h_t_REC_wRES_cut_wE[i]->SetMarkerColor(kOrange-3);
		//h_t_REC_wRES_cut_wE[i]->Draw("P same");

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
		TLegend* w7 = new TLegend(0.26, 0.64, 0.93, 0.76);
		w7->SetLineColor(kWhite);
		w7->SetFillColor(0);
		w7->SetTextSize(17);
		w7->SetTextFont(45);
		w7->AddEntry(h_t_MC, Form("Sartre %s MC %.f evts", vm_label.Data(),integral_MC), "L");
		w7->AddEntry(h_t_REC, Form("Sartre %s RECO w. EEMC %.f evts", vm_label.Data(),integral_REC), "P");
		w7->AddEntry(h_t_REC_new_method, Form("Sartre %s RECO' new meth. %.f evts", vm_label.Data(),integral_REC_new_method), "P");
		//w7->AddEntry(h_t_REC_new_method_wZ, Form("Sartre %s RECO' new meth. w. z %.f evts", vm_label.Data(),integral_REC_new_method_wZ), "P");
		//w7->AddEntry(h_t_REC_new_method_wE, Form("Sartre %s RECO' new meth. w. E %.f evts", vm_label.Data(),integral_REC_new_method_wE), "P");
		w7->Draw("same");

		TLegend *w8 = new TLegend(0.48,0.61,0.93,0.56);
		w8->SetLineColor(kWhite);
		w8->SetFillColor(0);
		w8->SetTextSize(17);
		w8->SetTextFont(45);
		w8->AddEntry(h_t_REC_new_method,"normalization: #frac{#pi}{2} / #frac{#pi}{2}","P");
		//w8->AddEntry(h_t_REC_new_method_wZ,"","P");
		//w8->AddEntry(h_t_REC_new_method_wE,"","P");
		w8->Draw("same");
		// Save figure
		Lnorm_newMethod_canvases[i]->Print("./figures/Lnorm_new_method.pdf");
	}
*/

/*	//plots with resolution and cut, L norm
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
		// get counts
		double integral_MC = h_t_MC->Integral();
		double integral_REC = h_t_REC->Integral();
		double integral_REC_wRES_cut[numCanvas];
		integral_REC_wRES_cut[i] = h_t_REC_wRES_cut[i]->Integral();
		double integral_REC_wRES_cut_wZ[numCanvas];
		integral_REC_wRES_cut_wZ[i] = h_t_REC_wRES_cut_wZ[i]->Integral();
		double integral_REC_wRES_cut_wE[numCanvas];
		integral_REC_wRES_cut_wE[i] = h_t_REC_wRES_cut_wE[i]->Integral();
		// normalize
		h_t_REC_wRES_cut[i]->Scale((M_PI/2)/(angles[i]));
		h_t_REC_wRES_cut_wZ[i]->Scale((M_PI/2)/(angles[i]));
		h_t_REC_wRES_cut_wE[i]->Scale((M_PI/2)/(angles[i]));
		
		// Draw histograms with different angle cuts
		h_t_MC->Draw("same");

		h_t_REC->SetMarkerStyle(20);
		h_t_REC->Draw("PEsame");

		h_t_REC_wRES_cut[i]->SetMarkerStyle(30);
		h_t_REC_wRES_cut[i]->SetMarkerColor(kRed);
		h_t_REC_wRES_cut[i]->Draw("P same");

		h_t_REC_wRES_cut_wZ[i]->SetMarkerStyle(25);
		h_t_REC_wRES_cut_wZ[i]->SetMarkerColor(kGreen+1);
		//h_t_REC_wRES_cut_wZ[i]->Draw("P same");

		h_t_REC_wRES_cut_wE[i]->SetMarkerStyle(5);
		h_t_REC_wRES_cut_wE[i]->SetMarkerColor(kOrange-3);
		//h_t_REC_wRES_cut_wE[i]->Draw("P same");

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
		TLegend* w7 = new TLegend(0.26, 0.64, 0.93, 0.76);
		w7->SetLineColor(kWhite);
		w7->SetFillColor(0);
		w7->SetTextSize(17);
		w7->SetTextFont(45);
		w7->AddEntry(h_t_MC, Form("Sartre %s MC %.f evts", vm_label.Data(),integral_MC), "L");
		w7->AddEntry(h_t_REC, Form("Sartre %s RECO w. EEMC %.f evts", vm_label.Data(),integral_REC), "P");
		w7->AddEntry(h_t_REC_wRES_cut[i], Form("Sartre %s RECO' #theta = %s %.f evts", vm_label.Data(), angleLabels[i].Data(),integral_REC_wRES_cut[i]), "P");
		//w7->AddEntry(h_t_REC_wRES_cut_wZ[i], Form("Sartre %s RECO' w. z #theta = %s: %.f evts", vm_label.Data(), angleLabels[i].Data(),integral_REC_wRES_cut_wZ[i]), "P");
		//w7->AddEntry(h_t_REC_wRES_cut_wE[i], Form("Sartre %s RECO' w. E #theta = %s: %.f evts", vm_label.Data(), angleLabels[i].Data(),integral_REC_wRES_cut_wE[i]), "P");
		w7->Draw("same");

		TLegend *w8 = new TLegend(0.48,0.61,0.93,0.56);
		w8->SetLineColor(kWhite);
		w8->SetFillColor(0);
		w8->SetTextSize(17);
		w8->SetTextFont(45);
		w8->AddEntry(h_t_REC_wRES_cut[i],"normalization: #frac{#pi}{2} / #theta","P");
		//w8->AddEntry(h_t_REC_wRES_cut_wZ[i],"","P");
		//w8->AddEntry(h_t_REC_wRES_cut_wE[i],"","P");
		w8->Draw("same");
	
		// Save figure
		TString cleanLabel = angleLabels[i];
		cleanLabel.ReplaceAll("#", "");
		cleanLabel.ReplaceAll("/", "");
		Lnorm_canvases_resCut[i]->Print(Form("./figures/Lnorm_angle%s_wRES_cut.pdf",cleanLabel.Data()));
	}
*/

	//plots with resolution and cut, total norm
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
		// get counts
		double integral_MC = h_t_MC->Integral();
		double integral_REC = h_t_REC->Integral();
		double integral_REC_wRES_cut[numCanvas];
		integral_REC_wRES_cut[i] = h_t_REC_wRES_cut[i]->Integral();
		double integral_REC_wRES_cut_wZ[numCanvas];
		integral_REC_wRES_cut_wZ[i] = h_t_REC_wRES_cut_wZ[i]->Integral();
		double integral_REC_wRES_cut_wE[numCanvas];
		integral_REC_wRES_cut_wE[i] = h_t_REC_wRES_cut_wE[i]->Integral();
		// normalize
		if(integral_MC>0 && integral_REC_wRES_cut[i]>0 && integral_REC_wRES_cut_wZ[i]>0 && integral_REC_wRES_cut_wE[i]>0) 
		{
    		h_t_REC_wRES_cut[i]->Scale(integral_MC/integral_REC_wRES_cut[i]);
			h_t_REC_wRES_cut_wZ[i]->Scale(integral_MC/integral_REC_wRES_cut_wZ[i]);
			h_t_REC_wRES_cut_wE[i]->Scale(integral_MC/integral_REC_wRES_cut_wE[i]);
		}
			
		// Draw histograms with different angle cuts
		h_t_MC->Draw("same");

		h_t_REC->SetMarkerStyle(20);
		h_t_REC->Draw("PEsame");

		h_t_REC_wRES_cut[i]->SetMarkerStyle(30);
		h_t_REC_wRES_cut[i]->SetMarkerColor(kRed);
		h_t_REC_wRES_cut[i]->Draw("P same");

		h_t_REC_wRES_cut_wZ[i]->SetMarkerStyle(25);
		h_t_REC_wRES_cut_wZ[i]->SetMarkerColor(kGreen+1);
		//h_t_REC_wRES_cut_wZ[i]->Draw("P same");

		h_t_REC_wRES_cut_wE[i]->SetMarkerStyle(5);
		h_t_REC_wRES_cut_wE[i]->SetMarkerColor(kOrange-3);
		//h_t_REC_wRES_cut_wE[i]->Draw("P same");

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
		TLegend* w7 = new TLegend(0.26, 0.64, 0.93, 0.76);
		w7->SetLineColor(kWhite);
		w7->SetFillColor(0);
		w7->SetTextSize(17);
		w7->SetTextFont(45);
		w7->AddEntry(h_t_MC, Form("Sartre %s MC %.f evts", vm_label.Data(),integral_MC), "L");
		w7->AddEntry(h_t_REC, Form("Sartre %s RECO w. EEMC %.f evts", vm_label.Data(),integral_REC), "P");
		w7->AddEntry(h_t_REC_wRES_cut[i], Form("Sartre %s RECO' #theta = %s %.f evts", vm_label.Data(), angleLabels[i].Data(),integral_REC_wRES_cut[i]), "P");
		//w7->AddEntry(h_t_REC_wRES_cut_wZ[i], Form("Sartre %s RECO' w. z #theta = %s %.f evts", vm_label.Data(), angleLabels[i].Data(),integral_REC_wRES_cut_wZ[i]), "P");
		//w7->AddEntry(h_t_REC_wRES_cut_wE[i], Form("Sartre %s RECO' w. E #theta = %s %.f evts", vm_label.Data(), angleLabels[i].Data(),integral_REC_wRES_cut_wE[i]), "P");
		w7->Draw("same");

		TLegend *w8 = new TLegend(0.48,0.61,0.93,0.56);
		w8->SetLineColor(kWhite);
		w8->SetFillColor(0);
		w8->SetTextSize(17);
		w8->SetTextFont(45);
		w8->AddEntry(h_t_REC_wRES_cut[i],"normalization: #int|#it{t}|_{MC}/#int|#it{t}|_{RECO'}","P");
		//w8->AddEntry(h_t_REC_wRES_cut_wZ[i],"","P");
		//w8->AddEntry(h_t_REC_wRES_cut_wE[i],"","P");
		w8->Draw("same");
	
		// Save figure
		TString cleanLabel = angleLabels[i];
		cleanLabel.ReplaceAll("#", "");
		cleanLabel.ReplaceAll("/", "");
		totalNorm_canvases_resCut[i]->Print(Form("./figures/totalNorm_angle%s_wRES_cut.pdf",cleanLabel.Data()));
	}

	
/*	//plots with resolution and cut, no norm
	TCanvas* noNorm_canvases_resCut[numCanvas];
	for (int i=0; i<numCanvas; i++) 
	{
    	TString canvasName = Form("noNorm_resCut_c%d", i+1);
    	noNorm_canvases_resCut[i] = new TCanvas(canvasName, canvasName, 1, 1, 600, 600);
    	gPad->SetLogy(1);
    	gPad->SetTicks();
    	gPad->SetLeftMargin(0.15);
    	gPad->SetBottomMargin(0.15);
    	gPad->SetRightMargin(0.01);
		// Generate histogram for each canvas
    	TString noNorm_histName_resCut = Form("noNorm_resCut_base%d", i+1);
    	TH1D* noNorm_baseHist_resCut = makeHist(noNorm_histName_resCut, "", "|#it{t} | (GeV^{2})", "dN/d|#it{t} | (GeV^{-2}) ", 100, 0, 0.18, kBlack);
    	noNorm_baseHist_resCut->GetYaxis()->SetRangeUser(8e-2, 8e8);
    	noNorm_baseHist_resCut->GetXaxis()->SetTitleColor(kBlack);
    	fixedFontHist1D(noNorm_baseHist_resCut, 1., 1.2);
    	noNorm_baseHist_resCut->GetYaxis()->SetTitleSize(noNorm_baseHist_resCut->GetYaxis()->GetTitleSize()*1.5);
    	noNorm_baseHist_resCut->GetXaxis()->SetTitleSize(noNorm_baseHist_resCut->GetXaxis()->GetTitleSize()*1.5);
    	noNorm_baseHist_resCut->GetYaxis()->SetLabelSize(noNorm_baseHist_resCut->GetYaxis()->GetLabelSize()*1.5);
    	noNorm_baseHist_resCut->GetXaxis()->SetLabelSize(noNorm_baseHist_resCut->GetXaxis()->GetLabelSize()*1.5);
    	noNorm_baseHist_resCut->GetXaxis()->SetNdivisions(4,4,0);
    	noNorm_baseHist_resCut->GetYaxis()->SetNdivisions(5,5,0);
    	noNorm_baseHist_resCut->Draw();
		// get counts
		double integral_MC = h_t_MC->Integral();
		double integral_REC = h_t_REC->Integral();
		double integral_REC_wRES_cut[numCanvas];
		integral_REC_wRES_cut[i] = h_t_REC_wRES_cut[i]->Integral();
		double integral_REC_wRES_cut_wZ[numCanvas];
		integral_REC_wRES_cut_wZ[i] = h_t_REC_wRES_cut_wZ[i]->Integral();
		double integral_REC_wRES_cut_wE[numCanvas];
		integral_REC_wRES_cut_wE[i] = h_t_REC_wRES_cut_wE[i]->Integral();

		// Draw histograms with different angle cuts
		h_t_MC->Draw("same");

		h_t_REC->SetMarkerStyle(20);
		h_t_REC->Draw("PEsame");

		h_t_REC_wRES_cut[i]->SetMarkerStyle(30);
		h_t_REC_wRES_cut[i]->SetMarkerColor(kRed);
		h_t_REC_wRES_cut[i]->Draw("P same");

		h_t_REC_wRES_cut_wZ[i]->SetMarkerStyle(25);
		h_t_REC_wRES_cut_wZ[i]->SetMarkerColor(kGreen+1);
		//h_t_REC_wRES_cut_wZ[i]->Draw("P same");

		h_t_REC_wRES_cut_wE[i]->SetMarkerStyle(5);
		h_t_REC_wRES_cut_wE[i]->SetMarkerColor(kOrange-3);
		//h_t_REC_wRES_cut_wE[i]->Draw("P same");

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
		TLegend* w7 = new TLegend(0.26, 0.64, 0.93, 0.76);
		w7->SetLineColor(kWhite);
		w7->SetFillColor(0);
		w7->SetTextSize(17);
		w7->SetTextFont(45);
		w7->AddEntry(h_t_MC, Form("Sartre %s MC %.f evts", vm_label.Data(),integral_MC), "L");
		w7->AddEntry(h_t_REC, Form("Sartre %s RECO w. EEMC %.f evts", vm_label.Data(),integral_REC), "P");
		w7->AddEntry(h_t_REC_wRES_cut[i], Form("Sartre %s RECO' #theta = %s, %.f evts", vm_label.Data(), angleLabels[i].Data(),integral_REC_wRES_cut[i]), "P");
		//w7->AddEntry(h_t_REC_wRES_cut_wZ[i], Form("Sartre %s RECO' w. z #theta = %s, %.f evts", vm_label.Data(), angleLabels[i].Data(),integral_REC_wRES_cut_wZ[i]), "P");
		//w7->AddEntry(h_t_REC_wRES_cut_wE[i], Form("Sartre %s RECO' w. E #theta = %s, %.f evts", vm_label.Data(), angleLabels[i].Data(),integral_REC_wRES_cut_wE[i]), "P");
		w7->Draw("same");
		// Save figure
		TString cleanLabel = angleLabels[i];
		cleanLabel.ReplaceAll("#", "");
		cleanLabel.ReplaceAll("/", "");
		noNorm_canvases_resCut[i]->Print(Form("./figures/noNorm_angle%s_wRES_cut.pdf",cleanLabel.Data()));
	}
*/

/*	//plot with resolution, L norm
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
		// get counts
		double integral_MC = h_t_MC->Integral();
		double integral_REC = h_t_REC->Integral();
		double integral_REC_wRES = h_t_REC_wRES->Integral();
		double integral_REC_wRES_wZ = h_t_REC_wRES_wZ->Integral();
		double integral_REC_wRES_wE = h_t_REC_wRES_wE->Integral();
		// normalize 
		h_t_REC_wRES->Scale((M_PI/2)/(M_PI/2));
		h_t_REC_wRES_wZ->Scale((M_PI/2)/(M_PI/2));
		h_t_REC_wRES_wE->Scale((M_PI/2)/(M_PI/2));

		// Draw histograms
		h_t_MC->Draw("same");

		h_t_REC->SetMarkerStyle(20);
		h_t_REC->Draw("PEsame");

		h_t_REC_wRES->SetMarkerStyle(30);
		h_t_REC_wRES->SetMarkerColor(kRed);
		h_t_REC_wRES->Draw("P same");

		h_t_REC_wRES_wZ->SetMarkerStyle(25);
		h_t_REC_wRES_wZ->SetMarkerColor(kGreen+1);
		//h_t_REC_wRES_wZ->Draw("P same");

		h_t_REC_wRES_wE->SetMarkerStyle(5);
		h_t_REC_wRES_wE->SetMarkerColor(kOrange-3);
		//h_t_REC_wRES_wE->Draw("P same");

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
		TLegend* w7 = new TLegend(0.26, 0.64, 0.93, 0.76);
		w7->SetLineColor(kWhite);
		w7->SetFillColor(0);
		w7->SetTextSize(17);
		w7->SetTextFont(45);
		w7->AddEntry(h_t_MC, Form("Sartre %s MC %.f evts", vm_label.Data(),integral_MC), "L");
		w7->AddEntry(h_t_REC, Form("Sartre %s RECO w. EEMC %.f evts", vm_label.Data(),integral_REC), "P");
		w7->AddEntry(h_t_REC_wRES, Form("Sartre %s RECO' new meth. %.f evts", vm_label.Data(),integral_REC_wRES), "P");
		//w7->AddEntry(h_t_REC_wRES_wZ, Form("Sartre %s RECO' new meth. w. z %.f evts", vm_label.Data(),integral_REC_wRES_wZ), "P");
		//w7->AddEntry(h_t_REC_wRES_wE, Form("Sartre %s RECO' new meth. w. E %.f evts", vm_label.Data(),integral_REC_wRES_wE), "P");
		w7->Draw("same");

		TLegend *w8 = new TLegend(0.48,0.61,0.93,0.56);
		w8->SetLineColor(kWhite);
		w8->SetFillColor(0);
		w8->SetTextSize(17);
		w8->SetTextFont(45);
		w8->AddEntry(h_t_REC_wRES,"normalization: #frac{#pi}{2} / #frac{#pi}{2}","P");
		//w8->AddEntry(h_t_REC_wRES_wZ,"","P");
		//w8->AddEntry(h_t_REC_wRES_wE,"","P");
		w8->Draw("same");

		//save figure
		Lnorm_canvases_res[i]->Print("./figures/Lnorm_wRES.pdf");
	}
*/

	//plot with resolution, total norm
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
		// get counts
		double integral_MC = h_t_MC->Integral();
		double integral_REC = h_t_REC->Integral();
		double integral_REC_wRES = h_t_REC_wRES->Integral();
		double integral_REC_wRES_wZ = h_t_REC_wRES_wZ->Integral();
		double integral_REC_wRES_wE = h_t_REC_wRES_wE->Integral();
		// normalize
		if(integral_MC>0 && integral_REC_wRES>0 && integral_REC_wRES_wZ>0 && integral_REC_wRES_wE>0) 
		{
    		h_t_REC_wRES->Scale(integral_MC/integral_REC_wRES);
			h_t_REC_wRES_wZ->Scale(integral_MC/integral_REC_wRES_wZ);
			h_t_REC_wRES_wE->Scale(integral_MC/integral_REC_wRES_wE);
		}

		// Draw histograms 
		h_t_MC->Draw("same");

		h_t_REC->SetMarkerStyle(20);
		h_t_REC->Draw("PEsame");

		h_t_REC_wRES->SetMarkerStyle(30);
		h_t_REC_wRES->SetMarkerColor(kRed);
		h_t_REC_wRES->Draw("P same");

		h_t_REC_wRES_wZ->SetMarkerStyle(25);
		h_t_REC_wRES_wZ->SetMarkerColor(kGreen+1);
		//h_t_REC_wRES_wZ->Draw("P same");

		h_t_REC_wRES_wE->SetMarkerStyle(5);
		h_t_REC_wRES_wE->SetMarkerColor(kOrange-3);
		//h_t_REC_wRES_wE->Draw("P same");

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
		TLegend* w7 = new TLegend(0.26, 0.64, 0.93, 0.76);
		w7->SetLineColor(kWhite);
		w7->SetFillColor(0);
		w7->SetTextSize(17);
		w7->SetTextFont(45);
		w7->AddEntry(h_t_MC, Form("Sartre %s MC %.f evts", vm_label.Data(),integral_MC), "L");
		w7->AddEntry(h_t_REC, Form("Sartre %s RECO w. EEMC %.f evts", vm_label.Data(),integral_REC), "P");
		w7->AddEntry(h_t_REC_wRES, Form("Sartre %s RECO' new meth. %.f evts", vm_label.Data(),integral_REC_wRES), "P");
		//w7->AddEntry(h_t_REC_wRES_wZ, Form("Sartre %s RECO' new meth. w. z %.f evts", vm_label.Data(),integral_REC_wRES_wZ), "P");
		//w7->AddEntry(h_t_REC_wRES_wE, Form("Sartre %s RECO' new meth. w. E %.f evts", vm_label.Data(),integral_REC_wRES_wE), "P");
		w7->Draw("same");

		TLegend *w8 = new TLegend(0.48,0.61,0.93,0.56);
		w8->SetLineColor(kWhite);
		w8->SetFillColor(0);
		w8->SetTextSize(17);
		w8->SetTextFont(45);
		w8->AddEntry(h_t_REC_wRES,"normalization: #int|#it{t}|_{MC}/#int|#it{t}|_{RECO'}","P");
		//w8->AddEntry(h_t_REC_wRES_wZ,"","P");
		//w8->AddEntry(h_t_REC_wRES_wE,"","P");
		w8->Draw("same");
	
		// Save figure
		totalNorm_canvases_res[i]->Print("./figures/totalNorm_wRES.pdf");
	}


/*	//plot with resolution, no norm
	TCanvas* noNorm_canvases_res[numCanvas];
	for (int i=0; i<1; i++) 
	{
    	TString canvasName = Form("noNorm_res_c%d", i+1);
    	noNorm_canvases_res[i] = new TCanvas(canvasName, canvasName, 1, 1, 600, 600);
    	gPad->SetLogy(1);
    	gPad->SetTicks();
    	gPad->SetLeftMargin(0.15);
    	gPad->SetBottomMargin(0.15);
    	gPad->SetRightMargin(0.01);
		// Generate histogram for each canvas
    	TString noNorm_histName_res = Form("noNorm_res_base%d", i+1);
    	TH1D* noNorm_baseHist_res = makeHist(noNorm_histName_res, "", "|#it{t} | (GeV^{2})", "dN/d|#it{t} | (GeV^{-2}) ", 100, 0, 0.18, kBlack);
    	noNorm_baseHist_res->GetYaxis()->SetRangeUser(8e-2, 8e8);
    	noNorm_baseHist_res->GetXaxis()->SetTitleColor(kBlack);
    	fixedFontHist1D(noNorm_baseHist_res, 1., 1.2);
    	noNorm_baseHist_res->GetYaxis()->SetTitleSize(noNorm_baseHist_res->GetYaxis()->GetTitleSize()*1.5);
    	noNorm_baseHist_res->GetXaxis()->SetTitleSize(noNorm_baseHist_res->GetXaxis()->GetTitleSize()*1.5);
    	noNorm_baseHist_res->GetYaxis()->SetLabelSize(noNorm_baseHist_res->GetYaxis()->GetLabelSize()*1.5);
    	noNorm_baseHist_res->GetXaxis()->SetLabelSize(noNorm_baseHist_res->GetXaxis()->GetLabelSize()*1.5);
    	noNorm_baseHist_res->GetXaxis()->SetNdivisions(4,4,0);
    	noNorm_baseHist_res->GetYaxis()->SetNdivisions(5,5,0);
    	noNorm_baseHist_res->Draw();
		// get counts
		double integral_MC = h_t_MC->Integral();
		double integral_REC = h_t_REC->Integral();
		double integral_REC_wRES = h_t_REC_wRES->Integral();
		double integral_REC_wRES_wZ = h_t_REC_wRES_wZ->Integral();
		double integral_REC_wRES_wE = h_t_REC_wRES_wE->Integral();
		
		// Draw histograms 
		h_t_MC->Draw("same");

		h_t_REC->SetMarkerStyle(20);
		h_t_REC->Draw("PEsame");

		h_t_REC_wRES->SetMarkerStyle(30);
		h_t_REC_wRES->SetMarkerColor(kRed);
		h_t_REC_wRES->Draw("P same");

		h_t_REC_wRES_wZ->SetMarkerStyle(25);
		h_t_REC_wRES_wZ->SetMarkerColor(kGreen+1);
		//h_t_REC_wRES_wZ->Draw("P same");

		h_t_REC_wRES_wE->SetMarkerStyle(5);
		h_t_REC_wRES_wE->SetMarkerColor(kOrange-3);
		//h_t_REC_wRES_wE->Draw("P same");

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
		TLegend* w7 = new TLegend(0.26, 0.64, 0.93, 0.76);
		w7->SetLineColor(kWhite);
		w7->SetFillColor(0);
		w7->SetTextSize(17);
		w7->SetTextFont(45);
		w7->AddEntry(h_t_MC, Form("Sartre %s MC: %.f evts", vm_label.Data(),integral_MC), "L");
		w7->AddEntry(h_t_REC, Form("Sartre %s RECO w. EEMC: %.f evts", vm_label.Data(),integral_REC), "P");
		w7->AddEntry(h_t_REC_wRES, Form("Sartre %s RECO' new method: %.f evts", vm_label.Data(),integral_REC_wRES), "P");
		//w7->AddEntry(h_t_REC_wRES_wZ, Form("Sartre %s RECO' new method w. z: %.f evts", vm_label.Data(),integral_REC_wRES_wZ), "P");
		//w7->AddEntry(h_t_REC_wRES_wE, Form("Sartre %s RECO' new method w. E: %.f evts", vm_label.Data(),integral_REC_wRES_wE), "P");
		w7->Draw("same");

		// Save figure
		noNorm_canvases_res[i]->Print("./figures/noNorm_wRES.pdf");
	}
*/

	//plots with cut, total norm
	TCanvas* totalNorm_canvases_Cut[numCanvas];
	for (int i=0; i<1; i++) 
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
		// get counts
		double integral_MC = h_t_MC->Integral();
		double integral_REC = h_t_REC->Integral();
		double integral_REC_wCut = h_t_REC_wCUT->Integral();
		double integral_REC_wCut_wZ = h_t_REC_wCUT_wZ->Integral();
		double integral_REC_wCut_wE = h_t_REC_wCUT_wZ->Integral();
		// normalize
		if(integral_MC>0 && integral_REC_wCut>0 && integral_REC_wCut_wZ>0 && integral_REC_wCut_wE>0) 
		{
    		h_t_REC_wCUT->Scale(integral_MC/integral_REC_wCut);
			h_t_REC_wCUT_wZ->Scale(integral_MC/integral_REC_wCut_wZ);
			h_t_REC_wCUT_wE->Scale(integral_MC/integral_REC_wCut_wE);
		}

		// Draw histograms with different angle cuts
		h_t_MC->Draw("same");

		h_t_REC->SetMarkerStyle(20);
		h_t_REC->Draw("PEsame");

		h_t_REC_wCUT->SetMarkerStyle(30);
		h_t_REC_wCUT->SetMarkerColor(kRed);
		h_t_REC_wCUT->Draw("P same");

		h_t_REC_wCUT_wZ->SetMarkerStyle(25);
		h_t_REC_wCUT_wZ->SetMarkerColor(kGreen+1);
		//h_t_REC_wCUT_wZ->Draw("P same");

		h_t_REC_wCUT_wE->SetMarkerStyle(5);
		h_t_REC_wCUT_wE->SetMarkerColor(kOrange-3);
		//h_t_REC_wCUT_wE->Draw("P same");

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
		TLegend* w7 = new TLegend(0.26, 0.64, 0.93, 0.76);
		w7->SetLineColor(kWhite);
		w7->SetFillColor(0);
		w7->SetTextSize(17);
		w7->SetTextFont(45);
		w7->AddEntry(h_t_MC, Form("Sartre %s MC %.f evts", vm_label.Data(),integral_MC), "L");
		w7->AddEntry(h_t_REC, Form("Sartre %s RECO w. EEMC %.f evts", vm_label.Data(),integral_REC), "P");
		w7->AddEntry(h_t_REC_wCUT, Form("Sartre %s RECO' #theta = #pi/9 %.f evts", vm_label.Data(),integral_REC_wCut), "P");
		//w7->AddEntry(h_t_REC_wCUT_wZ, Form("Sartre %s RECO' w. z #theta = #pi/9 %.f evts", vm_label.Data(),integral_REC_wCut_wZ), "P");
		//w7->AddEntry(h_t_REC_wCUT_wE, Form("Sartre %s RECO' w. E #theta = #pi/9 %.f evts", vm_label.Data(),integral_REC_wCut_wE), "P");
		w7->Draw("same");

		TLegend *w8 = new TLegend(0.48,0.61,0.93,0.56);
		w8->SetLineColor(kWhite);
		w8->SetFillColor(0);
		w8->SetTextSize(17);
		w8->SetTextFont(45);
		w8->AddEntry(h_t_REC_wCUT,"normalization: #int|#it{t}|_{MC}/#int|#it{t}|_{RECO'}","P");
		//w8->AddEntry(h_t_REC_wCUT_wZ,"","P");
		//w8->AddEntry(h_t_REC_wCUT_wE,"","P");
		w8->Draw("same");

		// Save figure
		totalNorm_canvases_Cut[i]->Print("./figures/totalNorm_wCUT.pdf");
	}


/*	//plots with cut, L norm
	TCanvas* Lnorm_canvases_Cut[numCanvas];
	for (int i=0; i<1; i++) 
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
		// get counts
		double integral_MC = h_t_MC->Integral();
		double integral_REC = h_t_REC->Integral();
		double integral_REC_wCut = h_t_REC_wCUT->Integral();
		double integral_REC_wCut_wZ = h_t_REC_wCUT_wZ->Integral();
		double integral_REC_wCut_wE = h_t_REC_wCUT_wE->Integral();
		// normalize
		h_t_REC_wCUT->Scale((M_PI/2)/(M_PI/9));
		h_t_REC_wCUT_wZ->Scale((M_PI/2)/(M_PI/9));
		h_t_REC_wCUT_wE->Scale((M_PI/2)/(M_PI/9));

		// Draw histograms with different angle cuts
		h_t_MC->Draw("same");

		h_t_REC->SetMarkerStyle(20);
		h_t_REC->Draw("PEsame");

		h_t_REC_wCUT->SetMarkerStyle(30);
		h_t_REC_wCUT->SetMarkerColor(kRed);
		h_t_REC_wCUT->Draw("P same");

		h_t_REC_wCUT_wZ->SetMarkerStyle(25);
		h_t_REC_wCUT_wZ->SetMarkerColor(kGreen+1);
		//h_t_REC_wCUT_wZ->Draw("P same");

		h_t_REC_wCUT_wE->SetMarkerStyle(5);
		h_t_REC_wCUT_wE->SetMarkerColor(kOrange-3);
		//h_t_REC_wCUT_wE->Draw("P same");

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
		TLegend* w7 = new TLegend(0.26, 0.64, 0.93, 0.76);
		w7->SetLineColor(kWhite);
		w7->SetFillColor(0);
		w7->SetTextSize(17);
		w7->SetTextFont(45);
		w7->AddEntry(h_t_MC, Form("Sartre %s MC %.f evts", vm_label.Data(),integral_MC), "L");
		w7->AddEntry(h_t_REC, Form("Sartre %s RECO w. EEMC %.f evts", vm_label.Data(),integral_REC), "P");
		w7->AddEntry(h_t_REC_wCUT, Form("Sartre %s RECO' #theta = #pi/9 %.f evts", vm_label.Data(),integral_REC_wCut), "P");
		//w7->AddEntry(h_t_REC_wCUT_wZ, Form("Sartre %s RECO' w. z #theta = #pi/9 %.f evts", vm_label.Data(),integral_REC_wCut_wZ), "P");
		//w7->AddEntry(h_t_REC_wCUT_wE, Form("Sartre %s RECO' w. E #theta = #pi/9 %.f evts", vm_label.Data(),integral_REC_wCut_wE), "P");
		w7->Draw("same");

		TLegend *w8 = new TLegend(0.48,0.61,0.93,0.56);
		w8->SetLineColor(kWhite);
		w8->SetFillColor(0);
		w8->SetTextSize(17);
		w8->SetTextFont(45);
		w8->AddEntry(h_t_REC_wCUT,"normalization: #frac{#pi}{2} / #theta_{max}","P");
		//w8->AddEntry(h_t_REC_wCUT_wZ,"","P");
		//w8->AddEntry(h_t_REC_wCUT_wE,"","P");
		w8->Draw("same");

		// Save figure
		Lnorm_canvases_Cut[i]->Print("./figures/Lnorm_wCUT.pdf");
	}
*/

/*	//plots with cut, no norm
	TCanvas* noNorm_canvases_Cut[numCanvas];
	for (int i=0; i<1; i++) 
	{
		TString canvasName = Form("noNorm_Cut_c%d", i+1);
		noNorm_canvases_Cut[i] = new TCanvas(canvasName, canvasName, 1, 1, 600, 600);
		gPad->SetLogy(1);
		gPad->SetTicks();
		gPad->SetLeftMargin(0.15);
		gPad->SetBottomMargin(0.15);
		gPad->SetRightMargin(0.01);
		// Generate histogram for each canvas
		TString noNorm_histName_Cut = Form("noNorm_Cut_base%d", i+1);
		TH1D* noNorm_baseHist_Cut = makeHist(noNorm_histName_Cut, "", "|#it{t} | (GeV^{2})", "dN/d|#it{t} | (GeV^{-2}) ", 100, 0, 0.18, kBlack);
		noNorm_baseHist_Cut->GetYaxis()->SetRangeUser(8e-2, 8e8);
		noNorm_baseHist_Cut->GetXaxis()->SetTitleColor(kBlack);
		fixedFontHist1D(noNorm_baseHist_Cut, 1., 1.2);
		noNorm_baseHist_Cut->GetYaxis()->SetTitleSize(noNorm_baseHist_Cut->GetYaxis()->GetTitleSize()*1.5);
		noNorm_baseHist_Cut->GetXaxis()->SetTitleSize(noNorm_baseHist_Cut->GetXaxis()->GetTitleSize()*1.5);
		noNorm_baseHist_Cut->GetYaxis()->SetLabelSize(noNorm_baseHist_Cut->GetYaxis()->GetLabelSize()*1.5);
		noNorm_baseHist_Cut->GetXaxis()->SetLabelSize(noNorm_baseHist_Cut->GetXaxis()->GetLabelSize()*1.5);
		noNorm_baseHist_Cut->GetXaxis()->SetNdivisions(4,4,0);
		noNorm_baseHist_Cut->GetYaxis()->SetNdivisions(5,5,0);
		noNorm_baseHist_Cut->Draw();
		// get counts
		double integral_MC = h_t_MC->Integral();
		double integral_REC = h_t_REC->Integral();
		double integral_REC_wCut = h_t_REC_wCUT->Integral();
		double integral_REC_wCut_wZ = h_t_REC_wCUT_wZ->Integral();
		double integral_REC_wCut_wE = h_t_REC_wCUT_wE->Integral();

		// Draw histograms with different angle cuts
		h_t_MC->Draw("same");

		h_t_REC->SetMarkerStyle(20);
		h_t_REC->Draw("PEsame");

		h_t_REC_wCUT->SetMarkerStyle(30);
		h_t_REC_wCUT->SetMarkerColor(kRed);
		h_t_REC_wCUT->Draw("P same");

		h_t_REC_wCUT_wZ->SetMarkerStyle(25);
		h_t_REC_wCUT_wZ->SetMarkerColor(kGreen+1);
		//h_t_REC_wCUT_wZ->Draw("P same");

		h_t_REC_wCUT_wE->SetMarkerStyle(5);
		h_t_REC_wCUT_wE->SetMarkerColor(kOrange-3);
		//h_t_REC_wCUT_wE->Draw("P same");

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
		TLegend* w7 = new TLegend(0.26, 0.64, 0.93, 0.76);
		w7->SetLineColor(kWhite);
		w7->SetFillColor(0);
		w7->SetTextSize(17);
		w7->SetTextFont(45);
		w7->AddEntry(h_t_MC, Form("Sartre %s MC %.f evts", vm_label.Data(),integral_MC), "L");
		w7->AddEntry(h_t_REC, Form("Sartre %s RECO w. EEMC %.f evts", vm_label.Data(),integral_REC), "P");
		w7->AddEntry(h_t_REC_wCUT, Form("Sartre %s RECO' #theta = #pi/9 %.f evts", vm_label.Data(),integral_REC_wCut), "P");
		//w7->AddEntry(h_t_REC_wCUT_wZ, Form("Sartre %s RECO' w. z #theta = #pi/9 %.f evts", vm_label.Data(),integral_REC_wCut_wZ), "P");
		//w7->AddEntry(h_t_REC_wCUT_wE, Form("Sartre %s RECO' w. E #theta = #pi/9 %.f evts", vm_label.Data(),integral_REC_wCut_wE), "P");
		w7->Draw("same");

		// Save figure
		noNorm_canvases_Cut[i]->Print("./figures/noNorm_wCUT.pdf");
	}
*/

/*	//2d plot truth
	TCanvas* canvases_2d[numCanvas];
	for (int i=0; i<1; i++) 
	{
    	TString canvasName = Form("2d_c%d", i+1);
    	canvases_2d[i] = new TCanvas(canvasName, canvasName, 1, 1, 600, 600);
    	gPad->SetLogz(1);
		gPad->SetLeftMargin(0.15);
		// Draw histograms with resolution
		h_t_REC_2d->Draw();
		// Save figure
		canvases_2d[i]->Print("./figures/2d.pdf");
	}

	//2d plot truth with z component in qx
	TCanvas* canvases_2d_wZ[numCanvas];
	for (int i=0; i<1; i++) 
	{
    	TString canvasName = Form("2d_wZ_c%d", i+1);
    	canvases_2d_wZ[i] = new TCanvas(canvasName, canvasName, 1, 1, 600, 600);
    	gPad->SetLogz(1);
		gPad->SetLeftMargin(0.15);
		// Draw histograms with resolution
		h_t_REC_2d_wZ->Draw();
		// Save figure
		canvases_2d_wZ[i]->Print("./figures/2d_wZ.pdf");
	}

	//2d plot truth with E component in qx
	TCanvas* canvases_2d_wE[numCanvas];
	for (int i=0; i<1; i++) 
	{
    	TString canvasName = Form("2d_wE_c%d", i+1);
    	canvases_2d_wE[i] = new TCanvas(canvasName, canvasName, 1, 1, 600, 600);
    	gPad->SetLogz(1);
		gPad->SetLeftMargin(0.15);
		// Draw histograms with resolution
		h_t_REC_2d_wE->Draw();
		// Save figure
		canvases_2d_wE[i]->Print("./figures/2d_wE.pdf");
	}


	//2d plots with resolution and cut
	TCanvas* canvases_2dresCut[numCanvas];
	for (int i=0; i<numCanvas; i++) 
	{
    	TString canvasName = Form("2dresCut_c%d", i+1);
    	canvases_2dresCut[i] = new TCanvas(canvasName, canvasName, 1, 1, 600, 600);
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
		canvases_2dresCut[i]->Print(Form("./figures/2d_angle%s_wRES_cut.pdf",cleanLabel.Data()));
	}

	//2d plots with resolution and cut with z absorbed into qx
	TCanvas* canvases_2dresCut_wZ[numCanvas];
	for (int i=0; i<numCanvas; i++) 
	{
    	TString canvasName = Form("2dresCut_wZ_c%d", i+1);
    	canvases_2dresCut_wZ[i] = new TCanvas(canvasName, canvasName, 1, 1, 600, 600);
    	gPad->SetLogz(1);
		gPad->SetLeftMargin(0.15);
		// Draw histograms with different angle cuts
		h_t_REC_2d_wRES_cut_wZ[i]->Draw();
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
		canvases_2dresCut_wZ[i]->Print(Form("./figures/2d_angle%s_wZ_wRES_cut.pdf",cleanLabel.Data()));
	}

	//2d plots with resolution and cut with E absorbed into qx
	TCanvas* canvases_2dresCut_wE[numCanvas];
	for (int i=0; i<numCanvas; i++) 
	{
    	TString canvasName = Form("2dresCut_wE_c%d", i+1);
    	canvases_2dresCut_wE[i] = new TCanvas(canvasName, canvasName, 1, 1, 600, 600);
    	gPad->SetLogz(1);
		gPad->SetLeftMargin(0.15);
		// Draw histograms with different angle cuts
		h_t_REC_2d_wRES_cut_wE[i]->Draw();
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
		canvases_2dresCut_wE[i]->Print(Form("./figures/2d_angle%s_wE_wRES_cut.pdf",cleanLabel.Data()));
	}

	//2d plot with resolution 
	TCanvas* canvases_2dres[numCanvas];
	for (int i=0; i<1; i++) 
	{
    	TString canvasName = Form("2dres_c%d", i+1);
    	canvases_2dres[i] = new TCanvas(canvasName, canvasName, 1, 1, 600, 600);
    	gPad->SetLogz(1);
		gPad->SetLeftMargin(0.15);
		// Draw histograms with resolution
		h_t_REC_2d_wRES->Draw();
		// Save figure
		canvases_2dres[i]->Print("./figures/2d_wRES.pdf");
	}
	
	//2d plot with resolution with z absorbed into qx
	TCanvas* canvases_2dres_wZ[numCanvas];
	for (int i=0; i<1; i++) 
	{
    	TString canvasName = Form("2dres_wZ_c%d", i+1);
    	canvases_2dres_wZ[i] = new TCanvas(canvasName, canvasName, 1, 1, 600, 600);
    	gPad->SetLogz(1);
		gPad->SetLeftMargin(0.15);
		// Draw histograms with resolution
		h_t_REC_2d_wRES_wZ->Draw();
		// Save figure
		canvases_2dres_wZ[i]->Print("./figures/2d_wZ_wRES.pdf");
	}

	//2d plot with resolution with E absorbed into qx
	TCanvas* canvases_2dres_wE[numCanvas];
	for (int i=0; i<1; i++) 
	{
    	TString canvasName = Form("2dres_wE_c%d", i+1);
    	canvases_2dres_wE[i] = new TCanvas(canvasName, canvasName, 1, 1, 600, 600);
    	gPad->SetLogz(1);
		gPad->SetLeftMargin(0.15);
		// Draw histograms with resolution
		h_t_REC_2d_wRES_wE->Draw();
		// Save figure
		canvases_2dres_wE[i]->Print("./figures/2d_wE_wRES.pdf");
	}

	//2d plots with cut
	TCanvas* canvases_2dCut[numCanvas];
	for (int i=0; i<numCanvas; i++) 
	{
    	TString canvasName = Form("2dCut_c%d", i+1);
    	canvases_2dCut[i] = new TCanvas(canvasName, canvasName, 1, 1, 600, 600);
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
		canvases_2dCut[i]->Print(Form("./figures/2d_angle%s_wCUT.pdf",cleanLabel.Data()));
	}

	//2d plots with cut with z absorbed into qx
	TCanvas* canvases_2dCut_wZ[numCanvas];
	for (int i=0; i<numCanvas; i++) 
	{
    	TString canvasName = Form("2dCut_wZ_c%d", i+1);
    	canvases_2dCut_wZ[i] = new TCanvas(canvasName, canvasName, 1, 1, 600, 600);
    	gPad->SetLogz(1);
		gPad->SetLeftMargin(0.15);
		// Draw histograms with different angle cuts
		h_t_REC_2d_wCUT_wZ[i]->Draw();

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
		canvases_2dCut_wZ[i]->Print(Form("./figures/2d_angle%s_wZ_wCUT.pdf",cleanLabel.Data()));
	}

	//2d plots with cut with E absorbed into qx
	TCanvas* canvases_2dCut_wE[numCanvas];
	for (int i=0; i<numCanvas; i++) 
	{
    	TString canvasName = Form("2dCut_wE_c%d", i+1);
    	canvases_2dCut_wE[i] = new TCanvas(canvasName, canvasName, 1, 1, 600, 600);
    	gPad->SetLogz(1);
		gPad->SetLeftMargin(0.15);
		// Draw histograms with different angle cuts
		h_t_REC_2d_wCUT_wE[i]->Draw();

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
		canvases_2dCut_wE[i]->Print(Form("./figures/2d_angle%s_wE_wCUT.pdf",cleanLabel.Data()));
	}
*/
}