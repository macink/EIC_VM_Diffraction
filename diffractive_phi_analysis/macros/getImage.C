#include "ePICStyle.h"
gROOT->ProcessLine("SetePICStyle()");

/*-----------------------------------------------------
    Fourier Transformation of t Distribution
- Purpose: transforms the momentum transfer distribution 
            (t) to obtain the spatial profile of the nucleus
- Output: root file and histogram
------------------------------------------------------*/

using namespace std;

void getImage()
{
    const char* file = "merged_sartre_18x110_cut.root";
    TFile* input = new TFile(file);
    
    TH1D* hdsigmadt_MC = (TH1D*)input->Get("h_t_MC");
    TH1D* hdsigmadt_REC = (TH1D*)input->Get("h_t_REC_EEMC"); // For method L reco
    TH1D* hdsigmadt_REC_new = (TH1D*)input->Get("h_t_REC_wRES_cut_pi12"); // For new method reco

    int nbins = hdsigmadt_MC->GetNbinsX();
    double dsigmadt_MC,dsigmadt_REC,dsigmadt_REC_new,tBinWidth,t,b,delta,F_b_MC,F_b_REC,F_b_REC_new,result1=0,result2=0,result3=0;

    double t_cut = 0.2;
    double bmin= -12;
    double bmax= 12;
    double noOfBins = 300;
    double hbarc = 0.197;
    
    TH1D* hF_b_MC = new TH1D("hF_b_MC", "F_b_MC", noOfBins, bmin, bmax);
    TH1D* hF_b_REC = new TH1D("hF_b_REC", "F_b_MC", noOfBins, bmin, bmax);
    TH1D* hF_b_REC_new = new TH1D("hF_b_REC_new", "F_b_MC", noOfBins, bmin, bmax);

    // loop to perform the transformation
    for (int j=1; j<=noOfBins; j++) 
    {
        F_b_MC= 0;F_b_REC=0,F_b_REC_new=0;
        b = hF_b_MC->GetBinCenter(j); 
        double prefactor = 1.0/(2*TMath::Pi());

        for (int i= 1; i<=nbins; i++) 
        {

            tBinWidth = hdsigmadt_MC->GetBinWidth(i);
            t = hdsigmadt_MC->GetBinCenter(i); // [GeV^2]
            delta = sqrt(fabs(t)); // = q [GeV]

            dsigmadt_MC = hdsigmadt_MC->GetBinContent(i); // [nb/GeV^2]
            dsigmadt_REC = hdsigmadt_REC->GetBinContent(i); // [nb/GeV^2]
            dsigmadt_REC_new = hdsigmadt_REC_new->GetBinContent(i); // [nb/GeV^2]

            // Convert to fm^2 / GeV^2
            dsigmadt_MC /= 1e7;
            dsigmadt_REC /= 1e7;
            dsigmadt_REC_new /= 1e7;

            //double bessel = TMath::BesselJ0(b*delta/hbarc); // for 2D transformation
            double bessel = 0.0; // for 3D transformation
            double x = b*delta/hbarc;

            if (dsigmadt_MC>0 || dsigmadt_REC>0 || dsigmadt_REC_new>0)
            {
                // for 3D transformation (comment out "if else" statement if want 2D transformation)
                if (fabs(x)<1e-6 || fabs(x-TMath::Pi())<1e-6) 
                {
                    bessel = 0.0;
                } 
                else
                {
                    bessel = sin(x)/x; // j0: spherical bessel function
                }

                double amp_MC = sqrt(dsigmadt_MC);
                double amp_REC = sqrt(dsigmadt_REC);
                double amp_REC_new = sqrt(dsigmadt_REC_new);

                if (t>t_cut) continue;
                // Apply minima sign flips as needed
                if (t>0.02) { amp_MC *= -1; amp_REC *= -1; amp_REC_new *= -1; }
                if (t>0.055) { amp_MC *= -1; amp_REC *= -1; amp_REC_new *= -1;}
                if (t>0.11) { amp_MC *= -1; amp_REC *= -1; amp_REC_new *= -1;}

                result1 = amp_MC * bessel *tBinWidth/2;
                result2 = amp_REC * bessel  *tBinWidth/2;
                result3 = amp_REC_new * bessel * tBinWidth/2;

                F_b_MC += result1;
                F_b_REC += result2;
                F_b_REC_new += result3;
            }
        }
        // for 2d
        //F_b_MC*=prefactor;F_b_MC/=hbarc;
        //F_b_REC*=prefactor;F_b_REC/=hbarc;
        //F_b_REC_new*=prefactor;F_b_REC_new/=hbarc;

        // for 3d
        F_b_MC *= prefactor;F_b_MC /= (2*hbarc);
        F_b_REC *= prefactor;F_b_REC /= (2*hbarc);
        F_b_REC_new *= prefactor;F_b_REC_new /= (2*hbarc);

        hF_b_MC->SetBinContent(j, F_b_MC);
        hF_b_MC->SetBinError(j, F_b_MC * 0.001);
        hF_b_REC->SetBinContent(j, F_b_REC);
        hF_b_REC->SetBinError(j, F_b_REC * 0.001);
        hF_b_REC_new->SetBinContent(j, F_b_REC_new);
        hF_b_REC_new->SetBinError(j, F_b_REC_new * 0.001);
    }

    // make histogram
    TCanvas* c1 = new TCanvas("c1","c1",600,600);
    gPad->SetTicks();
    gPad->SetLeftMargin(0.18);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.12);  
    gPad->SetLogx(0);

    // normalize
    hF_b_MC->Scale(1.0 / hF_b_MC->Integral("width"));
    hF_b_REC->Scale(1.0 / hF_b_REC->Integral("width"));
    hF_b_REC_new->Scale(1.0 / hF_b_REC_new->Integral("width"));

    // calculate ion size based off of the spatial distribution
    double max_MC=hF_b_MC->GetMaximum();
    double max_REC=hF_b_REC->GetMaximum();
    double max_REC_new=hF_b_REC_new->GetMaximum();
    double FWHM_MC = 0.;
    double FWHM_REC = 0.;
    double FWHM_REC_new = 0.;
    for(int i=0;i<hF_b_MC->GetNbinsX();i++)
    {
        if(hF_b_MC->GetBinContent(i+1)>=max_MC/2)
        {
            FWHM_MC=fabs(hF_b_MC->GetBinCenter(i+1));
            break;
        }
    }
    for(int i=0;i<hF_b_REC->GetNbinsX();i++)
    {
        if(hF_b_REC->GetBinContent(i+1)>=max_REC/2)
        {
            FWHM_REC=fabs(hF_b_REC->GetBinCenter(i+1));
            break;
        }
    }
    for(int i=0;i<hF_b_REC_new->GetNbinsX();i++)
    {
        if(hF_b_REC_new->GetBinContent(i+1)>=max_REC_new/2)
        {
            FWHM_REC_new=fabs(hF_b_REC_new->GetBinCenter(i+1));
            break;
        }
    }
    
   /* gStyle->SetOptStat(0); 
    gStyle->SetOptTitle(1);
    gStyle->SetTitleFontSize(.043); 
    gStyle->SetTitleX(.23);
    gStyle->SetTitleY(0.96);*/
    
    //hF_b_MC->SetTitle("3D Fourier-Bessel Transform of |t| Distribution");
    //hF_b_MC->SetTitle("2D Fourier-Bessel Transform of |t| Distribution");
    hF_b_MC->GetXaxis()->SetTitle("b [fm]");
    hF_b_MC->GetYaxis()->SetTitle("F(b)/#scale[0.6]{#int} F(b) db");
    hF_b_MC->GetXaxis()->SetTitleOffset(1.2);
    hF_b_MC->GetYaxis()->SetTitleOffset(1.5);
    hF_b_MC->GetYaxis()->SetRangeUser(-0.05, 0.17);
    //hF_b_MC->SetMarkerStyle(21); 
    //hF_b_MC->SetMarkerColor(kBlack);
    hF_b_MC->SetLineColor(kBlack);
    //hF_b_MC->SetMarkerSize(0.8);
    hF_b_MC->SetLineWidth(2);
    hF_b_MC->GetXaxis()->SetLabelSize(0.04);  
    hF_b_MC->GetYaxis()->SetLabelSize(0.04);

    hF_b_REC->SetMarkerStyle(24); 
    hF_b_REC->SetMarkerColor(kBlue);
    hF_b_REC->SetLineColor(kBlue);
    hF_b_REC->SetMarkerSize(0.8);

    hF_b_REC_new->SetMarkerStyle(30); 
    hF_b_REC_new->SetMarkerColor(kRed);
    hF_b_REC_new->SetLineColor(kRed);
    hF_b_REC_new->SetMarkerSize(0.8);

    hF_b_MC->Draw("same");
    hF_b_REC->Draw("Psame");
    hF_b_REC_new->Draw("Psame");

    TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

    // Add labels
	// Add labels
		TLatex* r42 = new TLatex(0.22, 0.9, "eAu 10x100 GeV");
		r42->SetNDC();
		r42->Draw("same");

    TLatex* label_2 = new TLatex(0.22,0.81, Form("True ion size = %.2f fm",FWHM_MC));
    label_2->SetNDC();
    label_2->SetTextSize(0.03);
    label_2->Draw("same");

    TLatex* label_3 = new TLatex(0.22,0.77, Form("Method L ion size = %.2f fm",FWHM_REC));
    label_3->SetNDC();
    label_3->SetTextSize(0.03);
    label_3->Draw("same");

    TLatex* label_4 = new TLatex(0.22,0.73, Form("Projection method ion size = %.2f fm",FWHM_REC_new));
    label_4->SetNDC();
    label_4->SetTextSize(0.03);
    label_4->Draw("same");

	TLatex* label_5 = new TLatex(0.2, 0.2, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->Draw("same");

    TLegend *leg = new TLegend(0.67, 0.75, 0.85, 0.85); 
    leg->AddEntry(hF_b_MC, " MC", "L");
    leg->AddEntry(hF_b_REC, " Method L", "P"); 
    leg->AddEntry(hF_b_REC_new, " Projection method", "P"); 
    leg->SetBorderSize(0);      
    leg->SetFillStyle(0);       
    leg->SetTextFont(42);       
    leg->SetTextSize(0.03);    
    leg->Draw("same");

    c1->Print("F_b.3d.pdf");

    // save into a root file
    string rootfile= "spaitialProfile.root";
    TFile *hfile = 0;
    hfile = new TFile(rootfile.c_str(), "RECREATE");
    hF_b_MC->Write();
    hF_b_REC->Write();
    hF_b_REC_new->Write();
    
    hfile->Close();

    cout << "True ion size = " << FWHM_MC << " fm " << endl;
    cout << "Method L ion size = " << FWHM_REC  << " fm " << endl;
    cout << "New method ion size = " << FWHM_REC_new  << " fm " << endl;

    cout << rootfile.c_str() <<" written." << endl;
    cout << "All done. Bye." << endl;
}




<<<<<<< HEAD
#include "ePICStyle.h"
gROOT->ProcessLine("SetePICStyle()");

/*-----------------------------------------------------
    Fourier Transformation of t Distribution
- Purpose: transforms the momentum transfer distribution 
            (t) to obtain the spatial profile of the nucleus
- Output: root file and histogram
------------------------------------------------------*/

using namespace std;

void getImage()
{
    const char* file = "merged_sartre_18x110_cut.root";
    TFile* input = new TFile(file);
    
    TH1D* hdsigmadt_MC = (TH1D*)input->Get("h_t_MC");
    TH1D* hdsigmadt_REC = (TH1D*)input->Get("h_t_REC_EEMC"); // For method L reco
    TH1D* hdsigmadt_REC_new = (TH1D*)input->Get("h_t_REC_wRES_cut_pi12"); // For new method reco

    int nbins = hdsigmadt_MC->GetNbinsX();
    double dsigmadt_MC,dsigmadt_REC,dsigmadt_REC_new,tBinWidth,t,b,delta,F_b_MC,F_b_REC,F_b_REC_new,result1=0,result2=0,result3=0;

    double t_cut = 0.2;
    double bmin= -12;
    double bmax= 12;
    double noOfBins = 300;
    double hbarc = 0.197;
    
    TH1D* hF_b_MC = new TH1D("hF_b_MC", "F_b_MC", noOfBins, bmin, bmax);
    TH1D* hF_b_REC = new TH1D("hF_b_REC", "F_b_MC", noOfBins, bmin, bmax);
    TH1D* hF_b_REC_new = new TH1D("hF_b_REC_new", "F_b_MC", noOfBins, bmin, bmax);

    // loop to perform the transformation
    for (int j=1; j<=noOfBins; j++) 
    {
        F_b_MC= 0;F_b_REC=0,F_b_REC_new=0;
        b = hF_b_MC->GetBinCenter(j); 
        double prefactor = 1.0/(2*TMath::Pi());

        for (int i= 1; i<=nbins; i++) 
        {

            tBinWidth = hdsigmadt_MC->GetBinWidth(i);
            t = hdsigmadt_MC->GetBinCenter(i); // [GeV^2]
            delta = sqrt(fabs(t)); // = q [GeV]

            dsigmadt_MC = hdsigmadt_MC->GetBinContent(i); // [nb/GeV^2]
            dsigmadt_REC = hdsigmadt_REC->GetBinContent(i); // [nb/GeV^2]
            dsigmadt_REC_new = hdsigmadt_REC_new->GetBinContent(i); // [nb/GeV^2]

            // Convert to fm^2 / GeV^2
            dsigmadt_MC /= 1e7;
            dsigmadt_REC /= 1e7;
            dsigmadt_REC_new /= 1e7;

            double bessel = TMath::BesselJ0(b*delta/hbarc); // for 2D transformation
            //double bessel = 0.0; // for 3D transformation
            double x = b*delta/hbarc;

            if (dsigmadt_MC>0 || dsigmadt_REC>0 || dsigmadt_REC_new>0)
            {
                // for 3D transformation (comment out "if else" statement if want 2D transformation)
                if (fabs(x)<1e-6 || fabs(x-TMath::Pi())<1e-6) 
                {
                    bessel = 0.0;
                } 
                else
                {
                    bessel = sin(x)/x; // j0: spherical bessel function
                }

                double amp_MC = sqrt(dsigmadt_MC);
                double amp_REC = sqrt(dsigmadt_REC);
                double amp_REC_new = sqrt(dsigmadt_REC_new);

                if (t>t_cut) continue;
                // Apply minima sign flips as needed
                if (t>0.02) { amp_MC *= -1; amp_REC *= -1; amp_REC_new *= -1; }
                if (t>0.055) { amp_MC *= -1; amp_REC *= -1; amp_REC_new *= -1;}
                if (t>0.11) { amp_MC *= -1; amp_REC *= -1; amp_REC_new *= -1;}

                result1 = amp_MC * bessel *tBinWidth/2;
                result2 = amp_REC * bessel  *tBinWidth/2;
                result3 = amp_REC_new * bessel * tBinWidth/2;

                F_b_MC += result1;
                F_b_REC += result2;
                F_b_REC_new += result3;
            }
        }
        // for 2d
        //F_b_MC*=prefactor;F_b_MC/=hbarc;
        //F_b_REC*=prefactor;F_b_REC/=hbarc;
        //F_b_REC_new*=prefactor;F_b_REC_new/=hbarc;

        // for 3d
        F_b_MC *= prefactor;F_b_MC /= (2*hbarc);
        F_b_REC *= prefactor;F_b_REC /= (2*hbarc);
        F_b_REC_new *= prefactor;F_b_REC_new /= (2*hbarc);

        hF_b_MC->SetBinContent(j, F_b_MC);
        hF_b_MC->SetBinError(j, F_b_MC * 0.001);
        hF_b_REC->SetBinContent(j, F_b_REC);
        hF_b_REC->SetBinError(j, F_b_REC * 0.001);
        hF_b_REC_new->SetBinContent(j, F_b_REC_new);
        hF_b_REC_new->SetBinError(j, F_b_REC_new * 0.001);
    }

    // make histogram
    TCanvas* c1 = new TCanvas("c1","c1",600,600);
    gPad->SetTicks();
    gPad->SetLeftMargin(0.18);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.12);  
    gPad->SetLogx(0);

    // normalize
    hF_b_MC->Scale(1.0 / hF_b_MC->Integral("width"));
    hF_b_REC->Scale(1.0 / hF_b_REC->Integral("width"));
    hF_b_REC_new->Scale(1.0 / hF_b_REC_new->Integral("width"));

    // calculate ion size based off of the spatial distribution
    double max_MC=hF_b_MC->GetMaximum();
    double max_REC=hF_b_REC->GetMaximum();
    double max_REC_new=hF_b_REC_new->GetMaximum();
    double FWHM_MC = 0.;
    double FWHM_REC = 0.;
    double FWHM_REC_new = 0.;
    for(int i=0;i<hF_b_MC->GetNbinsX();i++)
    {
        if(hF_b_MC->GetBinContent(i+1)>=max_MC/2)
        {
            FWHM_MC=fabs(hF_b_MC->GetBinCenter(i+1));
            break;
        }
    }
    for(int i=0;i<hF_b_REC->GetNbinsX();i++)
    {
        if(hF_b_REC->GetBinContent(i+1)>=max_REC/2)
        {
            FWHM_REC=fabs(hF_b_REC->GetBinCenter(i+1));
            break;
        }
    }
    for(int i=0;i<hF_b_REC_new->GetNbinsX();i++)
    {
        if(hF_b_REC_new->GetBinContent(i+1)>=max_REC_new/2)
        {
            FWHM_REC_new=fabs(hF_b_REC_new->GetBinCenter(i+1));
            break;
        }
    }
    
    gStyle->SetOptStat(0); 
    gStyle->SetOptTitle(1);
    gStyle->SetTitleFontSize(.043); 
    gStyle->SetTitleX(.23);
    gStyle->SetTitleY(0.96);
    
    hF_b_MC->SetTitle("3D Fourier-Bessel Transform of |t| Distribution");
    //hF_b_MC->SetTitle("2D Fourier-Bessel Transform of |t| Distribution");
    hF_b_MC->GetXaxis()->SetTitle("b [fm]");
    hF_b_MC->GetYaxis()->SetTitle("F(b)/#scale[0.6]{#int} F(b) db");
    hF_b_MC->GetXaxis()->SetTitleOffset(1.2);
    hF_b_MC->GetYaxis()->SetTitleOffset(1.5);
    hF_b_MC->GetYaxis()->SetRangeUser(-0.05, 0.17);
    hF_b_MC->SetMarkerStyle(21); 
    hF_b_MC->SetMarkerColor(kBlack);
    hF_b_MC->SetLineColor(kBlack);
    hF_b_MC->SetMarkerSize(0.8);
    hF_b_MC->SetLineWidth(2);
    hF_b_MC->GetXaxis()->SetLabelSize(0.04);  
    hF_b_MC->GetYaxis()->SetLabelSize(0.04);

    hF_b_REC->SetMarkerStyle(24); 
    hF_b_REC->SetMarkerColor(kBlue);
    hF_b_REC->SetLineColor(kBlue);
    hF_b_REC->SetMarkerSize(0.8);

    hF_b_REC_new->SetMarkerStyle(29); 
    hF_b_REC_new->SetMarkerColor(kRed);
    hF_b_REC_new->SetLineColor(kRed);
    hF_b_REC_new->SetMarkerSize(0.8);

    hF_b_MC->Draw("same");
    hF_b_REC->Draw("Psame");
    hF_b_REC_new->Draw("Psame");

    TLegend *leg = new TLegend(0.65, 0.75, 0.8, 0.85); 
    leg->AddEntry(hF_b_MC, "MC t_{max} = 0.2", "P");
    leg->AddEntry(hF_b_REC, "Method L t_{max} = 0.2", "P"); 
    leg->AddEntry(hF_b_REC_new, "New method t_{max} = 0.2", "P"); 
    leg->SetBorderSize(0);      
    leg->SetFillStyle(0);       
    leg->SetTextFont(42);       
    leg->SetTextSize(0.03);    
    leg->Draw("same");

    TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

    // Add labels
	TLatex* label_1 = new TLatex(0.48, 0.2, "EPIC: eAu 10x100 GeV");
	label_1->SetNDC();
	label_1->Draw("same");

    TLatex* label_2 = new TLatex(0.21,0.81, Form("True ion size = %.2f fm",FWHM_MC));
    label_2->SetNDC();
    label_2->SetTextSize(0.03);
    label_2->Draw("same");

    TLatex* label_3 = new TLatex(0.21,0.78, Form("Method L ion size = %.2f fm",FWHM_REC));
    label_3->SetNDC();
    label_3->SetTextSize(0.03);
    label_3->Draw("same");

    TLatex* label_4 = new TLatex(0.21,0.75, Form("New method ion size = %.2f fm",FWHM_REC_new));
    label_4->SetNDC();
    label_4->SetTextSize(0.03);
    label_4->Draw("same");

	TLatex* label_5 = new TLatex(0.2, 0.2, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->Draw("same");
    c1->Print("F_b.3d.pdf");

    // save into a root file
    string rootfile= "spaitialProfile.root";
    TFile *hfile = 0;
    hfile = new TFile(rootfile.c_str(), "RECREATE");
    hF_b_MC->Write();
    hF_b_REC->Write();
    hF_b_REC_new->Write();
    
    hfile->Close();

    cout << "True ion size = " << FWHM_MC << " fm " << endl;
    cout << "Method L ion size = " << FWHM_REC  << " fm " << endl;
    cout << "New method ion size = " << FWHM_REC_new  << " fm " << endl;

    cout << rootfile.c_str() <<" written." << endl;
    cout << "All done. Bye." << endl;
}



=======
#include "ePICStyle.h"
gROOT->ProcessLine("SetePICStyle()");

/*-----------------------------------------------------
    Fourier Transformation of t Distribution
- Purpose: transforms the momentum transfer distribution 
            (t) to obtain the spatial profile of the nucleus
- Output: root file and histogram
------------------------------------------------------*/

using namespace std;

void getImage()
{
    const char* file = "merged_sartre_18x110_cut.root";
    TFile* input = new TFile(file);
    
    TH1D* hdsigmadt_MC = (TH1D*)input->Get("h_t_MC");
    TH1D* hdsigmadt_REC = (TH1D*)input->Get("h_t_REC_EEMC"); // For method L reco
    TH1D* hdsigmadt_REC_new = (TH1D*)input->Get("h_t_REC_wRES_cut_pi12"); // For new method reco

    int nbins = hdsigmadt_MC->GetNbinsX();
    double dsigmadt_MC,dsigmadt_REC,dsigmadt_REC_new,tBinWidth,t,b,delta,F_b_MC,F_b_REC,F_b_REC_new,result1=0,result2=0,result3=0;

    double t_cut = 0.2;
    double bmin= -12;
    double bmax= 12;
    double noOfBins = 300;
    double hbarc = 0.197;
    
    TH1D* hF_b_MC = new TH1D("hF_b_MC", "F_b_MC", noOfBins, bmin, bmax);
    TH1D* hF_b_REC = new TH1D("hF_b_REC", "F_b_MC", noOfBins, bmin, bmax);
    TH1D* hF_b_REC_new = new TH1D("hF_b_REC_new", "F_b_MC", noOfBins, bmin, bmax);

    // loop to perform the transformation
    for (int j=1; j<=noOfBins; j++) 
    {
        F_b_MC= 0;F_b_REC=0,F_b_REC_new=0;
        b = hF_b_MC->GetBinCenter(j); 
        double prefactor = 1.0/(2*TMath::Pi());

        for (int i= 1; i<=nbins; i++) 
        {

            tBinWidth = hdsigmadt_MC->GetBinWidth(i);
            t = hdsigmadt_MC->GetBinCenter(i); // [GeV^2]
            delta = sqrt(fabs(t)); // = q [GeV]

            dsigmadt_MC = hdsigmadt_MC->GetBinContent(i); // [nb/GeV^2]
            dsigmadt_REC = hdsigmadt_REC->GetBinContent(i); // [nb/GeV^2]
            dsigmadt_REC_new = hdsigmadt_REC_new->GetBinContent(i); // [nb/GeV^2]

            // Convert to fm^2 / GeV^2
            dsigmadt_MC /= 1e7;
            dsigmadt_REC /= 1e7;
            dsigmadt_REC_new /= 1e7;

            double bessel = TMath::BesselJ0(b*delta/hbarc); // for 2D transformation
            //double bessel = 0.0; // for 3D transformation
            double x = b*delta/hbarc;

            if (dsigmadt_MC>0 || dsigmadt_REC>0 || dsigmadt_REC_new>0)
            {
                // for 3D transformation (comment out "if else" statement if want 2D transformation)
                if (fabs(x)<1e-6 || fabs(x-TMath::Pi())<1e-6) 
                {
                    bessel = 0.0;
                } 
                else
                {
                    bessel = sin(x)/x; // j0: spherical bessel function
                }

                double amp_MC = sqrt(dsigmadt_MC);
                double amp_REC = sqrt(dsigmadt_REC);
                double amp_REC_new = sqrt(dsigmadt_REC_new);

                if (t>t_cut) continue;
                // Apply minima sign flips as needed
                if (t>0.02) { amp_MC *= -1; amp_REC *= -1; amp_REC_new *= -1; }
                if (t>0.055) { amp_MC *= -1; amp_REC *= -1; amp_REC_new *= -1;}
                if (t>0.11) { amp_MC *= -1; amp_REC *= -1; amp_REC_new *= -1;}

                result1 = amp_MC * bessel *tBinWidth/2;
                result2 = amp_REC * bessel  *tBinWidth/2;
                result3 = amp_REC_new * bessel * tBinWidth/2;

                F_b_MC += result1;
                F_b_REC += result2;
                F_b_REC_new += result3;
            }
        }
        // for 2d
        //F_b_MC*=prefactor;F_b_MC/=hbarc;
        //F_b_REC*=prefactor;F_b_REC/=hbarc;
        //F_b_REC_new*=prefactor;F_b_REC_new/=hbarc;

        // for 3d
        F_b_MC *= prefactor;F_b_MC /= (2*hbarc);
        F_b_REC *= prefactor;F_b_REC /= (2*hbarc);
        F_b_REC_new *= prefactor;F_b_REC_new /= (2*hbarc);

        hF_b_MC->SetBinContent(j, F_b_MC);
        hF_b_MC->SetBinError(j, F_b_MC * 0.001);
        hF_b_REC->SetBinContent(j, F_b_REC);
        hF_b_REC->SetBinError(j, F_b_REC * 0.001);
        hF_b_REC_new->SetBinContent(j, F_b_REC_new);
        hF_b_REC_new->SetBinError(j, F_b_REC_new * 0.001);
    }

    // make histogram
    TCanvas* c1 = new TCanvas("c1","c1",600,600);
    gPad->SetTicks();
    gPad->SetLeftMargin(0.18);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.12);  
    gPad->SetLogx(0);

    // normalize
    hF_b_MC->Scale(1.0 / hF_b_MC->Integral("width"));
    hF_b_REC->Scale(1.0 / hF_b_REC->Integral("width"));
    hF_b_REC_new->Scale(1.0 / hF_b_REC_new->Integral("width"));

    // calculate ion size based off of the spatial distribution
    double max_MC=hF_b_MC->GetMaximum();
    double max_REC=hF_b_REC->GetMaximum();
    double max_REC_new=hF_b_REC_new->GetMaximum();
    double FWHM_MC = 0.;
    double FWHM_REC = 0.;
    double FWHM_REC_new = 0.;
    for(int i=0;i<hF_b_MC->GetNbinsX();i++)
    {
        if(hF_b_MC->GetBinContent(i+1)>=max_MC/2)
        {
            FWHM_MC=fabs(hF_b_MC->GetBinCenter(i+1));
            break;
        }
    }
    for(int i=0;i<hF_b_REC->GetNbinsX();i++)
    {
        if(hF_b_REC->GetBinContent(i+1)>=max_REC/2)
        {
            FWHM_REC=fabs(hF_b_REC->GetBinCenter(i+1));
            break;
        }
    }
    for(int i=0;i<hF_b_REC_new->GetNbinsX();i++)
    {
        if(hF_b_REC_new->GetBinContent(i+1)>=max_REC_new/2)
        {
            FWHM_REC_new=fabs(hF_b_REC_new->GetBinCenter(i+1));
            break;
        }
    }
    
    gStyle->SetOptStat(0); 
    gStyle->SetOptTitle(1);
    gStyle->SetTitleFontSize(.043); 
    gStyle->SetTitleX(.23);
    gStyle->SetTitleY(0.96);
    
    hF_b_MC->SetTitle("3D Fourier-Bessel Transform of |t| Distribution");
    //hF_b_MC->SetTitle("2D Fourier-Bessel Transform of |t| Distribution");
    hF_b_MC->GetXaxis()->SetTitle("b [fm]");
    hF_b_MC->GetYaxis()->SetTitle("F(b)/#scale[0.6]{#int} F(b) db");
    hF_b_MC->GetXaxis()->SetTitleOffset(1.2);
    hF_b_MC->GetYaxis()->SetTitleOffset(1.5);
    hF_b_MC->GetYaxis()->SetRangeUser(-0.05, 0.17);
    hF_b_MC->SetMarkerStyle(21); 
    hF_b_MC->SetMarkerColor(kBlack);
    hF_b_MC->SetLineColor(kBlack);
    hF_b_MC->SetMarkerSize(0.8);
    hF_b_MC->SetLineWidth(2);
    hF_b_MC->GetXaxis()->SetLabelSize(0.04);  
    hF_b_MC->GetYaxis()->SetLabelSize(0.04);

    hF_b_REC->SetMarkerStyle(24); 
    hF_b_REC->SetMarkerColor(kBlue);
    hF_b_REC->SetLineColor(kBlue);
    hF_b_REC->SetMarkerSize(0.8);

    hF_b_REC_new->SetMarkerStyle(29); 
    hF_b_REC_new->SetMarkerColor(kRed);
    hF_b_REC_new->SetLineColor(kRed);
    hF_b_REC_new->SetMarkerSize(0.8);

    hF_b_MC->Draw("same");
    hF_b_REC->Draw("Psame");
    hF_b_REC_new->Draw("Psame");

    TLegend *leg = new TLegend(0.65, 0.75, 0.8, 0.85); 
    leg->AddEntry(hF_b_MC, "MC t_{max} = 0.2", "P");
    leg->AddEntry(hF_b_REC, "Method L t_{max} = 0.2", "P"); 
    leg->AddEntry(hF_b_REC_new, "New method t_{max} = 0.2", "P"); 
    leg->SetBorderSize(0);      
    leg->SetFillStyle(0);       
    leg->SetTextFont(42);       
    leg->SetTextSize(0.03);    
    leg->Draw("same");

    TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

    // Add labels
	TLatex* label_1 = new TLatex(0.48, 0.2, "EPIC: eAu 10x100 GeV");
	label_1->SetNDC();
	label_1->Draw("same");

    TLatex* label_2 = new TLatex(0.21,0.81, Form("True ion size = %.2f fm",FWHM_MC));
    label_2->SetNDC();
    label_2->SetTextSize(0.03);
    label_2->Draw("same");

    TLatex* label_3 = new TLatex(0.21,0.78, Form("Method L ion size = %.2f fm",FWHM_REC));
    label_3->SetNDC();
    label_3->SetTextSize(0.03);
    label_3->Draw("same");

    TLatex* label_4 = new TLatex(0.21,0.75, Form("New method ion size = %.2f fm",FWHM_REC_new));
    label_4->SetNDC();
    label_4->SetTextSize(0.03);
    label_4->Draw("same");

	TLatex* label_5 = new TLatex(0.2, 0.2, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->Draw("same");
    c1->Print("F_b.3d.pdf");

    // save into a root file
    string rootfile= "spaitialProfile.root";
    TFile *hfile = 0;
    hfile = new TFile(rootfile.c_str(), "RECREATE");
    hF_b_MC->Write();
    hF_b_REC->Write();
    hF_b_REC_new->Write();
    
    hfile->Close();

    cout << "True ion size = " << FWHM_MC << " fm " << endl;
    cout << "Method L ion size = " << FWHM_REC  << " fm " << endl;
    cout << "New method ion size = " << FWHM_REC_new  << " fm " << endl;

    cout << rootfile.c_str() <<" written." << endl;
    cout << "All done. Bye." << endl;
}



>>>>>>> 3e7449a (Reorganized analysis macros)
