#include <TFile.h>
#include <TH1D.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>
#include <iostream>
#include <vector>

using namespace std;

void combine_wRES(TString inputDir, TString outputFile) 
{
    TSystemDirectory dir("inputDir", inputDir);
    TList* files = dir.GetListOfFiles();
    if(!files) 
    {
        cout << "Error: No files found in directory " << inputDir << endl;
        return;
    }
    vector<TString> fileNames;
    TSystemFile* file;
    TIter next(files);
    while((file = (TSystemFile*)next())) 
    {
        TString fileName = file->GetName();
        if(fileName.EndsWith("_output.root")) 
        {
            fileNames.push_back(inputDir + "/" + fileName);
        }
    }

    TH1D* h_t_MC_combined = nullptr;
    TH1D* h_t_REC_combined = nullptr;
    TH1D* h_t_REC_wRES_combined = nullptr;

    for(const auto& fileName : fileNames) 
    {
        TFile* inputFile = TFile::Open(fileName);
        if(!inputFile || inputFile->IsZombie()) 
        {
            cout << "Error: Cannot open file " << fileName << endl;
            continue;
        }

        TH1D* h_t_MC = dynamic_cast<TH1D*>(inputFile->Get("h_t_MC"));
        TH1D* h_t_REC = dynamic_cast<TH1D*>(inputFile->Get("h_t_REC"));
        TH1D* h_t_REC_wRES = dynamic_cast<TH1D*>(inputFile->Get("h_t_REC_wRES"));

        if(h_t_MC) 
        {
            if(!h_t_MC_combined) 
            {
                h_t_MC_combined = dynamic_cast<TH1D*>(h_t_MC->Clone("h_t_MC_combined"));
                h_t_MC_combined->SetDirectory(nullptr);
            } 
            else 
            {
                h_t_MC_combined->Add(h_t_MC);
            }
        }

        if(h_t_REC) 
        {
            if(!h_t_REC_combined) 
            {
                h_t_REC_combined = dynamic_cast<TH1D*>(h_t_REC->Clone("h_t_REC_combined"));
                h_t_REC_combined->SetDirectory(nullptr);
            } 
            else 
            {
                h_t_REC_combined->Add(h_t_REC);
            }
        }

        if(h_t_REC_wRES) 
        {
            if(!h_t_REC_wRES_combined) 
            {
                h_t_REC_wRES_combined = dynamic_cast<TH1D*>(h_t_REC_wRES->Clone("h_t_REC_wRES_combined"));
                h_t_REC_wRES_combined->SetDirectory(nullptr);
            } 
            else 
            {
                h_t_REC_wRES_combined->Add(h_t_REC_wRES);
            }
        }
        delete inputFile;
    }

    TFile* outputFileHandle = new TFile(outputFile, "RECREATE");
    if(!outputFileHandle || outputFileHandle->IsZombie()) 
    {
        cout << "Error: Unable to create output file " << outputFile << endl;
        return;
    }
    outputFileHandle->cd();

    if(h_t_MC_combined) 
    {
        cout << "Saving combined histogram: h_t_MC" << endl;
        h_t_MC_combined->Write();
    } 
    else 
    {
        cout << "Error: Combined histogram h_t_MC is null!" << endl;
    }
    
    if(h_t_REC_combined) 
    {
        cout << "Saving combined histogram: h_t_REC" << endl;
        h_t_REC_combined->Write();
    } 
    else 
    {
        cout << "Error: Combined histogram h_t_REC is null!" << endl;
    }

    if(h_t_REC_wRES_combined) 
    {
        cout << "Saving combined histogram: h_t_REC_wRES" << endl;
        h_t_REC_wRES_combined->Write();
    } 
    else 
    {
        cout << "Error: Combined histogram h_t_REC_wRES is null!" << endl;
    }
    
    if(h_t_MC_combined) h_t_MC_combined->Write();
    if(h_t_REC_combined) h_t_REC_combined->Write();
    if(h_t_REC_wRES_combined) h_t_REC_wRES_combined->Write();
    outputFileHandle->Close();

    cout << "Histograms combined and saved to: " << outputFile << endl;
}
