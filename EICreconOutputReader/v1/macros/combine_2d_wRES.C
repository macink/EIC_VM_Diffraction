#include <TFile.h>
#include <TH2D.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>
#include <iostream>
#include <vector>

using namespace std;

void combine_2d_wRES(TString inputDir, TString outputFile) 
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

    TH2D* h_t_REC_2d_wRES_combined = nullptr;

    for(const auto& fileName : fileNames) 
    {
        TFile* inputFile = TFile::Open(fileName);
        if(!inputFile || inputFile->IsZombie()) 
        {
            cout << "Error: Cannot open file " << fileName << endl;
            continue;
        }

        TH2D* h_t_REC_2d_wRES = dynamic_cast<TH2D*>(inputFile->Get("h_t_REC_2d_wRES"));

        if(h_t_REC_2d_wRES) 
        {
            if(!h_t_REC_2d_wRES_combined) 
            {
                h_t_REC_2d_wRES_combined = dynamic_cast<TH2D*>(h_t_REC_2d_wRES->Clone("h_t_REC_2d_wRES_combined"));
                h_t_REC_2d_wRES_combined->SetDirectory(nullptr);
            } 
            else 
            {
                h_t_REC_2d_wRES_combined->Add(h_t_REC_2d_wRES);
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

    if(h_t_REC_2d_wRES_combined) 
    {
        cout << "Saving combined histogram: h_t_REC_2d_wRES" << endl;
        h_t_REC_2d_wRES_combined->Write();
    } 
    else 
    {
        cout << "Error: Combined histogram h_t_REC_2d_wRES is null!" << endl;
    }

    
    if(h_t_REC_2d_wRES_combined) h_t_REC_2d_wRES_combined->Write();
    outputFileHandle->Close();

    cout << "Histograms combined and saved to: " << outputFile << endl;
}

