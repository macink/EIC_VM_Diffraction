#ifndef FormFactor_saturation_data_H
#define FormFactor_saturation_data_H

#include "EICvalueconst.h"

using namespace std;


/********************************************** 
 
 * Form Factor saturation model from scanned data points
 
 * Parameters:
    Vo: Nuclear density (hard sphere approximation) [fm^-3]
    A: Atomic mass number
    R: Radius of nucleus [fm]
    a0: Range of Yukawa potential [fm]
 * Variables:
    t: Four-momentum transfer squared >> t = tx+ty = q^2 [GeV^2]

**********************************************/
class FormFactor_saturation_data 
{
public:
    FormFactor_saturation_data(int bins, double* x_vals, double* y_vals) 
    {
        // Create TGraph from data points
        graph = new TGraph(bins, x_vals, y_vals);

        // Create TH1D from data points
        hist = new TH1D("", "", bins, x_vals[0], x_vals[bins-1]);
        for(int i=0;i<bins;i++) 
        {
            hist->Fill(x_vals[i], y_vals[i]);
        }
    }


    // Functions
    TH1D *getHist() const {return hist;}
    TGraph *getGraph() const {return graph;}


private:
    TGraph *graph;
    TH1D *hist;
};



#endif


<<<<<<< HEAD:FormFactor_saturation_data.h
=======
#endif
>>>>>>> 8fe1b685e5f0f826f3c9fb07ad42e59f8355018e:old_code/FormFactor_saturation_data.h
