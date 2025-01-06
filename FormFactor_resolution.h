#ifndef FormFactor_resolution_H
#define FormFactor_resolution_H

#include "EICvalueconst.h"

using namespace std;




/********************************************** 
 
 * Form Factor with detector resolution: F(q) and F(t) [-]
 
 * Parameters:
    Vo: Nuclear density (hard sphere approximation) [fm^-3]
    A: Atomic mass number
    R: Radius of nucleus [fm]
    a0: Range of Yukawa potential [fm]
 * Variables:
    q: Transverse momentum >> q = sqrt(qx^2+qy^2) [GeV]
    t: Four-momentum transfer squared >> t = tx+ty = q^2 [GeV^2]

**********************************************/
class FormFactor_resolution
{
public:
    FormFactor_resolution(double A_init, double Vo_init, double R_init, double a0_init, double q_min, double q_max, double t_min, double t_max, 
        double phi_min, double phi_max, double qy_min, double qy_max, double qx_prime_min, double qx_prime_max, double ty_min, double ty_max, 
        double tx_prime_min, double tx_prime_max, double bins, double qxBins, double qyBins, double txBins, double tyBins, double x_min,
        double x_max, double y_min, double y_max, double r_min, double r_max, double sigma): A(A_init), Vo(Vo_init), R(R_init), a0(a0_init)
    {
        // Initialize |F(q)|^2 with resolution TF2
        smearedFFq = new TF2("Smeared Form Factor: |F(q)|^{2}",FormFactorq_wResolution_squared,qx_prime_min,qx_prime_max,qy_min,qy_max,5);
        smearedFFq->SetParameters(A,Vo,R,a0,sigma);

        // Initialize |F(t)|^2 with resolution TF2
        smearedFFt = new TF2("Smeared Form Factor: |F(t)|^{2}",FormFactort_wResolution_squared,tx_prime_min,tx_prime_max,ty_min,ty_max,5);
        smearedFFt->SetParameters(A,Vo,R,a0,sigma);

        // Initialize smeared |F(q)|^2 with wedge cut TF1
        FFsmearWedge_q = new TF1("Form Factor |F(q)|^{2} with Resolution and Cut", calc_FormFactor_wResolution_WedgeCut_q, q_min, q_max, 11);
        FFsmearWedge_q->SetParameters(A, Vo, R, a0, phi_min, phi_max, qy_min, qy_max, qx_prime_min, qx_prime_max,sigma);

        // Initialize smeared |F(t)|^2 with wedge cut TF1
        FFsmearWedge_t = new TF1("Form Factor |F(t)|^{2} with Resolution and Cut", calc_FormFactor_wResolution_WedgeCut_t, t_min, t_max, 11);
        FFsmearWedge_t->SetParameters(A, Vo, R, a0, phi_min, phi_max, ty_min, ty_max, tx_prime_min, tx_prime_max,sigma);

        test = new TF1("Cut Form Factor with Resolution: |F(t)|^{2}",[phi_min,phi_max,sigma](double *var, double *par){double t = var[0]; 
            TF1 ffq("ffq",calc_FormFactor_wResolution_WedgeCut_t,phi_min,phi_max,2);
            ffq.SetParameters(t,sigma);
            return ffq.Integral(phi_min,phi_max,1e-12);
        },t_min,t_max,1);
        test->SetParameters(sigma);
    


    }
    


    // Calculations to add resolution to form factor
    static double calc_FormFactor(double q)
    {
	    if(q==0){return 0;}
        else{
        double A = 197, rho0 = 2.12, R = 6.38, a0 = 0.7;
	    const double arg1 = q * R / hbarc;  
	    const double arg2 = hbarc /q;
	    const double sph  = (sin(arg1) - arg1 * cos(arg1)) * 4*pi*rho0 * arg2 * arg2 *arg2/ double(A); 
        const double result = sph * 1/(1. + (a0 * a0 * q*q/(hbarc*hbarc)));  // [-]
	    return result*result;
        }
    }
    static double combine_guassian_formFactor(double *var, double *par)
    {
        double sigma = par[2];  // [GeV]
        double qy = par[0], qx_prime = par[1];
        double qx = var[0];  // integration variable
        double q = sqrt(qx*qx+qy*qy);
        double combined_funcs = calc_FormFactor(q)*exp(-(qx-qx_prime)*(qx-qx_prime)/(2*sigma*sigma));
	    return combined_funcs;
    }
    static double integrate_combined_funcs(double qy, double qx_prime, double sigma)
    {
        double qx_min = 0, qx_max = 5;
        TF1 *new_ff = new TF1("new_ff",combine_guassian_formFactor,qx_min,qx_max,3);  
        new_ff->SetParameters(qy,qx_prime,sigma);
        return new_ff->Integral(qx_min,qx_max,1e-12); // F(qx',qy) [GeV]
    }
    static double FormFactorq_wResolution_squared(double *var, double *par)
    {
        //double qy = var[0], qx_prime = var[1]; 
        double q = par[0], sigma = par[1];
        //double sigma = par[0];
        double phi = var[0];
        double qx_prime = q*sin(phi), qy = q*cos(phi);
        double smeared_formFactor = integrate_combined_funcs(qy,qx_prime,sigma);//*integrate_combined_funcs(qy,qx_prime,sigma);
        return smeared_formFactor;
    }
    static double FormFactort_wResolution_squared(double *var, double *par)
    {
        double t = par[0], sigma = par[1];
        double phi = var[0];
        double q = sqrt(t);
        double qx_prime = q*sin(phi), qy = q*cos(phi);
        //double qy = sqrt(ty), qx_prime = sqrt(tx_prime);
        double smeared_formFactor = integrate_combined_funcs(qy, qx_prime,sigma);//*integrate_combined_funcs(qy, qx_prime,sigma);
        return smeared_formFactor;
    }
    

    // Calculate smeared form factor with wedge cut
    static double calc_FormFactor2_qWedge(double *var, double *par)
    {
        double A = par[0], Vo = par[1], R = par[2], a0 = par[3], q = par[4];
        double phi = var[0];
        double qx = q*sin(phi), qy = q*cos(phi), qq = sqrt(qx*qx+qy*qy);
	    const double arg1 = qq*R / hbarc;  
	    const double arg2 = hbarc / qq;
	    const double arg3  = (sin(arg1)-arg1*cos(arg1)) * 4*pi*Vo*arg2*arg2*arg2 / A; 
        double result = 1/(2*pi)*arg3 / (1. + (a0 * a0 * qq*qq/(hbarc*hbarc)));
        return result*result;  // [-]
    }
    static double calc_FormFactor_wResolution_WedgeCut_q(double* var, double* par)
    {
        double q = var[0];
        double phi_min = par[4], phi_max = par[5], sigma = par[10];
        // Initialize |F(q)|^2 with resolution TF2
        TF2 smearedFFq("Smeared Form Factor: |F(q)|^{2}", FormFactorq_wResolution_squared, par[6], par[7], par[8], par[9], 5);

        smearedFFq.SetParameters(par[0], par[1], par[2], par[3],sigma);

        // Evaluate smeared form factor at t
        double smearedValue = smearedFFq.Eval(q);

        // Initialize |F(q)|^2 wedge cut TF1
        TF1 cutFFq2("Cut Form Factor: |F(q)|^{2}", [phi_min, phi_max,smearedValue](double* var, double* par) {double q = var[0];  // q is variable
            TF1 ffq2("ffq", calc_FormFactor2_qWedge, phi_min, phi_max, 5);
            ffq2.SetParameters(par[0], par[1], par[2], par[3], q);
            double integral = ffq2.Integral(phi_min, phi_max, 1e-12);
            return integral;
        }, par[0], par[1], 4);

        cutFFq2.SetParameters(par[0], par[1], par[2], par[3]);
        cutFFq2.SetParNames("A", "Vo", "R", "a0");

        // Calculate the form factor with resolution and wedge cut
        double result = cutFFq2.Eval(q);
        return result;
    }
    static double calc_FormFactor2_tWedge(double *var, double *par)
    {
        double A = par[0], Vo = par[1], R = par[2], a0 = par[3], t = par[4];
        double phi = var[0];
        double q = sqrt(t), qx = q*sin(phi), qy = q*cos(phi), qq = sqrt(qx*qx+qy*qy);
	    const double arg1 = q*R / hbarc;  
	    const double arg2 = hbarc / q;
	    const double arg3  = (sin(arg1)-arg1*cos(arg1)) * 4*pi*Vo*arg2*arg2*arg2 / A; 
	    double result = 1/(2*pi)*arg3 / (1. + (a0 * a0 * q*q/(hbarc*hbarc))); 
        return result*result; // [-]
    }
    static double calc_FormFactor_wResolution_WedgeCut_t(double* var, double* par)
    {
        double t = var[0];
        double phi_min = par[4], phi_max = par[5], sigma = par[10];
        // Initialize |F(t)|^2 with resolution TF2
        TF2 smearedFFt("Smeared Form Factor: |F(t)|^{2}", FormFactort_wResolution_squared, par[6], par[7], par[8], par[9], 5);

        smearedFFt.SetParameters(par[0], par[1], par[2], par[3],sigma);

        // Evaluate smeared form factor at t
        double smearedValue = smearedFFt.Eval(t);

        // Initialize |F(t)|^2 wedge cut TF1
        TF1 cutFFq2("Cut Form Factor: |F(t)|^{2}", [phi_min, phi_max,smearedValue](double* var, double* par) {double t = var[0];  // t is variable
            TF1 ffq2("ffq", calc_FormFactor2_tWedge, phi_min, phi_max, 5);
            ffq2.SetParameters(par[0], par[1], par[2], par[3], t);
            double integral = ffq2.Integral(phi_min, phi_max, 1e-12);
            return integral;
        }, par[0], par[1], 4);

        cutFFq2.SetParameters(par[0], par[1], par[2], par[3]);
        cutFFq2.SetParNames("A", "Vo", "R", "a0");

        // Calculate the form factor with resolution and wedge cut
        double result = cutFFq2.Eval(t);
        return result;
    }
    

    // Functions
    TF2 *getSmearedFormFactor_q() const {return smearedFFq;}
    TF2 *getSmearedFormFactor_t() const {return smearedFFt;}
    TF1 *getFormFactor_wResWedge_q() const {return FFsmearWedge_q;}
    TF1 *getFormFactor_wResWedge_t() const {return FFsmearWedge_t;}
    TH2D *getTransformFFhist() const {return hankelTransformFF;}
    TF1 *getTest() const {return test;}
    TF1 *getTest2() const {return f1;}
    
    TH2D *create_FormFactorq_2Dcontour(const char *name, const char *title, int qxBins, double qx_min, double qx_max, int qyBins, double qy_min, double qy_max)
    {
        TF2 *FF2 = getSmearedFormFactor_q();
        // Create a 2D histogram to represent the function
        TH2D *hist2D = new TH2D(name, title, qxBins, qx_min, qx_max, qyBins, qy_min, qy_max);
        hist2D->Sumw2();
        for (int i=1;i<=qxBins;i++) 
        {
            for (int j=1;j<=qyBins;j++) 
            {
                double qx = hist2D->GetXaxis()->GetBinCenter(i);
                double qy = hist2D->GetYaxis()->GetBinCenter(j);
                hist2D->SetBinContent(i, j, FF2->Eval(qx, qy));
            }
            //cout << "bin: " << i << endl;
        }
        //hist2D->SetTitle("Form Factor 2D: |F(qx,qy)|^{2}");
        hist2D->GetYaxis()->SetTitle("qy [GeV]");
	    hist2D->GetXaxis()->SetTitle("qx [GeV]");
        hist2D->GetZaxis()->SetTitle("|F(qx,qy)|^{2}");
        return hist2D;
    }
    TH2D *create_FormFactort_2Dcontour(const char *name, const char *title, int txBins, double tx_min, double tx_max, int tyBins, double ty_min, double ty_max)
    {
        TF2 *FF2 = getSmearedFormFactor_t();
        // Create a 2D histogram to represent the function
        TH2D *hist2D = new TH2D(name, title, txBins, tx_min, tx_max, tyBins, ty_min, ty_max);
        hist2D->Sumw2();
        for (int i=1;i<=txBins;i++) 
        {
            for (int j=1;j<=tyBins;j++) 
            {
                double tx = hist2D->GetXaxis()->GetBinCenter(i);
                double ty = hist2D->GetYaxis()->GetBinCenter(j);
                hist2D->SetBinContent(i, j, FF2->Eval(tx, ty));
            }
            //cout << "bin: " << i << endl;
        }
        //hist2D->SetTitle("Form Factor 2D: |F(tx,ty)|^{2}");
        hist2D->GetYaxis()->SetTitle("ty [GeV]");
	    hist2D->GetXaxis()->SetTitle("tx [GeV]");
        hist2D->GetZaxis()->SetTitle("|F(tx,ty)|^{2}");
        return hist2D;
    }


private:
    double A, Vo, R, a0;
    TF1 *FFsmearWedge_q, *FFsmearWedge_t, *test, *f1;
    TF2 *smearedFFq, *smearedFFt;
    TH2D *hankelTransformFF;
};


#endif
