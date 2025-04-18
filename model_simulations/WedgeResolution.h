#ifndef WedgeResolution_H
#define WedgeResolution_H

#include "EICvalueconst.h"

using namespace std;



/************************************************ 
 
 *  Woods-Saxon Dentsity Distribution 1D: G(r) [1/fm^-3]
 
 * Parameters:
    Vo: Nuclear density (hard sphere approximation) [fm^-3]
    R: Radius of nucleus [fm]
    a: Skin-depth [fm]
 * Variables:
    r: Radial distance from center of nucleus [fm]
        r = sqrt(x^2+y^2)
 
**************************************************/ 
class WedgeResolution {
    
    public:
        WedgeResolution(TF1* inputTF1, double A_init, double Vo_init, double R_init, double a0_init, 
                        double t_min, double t_max, double bins, double phi_min, double phi_max, double sigma, double r_min, double r_max) {
            
    // Retrieve TF1 and TH1D from Class1
    //TF1* dsigma_dt = class1_obj.get_fun_dsigma(); // Get TF1 from Class1
    //hist1D = class1_obj.get_hist_dsigma();        // Get TH1D from Class1
                            // Initialize |F(t)|^2 wedge cut TF1

                    
                    // Validate the input TF1
                    if (!inputTF1) {
                        throw std::invalid_argument("inputTF1 cannot be nullptr");
                    }
                
                    cutFFt2 = new TF1("Cut Form Factor: |F(t)|^{2}", [inputTF1, phi_min, phi_max, sigma](double* var, double* /*par*/) -> double {
                        double t = var[0];
                        if (t < 0) {
                            std::cerr << "Error: Negative t encountered!" << std::endl;
                            return 0.0;
                        }
                    
                        auto FF_cut_wRes_lambda = [inputTF1, sigma](double* x, double* par) -> double {
                            return FF_cut_wRes(x, par, inputTF1);
                        };
                    
                        TF1 fft2("temp", FF_cut_wRes_lambda, phi_min, phi_max, 2);
                        fft2.SetParameters(t, sigma);
                    
                        double integral = 0.0;
                        try {
                            integral = fft2.Integral(phi_min, phi_max, 1e-12);
                        } catch (const std::exception& e) {
                            std::cerr << "Error during phi integration: " << e.what() << std::endl;
                            return 0.0;
                        }
                    
                        double dsigma_value = inputTF1->Eval(t);
                        std::cout << "t: " << t << ", Integral: " << integral << ", dsigma_value: " << dsigma_value << std::endl;
                    
                        return integral * dsigma_value;
                    }, t_min, t_max, 0);
                    
                
 /*       cutFFt2 = new TF1("Cut Form Factor: |F(t)|^{2}", [this, phi_min, phi_max, sigma] (double *var, double *par)
        {
            double t = var[0];  // t is variable
            TF1 fft2("", FF_cut_wRes, phi_min, phi_max, 2);
            fft2.SetParameters(t, sigma);
            double integral = fft2.Integral(phi_min, phi_max, 1e-12);
            return integral;
        }, t_min, t_max, 0);

    cutFFt2 = new TF1("Cut Form Factor: |F(t)|^{2}", [dsigma_dt, phi_min, phi_max, sigma](double* var, double* //par) -> double {
        double t = var[0];  // Independent variable t
    
        // Call wedge_FF directly to integrate over phi
        double integral = wedge_FF(t, sigma, phi_min, phi_max, dsigma_dt);
    
        double dsigma_value = inputTF1->Eval(t);  // Evaluate dsigma_dt
        cout << "dsigma" << dsigma_value << endl;
        return integral * dsigma_value;           // Combine wedge_FF result with dsigma_dt
    }, t_min, t_max, 0);
    
 
    // Initialize |F(t)|^2 wedge cut TF1, incorporating the Class1 TF1
    cutFFt2 = new TF1("Cut Form Factor: |F(t)|^{2}", [dsigma_dt, phi_min, phi_max, sigma](double* var, double* par) {
        double t = var[0];  // Independent variable t
        TF1 fft2("temp", FF_cut_wRes, phi_min, phi_max, 2);
        fft2.SetParameters(t, sigma);
        // Use dsigma_dt->Eval(t) from Class1 TF1 in the calculation
        double dsigma_value = dsigma_dt->Eval(t);
        double integral = fft2.Integral(phi_min, phi_max, 1e-12);
        return integral * dsigma_value;  // Combine with Class1's TF1 evaluation
    }, t_min, t_max, 0);*/

    // Initialize and update the histogram with values from cutFFt2
    hist1D->Sumw2();
    for (int i = 0; i < bins; i++) {
        double t = hist1D->GetXaxis()->GetBinCenter(i + 1);
        double value = cutFFt2->Eval(t);
        hist1D->SetBinContent(i + 1, value); 
    }
        }
    /*
    // Integrand to calcuate smearing
    static double guassian(double *x, double *par, TF1* inputTF1=nullptr)
    {
        double qy = par[0], qx_prime = par[1], sigma = par[2];
        //double qy = sqrt(ty), qx_prime = sqrt(tx_prime);
        double qx = x[0];  // integration variable
        double q = sqrt(qx*qx+qy*qy);
        double t = q*q;
        //cout << "q: " << q << endl;
        // Replace inte_y(t) with values from TF1 or TH1D
    double value = 0.0;
    if (inputTF1) {
        value = inputTF1->Eval(t); // Use TF1 to evaluate at t
    }else {
        throw std::runtime_error("Neither TF1 nor TH1D is provided!");
    }
	    return value*exp(-(qx-qx_prime)*(qx-qx_prime)/(2*sigma*sigma));
    }
    static TF1* createGaussianWithCapture(TF1* dsigma_dt, double qy, double qx_prime, double sigma, double qqx_min, double qqx_max) {
        // Capturing lambda that wraps the logic of guassian()
        auto gaussianLambda = [dsigma_dt, qy, qx_prime, sigma](double* x, double* //par) -> double {
            double qx = x[0];  // Integration variable
            double q = sqrt(qx * qx + qy * qy);
            double t = q * q;
    
            // Use TF1 or TH1D to evaluate values at t
            double value = 0.0;
            if (dsigma_dt) {
                value = dsigma_dt->Eval(t);  // Use TF1 to evaluate at t
            } else {
                throw std::runtime_error("Neither TF1 nor TH1D is provided!");
            }
    
            // Perform the Gaussian calculation
            cout << "value: " << value << endl;
            return value * exp(-(qx - qx_prime) * (qx - qx_prime) / (2 * sigma * sigma));
        };
    
        // Create a TF1 object with the capturing lambda
        return new TF1("gaussian", gaussianLambda, qqx_min, qqx_max, 0);
    }
    
    static double FF_wRes(TF1* dsigma_dt, double qy, double qx_prime, double sigma) {
        double qqx_min = 0, qqx_max = 20;
    
        // Use a capturing lambda to evaluate the Gaussian
        auto gaussianLambda = [dsigma_dt, qy, qx_prime, sigma](double* x, double* //par) -> double {
            double qx = x[0];
            double q = sqrt(qx * qx + qy * qy);
            double t = q * q;
    
            // Evaluate dsigma_dt at t
            if (!dsigma_dt) {
                throw std::runtime_error("FF_wRes requires a valid input TF1.");
            }
            double value = dsigma_dt->Eval(t);
    
            // Return Gaussian calculation
            return value * exp(-(qx - qx_prime) * (qx - qx_prime) / (2 * sigma * sigma));
        };
    
        // Create TF1 to perform the integration
        TF1* f = new TF1("gaussian", gaussianLambda, qqx_min, qqx_max, 0);
    
        double result = f->Integral(qqx_min, qqx_max, 1e-12);
        cout << "result: " << result << endl;
        delete f; // Clean up dynamically allocated TF1
        return result;
    }
    
    static double wedge_FF(double t, double sigma, double phi_min, double phi_max, TF1* dsigma_dt) {
        // Create a lambda function that captures dsigma_dt and performs the calculation
        auto FF_cut_wRes_lambda = [dsigma_dt, sigma](double* x, double* par) -> double {
            double t_local = par[0];      // Parameter t
            double sigma_local = par[1]; // Parameter sigma
            double phi = x[0];           // Integration variable
            double q = sqrt(t_local);    // Compute q from t
            double qx_prime = q * sin(phi);
            double qy = q * cos(phi);
    
            // Call FF_wRes with captured TF1
            return FF_wRes(dsigma_dt, qy, qx_prime, sigma_local);
        };
    
        // Create a TF1 to perform the integral
        TF1* f = new TF1("", FF_cut_wRes_lambda, phi_min, phi_max, 2);
        f->SetParameters(t, sigma);
    
        // Perform the integral over phi
        double result = f->Integral(phi_min, phi_max, 1e-12);
        delete f; // Clean up TF1
        return result;
    }
    */
/*    
    // Integral over qx -> new form factor with smearing
    static double FF_wRes(TF1* dsigma_dt,double qy, double qx_prime, double sigma) 
    {
        double qqx_min = 0, qqx_max = 20;
        //TF1 *f = new TF1("f", guassian, qqx_min, qqx_max, 3);  
        TF1 *f = createGaussianWithCapture(dsigma_dt,qy, qx_prime, sigma, qqx_min, qqx_max);
        //f->SetParameters(qy, qx_prime, sigma);
        double result = f->Integral(qqx_min, qqx_max, 1e-12); // F(qx',qy) [GeV]
        cout << "result: " << result << endl;
        return result; 
    }
    // Integrand for wedge cut
    static double FF_cut_wRes(double *x, double *par,TF1* dsigma_dt)
    {
        double t = par[0], sigma = par[1];
        double phi = x[0];  // integration variable
        double q = sqrt(t), qx_prime = q*sin(phi), qy = q*cos(phi); 
        // call form factor with resolution
	    return FF_wRes(dsigma_dt,qy, qx_prime, sigma); // [GeV^2]
    }
    // Integral to integrate over theta
    static double wedge_FF(double t, double sigma, double phi_min, double phi_max)
    {
        TF1 *f = new TF1("", FF_cut_wRes, phi_min, phi_max, 2);  
        f->SetParameters(t, sigma);
        return f->Integral(phi_min, phi_max, 1e-12); //[GeV]
    }*/
   
    static double guassian(double* x, double* par, TF1* inputTF1) {
        double qy = par[0], qx_prime = par[1], sigma = par[2];
        double qx = x[0];  // Integration variable
        double q = sqrt(qx * qx + qy * qy);
    
        if (!inputTF1) {
            std::cerr << "Error: inputTF1 is nullptr!" << std::endl;
            return 0.0;
        }
    
        double FF_value = inputTF1->Eval(q);  // Evaluate form factor
        std::cout << "qx: " << qx << ", q: " << q << ", FF_value: " << FF_value << std::endl;
    
        return FF_value * exp(-(qx - qx_prime) * (qx - qx_prime) / (2 * sigma * sigma));
    }
    
    static double FF_wRes(double qy, double qx_prime, double sigma, TF1* inputTF1) {
        double qqx_min = 0, qqx_max = 20;
    
        // Validate inputTF1
        if (!inputTF1) {
            std::cerr << "Error: inputTF1 is nullptr!" << std::endl;
            return 0.0;
        }
    
        // Use a lambda for the TF1 evaluation
        auto gaussianLambda = [inputTF1, qy, qx_prime, sigma](double* x, double* /*par*/) -> double {
            double qx = x[0];
            double q = sqrt(qx * qx + qy * qy);
            if (q < 0) {
                std::cerr << "Error: Negative q encountered!" << std::endl;
                return 0.0;
            }
            double FF_value = inputTF1->Eval(q);
            return FF_value * exp(-(qx - qx_prime) * (qx - qx_prime) / (2 * sigma * sigma));
        };
    
        TF1* f = new TF1("gaussian", gaussianLambda, qqx_min, qqx_max, 0);
        double result = 0.0;
    
        try {
            result = f->Integral(qqx_min, qqx_max, 1e-12);
        } catch (const std::exception& e) {
            std::cerr << "Error during integration: " << e.what() << std::endl;
        }
    
        delete f; // Avoid memory leaks
        std::cout << "Result of FF_wRes: " << result << std::endl;
        return result;
    }
    
/*    
    // Integral over qx -> new form factor with smearing
    static double FF_wRes(double qy, double qx_prime, double sigma, TF1* inputTF1) 
    {
        double qqx_min = 0, qqx_max = 20;
        TF1 *f = new TF1("f", guassian, qqx_min, qqx_max, 3);  
        double q = qx_prime*qx_prime + qy*qy;
        double FF_val = inputTF1->Eval(q);
        f->SetParameters(qy, qx_prime, sigma, FF_val);
        return f->Integral(qqx_min, qqx_max, 1e-12); // F(qx',qy) [GeV]
    }*/
    // Integrand for wedge cut
    static double FF_cut_wRes(double *x, double *par, TF1* inputTF1)
    {
        double t = par[0], sigma = par[1];
        double phi = x[0];  // integration variable
        double q = sqrt(t), qx_prime = q*sin(phi), qy = q*cos(phi); 
        // call form factor with resolution
	    return FF_wRes(qy, qx_prime, sigma, inputTF1); // [GeV^2]
    }
    
  




// Functions
TF1 *getWedgeRes_fun_1D() const {return cutFFt2;}
TH1D *getWedgeRes_hist_1D() const {return hist1D;}

private:
        TF1* cutFFt2;
        TH1D* hist1D;

    };
    

#endif
