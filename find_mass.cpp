#include <iostream>
#include <cmath>

#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"


int main(int argc, char** argv){
  if (argc != 3){
    std::cout << "wrong number of arguments\n";
    std::cout << "Try find_mass /input/file /output/file\n";
    return -1;
  }
  
//Get D and A trees. Open output file.
  TFile * input_file = new TFile(argv[1]);
  TFile * output_file = new TFile(argv[2],"RECREATE");
    
//Make trees and histograms for the nuclei
  TTree * generator_tree = (TTree*)input_file->Get("genT");
  TH2D * his_P1_Mtf = new TH2D("P1_VS_Mtf","P1_VS_Mtf;P1;Mtf;counts", 800, 0, 800, 800, 400, 1200); // bin #, min1, min2, max
  TH2D * his_P1_Mtf_mean = new TH2D("P1_VS_Mtf_mean","P1_VS_Mtf_mean;P1;Mtf;counts", 1000, 0, 1000, 1000, 0, 1000); // bin #, min1, min2, max
  TH2D * his_theta_P1prime_q = new TH2D("theta_VS_P1prime_q","theta_VS_P1prime_q;P1prime/q;theta;counts", 20, 0, 1.3, 30, 0, 1.5); // bin #, min1, min2, max

// Define Variables
  double X_b, Q_2, q_vec[3], P1_prime[3], P1_prime_mag, q_mag, weight, P1;  // Define variables from generator.cpp
  double M_n = 0.93827231;	                                // Mass of exiting proton; Gev
  
// Get branch (values) from generator.cpp   
  generator_tree->SetBranchAddress("X_b",&X_b);
  generator_tree->SetBranchAddress("Q_2",&Q_2);
  generator_tree->SetBranchAddress("P1_prime_mag",&P1_prime_mag);
  generator_tree->SetBranchAddress("q_mag",&q_mag);
  generator_tree->SetBranchAddress("weight",&weight);
  generator_tree->SetBranchAddress("P1",&P1);
  generator_tree->SetBranchAddress("P1_prime",P1_prime);
  generator_tree->SetBranchAddress("q_vec",q_vec);

// Loop over all entries
  for (int i = 0; i < generator_tree->GetEntries(); i++){
    generator_tree->GetEvent(i);
    if (X_b < 1.) continue;
    if ((P1_prime_mag/q_mag) > 0.96) continue;
    if ((P1_prime_mag/q_mag) < 0.62) continue;    
    double q_dot_pprime = (P1_prime[0]*q_vec[0]+P1_prime[1]*q_vec[1]+P1_prime[2]*q_vec[2]);
    double cos_theta_P1_prime_q = (q_dot_pprime/(P1_prime_mag*q_mag)); // theta (angle) between P1_prime and q
    if (fabs(cos_theta_P1_prime_q) > 1) continue;
    double theta_P1_prime_q = acos(cos_theta_P1_prime_q);
    if (theta_P1_prime_q > 25.*(2.*M_PI/360.)) continue;
    double w = Q_2 / (2 * M_n * X_b);	           // Energy transfered to nucleon from photon; Gev
    double E1_prime = sqrt(P1_prime_mag*P1_prime_mag + M_n*M_n);  // final Energy of ejected nucleon; Gev
    double Mtf = sqrt(-Q_2 + 4*M_n*M_n - P1_prime_mag*P1_prime_mag + E1_prime*E1_prime + 4*M_n*(w - E1_prime) + 2*q_dot_pprime - 2*E1_prime*w);
    his_P1_Mtf->Fill(P1*1000,Mtf*1000, weight);
    his_theta_P1prime_q->Fill((P1_prime_mag/q_mag),theta_P1_prime_q, weight);
  }
  
// Sum of squares (error)
  his_P1_Mtf->Sumw2();
  his_theta_P1prime_q->Sumw2();

// Find mean for every x point and replot
  double epsilon = 50;
  double init_point = 0;
  double x[20];
  double y[20];
  double ex[20];
  double ey[20];
  


for (int round = 0; round < 19; round++){ 
  TH1D * proj_Mtf = his_P1_Mtf->ProjectionY(" ", init_point + epsilon*round, init_point + epsilon*(round+1.));
  double Mtf_bin_mean = proj_Mtf->GetMean(); // find mean given in bin graph
  double Mtf_bin_std = proj_Mtf->GetStdDev();
  //proj_Mtf->Fit("gaus");
  x[round] = round*epsilon + 25;
  y[round] = Mtf_bin_mean;
  ex[round] = 0.01;
  ey[round] = Mtf_bin_std;
  his_P1_Mtf_mean->Fill(round*epsilon+25,Mtf_bin_mean);
 }

 TCanvas *c1 = new TCanvas("c1","A Simple Graph Example",900,900);

  // TGraphErrors * error = new TGraphErrors(ADC_3_x.size(),&ADC_3_x[0],&ADC_3_y[0]);

  
  TGraphErrors* Mtf_error = new TGraphErrors(20, x, y, ex, ey);
   Mtf_error->SetTitle("");      // graph title
   Mtf_error->SetLineColor(kBlue);                                
   Mtf_error->SetLineWidth(1);                                   
   Mtf_error->SetMarkerStyle(48);                                // rad style for point
   Mtf_error->SetMarkerSize(1);                                // point size
   Mtf_error->SetMarkerColor(46);                                // cool color points
   Mtf_error->GetYaxis()->SetTitle("M_tf");  // y-axis title
   Mtf_error->GetYaxis()->SetTitleOffset(1.3);                   // y-axis title offset from axis
   Mtf_error->GetYaxis()->CenterTitle(true);                     // center on axis
   Mtf_error->GetXaxis()->SetTitle("P1");                     // X-axis title
   Mtf_error->GetXaxis()->CenterTitle(true);                     // Center on axis
   Mtf_error->SetMinimum(0);                                     // y-axis minimum
   Mtf_error->SetMaximum(1100);                                  // y-axis maximum
   Mtf_error->Draw("AP");                                        // Draw on canvas
   


  
    input_file->Close();
    his_P1_Mtf->Write();
    his_P1_Mtf_mean->Write();
    his_theta_P1prime_q->Write();
    Mtf_error->Write();
    // proj_Mtf->Write();
    output_file->Close();
  
  return 0;
}
