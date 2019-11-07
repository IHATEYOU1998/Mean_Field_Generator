#include <iostream>
#include <cmath>

#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TF1.h"


int main(int argc, char** argv){
  if (argc != 3){
    std::cout << "wrong number of arguments\n";
    std::cout << "Try find_mass /input/file /output/file\n";
    return -1;
  }

  
//Get D and A trees. Open output file.
  TFile * input_file = new TFile(argv[1]);
  TFile * output_file = new TFile(argv[2],"RECREATE");


  double total_bins = 300.;  
//Make trees and histograms for the nuclei
  TTree * generator_tree = (TTree*)input_file->Get("genT");
  TH2D * his_P1_Mtf = new TH2D("P1_VS_Mtf","P1_VS_Mtf;P1;Mtf", total_bins, 0, 300, total_bins, 860, 1040); // bin #, min1, min2, max
  TH2D * his_P1_Mtf_mean = new TH2D("P1_VS_Mtf_mean","P1_VS_Mtf_mean;P1;Mtf", total_bins, 0, 300, total_bins, 860, 1040); // bin #, min1, min2, max
  TH2D * his_theta_P1prime_q = new TH2D("theta_VS_P1prime_q","theta_VS_P1prime_q;P1prime/q;theta", 1000, 0, 1.5, 1000, 0, 200); // bin #, min1, min2, max
  TH2D * his_E1_theta_P1q = new TH2D("E1_VS_theta_P1q","graph;Theta Between P1 and q: Radians;E1: GeV", total_bins, 0, 3.5, total_bins, 0.8, 1.); // bin #, min1, min2, max
  TH2D * his_P1_theta_P1prime_Pk = new TH2D("P1_VS_theta_P1_prime_Pk","graph2;P1(phi)-Pk(phi): Radians; P1: GeV", total_bins, 0, 360, total_bins, 0., 0.5); // bin #, min1, min2, max

  
// Define Variables
  double X_b, Q_2, q_vec[3], P1_prime_vec[3], P1_prime_mag, q_mag, weight, P1_mag, w;  // Define variables from generator.cpp
  double P1_vec[3];  // Define variables from generator.cpp
  double M_n = 0.93827231;	                                                       // Mass of exiting proton (close to neutron); Gev
  double Phi_k;
  
// Get branch (values) from generator.cpp   
  generator_tree->SetBranchAddress("X_b",&X_b);
  generator_tree->SetBranchAddress("Q_2",&Q_2);
  generator_tree->SetBranchAddress("P1_prime_mag",&P1_prime_mag);
  generator_tree->SetBranchAddress("q_mag",&q_mag);
  generator_tree->SetBranchAddress("weight",&weight);
  generator_tree->SetBranchAddress("P1_mag",&P1_mag);
  generator_tree->SetBranchAddress("w",&w);
  generator_tree->SetBranchAddress("P1_prime_vec",P1_prime_vec);
  generator_tree->SetBranchAddress("q_vec",q_vec);
  generator_tree->SetBranchAddress("P1_vec",P1_vec);
  generator_tree->SetBranchAddress("Phi_k",&Phi_k);

  
// Loop over all entries
  for (int i = 0; i < generator_tree->GetEntries(); i++){   
    generator_tree->GetEvent(i);
    
    TVector3 P1_prime_TVec(P1_prime_vec[0], P1_prime_vec[1], P1_prime_vec[2]); //turn to TVector
    TVector3 q_TVec(q_vec[0], q_vec[1], q_vec[2]);                             //turn to TVector
    TVector3 P1_TVec(P1_vec[0], P1_vec[1], P1_vec[2]);                         //turn to TVector

    his_P1_theta_P1prime_Pk->Fill((Phi_k - P1_prime_TVec.Phi())*360/(2*M_PI), P1_mag);
    
    TLorentzVector test_P1_prime;
    test_P1_prime.SetPxPyPzE(P1_prime_vec[0],P1_prime_vec[1],P1_prime_vec[2],sqrt(P1_prime_TVec.Mag2() + M_n*M_n));
    TLorentzVector test_q;
    test_q.SetPxPyPzE(q_vec[0], q_vec[1], q_vec[2], w);
    TLorentzVector test_pair;
    test_pair.SetPxPyPzE(0, 0, 0, 2*M_n);
        
    Double_t theta_P1_prime_q = P1_prime_TVec.Angle(q_TVec);                   //give angle between vectors, Radians
    double E1_prime = sqrt(P1_prime_mag*P1_prime_mag + M_n*M_n);  // final Energy of ejected nucleon; Gev
    Double_t theta_P1_q = P1_TVec.Angle(q_TVec);                   //give angle between vectors, Radians
    his_E1_theta_P1q->Fill(theta_P1_q, (E1_prime - w));
    if (theta_P1_prime_q > 25.*(2.*M_PI/360.)) continue;
    double Mtf = sqrt((test_q+test_pair-test_P1_prime).Mag2());
    
    if (X_b < 1.15) continue;   
    if ((P1_prime_mag/q_mag) > 0.96) continue;
    if ((P1_prime_mag/q_mag) < 0.62) continue;
    if (Mtf < 0.) continue;
    his_P1_Mtf->Fill(P1_mag*1000,Mtf*1000, weight);
    his_theta_P1prime_q->Fill((P1_prime_mag/q_mag), theta_P1_prime_q*(180/M_PI), weight);   
  }
  
  
// Sum of squares (error)
  his_P1_Mtf->Sumw2();
  his_theta_P1prime_q->Sumw2();
  his_E1_theta_P1q->Sumw2();
  his_P1_theta_P1prime_Pk->Sumw2();

// Find mean for every x point and replot
  double section_width = 15.;
  double total_sections = total_bins/section_width;
  double x[(int) total_sections];
  double y[(int) total_sections];
  double ex[(int) total_sections];
  double ey[(int) total_sections];
  
// Loop through to find mean of sections in graph
for (int round = 0.; round < total_sections; round++){ 
  TH1D * proj_Mtf = his_P1_Mtf->ProjectionY(":)", (total_bins/total_sections)*round, (total_bins/total_sections)*(round+1.));
  double Mtf_bin_mean = proj_Mtf->GetMean(); // find mean given in bin graph
  double Mtf_bin_std = proj_Mtf->GetStdDev();
  x[round] = (round+0.5)*section_width;  
  y[round] = Mtf_bin_mean;
  ex[round] = 10;//section_width/2.;
  ey[round] = 10;//Mtf_bin_std;
  his_P1_Mtf_mean->Fill((round+0.5)*section_width, Mtf_bin_mean);
  his_P1_Mtf_mean->SetBinError(x[round], 10);
  std::cout << "std:  " <<  Mtf_bin_std << "        mean:  " << Mtf_bin_mean << "\n"; 
  proj_Mtf->Reset("ICESM");
 }
his_P1_Mtf_mean->Sumw2();
his_P1_Mtf_mean->Draw();
his_P1_Mtf_mean->SetMarkerSize(1);

 
// Error Bar Graph of Mtf Mean
//c1 = new TCanvas("c1","A Simple Graph with error bars",200,10,700,500);
TGraphErrors* Mtf_error = new TGraphErrors(total_sections, x, y, ex, ey);
// Make Pretty and Draw   
   Mtf_error->SetTitle("Mean Field Contribution");            // graph title
   Mtf_error->SetLineColor(4);//kBlue);                                
   Mtf_error->SetLineWidth(0);                                   
   Mtf_error->SetMarkerStyle(8);//48);                        // rad style for point
   Mtf_error->SetMarkerSize(1);                               // point size
   Mtf_error->SetMarkerColor(4);//46);                        // cool color points
   Mtf_error->GetYaxis()->SetTitle("Missing Mass [MeV]");
   Mtf_error->GetYaxis()->SetLabelSize(0.03);     // y-axis title
   Mtf_error->GetYaxis()->SetTitleOffset(1.3);                // y-axis title offset from axis
   Mtf_error->GetYaxis()->CenterTitle(true);                  // center on axis
   Mtf_error->GetXaxis()->SetTitle("Missing Momentum [MeV]"); // X-axis title
   Mtf_error->GetXaxis()->SetLabelSize(0.03);     // y-axis title
   Mtf_error->GetXaxis()->CenterTitle(true);                  // Center on axis
   Mtf_error->SetMinimum(860);                                  // y-axis minimum
   Mtf_error->SetMaximum(1040);                               // y-axis maximum
   Mtf_error->Draw("APE");                                     // Draw on canvas

  
// Draw a horizontal line (on error bar graph) for the mass of a proton
TF1 function_3_2_1("Linear law","[0]", 1);// Let's make the funcion line nicer//#7
function_3_2_1.SetParLimits(0,1000.*M_n,1000.*M_n);
function_3_2_1.SetLineColor(kRed);  function_3_2_1.SetLineStyle(10);// Fit it to the graph and draw it//#8
Mtf_error->Fit(&function_3_2_1);
   
 
// Cleanup and Close
    input_file->Close();
    his_P1_Mtf->Write();
    his_P1_Mtf_mean->Write();
    his_theta_P1prime_q->Write();
    Mtf_error->Write();
    his_E1_theta_P1q->Write();
    his_P1_theta_P1prime_Pk->Write();
    output_file->Close();
  return 0;
}
