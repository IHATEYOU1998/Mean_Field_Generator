#include <iostream>
#include <cmath>

#include <string>
#include <iostream>
#include<vector>

#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TGraphAsymmErrors.h"
#include <math.h>


int main(int argc, char** argv){
  if (argc != 3){
    std::cout << "wrong number of arguments\n";
    std::cout << "Try find_mass /input/file /output/file\n";
    return -1;
  }

  
//Get D and A trees. Open output file.
  TFile * input_file = new TFile(argv[1]);
  TFile * output_file = new TFile(argv[2],"RECREATE");


// Label constant variables to use in program:
  const double total_bins = 20.; // Histogram bins. Later on, you want to make sure you can divide by section width
  const int Mass_y_min = 860.;   //  
  const int Mass_y_max = 1040.;
  const int Mass_x_min = 0.;
  const int Mass_x_max = 1000.; // make sure it can work with total sections
  
//Make trees and histograms for the nuclei
  TTree * generator_tree = (TTree*)input_file->Get("genT");
  TH2D * his_P1_Mtf = new TH2D("P1_VS_Mtf","P1_VS_Mtf;P1;Mtf", total_bins, Mass_x_min, Mass_x_max, total_bins, Mass_y_min, Mass_y_max); // bin #, min1, min2, max
  TH2D * his_P1_Mtf_mean = new TH2D("P1_VS_Mtf_mean","P1_VS_Mtf_mean;P1;Mtf", total_bins, Mass_x_min, Mass_x_max, total_bins, Mass_y_min, Mass_y_max); // bin #, min1, min2, max
  TH2D * his_theta_P1prime_q = new TH2D("theta_VS_P1prime_q","theta_VS_P1prime_q;P1prime/q;theta", 1000, 0, 1.5, 1000, 0, 200); // bin #, min1, min2, max
  TH1D * his_Q2_weight = new TH1D("Q2_weight","Q2 [GeV^2];Counts", 20, 0., 4.1);
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
  
// Find mean for every x point and replot
  int section_width = 50.;
  int total_sections = Mass_x_max/section_width;
  double bin_per_section = total_bins/total_sections; // check to make sure divides
  double x[(int) total_sections];
  double y[(int) total_sections];
  double ex[(int) total_sections];
  double ey[(int) total_sections];
  TH1F *proj_histo[total_sections] = { NULL };

// Create Histograms for every section (projection section) of the graph
for (int round = 0.; round < total_sections; round++){
  std::ostringstream histogramNameStream;
  histogramNameStream << "projection_from_" << section_width*round << "_to_" << section_width*(round+1.);
  proj_histo[round] = new TH1F( histogramNameStream.str().c_str() , ":(" ,  total_bins, Mass_y_min, Mass_y_max);
 }

// Recreate graph with mean plotted
TCanvas *c1 = new TCanvas("c1","c1",900,900);// Mass_x_min, Mass_x_max, Mass_y_min, Mass_y_max);
for (int round = 0.; round < total_sections; round++){
  // project the Y axis
  TH1D * proj_Mtf = his_P1_Mtf->ProjectionY(":)", bin_per_section*round + 1, bin_per_section*(round+1.)); // name, first bin, last bin

  // fit the function and save
  proj_histo[round]->Add(proj_Mtf);
  proj_histo[round]->Fit("gaus", "QEM"); // fit the function to a gaussian

  // Write each fit    
  proj_histo[round]->Write();

  // Get values for Mass Graph
  double real_Mtf_bin_mean = proj_histo[round]->GetMean(); // find mean given in bin graph
  double real_Mtf_bin_std = proj_histo[round]->GetStdDev();
  int num_points = proj_histo[round]->GetEntries(); // gives number of entries as double -> turn to int
  x[round] = (round+0.5)*section_width;  
  y[round] = real_Mtf_bin_mean;
  ex[round] = 0.;
  ey[round] = real_Mtf_bin_std/sqrt(num_points - 1);
  his_P1_Mtf_mean->Fill((round+0.5)*section_width, real_Mtf_bin_mean);

  
  std::cout << "P1 Range: [" << section_width*round << ", " << section_width*(round+1.) << "]          "
	    << "no fit std:  " <<  real_Mtf_bin_std << "         no fit mean:  " <<  real_Mtf_bin_mean
	    << "          number of points in sample:  " << num_points << "  size: " << sizeof(proj_histo[round])/sizeof(proj_histo[round][0]) << "\n";
  proj_Mtf->Reset("ICESM");
 }
 
// Make mean graph pretty
  his_P1_Mtf_mean->SetMarkerStyle(8);  
  his_P1_Mtf_mean->SetMarkerSize(1);   
  his_P1_Mtf_mean->SetMarkerColor(4);  

 
// Error Bar Graph of Mtf Mean
   TCanvas *c2 = new TCanvas("c2","c2",900,900);// Mass_x_min, Mass_x_max, Mass_y_min, Mass_y_max);
   TGraphAsymmErrors *Mtf_error = new TGraphAsymmErrors(total_sections, x, y, ex, ex, ey, ey);
// Make Pretty and Draw   
   Mtf_error->SetTitle("Mean Field Contribution");            
   Mtf_error->SetLineColor(4);
   Mtf_error->SetMarkerStyle(8);
   Mtf_error->SetMarkerSize(1);                               
   Mtf_error->SetMarkerColor(4);
   Mtf_error->GetYaxis()->SetTitle("Missing Mass [MeV]");
   Mtf_error->GetYaxis()->SetLabelSize(0.03);     
   Mtf_error->GetYaxis()->SetTitleOffset(1.3);                
   Mtf_error->GetYaxis()->CenterTitle(true);                  
   Mtf_error->GetXaxis()->SetTitle("Missing Momentum [MeV]"); 
   Mtf_error->GetXaxis()->SetLabelSize(0.03);     
   Mtf_error->GetXaxis()->CenterTitle(true);                  
   Mtf_error->SetMinimum(860);                                
   Mtf_error->SetMaximum(1040);                               
   Mtf_error->SetFillColor(6);
   Mtf_error->SetFillStyle(3005);
   Mtf_error->Draw("AP");                                     
   c2->Update();
   

// Produce graphs and close
    input_file->Close();
    output_file->cd();
    his_P1_Mtf->Write();
    his_P1_Mtf_mean->Write();
    his_theta_P1prime_q->Write();
    Mtf_error->Write();
    c2->Update();
    output_file->Close();
    return 0;
}


// Add Gaussian Fit:
  /** get gaussian fit
  #include "TF1.h"
  double fit_Mtf_bin_mean = 0;
  double fit_Mtf_bin_std = 0.;
  TF1 *fit = (TF1 *)proj_histo[round]->GetFunction("gaus");
  if(fit != NULL){
    fit_Mtf_bin_mean = fit->GetParameter(1);
    fit_Mtf_bin_std = fit->GetParameter(2);}
  else{
    std:: cout << "No Data to Fit. Set: Mean = STD = 0 \n";} **/
