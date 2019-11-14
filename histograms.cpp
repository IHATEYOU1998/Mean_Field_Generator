
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TVector3.h"

using namespace std;

// given by Andrew on Friday April 12th 2019; modified after

int main(int argc, char ** argv){

  if( argc != 3){
    
    cerr<<"Wrong number of arguments. Instead try:\n\t"
	<< "src_reader /path/to/input/nucleus/file /path/to/output/file \n";

    return -1;
    
  }
  
  //Get D and A trees. Open output file.
  TFile * DataH = new TFile(argv[1]);
  TFile * fo = new TFile(argv[2],"RECREATE");
  
  cerr<<"File has been open from: "<<argv[1]<<"\n";
  
//Make trees and histograms for the nuclei
  TTree * TreeH = (TTree*)DataH->Get("genT");
  TH1D * his_Q_2 = new TH1D("Q_2","Q_2 [GeV^2];Counts",36,1.,5.);
  TH1D * his_X_b =  new TH1D("Xb" ,"X_b [unitless]",40,0.5,1.5);
  TH1D * his_E_miss = new TH1D("E_star","E_star [GeV]",36,0,1);
  TH1D * his_P_k_mag = new TH1D("P_k_mag","P_k_mag [GeV]",36,0,12);
  TH2D * his_P1_Xb = new TH2D("P1_VS_Xb","P1_VS_Xb; X_b; P1_mag; counts", 300, 0.5, 1.5, 300, 0., 0.25); // bin #, min1, min2, ] x,y
  TH2D * his_Phi_k_Phi_P1_prime = new TH2D("Phi(scattered_electron)_VS_Phi(ejected_nucleon)", "Phi(scattered_electron)_VS_Phi(ejected_nucleon); Phi_k; Phi_P1_prime; counts", 300, 0, 2*M_PI, 300, -3, 3);
  TH2D * his_P_k_P1_prime = new TH2D("P(scattered_electron)_VS_P(ejected_nucleon)", "P(scattered_electron)_VS_P(ejected_nucleon); P_k; P1_prime; counts", 300, 1, 3.5, 300, 1.5, 4);
   

// Sum of squares (error)
  his_Q_2 ->Sumw2();
  his_X_b ->Sumw2();
  his_E_miss ->Sumw2();
  his_P_k_mag ->Sumw2();
  his_P1_Xb ->Sumw2();
  his_Phi_k_Phi_P1_prime->Sumw2();
  his_P_k_P1_prime->Sumw2();

  cerr<<"Histograms and Trees successfully created\n";

//Define variables needed for histograms
  Double_t Q_2, X_b, E_miss, P_k_vec[3], weight, P1_mag, Phi_P1_prime, Phi_k, P1_prime_mag;
  
//Set proton and neutron numbers
  Int_t neunum = 2112, pronum = 2112;

//Set addresses for D
  TreeH->SetBranchAddress("Q_2",&Q_2);
  TreeH->SetBranchAddress("X_b",&X_b);
  TreeH->SetBranchAddress("E_miss",&E_miss);
  TreeH->SetBranchAddress("P_k_vec",P_k_vec);
  TreeH->SetBranchAddress("weight",&weight);
  TreeH->SetBranchAddress("P1_mag",&P1_mag);
  TreeH->SetBranchAddress("Phi_P1_prime",&Phi_P1_prime);   // Angle of ejected nucleon; Spherical: to Z-axis
  TreeH->SetBranchAddress("Phi_k",&Phi_k);   // Angle of ejected nucleon; Spherical: to Z-axis
  TreeH->SetBranchAddress("P1_prime_mag",&P1_prime_mag);   // Angle of ejected nucleon; Spherical: to Z-axis

    
//Loop over TTree
  for(int i = 0; i < TreeH->GetEntries(); i++){
    TreeH->GetEntry(i);
    TVector3 P_k_vecT(P_k_vec[0],P_k_vec[1],P_k_vec[2]);
    double theta_k = P_k_vecT.Theta();
    double phi_k = P_k_vecT.Phi();
    double P_k_mag = P_k_vecT.Mag();
    his_Q_2->Fill(Q_2, weight);
    his_X_b->Fill(X_b, weight);
    his_E_miss->Fill(E_miss, weight);
    his_P_k_mag->Fill(P_k_mag, weight);
    his_P1_Xb->Fill(X_b, P1_mag, weight);
    his_Phi_k_Phi_P1_prime->Fill(Phi_k, Phi_P1_prime,  weight);
    his_P_k_P1_prime->Fill(P1_prime_mag, P_k_mag);
  }
    cerr<<"Finished filling histogram\n";


    DataH->Close();
    his_Q_2->Write();
    his_X_b->Write();
    his_E_miss->Write();
    his_P_k_mag->Write();
    his_P1_Xb->Write();
    his_Phi_k_Phi_P1_prime->Write();
    his_P_k_P1_prime->Write();
    fo->Close();
  cerr<< argv[3]<<" has been completed. \n\n\n";
  
  return 0;
}
