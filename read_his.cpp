
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <math.h> 

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
  TH1D * his_X_b =  new TH1D("X_b" ,"X_b [unitless]",40,0.5,1.5);
  TH1D * his_E_miss = new TH1D("E_star","E_star [GeV]",36,0,1);
  TH1D * his_P_k_mag = new TH1D("P_k_mag","P_k_mag [GeV]",36,0,12);
  
  TH2D * his_P1_Xb = new TH2D("P1_Vs_X_b","P1_V_S_X_b;;X_b;P1_mag;counts", 36, 0.5, 0, 36, 1.5, 12); // bin #, min1, min2, max
 

// Sum of squares (error)
  his_Q_2 ->Sumw2();
  his_X_b ->Sumw2();
  his_E_miss ->Sumw2();
  his_P_k_mag ->Sumw2();
  his_P1_Xb ->Sumw2();

  cerr<<"Histograms and Trees successfully created\n";

//Define variables needed for histograms
  Double_t Q_2, X_b, E_miss, P_k_vec[3], weight, P1_mag;
  
//Set proton and neutron numbers
  Int_t neunum = 2112, pronum = 2112;

  //Set addresses for D
  TreeH->SetBranchAddress("Q_2",&Q_2);
  TreeH->SetBranchAddress("X_b",&X_b);
  TreeH->SetBranchAddress("E_miss",&E_miss);
  TreeH->SetBranchAddress("P_k_vec",P_k_vec);
  TreeH->SetBranchAddress("weight",&weight);
  TreeH->SetBranchAddress("P1_mag",&weight);

  
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
  }
    cerr<<"Finished filling histogram\n";


    DataH->Close();
    his_Q_2->Write();
    his_X_b->Write();
    his_E_miss->Write();
    his_P_k_mag->Write();
    his_P1_Xb->Write();
    fo->Close();
  cerr<< argv[3]<<" has been completed. \n\n\n";
  
  return 0;
}
