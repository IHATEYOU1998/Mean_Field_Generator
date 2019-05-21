
#include <iostream>

#include "TVector3.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"

#include "spec_info.h"
#include "Cross_Sections.h"


// Physical Constants
  const double M_e = 0.0005109906;                 // Mass of electron; GeV
  const double M_A = 4.0026032 * 0.93149 - M_e*2;  // Mass of initial nucleus; Helium-4  Gev
  const double BE_A = 0.028295673;                 // Binding energy of initial nucleus; Helium-4; from kaeri; Gev
  const double P_b_x = 0;	                   // initial electron momentum 'x'; Gev  (if non-zero, need to change code)
  const double P_b_y = 0;      	                   // initial electron momentum 'y'; Gev (if non-zero, need to change code)

// Randomly generated variabe's range
  const double xb_max = 1.5;      // Maximum in X_b; unitless
  const double xb_min = 0.5;      // Minimum in X_b; unitless
  const double Q2_max = 4;        // Maximum in Q^2; GeV^2
  const double Q2_min = 1;        // Minimum in Q^2; GeV^2
  const double P1x_max = 0.4;     // Maximum of initial nucleon's x-momentum; GeV
  const double P1x_min = -0.4;    // Minimum of initial nucleon's x-momentum; GeV
  const double P1y_max = 0.4;     // Maximum of initial nucleon's y-momentum; GeV
  const double P1y_min = -0.4;    // Minimum of initial nucleon's y-momentum; GeV
  const double P1z_max = 0.4;     // Maximum of initial nucleon's z-momentum; GeV
  const double P1z_min = -0.4;    // Minimum of initial nucleon's z-momentum; GeV
  const double phi_max = 2*M_PI;  // Maximum in X_b; radians
  const double phi_min = 0;       // Minimum in X_b; radians


int main(int argc, char ** argv){
  
// Print error message if not enough / too much arguments inputted. Note: for 3 input variables, argc != 4.
// arguments: 0: unimportant; 1: initial electron z momentum (GeV); 2: number of events to record; 3: where to store data (../DATA/gen_file.root);
  if( argc != 4)
    {
      std::cerr << "Wrong number of arguments. Instead try: \n\t"
	   << "./gen_root Energy_of_beam number_of_events output_file ..DATA/gen_file.root\n";
      
      return -1;
    }
  
   
// Input variables
  double P_b_z = atof(argv[1]);      // initial electron momentum 'z'; usually between 4 - 12 Gev?; P_z = E_b;
  int nEvents = atoi(argv[2]);       // number of times to loop and find differential cross section

// Define variables for later use
  double M_A_1;                      // Mass of final nucleus; Tritium or He-3; Gev
  string Spec_CSV_file;              // Define the spectral function csv file for later on
  double M_n;                        // Define mass of exiting nucleon used later on
  bool isProton;                     // Define boolean of whether we have a proton or neutron
  double BE_A_1;                     // Binding energy of recoil nucleus

  
// TTree to save values during the loop
  TFile * outfile = new TFile(argv[3], "RECREATE");
  TTree * outtree = new TTree("genT", "Generator Tree");
  double weight, Q_2, X_b, E_miss, P1_prime[3], P_k[3], nucleon_type, q_mag, q_vec[3], P1_prime_mag, P1;
  double theta_k, theta_P1_prime, Phi_P1_prime, P_k_mag, E_b, P_b_vec[3], P_A_1, P_A_1_vec[3], theta_A_1, Phi_A_1, P1_vec[3], E_A_1;
  
//Variables we are saving: "name", &variable, "name with elements"/double
  outtree->Branch("weight",&weight,"weight/D");                    // The weighed differential cross section; proton or neutron
  outtree->Branch("X_b",&X_b,"X_b/D");                             // Bjorken scaling variable X_b, randomly selected
  outtree->Branch("Q_2",&Q_2,"Q_2/D");                             // Q^2, randomly selected
  outtree->Branch("E_miss",&E_miss,"E_miss/D");                    // Excited energy of nucleus, goes in spectral function, E*
  outtree->Branch("P1_prime",P1_prime,"P1_prime[3]/D");         // Momentum of exiting nucleus <x, y, z>
  outtree->Branch("P_k",P_k,"P_k[3]/D");                           // Momentum of scattered electron <x, y, z>
  outtree->Branch("nucleon_type",&nucleon_type,"nucleon_type/I");  // True if proton, false if neutron

  outtree->Branch("q_mag",&q_mag,"q_mag/D");                              // Momentum of scattered electron <x, y, z>
  outtree->Branch("q_vec",q_vec,"q_vec[3]/D");         // Momentum of exiting nucleus <x, y, z>
  
  outtree->Branch("E_A_1",&E_A_1,"E_A_1k/D");                             // Momentum of scattered electron <x, y, z>
  outtree->Branch("theta_k",&theta_k,"theta_k/D");                        // Momentum of scattered electron <x, y, z>
  outtree->Branch("theta_P1_prime",&theta_P1_prime,"theta_P1_prime/D");   // Momentum of scattered electron <x, y, z>
  outtree->Branch("Phi_P1_prime",&Phi_P1_prime,"Phi_P1_prime/D");         // Momentum of scattered electron <x, y, z>
  outtree->Branch("P_k_mag",&P_k_mag,"P_k_mag/D");                        // Momentum of scattered electron <x, y, z>
  outtree->Branch("P1_prime_mag",&P1_prime_mag,"P1_prime_mag/D");      // Momentum of scattered electron <x, y, z>
  outtree->Branch("E_b",&E_b,"E_b/D");                                    // Momentum of scattered electron <x, y, z>
  outtree->Branch("P_b_vec",P_b_vec,"P_b_vec[3]/D");                      // Momentum of scattered electron <x, y, z>
  outtree->Branch("P_A_1",&P_A_1,"P_A_1/D");                              // Momentum of scattered electron <x, y, z>
  outtree->Branch("P_A_1_vec",P_A_1_vec,"P_A_1_vec[3]/D");                // Momentum of scattered electron <x, y, z>
  outtree->Branch("theta_A_1",&theta_A_1,"theta_A_1/D");                  // Momentum of scattered electron <x, y, z>
  outtree->Branch("Phi_A_1",&Phi_A_1,"Phi_A_1/D");                        // Momentum of scattered electron <x, y, z>
  outtree->Branch("P1",&P1,"P1/D");                                    // Momentum of scattered electron <x, y, z>
  outtree->Branch("P1_vec",P1_vec,"P1_vec[3]/D");                      // Momentum of scattered electron <x, y, z>

  
// Initialize classes
  spec_info Spec_func;                  // Access the spectral function interpolator spec_find from the spec_info class
  Cross_Sections Sig(onshell, dipole);  // Access cross section sigma_en from the Cross_Sections class (given by andrew denniston)
  TRandom3 myRand(0);                   // Get random number

  
// loop to find the differential cross section "nEvents" amount of times.
for (int event = 0 ; event < nEvents ; event++){
  
// Reset tree variable to zero after every loop
  weight = Q_2 = X_b = E_miss = nucleon_type = q_mag = 0;
  theta_k = theta_P1_prime = Phi_P1_prime = P_k_mag = P1_prime_mag = E_b = P_A_1 = theta_A_1 = Phi_A_1 = P1 = 0;
  memset(P1_prime, 0, sizeof(P1_prime));
  memset(P_k, 0, sizeof(P_k));
  memset(P1_vec,0, sizeof(P1_vec));
  memset( P_b_vec, 0, sizeof(P_b_vec));
  memset(P_A_1_vec, 0, sizeof(P_A_1_vec));
  memset(q_vec, 0, sizeof(q_vec));

// Define Random variables
  double Phi_k = phi_min + myRand.Rndm()*(phi_max-phi_min);         // electron scatter angle: radians. Random value between 0 and 2*pi
  Q_2 = Q2_min + myRand.Rndm()*(Q2_max-Q2_min);              // Mev^2. Random value between 0 Gev^2 and 5 Gev^2
  X_b = xb_min + myRand.Rndm()*(xb_max-xb_min);              // Bjorken scaling variable: unitless; Random  between 0.5 and 1.5
  double P1_x = P1x_min + myRand.Rndm()*(P1x_max-P1x_min);  // initial nucleon momentum 'x'; Mev; Random between [-0.4,0.4] GeV
  double P1_y = P1y_min + myRand.Rndm()*(P1y_max-P1y_min);  // initial nucleon momentum 'y'; Mev; Random between [-0.4,0.4] GeV
  double P1_z = P1z_min + myRand.Rndm()*(P1z_max-P1z_min);  // initial nucleon momentum 'z'; Mev; Random between [-0.4,0.4] GeV
  int n_or_p = 0;//(int)(myRand.Rndm() + .5);                    // A random number 0 (proton) or 1 (neutron)

  
// Choose the proton or neutral spectral function set (50:50 chance)
  if (n_or_p == 0){
    Spec_CSV_file = "SSP4.csv";         // experimental proton spectral function values
    isProton = true;                    // input needed in sigma_eN
    M_n = 0.93827231;	                // Mass of exiting proton; Gev
    nucleon_type = 2212;                // As specified by some random institute for a proton
    M_A_1 = 3.0160493 * 0.93149 - M_e;  // Mass of final nucleus; Tritium; Gev
    BE_A_1 = 0.008481821;}              // Binding energy of H-3; from kaeri; Gev
  else{
    Spec_CSV_file = "SSN4.csv";           // experimental neutron spectral function values
    isProton = false;                     // input needed in sigma_eN
    M_n = 0.93957;                        // Mass of exiting neutron; Gev
    nucleon_type = 2112;                  // As specified by some random institute for a neutron
    M_A_1 = 3.0160293 * 0.93149 - 2*M_e;  // Mass of final nucleus; He-3; Gev
    BE_A_1 = 0.007718058;}                // Binding energy of He-3; from kaeri; Gev

  
// Define variables of incoming electron
  P_b_vec[0] = P_b_x;
  P_b_vec[1] = P_b_y;
  P_b_vec[2] = P_b_z;
  double P_b = sqrt(P_b_x*P_b_x + P_b_y*P_b_y + P_b_z*P_b_z);  // Magnitude of incoming elecron momentum; Gev
  E_b = sqrt(P_b*P_b + M_e*M_e);                               // initial electron beam energy (about = P_z); Gev

  
// Define variables from electron scatter
  double w = Q_2 / (2 * M_n * X_b);	           // Energy transfered to nucleon from photon; Gev
  if (w > E_b) continue;                           // discontinue if unphysical
  double E_k = E_b - w;	                           // Energy of scattered electron; Gev
  P_k_mag = sqrt(E_k*E_k - M_e*M_e);               // Momentum of Scattered Electron; Gev
  q_mag = sqrt(Q_2 + w*w);                  // Momentum transfered to nucleon from photon; Gev
  double cos_theta_k = (P_b*P_b + P_k_mag*P_k_mag - q_mag*q_mag)/(2*P_b*P_k_mag);  // cos(angle) of scattered electron with Z-axis
  if (fabs(cos_theta_k) > 1.) continue;            // discontinue if unphysical
  theta_k = acos(cos_theta_k);                     // angle of scattered electron with z axis
  P_k[0] = P_k_mag*sin(theta_k)*cos(Phi_k);        // Scattered electron's X-momentum; Gev
  P_k[1] = P_k_mag*sin(theta_k)*sin(Phi_k);        // Scattered electron's Y-momentum; Gev
  P_k[2] = P_k_mag*cos(theta_k);                   // Scattered electron's Z-momentum; Gev
  q_vec[0] = -P_k[0];                            // X-momentum transfered to nucleon from photon; Gev
  q_vec[1] = -P_k[1];                            // Y-momentum transfered to nucleon from photon; Gev
  q_vec[2] = P_b_z - P_k[2];                     // Z-momentum transfered to nucleon from photon Gev

  
// Define nucleus momentum and enery variables
  P1_vec[0] = P1_x;
  P1_vec[1] = P1_y;
  P1_vec[2] = P1_z;
  P1 = sqrt(P1_x*P1_x + P1_y*P1_y + P1_z*P1_z);             // Momentum of initial nucleon to be ejected; Gev
  P1_prime[0] = P1_x + q_vec[0];                                      // final X-momentum of ejected nucleon; Gev   
  P1_prime[1] = P1_y + q_vec[1];                                      // final Y-momentum of ejected nucleon; Gev
  P1_prime[2] = P1_z + q_vec[2];                                      // final Z-momentum of ejected nucleon; Gev
  P1_prime_mag = sqrt(P1_prime[0]*P1_prime[0] + P1_prime[1]*P1_prime[1] + P1_prime[2]*P1_prime[2]);  // final Momentum of ejected nucleon; Gev
  double E1_prime = sqrt(P1_prime_mag*P1_prime_mag + M_n*M_n);  // final Energy of ejected nucleon; Gev
  P_A_1 = P1;                                                     // Momentum of final recoil nucleus; Gev
  P_A_1_vec[0] = -P1_x;
  P_A_1_vec[1] = -P1_y;
  P_A_1_vec[2] = -P1_z;
  double P_A = 0;                                                  // Momentum of initial nucleus; approx as zero; Gev
  double E_A = sqrt(P_A*P_A + M_A*M_A);                            // Energy of initial nucleus; Gev
  E_A_1 = E_A + w - E1_prime;                                     // Energy of final nucleus; Gev
  
  E_miss = M_n - (E1_prime - w);
  
  if ((E_A_1*E_A_1 - P_A_1*P_A_1) < M_A_1*M_A_1) continue;
  theta_P1_prime = acos(P1_prime[2]/P1_prime_mag);
  Phi_P1_prime = atan2(P1_prime[1],P1_prime[0]);
  theta_A_1 = acos(P_A_1_vec[2]/P_A_1);
  Phi_A_1 = atan2(P_A_1_vec[1],P_A_1_vec[0]);
  
  
// Variables that require outside calculation
  TVector3 P_k_TVec(P_k[0], P_k[1], P_k[2]);                             //turn to TVector
  TVector3 P1_prime_TVec(P1_prime[0], P1_prime[1], P1_prime[2]);     //turn to TVector
  double sigma = Sig.sigma_eN(E_b, P_k_TVec, P1_prime_TVec, isProton);  // cross section sigma_eN from Andrew and Jackson: Takes in (double Ebeam, TVector3 k, TVector3 p, bool isProton) returns in cm^2
  double spec = (pow(0.197345,3))*(0.001)*Spec_func.spec_find((double)P1/0.197345, (double)(E_miss)*1000., nucleon_type);  // Spectral function for (P_1,E*); Uses spec_find: P_1: inverse femtometers, E*: MeV; need spec in units 1/(geV^4);

// Find differential Cross Section
  double diff_cross = spec * sigma * (w / (2 * E_b * E_k * X_b));   

// Creating weighting function: uniform
  double normalize_range = 1./((xb_max-xb_min)*(Q2_max-Q2_min)*(phi_max-phi_min)*(P1x_max-P1x_min)*(P1y_max-P1y_min)*(P1z_max-P1z_min));  // (range of data points for spectral function)

// Divide Differential Cross Section by associated weighting
  weight = diff_cross/normalize_range;  // scaled by weighting factor;
  
// Only record physical values
  if (weight > 0){// and (M_A - E_A_1) < epsilon){   
    outtree->Fill();}  
  }
  
// Clean up
  outtree->Write();
  outfile->Close();
  return 0;
  
}
