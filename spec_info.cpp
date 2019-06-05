
#include "spec_info.h"
#include "TH2D.h"

#include <fstream>
#include <iostream>
#include <stdio.h>
#include <cmath>
using namespace std;

// k = momentum; es = energy (for the corresponding spectral value you want to generate)

spec_info::spec_info(){read_file();}
spec_info::~spec_info(){}
void spec_info::read_file(){
// Parse the CSV of the experimental proton or neutron spectral function data
  neutron_spec = new TH2D("neutron spec", "neutron spec", 200, 0.025, 9.975, 1000, 0.5, 999.5);
  proton_spec = new TH2D("proton spec", "proton spec", 200, 0.025, 9.975, 1000, 0.5, 999.5);

  std::ifstream  data_p("SSP4.csv");
// check if file opens
    if ( ! data_p.is_open() ) {                 
      cout <<" Failed to open" << endl;}     
    std::string line_p;
    std::vector<std::vector<std::string> > parsedCsv_p;
    while(std::getline(data_p,line_p))
    {
        std::stringstream lineStream(line_p);
        std::string cell;
        std::vector<std::string> parsedRow_p;
        while(std::getline(lineStream,cell,','))
        {
            parsedRow_p.push_back(cell);
        }
        parsedCsv_p.push_back(parsedRow_p);
    }
    std::ifstream  data_n("SSN4.csv");
// check if file opens
    if ( ! data_n.is_open() ) {                 
      cout <<" Failed to open" << endl;}     
    std::string line_n;
    std::vector<std::vector<std::string> > parsedCsv_n;
    while(std::getline(data_n,line_n))
    {
        std::stringstream lineStream(line_n);
        std::string cell;
        std::vector<std::string> parsedRow_n;
        while(std::getline(lineStream,cell,','))
        {
            parsedRow_n.push_back(cell);
        }
        parsedCsv_n.push_back(parsedRow_n);
    }


  for(int i = 0; i < 200; i++){
    for(int j = 0; j < 1000; j++){
      neutron_spec->Fill(0.025 + (0.05*i), 0.5 + j); // fill in the x, y values at the i,j index (formulas specific to file given)
      neutron_spec->SetBinContent(neutron_spec->GetXaxis()->FindBin(0.025 + (0.05*i)),neutron_spec->GetYaxis()->FindBin(0.5 + j),atof(parsedCsv_n[i][j].c_str()));
      
      proton_spec->Fill(0.025 + (0.05*i), 0.5 + j); // fill in the x, y values at the i,j index (formulas specific to file given)
      proton_spec->SetBinContent(proton_spec->GetXaxis()->FindBin(0.025 + (0.05*i)),proton_spec->GetYaxis()->FindBin(0.5 + j),atof(parsedCsv_p[i][j].c_str()));}}
  }

 
double spec_info::spec_find(double k_test, double es_test, int code)
{
// If inputs are out of bounds, spec_fun = 0. Energy and momentum are in units inverse fermi
  if (es_test > 999.5 or es_test < 0.5){
    return 0;}
  if (k_test > 9.975 or k_test < 0.025){
    return 0;}
  
//Interpolate based on neutron vs proton
   if (code == 2112){
      return neutron_spec->Interpolate(k_test, es_test);}
   if (code == 2212) {
      return proton_spec->Interpolate(k_test, es_test);}
}

double spec_info::spec_find_fake(double k_test, double es_test, int code)
{
  double k_in_MeV = k_test * 197.3;
  // Return prob(k_test) * Prob(es_test)
  double E_exponent = (es_test - 20. - k_in_MeV*k_in_MeV/(2*3.016*931.49))/0.3; // 0.3 MeV width
  double k_exponent = k_in_MeV/100.; // 100 MeV corresponds to kF approx 200

  // std::cerr << k_exponent << " " << E_exponent << " " << exp(-0.5*(k_exponent * k_exponent + E_exponent*E_exponent)) << "\n";
  return exp(-0.5*(k_exponent * k_exponent + E_exponent*E_exponent));
}

