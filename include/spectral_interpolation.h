#ifndef __SPECTRAL_INTERPOLATION_H__
#define __SPECTRAL_INTERPOLATION_H__

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>
#include <cstdlib>
#include <list>
#include <stdio.h>
using namespace std;


#include "TH2D.h"

class spec_info{
 private:
  void read_file();
  TH2D *neutron_spec; //size of CSV file 200 rows, 1000 columns
  TH2D *proton_spec;  //size of CSV file 200 rows, 1000 columns
  
 public:
  double spec_find(double k_test, double es_test, int code);
  double spec_find_fake(double k_test, double es_test, int code);

  spec_info();  
  ~spec_info();

};

#endif
