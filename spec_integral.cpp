#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
//#include <math.h>

#include <algorithm>
#include <string>
using namespace std;



int main(){

// Parse the CSV of the experimental proton or neutron spectral function data 
    std::ifstream  data("SSN4.csv");
    // check if file opens
    if ( ! data.is_open() ) {                 
      cout <<" Failed to open" << endl;
      }     
    std::string line;
    std::vector<std::vector<std::string> > parsedCsv;
    while(std::getline(data,line))
    {
        std::stringstream lineStream(line);
        std::string cell;
        std::vector<std::string> parsedRow;
        while(std::getline(lineStream,cell,','))
        {
            parsedRow.push_back(cell);
        }

        parsedCsv.push_back(parsedRow);
    }



double spec_integral = 0;

 for (int i=0; i<199; i++){
   for (int j=0; j<999; j++){
     double spec_av = (atof(parsedCsv[i][j].c_str())+atof(parsedCsv[i+1][j].c_str())+atof(parsedCsv[i][j+1].c_str())+atof(parsedCsv[i+1][j+1].c_str()))/4;
     spec_integral = spec_integral + spec_av*pow(0.025+0.05*i+0.025,2);
  }}


  spec_integral = spec_integral*0.05*1*(4*M_PI)/(pow(2*M_PI,3));

   std::cout <<spec_integral << "\n";

  
   spec_integral = spec_integral /(pow(0.197345,1) * 1000);

  std::cout << spec_integral;

return 0;

}

