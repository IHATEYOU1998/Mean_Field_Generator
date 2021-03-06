
#include <stdio.h>
#include <iostream>
#include <math.h>

#include "echelon.cpp"

int main()
{
// Values to interpolate
  long double k_test = 0.5;   // x value
  long double es_test = 1.5;  // y value

  
// Set up known values needed for matrix
  long double x2 = 0;
  long double x3 = 0;
  long double x6 = 1;
  long double x7 = 1;
  long double x10 = 2;
  long double x11 = 2;
  long double x14 = 3;
  long double x15 = 3;

  long double y5 = 0;
  long double y6 = 1;
  long double y7 = 2;
  long double y8 = 3;
  long double y9 = 0;
  long double y10 = 1;
  long double y11 = 2;
  long double y12 = 3;

// Set up values in using equation for Z-vales
  long double x5 = 1;
  long double x8 = 1;
  long double x9 = 2;
  long double x12 = 2;

  long double y2 = 1;
  long double y3 = 2;
  long double y14 = 1;
  long double y15 = 2;

// Initializing the Z-vales from function
  long double z2 = pow(y2,1)*x2;
  long double z3 = pow(y3,1)*x3;
  long double z5 = pow(y5,1)*x5;  
  long double z6 = pow(y6,1)*x6;   
  long double z7 = pow(y7,1)*x7;
  long double z8 = pow(y8,1)*x8; 
  long double z9 = pow(y9,1)*x9;
  long double z10 = pow(y10,1)*x10;
  long double z11 = pow(y11,1)*x11;
  long double z12 = pow(y12,1)*x12;
  long double z14 = pow(y14,1)*x14;
  long double z15 = pow(y15,1)*x15;
    
// Solve the Matrix for power series' variables using row echelon form
// all equations solved on paper. In the form of Ax = d (Last column is d, 0-39 is A)
long double M[40][41] = {{1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,z6},
                         {1,0,(y5-y6),0,0,pow(y5-y6,2),0,0,0,pow(y5-y6,3),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,z5},
                         {1,0,(y7-y6),0,0,pow(y7-y6,2),0,0,0,pow(y7-y6,3),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,z7},
                         {1,(x10-x6),0,pow(x10-x6,2),0,0,pow(x10-x6,3),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,z10},
                         {1,(x2-x6),0,pow(x2-x6,2),0,0,pow(x2-x6,3),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,z2},
  
                         {0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,z7},
                         {0,0,0,0,0,0,0,0,0,0,1,0,(y8-y7),0,0,pow(y8-y7,2),0,0,0,pow(y8-y7,3),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,z8},
                         {0,0,0,0,0,0,0,0,0,0,1,0,(y6-y7),0,0,pow(y6-y7,2),0,0,0,pow(y6-y7,3),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,z6},
                         {0,0,0,0,0,0,0,0,0,0,1,(x3-x7),0,pow(x3-x7,2),0,0,pow(x3-x7,3),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,z3},
                         {0,0,0,0,0,0,0,0,0,0,1,(x11-x7),0,pow(x11-x7,2),0,0,pow(x11-x7,3),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,z11},

                         {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,z10},
                         {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,(x14-x10),0,pow(x14-x10,2),0,0,pow(x14-x10,3),0,0,0,0,0,0,0,0,0,0,0,0,0,z14},
                         {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,(x6-x10),0,pow(x6-x10,2),0,0,pow(x6-x10,3),0,0,0,0,0,0,0,0,0,0,0,0,0,z6},
                         {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,(y11-y10),0,0,pow(y11-y10,2),0,0,0,pow(y11-y10,3),0,0,0,0,0,0,0,0,0,0,z11},
                         {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,(y9-y10),0,0,pow(y9-y10,2),0,0,0,pow(y9-y10,3),0,0,0,0,0,0,0,0,0,0,z9},
                         
                         {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,z11},
                         {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,(x15-x11),0,pow(x15-x11,2),0,0,pow(x15-x11,3),0,0,0,z15},
                         {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,(x7-x11),0,pow(x7-x11,2),0,0,pow(x7-x11,3),0,0,0,z7},
                         {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,(y12-y11),0,0,pow(y12-y11,2),0,0,0,pow(y12-y11,2),z12},
                         {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,(y10-y11),0,0,pow(y10-y11,2),0,0,0,pow(y10-y11,3),z10},
                        
                         {0,-1,0,0,0,0,0,0,0,0,0,1,0,0,(y6-y7),0,0,0,pow(y6-y7,2),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                         {0,1,0,0,(y7-y6),0,0,0,pow(y7-y6,2),0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                         {0,0,0,-1,0,0,0,0,0,0,0,0,0,1,0,0,0,(y6-y7),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                         {0,0,0,1,0,0,0,(y7-y6),0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                        
                         {0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,(x6-x10),0,0,pow(x6-x10,2),0,0,0,0,0,0,0,0,0,0,0,0,0},
                         {0,0,1,0,(x10-x6),0,0,pow(x10-x6,2),0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                         {0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,(x6-x10),0,0,0,0,0,0,0,0,0,0,0,0},
                         {0,0,0,0,0,1,0,0,(x10-x6),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                        
                         {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,1,0,0,(y10-y11),0,0,0,pow(y10-y11,2),0,0},
                         {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,(y11-y10),0,0,0,pow(y11-y10,2),0,0,-1,0,0,0,0,0,0,0,0,0},
                         {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,1,0,0,0,(y10-y11),0,0,0},
                         {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,(y11-y10),0,0,0,0,0,-1,0,0,0,0,0,0,0},
                        
                         {0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,(x7-x11),0,0,pow(x7-x11,2),0,0,0},
                         {0,0,0,0,0,0,0,0,0,0,0,0,1,0,(x11-x7),0,0,pow(x11-x7,2),0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0},
                         {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,(x7-x11),0,0},
                         {0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,3*(x6-x10),0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                        
                         {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,1,0,0,0,3*(y10-y11),0},
                         {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,3*(y11-y10),0,0,0,0,0,-1,0,0,0,0,0},
                         {0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,1,0,0,0,3*(y6-y7),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                         {0,0,0,1,0,0,3*(x10-x6),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0} };

   to_reduced_row_echelon_form(M);

// EDGE CASE: IF ALL 4 VALUES AROUND IT ARE 0 (should not interpolate: its zero)
// Interperlator will try and give it a super small number.
   //  if (z6 == 0 and z7 == 0  and z10 == 0 and z11 == 0){
   // return 0;}
   
// Find which point the function is closest to: 6, 7, 10, or 11
   double dist_6 = sqrt(pow(x6 - k_test,2) + pow(y6 - es_test,2));
   double dist_7 = sqrt(pow(x7 - k_test,2) + pow(y7 - es_test,2));
   double dist_10 = sqrt(pow(x10 - k_test,2) + pow(y10 - es_test,2));
   double dist_11 = sqrt(pow(x11 - k_test,2) + pow(y11 - es_test,2));



   long double spec_fun_6 = (M[0][40]+M[1][40]*(k_test-x6)+M[2][40]*(es_test-y6)+M[3][40]*pow(k_test-x6,2)+M[4][40]*(k_test-x6)*(es_test-y6)+M[5][40]*pow(es_test-y6,2)+M[6][40]*pow(k_test-x6,3)+M[7][40]*pow(k_test-x6,2)*(es_test-y6)+M[8][40]*(k_test-x6)*pow(es_test-y6,2)+M[9][40]*pow(es_test-y6,3));
   long double spec_fun_7 = (M[10][40]+M[11][40]*(k_test-x7)+M[12][40]*(es_test-y7)+M[13][40]*pow(k_test-x7,2)+M[14][40]*(k_test-x7)*(es_test-y7)+M[15][40]*pow(es_test-y7,2)+M[16][40]*pow(k_test-x7,3)+M[17][40]*pow(k_test-x7,2)*(es_test-y7)+M[18][40]*(k_test-x7)*pow(es_test-y7,2)+M[19][40]*pow(es_test-y7,3));
   long double spec_fun_10 = (M[20][40]+M[21][40]*(k_test-x10)+M[22][40]*(es_test-y10)+M[23][40]*pow(k_test-x10,2)+M[24][40]*(k_test-x10)*(es_test-y10)+M[25][40]*pow(es_test-y10,2)+M[26][40]*pow(k_test-x10,3)+M[27][40]*pow(k_test-x10,2)*(es_test-y10)+M[28][40]*(k_test-x10)*pow(es_test-y10,2)+M[29][40]*pow(es_test-y10,3));
  long double spec_fun_11 = (M[30][40]+M[31][40]*(k_test-x11)+M[32][40]*(es_test-y11)+M[33][40]*pow(k_test-x11,2)+M[34][40]*(k_test-x11)*(es_test-y11)+M[35][40]*pow(es_test-y11,2)+M[36][40]*pow(k_test-x11,3)+M[37][40]*pow(k_test-x11,2)*(es_test-y11)+M[38][40]*(k_test-x11)*pow(es_test-y11,2)+M[39][40]*pow(es_test-y11,3));


   std::cout << spec_fun_6 << "\n dist6: " << dist_6 << "\n";
   std::cout << spec_fun_7 << "\n dist7: " << dist_7 << "\n";
   std::cout << spec_fun_10 << "\n dist10: " << dist_10 << "\n";
   std::cout << spec_fun_11 << "\n dist11: " << dist_11 << "\n";

// if you want to print the matrix in reduced row echelon form 
//  for (int i = 0; i < 40; ++i)
//  {
//    for (int j = 0; j < 41; ++j)
//      std::cout << M[i][j] << '\t';
//      std::cout << "\n";
//      }

// if you want to print out the power series coefficient: zx6(0-9) , z7(10-19) , z10(21-29) , z11(30-39)  
  for (int i = 0; i <40; ++i)
    {
      std::cout << M[i][40];
      std::cout << "\n";
      }

     return 0;
}

