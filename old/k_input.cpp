#include <iostream>
#include <cstdlib>
#include <math.h>
#include <stdlib.h>
#include <string>
#include "k_input.h"

// Takes in an impossible es_test and whether the user wants a random one or not
int k_input(long double k_test, std::string want_rand)
{
/** generate a random momentum and energy value (if rand = Y) or ask the user for an input**/
  if (want_rand == "Y" or want_rand == "y")
  {
    srand(time(NULL));
    k_test = fmod(rand(), 9.925) + 0.075; // Random valye between 0.075 and 9.925
//    es_test = fmod(rand(), 999) + 0.5;
  }
  else
  {
    while (k_test > 9.925 or k_test < 0.075) // loops for user to input value between 0.075 and 9.925
    {
      std::cout << "Please enter a Momentum value in units of inverse fermi between [0.075,9.925]: \n";
      std::cin >> k_test;
    }
//    while (es_test > 999.5 or es_test < 0.5)
  //  {
    //  std::cout << "Please enter an Energy value in units of inverse fermi between [0.5,999.5]: \n";
      //std::cin >> es_test;
   // }
  }
  return k_test*10000; // *10000 because return only returns an integer
  //std::make_tuple(k_test, es_test);
}

// use typeid(k_test).name() to check type for string error
