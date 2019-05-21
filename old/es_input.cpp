#include <iostream>
#include <cstdlib>
#include <math.h>
#include <stdlib.h>
#include <string>
#include "es_input.h"

// Takes in an impossible es_test and whether the user wants a random one or not
long double es_input(long double es_test, std::string want_rand)
{
/** generate a random momentum and energy value (if rand = Y) or ask the user for an input**/
  if (want_rand == "Y" or want_rand == "y")
  {
    srand(time(NULL));
    es_test = fmod(rand(), 998.5) + 1.5; // Generate random value between 1.5 and 998.5
  }
  else
  {
    while (es_test > 998.5 or es_test < 1.5) // Keeps user in the loop until they output a value between 1.5 and 998.5
    {
      std::cout << "Please enter an Energy value in units of inverse fermi between [1.5,998.5]: \n";
      std::cin >> es_test;
    }
  }
  return es_test*100; // *100 because the return function only returns integers
}

// use typeid(k_test).name() to check type for string error
