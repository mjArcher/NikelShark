#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdio.h>

using namespace std;

/* extern "C" */ 
/* { */
/* } */

//
//write a macro to expand G matrix into G11, G12, G13 etc (this should be fairly easy to do!) and include this macro whenever you need these variables. Do the same for all other variables. 
//G, F, 

void aFunction()
{
  double B0(1.), G11(10.), G22(1.), G33(2), G23(3), G32(1), G12(2), G13(1), S(2), T0(3), K0(3), cv(3), alpha(4), beta(3), gamma(5);
  printf("%f", G12);
  double result = 0;;
  #include "./ots22/romenskii_energy.inc"
  cout << result << endl;
}

int main(int argc, char ** argv)
{
  aFunction();
}




