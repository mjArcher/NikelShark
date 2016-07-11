#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdio.h>

#include <Eigen/Dense>
#include <Eigen/Core>
#include <string>

#include "ElasticPrimState.h"
#include "SolidSystem.h"
#include "ElasticEOS/ElasticEOS.h"

using namespace std;
using namespace Eigen;

/* extern "C" */ 
/* { */
/* } */

#define strain_comps(G) double G11 = G(0,0), G12 = G(0,1), G13 = G(0,2), G22 = G(1,1), G23 = G(1,2), G33 = G(2,2)
/* #define eos_comps(eos) double B0 = eos.B0 */

double energy(Matrix3d G, ElasticEOS eos)
{
  double result;
  strain_comps(G);
  cout << G11 << " " << G(0,0) << endl;
  double B0, beta, K0, alpha, T0, cv, S, gamma;
  #include "./codegen/romenskii_energy.inc"
  /* eos_comps(eos) */
  cout << B0 << endl;

  return result;
}

//compute A using codegen functions
//


int main(int argc, char ** argv)
{
  printf("%s\n", "Testing SymPy energy");
  Matrix3d FL, FR;
	Vector3d uL, uR;
	double SL, SR;
	FL <<	 1.01, 0.01, 0, -0.1, 0.95, 0, 0, 0, 0.99;
	uL << 0, 0, 0;
	SL = 10.;

	FR << 1., 0, 0, 0, 1., 0.1, 0, 0, 1.;
	uR << 0, 0, 0;
	SR = 0.;

  ElasticEOS Eos("copper");
  ElasticPrimState primStateL(uL, FL, SL);
  System sys(Eos);
	const Matrix3d G = sys.strainTensor(primStateL.F_());
  double en = energy(G, Eos);
  /* printf("%f\n", en); */
}




