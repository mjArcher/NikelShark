#include "SolidSystem.h"

//equation of state and system tests
//
//

// System Derivatives to test 
//dsigma/dG
//drho/dF
//dStress/dFe

//Equation of state derivatives to test
//


int main(void) 
{
	//initialize example states
	//
	
  /* ElasticPrimState pri(SquareMatrix( 1.01, 0.01, 0, -0.1, 0.95, 0, 0, 0, 0.99), */
		       /* TinyVector<double, 3>(0, 0, 0), */
		       /* 10, 0.6); */
  /* SquareMatrix Fp(1.2,0.1,0.0,0.02,1.02,0,0,0,0.9); */
  /* Fp/=cbrt(Fp.Determinant()); */

  /* ElasticPrimState priB(SquareMatrix( 1, 0.1, 0, -0.2, 1, 0, 0, 0, 1), */
		       /* TinyVector<double, 3>(0, 0, 0), */
		       /* 0, 0); */


  Eigen::Matrix3d FL, FR;
  Eigen::Vector3d uL, uR;
	double SL, SR;

	FL <<	 1.01, 0.01, 0, -0.1, 0.95, 0, 0, 0, 0.99;
	uL << 0, 0.5, 1;
	SL = 1.e-3;

	FR << 1., 0, 0, 0, 1., 0.1, 0, 0, 1.;
	uR << 0, 0, 0;
	SR = 0.;

  ElasticPrimState primStateL(uL, FL, SL);  
  ElasticEOS eos("copper"); // the romenski eos by default (to change)
  System sys(eos);
  ElasticState consStateL = sys.primitiveToConservative(primStateL);
  ElasticPrimState primStateCheck = sys.conservativeToPrimitive(consStateL);
  std::cout << "Check that conversion works" << std::endl;
  if(primStateL==primStateCheck)
    std::cout << "PrimState check passes" << std::endl;

  //check eigendecomposition and acoustic tensor are computed correctly (check against Kevin's implementation)


}
