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

	FL <<	 0.98, 0, 0, 0.02, 1, 0.1, 0, 0, 1;
	uL << 0, 0.5, 1;
	SL = 1000;

	FR << 1., 0, 0, 0, 1., 0.1, 0, 0, 1.;
	uR << 0, 0, 0;
	SR = 0.;

  ElasticPrimState primStateL(uL, FL, SL);  
  ElasticPrimState primStateR(uR, FR, SR);  
  ElasticEOS eos("copper"); // the romenski eos by default (to change)
  System sys(eos);
  ElasticState consStateL = sys.primitiveToConservative(primStateL);
  ElasticPrimState primStateCheck = sys.conservativeToPrimitive(consStateL);
  //std::cout << "Check that conversion works" << std::endl;
  assert(primStateL==primStateCheck);
  assert(primStateR!=primStateCheck);
  std::cout << "\nPrimState check passes" << std::endl;
  //use left and right states: obtain godunov state.
  //check that acoustic tensor is correct.
  //compute strain tensor and invariants:
	const Eigen::Matrix3d G = sys.strainTensor(primStateL.F_());
	const Eigen::Vector3d I = sys.getInvariants(primStateL.F_());
  std::cout << "Compute acoustic tensor test" << std::endl;
  
  //test the contraction of 3rd order tensor by 
  ElasticPrimState pW = 0.5*(primStateL + primStateR);
  std::cout << pW << std::endl;
  std::cout << sys.AcousticTensor(pW) << std::endl;
  //this is not initially correct: so compute derivatives
  //test (rhoR - rhoL)/(FR - FL)
  double rhoL = sys.Density(primStateL);
  Eigen::Matrix3d FD_approx;
  double h = 1e-8; // for this test h = 1e-08 is the best, why is this? 
  // can we be more intelligent about how to implement this? 
  // the variable will be a function of several variables?
  // 
  for(int i = 0; i < 3; i++)
    for(int j = 0; j < 3; j++)
    {
      //apply to the left 
      Eigen::Matrix3d FL_temp = FL;
      FL_temp(i,j) += h;
      ElasticPrimState primStateLTemp(uL, FL_temp, SL);
      double rhoph = sys.Density(primStateLTemp);
      FD_approx(i,j) = (rhoph-rhoL)/h;
    }
  
  std::cout.precision(12);
  std::cout << FD_approx << std::endl;
  std::cout << -rhoL*FL.transpose().inverse() << std::endl;
  
  /* std::cout << drho_dF << std::endl; */

  /* //density */
  /* double rhoW = sys.Density(pW); */ 
  /* Eigen::Matrix3d FW = pW.F_(); */
    /* = (rhoR - rhoL)/(FR - FL); */

   
  //check eigendecomposition and acoustic tensor are computed correctly (check against Kevin's implementation)
  //
  //


}
