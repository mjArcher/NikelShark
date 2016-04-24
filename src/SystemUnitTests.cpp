#include "SolidSystem.h"
#include <vector> 

//equation of state and system tests

void check_dI_dF(ElasticPrimState pri, System sys);
void check_depsi_dI_dI(ElasticPrimState pri, System sys, ElasticEOS eos);
void check_depsi_dI_dS(ElasticPrimState pri, System sys, ElasticEOS eos);
void check_deps_dF(ElasticPrimState pri, System sys, ElasticEOS eos);
void check_drho_dF(ElasticPrimState, System sys);
void check_dsigma_dG(ElasticPrimState pL, System sys, ElasticEOS eos);
void check_dsigma_deps(ElasticPrimState prim, System sys);
void check_B(ElasticPrimState prim, System sys, ElasticEOS eos);
void construct_Eigenvectors(System sys, ElasticPrimState pW, ElasticEOS eos);


int main(void) 
{
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
  std::cout << "\nPRIMSTATE CHECK PASSED" << std::endl;
  //use left and right states: obtain godunov state.
  //check that acoustic tensor is correct.
  //compute strain tensor and invariants:
  std::cout.precision(10); 
  ElasticPrimState pW = 0.5*(primStateL + primStateR);
	const Eigen::Matrix3d G = sys.strainTensor(pW.F_());
	const Eigen::Vector3d I = sys.getInvariants(pW.F_());

  SquareTensor3 dstressdF = sys.dstress_dF(pW, G, I);
  /* std::cout << dstressdF << std::endl; */
  /* check_B(pW, sys, eos); */
  // therefore Q-1 = V 
  //construct eigenvectors
  construct_Eigenvectors(sys, pW, eos);

  //we can also do tests of the eigendecomposition

  //also experiment with single iteration, is this more efficient?
}

void construct_Eigenvectors(System sys, ElasticPrimState pW, ElasticEOS eos)
{
  std::cout << "CONSTRUCT EIGENVECTORS\n" << std::endl;
  int dirn = 0;
  /* Eigen::MatrixXcd Q = es.eigenvectors(); */
  Eigen::Matrix3d omega = sys.AcousticTensor(pW);
  Eigen::EigenSolver<Eigen::Matrix3d> es(omega);
  Eigen::Matrix3d V = es.pseudoEigenvectors();
  Eigen::Matrix3d eigs = es.pseudoEigenvalueMatrix();

  //we can check once they have been constructed
  //Q-1 = V
  //L = D^2
  Eigen::Matrix3d D;
  /* std::cout << L << std::endl; */
  D(0,0) = pow(eigs(0,0),0.5);
  D(1,1) = pow(eigs(1,1),0.5);
  D(2,2) = pow(eigs(2,2),0.5);
  Eigen::Matrix3d Q = V.inverse();
  Eigen::Matrix3d Qinv = V;
  /* std::cout << D << std::endl; */
 
  Eigen::Matrix3d perm = Eigen::Matrix3d::Zero();
  //permutation matrix
  perm(2,0) = 1;
  perm(1,1) = 1;
  perm(0,2) = 1;

  const Eigen::Matrix3d F = pW.F_();
  const Eigen::Matrix3d Finv = F.inverse();
	const Eigen::Matrix3d G = sys.strainTensor(F);
	const Eigen::Vector3d I = sys.getInvariants(F);
  double rho = sys.Density(pW);
  const Eigen::Matrix3d m2rho = -2.*rho*G;
	const SquareTensor3 dsdeps = m2rho * pW.dI_dG(G,I); //cannot remember this derivation?
  Eigen::Matrix3d dsdeps_ = dsdeps[2];
  Eigen::Matrix3d dsdS = eos.depsi_dI_dS(I, pW.S_())[2] * dsdeps_; //only non-zero component
  Eigen::Vector3d B;
  for(int i = 0; i < 3; i++) {
    B[i] = dsdS(dirn, i); 
  }
  B *= 1./rho;
  
  SquareTensor3 dstressdF = sys.dstress_dF(pW, G, I);
  Eigen::Matrix3d DQ = D*Q;
  Eigen::Matrix3d QA_11 = Q*dstressdF[0];
  Eigen::Matrix3d QA_12 = Q*dstressdF[1];
  Eigen::Matrix3d QA_13 = Q*dstressdF[2];
  Eigen::Vector3d QB1 = Q*B;
  //the first row is B when dir = 0, self explanatory
 
  Eigen::Matrix3d mpiQA_11 = -perm*QA_11;
  Eigen::Matrix3d mpiQA_12 = -perm*QA_12;
  Eigen::Matrix3d mpiQA_13 = -perm*QA_13;
  Eigen::Vector3d mpiQB1 = -perm*QB1;
  //representation of unit dyads 
  //
  //
  //L = 13 * 13 

  //L = 13 * 13
  //This is quite cool: specify compile time sizes
  //do these need to be know at compile time?
  //we can make our own tinyvector this way
  typedef Eigen::Matrix<double, 13, 1> Vector13d; //each row corresponds to eigenvector
  Vector13d testtype = Vector13d::Zero();
  //may need to change this representation: think of Leveque's linearised Riemann solver
 
  typedef Eigen::Matrix<Vector13d, 13, 1> Vector13V;
  Vector13V L, R;
  /* L[0] = Vector13d::Ones(); */
  /* L[1] = Vector13d::Ones()*4; */
  /* std::cout << L[0].dot(L[1]) << std::endl; */
  /* //we can also get row vector and column vector using row and col functions */
  /* //we can also concatenate vectors using this approach */
  Eigen::Vector3d zeros = Eigen::Vector3d::Zero();
  Eigen::Vector3d testV1 = Eigen::Vector3d::Ones()*2;
  Eigen::Vector3d testV2 = Eigen::Vector3d::Ones()*3;
  Eigen::Matrix3d onesM = Eigen::Matrix3d::Ones();

  //could either use eigen or ElasticPrimState.
  //Eigen does give some flexibility

  for(int i = 0; i < 3; i++)
  {
    L(i) << onesM<Vector3d>.row(i), onesM.row(i), onesM.row(i), onesM.row(i), 0;
    /* L(i) << DQ.row(i), QA_11.row(i), QA_12.row(i), QA_13.row(i), 0; */
    /* L[12-i] << DQ.row(i), -QA_11.row(i), -QA_12.row(i), -QA_13.row(i), QB1(i); */
  }


  /* L.transpose(); */
  /* std::cout << L[0] << std::endl; */

  // the r representation might not be correct here?
  
  /* std::cout << (Eigen::Matrix3d)perm << std::endl; */
  /* Eigen::Matrix3d permMatrix = (Eigen::Matrix3d)perm*Qinv; */
  /* std::cout << permMatrix << std::endl; */
  /* Matrix3d row1_3 = */ 
}
//problem with this approach is that I cannot assign by concatenating to a row or column of a matrix
//but it would be good to do it this way cause I could take the transpose at will.
  /* typedef Eigen::Matrix<double, 13, 13> Matrix13d; */ 
  /* Matrix13d L, R; */
  /* L = Matrix13d::Zero(); */
  /* std::cout << L.row(0) << std::endl; */
  /* typedef Eigen::Matrix<double, 13, 1> Vector13d; */
  /* Vector13d testtype = Vector13d::Zero(); */
  /* Vector13d::Ones(); */
  /* Eigen::Vector3d vec1; */
  /* vec1 = Eigen::Vector3d::Ones(); */

void compute_AcousticTensor(System sys, ElasticPrimState pW)
{
  Eigen::Matrix3d omega = sys.AcousticTensor(pW);
  std::cout << "Compute acoustic tensor test\n" << omega  << std::endl;
  printf("\nEigenvalues and orthogonal (eigenvectors) matrix Q\n"); 
  Eigen::EigenSolver<Eigen::Matrix3d> es(omega);
  std::cout << es.eigenvectors() << std::endl; 
  Eigen::MatrixXcd D = es.eigenvalues().asDiagonal();
  std::cout << D << std::endl;
  /* std::cout << "\n eigenvalues " << es.eigenvalues().asDiagonal() << std::endl; */ 
  /* std::cout <<  es.eigenvectors().col(0)[0] << std::endl; */ 
  printf("\nCheck that orthonormal\n");
  Eigen::MatrixXcd Q = es.eigenvectors();
  //how to only get non-trivial eigenvalues and eigenvectors from eigensolver
  Eigen::Matrix3d V = es.pseudoEigenvectors();
  Eigen::Matrix3d L = es.pseudoEigenvalueMatrix();

  printf("\nCheck that we can reconstruct the acoustic tensor: omega = Q*L*Q^-1");
  //this checks out
  Eigen::Matrix3d omegaR = V * L * V.inverse();
  std::cout << "\n" << omegaR  << std::endl;
}

void check_B(ElasticPrimState prim, System sys, ElasticEOS eos)
{
  std::cout << "Compute numerical result B " << std::endl;
  Eigen::Matrix3d sigma = sys.stress(prim);
  double h = 1e-7;
  ElasticPrimState prim_T(prim.u_(), prim.F_(), prim.S_() + h);
  Eigen::Matrix3d sigma_T = sys.stress(prim_T);

  std::cout << (sigma_T - sigma)/h << std::endl;
  printf("\n compute analytic result\n");

	const Eigen::Matrix3d G = sys.strainTensor(prim.F_());
	const Eigen::Vector3d I = sys.getInvariants(prim.F_());
  double rho = sys.Density(prim);
	const Eigen::Matrix3d m2rho = -2.*rho*G;
	const SquareTensor3 dsdeps = m2rho * prim.dI_dG(G,I); //cannot remember this derivation?
  /* std::cout << eos.depsi_dI_dS(I, prim.S_()) << std::endl; */
  check_depsi_dI_dS(prim, sys, eos);


  Eigen::Matrix3d dsdeps_ = dsdeps[2];
  Eigen::Matrix3d dsdS = eos.depsi_dI_dS(I, prim.S_())[2] * dsdeps_;
  //the first row is B when dir = 1, self explanatory
  std::cout << dsdS << std::endl;
  int dirn = 0;
  Eigen::Vector3d B;
  std::cout << " B " << std::endl;
  for(int i = 0; i < 3; i++)
  {
    B[i] = dsdS(dirn, i); 
  }
  B *= 1./rho;
  std::cout << B << std::endl;
}

//general matrix by matrix derivative
//pass in function pointers 
//try lambda functions C++11
void fd_dscalar_dmatrix(double scalar, Eigen::Matrix3d matrix)
{
  
}

//likely to be incorrect:
void check_dsigma_deps(ElasticPrimState prim, System sys)
{
  double rho = sys.Density(prim);
  Eigen::Matrix3d F = prim.F_();
	const Eigen::Vector3d I = sys.getInvariants(F);
	const Eigen::Matrix3d G = sys.strainTensor(prim.F_());
	const Eigen::Matrix3d m2rho = -2.*rho*G;
  const SquareTensor3 dsdeps = m2rho * prim.dI_dG(G,I);
  std::cout << "dsigma_deps analytical result: \n " <<  dsdeps << std::endl;
  for(int i = 0; i < 3; i++)
    for(int j = 0; j < 3; j++)
    {
       
    }
}

void check_drho_dF(ElasticPrimState pL, System sys)
{
  double rhoL = sys.Density(pL);
  Eigen::Matrix3d FL = pL.F_();
  std::cout << "drho_dF analytical result: \n " << -rhoL*FL.transpose().inverse() << std::endl;
  Eigen::Vector3d uL = pL.u_();
  double SL = pL.S_();
  Eigen::Matrix3d FD_approx;
  double h = 1e-8; // for this test h = 1e-08 is the best, why is this? 
  for(int i = 0; i < 3; i++)
    for(int j = 0; j < 3; j++)
    {
      Eigen::Matrix3d FL_temp = FL;
      FL_temp(i,j) += h;
      ElasticPrimState primStateLTemp(uL, FL_temp, SL);
      double rhoph = sys.Density(primStateLTemp);
      FD_approx(i,j) = (rhoph-rhoL)/h;
    }
  std::cout << "drho_dF numerical approximation: \n" << FD_approx << std::endl;
}

//checked against kevin
void check_dsigma_dG(ElasticPrimState pL, System sys, ElasticEOS eos)
{
  double rhoL = sys.Density(pL);
  Eigen::Matrix3d FL = pL.F_();
	const Eigen::Vector3d IL = sys.getInvariants(FL);
	const Eigen::Matrix3d GL = sys.strainTensor(pL.F_());
	const Eigen::Vector3d de_dI = eos.depsi_dI(IL, pL.S_());
  std::cout << "Analytical result: \n " << std::endl;
  
  for(int i = 0; i < 3; i++)
  {
    const Eigen::Matrix3d ds_dG = pL.dsigma_dG(i, de_dI, GL, rhoL);
    std::cout << ds_dG << "\n " << std::endl;
  }
}

void check_dG_dF(ElasticPrimState pL, System sys)
{
  std::cout << "Analytical: dG_dF " << std::endl;
  Eigen::Matrix3d FL = pL.F_();
	const Eigen::Matrix3d GL = sys.strainTensor(pL.F_());
  for (int j = 0; j < 3; j++)
    for (int m = 0; m < 3; m++){
      Eigen::Matrix3d dGdF = pL.dG_dF(GL,FL,j,m);
      int k1 = j*3 + m;
      std::cout << "\n" << k1 << ".\n" << dGdF << std::endl;
    }

}

//from Barton et al
//checked against kevin
void check_dI_dF(ElasticPrimState pri, System sys)
{
  Eigen::Matrix3d F = pri.F_();
  const Eigen::Matrix3d G = sys.strainTensor(F);
  const Eigen::Vector3d I = sys.getInvariants(F);
  std::cout << "Analytical: dI_dF\n" << pri.dI_dF(G,I) << std::endl;
}

void check_deps_dF(ElasticPrimState pri, System sys, ElasticEOS eos)
{
  Eigen::Matrix3d F = pri.F_();
  const Eigen::Matrix3d GO = sys.strainTensor(F);
  const Eigen::Vector3d IO = sys.getInvariants(F);
	const std::vector<Eigen::Matrix3d> depsdF = sys.dep_dF(pri.dI_dF(GO,IO), eos.depsi_dI_dI(IO, pri.S_())); // this is the analytical solution
  std::cout << "Analytical: deps_dF\n" << depsdF << std::endl;

  std::cout << "Finite-difference approximation (correct): deps_dF\n" << std::endl;
  double h = 1e-7;
  ElasticPrimState priT = pri;
  for(int s = 0; s < 3; s++)
  {
    for(int i = 0; i < 3; i++){
      for(int j = 0; j < 3; j++) {
        Eigen::Matrix3d FTemp = pri.F_();
        FTemp(i,j) += h;
        const Eigen::Matrix3d G = sys.strainTensor(FTemp);
        const Eigen::Vector3d I = sys.getInvariants(FTemp);
        //wrap in primstate to calculate entropy
        Eigen::Vector3d Id = eos.depsi_dI(I, pri.S_());
        std::cout << (Id[s] - eos.depsi_dI(sys.getInvariants(F), pri.S_())[s])/h << "\t";
      }   
      std::cout << "\n";
    } 
    std::cout << "\n";
  }
}

//this checks out
void check_depsi_dI_dI(ElasticPrimState pri, System sys, ElasticEOS eos)
{
  Eigen::Matrix3d F = pri.F_();
  const Eigen::Vector3d I = sys.getInvariants(F);
  Eigen::Matrix3d depsi_dI_dI = eos.depsi_dI_dI(I, pri.S_());
  printf("\ndepsi_dI_dI\n");
  std::cout << depsi_dI_dI << std::endl;
  //numerical bit
  double h = 1e-08;
  double S = pri.S_();
  Eigen::Vector3d eps_dI = eos.depsi_dI(I, S);
  //change each individual invariant:
  printf("\nFinite difference approximation: Correct \n");
  for(int i = 0; i < 3; i++)
  {
    Eigen::Vector3d IC_temp = I;
    IC_temp[i] += h;
    std::cout << (eos.depsi_dI(IC_temp, S) - eos.depsi_dI(I, S))/h << '\n' <<  std::endl;
  }

}

void check_depsi_dI_dS(ElasticPrimState pri, System sys, ElasticEOS eos)
{
  Eigen::Matrix3d F = pri.F_();
  const Eigen::Vector3d I = sys.getInvariants(F);
  Eigen::Vector3d depsi_dI_dS = eos.depsi_dI_dS(I, pri.S_());
  printf("\ndepsi_dI_dS\n");
  std::cout << depsi_dI_dS << std::endl;
  //numerical bit
  double h = 1e-08;
  double S = pri.S_();
  Eigen::Vector3d eps_dI = eos.depsi_dI(I, S);
  //change each individual invariant:
  printf("\nFinite difference approximation: Correct \n");
  std::cout << (eos.depsi_dI(I, S+=h) - eos.depsi_dI(I, S))/h << '\n' <<  std::endl;
}

