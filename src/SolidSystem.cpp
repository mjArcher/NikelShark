#include "SolidSystem.h"

// getInvariants
using namespace Eigen;
using namespace std; 

System::System(const ElasticEOS eos):Eos(eos){}

System::~System(){}

//solid problem worth addressing for efficiency purposes is that when convert from 

ElasticPrimState System::conservativeToPrimitive(const ElasticState& consState) const 
{
	//calculate rho from inside class or calculate here
	const double rho = Density(consState);
	const Vector3d u = consState.mom()/rho;
	const Matrix3d F = consState.rhoF()/rho;
	double E = consState.rhoE()/rho;
	Vector3d inv = getInvariants(F);
	//! calculate internl energy from consState
	double ie =	E - u.dot(u)/2.; 	
	/* cout << "Cons INTERNAL ENERGY " << ie << endl; */
	double S = Eos.entropy(inv, ie);
	ElasticPrimState primState(u, F, S);
	return primState;
}

ElasticState System::primitiveToConservative(const ElasticPrimState& primState) const
{
	//finger tensor, invariants, internal energy, Energy (using eigen dot), 
	const Matrix3d F = primState.F_();
	const Vector3d invariants = getInvariants(F);
	/* cout << "Primitive invariants " << invariants << " G " << G << endl; */
	const double ie =	Eos.internalEnergy(invariants, primState.S_());
	Vector3d u = primState.u_();
	double E = ie + u.dot(u)/2.;
	/* cout << " Prim INTERNAL ENERGY " << ie << endl; */
	double rho = Density(primState);
	ElasticState consState(rho*primState.u_(), rho*F, rho*E);
	return consState;
}

double System::Density(const Eigen::Matrix3d& F) const 
{
	double rho = Eos.rho0/(F.determinant());	
	return rho;
}

double System::Density(const ElasticPrimState& primState) const 
{
	Matrix3d F = primState.F_();
	double rho = Eos.rho0/(F.determinant());	
	return rho;
}

double System::Density(const ElasticState& consState) const
{
	Matrix3d rhoF = consState.rhoF();	
	double rho = sqrt(rhoF.determinant()/Eos.rho0);	
	return rho;
}

Matrix3d System::strainTensor(const Matrix3d& F) const
{
	Matrix3d G = (F.inverse().transpose())*(F.inverse());	
	/* Matrix3d G = F.inverse()*(F.inverse().transpose()); */	
	return G;
}

Vector3d System::getInvariants(const Matrix3d& F) const
{
	const Matrix3d G = strainTensor(F);
	// Calculate invariants
	Vector3d I;
	I[0] = G.trace();
	I[1] = 0.5*(G.trace()*G.trace()-(G*G).trace());
	I[2] = G.determinant();
	return I;
}

ElasticState System::flux(const ElasticPrimState& primState) const 
{
  return flux(primitiveToConservative(primState), primState);
}

ElasticState System::flux(const ElasticState& consState) const
{
	return flux(consState, conservativeToPrimitive(consState));
}

ElasticState System::flux(const ElasticState& consState, const ElasticPrimState& primState) const
{
	ElasticState Fl; //flux
	Matrix3d rhoF = consState.rhoF();	
	double rhoE = consState.rhoE();
	Vector3d m = consState.mom();
	Vector3d u = primState.u_();
	Matrix3d sigma = stress(primState);
	//but we can use both the primitive state and conservative state to compute the flux
	// write for loops here instead

	Fl[0]  = m(0)*u(0) - sigma(0,0);
	Fl[1]  = m(1)*u(0) - sigma(0,1);
	Fl[2]  = m(2)*u(0) - sigma(0,2);
	Fl[3]  = 0; 
	Fl[4]  = 0; 
	Fl[5]  = 0; 
	Fl[6]  = rhoF(1,0)*u(0) - rhoF(0,0)*u(1);
	Fl[7]  = rhoF(1,1)*u(0) - rhoF(0,1)*u(1);
	Fl[8]  = rhoF(1,2)*u(0) - rhoF(0,2)*u(1);
	Fl[9]  = rhoF(2,0)*u(0) - rhoF(0,0)*u(2);
	Fl[10] = rhoF(2,1)*u(0) - rhoF(0,1)*u(2);
	Fl[11] = rhoF(2,2)*u(0) - rhoF(0,2)*u(2);
	Fl[12] = rhoE*u(0) - u(0)*sigma(0,0)-u(1)*sigma(0,1)-u(2)*sigma(0,2);
	return Fl;

	/* Fl[0]  = v[0]*u[0] - sigma(0,0); */
	/* Fl[1]  = v[1]*u[0] - sigma(0,1); */
	/* Fl[2]  = v[2]*u[0] - sigma(0,2); */
	/* Fl[3]  = 0; */ 
	/* Fl[4]  = 0; */ 
	/* Fl[5]  = 0; */ 
	/* Fl[6]  = rho*F(1,0)*u(0) - rho*F(0,0)*u(1); */
	/* Fl[7]  = rho*F(1,1)*u(0) - rho*F(0,1)*u(1); */
	/* Fl[8]  = rho*F(1,2)*u(0) - rho*F(0,2)*u(1); */
	/* Fl[9]  = rho*F(2,0)*u(0) - rho*F(0,0)*u(2); */
	/* Fl[10] = rho*F(2,1)*u(0) - rho*F(0,1)*u(2); */
	/* Fl[11] = rho*F(2,2)*u(0) - rho*F(0,2)*u(2); */
	/* Fl[12] = rho*u(0)*E_() - u(0)*sigma(0,0)-u(1)*sigma(0,1)-u(1)*sigma(0,2); */
}


//kronecker in utility functions
Matrix3d System::stress(const ElasticPrimState& primState) const
{
	Matrix3d F = primState.F_();
	Vector3d I = getInvariants(F);
	double rho = Density(primState);
	// can we just get the entropy from e primitive state?
	// yes but we still need to calculate the stress
	
	/* double S = Eos.entropy(inv, EosVgc */
	double S = primState.S_();
	Vector3d Id = Eos.depsi_dI(I, S);
	Matrix3d G = strainTensor(F);
	Matrix3d GI = G.inverse();
	double sigma_t;
	double term1 = 0., term2 = 0., term3 = 0.;
	Matrix3d sigma;
	//this needs to be reimplemented to take advantage of library matrix multiplication
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
		{	
			sigma_t = 0;
			for(int k = 0; k < 3; k++)
			{
				term1 = Id[0]*kron(j,k);
				term2 =	Id[1]*(I[0]*kron(j,k) - G(j,k));
				term3 = Id[2]*(I[2]*GI(j,k));
				sigma_t += -2.*rho*G(i,k)*(term1 + term2 + term3);
			}	
			sigma(i,j) = sigma_t;	
		}
	return sigma;
}

Matrix3d System::AcousticTensor(const ElasticPrimState& primState) const
{
  double rho = Density(primState);
	const Matrix3d F = primState.F_();
	const Matrix3d G = strainTensor(F);
	const Vector3d I = getInvariants(F);
  Matrix3d omega;
	const SquareTensor3 dsdFe = dstress_dF(primState, G, I); //need to decide on a return type here
	//construct acoustic tensor	i, j, k, m : i does not change for 1D 
	//population of acoustic tensor
	//dirn = 0 for 1D problem
  for(int i=0; i<3; ++i){
    for(int j=0; j<3; ++j){
      omega(i,j) = 0.0;
      for(int k=0; k<3; ++k){	
				omega(i,j) += dsdFe(i,j,k) * F(0,k); //symmetric 
      }
      omega(i,j) /= rho;
    }
  }
  return omega;
}

// reconsider the implementation of this
vector<Matrix3d> System::dep_dF(const SquareTensor3 dI_dF, const Matrix3d depsi_dI_dI) const
{
	vector<Matrix3d> depsi_dF(3);	

	Matrix3d dEps_pdF;

	for(int p = 0; p < 3; ++p)//row in depsi_dI_dI
	{
		dEps_pdF =  depsi_dI_dI(p, 0)*dI_dF[0] + depsi_dI_dI(p, 1)*dI_dF[1] + depsi_dI_dI(p, 2)*dI_dF[2]; //
		depsi_dF[p] = dEps_pdF;
	}

	/* depsi_dF[0] = depsi_dI_dI[0] * dI_dF[0]; //row */ 
	/* depsi_dF[1] = depsi_dI_dI[1] * dI_dF[1]; */
	/* depsi_dF[2] = depsi_dI_dI[2] * dI_dF[2]; */
	return depsi_dF;
}

// Need only store 
// Aijk 
// components of the acoustic tensor
/* SquareTensor3 System::dstress_dF(const ElasticPrimState& primState, const Matrix3d& G, const Vector3d& I) const */
/* { */
/* 	const double rho = Density(primState); */
/* 	const Matrix3d F = primState.F_(); */	
/* 	const double ie =	Eos.internalEnergy(I, primState.S_()); */
/* 	const Matrix3d dstress_drho = stress(primState)/rho; */
/* 	const Matrix3d drho_dF = - rho * (F.inverse()).transpose(); */
/* 	const Matrix3d m2rho = -2. * rho * G; */
/* 	const vector<Eigen::Matrix3d> depsdF = dep_dF(primState.dI_dF(G,I), Eos.depsi_dI_dI(I, primState.S_())); */
/* 	const Vector3d de_dI = Eos.depsi_dI(I, primState.S_()); */
/* 	const SquareTensor3 dsdeps = m2rho * primState.dI_dG(G,I); //cannot remember this derivation? */
/* 	double sigma_rho; */

/* 	vector<Matrix3d> A(3); */
/* 	const Matrix3d sigma = stress(primState); */

/* 	for (int k = 0; k < 3; k++) */
/*   { */
/* 		sigma_rho = dstress_drho(0,k); */
/* 		Matrix3d ds_dG = primState.dsigma_dG(k, de_dI, G, rho); */ 
/*     for(int j = 0; j < 3; j++) */
/*     { */
/*       for(int m = 0; m < 3; m++){ */
/* 				/1* A(k, m, j) -- row, column, depth -- index using A[j](k,m) - this is a std array *1/ */
/* 				Matrix3d dGdF = primState.dG_dF(G,F,j,m); //returns 2d slice of tensor at j, m (denominator constant) */
/* 				A[k](j,m) = sigma_rho * drho_dF(j,m) */ 
/* 					 + depsdF[0](j,m) * dsdeps[0](0,k) */ 
/*            + depsdF[1](j,m) * dsdeps[1](0,k) */ 
/*            + depsdF[2](j,m) * dsdeps[2](0,k) */
/*         + (dGdF * ds_dG).trace(); */
/* 			} */
/* 		} */
/* 	} */ 
/* 	SquareTensor3 dstressdF(A); */
/* 	return dstressdF; */
/* } */

SquareTensor3 System::dstress_dF(const ElasticPrimState& primState, const Matrix3d& G, const Vector3d& I) const
{
	const double rho = Density(primState);
	const Matrix3d F = primState.F_();	
	const double ie =	Eos.internalEnergy(I, primState.S_());
	const Matrix3d dstress_drho = stress(primState)/rho;
	const Matrix3d drho_dF = - rho * (F.inverse()).transpose();
	const Matrix3d m2rho = -2. * rho * G;
	const vector<Eigen::Matrix3d> depsdF = dep_dF(primState.dI_dF(G,I), Eos.depsi_dI_dI(I, primState.S_()));
	const Vector3d de_dI = Eos.depsi_dI(I, primState.S_());
	const SquareTensor3 dsdeps = m2rho * primState.dI_dG(G,I); //cannot remember this derivation?
	double sigma_rho;

	vector<Matrix3d> A(3);
	const Matrix3d sigma = stress(primState);

	for (int k = 0; k < 3; k++)
  {
		sigma_rho = dstress_drho(0,k);
		Matrix3d ds_dG = primState.dsigma_dG(k, de_dI, G, rho); 
    for(int j = 0; j < 3; j++)
    {
      for(int m = 0; m < 3; m++){
				/* A(k, m, j) -- row, column, depth -- index using A[j](k,m) - this is a std array */
				Matrix3d dGdF = primState.dG_dF(G,F,j,m); //returns 2d slice of tensor at j, m (denominator constant)
				A[j](k,m) = sigma_rho * drho_dF(j,m) 
					 + depsdF[0](j,m) * dsdeps[0](0,k) 
           + depsdF[1](j,m) * dsdeps[1](0,k) 
           + depsdF[2](j,m) * dsdeps[2](0,k)
        + (dGdF * ds_dG).trace();
			}
		}
	} 

	SquareTensor3 dstressdF(A);
	return dstressdF;

}

//construct and then decompose acoustic tensor: get Q, D, A and then construct Le, Re and 
VectorXd System::stateEigenDecompose(const ElasticPrimState& pW,
		const int dirn,
		vector<ElasticPrimState>& Le,
		vector<ElasticPrimState>& Re) const 
{
	const double rho = Density(pW);
	const Matrix3d G = strainTensor(pW.F_());
	/* const Invariants inv = G.getInvariants(); */
	const ElasticState consState = primitiveToConservative(pW);
	
	const Matrix3d omega = AcousticTensor(pW);

}

ElasticState System::godunovFlux(const ElasticState& qL, const ElasticState& qR) const
{
  return flux(godunovState(conservativeToPrimitive(qL), conservativeToPrimitive(qR)));
}

/* ElasticState System::godunovFlux(const ElasticPrimState& pL, const ElasticPrimState& pR) */
/* { */
/*   return flux(primitiveToConservative(godunovState(pL, pR))); */
/* } */

//compute state at interface between L and R states 
ElasticPrimState System::godunovState(ElasticPrimState pL, ElasticPrimState pR) const
{
  if (pL == pR)
    return pL;
  // compute the eigenvalues, L and R eigenvectors.
  // decompose the acoustic tensor to get Q, D and Q
  // construct the diagonal matrix of eigenvalues (w).
  // construct the godunov state based on this 
  //
  /* const ElasticPrimState pW = interpolate(0.5, Q1, Q2);  // this is a very complicated function (can we just take the avg?) */

  const ElasticPrimState pW = 0.5*(pL + pR);
  vector<ElasticPrimState> Le(14), Re(14); //test if this has actually initialised correctly.
  const VectorXd lambda = stateEigenDecompose(pW, 0, Le, Re); // within this construct the acoustic tensor

  /* VectorXd VectorXd(ElasticPrimState::e_size); */
  //compute the primitive variable jump across the Riemann problem:
  const ElasticPrimState dp = pR - pL; 

  ElasticPrimState god = pL;
  for(int w=0; w<14; ++w) {
    if(lambda[w]<0.0) {
      /* const double strength = dot(Le[w], dp); //need to rewrite this line */ 
      /* if(strength!=0.0) god+=strength*Re[w]; */
    }
  }
  return god;
}
	
/* TinyVector<double, 14> System::StateEigenDecompose(const ElasticPrimState& pW, */
/* 		const int dirn, */
/* 		TinyVector<ElasticPrimState, 14>& Le, */
/* 		TinyVector<ElasticPrimState, 14>& Re) const */ 
/* { */
//in equation of state
//deps_dI
//deps_dI_dI

//Primstate
/* dI_dF */
//calculate in dstress_dF function
/* depsi_dF */

// Aij = dsigma_dF


