#include "ConsState.h"
#include "SolidSystem.h"

using namespace Eigen;
using namespace std;

//1. Calc U,
ConsState::ConsState(Vector3d u_i, Matrix3d F_i, double rho_i, double S_0, double e_0):rho_0(rho_i), S(S_0)
{
	u_0 = u_i;
	F_0 = F_i;
	VectorXd temp(e_size);
	v = temp;
	rho = rho_0/(F_0.determinant()); 
	/* cout << "rho Calculated " << rho << " mod u " << mod_u() << endl; */
	E_0 = e_0 + (u_0(0)*u_0(0)+u_0(1)*u_0(1)+u_0(2)*u_0(2))/2.;
	E = E_0;
	Vector3d mom = rho*u_0;
	Matrix3d rhoF = rho*F_0;
	v[0] = mom(0), v[1] = mom(1), v[2] = mom(2),
		v[3] = rhoF(0,0), v[4] = rhoF(0,1), v[5] = rhoF(0,2),
		v[6] = rhoF(1,0), v[7] = rhoF(1,1), v[8] = rhoF(1,2),
		v[9] = rhoF(2,0), v[10] = rhoF(2,1), v[11] = rhoF(2,2),
		v[12] = rho*E_0;
	/* cout << "calc Rho" << rhoNew << endl; */
}

double ConsState::mod_u()
{
	return u_(0)*u_(0)+u_(1)*u_(1)+u_(2)*u_(2);
}

void ConsState::initialStates(std::ostream& os)
{
	os <<   "Initial States\n " << 
		"Density Initial\t" << rho_0 <<"\n" 
		"Energy Initial\t"	<< E_0 << "\n"
		"Velocity \t" 		<< u_0 << "\n";
}
// calc /* //default constructor */
ConsState::ConsState(){
	VectorXd temp(e_size);
	v = temp;
}
//destructor
ConsState::~ConsState(){}


Matrix3d ConsState::strainTensor(Matrix3d F)
{
	Matrix3d G = (F.inverse().transpose())*(F.inverse());	
	return G;
}


double ConsState::e_()
{
	/* cout << "internal from Energy "<< E_() - mod_u()/2. << endl; */
	return E_() - mod_u()/2.;
}


const string ConsState::name[] = {"rhou_1", "rhou_2", "rhou_3", "rhoF_11", "rhoF_12", "rhoF_13", "rhoF_21", "rhoF_22", "rhoF_23", "rho_31", "rho_32", "rho_33", "rhoE"};

/* // derivative incorrect */ 
/*  void ConsState::stress(Matrix3d G) */
/* { */
/* 	Matrix3d GI = G.inverse(); */
/* 	double gProd1, gProd2; */
/* 	for(int i = 0; i < 3; i++) */	
/* 		for(int j = 0; j < 3; j++) */
/* 		{ */
/* 			for(int k = 0; k < 3; k++) */
/* 			{ */
/* 				gProd1 = G(i,k)*G(j,k); */
/* 				gProd2 = G(i,k)*GI(j,k); */
/* 			} */
/* 			sigma(i,j) = -2*rho*(G(i,j)*I[3] + I[4]*(I[0]*G(i,j)- */
/* 							gProd1) + I[5]*I[2]*gProd2); */
/* 		} */
/* 	stressi(G); */
/* 	cout << "first\n" << sigma << endl; */
/* } */

//Will eventaully remove this function
//Requires the first few invariants and entropy

void ConsState::prepareOutputs()
{
	Matrix3d rhoF;
	ConsState temp = *this;	
	rhoF <<  v[3], v[4], v[5],
			 v[6], v[7], v[8],
			 v[9], v[10], v[11];
	/* rho = pow(rhoF.determinant()/2., 2.); */
	rho = sqrt(rhoF.determinant()/rho_0);
	//deformation gradient
	Matrix3d F = rhoF/rho;
	// Calculate finger tensor
	Matrix3d G = strainTensor(F);

	//calculate invariants
	/* cout << B_0 << endl; */
	invariants(G);
	dInvariants_G();

	//Stress
	stress(G);
}

double ConsState::rho_()
{
	/* rho = sqrt(rhoF.determinant()/rho_0); */
	/* cout << rho << endl; */
	return rho;
}

// throw error if index is greater than certain value
double ConsState::u_(int index)
{
	return v[index]/rho_();
}

double ConsState::sigma_(int index_i, int index_j)
{
	return sigma(index_i, index_j);
}

double ConsState::E_()
{
	// This must be calculated after the flux
	return v[12]/rho_();
}
//explicit expression for entropy calculated for stress

//overload operators
ConsState ConsState::operator+(ConsState param)
{
	ConsState temp = *this;	
	for (int i = 0; i < e_size; i++)
	{
		temp[i] += param[i];
	}	
	return temp; 
}

ConsState ConsState::operator-(ConsState param)
{       
	ConsState temp = *this;					
	for (int i = 0; i < e_size; i++)
	{
		temp[i] -= param[i];     
	}       
	return temp;
}

ConsState ConsState::operator*(double scalar)
{
	ConsState temp = *this;
	for (int i = 0; i < e_size; i++)
	{
		temp[i] *= scalar;        
	}                                    
	return temp;
}

ConsState operator*(double k, ConsState v)
{
	return v*k;
}

ConsState ConsState::operator/(ConsState rhs)
{
	ConsState temp = *this;
	for(int i = 0; i < e_size; i++)
	{
		temp[i] /= rhs[i];
	}
	return temp;
}	

/*  ConsState ConsState::operator=(ConsState rhs) */
/* { */
/* 	/1* for(int i = 0; i < e_size; i++) *1/ */
/* 	/1* { *1/ */
/* 	/1* 	(*this)[i] = rhs[i]; *1/ */
/* 	/1* } *1/ */
/* 	// i may decide to change v back to a vector */ 
/* 	// */
/* 	(*this).v = rhs.v; */
/* } */

std::ostream& operator<<(std::ostream& os, ConsState& param)
{
	/* param.initialstates(os); */
	/* os << "current states \n" << */
	/* 	  "density\t"    << param.rho_() << "\n" */
	/* 	  "energy\t"     << param.e_() << "\n" */
	/* 	  "u or f vector\t " << "\n"; */
	os.precision(3);
	for (int i = 0; i < ConsState::e_size; i++)
	{
		os << setw(3) << " |" << setw(7) << param[i];
	}	
	os << setw(7) << " | ";
	return os;
}

double ConsState::getInvariant(int index)
{
	return I[index];
}

