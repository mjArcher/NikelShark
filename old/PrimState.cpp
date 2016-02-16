#include "ElasticPrimState.h"
#include "SolidSystem.h"

using namespace Eigen;
using namespace std;

//1. Calc U,
ElasticPrimState::ElasticPrimState(Vector3d u, Matrix3d F, double S)
{
	v = VectorXd(ElasticPrimState::e_size);

	v[0] << u(0), u(1), u(2),
		F(0,0), F(0,1), F(0,2),
		F(1,0), F(1,1), F(1,2),
		F(2,0), F(2,1), F(2,2),	S; /* cout << "calc Rho" << rhoNew << endl; */
}

// default constructor

ElasticPrimState::ElasticPrimState(){
	v = VectorXd(ElasticPrimState::e_size);
}

ElasticPrimState::~ElasticPrimState(){}

//! velocity vector
Vector3d ElasticPrimState::u_()
{
	Vector3d vel;
	vel << v(0), v(1), v(2);
	return vel;
}

//! Rho * deformation gradient matrix
Vector3d ElasticPrimState::F_()
{
	Matrix3d rhoF;
	rhoF << v(3), v(4), v(5), v(6), v(7), v(8), v(9), v(10), v(11);
	return rhoF;
}

//! rhoE double
double ElasticPrimState::S_()
{
	return v[12];
}

std::ostream& operator<<(std::ostream& os, const ElasticPrimState& param){
	os.precision(3);
	for (int i = 0; i < ElasticPrimState::e_size; i++)
	{
		os << setw(3) << " |" << setw(7) << param[i];
	}	
	os << setw(7) << " | ";
	return os;
}

ElasticPrimState& ElasticPrimState::operator+=(const ElasticPrimState& consState)
{
	// actual addition of rhs to *this
	v += consState.getStateVector();
	return *this;
}

ElasticPrimState& ElasticPrimState::operator-=(const ElasticPrimState& consState)
{
	// actual addition of rhs to *this
	v -= consState.getStateVector();
	return *this;
}

ElasticPrimState& ElasticPrimState::operator*=(const double& scalar)
{
	// actual addition of rhs to *this
	v *= scalar;
	return *this;
}

ElasticPrimState& ElasticPrimState::operator/=(const double& scalar)
{
	// actual addition of rhs to *this
	v /= scalar;
	return *this;
}

// const safety
ElasticPrimState operator+(const ElasticPrimState& lhs, const ElasticPrimState& rhs)
{
	ElasticPrimState C(lhs);
  C += rhs;
  return C;
}

//const safety
ElasticPrimState operator-(const ElasticPrimState& lhs, const ElasticPrimState& rhs)
{
	ElasticPrimState C(lhs);
  C -= rhs;
  return C;
}

ElasticPrimState operator*(const ElasticPrimState& lhs, const double& s)
{
  return lhs/s;
}

ElasticPrimState operator*(const double& s, ElasticPrimState& rhs)
{
	rhs /= s;
  return rhs;
}

ElasticPrimState operator/(const ElasticPrimState& lhs, const double& s)
{
  return lhs*s;
}

//test functions - eigen library
Tensor<double, 4> dsigma_dG()
{
	Tensor<double, 4>	t(10, 10, 10, 10);	
	
	return t;
}

Tensor<double, 4> dG_dF()
{
	Tensor<double, 4>	t(10, 10, 10, 10);	
	return t;
}
