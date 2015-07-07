#ifndef ELASTICPRIMSTATE_H
#define ELASTICPRIMSTATE_H

#include <math.h>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <tvmet/Matrix.h>
#include <fstream>
#include <vector>
#include <iostream>
 
class ElasticPrimState
{
	//velocity, deformation gradient, entropy
	public: 
		ElasticPrimState(Eigen::Vector3d u, Eigen::Matrix3d F, double S);

		ElasticPrimState();

		~ElasticPrimState();	

		Eigen::Vector3d u_();		

		Eigen::Matrix3d F();
		
		double S_();

		Eigen::VectorXd getStateVector(){return v;};
		//! Acoustic tensor derivatives
			
		Eigen::Tensor<double, 4> dsigma_dG();

		Eigen::Tensor<double, 4> dG_dF();	
	
		//! Overloaded operators
		
		double& operator[] (int x){return v[x];}

		ElasticPrimState& ElasticPrimState::operator+=(const ElasticPrimState& consState);

		ElasticPrimState& ElasticPrimState::operator-=(const ElasticPrimState& consState);

		ElasticPrimState& ElasticPrimState::operator*=(const double& scalar);

		ElasticPrimState& ElasticPrimState::operator/=(const double& scalar);
		
/* ElasticPrimState operator=(ElasticPrimState); */	

		friend std::ostream& operator<<(const std::ostream&, const ElasticPrimState&);

	private:
		/* Initial states */
		static const int e_size = 13;
		
	//! Solution vector
		Eigen::VectorXd v;
};

		ElasticPrimState operator+(const ElasticPrimState& lhs, const ElasticPrimState& rhs);

		ElasticPrimState operator-(const ElasticPrimState& lhs, const ElasticPrimState& rhs);

		ElasticPrimState operator*(const ElasticPrimState& lhs, const double& s);

		ElasticPrimState operator*(const double& s, ElasticPrimState& rhs);

		ElasticPrimState operator/(const ElasticPrimState& lhs, const double& s);

#endif
