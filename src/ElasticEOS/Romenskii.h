#ifndef ROMENSKII_H
#define ROMENSKII_H

#include "ElasticEOS.h"
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <Eigen/Dense>
#include <Eigen/Core>

class Romenskii : public ElasticEOS 
{

	public:

		Romenskii(const std::string&);
		
		//! Default constuctor. Careful!
		virtual double entropy(const Eigen::Vector3d invariants, double internalEnergy) const;	

		virtual double internalEnergy(const Eigen::Vector3d invariants, double Entropy) const; //invariants, 

		virtual void checkEosConstants() const;

		virtual double soundSpeed() const;

		virtual Eigen::Vector3d depsi_dI(const Eigen::Vector3d I, double S) const;

		virtual Eigen::Matrix3d depsi_dI_dI(const Eigen::Vector3d& I, double S) const;

    virtual Eigen::Vector3d depsi_dI_dS(const Eigen::Vector3d& I, double S) const;

	private:
		//! Equation of state parameters 	
		//! Density
		double rho0; //this is public because it is required by  the density function in the system class 
		//! Heat capacity at constant volume
		double cv;

		//! Temperature
		double T0;

		//! squared bulk shear wave
		double B0;	

		//! squared bulk speed of sound
		double K0;

		//! Constants characterising the non-linear dependence of sound speeds and
		//! temperature on the mass density
		double alpha;

		double beta;

		double gamma;
};

#endif
