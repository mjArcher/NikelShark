#ifndef ELASTICEOS_H
#define ELASTICEOS_H

#include <cstdlib>
#include <iostream>
#include <sstream>


#include <Eigen/Dense>
#include <Eigen/Core>
// this is the Romenski equation of state
//

class ElasticEOS 
{
	public:

    //default 
		ElasticEOS() {}

		//! Default constuctor. Careful!
		virtual ~ElasticEOS() {}

		virtual double entropy(const Eigen::Vector3d invariants, double internalEnergy) const=0;	

		virtual double internalEnergy(const Eigen::Vector3d invariants, double Entropy) const=0; //invariants, 

		virtual void checkEosConstants() const=0;

		virtual double soundSpeed() const=0;

		virtual Eigen::Vector3d depsi_dI(const Eigen::Vector3d I, double S) const=0;

		virtual Eigen::Matrix3d depsi_dI_dI(const Eigen::Vector3d& I, double S) const=0;

    virtual Eigen::Vector3d depsi_dI_dS(const Eigen::Vector3d& I, double S) const=0;

    //! Reference density
    virtual double rho0_() const
    {
      return rho0;
    }

    //additional derivatives for plasticity implementation.

  protected:

    double rho0;

};
#endif
