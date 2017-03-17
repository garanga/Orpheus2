/*
 * Isotropic.hpp
 *
 *  Created on: Jan 17, 2017
 *      Author: pavel
 */

#ifndef ISOTROPIC_HPP_
#define ISOTROPIC_HPP_

#include "Material.hpp"

#include <string>

class Isotropic : public Material
{

public:

	Isotropic(std::string name, double young, double poisson);

   ~Isotropic();

//	std::string
//	getName() const override;

	double
	getYoung() const override;

	double
	getPoisson() const override;

private:

	double      young_;
	double      poisson_;

};

#endif /* ISOTROPIC_HPP_ */
