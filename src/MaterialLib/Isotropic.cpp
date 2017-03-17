/*
 * Isotropic.cpp
 *
 *  Created on: Jan 17, 2017
 *      Author: pavel
 */

#include "Isotropic.hpp"

#include <iostream>
#include <string>


//std::string
//Material::getName() const
//{
//	return name_;
//}

Isotropic::Isotropic(std::string name, double young, double poisson)
	: Material(name),
	  young_(young), poisson_(poisson)
{
	std::cout << "The isotropic material '" << name << "' is created" << std::endl;
}


Isotropic::~Isotropic()
{

}





double
Isotropic::getYoung() const
{
	return young_;
}


double Isotropic::getPoisson() const
{
	return poisson_;
}
