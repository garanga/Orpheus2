/*
 * Constraint.cpp
 *
 *  Created on: Feb 1, 2017
 *      Author: pavel
 */

#include "Constraint.hpp"

Constraint::Constraint(std::string myName, std::vector<int> myRegion, DisplacementConstraintType myType, double* myValue)
{
	mName = myName;
	type = myType;
	region = myRegion;
	value = myValue;

//	std::cout << std::endl;
//	std::cout << mValue[0] << "   " << mValue[1] << std::endl;
//	std::cout << mRegion[0] << std::endl;
//	std::cout << std::endl;

}

Constraint::~Constraint()
{

}


