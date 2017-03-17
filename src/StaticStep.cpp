/*
 * StaticStep.cpp
 *
 *  Created on: Jan 27, 2017
 *      Author: pavel
 */

#include "StaticStep.hpp"

#include "includes.hpp"
#include <string>


StaticStep::StaticStep(std::string myName,
                       double myTimeBegin,
					   double myTimeEnd,
					   double myTimeIncrement,
					   double myLoadFactorBegin,
					   double myLoadFactorEnd)
	: Step(myName),
	  timeBegin(myTimeBegin), timeEnd(myTimeEnd), timeIncrement(myTimeIncrement), loadFactorBegin(myLoadFactorBegin), loadFactorEnd(myLoadFactorEnd)
{

}

StaticStep::~StaticStep()
{

}

DisplacementConstraint* StaticStep::CreateDispacementConstraint(std::string myName, std::vector<int> myRegion, DisplacementConstraintType myType, double* myValue)
{
	DisplacementConstraint* displacementConstraint = new DisplacementConstraint(myName, myRegion, myType, myValue);
	constraints.push_back(displacementConstraint);
	return displacementConstraint;
}

ConcentratedLoad* StaticStep::CreateConcentratedLoad(std::string myName, std::vector<int> myRegion, ConcentratedLoadType myType, double* myValue)
{
	ConcentratedLoad* concentratedLoad = new ConcentratedLoad(myName, myRegion, myType, myValue);
	loads.push_back(concentratedLoad);
	return concentratedLoad;
}
