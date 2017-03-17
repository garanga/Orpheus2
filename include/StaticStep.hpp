/*
 * StaticStep.hpp
 *
 *  Created on: Jan 27, 2017
 *      Author: pavel
 */

#ifndef STATICSTEP_HPP_
#define STATICSTEP_HPP_


#include "Step.hpp"
#include <string>
#include <vector>

class Constraint;
class DisplacementConstraint;
class Load;
class ConcentratedLoad;

enum class DisplacementConstraintType;
enum class ConcentratedLoadType;

class StaticStep  : public Step
{
private:

public:

	double timeBegin;
	double timeEnd;
	double timeIncrement;
	double loadFactorBegin;
	double loadFactorEnd;

	std::vector<Constraint*> constraints;
	std::vector<Load*> loads;

	// Specialized constructor
	StaticStep(std::string myName,
               double myTimeBegin,
			   double myTimeEnd,
			   double myTimeIncrement,
			   double myLoadFactorBegin,
			   double myLoadFactorEnd);

	~StaticStep();

	// A method creating a DisplacementConstraint object
	DisplacementConstraint*
	CreateDispacementConstraint(std::string myName, std::vector<int> myRegion, DisplacementConstraintType myType, double* myValue);

	// A method creating a ConcentratedLoad object
	ConcentratedLoad*
	CreateConcentratedLoad(std::string myName, std::vector<int> myRegion, ConcentratedLoadType myType, double* myValue);

};


#endif /* STATICSTEP_HPP_ */
