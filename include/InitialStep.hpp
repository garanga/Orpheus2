/*
 * InitialStep.hpp
 *
 *  Created on: Feb 9, 2017
 *      Author: pavel
 */

#ifndef INITIALSTEP_HPP_
#define INITIALSTEP_HPP_


#include "Step.hpp"

#include <iostream>
#include <vector>

class DisplacementConstraint;
class ConcentratedLoad;

enum class DisplacementConstraintType;
enum class ConcentratedLoadType;

class InitialStep : public Step
{
private:

public:
	// Specialized constructor
	InitialStep(std::string myName);
	~InitialStep();
};


#endif /* INITIALSTEP_HPP_ */
