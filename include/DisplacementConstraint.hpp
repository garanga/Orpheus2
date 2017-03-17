/*
 * DisplacementConstraint.hpp
 *
 *  Created on: Jan 30, 2017
 *      Author: pavel
 */

#ifndef DISPLACEMENTCONSTRAINT_HPP_
#define DISPLACEMENTCONSTRAINT_HPP_

#include <string>
#include <vector>

#include "Constraint.hpp"

enum class DisplacementConstraintType;

class DisplacementConstraint : public Constraint
{

public:
	// Specialized constructor
	DisplacementConstraint(std::string myName, std::vector<int> myRegion, DisplacementConstraintType myType, double* myValue);

	// Virtual destructor
	virtual ~DisplacementConstraint();

};


#endif /* DISPLACEMENTCONSTRAINT_HPP_ */
