/*
 * Constraint.hpp
 *
 *  Created on: Feb 1, 2017
 *      Author: pavel
 */

#ifndef CONSTRAINT_HPP_
#define CONSTRAINT_HPP_


#include <iostream>
#include <vector>

class Job;
enum class DisplacementConstraintType;

class Constraint
{
private:
	std::string mName;
public:

	DisplacementConstraintType type;
	std::vector<int> region;
	double* value;

	Constraint(std::string myName, std::vector<int> myRegion, DisplacementConstraintType myType, double* myValue);
	virtual ~Constraint();
	friend class Job;
};



#endif /* CONSTRAINT_HPP_ */
