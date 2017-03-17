/*
 * Step.cpp
 *
 *  Created on: Feb 9, 2017
 *      Author: pavel
 */

#include "Step.hpp"

#include <vector>

Step::Step(std::string name_)
	: name(name_)
{

}

Step::~Step()
{

}

void
Step::addFieldOutputRequest(FieldType fieldType_)
{
	fieldOutputRequest.push_back(fieldType_);
}

void
Step::addFieldOutputRequest(std::vector<FieldType> fieldTypes_)
{
	fieldOutputRequest.insert(fieldOutputRequest.end(),
			                  fieldTypes_.begin(), fieldTypes_.end());
}

std::string
Step::getName() const
{
	return name;
}

std::vector<FieldType> Step::getFieldOutputRequest() const
{
	return fieldOutputRequest;
}


