/*
 * FieldOutput.cpp
 *
 *  Created on: Mar 2, 2017
 *      Author: pavel
 */

#include "FieldOutput.hpp"

#include "includes.hpp"
#include <vector>

FieldOutput::FieldOutput(FieldType fieldType, std::vector<double> data)
	: fieldType_(fieldType), data_(data)
{

}

FieldOutput::~FieldOutput()
{

}

FieldType FieldOutput::getFieldType() const
{
	return fieldType_;
}

std::vector<double>
FieldOutput::getData() const
{
	return data_;
}


