/*
 * ElementType.cpp
 *
 *  Created on: Jan 18, 2017
 *      Author: pavel
 */

#include "ElementType.hpp"

#include <string>

ElementType::ElementType()
{

}

ElementType::~ElementType()
{

}

std::string
ElementType::getName() const
{
	return name_;
}

int
ElementType::getNodesNum() const
{
	return nodesNum_;
}
