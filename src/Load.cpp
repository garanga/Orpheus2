/*
 * Load.cpp
 *
 *  Created on: Feb 3, 2017
 *      Author: pavel
 */

#include "Load.hpp"

Load::Load(std::string myName, std::vector<int> myRegion, ConcentratedLoadType myType, double* myValue)
{
	mName = myName;
	type = myType;
	region = myRegion;
	value = myValue;
}

Load::~Load()
{

}


