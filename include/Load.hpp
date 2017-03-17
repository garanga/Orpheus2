/*
 * Load.hpp
 *
 *  Created on: Feb 3, 2017
 *      Author: pavel
 */

#ifndef LOAD_HPP_
#define LOAD_HPP_

#include <iostream>
#include <vector>

class Job;
enum class ConcentratedLoadType;

class Load
{
private:
	std::string mName;
public:

	std::vector<int> region;
	ConcentratedLoadType type;
	double* value;

	Load(std::string myName, std::vector<int> myRegion, ConcentratedLoadType myType, double* myValue);
	virtual ~Load();
	friend class Job;
};



#endif /* LOAD_HPP_ */
