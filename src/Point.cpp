/*
 * Point.cpp
 *
 *  Created on: Feb 27, 2017
 *      Author: pavel
 */

#include "Point.hpp"

#include <vector>

Point::Point(double x1, double x2)
{
	coord_.push_back(x1);
	coord_.push_back(x2);
}

Point::~Point()
{

}

std::vector<double>
Point::getCoord() const
{
	return coord_;
}


