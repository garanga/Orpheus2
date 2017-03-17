/*
 * Line.cpp
 *
 *  Created on: Feb 27, 2017
 *      Author: pavel
 */

#include "Line.hpp"

#include <vector>

Line::Line(int point1No, int point2No)
{
	connect_.push_back(point1No);
	connect_.push_back(point2No);
}

Line::~Line()
{

}

void
Line::setDivisions(int divisions)
{
	divisions_ = divisions;
}

std::vector<int>
Line::getConnect() const
{
	return connect_;
}

int
Line::getDivisions() const
{
	return divisions_;
}



