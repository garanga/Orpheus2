/*
 * Point.hpp
 *
 *  Created on: Feb 27, 2017
 *      Author: pavel
 */

#ifndef POINT_HPP_
#define POINT_HPP_

#include <vector>

class Point
{
public:

	Point(double x1, double x2);
   ~Point();

	std::vector<double>
	getCoord() const;

private:

	std::vector<double> coord_;

};



#endif /* POINT_HPP_ */
