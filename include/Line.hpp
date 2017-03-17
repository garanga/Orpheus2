/*
 * Line.hpp
 *
 *  Created on: Feb 27, 2017
 *      Author: pavel
 */

#ifndef LINE_HPP_
#define LINE_HPP_

#include <vector>

class Line
{

public:

	Line(int point1No, int point2No);
   ~Line();

    void
	setDivisions(int divisions);

	std::vector<int>
	getConnect() const;

	int
	getDivisions() const;

private:

	std::vector<int> connect_;
	int divisions_ = 1;

};



#endif /* LINE_HPP_ */
