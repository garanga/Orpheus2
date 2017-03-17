/*
 * ElementEdge.hpp
 *
 *  Created on: Feb 8, 2017
 *      Author: pavel
 */

#ifndef ELEMENTEDGE_HPP_
#define ELEMENTEDGE_HPP_


#include <vector>

class Element;

class ElementEdge
{
private:
	int mId;
	int* mConnectivity;
	std::vector<Element*> mElements;
public:
	ElementEdge(int myId, int* myConnectivity, Element* myElement);
	int id() const;
	int* connectivity() const;
	std::vector<Element*> elements() const;
};


#endif /* ELEMENTEDGE_HPP_ */
