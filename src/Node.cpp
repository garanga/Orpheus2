/*
 * Node.cpp
 *
 *  Created on: Jan 18, 2017
 *      Author: pavel
 */

#ifndef NODE_CPP_
#define NODE_CPP_

#include "Node.h"

Node::Node(int myId, double* myCoord)
{
	mId = myId;
	mCoord = myCoord;
}

int Node::getId() const
{
	return mId;
}

double* Node::getCoord() const
{
	return mCoord;
}



#endif /* NODE_CPP_ */
