/*
 * ElementEdge.cpp
 *
 *  Created on: Feb 8, 2017
 *      Author: pavel
 */

#include "ElementEdge.hpp"

#include <algorithm>
#include <vector>
#include "Element.h"

ElementEdge::ElementEdge(int myId, int* myConnectivity, Element* myElement)
	: mId(myId),
	  mConnectivity(myConnectivity)
{

	if (std::find(mElements.begin(), mElements.end(), myElement) == mElements.end())
	{
		mElements.push_back(myElement);
	}

}

int ElementEdge::id() const
{
	return mId;
}

int* ElementEdge::connectivity() const
{
	return mConnectivity;
}

std::vector<Element*> ElementEdge::elements() const
{
	return mElements;
}


