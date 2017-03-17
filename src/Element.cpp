/*
 * Element.cpp
 *
 *  Created on: Jan 18, 2017
 *      Author: pavel
 */

#include "Element.h"

#include "ElementLib/ElementType.hpp"

Element::Element(int myId, int* myConnectivity, ElementType* myType)
{
	mId = myId;
	mConnectivity = myConnectivity;
	mType = myType;
}

int Element::id() const
{
	return mId;
}

int* Element::getConnect() const
{
	return mConnectivity;
}

ElementType* Element::type() const
{
	return mType;
}


