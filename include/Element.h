/*
 * Element.hpp
 *
 *  Created on: Jan 18, 2017
 *      Author: pavel
 */

#ifndef ELEMENT_H_
#define ELEMENT_H_

class ElementType;

class Element
{

public:

	Element(int, int*, ElementType*);
	int id() const;
	int* getConnect() const;
	ElementType* type() const;

private:

	int mId;
	int* mConnectivity;
	ElementType* mType;

};

#endif /* ELEMENT_H_ */
