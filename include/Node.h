/*
 * Node.h
 *
 *  Created on: Jan 18, 2017
 *      Author: pavel
 */

#ifndef NODE_H_
#define NODE_H_

class Node
{
private:
	int mId;
	double* mCoord;
public:
	Node(int myId, double* myCoord);
	int getId() const;
	double* getCoord() const;
};

#endif /* NODE_H_ */
