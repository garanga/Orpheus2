/*
 * Part.hpp
 *
 *  Created on: Jan 17, 2017
 *      Author: pavel
 */

#ifndef PART_HPP_
#define PART_HPP_


#include <string>
#include <vector>

class Material;
class ElementType;
struct Mesh;

//! A test class1
/*!
 *
 * A more detailed description1 \f$ \frac{1}{2}  \f$
 *
 */

class Part
{

public:

	Part(std::string name);

	void
	setSizes(double* sizes);

	void
	setDivisions(int* divisions);

	// A method creating a mesh on part
	Mesh*
	CreateMesh(ElementType* type, bool writeMesh = true);

	std::string
	getName() const;

	double*
	getSizes() const;

	int*
	getDivisions() const;

	Mesh* mesh = nullptr;

private:

	std::string name_;
	double*     sizes_ = nullptr;
	int*        divisions_ = nullptr;

};



#endif /* PART_HPP_ */
