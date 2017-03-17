/*
 * Mesh.hpp
 *
 *  Created on: Jan 18, 2017
 *      Author: pavel
 */

#ifndef MESH_HPP_
#define MESH_HPP_


#include <vector>

class Node;
class Element;
class ElementEdge;

struct Mesh
{
    int nodesNum;
	int elementsNum;
	int elementEdgesNum;

	std::vector <Node*> nodes;
	std::vector <Element*> elements;
	std::vector <ElementEdge*> elementEdges;

};


#endif /* MESH_HPP_ */
