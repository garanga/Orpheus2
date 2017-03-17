/*
 * Body.cpp
 *
 *  Created on: Feb 27, 2017
 *      Author: pavel
 */

#include "Body.hpp"

#include "includes.hpp"

#include "triangle.h"
#include "algorithm"
#include <string>

Body::Body(std::string name)
	: name_(name)
{

}

Body::~Body()
{

}

void
Body::addPoint(double x1, double x2)
{
	Point* point = new Point(x1, x2);
	points_.push_back(point);
}

void
Body::addLine(int point1No, int point2No)
{
	Line* line = new Line(point1No, point2No);
	lines_.push_back(line);
}

void
Body::setLineDivisions(int lineNo, int divisions)
{
	lines_[lineNo]->setDivisions(divisions);
}

Mesh*
Body::createMesh(ElementType* type)
{
	triangulateio* in     = new(triangulateio);
	triangulateio* out    = new(triangulateio);
	triangulateio* vorout = new(triangulateio);

	// In

	for (auto line = lines_.begin(); line < lines_.end(); ++line)
	{
		int divisions = (*line)->getDivisions();

		if (divisions > 1)
		{

			std::vector<int> connect = (*line)->getConnect();

			Point* point1 = points_[connect[0]];
			Point* point2 = points_[connect[1]];

			std::vector<double> coord1 = point1->getCoord();
			std::vector<double> coord2 = point2->getCoord();

			double length = std::sqrt(std::pow(coord2[0]-coord1[0],2)+std::pow(coord2[1]-coord1[1],2));

			std::vector<double> direction {(coord2[0]-coord1[0])/length, (coord2[1]-coord1[1])/length};

			for (int i = 0; i < divisions - 1; ++i)
			{
				addPoint(coord1[0] + direction[0]*length/divisions*(i + 1), coord1[1] + direction[1]*length/divisions*(i + 1));
			}

		}

	}

	int numberofpoints = points_.size();

	std::vector<double> v1;
	for (auto point = points_.begin(); point < points_.end(); ++point)
	{
		std::vector<double> coord = (*point)->getCoord();
		v1.insert(v1.end(), coord.begin(), coord.end());
	}

	double* pointlist = &v1[0];

	int numberofsegments = lines_.size();

	std::vector<int> v2;
	for (auto segment = lines_.begin(); segment < lines_.end(); ++segment)
	{
		std::vector<int> connect = (*segment)->getConnect();
		std::for_each(connect.begin(), connect.end(), [](int& n) {n += 1.0;});
		v2.insert(v2.end(), connect.begin(), connect.end());
	}

	int* segmentlist = &v2[0];

	int numberofpointattributes = 0;
	double* pointattributelist;
	pointattributelist = new double[numberofpoints*numberofpointattributes];

	int* pointmarkerlist;
	pointmarkerlist = nullptr;

	int* segmentmarkerlist;
	segmentmarkerlist = nullptr;

	int numberofholes = 0;
	double* holelist;
	holelist = new double[numberofholes*2];

	int numberofregions = 0;
	double* regionlist;
	regionlist = new double[numberofregions*4];

	in->numberofpoints = numberofpoints;
	in->pointlist = pointlist;
	in->numberofsegments = numberofsegments;
	in->segmentlist = segmentlist;

	in->pointattributelist = pointattributelist;
	in->pointmarkerlist = pointmarkerlist;
	in->numberofpointattributes = numberofpointattributes;
	in->segmentmarkerlist = segmentmarkerlist;
	in->holelist = holelist;
	in->numberofholes = 0;
	in->regionlist = regionlist;
	in->numberofregions = numberofregions;

	// Out
	out->pointlist = nullptr;
	out->pointmarkerlist = nullptr;
	out->trianglelist = nullptr;
	out->neighborlist = nullptr;
	out->segmentlist = nullptr;
	out->segmentmarkerlist = nullptr;

	triangulate((char*)"a0.05",in,out,vorout);

	mesh_ = new Mesh;

	mesh_->nodesNum = out->numberofpoints;

	for (int i = 0; i < out->numberofpoints; ++i)
	{
		double* coord = new double[2];
		coord[0] = out->pointlist[2*i+0];
		coord[1] = out->pointlist[2*i+1];
		Node* node = new Node(i,coord);
		mesh_->nodes.push_back(node);
		coord = nullptr;
		node = nullptr;
	}


	std::cout << out->numberoftriangles << std::endl;
	std::cout << out->numberofpoints << std::endl;

}

int
Body::getLinesNum() const
{
	return lines_.size();
}

std::vector<Line*>
Body::getLines() const
{
	return lines_;
}

Mesh*
Body::getMesh() const
{
	return mesh_;
}
