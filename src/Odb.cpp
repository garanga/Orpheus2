/*
 * Odb.cpp
 *
 *  Created on: Jan 31, 2017
 *      Author: pavel
 */

#include "Odb.hpp"

#include "includes.hpp"

#include <iostream>
#include <fstream>
#include <map>

Odb::Odb(Model* model)
{
	parts_  = model->parts;
	bodies_ = model->bodies;

	// A number of frames is unknown so they are created during the simulation

	std::cout << "An output database is created" << std::endl;
}

Odb::~Odb()
{

}

OdbStep*
Odb::createStep(Step* step)
{
	OdbStep* odbStep = new OdbStep(step);
	steps_.push_back(odbStep);
	return odbStep;
}

int
Odb::getStepNum() const
{
	return steps_.size();
}

OdbStep*
Odb::getStep(int stepNo) const
{
	return steps_[stepNo];
}

OdbStep*
Odb::getStep(Constant constant) const
{
	switch(constant)
	{
		case Constant::LAST:
		{
			return steps_.back();
		}
	}
}

void
saveOdb(const Odb* odb, std::string filePath)
{

	std::ofstream file(filePath);

	file << odb->parts_[0]->mesh->nodesNum << "\n";

	file.precision(3);
	file.setf(std::ios::scientific);

	for (auto node = (odb->parts_[0]->mesh->nodes).begin();
	          node < (odb->parts_[0]->mesh->nodes).end(); ++node)
	{
		double* coord = (*node)->getCoord();
		file << std::showpos << coord[0] << " " << coord[1] << "\n";
		delete[] coord;
	}

	file << std::noshowpos << odb->parts_[0]->mesh->elementsNum << "\n";

	for (auto element = (odb->parts_[0]->mesh->elements).begin();
	          element < (odb->parts_[0]->mesh->elements).end(); ++element)
	{
		int* connect = (*element)->getConnect();
		file << connect[0] << " "
			 << connect[1] << " "
			 << connect[2] << " "
			 << connect[3] << "\n";
		delete[] connect;
	}

	file << (odb->steps_).size() << "\n";

	for (auto step = odb->steps_.begin(); step < odb->steps_.end(); ++step)
	{
		file << (*step)->getFieldsOutputRequestNum() << "\n";
		file << (*step)->getFrameNum() << "\n";
		for (int i = 0; i < (*step)->getFrameNum(); ++i)
		{
			OdbFrame* frame = (*step)->getFrame(i);
			FieldOutput* fieldOutput = frame->getFieldOutput(0);

			FieldType fieldType = fieldOutput->getFieldType();
			file << static_cast<int>(fieldType) << "\n";

			std::vector<double> data;
			data = fieldOutput->getData();

			for (auto it = data.begin(); it < data.end(); ++it)
			{
				file << std::showpos << *it << "\n";
			}







//			delete fieldOutput;
//			delete frame;

		}
	}



//	file << bodies_[0]->getMesh()->nodesNum << "\n";
//
//	file.precision(3);
//	file.setf(std::ios::showpos);
//	file.setf(std::ios::scientific);
//
//	for (auto node = (bodies_[0]->getMesh()->nodes).begin(); node < (bodies_[0]->getMesh()->nodes).end(); ++node)
//	{
//		double* coord = (*node)->getCoord();
//		file << coord[0] << " " << coord[1] << "\n";
//		coord = nullptr;
//	}

	file.close();
}


