/*
 * Model.cpp
 *
 *  Created on: Jan 17, 2017
 *      Author: pavel
 */

#include "Model.hpp"

#include "includes.hpp"

#include <iostream>

Model::Model(std::string name)
	: name_(name)
{

	std::cout << "The model '" << name_ << "' is created" << std::endl;

	// Creating an initial step
//	Step* step = new InitialStep("Initial");
//	steps.push_back(step);

}

Model::~Model()
{
	for (auto body = bodies.begin(); body < bodies.end(); ++body)
	{
		delete *body;
	}

	for (auto part = parts.begin(); part < parts.end(); ++part)
	{
		delete *part;
	}

	for (auto material = materials.begin(); material < materials.end(); ++material)
	{
		delete *material;
	}

	for (auto step = steps.begin(); step < steps.end(); ++step)
	{
		delete *step;
	}
}

std::string
Model::getName() const
{
	return name_;
}

Body*
Model::createBody(std::string name)
{
	Body* body = new Body(name);
	bodies.push_back(body);
	return body;
}

Part*
Model::createPart(std::string name)
{
	Part* part = new Part(name);
	parts.push_back(part);
	return part;
}

Material*
Model::createIsotropic(std::string name, double young, double poisson)
{
	Isotropic* material = new Isotropic(name, young, poisson);
	materials.push_back(material);
	return material;
}

//Step* Model::CreateStaticStep(std::string myName,
//                              double myTimeBegin,
//							  double myTimeEnd,
//							  double myTimeIncrement,
//							  double myLoadFactorBegin,
//							  double myLoadFactorEnd)
//{
//	Step* step = new StaticStep(myName, myTimeBegin, myTimeEnd, myTimeIncrement, myLoadFactorBegin, myLoadFactorEnd);
//	steps.push_back(step);
//	return step;
//}

StaticStep*
Model::createStaticStep(std::string name,
                        double timeBegin,
                        double timeEnd,
						double timeIncrement,
						double loadFactorBegin,
						double loadFactorEnd)
{
	StaticStep* step = new StaticStep(name, timeBegin, timeEnd, timeIncrement, loadFactorBegin, loadFactorEnd);
	steps.push_back(step);
	return step;
}
