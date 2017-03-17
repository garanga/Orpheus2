/*
 * Mdb.cpp
 *
 *  Created on: Jan 17, 2017
 *      Author: pavel
 */

#include "Mdb.hpp"

#include "Model.hpp"
#include "Job.hpp"

#include <algorithm>
#include <iostream>

Mdb::Mdb()
{
	std::cout << "An empty models database is created" << std::endl;
}

Mdb::~Mdb()
{
	for (auto model = models_.begin(); model < models_.end(); ++model)
	{
		delete *model;
	}

	for (auto job = jobs_.begin(); job < jobs_.end(); ++job)
	{
		delete *job;
	}
}

Model*
Mdb::getModel(int No) const
{
	return models_[No];
}

Model*
Mdb::getModel(std::string name) const
{
	auto it = std::find_if(models_.begin(),models_.end(),[name](Model* model)
	          {
		 	  	  return model->getName() == name;
	          });
	return *it;
}

Job*
Mdb::getJob(int No) const
{
	return jobs_[No];
}

Job*
Mdb::getJob(std::string name) const
{
	auto it = std::find_if(jobs_.begin(),jobs_.end(),[name](Job* job)
	          {
		 	  	  return job->getName() == name;
	          });
	return *it;
}

Model*
Mdb::createModel(std::string myName)
{
	Model* model = new Model(myName);
	models_.push_back(model);
	return model;
}

Job*
Mdb::createJob(std::string name, Model* model)
{
	Job* job = new Job(name, model);
	jobs_.push_back(job);
	return job;
}



