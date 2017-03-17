/*
 * Job.hpp
 *
 *  Created on: Jan 31, 2017
 *      Author: pavel
 */

#ifndef JOB_HPP_
#define JOB_HPP_

#include <string>

class Model;
class Odb;

class Job
{

public:

	Job(std::string name, Model* model);
   ~Job();

    std::string
	getName() const;

	Odb*
	Submit();

private:

	std::string name_;
	Model*      model_;
	Odb*        odb_;

};


#endif /* JOB_HPP_ */
