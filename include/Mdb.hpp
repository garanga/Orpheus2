/*
 * Mdb.hpp
 *
 *  Created on: Jan 17, 2017
 *      Author: pavel
 */

#ifndef MDB_HPP_
#define MDB_HPP_

#include <string>
#include <vector>


class Model;
class Job;

class Mdb
{

public:

	Mdb();

   ~Mdb();

	Model*
	getModel(int No) const;

	Model*
	getModel(std::string name) const;

	Job*
	getJob(int No) const;

	Job*
	getJob(std::string name) const;

	Model*
	createModel(std::string myName);

	Job*
	createJob(std::string name, Model* model);

private:

	std::vector<Model*>	models_;
	std::vector<Job*>   jobs_;

};

#endif /* MDB_HPP_ */
