/*
 * Step.hpp
 *
 *  Created on: Feb 9, 2017
 *      Author: pavel
 */

#ifndef STEP_HPP_
#define STEP_HPP_


#include <string>
#include <vector>

enum class FieldType;

class Step
{

public:

	Step(std::string);

	virtual ~Step();

	void addFieldOutputRequest(FieldType);
	void addFieldOutputRequest(std::vector<FieldType>);

	std::string getName() const;
	std::vector<FieldType> getFieldOutputRequest() const;

private:

	std::string name;
	std::vector<FieldType> fieldOutputRequest;

};


#endif /* STEP_HPP_ */
