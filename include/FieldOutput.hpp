/*
 * FieldOutput.hpp
 *
 *  Created on: Mar 2, 2017
 *      Author: pavel
 */

#ifndef FIELDOUTPUT_HPP_
#define FIELDOUTPUT_HPP_

#include <vector>

enum class FieldType;

class FieldOutput
{

public:

	FieldOutput(FieldType, std::vector<double>);
   ~FieldOutput();

    FieldType getFieldType() const;

    std::vector<double> getData() const;

private:

   FieldType fieldType_;
   std::vector<double> data_;



};



#endif /* FIELDOUTPUT_HPP_ */
