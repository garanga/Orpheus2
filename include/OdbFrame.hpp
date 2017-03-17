/*
 * OdbFrame.hpp
 *
 *  Created on: Mar 2, 2017
 *      Author: pavel
 */

#ifndef ODBFRAME_HPP_
#define ODBFRAME_HPP_

#include <vector>

class FieldOutput;

class OdbFrame
{

public:

	OdbFrame();
   ~OdbFrame();

    void
	addFieldOutput(FieldOutput*);

    FieldOutput*
	getFieldOutput(int fieldOutputNo) const;

private:

   std::vector<FieldOutput*> fieldOutputs_;

};



#endif /* ODBFRAME_HPP_ */
