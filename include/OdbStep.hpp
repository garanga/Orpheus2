/*
 * OdbStep.hpp
 *
 *  Created on: Mar 2, 2017
 *      Author: pavel
 */

#ifndef ODBSTEP_HPP_
#define ODBSTEP_HPP_

#include <vector>

class Step;
class OdbFrame;
enum class FieldType;

class OdbStep
{

public:

	OdbStep(Step*);
   ~OdbStep();

    OdbFrame* createFrame();

    int getFieldsOutputRequestNum() const;

    int getFrameNum() const;

    OdbFrame* getFrame(int) const;

private:

    std::vector<FieldType> fieldOutputRequest;
	std::vector<OdbFrame*> frames_;

};



#endif /* ODBSTEP_HPP_ */
