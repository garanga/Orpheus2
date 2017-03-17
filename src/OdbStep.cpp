/*
 * OdbStep.cpp
 *
 *  Created on: Mar 2, 2017
 *      Author: pavel
 */

#include "OdbStep.hpp"

#include "includes.hpp"

OdbStep::OdbStep(Step* step)
{
	fieldOutputRequest = step->getFieldOutputRequest();
}

OdbStep::~OdbStep()
{

}

OdbFrame*
OdbStep::createFrame()
{
	OdbFrame* frame = new OdbFrame();
	frames_.push_back(frame);
	return frame;
}

int OdbStep::getFieldsOutputRequestNum() const
{
	return fieldOutputRequest.size();
}

int
OdbStep::getFrameNum() const
{
	return frames_.size();
}

OdbFrame*
OdbStep::getFrame(int frameNo) const
{
	return frames_[frameNo];
}


