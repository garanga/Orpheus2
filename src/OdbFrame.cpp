/*
 * OdbFrame.cpp
 *
 *  Created on: Mar 2, 2017
 *      Author: pavel
 */

#include "OdbFrame.hpp"

#include "includes.hpp"

OdbFrame::OdbFrame()
{

}

OdbFrame::~OdbFrame()
{
	for (auto fieldOutput = fieldOutputs_.begin(); fieldOutput < fieldOutputs_.end(); ++fieldOutput)
	{
		delete *fieldOutput;
	}
}

void
OdbFrame::addFieldOutput(FieldOutput* fieldOutput)
{
	fieldOutputs_.push_back(fieldOutput);
}

FieldOutput*
OdbFrame::getFieldOutput(int fieldOutputNo) const
{
	return fieldOutputs_[fieldOutputNo];
}


