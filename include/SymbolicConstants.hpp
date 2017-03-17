/*
 * SymbolicConstants.hpp
 *
 *  Created on: Feb 1, 2017
 *      Author: pavel
 */

#ifndef SYMBOLICCONSTANTS_HPP_
#define SYMBOLICCONSTANTS_HPP_

enum class DisplacementConstraintType
{
	UX = 1 << 0,
	UY = 1 << 1,
	UXY = UX | UY
};

enum class ConcentratedLoadType
{
	FX = 1 << 0,
	FY = 1 << 1,
	FXY = FX | FY
};


enum class StaticStepConstants
{
	// matrixSolver
	DIRECT,
	ITERATIVE,
	//

};

enum class  Constant
{
	LAST
};

enum class FieldType
{
	UX	=  1,
	UY	=  2,
	UZ	=  3,
	SX	=  4,
	SY	=  5,
	SZ	=  6,
	SYZ	=  7,
	SXZ	=  8,
	SXY	=  9,
	U 	= 10
};


#endif /* SYMBOLICCONSTANTS_HPP_ */
