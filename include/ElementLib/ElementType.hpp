/*
 * ElementType.hpp
 *
 *  Created on: Jan 18, 2017
 *      Author: pavel
 */

#ifndef ELEMENTTYPE_HPP_
#define ELEMENTTYPE_HPP_

#include <string>

#include <Eigen/Dense>

//! A base class for all the classes representing finite elements
/*
 * Detailed description
 */
class ElementType
{

public:

	ElementType();

	virtual
   ~ElementType();

	virtual Eigen::Matrix<double,2,2>
	calcLocK(int i, int j, Eigen::Matrix<double,2,4> &nodesCoordGlob) const = 0;

	std::string
	getName() const;

	int
	getNodesNum() const;



	virtual void
	update(Eigen::Matrix<double,2,4> &nodesCoordGlob, Eigen::VectorXd &u, Eigen::MatrixXd &sigma, Eigen::VectorXd &force) const = 0;



	virtual void
	updateForce(Eigen::Matrix<double,2,4> &nodesCoordGlob, Eigen::VectorXd &u, Eigen::VectorXd &force) const = 0;

protected:

	std::string name_ = "";
	int dimension_ = 0;
	int nodesNum_  = 0;

};




#endif /* ELEMENTTYPE_HPP_ */
