/*
 * LinearSolver.hpp
 *
 *  Created on: Jan 27, 2017
 *      Author: pavel
 */

#ifndef LINEARSOLVER_HPP_
#define LINEARSOLVER_HPP_


#include <string>

class Model;
class Odb;

#include "Eigen/Sparse"

class LinearSolver
{

public:

	LinearSolver(Model*, Odb*);

	void solve();

private:

	Model* model_;
	Odb*   odb_;

	void update(Eigen::VectorXd &u, Eigen::MatrixXd &sigma, Eigen::VectorXd & force) const;
	Eigen::SparseMatrix<double> calcGlobK() const;
	void updateForce(Eigen::VectorXd &u, Eigen::VectorXd &force) const;

};


#endif /* LINEARSOLVER_HPP_ */
