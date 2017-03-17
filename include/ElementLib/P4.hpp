/*
 * P4.hpp
 *
 *  Created on: Jan 19, 2017
 *      Author: pavel
 */

#ifndef ELEMENTLIB_P4_HPP_
#define ELEMENTLIB_P4_HPP_

#include "ElementType.hpp"

#include <Eigen/Dense>

class Material;

//! A class representing a 4-node plane isoparametric element
/*! The shape functions for the element are defined in the
 *  reference element so that different elements have the same shape functions.
 *  Let \f$ \bm{u}_I = \{u_{I1},u_{I2}\}^T \f$
 *  be the displacement vector at node \f$ I, I=1,...,4 \f$ and
 *  let \f$ \bm{\xi}=\{\xi,\eta \}^T \f$ be the local coordinates
 *  retated to the reference element.
 *  Then for isoparametric element the following mapping is valid:
 *  \f[
 *  	\bm{x}(\bm{\xi})=\sum_{I=1}^{4} N_I(\bm{\xi}) \bm{x}_I
 *  \f]
 */

class P4 : public ElementType
{

public:

	P4(Material* material, bool isPlaneStrain=false);
   ~P4();

	// A method returns a stiffness matrix K: the integral over the domain occupied by an element
	// of a quantity B(i,...).transpose()*D*B(j,...)
	Eigen::Matrix<double,2,2>
	calcLocK(int i, int j, Eigen::Matrix<double,2,4> &nodesCoordGlob) const override;

	// A method returns an elastic matrix D
	Eigen::Matrix<double,3,3>
	getD() const;

	// A method updates the output variables and residual force
	void
	update(Eigen::Matrix<double,2,4> &nodesCoordGlob, Eigen::VectorXd &u, Eigen::MatrixXd &sigma, Eigen::VectorXd &force) const;

	void
	updateForce(Eigen::Matrix<double,2,4> &nodesCoordGlob, Eigen::VectorXd &u, Eigen::VectorXd &force) const;



private:

	double* mAlpha;
	double* mBeta;
	double* mGamma;
	double* mDelta;

	Eigen::Matrix<double,3,3> mD;

	Eigen::Matrix<double,2,4> mQuadraturePoints;
	Eigen::Matrix<double,1,4> mQuadratureWeights;

	//! A method returns a shape function corresponded to node \f$ I \f$
	//! calculated in a point \f$ (\xi,\eta) \f$ in the local coordinates
	double calcShapeFuncLoc(int i, double xi, double eta) const;

	// A method returns a shape function of node <<i>> local derivatives
	// calculated in a local point (xi,eta): ||dNi/d(xi),dNi/d(eta)||(xi,eta)
	Eigen::Matrix<double,1,2> calcShapeFuncLocDerLoc(int i, double xi, double eta) const;

	// A method returns the Jacobian matrix d(x,y)/d(xi,eta)
	// calculated in a local point (xi,eta): ||dx/dxi,dx/deta||
	//                                       ||dy/dxi,dy/deta||(xi,eta)
	// nodesCoordGlob is a matrix containing element nodes global coordinates:
	// || x1,x2,x3,x4
	//    y1,y2,y3,y4 ||
	Eigen::Matrix<double,2,2> calcJacobianMatrix(double xi, double eta, Eigen::Matrix<double,2,4>& nodesCoordGlob) const;

	// A method returns a shape function of node <<i>> global derivatives
	// calculated in a local point (xi,eta): ||dNi/d(x),dNi/d(y)||(xi,eta)
	Eigen::Matrix<double,1,2> calcShapeFuncGlobDerLoc(int i, double xi, double eta, Eigen::Matrix<double,2,4>& nodesCoordGlob) const;

	//! A method returns a matrix \f$ \bm{B}_I \f$ composed of
	//! global derivatives of a shape function \f$ N_I \f$
	Eigen::Matrix<double,3,2> calcBi(int i, double xi, double eta, Eigen::Matrix<double,2,4>& nodesCoordGlob) const;

	//! A method returns a matrix \f$ \bm{B} \f$ composed of
	//! global derivatives of shape functions \f$ N_I, I=1,\ldots,4 \f$
	Eigen::Matrix<double,3,8> calcB(double xi, double eta, Eigen::Matrix<double,2,4>& nodesCoordGlob) const;

};

#endif /* ELEMENTLIB_P4_HPP_ */
