/*
 * P4.cpp
 *
 *  Created on: Jan 19, 2017
 *      Author: pavel
 */

#include "P4.hpp"

#include <cmath>
using std::sqrt;
#include <iostream>

#include "../MaterialLib/Material.hpp"

P4:: P4(Material* myMaterial, bool isPlainStrain)
{
	name_      = "P4";
	dimension_ = 2;
	nodesNum_  = 4;

	mAlpha = new double[4];
	mAlpha[0] =  1.0; mAlpha[1] =  1.0; mAlpha[2] = 1.0; mAlpha[3] =  1.0;
	mBeta  = new double[4];
	mBeta[0]  = -1.0; mBeta[1]  =  1.0; mBeta[2]  = 1.0; mBeta[3]  = -1.0;
	mGamma = new double[4];
	mGamma[0] = -1.0; mGamma[1] = -1.0; mGamma[2] = 1.0; mGamma[3] =  1.0;
	mDelta = new double[4];
	mDelta[0] =  1.0; mDelta[1] = -1.0; mDelta[2] = 1.0; mDelta[3] = -1.0;

	double young = myMaterial->getYoung();
	double poisson = myMaterial->getPoisson();

	// exy is doubled
	if (!isPlainStrain)
	{
		mD(0,0) = 1.0    ; mD(0,1) = poisson; mD(0,2) = 0.0;
		mD(1,0) = poisson; mD(1,1) = 1.0    ; mD(1,2) = 0.0;
		mD(2,0) = 0.0    ; mD(2,1) = 0.0    ; mD(2,2) = (1.0-poisson)/2.0;
		mD *= young/(1.0-pow(poisson,2));
	}
	else
	{
		mD(0,0) = 1.0-poisson; mD(0,1) = poisson    ; mD(0,2) = 0.0;
		mD(1,0) = poisson    ; mD(1,1) = 1.0-poisson; mD(1,2) = 0.0;
		mD(2,0) = 0.0        ; mD(2,1) = 0.0        ; mD(2,2) = (1.0-2.0*poisson)/2.0;
		mD *= young/((1.0+poisson)*(1.0-2.0*poisson));
	} /* plainStrain */


	mQuadraturePoints << -1.0/sqrt(3.0), 1.0/sqrt(3.0), 1.0/sqrt(3.0), -1.0/sqrt(3.0),
			             -1.0/sqrt(3.0),-1.0/sqrt(3.0), 1.0/sqrt(3.0),  1.0/sqrt(3.0);

	mQuadratureWeights << 1.0, 1.0, 1.0, 1.0;

	Eigen::Matrix<double,2,4> x;
	x <<-1.0, 1.0, 1.0,-1.0,
	    -1.0,-1.0, 1.0, 1.0;

//	std::cout << x << std::endl;
//
//	std::cout << calcBi(0,0.25,0.25,x) << std::endl;
//	std::cout << calcBi(1,0.25,0.25,x) << std::endl;
//	std::cout << calcBi(2,0.25,0.25,x) << std::endl;
//	std::cout << calcBi(3,0.25,0.25,x) << std::endl;
//
//	std::cout << calcB(0.25,0.25,x) << std::endl;

//	Eigen::VectorXd u(8);
//	u << 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0;
//
//	Eigen::MatrixXd sigma;
//	Eigen::VectorXd force;
//
//	update(x,u,sigma,force);
//
//	std::cout << sigma << std::endl;
//
//	std::cout << force << std::endl;


}

P4::~P4()
{
	delete[] mAlpha;
	delete[] mBeta;
	delete[] mGamma;
	delete[] mDelta;
}

/*!
 *  Detailed description will be given later.
 *  @param[in]	i	The node number
 *  @param[in]	xi	The \f$ \xi \f$ coordinate in a reference frame
 *  @param[in]	eta	The \f$ \eta \f$ coordinate in a reference frame
 */
double P4::calcShapeFuncLoc(int i, double xi, double eta) const
{
	return (mAlpha[i] + mBeta[i]*xi + mGamma[i]*eta + mDelta[i]*xi*eta)/4.0;
}

Eigen::Matrix<double,1,2> P4::calcShapeFuncLocDerLoc(int i, double xi, double eta) const
{
	Eigen::Matrix<double,2,1> v;
	v << (mBeta[i] + mDelta[i]*eta)/4.0,(mGamma[i] + mDelta[i]*xi)/4.0;
	return v;
}

Eigen::Matrix<double,2,2> P4::calcJacobianMatrix(double xi, double eta, Eigen::Matrix<double,2,4>& nodesCoordGlob) const
{
	// tmp is used for the following matrix
	// ||dN0/d(xi),dN0/d(eta)
	//   dN1/d(xi),dN1/d(eta)
	//   dN2/d(xi),dN2/d(eta)
	//   dN3/d(xi),dN3/d(eta)||
	Eigen::Matrix<double,4,2> tmp;

	Eigen::Matrix<double,1,2> v;
	for (int i=0; i<nodesNum_; ++i)
	{
		v = calcShapeFuncLocDerLoc(i, xi, eta);
		tmp.row(i) = v;
	}

	return nodesCoordGlob*tmp;
}

Eigen::Matrix<double,1,2> P4::calcShapeFuncGlobDerLoc(int i, double xi, double eta, Eigen::Matrix<double,2,4>& nodesCoordGlob) const
{
	Eigen::Matrix<double,1,2> shapeFuncLocDerLoc;
	shapeFuncLocDerLoc = calcShapeFuncLocDerLoc(i, xi, eta);

//	std::cout << std::endl;
//	std::cout << "shapeFuncLocDerLoc" << std::endl;
//	std::cout << xi << "," << eta << std::endl;
//	std::cout << shapeFuncLocDerLoc << std::endl;
//	std::cout << std::endl;
//	std::cout << std::endl;

	Eigen::Matrix<double,2,2> jacobianMatrix;
	jacobianMatrix = calcJacobianMatrix(xi, eta, nodesCoordGlob);

//	std::cout << std::endl;
//	std::cout << "jacobianMatrix" << std::endl;
//	std::cout << xi << "," << eta << std::endl;
//	std::cout << jacobianMatrix << std::endl;
//	std::cout << std::endl;
//	std::cout << std::endl;

	Eigen::Matrix<double,1,2> v;
	v = shapeFuncLocDerLoc*jacobianMatrix.inverse();

	return v;
}


/*!
 * \f[ \bm{B}_I \left( \bm{\xi} \right) =
 * 		\left[
 *
 * 			\begin{array}{cc}
 * 				N_{I,1} & 0       \\
 * 				0       & N_{I,2} \\
 * 				N_{I,2} & N_{I,1}
 * 			\end{array}
 * 		\right]
 *
 * \f]
 */
Eigen::Matrix<double,3,2> P4::calcBi(int i, double xi, double eta, Eigen::Matrix<double,2,4>& nodesCoordGlob) const
{
	Eigen::Matrix<double,1,2> shapeFuncGlobDerLoc;
	shapeFuncGlobDerLoc = calcShapeFuncGlobDerLoc(i, xi, eta, nodesCoordGlob);

	Eigen::Matrix<double,3,2> Bi;
	Bi << shapeFuncGlobDerLoc(0), 0.0,
		  0.0                   , shapeFuncGlobDerLoc(1),
		  shapeFuncGlobDerLoc(1), shapeFuncGlobDerLoc(0);
	return Bi;
}

/*!
 * \f[ \bm{B} \left( \bm{\xi} \right) =
 * 		\left[
 *
 * 			\begin{array}{cccccccc}
 * 				N_{1,1} & 0       & N_{2,1} & 0       & N_{3,1} & 0       & N_{4,1} & 0       \\
 * 				0       & N_{1,2} & 0       & N_{2,2} & 0       & N_{3,2} & 0       & N_{4,2} \\
 * 				N_{1,2} & N_{1,1} & N_{2,2} & N_{2,1} & N_{3,2} & N_{3,1} & N_{4,2} & N_{4,1}
 * 			\end{array}
 * 		\right]
 *
 * \f]
 */
Eigen::Matrix<double,3,8> P4::calcB(double xi, double eta, Eigen::Matrix<double,2,4>& nodesCoordGlob) const
{
	Eigen::Matrix<double,3,2> B0, B1, B2, B3;
	B0 = calcBi(0, xi, eta, nodesCoordGlob);
	B1 = calcBi(1, xi, eta, nodesCoordGlob);
	B2 = calcBi(2, xi, eta, nodesCoordGlob);
	B3 = calcBi(3, xi, eta, nodesCoordGlob);

	Eigen::Matrix<double,3,8> B;
	B << B0, B1, B2, B3;

	return B;
}

Eigen::Matrix<double,2,2> P4::calcLocK(int i, int j, Eigen::Matrix<double,2,4> &nodesCoordGlob) const
{

	Eigen::Matrix<double,2,2> locK;
	locK = Eigen::Matrix<double,2,2>::Zero();

	double xi, eta, weight;
	Eigen::Matrix<double,3,2> Bi, Bj;
	// !!!Jacobian matrix is calculated many times
	Eigen::Matrix<double,2,2> jacobianMatrix;

	for (int k=0; k<mQuadraturePoints.cols(); ++k)
	{

		xi  = mQuadraturePoints(0,k);
		eta = mQuadraturePoints(1,k);
		weight = mQuadratureWeights(0,k);

		Bi = calcBi(i, xi, eta, nodesCoordGlob);
		Bj = calcBi(j, xi, eta, nodesCoordGlob);

		jacobianMatrix = calcJacobianMatrix(xi, eta, nodesCoordGlob);

		locK += weight*Bi.transpose()*mD*Bj*jacobianMatrix.determinant();
	}

	return locK;
}

/*!
 * The vector of strains \f$ \bm{\varepsilon} = \left\{ \varepsilon_{11}, \varepsilon_{22}, 2 \varepsilon_{12} \right\}^T \f$
 * \f[
 * 		\bm{\varepsilon} \left( \bm{d} \right) = \bm{B} \bm{d}
 * \f]
 * where \f$ \bm{d} = \left\{ u_{1x}, u_{1y}, u_{2x}, u_{2y}, u_{3x}, u_{3y}, u_{4x}, u_{4y}  \right\}^T \f$
 * \f[
  		\bm{\sigma} \left( \bm{d} \right) = \bm{D} \epsilon \left( \bm{d} \right)
 * \f]
 */
void P4::update(Eigen::Matrix<double,2,4> &nodesCoordGlob, Eigen::VectorXd &u, Eigen::MatrixXd &sigma, Eigen::VectorXd &force) const
{
	// May be simplify to set mQuadPoints only in one direction
	sigma = Eigen::MatrixXd::Zero(3,mQuadraturePoints.cols());

	for (int k=0; k<mQuadraturePoints.cols(); ++k)
	{
		double xi  = mQuadraturePoints(0,k);
		double eta = mQuadraturePoints(1,k);

		Eigen::Matrix<double,3,8> B = calcB(xi, eta, nodesCoordGlob);
		sigma.col(k) = mD*B*u;

		double weight = mQuadratureWeights(0,k);
		Eigen::Matrix<double,2,2> jacobianMatrix = calcJacobianMatrix(xi, eta, nodesCoordGlob);

		force -= weight*B.transpose()*sigma.col(k)*jacobianMatrix.determinant();
	}

}

void P4::updateForce(Eigen::Matrix<double,2,4> &nodesCoordGlob, Eigen::VectorXd &u, Eigen::VectorXd &force) const
{
	for (int k=0; k<mQuadraturePoints.cols(); ++k)
	{
		double xi  = mQuadraturePoints(0,k);
		double eta = mQuadraturePoints(1,k);

		Eigen::Matrix<double,3,8> B = calcB(xi, eta, nodesCoordGlob);

		Eigen::VectorXd sigma = mD*B*u;


		double weight = mQuadratureWeights(0,k);
		Eigen::Matrix<double,2,2> jacobianMatrix = calcJacobianMatrix(xi, eta, nodesCoordGlob);

		force -= weight*B.transpose()*sigma*jacobianMatrix.determinant();
	}
}






Eigen::Matrix<double,3,3> P4::getD() const
{
	return mD;
}



