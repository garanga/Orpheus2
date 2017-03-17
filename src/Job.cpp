/*
 * Job.cpp
 *
 *  Created on: Jan 31, 2017
 *      Author: pavel
 */
#include "Job.hpp"

#include "includes.hpp"

#include <iostream>
#include <string>
#include <vector>


Job::Job(std::string name, Model* model)
	: name_(name), model_(model)
{
	std::cout << "The job '" << name_ << "' is created" << std::endl;
	odb_ = new Odb(model);
}


Job::~Job()
{
	delete odb_;
}


std::string
Job::getName() const
{
	return name_;
}


Odb*
Job::Submit()
{
	LinearSolver* solver = new LinearSolver(model_, odb_);

	solver->solve();

	std::cout << odb_->getStepNum() << std::endl;

	for (int i = 0; i < odb_->getStepNum(); ++i)
	{
		std::cout << odb_->getStep(i)->getFrameNum() << std::endl;
	}

	{
	std::vector<double> data = odb_->getStep(0)->getFrame(0)->getFieldOutput(0)->getData();

	for (auto it = data.begin(); it < data.end(); ++it)
	{
		std::cout << *it << "  ";
	}
	std::cout << std::endl;
	}

	{
	std::vector<double> data = odb_->getStep(0)->getFrame(1)->getFieldOutput(0)->getData();

	for (auto it = data.begin(); it < data.end(); ++it)
	{
		std::cout << *it << "  ";
	}
	std::cout << std::endl;
	}

	saveOdb(odb_,"example");

	return odb_;
}


/*!
 * \f[
 *		\bm{R}^k = \iint\limits_{\Gamma_t} \bm{N}^T \bm{t} \, d\Gamma +
 *		           \iiint\limits_{\Omega} \bm{N}^T \bm{f}^b \, d\Omega -
 *		           \iiint\limits_{\Omega} \bm{B}^T \bm{\sigma}^k \, d\Omega
 * \f]
 */
//void updateSR(Eigen::MatrixXd &sigma, Eigen::VectorXd &force)
//{

//	Mesh* mesh = mModel->parts[0]->mesh;
//
//	Eigen::VectorXd sigma(3,3*3*mesh->elementsNum)
//
//	int* connectivity;
//	double* coord;
//	Eigen::Matrix<double,2,4> nodesCoordGlob; //!!!!!!!
//
//	int k,l;

//	Eigen::Matrix<double,2,2> locK;

//	for (auto it = (mesh->elements).begin(); it < (mesh->elements).end(); it++)
//	{
//		connectivity = (*it)->connectivity();
//
//		for (int i=0; i<(*it)->type()->getNumNodes(); i++)
//		{
//			coord = (mesh->nodes)[connectivity[i]]->getCoord();
//			nodesCoordGlob(0,i) = coord[0];
//			nodesCoordGlob(1,i) = coord[1];
//		}
//
//		for (int i=0; i<(*it)->type()->getNumNodes(); i++)
//		{
//			k = connectivity[i];
//			for (int j=0; j<(*it)->type()->getNumNodes(); j++)
//			{
//				l = connectivity[j];
//				locK = (*it)->type()->calcLocK(i, j, nodesCoordGlob);
//
//				triplets.push_back(T(2*k  ,2*l  ,locK(0,0)));
//				triplets.push_back(T(2*k  ,2*l+1,locK(0,1)));
//				triplets.push_back(T(2*k+1,2*l  ,locK(1,0)));
//				triplets.push_back(T(2*k+1,2*l+1,locK(1,1)));
//
//			}
//		}
//	}
//
//	Eigen::SparseMatrix<double> globK(2*mesh->nodesNum,2*mesh->nodesNum);
//
//	globK.setFromTriplets(triplets.begin(), triplets.end());
//
//	return globK;
//}



//Eigen::SparseMatrix<double> Job::calcGlobK()
//{
//
//	Mesh* mesh = mModel->parts[0]->mesh;
//
//	std::cout.precision(2);
//	std::cout.setf(std::ios::showpos);
//	std::cout.setf(std::ios::scientific);
//
//	std::vector<T> triplets;
//
//	int* connectivity;
//	double* coord;
//	Eigen::Matrix<double,2,4> nodesCoordGlob;
//
//	int k,l;
//
//	Eigen::Matrix<double,2,2> locK;
//
//	for (auto it = (mesh->elements).begin(); it < (mesh->elements).end(); it++)
//	{
//		connectivity = (*it)->connectivity();
//
//		for (int i=0; i<(*it)->type()->getNumNodes(); i++)
//		{
//			coord = (mesh->nodes)[connectivity[i]]->getCoord();
//			nodesCoordGlob(0,i) = coord[0];
//			nodesCoordGlob(1,i) = coord[1];
//		}
//
//		for (int i=0; i<(*it)->type()->getNumNodes(); i++)
//		{
//			k = connectivity[i];
//			for (int j=0; j<(*it)->type()->getNumNodes(); j++)
//			{
//				l = connectivity[j];
//				locK = (*it)->type()->calcLocK(i, j, nodesCoordGlob);
//
//				triplets.push_back(T(2*k  ,2*l  ,locK(0,0)));
//				triplets.push_back(T(2*k  ,2*l+1,locK(0,1)));
//				triplets.push_back(T(2*k+1,2*l  ,locK(1,0)));
//				triplets.push_back(T(2*k+1,2*l+1,locK(1,1)));
//
//			}
//		}
//	}
//
//	Eigen::SparseMatrix<double> globK(2*mesh->nodesNum,2*mesh->nodesNum);
//
//	globK.setFromTriplets(triplets.begin(), triplets.end());
//
////	globK.makeCompressed();
//
////	std::cout << globK << std::endl;
//
//	return globK;
//}


//
//    std::vector< Eigen::Triplet<int> > triplets;
//
//    triplets.push_back(Eigen::Triplet<int>(0, 1, 3));
//    triplets.push_back(Eigen::Triplet<int>(1, 0, 22));
//    triplets.push_back(Eigen::Triplet<int>(2, 1, 5));
//    triplets.push_back(Eigen::Triplet<int>(2, 3, 1));
//    triplets.push_back(Eigen::Triplet<int>(4, 2, 14));
//    triplets.push_back(Eigen::Triplet<int>(4, 4, 8));
//
//    A.setFromTriplets(triplets.begin(), triplets.end());
//
////    A.insert(0, 0);
//    std::cout << A;
//
////    A.makeCompressed();
//
////    std::cout << std::endl << A;
//
//
//
