/*
 * LinearSolver.cpp
 *
 *  Created on: Jan 27, 2017
 *      Author: pavel
 */

#include "LinearSolver.hpp"

#include "includes.hpp"

using namespace Eigen;

typedef Triplet<double> T;

LinearSolver::LinearSolver(Model* model, Odb* odb)
	: model_(model), odb_(odb)
{

}


void
LinearSolver::solve()
{
	////////////////
	// Input data //
	////////////////

	// For now it is only one part exists
	Mesh* mesh = model_->parts[0]->mesh;

	VectorXd  u    = VectorXd::Zero(2*mesh->nodesNum);
//	VectorXd du    = VectorXd::Zero(2*mesh->nodesNum);
	VectorXd force = VectorXd::Zero(2*mesh->nodesNum);

	// Initializing the first step
	auto step = (model_->steps).begin();
	// Creating the step in the output database
	odb_->createStep(*step);


	std::cout << "//////////////////////////////////////////////////" << std::endl;
	std::cout << "// The step " << (*step)->getName() << " is initialized" << std::endl;
	std::cout << "//////////////////////////////////////////////////" << std::endl;

	double timeBegin = (*step)->timeBegin;					// Starting time
	std::cout << "Starting time: " << timeBegin << std::endl;
	double timeEnd   = (*step)->timeEnd;				// Ending time
	std::cout << "Ending time: " << timeEnd << std::endl;
	double stepTime = timeEnd-timeBegin;					// Time interval for the load step

	double timeIncrement   = (*step)->timeIncrement;		// Time increment
	std::cout << "Time increment: " << timeIncrement << std::endl;
//	double timeIncrementOld = timeIncrement;					// Time increment of the previous iteration
//	std::cout << "Time increment old: " << timeIncrementOld << std::endl;

	double loadFactorBegin   = (*step)->loadFactorBegin;	// Starting load factor
	std::cout << "Starting load factor: " << loadFactorBegin << std::endl;
	double loadFactorEnd   = (*step)->loadFactorEnd;		// Ending load factor
	std::cout << "Ending load factor: " << loadFactorEnd << std::endl;




	double time = timeBegin;								// Current time

	/////////////////////////
	// Load increment loop //
	/////////////////////////

//	int istep = -1;											// Load increment number

	int flag10 = 1; // Open the door to the load increment loop
	while (flag10)
	{
	    flag10 = 0; 									// Close the door to the load increment loop

//		VectorXd uOld = u;								// Saving the displacement vector of the previous iteration

//		double timeOld = time;							// Saving the current time

//		updateOutputs(u, force);

		time += timeIncrement;												// Increase time
//		istep += 1;															// Increase load increment number

		if ((time-timeEnd)>1.0e-10) /* Using chosen time increment we are jump over the ending time of the load step */
		{
			if ((timeIncrement-(time-timeEnd))>1.0e-10)	/* variant1: on the previous iteration we were far from the ending time */
			{
				timeIncrement -= time-timeEnd;			// Choosing correct time increment
				time = timeEnd;							//
			}
			else										/* variant2: on the previous iteration we are near with the ending time */
			{
				++step;									// Progressing to the next load step

				if (step==(model_->steps).end())		/* The previous load step was last load step */
				{
					flag10 = 0;							// Close the door to the load increment loop													??????? It is already closed
					break;								// Leave the bisection loop. Since the door to the load increment loop is closed this will lead to ending of the load loop also
				}
				else													// Next load step
				{
					odb_->createStep(*step);

					std::cout << "//////////////////////////////////////////////////" << std::endl;
					std::cout << "// The step " << (*step)->getName() << " is initialized" << std::endl;
					std::cout << "//////////////////////////////////////////////////" << std::endl;

					time -= timeIncrement;

					timeBegin = (*step)->timeBegin;					// Starting time
					timeEnd = (*step)->timeEnd;						// Ending time
					stepTime = timeEnd-timeBegin;					// Time interval for the load step

					timeIncrement   = (*step)->timeIncrement;		// Time increment
//					timeIncrementOld = timeIncrement;

					loadFactorBegin   = (*step)->loadFactorBegin;	// Starting load factor
					loadFactorEnd   = (*step)->loadFactorEnd;		// Starting load factor

					time += timeIncrement;							// Current time
				}
			}
		}

		// Load factor and (displacements variation) factor

		double loadFactor = loadFactorBegin+(time-timeBegin)/stepTime*(loadFactorEnd-loadFactorBegin);
		std::cout << "load factor: " << loadFactor << std::endl;
		double dispFactor = timeIncrement/stepTime*(loadFactorEnd-loadFactorBegin);
		std::cout << "displacements variation factor: " << dispFactor << std::endl;

		VectorXd du = VectorXd::Zero(2*mesh->nodesNum);

		// Initialize global stiffness K and residual vector F

		SparseMatrix<double> globK(2*mesh->nodesNum,2*mesh->nodesNum);
		globK = calcGlobK();

		VectorXd force = VectorXd::Zero(2*mesh->nodesNum);
		updateForce(u, force);

		// Prescribed loads

		for (auto load = ((*step)->loads).begin(); load < ((*step)->loads).end(); ++load)
		{
			std::vector<int> indices;
			std::vector<double> values;

			if (static_cast<int>((*load)->type) & static_cast<int>(ConcentratedLoadType::FX))
			{
				std::vector<int> region = ((*load)->region);
				for (unsigned i=0; i<region.size(); ++i)
				{
					region[i] *= 2;
					region[i] += 0;
				}
				indices.insert(indices.end(),region.begin(),region.end());

				std::vector<double> value(region.size(),((*load)->value[0])*loadFactor);
				values.insert(values.end(),value.begin(),value.end());
			}

			if (static_cast<int>((*load)->type) & static_cast<int>(ConcentratedLoadType::FY))
			{
				std::vector<int> region = ((*load)->region);
				for (unsigned i=0; i<region.size(); ++i)
				{
					region[i] *= 2;
					region[i] += 1;
				}
				indices.insert(indices.end(),region.begin(),region.end());

				std::vector<double> value(region.size(),((*load)->value[1])*loadFactor);
				values.insert(values.end(),value.begin(),value.end());
			}

			for (unsigned k=0; k<force.size(); ++k)
			{
				for (unsigned i=0; i<indices.size(); ++i)
				{
					if (k==indices[i])
					{
						force(k) += values[i];
					}
				}
			}
		}

		// Prescribed constraints

		for (auto constraint = ((*step)->constraints).begin(); constraint < ((*step)->constraints).end(); ++constraint)
		{
			std::vector<int> indices;
			std::vector<double> values;

			if (static_cast<int>((*constraint)->type) & static_cast<int>(DisplacementConstraintType::UX))
			{
				std::vector<int> region = ((*constraint)->region);
				for (unsigned int i=0; i<region.size(); ++i)
				{
					region[i] *= 2;
					region[i] += 0;
				}
				indices.insert(indices.end(),region.begin(),region.end());
				std::vector<double> value(region.size(),((*constraint)->value[0])*dispFactor);
				values.insert(values.end(),value.begin(),value.end());
			}

			if (static_cast<int>((*constraint)->type) & static_cast<int>(DisplacementConstraintType::UY))
			{
				std::vector<int> region = ((*constraint)->region);
				for (unsigned int i=0; i<region.size(); ++i)
				{
					region[i] *= 2;
					region[i] += 1;
				}
				indices.insert(indices.end(),region.begin(),region.end());
				std::vector<double> value(region.size(),((*constraint)->value[1])*dispFactor);
				values.insert(values.end(),value.begin(),value.end());
			}

			for (unsigned k=0; k<globK.outerSize(); ++k)
			{
				for (Eigen::SparseMatrix<double>::InnerIterator it(globK,k); it; ++it)
				{
					for (unsigned i=0; i<indices.size(); ++i)
					{
						if (it.row() == indices[i])
						{
							it.valueRef() = it.row() == it.col() ? model_->materials.back()->getYoung() : 0.0;
						}
					}
				}

			}

			for (unsigned i=0; i<indices.size(); ++i)
			{
				force(indices[i]) = values[i]*model_->materials.back()->getYoung();
			}

		}

		SparseLU<SparseMatrix<double>> linearEquationsSolver(globK);
		VectorXd sol = linearEquationsSolver.solve(force);

		du = sol;
		u += du;

//		std::cout << u << std::endl;
//		std::cout << std::endl;

		OdbFrame* frame = odb_->getStep(Constant::LAST)->createFrame();

		FieldType field = FieldType::U;
		std::vector<double> data(u.data(), u.data() + u.rows());

		FieldOutput* fieldOutput = new FieldOutput(FieldType::U, data);

		frame->addFieldOutput(fieldOutput);

		flag10=1;

	}

}


SparseMatrix<double> LinearSolver::calcGlobK() const
{

	Mesh* mesh = model_->parts[0]->mesh;

	std::cout.precision(2);
	std::cout.setf(std::ios::showpos);
	std::cout.setf(std::ios::scientific);

	std::vector<T> triplets;

	int* connect;
	double* coord;
	Matrix<double,2,4> nodesCoordGlob;

	int k,l;

	Matrix<double,2,2> locK;

	for (auto it = (mesh->elements).begin(); it < (mesh->elements).end(); it++)
	{
		connect = (*it)->getConnect();

		for (int i=0; i<(*it)->type()->getNodesNum(); i++)
		{
			coord = (mesh->nodes)[connect[i]]->getCoord();
			nodesCoordGlob(0,i) = coord[0];
			nodesCoordGlob(1,i) = coord[1];
		}

		for (int i=0; i<(*it)->type()->getNodesNum(); i++)
		{
			k = connect[i];
			for (int j=0; j<(*it)->type()->getNodesNum(); j++)
			{
				l = connect[j];
				locK = (*it)->type()->calcLocK(i, j, nodesCoordGlob);

				triplets.push_back(T(2*k  ,2*l  ,locK(0,0)));
				triplets.push_back(T(2*k  ,2*l+1,locK(0,1)));
				triplets.push_back(T(2*k+1,2*l  ,locK(1,0)));
				triplets.push_back(T(2*k+1,2*l+1,locK(1,1)));

			}
		}
	}

	SparseMatrix<double> globK(2*mesh->nodesNum,2*mesh->nodesNum);

	globK.setFromTriplets(triplets.begin(), triplets.end());

//	globK.makeCompressed();

//	std::cout << globK << std::endl;

	return globK;
}


void LinearSolver::updateForce(VectorXd &u, VectorXd &force) const
{
	Mesh* mesh = model_->parts[0]->mesh;

	for (auto it = (mesh->elements).begin(); it < (mesh->elements).end(); it++)
	{
		int* connect = (*it)->getConnect();

		Matrix<double,2,4> nodesCoordGlob;
		VectorXd uLoc(8);
		VectorXd forceLoc(8);
		for (int i=0; i<(*it)->type()->getNodesNum(); i++)
		{
			double* coord = (mesh->nodes)[connect[i]]->getCoord();
			nodesCoordGlob(0,i) = coord[0];
			nodesCoordGlob(1,i) = coord[1];

			uLoc(2*i  ) = u(2*connect[i]  );
			uLoc(2*i+1) = u(2*connect[i]+1);

			forceLoc(2*i  ) = force(2*connect[i]  );
			forceLoc(2*i+1) = force(2*connect[i]+1);
		}

		(*it)->type()->updateForce(nodesCoordGlob, uLoc, forceLoc);

		// Updating the global force
		for (int i=0; i<(*it)->type()->getNodesNum(); i++)
		{
			force(2*connect[i]  ) = forceLoc(2*i  );
			force(2*connect[i]+1) = forceLoc(2*i+1);
		}

	}

}







void LinearSolver::update(VectorXd &u, MatrixXd &sigma, VectorXd &force) const
{
	Mesh* mesh = model_->parts[0]->mesh;

	for (auto it = (mesh->elements).begin(); it < (mesh->elements).end(); it++)
	{
		int* connect = (*it)->getConnect();

		Matrix<double,2,4> nodesCoordGlob;
		VectorXd uLoc(8);
		VectorXd forceLoc(8);
		for (int i=0; i<(*it)->type()->getNodesNum(); i++)
		{
			double* coord = (mesh->nodes)[connect[i]]->getCoord();
			nodesCoordGlob(0,i) = coord[0];
			nodesCoordGlob(1,i) = coord[1];

			uLoc(2*i  ) = u(2*connect[i]  );
			uLoc(2*i+1) = u(2*connect[i]+1);

			forceLoc(2*i  ) = force(2*connect[i]  );
			forceLoc(2*i+1) = force(2*connect[i]+1);
		}

		Eigen::MatrixXd sigmaLoc(3,8);

		(*it)->type()->update(nodesCoordGlob, uLoc, sigmaLoc, forceLoc);

		// Updating the global force
		for (int i=0; i<(*it)->type()->getNodesNum(); i++)
		{
			force(2*connect[i]  ) = forceLoc(2*i  );
			force(2*connect[i]+1) = forceLoc(2*i+1);
		}

	}

}
