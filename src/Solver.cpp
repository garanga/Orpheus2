/*
 * Solver.cpp
 *
 *  Created on: Jan 27, 2017
 *      Author: pavel
 */

#include "Solver.hpp"

#include "includes.hpp"
#include <cassert>
#include <algorithm>

using namespace Eigen;

typedef Triplet<double> T;

Solver::Solver(Model* model, Odb* odb)
	: model_(model), odb_(odb)
{

}


void
Solver::solve()
{
	////////////////
	// Input data //
	////////////////

	// For now it is only one part exists
	Mesh* mesh = model_->parts[0]->mesh;

	VectorXd  u    = Eigen::VectorXd::Zero(2*mesh->nodesNum);
	VectorXd du    = Eigen::VectorXd::Zero(2*mesh->nodesNum);
	VectorXd force = Eigen::VectorXd::Zero(2*mesh->nodesNum);

	// Initialize the first step

	auto step = (model_->steps).begin();
	odb_->createStep(*step);


	std::cout << "//////////////////////////////////////////////////" << std::endl;
	std::cout << "// The step " << (*step)->getName() << " is initialized" << std::endl;
	std::cout << "//////////////////////////////////////////////////" << std::endl;

	double timeBegin = (*step)->timeBegin;					// Starting time
	std::cout << "Starting time: " << timeBegin << std::endl;
	double timeEnd   = (*step)->timeEnd;					// Ending time
	std::cout << "Ending time: " << timeEnd << std::endl;
	double stepTime = timeEnd-timeBegin;					// Time interval for the load step

	double timeIncrement   = (*step)->timeIncrement;		// Time increment
	std::cout << "Time increment: " << timeIncrement << std::endl;
	double timeIncrementOld = timeIncrement;					// Time increment of the previous iteration
//	std::cout << "Time increment old: " << timeIncrementOld << std::endl;

	double loadFactorBegin   = (*step)->loadFactorBegin;	// Starting load factor
	std::cout << "Starting load factor: " << loadFactorBegin << std::endl;
	double loadFactorEnd   = (*step)->loadFactorEnd;		// Ending load factor
	std::cout << "Ending load factor: " << loadFactorEnd << std::endl;

	int bisectionLevel = 1;									// Initial bisection level
	VectorXd tary = VectorXd::Zero(maxBisectionLevel);		// Time stamps for the bisection process


	double time = timeBegin;								// Current time

	/////////////////////////
	// Load increment loop //
	/////////////////////////

	int istep = -1;											// Load increment number

	int flag10 = 1; // Open the door to the load increment loop
	while (flag10)
	{
		    flag10 = 0; 									// Close the door to the load increment loop
		int flag11 = 1; 									// Open  the door to the bisection loop
		int flag20 = 1; 									// Open  the door to the convergence iterations loop

		VectorXd uOld = u;									// Saving the displacement vector of the previous iteration

		if (bisectionLevel==1) /* where were no bisections before */
		{
			timeIncrement = timeIncrementOld;								// Time increment equals to the time increment of the previous iteration
			tary(bisectionLevel) = time + timeIncrement;
		}
		else /* where were bisections before  */
		{
			bisectionLevel -= 1;											// Reduce the bisection level
			timeIncrement = tary(bisectionLevel) - tary(bisectionLevel+1);	// New time increment
			tary(bisectionLevel+1) = 0;										// Empty converged bisection level

			istep -= 1;														// Decrease load increment 														???????
		}

//		std::cout << "time: " << time << std::endl;
//		std::cout << "timeIncrement: " << timeIncrement << std::endl;
//		std::cout << "timeIncrementOld: " << timeIncrementOld << std::endl;

		double timeOld = time;								// Saving the current time

		//////////////////////////////
		// Update history variables //
		//////////////////////////////

//		updateForce(u, force);								//																								???????

		time += timeIncrement;												// Increase time
		istep += 1;															// Increase load increment number

		////////////////////
		// Bisection loop //
		////////////////////

		std::cout << "33333333333333333333333333333333333333333333333333333333333333333333333333" << std::endl;

		while (flag11 == 1)
		{
			flag11 = 0;										// Close the door to the bisection loop

			if ((time-timeEnd)>1.0e-10) /* Using chosen time increment we are jump over the ending time of the load step */
			{
				if ((timeIncrement-(time-timeEnd))>1.0e-10)	/* variant1: on the previous iteration we were far from the ending time */
				{
					timeIncrement -= time-timeEnd;			// Choosing correct time increment
					timeIncrementOld = timeIncrement;		// Saving time increment
					time = timeEnd;							//
				}
				else										/* variant2: on the previous iteration we are near with the ending time */
				{
					std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!33333333333333333333333333333333333333" << std::endl;

					++step;									// Progressing to the next load step

					std::cout << "4444444444444444444444444444444444444444444444444444444444444444444444444" << std::endl;

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
						timeIncrementOld = timeIncrement;

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

//			/////////////////////////////////
//			// Start convergence iteration //
//			/////////////////////////////////

			int iter = 0;
			du = Eigen::VectorXd::Zero(2*mesh->nodesNum);		// Doubled

			while (flag20==1)
			{
				flag20 = 0;
				iter += 1;

				// Check maximum iteration
				assert(iter<=itra);

				// Initialize global stiffness K and residual vector F

				SparseMatrix<double> globK(2*mesh->nodesNum,2*mesh->nodesNum);
				VectorXd force = VectorXd::Zero(2*mesh->nodesNum);

				// Assemble K and F (update)

				globK = calcGlobK();

//				std::cout << globK << std::endl;

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

				// All indexes of the constrained dofs
				std::vector<int> constraintsIndexes;

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
						force(indices[i]) = iter==1 ? values[i]*model_->materials.back()->getYoung() : 0.0;
						constraintsIndexes.push_back(indices[i]);
					}

				}

				// Check convergence
				if (iter>1)
				{
					// Converting Eigen::VectorXd to std::vector
					std::vector<double> forceVector(force.data(), force.data()+force.rows());

					// Deleting the components of the vector corresponding to the fixed nodes
					std::vector<double> forceVectorReduced(forceVector);

					// Sorting the indexes in a container of fixed dofs indexes
					std::sort(constraintsIndexes.begin(),constraintsIndexes.end());

					for (auto it=constraintsIndexes.end()-1; it>=constraintsIndexes.begin(); --it)
					{
						forceVectorReduced.erase(forceVectorReduced.begin()+(*it));
					}

					// Converting std::vector forceVectorReduced to Eigen::VectorXd forceVectorReducedE
					VectorXd forceVectorReducedE;
					forceVectorReducedE = VectorXd::Map(forceVectorReduced.data(),forceVectorReduced.size());

					double resn = forceVectorReducedE.lpNorm<Infinity>();

					std::cout << " Iteration: "      << iter
							  << " Error: "          << resn
							  << " Time: "           << time
							  << " Time increment: " << timeIncrement
							  << std::endl;

					if (resn<tol)
					{
						flag10=1;

						std::cout << "Convergence!!!" << std::endl;

						// The iteration is converged. So the results should be written to a frame

						OdbFrame* frame = odb_->getStep(Constant::LAST)->createFrame();


						FieldType field = FieldType::U;
						std::vector<double> data(u.data(), u.data() + u.rows());

						FieldOutput* fieldOutput = new FieldOutput(FieldType::U, data);

						frame->addFieldOutput(fieldOutput);












						break;
					}

					if ((resn>atol)||(iter>itra))	// Start bisection
					{
						bisectionLevel += 1;
						if (bisectionLevel < maxBisectionLevel)
						{
							timeIncrement *= 0.5;
							time = timeOld + timeIncrement;
							tary(bisectionLevel) = time;
							u = uOld;
							std::cout << "Not converged. Bisecting load increment " << bisectionLevel << std::endl;
						}
						else
						{
							std::cout << "Max No. of bisections" << std::endl;
							std::exit(1);
						}
						flag11 = 1;
						flag20 = 1;
						break;

					}
				}

//				std::cout << "WHIS IS THE FIRST ITERATION" << std::endl;

				// if (iter < 1) or ((iter>1)& rtol not big not small)

				// Solve the system of equations
				if (flag11==0)
				{

//					std::cout << globK << std::endl;
//					std::cout << std::endl;
//					std::cout << force << std::endl;
//					std::cout << std::endl;

					SparseLU<SparseMatrix<double>> linearEquationsSolver(globK);
					VectorXd soln = linearEquationsSolver.solve(force);

					du += soln;
					 u += soln;
					 flag20 = 1; // To the next iteration

//					 std::cout << soln << std::endl;
//					 std::cout << std::endl;
					 std::cout << u << std::endl;
					 std::cout << std::endl;

					 std::cout << "22222222222222222222222222222222222222222222222222222222222222222222222222222222222" << std::endl;

				}
				else
				{
					flag20 = 0;
				}

				if (flag10==1) {break;}

			}	// 20 Convergence iteration
		}  		// 11 Bisection
	}			// 10 Load increment
}





















//void Solver::solve()
//{
//	////////////////
//	// Input data //
//	////////////////
//
//	Eigen::MatrixXd sdispt(3,3);
//	sdispt << 1, 0, 0,
//	          1, 1, 0,
//			  3, 0, 0;
//
//	Eigen::MatrixXd extforce(2,3);
//	extforce << 0, 1, -200.0e6,
//	            2, 1,  200.0e6;
//
//
//	// For now it is only one part exists
//	Mesh* mesh = mModel->parts[0]->mesh;
//
//	Eigen::VectorXd  u = Eigen::VectorXd::Zero(2*mesh->nodesNum);
//	Eigen::VectorXd du = Eigen::VectorXd::Zero(2*mesh->nodesNum);
//
//	// Let's introduce the variable for residual force "force"
//	// The vector is initialized by zeroes before the main loop
//	Eigen::VectorXd force = Eigen::VectorXd::Zero(2*mesh->nodesNum);
//
//
//	// Let's introduce the variable for stresses "sigma"
//	// The stresses are stored such that in columns sigma11,sigma22,sigma12
//	//                                   in rows, for example, sigma11 for all integration
//	//                                   points and all elements bypass sequentially
//
//	Eigen::MatrixXd sigma(3,2*2*mesh->elementsNum);
//
//
//	Eigen::MatrixXd tims(5,2);
//
//	tims << 0.0, 0.5,
//			0.5, 1.0,
//			0.1, 0.1,
//			0.0, 0.5,
//			0.5, 1.0;
//
//	int nload = tims.cols(); // Total number of load steps
//	int iload = 0; // Initialize the first step
//
//	double timef = tims(0,iload); // Starting time
//	double timei = tims(1,iload); // Ending time
//	double delta = tims(2,iload); // Time increment
//	int cur1 = tims(3,iload);
//	int cur2 = tims(4,iload);
//
//	int delta0 = delta; // Initially saved time increment the same as desired time increment
//
//	int time = timef; // Starting time
//	double tdelta = timei-timef; // Time interval for load step
//
//	int itol = 1; // Initial bisection level
//
//	int ntol = 10; // Maximum bisection level
//	Eigen::VectorXd tary = Eigen::VectorXd::Zero(ntol); // Time stamps for bisections
//
//	/////////////////////////
//	// Load increment loop //
//	/////////////////////////
//
//	int istep = -1;
//
//	int flag10 = 1;
//	while (flag10) // Solution has been converged
//	{
//		flag10 = 0;
//		int flag11 = 1;
//		int flag20 = 1;
//
//		Eigen::VectorXd cu = u;
//
//		if (itol==1)
//		{
//			delta = delta0;
//			tary(itol) = time + delta;
//		}
//		else // Recover previous bisection
//		{
//			itol = itol - 1; // Reduce the bisection level
//			delta = tary(itol) - tary(itol+1); // New time increment
//			tary(itol+1) = 0; // Empty converged bisection level
//
//			istep = istep - 1; // Decrease load increment
//		}
//
////		std::cout << force << std::endl;
////		std::cout << std::endl;
//
//		///////////////////////////////////////////
//		// Update stresses and history variables //
//		///////////////////////////////////////////
//
//		// Update stresses and hystory variables
//		update(u, sigma, force);
//
////		std::cout << force << std::endl;
////		std::cout << std::endl;
//
//		time = time + delta;	// Increase time
//		istep += 1;
//
//		// Check time and control bisection
//		while (flag11 == 1) // Bisection loop start
//		{
//			flag11 = 0;
//			if ((time-timei)>1.0e-10)
//			{
//				if ((delta-(time-timei))>1.0e-10)
//				{
//					delta -= time-timei;
//					delta0 = delta;
//					time = timei;
//				}
//				else
//				{
//					iload += 1;				// Progress to the next load step
//					if (iload>nload)		// Finished final load step
//					{
//						flag10 = 0;			// Stop the program
//						break;
//					}
//					else					// Next load step
//					{
//						time -= delta;
//						delta = tims(2,iload);
//						delta0 = delta;
//						time = time + delta;
//						timef = tims(0,iload);
//						timei = tims(1,iload);
//						tdelta = timei-timef;
//						cur1 = tims(3,iload);
//						cur2 = tims(4,iload);
//					}
//
//				}
//			}
//
//			// Load factor and prescribed displacements
//
//			double factor = cur1+(time-timef)/tdelta*(cur2-cur1);
//			Eigen::VectorXd sdisp = delta*sdispt.col(2)/tdelta*(cur2-cur1);
//
//			/////////////////////////////////
//			// Start convergence iteration //
//			/////////////////////////////////
//
//			int iter = 0;
//			du = Eigen::VectorXd::Zero(2*mesh->nodesNum);
//
//			while (flag20==1)
//			{
//				flag20 = 0;
//				iter += 1;
//				// Check maximum iteration
//				assert(iter<=itra);
//
//				// Initialize global stiffness K and residual vector F
//
//				Eigen::SparseMatrix<double> globK(2*mesh->nodesNum,2*mesh->nodesNum);
//				Eigen::VectorXd force = Eigen::VectorXd::Zero(2*mesh->nodesNum);
//
//				// Assemble K and F (update)
//
//				globK = calcGlobK();
//				updateForce(u, force);
//
//				// Increase external force
//				Eigen::VectorXd locationF = 2*extforce.col(0) + extforce.col(1);
//
//
//				assert(extforce.rows>0);
//				for (int i=0; i<locationF.rows(); ++i)
//				{
//					force(locationF(i)) += factor*extforce(i,2);
//				}
//
//				// Prescribed displacement BC
//				int ndisp = sdispt.rows();
//				assert(ndisp>0);
//				Eigen::VectorXd locationU = 2*sdisp.col(0) + sdisp.col(1);
//
//				//		std::vector<int> indices;
//				//		std::vector<double> values;
//				//
//				//		if (static_cast<int>((*it)->mType) & static_cast<int>(DisplacementConstraintType::UX))
//				//		{
//				//			std::vector<int> region = ((*it)->mRegion);
//				//			for (unsigned int i=0; i<region.size(); ++i)
//				//			{
//				//				region[i] *= 2;
//				//				region[i] += 0;
//				//			}
//				//			indices.insert(indices.end(),region.begin(),region.end());
//				//			std::vector<double> value(region.size(),(*it)->mValue[0]);
//				//			values.insert(values.end(),value.begin(),value.end());
//				//		}
//				//
//				//		if (static_cast<int>((*it)->mType) & static_cast<int>(DisplacementConstraintType::UY))
//				//		{
//				//			std::vector<int> region = ((*it)->mRegion);
//				//			for (unsigned int i=0; i<region.size(); ++i)
//				//			{
//				//				region[i] *= 2;
//				//				region[i] += 1;
//				//			}
//				//			indices.insert(indices.end(),region.begin(),region.end());
//				//			std::vector<double> value(region.size(),(*it)->mValue[1]);
//				//			values.insert(values.end(),value.begin(),value.end());
//				//		}
//				//
//				//
//				//		for (int k=0; k<globK.outerSize(); ++k)
//				//		{
//				//			for (Eigen::SparseMatrix<double>::InnerIterator it1(globK,k); it1; ++it1)
//				//			{
//				//
//				//				for (unsigned int i=0; i<indices.size(); ++i)
//				//				{
//				//					if (it1.row() == indices[i] || it1.col() == indices[i])
//				//					{
//				//
//				////						it1.valueRef() = it1.row() == it1.col() ? it1.value() : 0.0;
//				//
//				//						if (it1.row() == it1.col())
//				//						{
//				//							it1.valueRef() = it1.value();
//				//						}
//				//						else
//				//						{
//				//							if (it1.col()==indices[i]) // this is the given column
//				//							{
//				//								globF(it1.row()) -= it1.value()*values[i];
//				//							}
//				//							it1.valueRef() = 0.0;
//				//						} /* it1.row() != it1.col() */
//				//
//				//					}
//				//				}
//				//
//				//
//				//
//				//			}
//				//		}
//
//
//
//
//
//
//
//
//
//			}
//
//
//
//
//
//
//		}
//
//
//
//
//	}
//
//
//}






void Solver::update(VectorXd &u, MatrixXd &sigma, VectorXd &force) const
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

SparseMatrix<double> Solver::calcGlobK() const
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


void Solver::updateForce(VectorXd &u, VectorXd &force) const
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

