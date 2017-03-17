#include "includes.hpp"

#include <ctime>
#include <iostream>
#include <vector>

#include <Eigen/Dense>

int main(int argc, char* argv[])
{
	// Model database
	Mdb* mdb;
	mdb = new Mdb();

	// Model
	Model* model;
	model = mdb->createModel("Model-1");

	// Creating a solid body
	Body* body;
	body = model->createBody("Body-1");

	body->addPoint(-5.0,-5.0);	// 0
	body->addPoint( 5.0,-5.0);	// 1
	body->addPoint( 5.0, 5.0);	// 2
	body->addPoint(-5.0, 5.0);	// 3

	body->addLine(0,1);			// 0
	body->addLine(1,2);			// 1
	body->addLine(2,3);			// 2
	body->addLine(3,0);			// 3

	for (int i = 0; i < body->getLinesNum(); ++i)
	{
		body->setLineDivisions(i, 10);
	}

	// Creating a solid geometrical model
	Part* part;
	part = model->createPart("Part-1");
	double* sizes = new double[2];
	sizes[0] = 10.0; sizes[1] = 10.0;
	part->setSizes(sizes);
	int* divisions = new int[2];
	divisions[0] =  1; divisions[1] =  1;
	part->setDivisions(divisions);

	Material* material;
	material = model->createIsotropic("Material-1", 200.0e9, 0.25);

	// Creating a discrete geometrical model

	P4* type;
	type = new P4(material);

	////////////////////////////////////
	Mesh* mesh1;
	mesh1 = body->createMesh(type);
	////////////////////////////////////

	Mesh* mesh;
	mesh = part->CreateMesh(type);

	// Creating a step for boundary conditions application

	StaticStep* step1;
//	step1 = dynamic_cast<StaticStep*>(model->CreateStaticStep("Step-1", 0.0, 0.5, 0.1, 0.0, 0.5));
	step1 = model->createStaticStep("Step-1", 0.0, 0.5, 0.1, 0.0, 0.5);

	step1->addFieldOutputRequest(FieldType::U);
	{
		std::vector<int> region;
		double* value;
		region.push_back(2  -1);
		value = new double[2];
		value[0] = 0.0; value[1] = 0.0;

		step1->CreateDispacementConstraint("Constraint-1",
			                          	   region,
										   DisplacementConstraintType::UXY,
										   value);
	}

	{
		std::vector<int> region;
		double* value;

		region.push_back(2*2-1);
		value = new double[2];
		value[0] = 0.0;
		step1->CreateDispacementConstraint("Constraint-2",
				                           region,
										   DisplacementConstraintType::UX,
										   value);
	}

	{
		std::vector<int> region;
		double* value;

		region.push_back(0*2);
		value = new double[2];
		value[1] = -200.0e6;
		step1->CreateConcentratedLoad("Load-1",
				                      region,
									  ConcentratedLoadType::FY,
									  value);
	}

	{
		std::vector<int> region;
		double* value;

		region.push_back(1*2);
		value = new double[2];
		value[1] = 200.0e6;
		step1->CreateConcentratedLoad("Load-2",
				                      region,
									  ConcentratedLoadType::FY,
									  value);
	}



	StaticStep* step2;
//	step2 = dynamic_cast<StaticStep*>(model->CreateStaticStep("Step-2", 0.5, 1.0, 0.1, 0.5, 1.0));
	step2 = model->createStaticStep("Step-2", 0.5, 1.0, 0.1, 0.5, 1.0);

	step2->addFieldOutputRequest(FieldType::U);

	{
		std::vector<int> region;
		double* value;
		region.push_back(2  -1);
		value = new double[2];
		value[0] = 0.0; value[1] = 0.0;

		step2->CreateDispacementConstraint("Constraint-1",
			                          	   region,
										   DisplacementConstraintType::UXY,
										   value);
	}

	{
		std::vector<int> region;
		double* value;

		region.push_back(2*2-1);
		value = new double[2];
		value[0] = 0.0;
		step2->CreateDispacementConstraint("Constraint-2",
				                           region,
										   DisplacementConstraintType::UX,
										   value);
	}

	{
		std::vector<int> region;
		double* value;

		region.push_back(0*2);
		value = new double[2];
		value[1] = -200.0e6;
		step2->CreateConcentratedLoad("Load-1",
				                      region,
									  ConcentratedLoadType::FY,
									  value);
	}

	{
		std::vector<int> region;
		double* value;

		region.push_back(1*2);
		value = new double[2];
		value[1] = 200.0e6;
		step2->CreateConcentratedLoad("Load-2",
				                      region,
									  ConcentratedLoadType::FY,
									  value);
	}















	Job* job;
	job = mdb->createJob("Job-1", model);

	time_t time1;
	time(&time1);





	job->Submit();





	time_t time2;
	time(&time2);

	std::cout << time2-time1 << "s" << std::endl;

	std::cout << mdb->getModel(0) << std::endl;
	std::cout << mdb->getModel("Model-1") << std::endl;



//	delete job;
//	delete step1;
//	delete step2;
//	delete mesh;
//	delete material;
//	delete part;
//	delete model;
	delete mdb;

	return 0;




//	Eigen::Matrix<double,2,4> xy;
//	xy << -1.0,  1.0, 1.0, -1.0,
//		  -1.0, -1.0, 1.0,  1.0;
//
//	Eigen::Matrix<double,2,2> ans;
//
//	ans = type->calcLocK(0,0,xy);
//	ans = type->calcLocK(0,1,xy);
//	ans = type->calcLocK(0,2,xy);
//	ans = type->calcLocK(0,3,xy);
//	ans = type->calcLocK(1,1,xy);
//	ans = type->calcLocK(1,2,xy);
//	ans = type->calcLocK(1,3,xy);
//	ans = type->calcLocK(2,2,xy);
//	ans = type->calcLocK(2,3,xy);
//	ans = type->calcLocK(3,3,xy);


























//	IsotropicMaterial* material;
//	material = model -> CreateIsotropicMaterial("Material-1", 200.0e9, 0.3);
//
//
//	cout << material -> name << endl;

//	int numberOfPoints = 4, numberOfSegments = 4;
//	double* points;
//	points = new double[2*numberOfPoints];
//	points[0] = 0.0; points[1] = 0.0;
//	points[2] = 1.0; points[3] = 0.0;
//	points[4] = 1.0; points[5] = 1.0;
//	points[6] = 0.0; points[7] = 1.0;
//	int* segments;
//	segments = new int[2*numberOfSegments];
//	segments[0] = 1; segments[1] = 2;
//	segments[2] = 2; segments[3] = 3;
//	segments[4] = 3; segments[5] = 4;
//	segments[6] = 4; segments[7] = 1;
//
//	Part part(numberOfPoints, numberOfSegments, points, segments);
//
//	std::cout << part.GetNumberOfPoints() << std::endl;
//
//	Mesh mesh(part);
//
//    Eigen::SparseMatrix<int> A(5,5);
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
//	delete[] points;
//	delete[] segments;


}
