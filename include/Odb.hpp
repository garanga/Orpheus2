/*
 * Odb.hpp
 *
 *  Created on: Jan 31, 2017
 *      Author: pavel
 */

#ifndef ODB_HPP_
#define ODB_HPP_

#include <map>
#include <string>
#include <vector>

class Model;
class Part;		// OdbPart in Abaqus
class Body;
class Step;
class OdbStep;

enum class Constant;
enum class FieldType;

class Odb
{

public:

	Odb(Model*);
   ~Odb();

    OdbStep* createStep(Step*);

    int
	getStepNum() const;

    OdbStep*
	getStep(int) const;

    OdbStep*
	getStep(Constant) const;

    friend void
	saveOdb(const Odb*, std::string);

private:

   std::vector<Part*> parts_;
   std::vector<Body*> bodies_;



   std::vector<OdbStep*> steps_;



//	std::vector<Node*>       nodes_;
//	std::vector<Element*>    elements_;
//	std::vector<StaticStep*> steps_;

};


#endif /* ODB_HPP_ */
