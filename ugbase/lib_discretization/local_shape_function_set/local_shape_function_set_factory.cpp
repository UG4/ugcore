/*
 * trialspacefactory.cpp
 *
 *  Created on: 17.02.2010
 *      Author: andreasvogel
 */

#include "local_shape_function_set_factory.h"

namespace ug{

LocalShapeFunctionSetFactory::
LocalShapeFunctionSetFactory()
{
	if(init_standard_local_shape_function_sets<ReferenceEdge>() != true) assert(0);
	if(init_standard_local_shape_function_sets<ReferenceTriangle>() != true) assert(0);
	if(init_standard_local_shape_function_sets<ReferenceQuadrilateral>() != true) assert(0);
	if(init_standard_local_shape_function_sets<ReferenceTetrahedron>() != true) assert(0);
	if(init_standard_local_shape_function_sets<ReferencePrism>() != true) assert(0);
	if(init_standard_local_shape_function_sets<ReferenceHexahedron>() != true) assert(0);
};

LocalShapeFunctionSetFactory&
LocalShapeFunctionSetFactory::
inst()
{
	static LocalShapeFunctionSetFactory myInst;
	return myInst;
};

} // namespace ug

