/*
 * trialspacefactory.cpp
 *
 *  Created on: 17.02.2010
 *      Author: andreasvogel
 */

#include "local_shape_function_set_provider.h"

namespace ug{

struct UG_ERROR_CannotRegisterLocalShapeFunctionSet{};

LocalShapeFunctionSetProvider::
LocalShapeFunctionSetProvider()
{
	if(!init_standard_local_shape_function_sets<ReferenceEdge>())
		throw(UG_ERROR_CannotRegisterLocalShapeFunctionSet());
	if(!init_standard_local_shape_function_sets<ReferenceTriangle>())
		throw(UG_ERROR_CannotRegisterLocalShapeFunctionSet());
	if(!init_standard_local_shape_function_sets<ReferenceQuadrilateral>())
		throw(UG_ERROR_CannotRegisterLocalShapeFunctionSet());
	if(!init_standard_local_shape_function_sets<ReferenceTetrahedron>())
		throw(UG_ERROR_CannotRegisterLocalShapeFunctionSet());
	if(!init_standard_local_shape_function_sets<ReferencePyramid>())
		throw(UG_ERROR_CannotRegisterLocalShapeFunctionSet());
	if(!init_standard_local_shape_function_sets<ReferencePrism>())
		throw(UG_ERROR_CannotRegisterLocalShapeFunctionSet());
	if(!init_standard_local_shape_function_sets<ReferenceHexahedron>())
		throw(UG_ERROR_CannotRegisterLocalShapeFunctionSet());
};

} // namespace ug

