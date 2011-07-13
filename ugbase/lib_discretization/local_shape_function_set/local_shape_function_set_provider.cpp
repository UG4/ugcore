/*
 * trialspacefactory.cpp
 *
 *  Created on: 17.02.2010
 *      Author: andreasvogel
 */

#include "local_shape_function_set_provider.h"

namespace ug{

LocalShapeFunctionSetProvider::
LocalShapeFunctionSetProvider()
{
	static bool init = false;

	if(!init)
	{
		init = true;
		if(!init_standard_local_shape_function_sets<ReferenceEdge>())
			throw(UGFatalError("Cannot register standard Edge trial spaces."));
		if(!init_standard_local_shape_function_sets<ReferenceTriangle>())
			throw(UGFatalError("Cannot register standard Triangle trial spaces."));
		if(!init_standard_local_shape_function_sets<ReferenceQuadrilateral>())
			throw(UGFatalError("Cannot register standard Quadrilateral trial spaces."));
		if(!init_standard_local_shape_function_sets<ReferenceTetrahedron>())
			throw(UGFatalError("Cannot register standard Tetrahedron trial spaces."));
		if(!init_standard_local_shape_function_sets<ReferencePyramid>())
			throw(UGFatalError("Cannot register standard Pyramid trial spaces."));
		if(!init_standard_local_shape_function_sets<ReferencePrism>())
			throw(UGFatalError("Cannot register standard Prism trial spaces."));
		if(!init_standard_local_shape_function_sets<ReferenceHexahedron>())
			throw(UGFatalError("Cannot register standard Hexahedron trial spaces."));
	}
};

std::map<LSFSID, std::vector<const LocalShapeFunctionSetBase*> >
LocalShapeFunctionSetProvider::m_baseMap =
		std::map<LSFSID, std::vector<const LocalShapeFunctionSetBase*> >();


} // namespace ug

