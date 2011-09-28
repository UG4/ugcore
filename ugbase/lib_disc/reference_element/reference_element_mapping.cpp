/*
 * reference_element_mapping.cpp
 *
 *  Created on: 21.07.2011
 *      Author: andreasvogel
 */

#include "common/util/provider.h"

#include "reference_element_mapping.h"

#include "reference_edge.h"
#include "reference_triangle.h"
#include "reference_quadrilateral.h"
#include "reference_tetrahedron.h"
#include "reference_prism.h"
#include "reference_pyramid.h"
#include "reference_hexahedron.h"

namespace ug{

ReferenceMappingProvider::
ReferenceMappingProvider()
{
//	clear mappings
	for(int d = 0; d < 4; ++d)
		for(int rd = 0; rd < 4; ++rd)
			for(int roid = 0; roid < NUM_REFERENCE_OBJECTS; ++roid)
				m_vvvMapping[d][rd][roid] = NULL;

//	set mappings

//	edge
	set_mapping<1,1>(ROID_EDGE, Provider<DimReferenceMappingWrapper<ReferenceMapping<ReferenceEdge, 1> > >::get());
	set_mapping<1,2>(ROID_EDGE, Provider<DimReferenceMappingWrapper<ReferenceMapping<ReferenceEdge, 2> > >::get());
	set_mapping<1,3>(ROID_EDGE, Provider<DimReferenceMappingWrapper<ReferenceMapping<ReferenceEdge, 3> > >::get());

//	triangle
	set_mapping<2,2>(ROID_TRIANGLE, Provider<DimReferenceMappingWrapper<ReferenceMapping<ReferenceTriangle, 2> > >::get());
	set_mapping<2,3>(ROID_TRIANGLE, Provider<DimReferenceMappingWrapper<ReferenceMapping<ReferenceTriangle, 3> > >::get());

//	quadrilateral
	set_mapping<2,2>(ROID_QUADRILATERAL, Provider<DimReferenceMappingWrapper<ReferenceMapping<ReferenceQuadrilateral, 2> > >::get());
	set_mapping<2,3>(ROID_QUADRILATERAL, Provider<DimReferenceMappingWrapper<ReferenceMapping<ReferenceQuadrilateral, 3> > >::get());

//	3d elements
	set_mapping<3,3>(ROID_TETRAHEDRON, Provider<DimReferenceMappingWrapper<ReferenceMapping<ReferenceTetrahedron, 3> > >::get());
	set_mapping<3,3>(ROID_PRISM, Provider<DimReferenceMappingWrapper<ReferenceMapping<ReferencePrism, 3> > >::get());
	set_mapping<3,3>(ROID_PYRAMID, Provider<DimReferenceMappingWrapper<ReferenceMapping<ReferencePyramid, 3> > >::get());
	set_mapping<3,3>(ROID_HEXAHEDRON, Provider<DimReferenceMappingWrapper<ReferenceMapping<ReferenceHexahedron, 3> > >::get());
}


} // end namespace ug
