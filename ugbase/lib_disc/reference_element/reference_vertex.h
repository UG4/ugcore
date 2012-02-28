/*
 * reference_vertex.h
 *
 *  Created on: 15.06.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__REFERENCE_ELEMENT__REFERENCE_VERTEX__
#define __H__UG__LIB_DISC__REFERENCE_ELEMENT__REFERENCE_VERTEX__

#include "common/math/ugmath.h"
#include "lib_grid/grid/geometric_base_objects.h"
#include "reference_element_mapping.h"

namespace ug{

class ReferenceVertex : public DimReferenceElement<1>
{
	public:
	///	type of reference element
		static const ReferenceObjectID REFERENCE_OBJECT_ID = ROID_VERTEX;

	///	dimension of reference element
		static const int dim = 0;

	///	number of corners
		static const int num_corners = 1;

	///	number of eges
		static const int num_edges = 0;

	///	number of faces
		static const int num_faces = 0;

	///	number of volumes
		static const int num_volumes = 0;

	public:
	///	Constructor filling the arrays
		ReferenceVertex();

	/// \copydoc ug::ReferenceElement::reference_object_id()
		ReferenceObjectID reference_object_id() const {return REFERENCE_OBJECT_ID;}

	/// \copydoc ug::ReferenceElement::dimension()
		int dimension() const {return dim;}

	/// \copydoc ug::ReferenceElement::size()
		number size() const	{return 1.0;}
};

} // end namespace ug

#endif /* __H__UG__LIB_DISC__REFERENCE_ELEMENT__REFERENCE_VERTEX__ */
