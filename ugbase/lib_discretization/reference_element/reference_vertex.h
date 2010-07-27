/*
 * reference_vertex.h
 *
 *  Created on: 15.06.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__REFERENCE_ELEMENT__REFERENCE_VERTEX__
#define __H__LIBDISCRETIZATION__REFERENCE_ELEMENT__REFERENCE_VERTEX__

namespace ug{

class ReferenceVertex
{
public:
	static const ReferenceObjectID REFERENCE_OBJECT_ID = ROID_VERTEX;

	static const int dim = 0;
	static const int num_corners = 1;
	static const int num_edges = 0;
	static const int num_faces = 0;
	static const int num_volumes = 0;

	// TODO: Implement
};

template <>
class reference_element_traits<VertexBase>
{
	public:
		typedef ReferenceVertex reference_element_type;
};

} // end namespace ug

#endif /* __H__LIBDISCRETIZATION__REFERENCE_ELEMENT__REFERENCE_VERTEX__ */
