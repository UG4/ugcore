/*
 * domain_traits.h
 *
 *  Created on: 05.03.2012
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__DOMAIN_TRAITS__
#define __H__UG__LIB_DISC__DOMAIN_TRAITS__

#include <boost/mpl/list.hpp>
#include <boost/mpl/for_each.hpp>
#include "lib_grid/grid_objects/grid_objects.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
//	Element types for different dimensions
////////////////////////////////////////////////////////////////////////////////

/**
 * This Class provides boost::mpl::lists storing the type of elements used in
 * for the Domain. It can be used to control dimension dependent builds, where
 * not all template instantiations are available (e.g. a Hexahedron in 1d,2d, etc)
 * While DimElemList returns the Element Types in the dimenion of the domain,
 * the list AllElemList returns all elements contained in the Domain-dimension
 * and the dimensions below.
 */
template <int dim> struct domain_traits;

// 0d
template <> struct domain_traits<0> {
typedef boost::mpl::list<RegularVertex> DimElemList;
typedef boost::mpl::list<RegularVertex> AllElemList;

typedef geometry_traits<Vertex>::const_iterator const_iterator;
typedef geometry_traits<Vertex>::iterator iterator;

typedef geometry_traits<Vertex>::grid_base_object grid_base_object;
typedef geometry_traits<Vertex>::grid_base_object element_type;

const static size_t MaxNumVerticesOfElem = 1;
};

// 1d
template <> struct domain_traits<1> {
typedef boost::mpl::list<RegularEdge> DimElemList;
typedef boost::mpl::list<RegularEdge> AllElemList;
typedef boost::mpl::list<> ManifoldElemList;

typedef geometry_traits<Edge>::const_iterator const_iterator;
typedef geometry_traits<Edge>::iterator iterator;

typedef geometry_traits<Edge>::grid_base_object grid_base_object;
typedef geometry_traits<Edge>::grid_base_object element_type;
typedef geometry_traits<Vertex>::grid_base_object side_type;

const static size_t MaxNumVerticesOfElem = 2;

typedef MathVector<1> position_type;
typedef Attachment<position_type> position_attachment_type;
typedef Grid::VertexAttachmentAccessor<position_attachment_type> position_accessor_type;

};

// 2d
template <> struct domain_traits<2> {
typedef boost::mpl::list<Triangle, Quadrilateral> DimElemList;
typedef boost::mpl::list<RegularEdge, Triangle, Quadrilateral> AllElemList;
typedef boost::mpl::list<RegularEdge> ManifoldElemList;

typedef geometry_traits<Face>::const_iterator const_iterator;
typedef geometry_traits<Face>::iterator iterator;

typedef geometry_traits<Face>::grid_base_object grid_base_object;
typedef geometry_traits<Face>::grid_base_object element_type;
typedef geometry_traits<Edge>::grid_base_object side_type;

const static size_t MaxNumVerticesOfElem = 4;

typedef MathVector<2> position_type;
typedef Attachment<position_type> position_attachment_type;
typedef Grid::VertexAttachmentAccessor<position_attachment_type> position_accessor_type;
};

// 3d
template <> struct domain_traits<3> {
typedef boost::mpl::list<Tetrahedron, Prism, Pyramid, Hexahedron> DimElemList;
typedef boost::mpl::list<RegularEdge, Triangle, Quadrilateral,
								 Tetrahedron, Prism, Pyramid, Hexahedron> AllElemList;
typedef boost::mpl::list<Triangle, Quadrilateral> ManifoldElemList;

typedef geometry_traits<Volume>::const_iterator const_iterator;
typedef geometry_traits<Volume>::iterator iterator;

typedef geometry_traits<Volume>::grid_base_object grid_base_object;
typedef geometry_traits<Volume>::grid_base_object element_type;
typedef geometry_traits<Face>::grid_base_object side_type;

const static size_t MaxNumVerticesOfElem = 8;

typedef MathVector<3> position_type;
typedef Attachment<position_type> position_attachment_type;
typedef Grid::VertexAttachmentAccessor<position_attachment_type> position_accessor_type;

};


} // end namespace ug

#endif /* __H__UG__LIB_DISC__DOMAIN_TRAITS__ */
