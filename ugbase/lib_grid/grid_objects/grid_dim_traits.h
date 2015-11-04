#ifndef __H__LIB_GRID__GRID_DIM_TRAITS__
#define __H__LIB_GRID__GRID_DIM_TRAITS__

#include <boost/mpl/list.hpp>
#include "lib_grid/grid_objects/grid_objects.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
//	Element types for different dimensions
////////////////////////////////////////////////////////////////////////////////

/**
 * This class template provides the basic element types and their topological
 * properties depending on the topological dimensionality. Besides that, it
 * provides the lists of the element types for every topological dimension.
 * While DimElemList returns the Element Types in the dimenion of the domain,
 * the list AllElemList returns all elements contained in the Domain-dimension
 * and the dimensions below.
 */
template <int dim> struct grid_dim_traits;

// 0d
template <> struct grid_dim_traits<0>
{
	typedef boost::mpl::list<RegularVertex> DimElemList;
	typedef boost::mpl::list<RegularVertex> AllElemList;

	typedef geometry_traits<Vertex>::const_iterator const_iterator;
	typedef geometry_traits<Vertex>::iterator iterator;

	typedef geometry_traits<Vertex>::grid_base_object grid_base_object;
	typedef geometry_traits<Vertex>::grid_base_object element_type;

	const static size_t MaxNumVerticesOfElem = 1;
};

// 1d
template <> struct grid_dim_traits<1>
{
	typedef boost::mpl::list<RegularEdge> DimElemList;
	typedef boost::mpl::list<RegularEdge> AllElemList;
	typedef boost::mpl::list<> ManifoldElemList;

	typedef geometry_traits<Edge>::const_iterator const_iterator;
	typedef geometry_traits<Edge>::iterator iterator;

	typedef geometry_traits<Edge>::grid_base_object grid_base_object;
	typedef geometry_traits<Edge>::grid_base_object element_type;
	typedef geometry_traits<Vertex>::grid_base_object side_type;

	const static size_t MaxNumVerticesOfElem = 2;
};

// 2d
template <> struct grid_dim_traits<2>
{
	typedef boost::mpl::list<Triangle, Quadrilateral> DimElemList;
	typedef boost::mpl::list<RegularEdge, Triangle, Quadrilateral> AllElemList;
	typedef boost::mpl::list<Edge> ManifoldElemList;

	typedef geometry_traits<Face>::const_iterator const_iterator;
	typedef geometry_traits<Face>::iterator iterator;

	typedef geometry_traits<Face>::grid_base_object grid_base_object;
	typedef geometry_traits<Face>::grid_base_object element_type;
	typedef geometry_traits<Edge>::grid_base_object side_type;

	const static size_t MaxNumVerticesOfElem = 4;
};

// 3d
template <> struct grid_dim_traits<3>
{
	typedef boost::mpl::list<Tetrahedron, Prism, Pyramid, Hexahedron, Octahedron> DimElemList;
	typedef boost::mpl::list<RegularEdge, Triangle, Quadrilateral,
								Tetrahedron, Prism, Pyramid, Hexahedron, Octahedron> AllElemList;
	typedef boost::mpl::list<Triangle, Quadrilateral> ManifoldElemList;

	typedef geometry_traits<Volume>::const_iterator const_iterator;
	typedef geometry_traits<Volume>::iterator iterator;

	typedef geometry_traits<Volume>::grid_base_object grid_base_object;
	typedef geometry_traits<Volume>::grid_base_object element_type;
	typedef geometry_traits<Face>::grid_base_object side_type;

	const static size_t MaxNumVerticesOfElem = 8;
};

} // end namespace ug

#endif // __H__LIB_GRID__GRID_DIM_TRAITS__

/* End of File */
