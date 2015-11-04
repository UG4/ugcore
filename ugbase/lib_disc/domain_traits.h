
#ifndef __H__UG__LIB_DISC__DOMAIN_TRAITS__
#define __H__UG__LIB_DISC__DOMAIN_TRAITS__

#include <boost/mpl/for_each.hpp>
#include "lib_grid/grid_objects/grid_dim_traits.h"

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
template <> struct domain_traits<0> : public grid_dim_traits<0> {};

// 1d
template <> struct domain_traits<1> : public grid_dim_traits<1> {

typedef MathVector<1> position_type;
typedef Attachment<position_type> position_attachment_type;
typedef Grid::VertexAttachmentAccessor<position_attachment_type> position_accessor_type;

};

// 2d
template <> struct domain_traits<2> : public grid_dim_traits<2> {

typedef MathVector<2> position_type;
typedef Attachment<position_type> position_attachment_type;
typedef Grid::VertexAttachmentAccessor<position_attachment_type> position_accessor_type;

};

// 3d
template <> struct domain_traits<3> : public grid_dim_traits<3> {

typedef MathVector<3> position_type;
typedef Attachment<position_type> position_attachment_type;
typedef Grid::VertexAttachmentAccessor<position_attachment_type> position_accessor_type;

};


} // end namespace ug

#endif /* __H__UG__LIB_DISC__DOMAIN_TRAITS__ */
