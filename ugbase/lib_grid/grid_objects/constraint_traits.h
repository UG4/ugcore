#ifndef __H__UG_LIB_GRID__CONSTRAINT_TRAITS__
#define __H__UG_LIB_GRID__CONSTRAINT_TRAITS__

#include "grid_objects_0d.h"
#include "grid_objects_1d.h"
#include "grid_objects_2d.h"

namespace ug{
/**	constraint traits provide the associated constrained and constraining grid
 * grid object types to a given grid object type.*/
template <class TElem>
struct constraint_traits{
	typedef TElem	elem_t;
	typedef void	constrained_t;
	typedef void	constraining_t;
};

template <>
struct constraint_traits<Vertex>{
	typedef Vertex				elem_t;
	typedef ConstrainedVertex	constrained_t;
	typedef void				constraining_t;
};

template <>
struct constraint_traits<Edge>{
	typedef Edge				elem_t;
	typedef ConstrainedEdge		constrained_t;
	typedef ConstrainingEdge	constraining_t;
};

template <>
struct constraint_traits<Face>{
	typedef Face				elem_t;
	typedef ConstrainedFace		constrained_t;
	typedef ConstrainingFace	constraining_t;
};

template <>
struct constraint_traits<Triangle>{
	typedef Triangle				elem_t;
	typedef ConstrainedTriangle		constrained_t;
	typedef ConstrainingTriangle	constraining_t;
};

template <>
struct constraint_traits<Quadrilateral>{
	typedef Quadrilateral				elem_t;
	typedef ConstrainedQuadrilateral	constrained_t;
	typedef ConstrainingQuadrilateral	constraining_t;
};

}//	end of namespace

#endif
