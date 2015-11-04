
#ifndef __H__UG__LIB_DISC__QUADRATURE__GAUSS_QUAD__GAUSS_QUAD__
#define __H__UG__LIB_DISC__QUADRATURE__GAUSS_QUAD__GAUSS_QUAD__

#include "common/common.h"
#include "../quadrature.h"
#include "lib_disc/reference_element/reference_element.h"

namespace ug{

/// fixed order gauss quadrature
template <typename TRefElem, int order>
class GaussQuadrature;

/// wrapper to ease implementation
template <typename TImpl, int TDim, int TOrder, int TNip>
class GaussQuadBase
{
	public:
	/// Dimension of integration domain
		static const size_t dim = TDim;

	/// Position Type in Reference Element Space
		typedef MathVector<dim> position_type;

	///	Type of weights
		typedef number weight_type;

	/// Order of quadrature rule
		static const size_t p = TOrder;

	/// Number of integration points
		static const size_t nip = TNip;

	public:
	/// number of integration points
		static size_t size() {return nip;}

	/// returns i'th integration point
		static const MathVector<dim>& point(size_t i)
			{UG_ASSERT(i < size(), "Wrong index"); return m_vPoint[i];}

	/// returns all positions in an array of size()
		static const MathVector<dim>* points() {return m_vPoint;}

	/// return the i'th weight
		static number weight(size_t i)
			{UG_ASSERT(i < size(), "Wrong index"); return m_vWeight[i];}

	/// returns all weights in an array of size()
		static const number* weights() {return m_vWeight;}

	/// returns the order
		static size_t order() {return p;}

	protected:
	/// integration points
		static MathVector<dim> m_vPoint[nip];

	/// weights
		static number m_vWeight[nip];
};

/// flexible order gauss quadrature
/**
 * Providing gauss quadrature for an reference element. This class wrapps a
 * static GaussQuadrature into the Quadrature interface.
 *
 * \tparam 		TRefElem		Reference Element Type
 */
template <typename TRefElem>
class FlexGaussQuadrature
	: public QuadratureRule<TRefElem::dim>
{
	public:
	///	Constructor
		FlexGaussQuadrature(int order);

	///	Destructor
		~FlexGaussQuadrature() {}
};

} // namespace ug

// include implementation
#include "gauss_quad_edge.h"
#include "gauss_quad_triangle.h"
#include "gauss_quad_quadrilateral.h"
#include "gauss_quad_tetrahedron.h"
#include "gauss_quad_pyramid.h"
#include "gauss_quad_prism.h"
#include "gauss_quad_hexahedron.h"
#include "gauss_quad_octahedron.h"


#endif /* __H__UG__LIB_DISC__QUADRATURE__GAUSS_QUAD__GAUSS_QUAD__ */
