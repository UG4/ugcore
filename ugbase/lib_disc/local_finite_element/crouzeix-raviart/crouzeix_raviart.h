/*
 * crouzeix_raviart.h
 */

#ifndef __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__CROUZEIX_RAVIART__CROUZEIX_RAVIART__
#define __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__CROUZEIX_RAVIART__CROUZEIX_RAVIART__

#include "../common/lagrange1d.h"
#include "../local_shape_function_set.h"
#include "../local_dof_set.h"
#include "lib_disc/common/multi_index.h"
#include "common/util/provider.h"
#include "common/util/metaprogramming_util.h"
#include "lib_grid/grid/geometric_base_objects.h"

namespace ug{

/// Lagrange Shape Function Set without virtual functions and fixed order
template <typename TRefElem>
class CrouzeixRaviartLSFS
	: public BaseLocalShapeFunctionSet<CrouzeixRaviartLSFS<TRefElem>, TRefElem::dim>
{
	private:
	///	abbreviation for order
		static const size_t p = 1;

	///	base class
		typedef BaseLocalShapeFunctionSet<CrouzeixRaviartLSFS<TRefElem>, TRefElem::dim> base_type;

	public:
	///	Domain position type
		typedef typename base_type::position_type position_type;

	///	Shape type
		typedef typename base_type::shape_type shape_type;

	///	Gradient type
		typedef typename base_type::grad_type grad_type;

	public:
	///	Reference Element type
		typedef TRefElem reference_element_type;

	///	Order of Shape functions
		static const size_t order = p;

	///	Dimension, where shape functions are defined
		static const int dim = reference_element_type::dim;

	/// Number of shape functions
		size_t nsh;

	public:
	///	Constructor
		CrouzeixRaviartLSFS()
		{
			const TRefElem& rRefElem = Provider<TRefElem>::get();

		//	set number of shapes
			if(dim > 0) nsh = rRefElem.num(dim-1);
			else nsh = 0;
		}

	///	\copydoc ug::LocalShapeFunctionSet::type()
		inline static LFEID type() {return LFEID(LFEID::CROUZEIX_RAVIART, 1);}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		inline size_t num_sh() const {return nsh;}

	///	\copydoc ug::LocalShapeFunctionSet::position()
		inline bool position(size_t i, position_type& pos) const
		{
			return false;
		}

	///	\copydoc ug::LocalShapeFunctionSet::shape()
		inline number shape(const size_t i, const MathVector<dim>& x) const
		{
			TRefElem::check_position(x);

		// \TODO: THIS IS SIMPLICES ONLY !!!!
		// \TODO: RETHINK, MIGHT BE WRONG !!!!
			if(i == 0) {
				number x0 = 1.;
				for(int d = 0; d < dim; ++d) x0 -= x[d];
				return dim*(1./dim - x0);
			}
			else return dim*(1./dim - x[i]);
		}

	///	\copydoc ug::LocalShapeFunctionSet::shape()
		inline void grad(grad_type& g, const size_t i,	const position_type& x) const
		{
			TRefElem::check_position(x);

			if(i==0) {
				VecSet(g, dim);
			}
			else
			{
				VecSet(g, 0.0);
				g[i] = -dim;
			}
		}
};



} //namespace ug

#endif /* __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__CROUZEIX_RAVIART__CROUZEIX_RAVIART__ */

