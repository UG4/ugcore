/*
 * piecewise_constant.h
 *
 * Created on: 19.06.2012
 * Author: Christian Wehner
 */

#ifndef __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__CROUZEIX_RAVIART__PIECEWISE_CONSTANT__
#define __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__CROUZEIX_RAVIART__PIECEWISE_CONSTANT__

#include "../common/lagrange1d.h"
#include "../local_shape_function_set.h"
#include "../local_dof_set.h"
#include "lib_disc/common/multi_index.h"
#include "common/util/provider.h"
#include "common/util/metaprogramming_util.h"
#include "lib_grid/grid/geometric_base_objects.h"

namespace ug{

/// Elementwise constant shape functions
template <typename TRefElem>
class PiecewiseConstantLSFS
	: public BaseLocalShapeFunctionSet<PiecewiseConstantLSFS<TRefElem>, TRefElem::dim>
{
	private:
	///	base class
		typedef BaseLocalShapeFunctionSet<PiecewiseConstantLSFS<TRefElem>, TRefElem::dim> base_type;

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
		static const size_t order = 0;

	///	Dimension, where shape functions are defined
		static const int dim = reference_element_type::dim;

	/// Number of shape functions
		static const size_t nsh = 1;

	protected:
		MathVector<dim> bary;

	public:
	///	Constructor
		PiecewiseConstantLSFS()
		{
			const TRefElem& rRef = Provider<TRefElem>::get();
			//	get corner position integer
			int num = rRef.num(0);

			bary = rRef.corner(0);
            for (int j=1;j<num;j++){
            	bary+=rRef.corner(j);
            }
            bary*=1./(number)num;
		}

	///	\copydoc ug::LocalShapeFunctionSet::type()
		inline static LFEID type() {return LFEID(LFEID::PIECEWISE_CONSTANT, dim, 0);}

	///	\copydoc ug::LocalShapeFunctionSet::continuous()
		inline static bool continuous() {return false;}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		inline size_t num_sh() const {return 1;}

	///	\copydoc ug::LocalShapeFunctionSet::position()
		inline bool position(size_t i, position_type& pos) const
		{
			for(int d = 0; d < dim; d++)
				pos[d]=bary[d];
			return true;
		}

	///	\copydoc ug::LocalShapeFunctionSet::shape()
		inline number shape(const size_t i, const MathVector<dim>& x) const
		{
			return 1;
		}

	///	\copydoc ug::LocalShapeFunctionSet::shape()
		inline void grad(grad_type& g, const size_t i,	const position_type& x) const
		{
			TRefElem::check_position(x);
			VecSet(g, 0.0);
		}
};

} //namespace ug

#endif /* __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__CROUZEIX_RAVIART__PIECEWISE_CONSTANT__ */

