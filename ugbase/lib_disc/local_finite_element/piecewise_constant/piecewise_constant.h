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
#include "piecewise_constant_local_dof.h"
#include "lib_disc/common/multi_index.h"
#include "common/util/provider.h"
#include "common/util/metaprogramming_util.h"
#include "lib_grid/grid/geometric_base_objects.h"

namespace ug{

/// Elementwise constant shape functions
template <typename TRefElem>
class PiecewiseConstantLSFS
	: public PiecewiseConstantLDS<TRefElem>,
	  public BaseLocalShapeFunctionSet<PiecewiseConstantLSFS<TRefElem>, TRefElem::dim>
{
	private:
	///	base class
		typedef BaseLocalShapeFunctionSet<PiecewiseConstantLSFS<TRefElem>, TRefElem::dim> base_type;

	public:
	///	Shape type
		typedef typename base_type::shape_type shape_type;

	///	Gradient type
		typedef typename base_type::grad_type grad_type;

	public:
	///	Dimension, where shape functions are defined
		static const int dim = TRefElem::dim;

	protected:
	///	barycenter
		MathVector<dim> bary;

	public:
	///	Constructor
		PiecewiseConstantLSFS()
		{
			const TRefElem& rRef = Provider<TRefElem>::get();

			bary = rRef.corner(0);
            for (size_t j=1; j < rRef.num(0); ++j){
            	bary += rRef.corner(j);
            }
            bary *= 1./rRef.num(0);
		}

	///	\copydoc ug::LocalShapeFunctionSet::continuous()
		inline bool continuous() const {return false;}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		inline size_t num_sh() const {return 1;}

	///	\copydoc ug::LocalShapeFunctionSet::position()
		inline bool position(size_t i, MathVector<dim>& pos) const
		{
			pos = bary; return true;
		}

	///	\copydoc ug::LocalShapeFunctionSet::shape()
		inline number shape(const size_t i, const MathVector<dim>& x) const
		{
			return 1;
		}

	///	\copydoc ug::LocalShapeFunctionSet::shape()
		inline void grad(MathVector<dim>& g, const size_t i, const MathVector<dim>& x) const
		{
			TRefElem::check_position(x);
			VecSet(g, 0.0);
		}
};

} //namespace ug

#endif /* __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__CROUZEIX_RAVIART__PIECEWISE_CONSTANT__ */

