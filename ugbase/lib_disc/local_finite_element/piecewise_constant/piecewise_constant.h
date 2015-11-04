
#ifndef __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__PIECEWISE_CONSTANT__
#define __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__PIECEWISE_CONSTANT__

#include "common/util/provider.h"
#include "lib_grid/grid/grid_base_objects.h"
#include "lib_disc/local_finite_element/local_dof_set.h"
#include "lib_disc/reference_element/reference_element_util.h"

namespace ug{

/// Elementwise constant shape functions
template <typename TRefElem>
class PiecewiseConstantLSFS
	: public BaseLSFS<PiecewiseConstantLSFS<TRefElem>, TRefElem::dim>
{
	public:
	///	Dimension, where shape functions are defined
		static const int dim = TRefElem::dim;

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

            m_vLocalDoF = LocalDoF(dim, 0, 0);
		}

	public:
	///	\copydoc ug::DimLocalDoFSet::roid()
		ReferenceObjectID roid() const {return TRefElem::REFERENCE_OBJECT_ID;}

	///	\copydoc ug::DimLocalDoFSet::num_dof()
		size_t num_dof() const {return 1;};

	///	\copydoc ug::DimLocalDoFSet::num_dof()
		size_t num_dof(ReferenceObjectID type) const
		{
			if (type == TRefElem::REFERENCE_OBJECT_ID)   return 1;
			else return 0;
		}

	///	\copydoc ug::DimLocalDoFSet::local_dof()
		const LocalDoF& local_dof(size_t dof) const {return m_vLocalDoF;}

	///	\copydoc ug::DimLocalDoFSet::position()
		inline bool position(size_t i, MathVector<dim>& pos) const
		{
			pos = bary; return true;
		}

	///	\copydoc ug::DimLocalDoFSet::exact_position_available()
		bool exact_position_available() const {return true;};

	public:
	///	\copydoc ug::LocalShapeFunctionSet::continuous()
		bool continuous() const {return false;}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		size_t num_sh() const {return 1;}

	///	\copydoc ug::LocalShapeFunctionSet::shape()
		number shape(const size_t i, const MathVector<dim>& x) const
		{
			return 1;
		}

	///	\copydoc ug::LocalShapeFunctionSet::shape()
		void grad(MathVector<dim>& g, const size_t i, const MathVector<dim>& x) const
		{
			TRefElem::check_position(x);
			VecSet(g, 0.0);
		}

	protected:
		MathVector<dim> bary; ///< Barycenter
		LocalDoF m_vLocalDoF; ///< association to elements
};

} //namespace ug

#endif /* __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__PIECEWISE_CONSTANT__ */

