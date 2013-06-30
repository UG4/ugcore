/*
 * piecewise_constant_local_dof.h
 *
 * Created on: 19.06.2012
 * Author: Christian Wehner
 */

#ifndef __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__LAGRANGE__CROUZEIX_RAVIART__PIECEWISE_CONSTANT_LOCAL_DOF__
#define __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__LAGRANGE__CROUZEIX_RAVIART__PIECEWISE_CONSTANT_LOCAL_DOF__

#include "common/util/provider.h"
#include "../local_dof_set.h"
#include "lib_disc/reference_element/reference_element_util.h"
#include "lib_disc/common/multi_index.h"

namespace ug{

///////////////////////////////////////////////////////////////////////////////
// Piecewise constant Set
///////////////////////////////////////////////////////////////////////////////

/// Piecewise constant Set
template <typename TRefElem>
class PiecewiseConstantLDS : public LocalDoFSet
{
	protected:
	///	dimension of reference element
		static const int dim = TRefElem::dim;

	public:
	///	constructor
		PiecewiseConstantLDS() {m_vLocalDoF = LocalDoF(dim, 0, 0);}

	///	\copydoc ug::LocalShapeFunctionSet::type()
		inline LFEID type() const {return LFEID(LFEID::PIECEWISE_CONSTANT, dim, 0);}

	///	returns the type of reference element
		ReferenceObjectID roid() const {return TRefElem::REFERENCE_OBJECT_ID;}

	///	returns the total number of DoFs on the finite element
		size_t num_dof() const {return 1;};

	///	returns the number of DoFs on a sub-geometric object type
		size_t num_dof(ReferenceObjectID type) const
		{
			if (type == TRefElem::REFERENCE_OBJECT_ID)   return 1;
			else return 0;
		}

	///	returns the dof storage
		const LocalDoF& local_dof(size_t dof) const {return m_vLocalDoF;}

	///	returns if the local dof position are exact
		bool exact_position_available() const {return true;};

	protected:
		LocalDoF m_vLocalDoF; ///< association to elements
};


} //namespace ug

#endif /* __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__LAGRANGE__CROUZEIX_RAVIART__PIECEWISE_CONSTANT_LOCAL_DOF__ */
