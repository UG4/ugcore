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
		static const int refDim = TRefElem::dim;

	public:
	///	constructor
		PiecewiseConstantLDS()
		{
			if(refDim > 0)
			{
				nsh = 1;

			//	set local DoFs (all located at dim objects)
				m_vLocalDoF.resize(nsh);
				m_vLocalDoF[0] = LocalDoF(refDim, 0, 0);
			}
			else
			{
				nsh = 0;
				m_vLocalDoF.clear();
			}
		}

	///	returns the type of reference element
		ReferenceObjectID roid() const {return TRefElem::REFERENCE_OBJECT_ID;}

	///	returns the total number of DoFs on the finite element
		size_t num_dof() const {return nsh;};

	///	returns the number of DoFs on a sub-geometric object type
		size_t num_dof(ReferenceObjectID type) const
		{
			if (ReferenceElementDimension(type) == refDim)   return 1;
			else return 0;
		}

	///	returns the dof storage
		const LocalDoF& local_dof(size_t dof) const {return m_vLocalDoF[dof];}

	protected:
	///	number of shapes
		size_t nsh;

	///	order
		static const size_t p=0;

	///	association to elements
		std::vector<LocalDoF> m_vLocalDoF;
};


} //namespace ug

#endif /* __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__LAGRANGE__CROUZEIX_RAVIART__PIECEWISE_CONSTANT_LOCAL_DOF__ */
