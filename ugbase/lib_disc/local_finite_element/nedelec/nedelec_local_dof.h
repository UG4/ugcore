/*
 * nedelec_local_dof.h
 *
 * This file contains implementations of the dof distribution patterns for
 * the so-called Nedelec (or Whitney-1) elements: The edge dofs.
 *
 * Created on: 28.08.2012
 * Author: Dmitry Logashenko
 */

#ifndef __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__NEDELEC__NEDELEC_LOCAL_DOF__
#define __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__NEDELEC__NEDELEC_LOCAL_DOF__

#include "common/util/provider.h"
#include "../local_dof_set.h"
#include "lib_disc/common/multi_index.h"

namespace ug{

///////////////////////////////////////////////////////////////////////////////
// Nedelec Set
///////////////////////////////////////////////////////////////////////////////

/// Nedelec, i.e. the edge local dof set
template <typename TRefElem>
class NedelecLDS : public LocalDoFSet
{
	protected:
	///	dimension of reference element
		static const int refDim = TRefElem::dim;

	public:
	///	constructor
		NedelecLDS()
		{
			if(refDim < 2)
			{
			//	No dofs if the dimension is less than 2:
				nsh = 0;
				m_vLocalDoF.clear();
				return;
			}
			
			const TRefElem& rRefElem = Provider<TRefElem>::get();
			
			nsh = rRefElem.num(1); // number of the edges
		//	set local DoFs (all located at the edges)
			m_vLocalDoF.resize(nsh);
			for(size_t i = 0; i < nsh; ++i)
				m_vLocalDoF[i] = LocalDoF(1, i, 0);
		}

	///	returns the type of reference element
		ReferenceObjectID roid() const {return TRefElem::REFERENCE_OBJECT_ID;}

	///	\copydoc ug::LocalShapeFunctionSet::type()
		inline LFEID type() const {return LFEID(LFEID::NEDELEC, refDim, 1);}

	///	returns the total number of DoFs on the finite element
		size_t num_dof() const {return nsh;};

	///	returns the number of DoFs on a sub-geometric object type
		size_t num_dof(ReferenceObjectID type) const
		{
			if(ReferenceElementDimension(type) == 1) return 1;
			else return 0;
		}

	///	returns the dof storage
		const LocalDoF& local_dof(size_t dof) const {return m_vLocalDoF[dof];}

	///	returns if the local dof position are exact
		bool exact_position_available() const {return true;};

	protected:
	///	number of shapes (== number of edges)
		size_t nsh;

	///	association to elements
		std::vector<LocalDoF> m_vLocalDoF;
};

} // namespace ug

#endif /* __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__WHITNEY__LAGRANGE_LOCAL_DOF__ */
