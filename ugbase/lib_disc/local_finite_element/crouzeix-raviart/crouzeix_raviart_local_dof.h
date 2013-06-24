/*
 * crouzeix_raviart_local_dof.h
 *
 * Created on: 18.06.2012
 * Author: Christian Wehner
 */

#ifndef __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__LAGRANGE__CROUZEIX_RAVIART_LOCAL_DOF__
#define __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__LAGRANGE__CROUZEIX_RAVIART_LOCAL_DOF__

#include "common/util/provider.h"
#include "../local_dof_set.h"
#include "lib_disc/reference_element/reference_element_util.h"
#include "lib_disc/common/multi_index.h"

namespace ug{

///////////////////////////////////////////////////////////////////////////////
// Crouzeix - Raviart Set
///////////////////////////////////////////////////////////////////////////////

/// Crouzeix - Raviart Set
template <typename TRefElem>
class CrouzeixRaviartLDS : public LocalDoFSet
{
	protected:
	///	dimension of reference element
		static const int refDim = TRefElem::dim;

	public:
	///	constructor
		CrouzeixRaviartLDS()
		{
			const TRefElem& rRefElem = Provider<TRefElem>::get();

			p = 1;
			if(refDim > 0)
			{
				nsh = rRefElem.num(refDim-1);

			//	set local DoFs (all located at dim-1 objects)
				m_vLocalDoF.resize(nsh);
				for(size_t i = 0; i< rRefElem.num(refDim-1); ++i)
					m_vLocalDoF[i] = LocalDoF(refDim-1, i, 0);
			}
			else
			{
				nsh = 0;
				m_vLocalDoF.clear();
			}
		}

	///	returns the reference dimension
		int dim() const {return refDim;}

	///	returns the type of reference element
		ReferenceObjectID roid() const {return TRefElem::REFERENCE_OBJECT_ID;}

	///	returns the total number of DoFs on the finite element
		size_t num_dof() const {return nsh;};

	///	returns the number of DoFs on a sub-geometric object type
		int num_dof(ReferenceObjectID type) const
		{
			const int d = ReferenceElementDimension(type);
			if(d == refDim-1)   return 1;
			else return 0;
		}

	///	returns the number of DoFs on sub-geometric object in dimension and id
		size_t num_dof(int d, size_t id) const
		{
			if(d == refDim-1) return 1;
			else return 0;
		}

	///	returns the dof storage
		const LocalDoF& local_dof(size_t dof) const {return m_vLocalDoF[dof];}

	protected:
	///	number of shapes
		size_t nsh;

	///	order
		size_t p;

	///	association to elements
		std::vector<LocalDoF> m_vLocalDoF;
};


} //namespace ug

#endif /* __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__LAGRANGE__CROUZEIX_RAVIART_LOCAL_DOF__ */
