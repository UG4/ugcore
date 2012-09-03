/*
 * crouzeix_raviart_local_dof.h
 *
 * Created on: 18.06.2012
 * Author: Christian Wehner
 */

#ifndef __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__MINI__MINI_BUBBLE_LOCAL_DOF__
#define __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__MINI__MINI_BUBBLE_LOCAL_DOF__

#include "common/util/provider.h"
#include "../local_dof_set.h"
#include "lib_disc/reference_element/reference_element_util.h"
#include "lib_disc/common/multi_index.h"

namespace ug{

///////////////////////////////////////////////////////////////////////////////
// MiniBubble Set
///////////////////////////////////////////////////////////////////////////////

/// MiniBubble Set (2D only!)
template <typename TRefElem>
class MiniBubbleLDS : public ILocalDoFSet
{
	protected:
	///	dimension of reference element
		static const int refDim = TRefElem::dim;

	public:
	///	constructor
		MiniBubbleLDS()
		{
			// get _the_ reference element
			const TRefElem& rRefElem = Provider<TRefElem>::get();

			if(refDim >= 2)
			{
				// face (or volume???)
				//	set local DoFs (located at vertices+bubble)
				nsh = rRefElem.num(0)+1;

				m_vLocalDoF.resize(nsh);
				for(size_t i = 0; i< nsh-1; ++i)
					m_vLocalDoF[i] = LocalDoF(0, i, 0);

				m_vLocalDoF[nsh-1] = LocalDoF(refDim, nsh-1, 0); // bubble located at element
			}
			else
			{
				// edge or vertex
				nsh = refDim+1;
				m_vLocalDoF.resize(nsh);
				for(size_t i = 0; i< nsh-1; ++i)
					m_vLocalDoF[i] = LocalDoF(0, i, 0);

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
			return max_num_dof(ReferenceElementDimension(type));
		}

	///	returns the number of DoFs on sub-geometric object in dimension and id
		size_t num_dof(int d, size_t id) const
		{
			// each element and vertex hold one dof
			if((d == refDim)||(d==0)) return 1;
			else return 0;
		}

	///	returns if the storage needs objects of a given dimension
		size_t max_num_dof(int d) const
		{
			if (d==0) return 1;         // vertices
			if (d == refDim)   return 1;    // element
			else if (d > refDim) return -1;
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
