/*
 * lagrange_local_dof.h
 *
 *  Created on: 27.06.2011
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISCRETIZATION__LOCAL_SHAPE_FUNCTION_SET__LAGRANGE__LAGRANGE_LOCAL_DOF__
#define __H__UG__LIB_DISCRETIZATION__LOCAL_SHAPE_FUNCTION_SET__LAGRANGE__LAGRANGE_LOCAL_DOF__

#include "lagrange.h"
#include "../local_dof.h"
#include "lib_discretization/common/multi_index.h"

namespace ug{

/// Lagrange DoF Pattern
template <typename TRefElem, int TOrder>
struct LagrangeLDP{};

/// specialization for Edges
/**
 * Lagrange shape function of any order for the Reference Edge
 * \tparam 	TOrder		requested order
 */
template <>
template <int TOrder>
class LagrangeLDP<ReferenceEdge, TOrder>
{
	protected:
	///	corresponding local shape function set
		typedef LagrangeLSFS<ReferenceEdge, TOrder> LSFS;

	///	number of shapes
		static const size_t nsh = LSFS::nsh;

	public:
	///	constructor
		LagrangeLDP()
		{
		//	on vertex
			m_vLocalDoF[0] = LocalDoF(0, 1, 0);
			m_vLocalDoF[nsh-1] = LocalDoF(0, 0, 0);

		//	on edge
			for(size_t sh = 1; sh < nsh-1; ++sh)
				m_vLocalDoF[sh] = LocalDoF(1, 0, sh-1);
		}

	///	returns the total number of DoFs on the finite element
		static inline size_t num_sh() {return nsh;};

	///	returns the dof storage
		inline const LocalDoF& storage(size_t sh) const
			{return m_vLocalDoF[sh];}

	///	returns if the storage needs objects of a given dimension
		static inline bool storage_use(int dim)
		{
				 if(dim == 1) return true;
			else if(dim == 2) return nsh > 2;
			else return false;
		}

		static inline size_t map_offset(size_t offset, std::vector<size_t>& vNodeOrder)
		{
			if(vNodeOrder[0] < vNodeOrder[1]) return offset;
			else return nsh-1 - offset;
		}

		static inline size_t num_sh(int dim, size_t i)
		{
			if(dim == 1) return 1;
			else if (dim == 2) return nsh - 2;
			else return 0;
		}

	protected:
	///	association to elements
		LocalDoF m_vLocalDoF[nsh];
};

} //namespace ug

#endif /* __H__UG__LIB_DISCRETIZATION__LOCAL_SHAPE_FUNCTION_SET__LAGRANGE__LAGRANGE_LOCAL_DOF__ */
