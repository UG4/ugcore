/*
 * trialspace_impl.h
 *
 *  Created on: 16.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__LOCAL_DOF_PATTERN_IMPL__
#define __H__LIBDISCRETIZATION__LOCAL_DOF_PATTERN_IMPL__

namespace ug{

///////////////////////////////////////////
// Local DoF Pattern

template <typename TRefElem>
inline
int
LocalDoFPattern<TRefElem>::
total_num_dofs() const
{
	return m_total_num_dofs;
}

template <typename TRefElem>
inline
int
LocalDoFPattern<TRefElem>::
num_dofs(ReferenceElementType type) const
{
	return m_num_dof[type];
}

// sets the number of dofs
template <typename TRefElem>
void
LocalDoFPattern<TRefElem>::
set_num_dofs(ReferenceElementType type, int num)
{
	m_num_dof[type] = num;
}


template <typename TRefElem>
LocalDoFPattern<TRefElem>::
LocalDoFPattern() : m_total_num_dofs(0)
{
	for(int i = 0; i < NUM_REFERENCE_ELEMENTS; ++i)
	{
		m_num_dof[i] = 0;
	}
}


///////////////////////////////////////////
// ContinuousDoFPattern

inline
int
ContinuousDoFPattern::
total_num_dofs(ReferenceElementType type) const
{
	return m_total_num_dofs[type];
}

// returns the number of dofs needed on a Reference Element
inline
int
ContinuousDoFPattern::
num_dofs(ReferenceElementType type) const
{
	return m_num_dofs[type];
}


template <typename TRefElem>
bool
ContinuousDoFPattern::
add_local_dof_pattern(const LocalDoFPattern<TRefElem>& p)
{
	TRefElem refElem;
	int dim = refElem.dimension();

	// return false if initialized, by wrong pattern added
	if(m_dim != dim)
	{
		if(m_dim != -1)
		{
			return false;
		}
		else
		{
			m_dim = dim;
			for(int d = 0; d < dim; ++d)
			{
				for(int i = 0; i < (int)refElem.num_obj(d); ++i)
				{
					ReferenceElementType id = refElem.ref_elem_type(d, i);
					m_num_dofs[id] = p.num_dofs(id);
				}
			}
		}
	}
	else
	{
		// already initiallised. Look if no conflict apparent, add dofs
		for(int d = 0; d < dim; ++d)
		{
			for(int i = 0; i < (int) refElem.num_obj(d); ++i)
			{
				ReferenceElementType id = refElem.ref_elem_type(d, i);
				if(m_num_dofs[id] != p.num_dofs(id))
				{
					if(m_num_dofs[id] != -1) return false;
					m_num_dofs[id] = p.num_dofs(id);
				}
			}
		}
	}

	// update total num dofs (for element that is added)
	ReferenceElementType id = reference_element_traits<TRefElem>::REFERENCE_ELEMENT_TYPE;
	m_total_num_dofs[id] = m_num_dofs[id];
	for(int d = 0; d < dim; ++d)
	{
		for(int i = 0; i < (int)refElem.num_obj(d); ++i)
		{
			ReferenceElementType id2 = refElem.ref_elem_type(d, i);
			m_total_num_dofs[id] += refElem.num_ref_elem(id2) * m_num_dofs[id2];
		}
	}

	return true;
};

} // namespace ug


#endif /* __H__LIBDISCRETIZATION__LOCAL_DOF_PATTERN_IMPL__ */
