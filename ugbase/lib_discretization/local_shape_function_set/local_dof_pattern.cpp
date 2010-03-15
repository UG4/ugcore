/*
 * local_dof_pattern.cpp
 *
 *  Created on: 12.05.2009
 *      Author: andreasvogel
 */

#include "local_dof_pattern.h"

namespace ug{

ContinuousDoFPattern&
ContinuousDoFPattern::
operator=(const ContinuousDoFPattern& item)
{
	m_dim = item.m_dim;
	for(int i = 0; i < NUM_REFERENCE_ELEMENTS; ++i)
	{
		m_num_dofs[i] = item.m_num_dofs[i];
		m_total_num_dofs[i] = item.m_total_num_dofs[i];
	}
	return *this;
}

ContinuousDoFPattern&
ContinuousDoFPattern::
operator+=(const ContinuousDoFPattern& item)
{
	assert(m_dim == item.m_dim);
	for(int i = 0; i < NUM_REFERENCE_ELEMENTS; ++i)
	{
		m_num_dofs[i] += item.m_num_dofs[i];
		m_total_num_dofs[i] += item.m_total_num_dofs[i];
	}
	return *this;
}

ContinuousDoFPattern&
ContinuousDoFPattern::
operator-=(const ContinuousDoFPattern& item)
{
	assert(m_dim == item.m_dim);
	for(int i = 0; i < NUM_REFERENCE_ELEMENTS; ++i)
	{
		m_num_dofs[i] -= item.m_num_dofs[i];
		m_total_num_dofs[i] -= item.m_total_num_dofs[i];
	}
	return *this;
}

ContinuousDoFPattern::
ContinuousDoFPattern(const ContinuousDoFPattern& item)
{
	m_dim = item.m_dim;
	for(int i = 0; i < NUM_REFERENCE_ELEMENTS; ++i)
	{
		m_num_dofs[i] = item.m_num_dofs[i];
		m_total_num_dofs[i] = item.m_total_num_dofs[i];
	}
}

ContinuousDoFPattern::ContinuousDoFPattern()
{
	m_dim = -1;
	for(int i = 0; i < NUM_REFERENCE_ELEMENTS; ++i)
	{
		m_num_dofs[i] = -1;
		m_total_num_dofs[i] = -1;
	}
}

bool operator==(const ContinuousDoFPattern& lhs, const ContinuousDoFPattern& rhs)
{
	if(lhs.m_dim != rhs.m_dim) return false;
	for(int i = 0; i < NUM_REFERENCE_ELEMENTS; ++i)
	{
		if(lhs.m_num_dofs[i] != rhs.m_num_dofs[i]) return false;
		if(lhs.m_total_num_dofs[i] != rhs.m_total_num_dofs[i]) return false;
	}
	return true;
}

}

