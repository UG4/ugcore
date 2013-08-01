/*
 * dof_distribution.cpp
 *
 *  Created on: 24.01.2012
 *      Author: andreasvogel
 */

#include "dof_distribution.h"
#include "lib_disc/function_spaces/grid_function.h"

namespace ug{


////////////////////////////////////////////////////////////////////////////////
// DoFDistribution
////////////////////////////////////////////////////////////////////////////////


void DoFDistribution::manage_grid_function(IGridFunction& gridFct)
{
//	if already registered, we're done
	if(std::find(m_vpGridFunction.begin(), m_vpGridFunction.end(), &gridFct)
		!= m_vpGridFunction.end())
		return;

//	add to managed functions
	m_vpGridFunction.push_back(&gridFct);
}

void DoFDistribution::unmanage_grid_function(IGridFunction& gridFct)
{
	m_vpGridFunction.erase(
	std::remove(m_vpGridFunction.begin(), m_vpGridFunction.end(), &gridFct)
	, m_vpGridFunction.end());
}

void DoFDistribution::permute_values(const std::vector<size_t>& vIndNew)
{
//	swap values of handled grid functions
	for(size_t i = 0; i < m_vpGridFunction.size(); ++i)
		m_vpGridFunction[i]->permute_values(vIndNew);
}

void DoFDistribution::copy_values(const std::vector<std::pair<size_t, size_t> >& vIndexMap,
                                         bool bDisjunct)
{
//	swap values of handled grid functions
	for(size_t i = 0; i < m_vpGridFunction.size(); ++i)
		m_vpGridFunction[i]->copy_values(vIndexMap, bDisjunct);
}

void DoFDistribution::resize_values(size_t newSize)
{
//	swap values of handled grid functions
	for(size_t i = 0; i < m_vpGridFunction.size(); ++i)
		m_vpGridFunction[i]->resize_values(newSize);
}

} // end namespace ug
