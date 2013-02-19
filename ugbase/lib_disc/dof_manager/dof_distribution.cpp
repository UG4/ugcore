/*
 * dof_distribution.cpp
 *
 *  Created on: 24.01.2012
 *      Author: andreasvogel
 */

#include "dof_distribution.h"

namespace ug{

void DoFDistribution::add_transfer(SmartPtr<ILocalTransfer> spTransfer)
{
	for(int gbo = 0; gbo < NUM_GEOMETRIC_BASE_OBJECTS; ++gbo)
	{
		if(spTransfer->prolongation_needed((GeometricBaseObject)gbo))
			m_vProlongation[gbo].push_back(spTransfer);

		if(spTransfer->restriction_needed((GeometricBaseObject)gbo))
			m_vRestriction[gbo].push_back(spTransfer);
	}
}

void DoFDistribution::remove_transfer(SmartPtr<ILocalTransfer> spTransfer)
{
	for(int gbo = 0; gbo < NUM_GEOMETRIC_BASE_OBJECTS; ++gbo)
	{
		m_vProlongation[gbo].erase(std::remove(m_vProlongation[gbo].begin(),
		                                       m_vProlongation[gbo].end(),
		                                       spTransfer),
		                                       m_vProlongation[gbo].end());
		m_vRestriction[gbo].erase(std::remove(m_vRestriction[gbo].begin(),
		                                      m_vRestriction[gbo].end(),
		                                       spTransfer),
		                                       m_vRestriction[gbo].end());
	}
}

void DoFDistribution::clear_transfers()
{
	for(int gbo = 0; gbo < NUM_GEOMETRIC_BASE_OBJECTS; ++gbo)
	{
		m_vProlongation[gbo].clear();
		m_vRestriction[gbo].clear();
	}
}

} // end namespace ug
