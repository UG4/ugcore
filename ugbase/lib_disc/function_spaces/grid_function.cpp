/*
 * grid_function.cpp
 *
 *  Created on: 08.03.2012
 *      Author: andreasvogel
 */

#include "grid_function.h"
#include "lib_disc/dof_manager/surface_dof_distribution.h"
#include "lib_disc/dof_manager/level_dof_distribution.h"

namespace ug{


template <>
void IDDGridFunction<SurfaceDoFDistribution>::add_transfer(SmartPtr<ILocalTransfer> spTransfer)
{
	spTransfer->set_dof_distribution(m_spDD);
	m_vTransfer.push_back(spTransfer);
}

template <>
void IDDGridFunction<SurfaceDoFDistribution>::remove_transfer(SmartPtr<ILocalTransfer> spTransfer)
{
	m_vTransfer.erase(std::remove(m_vTransfer.begin(),
	                              m_vTransfer.end(),
	                              spTransfer),
	                              m_vTransfer.end());
}

template <>
void IDDGridFunction<SurfaceDoFDistribution>::clear_transfers()
{
	m_vTransfer.clear();
}

template <>
void IDDGridFunction<LevelDoFDistribution>::add_transfer(SmartPtr<ILocalTransfer> transfer)
{
	UG_THROW_FATAL("No Transfer operations for Level DoF Distributions.");
}
template <>
void IDDGridFunction<LevelDoFDistribution>::remove_transfer(SmartPtr<ILocalTransfer> spTransfer)
{
	UG_THROW_FATAL("No Transfer operations for Level DoF Distributions.");
}

template <>
void IDDGridFunction<LevelDoFDistribution>::clear_transfers()
{
	UG_THROW_FATAL("No Transfer operations for Level DoF Distributions.");
}


} // end namespace ug
