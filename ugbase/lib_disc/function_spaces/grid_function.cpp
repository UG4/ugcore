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


void IGridFunction::add_transfer(SmartPtr<ILocalTransfer> spTransfer)
{
	m_spDD->add_transfer(spTransfer);
}

void IGridFunction::remove_transfer(SmartPtr<ILocalTransfer> spTransfer)
{
	m_spDD->remove_transfer(spTransfer);
}

void IGridFunction::clear_transfers()
{
	m_spDD->clear_transfers();
}

} // end namespace ug
