/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
 * Author: Markus Breit
 *
 * This file is part of UG4.
 *
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 *
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 *
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 *
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#include "mem_info.h"

#ifdef UG_PARALLEL
#include "pcl/pcl_base.h"  // for NumProcs
#include "pcl/pcl_process_communicator.h"  // for NumProcs
#endif

namespace ug {


void MemInfo::communicate_process_values()
{
#ifdef UG_PARALLEL
	if (pcl::NumProcs() > 1)
	{
		pcl::ProcessCommunicator pc;
		pc.allreduce(&m_locRes, &m_gloRes, 1, PCL_RO_SUM);
		pc.allreduce(&m_locVirt, &m_gloVirt, 1, PCL_RO_SUM);
		pc.allreduce(&m_locRes, &m_maxRes, 1, PCL_RO_MAX);
		pc.allreduce(&m_locVirt, &m_maxVirt, 1, PCL_RO_MAX);
		return;
	}
#endif

	m_gloRes = m_locRes;
	m_gloVirt = m_locVirt;
	m_maxRes = m_locRes;
	m_maxVirt = m_locVirt;
}

number MemInfo::local_resident_memory() const
{
	return m_locRes;
}

number MemInfo::local_virtual_memory() const
{
	return m_locVirt;
}

number MemInfo::global_resident_memory() const
{
	return m_gloRes;
}

number MemInfo::global_virtual_memory() const
{
	return m_gloVirt;
}

number MemInfo::max_resident_memory() const
{
	return m_maxRes;
}

number MemInfo::max_virtual_memory() const
{
	return m_maxVirt;
}


} // namespace ug
