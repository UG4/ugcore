/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
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

#ifndef GLOBAL_LAYOUT_H_
#define GLOBAL_LAYOUT_H_

#ifdef UG_PARALLEL

#include <vector>
#include <map>
#include "pcl/pcl.h"
#include "parallelization.h"
#include "lib_algebra/parallelization/parallelization.h"

namespace ug
{

typedef std::map<int, std::vector<AlgebraID> > GlobalLayout;

void ReceiveGlobalLayout(pcl::InterfaceCommunicator<IndexLayout> &comm, const std::vector<int> &srcprocs,
		GlobalLayout &globalMasterLayout, GlobalLayout &globalSlaveLayout);
void SendGlobalLayout(pcl::InterfaceCommunicator<IndexLayout> &comm,
		const GlobalLayout &globalMasterLayout, const GlobalLayout &globalSlaveLayout, int pid);
void MergeGlobalLayout(GlobalLayout &globalLayout, const std::map<int, int> &merge);


template<typename TGlobalToLocal>
void AddLayoutFromGlobalLayout(IndexLayout &layout, GlobalLayout &globalLayout, TGlobalToLocal &globalToLocal, bool bRemoveDoubles=true)
{
	PROFILE_FUNC();
	for(GlobalLayout::iterator it = globalLayout.begin(); it != globalLayout.end(); ++it)
	{
		int pid = it->first;
		std::vector<AlgebraID> &v = it->second;
		if(v.size() == 0) continue;
		sort(v.begin(), v.end());

		IndexLayout::Interface &interface = layout.interface(pid);

		if(bRemoveDoubles)
		{
			interface.push_back(globalToLocal[v[0]]);
			for(size_t i=1; i<v.size(); i++)
				if(v[i-1] != v[i])
					interface.push_back(globalToLocal[v[i]]);
		}
		else
		{
			for(size_t i=0; i<v.size(); i++)
				interface.push_back(globalToLocal[v[i]]);
		}
	}
}

template<typename TGlobalToLocal>
inline void CreateLayoutFromGlobalLayout(IndexLayout &layout, GlobalLayout &globalLayout, TGlobalToLocal &globalToLocal, bool bRemoveDoubles=true)
{
	layout.clear();
	AddLayoutFromGlobalLayout(layout, globalLayout, globalToLocal, bRemoveDoubles);
}


template<typename TLocalToGlobal>
void CreateGlobalLayout(GlobalLayout &globalLayout,
                        const IndexLayout &layout, const TLocalToGlobal &localToGlobal)
{
	PROFILE_FUNC();
	for(IndexLayout::const_iterator iter = layout.begin(); iter != layout.end(); ++iter)
	{
		const IndexLayout::Interface &interface = layout.interface(iter);
		int pid = layout.proc_id(iter);

		std::vector<AlgebraID> &v = globalLayout[pid];
		for(IndexLayout::Interface::const_iterator iter2 = interface.begin(); iter2 != interface.end(); ++iter2)
		{
			size_t localIndex = interface.get_element(iter2);
			v.push_back(localToGlobal[localIndex]);
		}
	}
}

void PrintGlobalLayout(const GlobalLayout &globalLayout, const char *name=NULL);
} // namespace ug

#endif

#endif /* GLOBAL_LAYOUT_H_ */
