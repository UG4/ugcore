/*
 * global_layout.h
 *
 *  Created on: 02.08.2011
 *      Author: mrupp
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
