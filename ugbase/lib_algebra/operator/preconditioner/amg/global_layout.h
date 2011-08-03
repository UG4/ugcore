/*
 * global_layout.h
 *
 *  Created on: 02.08.2011
 *      Author: mrupp
 */

#ifndef GLOBAL_LAYOUT_H_
#define GLOBAL_LAYOUT_H_

#include <vector>
#ifdef UG_PARALLEL
#include "pcl/pcl.h"
#include "parallelization.h"
#include "lib_algebra/parallelization/parallelization.h"
#endif


namespace ug
{

typedef std::map<int, std::vector<AlgebraID> > GlobalLayout;

void ReceiveGlobalLayout(pcl::ParallelCommunicator<IndexLayout> &comm, std::vector<int> &srcprocs, GlobalLayout &globalMasterLayout, GlobalLayout &globalSlaveLayout);
void SendGlobalLayout(pcl::ParallelCommunicator<IndexLayout> &comm, GlobalLayout &globalMasterLayout, GlobalLayout &globalSlaveLayout, int pid);
void MergeGlobalLayout(GlobalLayout &globalLayout, std::map<int, int> &merge);

template<typename TGlobalToLocal>
inline void CreateLayoutFromGlobalLayout(IndexLayout &layout, GlobalLayout &globalLayout, TGlobalToLocal &globalToLocal)
{
	layout.clear();
	for(GlobalLayout::iterator it = globalLayout.begin(); it != globalLayout.end(); ++it)
	{
		int pid = it->first;
		std::vector<AlgebraID> &v = it->second;
		sort(v.begin(), v.end());

		IndexLayout::Interface &interface = layout.interface(pid);
		for(size_t i=0; i<v.size(); i++)
			interface.push_back(globalToLocal[v[i]]);
	}
}

template<typename TLocalToGlobal>
inline void CreateGlobalLayout(GlobalLayout &globalLayout, IndexLayout layout, const TLocalToGlobal &localToGlobal)
{
	for(IndexLayout::iterator iter = layout.begin(); iter != layout.end(); ++iter)
	{
		IndexLayout::Interface &interface = layout.interface(iter);
		int pid = layout.proc_id(iter);

		std::vector<AlgebraID> &v = globalLayout[pid];
		for(IndexLayout::Interface::iterator iter2 = interface.begin(); iter2 != interface.end(); ++iter2)
		{
			size_t localIndex = interface.get_element(iter2);
			v.push_back(localToGlobal[localIndex]);
		}
	}

}




}
#endif /* GLOBAL_LAYOUT_H_ */
