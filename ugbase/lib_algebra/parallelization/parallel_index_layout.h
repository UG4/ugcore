/*
 * parallel_index_layout.h
 *
 *  Created on: 21.5.2010
 *      Author: A. Vogel, S.Reiter
 */

#ifndef __H__LIB_ALGEBRA__PARALLELIZATION__PARALLEL_INDEX_LAYOUT__
#define __H__LIB_ALGEBRA__PARALLELIZATION__PARALLEL_INDEX_LAYOUT__

#include <vector>
#include <iomanip>
#include "pcl/pcl.h"

namespace ug
{

///\ingroup lib_algebra_parallelization

///	Allows communication between distributed vectors and matrices.
/**	Note that indices are stored in an std::vector in the moment.
 *	This allows fast iteration and memory allocation, if dynamic
 *	interfaces are required this may however be slower than a
 *	std::list container.
 */
typedef pcl::SingleLevelLayout<pcl::OrderedInterface<size_t, std::vector> >
		IndexLayout;

///	Logs the internals of an index layout.
/**
 * Writes information about an index interface. If depth >= 1 is passed, then
 * also the current indices in the interfaces are printed.
 */
inline void LogIndexLayout(IndexLayout& layout, int depth = 0)
{
	using namespace std;

	typedef IndexLayout::Interface Interface;
	typedef IndexLayout::iterator  InterfaceIter;

	UG_LOG("-- IndexLayout Informations: Proc "<< pcl::GetOutputProcRank() << " --\n");

	UG_LOG(" interface | target proc id |   size    ");
	if(depth >= 1) UG_LOG(" | indices ")
	UG_LOG("\n");

	int i = 0;
	for(InterfaceIter iiter = layout.begin();
		iiter != layout.end(); ++iiter, ++i)
	{
		Interface& interface = layout.interface(iiter);
		UG_LOG(" " << std::setw(9) << i << " | " << std::setw(14) <<
		       layout.proc_id(iiter) << " | " << std::setw(9) << interface.size() << " ");
		if(depth >= 1)
		{
			UG_LOG(" | (");
			for(Interface::iterator indexIter = interface.begin();
					indexIter != interface.end(); ++indexIter)
			{
			//  get index
				const size_t index = interface.get_element(indexIter);

			//	add comma
				if(indexIter != interface.begin())
					UG_LOG(", ");

			//	log index
				UG_LOG(index);
			}
			UG_LOG(")");
		}
		UG_LOG("\n");
	}
	UG_LOG(endl);
}

/// logs index infos for all procs successively
inline void LogIndexLayoutOnAllProcs(IndexLayout& layout, int depth = 0)
{
//	remember current outproc
	int outproc = pcl::GetOutputProcRank();

//	loop all procs
	for(int p = 0; p < pcl::GetNumProcesses(); ++p)
	{
	//	synchronize, to prevent other procs to write before this one has finished.
		pcl::SynchronizeProcesses();

	//	write process p
		if(p == pcl::GetProcRank())
		{
		//	set output proc to proc p
			pcl::SetOutputProcRank(p);

		//	write
			LogIndexLayout(layout, depth);
		}
	}
	pcl::SynchronizeProcesses();
	UG_LOG(std::flush);

//	reset output proc
	pcl::SetOutputProcRank(outproc);
}

/// replaces the indices in the layout based on a passed mapping
/**
 * This function replaces the indices in the layout by new indices. The mapping
 * between old and new indices is passed as newIndex = vMap[oldIndex]. The size
 * of the new index set must be smaller or equal to the old index set. If the
 * entry in the map is negative, the old index is removed from the layout.
 */
inline void ReplaceIndicesInLayout(IndexLayout& layout, const std::vector<int>& vMap)
{
//	interface iterators
	IndexLayout::iterator interfaceIter = layout.begin();
	IndexLayout::iterator interfaceEnd = layout.end();

//	iterate over interfaces
	for(; interfaceIter != interfaceEnd; ++interfaceIter)
	{
	//	get interface
		IndexLayout::Interface& interface = layout.interface(interfaceIter);

	//	loop over indices
		for(IndexLayout::Interface::iterator iter = interface.begin(); iter != interface.end();)
		{
		//  get index
			size_t& index = interface.get_element(iter);

		//	get new index
			const int newIndex = vMap[index];

		//	erase index if negative
			if(newIndex < 0)
				iter = interface.erase(iter);
		//	else replace index
			else
			{
				index = newIndex;
				 ++iter;
			}

		}
	}
}

} // end namespace ug

#endif /* __H__LIB_ALGEBRA__PARALLELIZATION__PARALLEL_INDEX_LAYOUT__ */
