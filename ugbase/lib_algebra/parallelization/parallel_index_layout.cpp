/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
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

#include "common/common.h"
#include <vector>
#include <iomanip>
#include "pcl/pcl.h"
#include "parallel_index_layout.h"
#include "common/util/string_util.h"

namespace ug{

void LogIndexLayout(IndexLayout& layout, int depth)
{
	using namespace std;

	using Interface = IndexLayout::Interface;
	using InterfaceIter = IndexLayout::iterator;

	UG_LOG("-- IndexLayout Informations: Proc "<< GetLogAssistant().get_output_process() << " --\n");

	UG_LOG(" interface | target proc id |   size    ");
	if(depth >= 1) UG_LOG(" | indices ")
	UG_LOG("\n");

	int i = 0;
	for(auto iiter = layout.begin(); iiter != layout.end(); ++iiter, ++i)
	{
		Interface& interface = layout.interface(iiter);
		UG_LOG(" " << std::setw(9) << i << " | " << std::setw(14) <<
		       layout.proc_id(iiter) << " | " << std::setw(9) << interface.size() << " ");
		if(depth >= 1)
		{
			UG_LOG(" | (");
			for(auto indexIter = interface.begin(); indexIter != interface.end(); ++indexIter)
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


std::ostream &operator << (std::ostream &out, const IndexLayout &layout)
{
	out << "IndexLayout: ";
	for(auto iter = layout.begin(); iter != layout.end(); ++iter)
	{
		size_t pid = layout.proc_id(iter);
		const IndexLayout::Interface &interface = layout.interface(iter);
		out << "\n to processor " << pid << " (size " << interface.size() << ") :";
		std::stringstream ss;
		size_t k=0;
		int ndigit= NumberOfDigits(interface.size());
		for(auto iter2 = interface.begin(); iter2 != interface.end(); ++iter2)
		{
			if(k % 10 == 0)
			{
				if(k!=0) ss << "\n";
				ss << std::setw(ndigit) << k << "-" << std::setw(ndigit) << std::min(interface.size(), k+9) << ": ";
			}
			ss << std::setw(5) << interface.get_element(iter2) << "  ";
			k++;
		}
		out << ConfigShift(ss.str());
	}

	return out;
}

void LogIndexLayoutOnAllProcs(IndexLayout& layout, int depth)
{
//	remember current outproc
	int outproc = GetLogAssistant().get_output_process();

//	loop all procs
	for(int p = 0; p < pcl::NumProcs(); ++p)
	{
	//	synchronize, to prevent other procs to write before this one has finished.
		pcl::SynchronizeProcesses();

	//	write process p
		if(p == pcl::ProcRank())
		{
		//	set output proc to proc p
			GetLogAssistant().set_output_process(p);

		//	write
			LogIndexLayout(layout, depth);
		}
	}
	pcl::SynchronizeProcesses();
	UG_LOG(std::flush);

//	reset output proc
	GetLogAssistant().set_output_process(outproc);
}

void ReplaceIndicesInLayout(IndexLayout& layout, const std::vector<int>& vMap)
{
	PROFILE_FUNC_GROUP("algebra parallelization");
//	interface iterators
	auto interfaceIter = layout.begin();
	auto interfaceEnd = layout.end();

//	iterate over interfaces
	for(; interfaceIter != interfaceEnd; ++interfaceIter)
	{
	//	get interface
		IndexLayout::Interface& interface = layout.interface(interfaceIter);

	//	loop over indices
		for(auto iter = interface.begin(); iter != interface.end();)
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

void MarkAllFromInterface(std::vector<bool> &mark, const IndexLayout::Interface &interface)
{
	for(auto iter = interface.begin(); iter != interface.end(); ++iter)
		mark[ interface.get_element(iter) ] = true;
}

void MarkAllFromLayout(std::vector<bool> &mark, const IndexLayout &layout)
{
	for(auto iter = layout.begin(); iter != layout.end(); ++iter)
		MarkAllFromInterface(mark, layout.interface(iter));
}


void AddAllFromInterface(std::set<size_t> &s, const IndexLayout::Interface &interface)
{
	for(auto iter = interface.begin(); iter != interface.end(); ++iter)
		s.insert(interface.get_element(iter));
}

void AddAllFromLayout(std::set<size_t> &s, const IndexLayout &layout)
{
	for(auto iter = layout.begin(); iter != layout.end(); ++iter)
		AddAllFromInterface(s, layout.interface(iter));
}

}
