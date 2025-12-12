/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
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

#ifndef __H__UG__partition_map__
#define __H__UG__partition_map__

#include "lib_grid/file_io/file_io.h"
#include "lib_grid/grid/grid.h"
#include "lib_grid/tools/subset_handler_grid.h"
#include "common/util/smart_pointer.h"

namespace ug {


/** \ingroup lib_grid_tools
 *  \{ */

///	Used to describe how a domain shall be distributed in a parallel environment.
/**	A partition map holds a subset handler and a map, which specifies the
 * target process of each subset. Make sure to assign a grid before performing
 * the partitioning.
 *
 * \todo	a partition map should feature a constructor which takes a grid.
 */
class PartitionMap{
	public:
		PartitionMap();

		void clear();

		void assign_grid(Grid& grid);

		SmartPtr<SubsetHandler> get_partition_handler();

		void add_target_proc(int tarProcRank);

		void add_target_procs(int first, int num);

		[[nodiscard]] size_t num_target_procs() const;

		[[nodiscard]] int get_target_proc(size_t index) const;

		int* get_target_procs();

		std::vector<int>& get_target_proc_vec();

	///	changes an existing target process. Make sure that index < num_target_procs
		bool change_target_proc(size_t index, int newRank);

	///	returns the index at which the given process lies. -1 if it doesn't exist.
		[[nodiscard]] int find_target_proc(int procRank) const;

	///	adds the given offset to all target-proc-ranks
		void shift_target_procs(int offset);

	private:
		SmartPtr<SubsetHandler>	m_shPartitions;
		std::vector<int>		m_targetProcs;
};

using SPPartitionMap = SmartPtr<PartitionMap>;



///	Save the partition map to a file.
/**	The resulting file will contain the grid on which the partition-map operates,
 * together with subsets, each representing the process on which the subset will
 * be sent.
 *
 * \todo	currently only the .ugx format is supported.
 */
template <typename TAPos>
bool SavePartitionMapToFile(PartitionMap& pm, const char* filename,
							TAPos& aPos)
{
	SubsetHandler& partsh = *pm.get_partition_handler();

//	make sure that a grid exists
	if(!partsh.grid()){
		UG_LOG("WARNING IN SavePartitionMapToFile: a grid has to be assigned "
				"to the PartitionMap. Aborting.\n");
		return false;
	}

	Grid& grid = *partsh.grid();

//	we need a subset-handler, which will have a 1-1 subset-process relation.
	SubsetHandler sh(grid);

//	add all partitions to the handler
	for(int si = 0; si < partsh.num_subsets(); ++si){
		int newSI = pm.get_target_proc(si);
		sh.assign_subset(partsh.begin<Vertex>(si),
						 partsh.end<Vertex>(si), newSI);
		sh.assign_subset(partsh.begin<Edge>(si),
						 partsh.end<Edge>(si), newSI);
		sh.assign_subset(partsh.begin<Face>(si),
						 partsh.end<Face>(si), newSI);
		sh.assign_subset(partsh.begin<Volume>(si),
						 partsh.end<Volume>(si), newSI);
	}

//	now save the grid to file
	return SaveGridToFile(grid, sh, filename, aPos);
}

/** \} */

}//	end of namespace

#endif
