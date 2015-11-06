/*
 * Copyright (c) 2014:  G-CSC, Goethe University Frankfurt
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


#ifndef __H__LIB_ALGEBRA__PARALLELIZATION__PARALLELIZATION_UTIL_IMPL__
#define __H__LIB_ALGEBRA__PARALLELIZATION__PARALLELIZATION_UTIL_IMPL__

#include <vector>
#include "parallelization_util.h"

namespace ug{

template <class TIndVec>
void GenerateGlobalConsecutiveIndices(TIndVec& indsOut, size_t numLocalInds,
									  const AlgebraLayouts& layouts)
{
	typedef typename TIndVec::value_type index_t;

	PU_PROFILE_FUNC();
	indsOut.clear();
	indsOut.resize(numLocalInds, 0);

//	count local unique indices. Ignore vmasters and hslaves.
//	we're using uint64 here to stay compatible with mpi
	uint64 numLocalUniqueInds = 0;
	const index_t invalidIndex = -1;
	SetLayoutValues(&indsOut, layouts.slave(), invalidIndex);
	SetLayoutValues(&indsOut, layouts.vertical_master(), invalidIndex);

	for(size_t i = 0; i < numLocalInds; ++i){
		if(indsOut[i] != invalidIndex){
			++numLocalUniqueInds;
		}
	}

//	we'll now gather the numbers of local unique indices on proces 0 and create
//	an offset map for each process. This offset map will then be scattered to
//	the separate processes again.
	const pcl::ProcessCommunicator& procCom = layouts.proc_comm();
	std::vector<uint64>	p0_uniqueIndNumbers(procCom.size());

	procCom.gather(&numLocalUniqueInds, 1, PCL_DT_UNSIGNED_LONG_LONG,
					&p0_uniqueIndNumbers.front(), 1, PCL_DT_UNSIGNED_LONG_LONG, 0);

	if(procCom.get_local_proc_id() == 0){
	//	convert the number of unique indices of each process to offsets
		uint64 offset = 0;
		for(size_t i = 0; i < p0_uniqueIndNumbers.size(); ++i){
			size_t newOffset = offset + p0_uniqueIndNumbers[i];
			p0_uniqueIndNumbers[i] = offset;
			offset = newOffset;
		}
	}

	uint64 indOffset = 0;
	procCom.scatter(&p0_uniqueIndNumbers.front(), 1, PCL_DT_UNSIGNED_LONG_LONG,
					&indOffset, 1, PCL_DT_UNSIGNED_LONG_LONG, 0);

//	using the local offset we can now assign the local unique indices to indsOut
	index_t newInd = static_cast<index_t>(indOffset);
	for(size_t i = 0; i < numLocalInds; ++i){
		if(indsOut[i] != invalidIndex){
			indsOut[i] = newInd;
			++newInd;
		}
	}

//	finally we have to copy the unique indices from hmasters to hslaves and
//	afterwards from vslaves to vmasters.
	pcl::InterfaceCommunicator<IndexLayout>& icom = layouts.comm();
	ComPol_VecCopy<TIndVec>	compolCopy(&indsOut);
	
	icom.send_data(layouts.master(), compolCopy);
	icom.receive_data(layouts.slave(), compolCopy);
	icom.communicate();

//	Todo:	This method is mainly used on surface-layouts.
//			Currently a bug in those layouts prevents the execution of the following code
//			which would be required for level-matrices.
//			Since ghosts (pure v-masters) are ignored in surface-vectors and matrices,
//			but pure v-slaves are considered, v-master and v-slave interfaces do not
//			match. Fortunately they are not required for surface-index-generation.
//			The following block should however be used as soon as those layouts are fixed.
	// icom.send_data(layouts.vertical_slave(), compolCopy);
	// icom.receive_data(layouts.vertical_master(), compolCopy);
	// icom.communicate();
}

}//	end of namespace

#endif	//__H__LIB_ALGEBRA__PARALLELIZATION__PARALLELIZATION_UTIL_IMPL__
