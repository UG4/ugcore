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

#ifndef __H__LIB_ALGEBRA__PARALLELIZATION__ALGEBRA_ID__
#define __H__LIB_ALGEBRA__PARALLELIZATION__ALGEBRA_ID__

#include <utility>
#include <vector>
#include <map>
#include "common/util/hash.h"
#include "pcl/pcl.h"

namespace ug{

///	this type is used to identify distributed objects.
struct AlgebraID : public std::pair<int, size_t>
{
	AlgebraID() { first = -1; second = -1; }
	AlgebraID(int _masterProc, size_t _indexOnMaster)
	{
		first = _masterProc;
		second = _indexOnMaster;
	}

	int master_proc() const { return first; }
	size_t index_on_master() const { return second; }
	bool is_slave() const { return master_proc() != pcl::ProcRank(); }
	bool is_master() const { return master_proc() == pcl::ProcRank(); }
};

template<>
size_t hash_key<AlgebraID>(const AlgebraID& key);

typedef std::vector<AlgebraID>	AlgebraIDVec;
typedef Hash<AlgebraID, size_t>	AlgebraIDHashList;

std::ostream& operator<<(std::ostream &out, const AlgebraID &ID);

///	Creates a hash which allows a algebraID->localIndex mapping
inline void GenerateAlgebraIDHashList(AlgebraIDHashList &hash, AlgebraIDVec& algebraIDs)
{
	PROFILE_FUNC_GROUP("algebra");
//	clear and resize the hash
//	We'll resize the hash to the number of algebraIDs, in order to have a
//	hopefully good mapping...
	hash.clear();
	hash.resize_hash(algebraIDs.size());
	hash.reserve(algebraIDs.size());

//	now add all ids. We use the algebraID as key and the local index as value
	for(size_t i = 0; i < algebraIDs.size(); ++i)
		hash.insert(algebraIDs[i], i);
}

} // end namespace ug

#endif /* __H__LIB_ALGEBRA__PARALLELIZATION__ALGEBRA_ID__ */
