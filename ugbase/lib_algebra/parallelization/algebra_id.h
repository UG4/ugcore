/*
 * parallelization_util.h
 *
 *  Created on: 18.11.2011
 *      Author: S.Reiter
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
	bool is_slave() const { return master_proc() != pcl::GetProcRank(); }
	bool is_master() const { return master_proc() == pcl::GetProcRank(); }
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
