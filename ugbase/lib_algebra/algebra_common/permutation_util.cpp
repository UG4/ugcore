
#include "permutation_util.h"


namespace ug{


bool GetInversePermutation(const std::vector<size_t> &perm, std::vector<size_t> &invPerm)
{
	invPerm.resize(perm.size());
	bool bId = true;
	for(size_t i=0; i<perm.size(); i++) invPerm[i] = (size_t) (-1);

	for(size_t i=0; i<perm.size(); i++)
	{
		UG_COND_THROW(invPerm[perm[i]] != (size_t) (-1), "not a bijective permutation "
			"(double mapping to index " << perm[i] << " by indices " << invPerm[perm[i]] << " and " << i << ")!");
		bId = bId && perm[i] == i;
		invPerm[perm[i]] = i;
	}
	return bId;
}


}
