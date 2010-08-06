/*
 * subset_group.h
 *
 *  Created on: 13.07.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__COMMON__SUBSET_GROUP__
#define __H__LIB_DISCRETIZATION__COMMON__SUBSET_GROUP__

#include <vector>

namespace ug{

// Subset Group is just a group integers, representing some subsets
class SubsetGroup
{
	public:
		SubsetGroup() {clear();}

		/// adds a subset to this group
		bool add_subset(int si);

		/// removes a subset from this group
		bool remove_subset(int si);

		/// clear all subsets
		void clear() {m_vSubset.clear();}

		/// number of subsets in this group
		size_t num_subsets() const {return m_vSubset.size();}

		/// number of subsets in this group
		int operator[](size_t i) const {return m_vSubset[i];}

		/// returns true if subset is contained in this group
		bool containes_subset(int si) const;

	protected:
		std::vector<int> m_vSubset;
};

} // end namespace ug

#endif /*__H__LIB_DISCRETIZATION__COMMON__SUBSET_GROUP__ */
