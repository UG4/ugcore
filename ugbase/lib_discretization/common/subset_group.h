/*
 * subset_group.h
 *
 *  Created on: 13.07.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__COMMON__SUBSET_GROUP__
#define __H__LIB_DISCRETIZATION__COMMON__SUBSET_GROUP__

#include <vector>
#include "lib_grid/tools/subset_handler_interface.h"

namespace ug{

////////////////////////////////////////////////////////////////////////
//	ERROR_BadIndexInSubsetGroup
struct ERROR_BadIndexInSubsetGroup{
	ERROR_BadIndexInSubsetGroup(int subsetIndex) : m_subsetIndex(subsetIndex)	{}
	int m_subsetIndex;
};

////////////////////////////////////////////////////////////////////////
//	ERROR_SubsetGroupHasNoSubsetHandler
struct ERROR_SubsetGroupHasNoSubsetHandler{};


// Subset Group is just a group integers, representing some subsets
class SubsetGroup
{
	public:
		SubsetGroup() : m_pSH(NULL) {clear();}

		/// set an underlying subset handler
		void set_subset_handler(const ISubsetHandler& sh) {m_pSH = &sh; clear();}

		/// get underlying subset handler
		const ISubsetHandler* get_subset_handler() const {return m_pSH;}

		/// adds a subset by number to this group
		bool add_subset(int si);

		/// adds all subset with a given name to this group
		/** adds all subset with a given name to this group
		 *
		 * \param[in]	name	Name of Subset(s) to be added
		 * \return 		true	if at least one subset added
		 * 				false	if no subset found with this name
		 */
		bool add_subset(const char* name);

		/// removes a subset from this group
		bool remove_subset(int si);

		/// removes all subset with a given name from this group
		/** removes all subset with a given name to this group
		 *
		 * \param[in]	name	Name of Subset(s) to be removed
		 * \return 		true	if at least one subset removed
		 * 				false	if no subset found with this name
		 */
		bool remove_subset(const char* name);

		/// select all subsets of underlying subset
		void add_all_subsets();

		/// clear all subsets
		void clear() {m_vSubset.clear();}

		/// number of subsets in this group
		inline size_t num_subsets() const
		{
			UG_ASSERT(is_init(), "No SubsetHandler set.");
			return m_vSubset.size();
		}

		/// subset i in this group
		inline int operator[](size_t i) const
		{
			UG_ASSERT(is_init(), "No SubsetHandler set.");
			UG_ASSERT(i < num_subsets(), "requested subset does not exist.");
			return m_vSubset[i];
		}

		///	name of subset
		const char* get_subset_name(size_t i) const;

		/// dimension of subset (i.e. highest dimension of grid entity in the subset)
		int get_subset_dimension(size_t i) const;

		/// common dimension of all subset (i.e. highest dimension of grid entity in the subset)
		/** common dimension of all subset
		 *
		 * \return 		-1			if no common dimension available
		 * 				dim	>= 0	common dimension of all subsets in this subset group
		 */
		int get_subset_dimension() const;

		/// returns true if subset is contained in this group
		bool containes_subset(int si) const;

		/// returns true if at least one subset with the given name is contained in this group
		bool containes_subset(const char* name) const;

	protected:
		// returns if SubsetGroup is ready for use
		bool is_init() const {return m_pSH != NULL;}

	protected:
		const ISubsetHandler* m_pSH;

		std::vector<int> m_vSubset;
};

} // end namespace ug

#endif /*__H__LIB_DISCRETIZATION__COMMON__SUBSET_GROUP__ */
