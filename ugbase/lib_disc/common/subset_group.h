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
#ifdef UG_PARALLEL
#include "pcl/pcl_process_communicator.h"
#endif

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


/// Group of subsets
/**
 * A SubsetGroup is used to describe a group of Subsets. Therefore, it has an
 * underlying SubsetHandler from which the subsets can be chosen. It is very
 * light-weight, since internally it is just a vector of the subset indices. But
 * it provides several comfort functions. Note, that the subset indices are
 * always sorted by increasing index in this subset group.
 */
class SubsetGroup
{
	public:
	///	Default Constructor
		SubsetGroup() : m_pSH(NULL) {clear();}

	/// set an underlying subset handler
		void set_subset_handler(const ISubsetHandler& sh) {m_pSH = &sh; clear();}

	/// get underlying subset handler
		const ISubsetHandler* get_subset_handler() const {return m_pSH;}

	/// adds a subset by number to this group
		bool add(int si);

	/// adds all subset with by name to this group
	/**
	 * This function adds all subset with by name to this group.
	 *
	 * \param[in]	name	Name of Subset(s) to be added
	 * \return 		true	if at least one subset added
	 * 				false	if no subset found with this name
	 */
		bool add(const char* name);

	/// adds all subsets of another subset to the group
		bool add(const SubsetGroup& ssGroup);

	/// select all subsets of underlying subset
		bool add_all();

	/// removes a subset from this group
		bool remove(int si);

	/// removes all subset with a given name from this group
	/**
	 * This function removes all subsets by name from this group
	 *
	 * \param[in]	name	Name of Subset(s) to be removed
	 * \return 		true	if at least one subset removed
	 * 				false	if no subset found with this name
	 */
		bool remove(const char* name);

	/// removes all subsets of another subset from the group
		bool remove(const SubsetGroup& ssGroup);

	/// clear all subsets
		void clear() {m_vSubset.clear();}

	/// returns if function group is empty
		bool empty() {return m_vSubset.empty();}

	/// number of subsets in this group
		inline size_t num_subsets() const
		{
			if (m_pSH==NULL) return 0;
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
		const char* name(size_t i) const;

	///	returns if a subset is a regular grid
		bool regular_grid(size_t i) const;

	/// dimension of subset
	/**
	 * Returns the dimension of the subset. The dimension of the subset
	 * is defined as the highest dimension of geometric objects contained in
	 * the subset.
	 * If a ProcessCommunicator is passed, the maximum among all procs
	 * is returned.
	 */
		int dim(size_t i
#ifdef UG_PARALLEL
                     ,pcl::ProcessCommunicator* pProcCom = NULL
#endif
				) const;

	/// highest dimension of all subset
	/**
	 * Returns the highest dimension of all subset. The dimension of a subset
	 * is defined as the highest dimension of geometric objects contained in
	 * the subset.
	 * No check between different processes is performed in parallel.
	 *
	 * \return 		-1			if no dimension available
	 * 				dim	>= 0	highest dimension of all subsets in this group
	 */
		int get_local_highest_subset_dimension() const;

	/// returns true if subset is contained in this group
		bool contains(int si) const;

	/// returns true if at least one subset of a name is contained in this group
		bool contains(const char* name) const;

	protected:
	// returns if SubsetGroup is ready for use
		bool is_init() const {return m_pSH != NULL;}

	protected:
		const ISubsetHandler* m_pSH; ///< underlying SubsetHandler
		std::vector<int> m_vSubset; ///< selected Subset Indices
};

} // end namespace ug

#endif /*__H__LIB_DISCRETIZATION__COMMON__SUBSET_GROUP__ */
