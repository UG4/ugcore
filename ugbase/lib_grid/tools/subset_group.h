/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
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

#ifndef __H__UG__LIB_GRID__SUBSET_GROUP__
#define __H__UG__LIB_GRID__SUBSET_GROUP__

#include <vector>
#include <string>
#include "lib_grid/tools/subset_handler_interface.h"

namespace ug{

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
		SubsetGroup();

	///	Constructor setting subset handler
		explicit SubsetGroup(ConstSmartPtr<ISubsetHandler> sh);

	///	Constructor setting subset handler and subsets
		SubsetGroup(ConstSmartPtr<ISubsetHandler> sh, const char* names);

	///	Constructor setting subset handler and subsets
		SubsetGroup(ConstSmartPtr<ISubsetHandler> sh, const std::string& names);

	///	Constructor setting subset handler and subsets
		SubsetGroup(ConstSmartPtr<ISubsetHandler> sh, const std::vector<std::string>& vName);

	/// set an underlying subset handler
		void set_subset_handler(ConstSmartPtr<ISubsetHandler> sh) {m_pSH = sh; clear();}

	/// get underlying subset handler
		ConstSmartPtr<ISubsetHandler> subset_handler() const {return m_pSH;}

	/// adds a subset by number to this group
		void add(int si);

	/// adds subset with a name to this group
		void add(const char* name);

	/// adds subset with a name to this group
		void add(const std::string& name);

	/// adds all subset with by name to this group
	/**
	 * This function adds all subset with by name to this group.
	 * \param[in]	vName	Name of Subset(s) to be added
	 */
		void add(const std::vector<std::string>& vName);

	/// adds all subsets of another subset to the group
		void add(const SubsetGroup& ssGroup);

	/// select all subsets of underlying subset handler
		void add_all();

	/// removes a subset from this group
		void remove(int si);

	/// removes subset with a given name from this group
		void remove(const char* name);

	/// removes subset with a given name from this group
		void remove(const std::string& name);

	/// removes subsets with given names from this group
	/**
	 * This function removes all subsets by name from this group
	 * \param[in]	vName	Name of Subset(s) to be removed
	 */
		void remove(const std::vector<std::string>& vName);

	/// removes all subsets of another subset from the group
		void remove(const SubsetGroup& ssGroup);

	/// clear all subsets
		void clear() {m_vSubset.clear();}

	/// returns if function group is empty
		[[nodiscard]] bool empty() const {return m_vSubset.empty();}

	/// number of subsets in this group
		[[nodiscard]] inline size_t size() const
		{
			if (!m_pSH.valid()) return 0;
			return m_vSubset.size();
		}

	/// index of the subset # i in this group
		inline int operator [] (size_t i) const
		{
			UG_ASSERT(is_init(), "No SubsetHandler set.");
			UG_ASSERT(i < size(), "requested subset does not exist.");
			return m_vSubset[i];
		}
		
	///	vector of the subset indices in the group
		[[nodiscard]] inline const std::vector<int>& index_vector() const
		{
			return m_vSubset;
		}

	///	name of subset
		[[nodiscard]] const char* name(size_t i) const;

	///	returns if a subset is a regular grid
		[[nodiscard]] bool regular_grid(size_t i) const;

	/// dimension of subset
	/**
	 * Returns the dimension of the subset. The dimension of the subset
	 * is defined as the highest dimension of geometric objects contained in
	 * the subset. This maximum is taken over all procs in parallel.
	 */
		[[nodiscard]] int dim(size_t i) const;

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
		[[nodiscard]] int get_highest_subset_dimension() const;

	/// returns true if subset is contained in this group
		[[nodiscard]] bool contains(int si) const;

	/// returns true if at least one subset of a name is contained in this group
		bool contains(const char* name) const;

	protected:
	// returns if SubsetGroup is ready for use
		[[nodiscard]] bool is_init() const {return m_pSH.valid();}

	protected:
		ConstSmartPtr<ISubsetHandler> m_pSH; ///< underlying SubsetHandler
		std::vector<int> m_vSubset; ///< selected Subset Indices
};

/**
 * Returns if dimension is the same in all subsets of the subset group
 * @param subsetGroup	subset group that is checked
 * @returns true if dimension is the same in all subsets, else false
 */
bool SameDimensionsInAllSubsets(const SubsetGroup& subsetGroup);

/**
 * Removes all subsets from the subset group that have a lower dimension than the
 * highest dimension contained in the subset group.
 * @param subsetGroup 	subset group that is modified
 */
void RemoveLowerDimSubsets(SubsetGroup& subsetGroup);

} // end namespace ug

#endif