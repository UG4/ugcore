/*
 * Copyright (c) 2009-2021:  G-CSC, Goethe University Frankfurt
 * Author: Dmitry Logashenko
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

#ifndef __H__LIBGRID__GRID_DEBUG__
#define __H__LIBGRID__GRID_DEBUG__

#include <vector>
#include <memory>

// ug4 headers:
#include "grid/grid.h"
#include "subset_handler.h"

namespace ug
{

/// Debugging tool for function that do have no direct access to the grid
/**
 * This class provides access to grid data in functions that have no
 * direct access to the grid (e.g. get only a pointer to an element but
 * need the subsets of the corners etc.).
 * 
 * This is the class for the base Grid class.
 * 
 * REMARK: THIS IS A PURELY DEBUGGING TOOL! IT MAY NOT BE USED IN THE
 * NORMAL ROUTINES FOR THE "EVERY-DAY" USE! THIS CLASS CAN BE PATCHED AND
 * CHANGED ANY TIME, THERE NO STABLE IMPLEMENTATION MAY BE ASSUMED!
 * 
 * This class provides a global pointer (which is static in the class)
 * referencing its single object (only if it is created - otherwise the
 * pointer is nullptr).
 */
class grid_global_debug_info_provider
{
	using this_type = grid_global_debug_info_provider;
	using grid_type = Grid;
	
public:

///	creates the object (if it did not exist)
/**
 * This function creates the single object of the class.
 * If the object exists, this function throws and exception
 * (to prevent the undesired reinitialization during the debugging process).
 */
	static void create
	(
		grid_type& rGrid, ///< the grid to use
		ISubsetHandler& rSH ///< the subset handler to use
	)
	{
		if (the_object == nullptr)
			the_object.reset (new grid_global_debug_info_provider (rGrid, rSH));
		else
			UG_THROW ("Reinitialization of the grid debugging info provider is not allowed!");
	};

///	checks if an element is in a subset
	static bool elem_in_subset
	(
		GridObject * elem, ///< the element to check
		int si ///< subset index to check
	)
	{
		if (! the_object)
			return false;
		return the_object->m_pSH->get_subset_index (elem) == si;
	};
	
///	checks if an element is in subsets from a list
	static bool elem_in_subsets
	(
		GridObject * elem, ///< the element to check
		std::vector<int> si_ar ///< subset indices to check
	)
	{
		if (! the_object)
			return false;
		int si = the_object->m_pSH->get_subset_index (elem);
		for (size_t i = 0; i < si_ar.size (); i++) {if (si == si_ar[i]) {return true;}}
		return false;
	};

///	checks if one of the associated elements is in a given subset
	template <typename TAssElem>
	static bool ass_elem_in_subset
	(
		GridObject * elem, ///< the element to check
		int si ///< subset index to check
	)
	{
		if (! the_object)
			return false;
		
		typename Grid::traits<TAssElem>::secure_container ass_elem_list;
		the_object->m_pGrid->associated_elements (ass_elem_list, elem);
		for (size_t i = 0; i < ass_elem_list.size (); i++)
			if (the_object->m_pSH->get_subset_index (ass_elem_list [i]) == si)
				return true;
		return false;
	};
	
///	checks if one of the associated elements is in (all or some of the) given subsets
	template <typename TAssElem>
	static bool ass_elem_in_subsets
	(
		GridObject * elem, ///< the element to check
		std::vector<int> si_ar, ///< subset indices to check
		bool in_all = false ///< if to check all the subsets of the list
	)
	{
		if (! the_object)
			return false;
		if (si_ar.empty())
			return in_all; // dummy value ("sum or product over the empty set")

		using flags_t = unsigned long;
		
		flags_t flags = 0, all_flags = 1;
		if (in_all)
		{
			if (si_ar.size () > sizeof (flags_t))
				UG_THROW ("Grid debugging info provider: too many subsets for the in_all flag");
			all_flags = (all_flags << si_ar.size ()) - 1; // 1s at the positions of the subset indices
		}
		
		typename Grid::traits<TAssElem>::secure_container ass_elem_list;
		int si;
		the_object->m_pGrid->associated_elements (ass_elem_list, elem);
		for (size_t i = 0; i < ass_elem_list.size (); i++)
		{
			si = the_object->m_pSH->get_subset_index (ass_elem_list [i]);
			for (size_t j = 0; j < si_ar.size (); j++)
				if (si == si_ar [j])
				{
					if (! in_all) return true;
					flags |= static_cast<flags_t>(1) << j;
				}
		}
		return flags == all_flags;
	};
	
protected:
	
/// (protected) constructor
	grid_global_debug_info_provider
	(
		grid_type& rGrid,
		ISubsetHandler& rSH
	)
	: m_pGrid (&rGrid), m_pSH (&rSH)
	{
	};
	
private:
//	Data:

	static std::unique_ptr<grid_global_debug_info_provider> the_object;

	grid_type * m_pGrid; ///< current grid
	ISubsetHandler * m_pSH; ///< current SubsetHandler to use
};

} // namespace ug

#endif
