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

#ifndef __H__UG__LIB_DISC__DOF_MANAGER__FUNCTION_PATTERN__
#define __H__UG__LIB_DISC__DOF_MANAGER__FUNCTION_PATTERN__

#include <vector>
#include <string>

//#include "common/common.h"
#include "common/util/string_util.h"
#include "lib_grid/tools/subset_handler_interface.h"
#include "lib_grid/tools/subset_group.h"
#include "lib_grid/algorithms/subset_dim_util.h"
#include "lib_disc/local_finite_element/local_finite_element_id.h"

namespace ug {

// predeclaration
class FunctionGroup;

/// Describes the setup of discrete functions on a SubsetHandler
/**
 * A Function Pattern is used to describe the setup for discrete functions on
 * a domain. Therefore, it has the underlying SubsetHandler of the Domain.
 * The functions (sometimes also called grid functions) are defined with respect
 * to the subsets: A function can 'live' on parts of the subsets as well as on
 * the whole domain.
 */
class FunctionPattern
{
	public:
	///	Default Constructor
		explicit FunctionPattern(ConstSmartPtr<ISubsetHandler> spSH) :
			m_bLocked(false), m_spSH(spSH)
		{clear();}

	///	sets new subsets handler
		void set_subset_handler(ConstSmartPtr<ISubsetHandler> spSH);

	/// get underlying subset handler
		[[nodiscard]] ConstSmartPtr<ISubsetHandler> subset_handler() const {return m_spSH;}

	/// add single solutions of LocalShapeFunctionSetID to the entire domain
	/**
	 * \param[in] 	vName		name(s) of single solution
	 * \param[in] 	id			Shape Function set id
	 */
		void add(const std::vector<std::string>& vName,
				 LFEID id);

	/// add single solutions of LocalShapeFunctionSetID to selected subsets
	/**
	 * \param[in] vName			name(s) of single solution (comma separated)
	 * \param[in] id			Shape Function set id
	 * \param[in] vSubset		Subsetnames as a std::vector of strings
	 */
		void add(const std::vector<std::string>& vName,
		         LFEID id,
				 const std::vector<std::string>& vSubset);

	private:
	/// add single solutions of LocalShapeFunctionSetID to selected subsets
	/**
	 * \param[in] vName			name(s) of single solution
	 * \param[in] id			Shape Function set id
	 * \param[in] ssGrp			SubsetGroup, where solution lives
	 */
		void add(const std::vector<std::string>& vName,
				 LFEID id,
				 const SubsetGroup& ssGrp);

	public:
	///	lock pattern (i.e. can not be changed then)
		void lock()	{m_bLocked = true;}

	///	returns if pattern is locked
		[[nodiscard]] bool is_locked() const {return m_bLocked;}

	/// clear all functions
		void clear()
		{
			if(is_locked()) UG_THROW("Pattern locked.");
			m_vFunction.clear();
		}

	///	number of subsets
		[[nodiscard]] int num_subsets() const {return m_spSH->num_subsets();}

	///	returns a subset group consisting of the subsets specified by their names
		SubsetGroup subset_grp_by_name(const char* names) const
		{
			return SubsetGroup(subset_handler(), TokenizeString(names));
		}

	/// returns a subset group consisting of all the subsets in the domain except for the specified ones
		SubsetGroup all_subsets_grp_except_for(const char* names) const
		{
			SubsetGroup subset_grp(subset_handler());
			subset_grp.add_all();
			subset_grp.remove(TokenizeString(names));
			return subset_grp;
		}

	/// returns the subset id
		[[nodiscard]] int subset_id_by_name(const char* name) const;

	///	dimension of subset
		[[nodiscard]] int dim_subset(int si) const {return DimensionOfSubset(*m_spSH, si);}

	///	returns the name of a subset
		[[nodiscard]] const char* subset_name(int si) const {return m_spSH->subset_info(si).name.c_str();}

	/// number of discrete functions this dof distributor handles
		[[nodiscard]] size_t num_fct() const {return m_vFunction.size();}

	/// number of discrete functions on a subset
		[[nodiscard]] size_t num_fct(int si) const
		{
			size_t num = 0;
			for(size_t fct = 0; fct < num_fct(); ++fct)
			{
				if(m_vFunction[fct].is_def_in_subset(si)) num++;
			}
			return num;
		}

	/// returns the trial space of a discrete function
	/// \{
		[[nodiscard]] const LFEID& local_finite_element_id(size_t fct) const
		{
			UG_ASSERT(fct < num_fct(), "Invalid index.");
			return m_vFunction[fct].id;
		}
		[[nodiscard]] const LFEID& lfeid(size_t fct) const {return local_finite_element_id(fct);}
	/// \}

	/// returns the name of a discrete function
		[[nodiscard]] const char* name(size_t fct) const
		{
			UG_ASSERT(fct < num_fct(), "Invalid index.");
			return m_vFunction[fct].name.c_str();
		}

	/// returns the names of a discrete function
		[[nodiscard]] std::vector<std::string> names() const
		{
			std::vector<std::string> vName(num_fct());
			for(size_t fct = 0; fct < vName.size(); ++fct)
				vName[fct] = name(fct);
			return vName;
		}

	/// returns function id of the loc_fct's function on subset si
		[[nodiscard]] size_t fct_id(size_t loc_fct, int si) const
		{
			size_t fct = 0;
			for(size_t i = 0; i < loc_fct; ++i)
			{
				if(is_def_in_subset(fct, si)) fct++;
			}
			return fct;
		}

	/// returns the function id if function with given name found in pattern, exception else
		size_t fct_id_by_name(const char* name) const;

	/// returns the dimension in which solution lives
		[[nodiscard]] int dim(size_t fct) const
		{
			UG_ASSERT(fct < num_fct(), "Invalid index.");
			return m_vFunction[fct].id.dim();
		}

	/// returns true if the discrete function nr_fct is defined on subset si
		[[nodiscard]] bool is_def_in_subset(size_t fct, int si) const
		{
			UG_ASSERT(fct < num_fct(), "Invalid index.");
			return m_vFunction[fct].is_def_in_subset(si);
		}

	/// returns true if the discrete function nr_fct is defined on all subsets
		[[nodiscard]] bool is_def_everywhere(size_t fct) const
		{
			UG_ASSERT(fct < num_fct(), "Invalid index.");
			return m_vFunction[fct].is_def_everywhere();
		}

	/// virtual destructor
		virtual ~FunctionPattern() = default;


	protected:
	///	internal structure to hold all functions
		struct Function
		{
			Function(const char* name_, LFEID id_, bool everywhere_,
			         const SubsetGroup& subsetIndices_)
				: name(name_), id(id_),
				  everywhere(everywhere_), subsetIndices(subsetIndices_){};

			std::string name;
			LFEID id;
			bool everywhere;
			SubsetGroup subsetIndices;

			[[nodiscard]] bool is_def_in_subset(int si) const
			{
				if(everywhere) return true;
				if(subsetIndices.contains(si)) return true;
				return false;
			}

			[[nodiscard]] bool is_def_everywhere() const {return everywhere;}
		};

	protected:
	// 	flag if pattern is fixed
		bool m_bLocked;

	// 	underlying subset handler
		ConstSmartPtr<ISubsetHandler> m_spSH;

	// 	informations about Functions
		std::vector<Function> m_vFunction;
};

} // end namespace ug

#endif