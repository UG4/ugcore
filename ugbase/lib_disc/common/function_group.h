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

#ifndef __H__UG__LIB_DISC__COMMON__FUNCTION_GROUP__
#define __H__UG__LIB_DISC__COMMON__FUNCTION_GROUP__

#include <vector>
#include <string>

#include "common/common.h"
#include "lib_disc/local_finite_element/local_finite_element_id.h"
#include "lib_disc/dof_manager/function_pattern.h"

namespace ug{

/** Group of functions represented by integers
 * FunctionGroup is just a group of size_t, representing some functions. The
 * function group is based on a FunctionPattern and the integer represent the
 * position of the function in the function pattern. Selection of functions is
 * best via usage of symbolic names of the functions.
 */
class FunctionGroup
{
	public:
	///	Default Constructor
		FunctionGroup();

	///	Constructor setting function pattern
		FunctionGroup(ConstSmartPtr<FunctionPattern> spFuncPattern);

	///	Constructor setting function pattern and function
		FunctionGroup(ConstSmartPtr<FunctionPattern> spFuncPattern, const char* name);

	///	Constructor setting function pattern and function
		FunctionGroup(ConstSmartPtr<FunctionPattern> spFuncPattern, const std::string& name);

	///	Constructor setting function pattern and functions
		FunctionGroup(ConstSmartPtr<FunctionPattern> spFuncPattern, const std::vector<std::string>& vName);

	/// set underlying function pattern
		void set_function_pattern(ConstSmartPtr<FunctionPattern> spFuncPattern);

	/// get underlying function pattern
		ConstSmartPtr<FunctionPattern> function_pattern() const
			{return m_spFunctionPattern;}

	/// adds a function by id to this group
		void add(size_t fct);

	/// adds function with a given name to this group
		void add(const char* name);

	/// adds function with a given name to this group
		void add(const std::string& name);

	/// adds functions with a given names to this group
		void add(const std::vector<std::string>& name);

	/// adds all functions of a function group
		void add(const FunctionGroup& fctGroup);

	/// selects all subsets in the order of the underlying pattern
		void add_all();

	/// removes a function by id from this group
		void remove(size_t fct);

	/// removes function with a given name from this group
		void remove(const char* name);

	/// removes function with a given name from this group
		void remove(const std::string& name);

	/// removes functions with a given names from this group
		void remove(const std::vector<std::string>& vName);

	/// clear all subsets
		void clear() {m_vFunction.clear();}

	/// sorts the selected functions by increasing unique id
		void sort();

	/// returns if function group is empty
		bool empty() const {return m_vFunction.empty();}

	/// number of functions in this group
		size_t size() const {return m_vFunction.size();}

	/// returns the name of a function
		const char* name(size_t i) const;

	/// returns the comma-separted names of all functions
		std::string names() const;

	/// returns unique function id of a function
		size_t operator [] (size_t i) const {return unique_id(i);}

	/// returns unique function id of a function
		size_t unique_id(size_t i) const
		{
			UG_ASSERT(i < size(), "Invalid function access");
			return m_vFunction[i];
		}

	/// dimension of a function
	/**
	 * Returns the dimension of a function. The dimension is defined
	 * to be the highest dimension of grid entity the function lives on
	 */
		int dim(size_t i) const;

	/// common dimension of all functions
	/**
	 *	Returns the commen dimension of all functions of the group. This
	 *	dimension is defined to be the highest dimension of grid entity
	 *	the functions live on
	 *
	 * \return 		-1			if no common dimension available
	 * 				dim	>= 0	common dimension
	 */
		int dim() const;

	/// returns the trial space of the discrete function fct
	/// \{
		LFEID local_finite_element_id(size_t i) const;
		LFEID lfeid(size_t i) const;
	/// \}

	/// returns true if unique id is contained in this group
		bool contains(size_t uniqueID) const;

	/// returns true if all unique ids of another function group are contained
		bool contains(const FunctionGroup& fctGroup) const;

	/// return index in Function group for a function
		size_t local_index(size_t uniqueID) const;

	protected:
	/// returns if FunctionGroup is ready for use
		bool is_init() const {return m_spFunctionPattern.valid();}

	protected:
	/// underlying function pattern
		ConstSmartPtr<FunctionPattern> m_spFunctionPattern;

	/// vector holding all selected unique function ids
		std::vector<size_t> m_vFunction;
};

/// describes a mapping between two local index sets
/**
 * This class is used to define a mapping between two index sets. The domain index
 * set must contain indices in consecutive order [0, ..., N],
 * while the codomain index set can have any size.
 */
class FunctionIndexMapping
{
	public:
	/// removes all connections
		void clear() {m_vMapping.clear();}

	/// adds a mapping between indexFrom and indexTo
		void push_back(size_t indexTo){m_vMapping.push_back(indexTo);}

	/// returns the number of indices that are mapped
		size_t num_fct() const {return m_vMapping.size();}

	/// returns the mapped index
		size_t operator [] (size_t i) const
		{
			UG_ASSERT(i < num_fct(), "Invalid index.\n");
			return m_vMapping[i];
		}

	protected:
	/// vector holding the mapped indices
		std::vector<size_t> m_vMapping;
};


inline
std::ostream& operator << (std::ostream& outStream, const FunctionIndexMapping& map)
{
	outStream << '[';
	for(size_t i = 0; i < map.num_fct(); ++i)
	{
		outStream << map[i];
		if(i !=  map.num_fct()-1) outStream << ',';
	}
	outStream << ']';
	return outStream;
}

inline
std::ostream& operator << (std::ostream& outStream, const FunctionGroup& grp)
{
	outStream << '[';
	for(size_t i = 0; i < grp.size(); ++i)
	{
		outStream << grp[i];
		if(i !=  grp.size()-1) outStream << ',';
	}
	outStream << ']';
	return outStream;
}

} // end namespace ug

#endif