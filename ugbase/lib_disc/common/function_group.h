/*
 * function_group.h
 *
 *  Created on: 13.07.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__COMMON__FUNCTION_GROUP__
#define __H__UG__LIB_DISC__COMMON__FUNCTION_GROUP__

#include <vector>

#include "common/common.h"
#include "lib_disc/local_finite_element/local_finite_element_id.h"
#include "lib_disc/dof_manager/function_pattern.h"

namespace ug{

// FunctionGroup is just a group of size_t, representing some functions
class FunctionGroup
{
	public:
	///	Default Constructor
		FunctionGroup() : m_pFunctionPattern(NULL) {clear();}

	///	Constructor setting function pattern
		FunctionGroup(const FunctionPattern& funcPattern)
			: m_pFunctionPattern(&funcPattern) {clear();}

	/// set underlying function pattern
		void set_function_pattern(const FunctionPattern& funcPattern)
			{m_pFunctionPattern = &funcPattern; clear();}

	/// get underlying function pattern
		const FunctionPattern* function_pattern() const
			{return m_pFunctionPattern;}

	/// adds a function by id to this group
		void add(size_t fct);

	/// adds all functions with a given name to this group
	/** adds all functions with a given name to this group
	 *
	 * \param[in]	name	Name of Function(s) to be added
	 * \return 		true	if at least one function added
	 * 				false	if no function found with this name
	 */
		void add(const char* name);

	/// adds all functions of a function group
		void add(const FunctionGroup& fctGroup);

	/// selects all subsets in the order of the underlying pattern
		void add_all();

	/// removes a function by id from this group
		void remove(size_t fct);

	/// removes all functions with a given name from this group
	/** removes all functions with a given name to this group
	 *
	 * \param[in]	name	Name of Function(s) to be removed
	 * \return 		true	if at least one function removed
	 * 				false	if no function found with this name
	 */
		void remove(const char* name);

	/// clear all subsets
		void clear() {m_vFunction.clear();}

	/// sorts the selected functions by increasing unique id
		void sort();

	/// returns if function group is empty
		bool empty() const {return m_vFunction.empty();}

	/// number of functions in this group
		size_t size() const {return num_fct();}

	/// number of functions in this group
		size_t num_fct() const {return m_vFunction.size();}

	/// returns the name of a function
		const char* name(size_t i) const;

	/// returns unique function id of a function
		size_t operator[](size_t i) const {return unique_id(i);}

	/// returns unique function id of a function
		size_t unique_id(size_t i) const
		{
			UG_ASSERT(i < num_fct(), "Invalid function access");
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
		LFEID local_finite_element_id(size_t i) const;

	/// returns true if unique id is contained in this group
		bool contains(size_t uniqueID) const;

	/// returns true if all unique ids of another function group are contained
		bool contains(const FunctionGroup& fctGroup) const;

	/// return index in Function group for a function
		size_t local_index(size_t uniqueID) const;

	protected:
	/// returns if FunctionGroup is ready for use
		bool is_init() const {return m_pFunctionPattern != NULL;}

	protected:
	/// underlying function pattern
		const FunctionPattern* m_pFunctionPattern;

	/// vector holding all selected unique function ids
		std::vector<size_t> m_vFunction;
};

/// describes a mapping between local ids of two FunctionGroups
class SubFunctionMap
{
	public:
		SubFunctionMap() {clear();}

		/** select function
		 *	The selected function gets the next free index
		 *
		 * \param[in] 	fct		index of function to be selected
		 * \return 				true is selected
		 * 						false if unable to select
		 */
		void select(size_t fct) {m_vMap.push_back(fct);}

		/// clear map
		void clear() {m_vMap.clear();}

		/// number of selected functions
		size_t num_fct() const {return m_vMap.size();}

		/// return number of function
		size_t operator[](size_t i) const
		{
			UG_ASSERT(i < num_fct(), "Invalid function access.");
			return m_vMap[i];
		}

	protected:
		// mapping SubFunction -> Function
		std::vector<size_t> m_vMap;
};

/// describes a mapping between two local index sets
/**
 * This class is used to define mapping between to index sets. The domain index
 * set must contain indices in consecutive order [0, ..., N],
 * while the codomain index set can have any size.
 */
class FunctionIndexMapping
{
	public:
	/// Constructor
		FunctionIndexMapping() {clear();}

	/// removes all connections
		void clear() {m_vMapping.clear();}

	/// adds a mapping between indexFrom and indexTo
		void add(size_t indexFrom, size_t indexTo)
		{
		//	increase size if needed, setting all new entries to invalid
			if(indexFrom >= num_fct()) m_vMapping.resize(indexFrom+1, -1);

		// 	set the mapping
			m_vMapping[indexFrom] = indexTo;
		}

	/// returns the number of indices that are mapped
		size_t num_fct() const {return m_vMapping.size();}

	/// returns the mapped index
		size_t operator[](size_t i) const
		{
			UG_ASSERT(i < num_fct(), "Invalid index.\n");
			UG_ASSERT(m_vMapping[i] >= 0, "Mapping undefined.\n");
			return m_vMapping[i];
		}

	/// returns if the mapping is valid for an index
		bool mapping_valid(size_t i) const
		{
			return (i < num_fct() && m_vMapping[i] >= 0);
		}

	protected:
	/// vector holding the mapped indices
		std::vector<int> m_vMapping;
};


inline
std::ostream& operator<< (std::ostream& outStream, const ug::FunctionIndexMapping& map)
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
std::ostream& operator<< (std::ostream& outStream, const ug::FunctionGroup& grp)
{
	outStream << '[';
	for(size_t i = 0; i < grp.num_fct(); ++i)
	{
		outStream << grp[i];
		if(i !=  grp.num_fct()-1) outStream << ',';
	}
	outStream << ']';
	return outStream;
}

} // end namespace ug

#endif /*__H__UG__LIB_DISC__COMMON__FUNCTION_GROUP__ */
