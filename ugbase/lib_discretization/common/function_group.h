/*
 * function_group.h
 *
 *  Created on: 13.07.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__COMMON__FUNCTION_GROUP__
#define __H__LIB_DISCRETIZATION__COMMON__FUNCTION_GROUP__

#include <vector>

#include "common/common.h"
#include "lib_discretization/dof_manager/function_pattern.h"

namespace ug{

////////////////////////////////////////////////////////////////////////
//	ERROR_BadIndexInFunctionGroup
struct ERROR_BadIndexInFunctionGroup{
	ERROR_BadIndexInFunctionGroup(int index) : m_index(index)	{}
	int m_index;
};

// FunctionGroup is just a group of size_t, representing some functions
class FunctionGroup
{
	public:
		FunctionGroup() : m_pFunctionPattern(NULL) {clear();}

		/// set underlying function pattern
		void set_function_pattern(const FunctionPattern& funcPattern) {m_pFunctionPattern = &funcPattern; clear();}

		/// get underlying function pattern
		const FunctionPattern* get_function_pattern() const {return m_pFunctionPattern;}

		/// adds a function by id to this group
		bool add_function(size_t fct);

		/// adds all functions with a given name to this group
		/** adds all functions with a given name to this group
		 *
		 * \param[in]	name	Name of Function(s) to be added
		 * \return 		true	if at least one function added
		 * 				false	if no function found with this name
		 */
		bool add_function(const char* name);

		/// removes a function by id from this group
		bool remove_function(size_t fct);

		/// removes all functions with a given name from this group
		/** removes all functions with a given name to this group
		 *
		 * \param[in]	name	Name of Function(s) to be removed
		 * \return 		true	if at least one function removed
		 * 				false	if no function found with this name
		 */
		bool remove_function(const char* name);

		/// selects all subsets in the order of the underlying pattern
		void add_all_functions();

		/// clear all subsets
		void clear() {m_vFunction.clear();}

		/// returns if function group is empty
		bool empty() {return m_vFunction.empty();}

		/// number of functions in this group
		size_t num_fct() const
		{
			UG_ASSERT(is_init(), "No FunctionPattern set.");
			return m_vFunction.size();
		}

		/// get function id
		size_t operator[](size_t i) const
		{
			UG_ASSERT(is_init(), "No FunctionPattern set.");
			return fct_id(i);
		}

		/// name of function in this group
		const char* get_function_name(size_t i) const;

		/// get function id
		size_t fct_id(size_t i) const
		{
			UG_ASSERT(is_init(), "No FunctionPattern set.");
			UG_ASSERT(i < num_fct(), "Invalid function access");
			return m_vFunction[i];
		}

		/// dimension of function (i.e. highest dimension of grid entity the function lives on)
		int get_function_dimension(size_t i) const;

		/// common dimension of all functions (i.e. highest dimension of grid entity the functions live on)
		/** common dimension of all functions
		 *
		 * \return 		-1			if no common dimension available
		 * 				dim	>= 0	common dimension of all functions in this function group
		 */
		int get_function_dimension() const;

		/// returns the trial space of the discrete function fct
		LocalShapeFunctionSetID local_shape_function_set_id(size_t i) const;

		/// returns true if subset is contained in this group
		bool containes_function(size_t fct) const;

		/// return index in Function group for a function
		size_t local_index(size_t fct) const;

	protected:
		// returns if FunctionGroup is ready for use
		bool is_init() const {return m_pFunctionPattern != NULL;}

	protected:
		const FunctionPattern* m_pFunctionPattern;

		std::vector<size_t> m_vFunction;
};

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


} // end namespace ug

#endif /*__H__LIB_DISCRETIZATION__COMMON__FUNCTION_GROUP__ */
