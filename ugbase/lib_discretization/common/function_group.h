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


namespace ug{

// FunctionGroup is just a group of size_t, representing some functions
class FunctionGroup
{
	public:
		FunctionGroup() {clear();}

		/// adds a subset to this group
		bool add_function(size_t fct);

		/// removes a subset from this group
		bool remove_function(size_t fct);

		/// clear all subsets
		void clear() {m_vFunction.clear();}

		/// number of functions in this group
		size_t num_fct() const {return m_vFunction.size();}

		/// get function id
		size_t operator[](size_t i) const {return fct_id(i);}

		/// get function id
		size_t fct_id(size_t i) const
		{
			UG_ASSERT(i < num_fct(), "Invalid function access");
			return m_vFunction[i];
		}

		/// returns true if subset is contained in this group
		bool containes_function(size_t fct) const;

	protected:
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
