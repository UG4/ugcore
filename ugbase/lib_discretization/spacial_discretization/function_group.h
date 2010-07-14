/*
 * function_group.h
 *
 *  Created on: 13.07.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__FUNCTION_GROUP__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__FUNCTION_GROUP__

#include <vector>

namespace ug{

// FunctionGroup is just a group of fctze_t, representing some functions
class FunctionGroup
{
	public:
		FunctionGroup() {m_vFunction.clear();}

		/// adds a subset to this group
		bool add_function(size_t fct);

		/// removes a subset from this group
		bool remove_function(size_t fct);

		/// clear all subsets
		void clear() {m_vFunction.clear();}

		/// number of functions in this group
		size_t num_functions() const {return m_vFunction.size();}

		/// get function id
		size_t operator[](size_t i) const {return m_vFunction[i];}

		/// returns true if subset is contained in this group
		bool containes_function(size_t fct) const;

	protected:
		std::vector<size_t> m_vFunction;
};

} // end namespace ug

#endif /*__H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__FUNCTION_GROUP__ */
