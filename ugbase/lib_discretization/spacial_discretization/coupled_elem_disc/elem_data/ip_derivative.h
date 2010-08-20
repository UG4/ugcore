/*
 * ip_derivative.h
 *
 *  Created on: 13.08.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__COUPLED_ELEM_DISC__ELEM_DATA__IP_DERIVATIVE__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__COUPLED_ELEM_DISC__ELEM_DATA__IP_DERIVATIVE__

#include <vector>

namespace ug{

template <typename TEntry>
class IPDerivative
{
	public:
		// own type
		typedef IPDerivative<TEntry> this_type;

		// entry type used by algebra
		typedef TEntry entry_type;

	public:
		IPDerivative() {clear();}

		///////////////////////////
		// setup
		///////////////////////////

		/// clear all
		void clear()
		{
			m_vvEntries.clear();
		}

		/// set number of functions
		void set_num_fct(size_t numFct)
		{
			m_vvEntries.resize(numFct);
			for(size_t i=0; i < m_vvEntries.size(); ++i)
				m_vvEntries[i].clear();
		}

		/// set number of dofs for function
		void set_num_dofs(size_t fct, size_t numDoFs)
		{
			UG_ASSERT(fct < m_vvEntries.size(), "Index not valid.");
			m_vvEntries[fct].resize(numDoFs);
		}

		///////////////////////////
		// access
		///////////////////////////

		/// number of local functions
		size_t num_fct() const {return m_vvEntries.size();}

		/// number of dofs (sum)
		size_t num_dofs() const
		{
			size_t sum = 0;
			for(size_t fct = 0; fct < num_fct(); ++fct)
				sum += num_dofs(fct);
			return sum;
		}

		/// number of dofs per function
		size_t num_dofs(size_t fct) const
		{
			UG_ASSERT(fct < m_vvEntries.size(), "Index not valid.");
			return m_vvEntries[fct].size();
		}

		/// access to derivative of function fct
		entry_type& operator()(size_t fct, size_t dof)
		{
			UG_ASSERT(fct < m_vvEntries.size(), "Index not valid.");
			UG_ASSERT(dof < m_vvEntries[fct].size(), "Index not valid.");
			return m_vvEntries[fct][dof];
		}

		/// const access to derivative of function fct
		const entry_type& operator()(size_t fct, size_t dof) const
		{
			UG_ASSERT(fct < m_vvEntries.size(), "Index not valid.");
			UG_ASSERT(dof < m_vvEntries[fct].size(), "Index not valid.");
			return m_vvEntries[fct][dof];
		}

	protected:
		// entries
		std::vector<std::vector<entry_type> > m_vvEntries;
};


} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__COUPLED_ELEM_DISC__ELEM_DATA__IP_DERIVATIVE__ */
