/*
 * dof_manager.h
 *
 *  Created on: 05.03.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__DOF_MANAGER__DOF_MANAGER__
#define __H__LIB_DISCRETIZATION__DOF_MANAGER__DOF_MANAGER__

namespace ug{

enum DoFManagerGroupStrategy {
	DMGS_INVALID = -1,
	DMGS_NO_GROUPING = 0,
	DMGS_GROUP_ALL,
	DMGS_NUM_STRATEGIES
};

}

template <typename TMultiIndex>
class DoFValueContainer{
	public:
		typedef TMultiIndex index_type;
		typedef number value_type;
	public:
		DoFValueContainer(){};

		void resize(std::size_t n)
		{
			m_values.resize(n);
			m_indices.resize(n);
		}

		std::size_t size()
		{
			return m_values.size();
		}

		const TMultiIndex& index(std::size_t i)
		{
			assert(i < m_indices.size());
			return m_indices[i];
		}

		value_type value(std::size_t i)
		{
			assert(i < m_values.size());
			return m_values[i];
		}

	protected:
		std::vector<value_type> m_values;
		std::vector<TMultiIndex> m_indices;
};

#endif /* __H__LIB_DISCRETIZATION__DOF_MANAGER__DOF_MANAGER__ */
