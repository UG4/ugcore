/*
 * parallel_dof_manager.h
 *
 *  Created on: 21.5.2010
 *      Author: A. Vogel
 */

#ifndef __H__LIB_DISCRETIZATION__PARALLELIZATION__PARALLEL_DOF_MANAGER__
#define __H__LIB_DISCRETIZATION__PARALLELIZATION__PARALLEL_DOF_MANAGER__

#include "parallelization.h"

namespace ug
{

template <typename TDoFManager>
class ParallelDoFManager : public TDoFManager
{
	public:
		typedef typename TDoFManager::geom_obj_container_type geom_obj_container_type;
	public:
		ParallelDoFManager(geom_obj_container_type& geomObjCont, GridLayoutMap& layoutMap)
		: TDoFManager(geomObjCont), m_layoutMap(layoutMap)
		{}

		bool create_index_layouts()
		{
			bool bRetVal = true;

			uint num_levels = this->num_levels();
			m_slaveLayouts.clear();
			m_slaveLayouts.resize(num_levels);
			for(uint l = 0; l < num_levels; ++l)
			{
				bRetVal &= CreateIndexLayout(m_masterLayouts[l], *this, m_layoutMap, INT_MASTER,l);
				bRetVal &= CreateIndexLayout(m_slaveLayouts[l], *this, m_layoutMap, INT_SLAVE,l);
			}

			return bRetVal;
		}

	private:
		// Layout map of grid
		GridLayoutMap& m_layoutMap;

		// index layout for each grid level
		std::vector<IndexLayout> m_slaveLayouts;

		// index layout for each grid level
		std::vector<IndexLayout> m_masterLayouts;

};


} // end namespace ug

#endif
