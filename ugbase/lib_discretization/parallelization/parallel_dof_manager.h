/*
 * parallel_dof_manager.h
 *
 *  Created on: 21.5.2010
 *      Author: A. Vogel
 */

#ifndef __H__LIB_DISCRETIZATION__PARALLELIZATION__PARALLEL_DOF_MANAGER__
#define __H__LIB_DISCRETIZATION__PARALLELIZATION__PARALLEL_DOF_MANAGER__

#include "pcl/pcl.h"
#include "parallelization_util.h"

namespace ug
{

template <typename TDoFManager>
class ParallelDoFManager : public TDoFManager
{
	public:
		typedef typename TDoFManager::geom_obj_container_type geom_obj_container_type;

	public:
		ParallelDoFManager(geom_obj_container_type& geomObjCont)
		: TDoFManager(geomObjCont), m_pLayoutMap(NULL)
		{}

		ParallelDoFManager(geom_obj_container_type& geomObjCont, GridLayoutMap& layoutMap)
		: TDoFManager(geomObjCont), m_pLayoutMap(&layoutMap)
		{}

		void set_grid_layout_map(GridLayoutMap& layoutMap)
		{
			m_pLayoutMap = &layoutMap;
		}

		//bool create_index_layouts()
		virtual bool finalize()
		{
			if(!m_pLayoutMap){
				UG_LOG("  no layout map specified. aborting.\n");
				return false;
			}

			TDoFManager::finalize();

			bool bRetVal = true;

			uint num_levels = this->num_levels();
			m_slaveLayouts.clear();
			m_slaveLayouts.resize(num_levels);
			m_masterLayouts.clear();
			m_masterLayouts.resize(num_levels);

			for(uint l = 0; l < num_levels; ++l)
			{
/*
				bRetVal &= CreateIndexLayout(m_masterLayouts[l],
							*static_cast<TDoFManager*>(this), *m_pLayoutMap, INT_MASTER,l);
				bRetVal &= CreateIndexLayout(m_slaveLayouts[l],
							*static_cast<TDoFManager*>(this), *m_pLayoutMap, INT_SLAVE,l);
*/
				bRetVal &= CreateIndexLayout(m_masterLayouts[l],
							*this, *m_pLayoutMap, INT_MASTER,l);
				bRetVal &= CreateIndexLayout(m_slaveLayouts[l],
							*this, *m_pLayoutMap, INT_SLAVE,l);

			}

			return bRetVal;
		}

		inline IndexLayout& get_slave_layout(size_t level)	{return m_slaveLayouts[level];}
		inline IndexLayout& get_master_layout(size_t level)	{return m_masterLayouts[level];}
		inline pcl::ParallelCommunicator<IndexLayout>& get_communicator()	{return m_Communicator;}

	private:
		// Layout map of grid
		GridLayoutMap* m_pLayoutMap;

		// index layout for each grid level
		std::vector<IndexLayout> m_slaveLayouts;

		// index layout for each grid level
		std::vector<IndexLayout> m_masterLayouts;

		// communicator
		pcl::ParallelCommunicator<IndexLayout> m_Communicator;

};

} // end namespace ug

#endif
