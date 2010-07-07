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

			size_t num_levels = TDoFManager::num_levels();
			m_slaveLayouts.clear();
			m_slaveLayouts.resize(num_levels);
			m_masterLayouts.clear();
			m_masterLayouts.resize(num_levels);
			m_verticalSlaveLayouts.clear();
			m_verticalSlaveLayouts.resize(num_levels);
			m_verticalMasterLayouts.clear();
			m_verticalMasterLayouts.resize(num_levels);
			m_processCommunicators.resize(num_levels);

		//TODO:	this communicator should be specified from the application
			pcl::ProcessCommunicator commWorld;

			for(size_t l = 0; l < num_levels; ++l)
			{
				bRetVal &= CreateIndexLayout(m_masterLayouts[l],
							*this, *m_pLayoutMap, INT_MASTER,l);
				bRetVal &= CreateIndexLayout(m_slaveLayouts[l],
							*this, *m_pLayoutMap, INT_SLAVE,l);
				bRetVal &= CreateIndexLayout(m_verticalMasterLayouts[l],
							*this, *m_pLayoutMap, INT_VERTICAL_MASTER,l);
				bRetVal &= CreateIndexLayout(m_verticalSlaveLayouts[l],
							*this, *m_pLayoutMap, INT_VERTICAL_SLAVE,l);

			//	create local process communicator
			//	if a process has only vertical slaves, it is not involved in process communication.
			//TODO: perform a more precise check
				bool participate = m_verticalSlaveLayouts[l].empty()
								   && !commWorld.empty()
								   && (TDoFManager::num_dofs(l) > 0);
				m_processCommunicators[l] = commWorld.create_sub_communicator(participate);
			}

			return bRetVal;
		}

		inline IndexLayout& get_slave_layout(size_t level)	{return m_slaveLayouts[level];}
		inline IndexLayout& get_master_layout(size_t level)	{return m_masterLayouts[level];}
		inline IndexLayout& get_vertical_slave_layout(size_t level)		{return m_verticalSlaveLayouts[level];}
		inline IndexLayout& get_vertical_master_layout(size_t level)	{return m_verticalMasterLayouts[level];}

		inline pcl::ParallelCommunicator<IndexLayout>& get_communicator()	{return m_communicator;}

		inline pcl::ProcessCommunicator get_process_communicator(size_t level)	{return m_processCommunicators[level];}

	private:
		// Layout map of grid
		GridLayoutMap* m_pLayoutMap;

		// index layout for each grid level
		std::vector<IndexLayout> m_slaveLayouts;

		// index layout for each grid level
		std::vector<IndexLayout> m_masterLayouts;

		// index layout for each grid level
		std::vector<IndexLayout> m_verticalMasterLayouts;

		// index layout for each grid level
		std::vector<IndexLayout> m_verticalSlaveLayouts;

		// communicator
		pcl::ParallelCommunicator<IndexLayout> m_communicator;

		// process communicator
		std::vector<pcl::ProcessCommunicator> m_processCommunicators;

};

} // end namespace ug

#endif
