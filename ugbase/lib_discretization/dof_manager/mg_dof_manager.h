/*
 * mg_dof_manager.h
 *
 *  Created on: 12.07.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__DOF_MANAGER__MG_DOF_MANAGER__
#define __H__LIB_DISCRETIZATION__DOF_MANAGER__MG_DOF_MANAGER__

#include <vector>

#include "lib_grid/lib_grid.h"
#include "./function_pattern.h"
#include "lib_discretization/geom_object_container/surface_view.h"
#include "lib_grid/algorithms/subset_util.h"

namespace ug{

template <typename TDoFDistribution>
class MGDoFManager
{
	public:
		typedef TDoFDistribution dof_distribution_type;

	public:
		MGDoFManager()
			: m_pMGSubsetHandler(NULL), m_pMultiGrid(NULL), m_pSurfaceView(NULL),
				m_pFunctionPattern(NULL), m_pSurfaceDoFDistribution(NULL)
		{
			m_vLevelDoFDistribution.clear();
		};

		MGDoFManager(MultiGridSubsetHandler& mgsh, FunctionPattern& dp)
			: m_pMGSubsetHandler(NULL), m_pMultiGrid(NULL), m_pSurfaceView(NULL),
				m_pFunctionPattern(NULL), m_pSurfaceDoFDistribution(NULL)
		{
			m_vLevelDoFDistribution.clear();
			assign_multi_grid_subset_handler(mgsh);
			assign_function_pattern(dp);
		};

		// set multi grid subset handler
		bool assign_multi_grid_subset_handler(MultiGridSubsetHandler& mgsh)
		{
			// Remember Subsethandler and MultiGrid
			m_pMGSubsetHandler = &mgsh;
			m_pMultiGrid = m_pMGSubsetHandler->get_assigned_multi_grid();

			// Get StorageManager for levels
			m_levelStorageManager.set_subset_handler(mgsh);

			// set Function pattern if already assigned
			if(m_pFunctionPattern != NULL)
				if(!assign_function_pattern(*m_pFunctionPattern))
					return false;

			return true;
		}

		// set multi grid subset handler
		bool assign_function_pattern(FunctionPattern& dp)
		{
			m_pFunctionPattern = &dp;

			// if already subsethandler set
			if(m_pMGSubsetHandler != NULL)
			{
				// remember current levels
				size_t num_level = m_vLevelDoFDistribution.size();

				// free memory
				delete_distributions();

				// reallocate for new pattern
				if(num_level > 0)
					if(!level_distribution_required(num_level-1))
						return false;
			}
			return true;
		}

		// number of levels
		inline size_t num_levels() const {return m_pMGSubsetHandler->num_levels();}

		// distribute dofs on all levels + surface level
		bool distribute_dofs()
		{
			// require distributions on all levels
			if(!level_distribution_required(num_levels()-1))
				{UG_LOG("Cannot access distribution of level.\n"); return false;}

			// distribute on level grids
			for(size_t l = 0; l < num_levels(); ++l)
			{
				if(!m_vLevelDoFDistribution[l]->distribute_dofs())
					{UG_LOG("Cannot distribute dofs on level "<<l<<".\n"); return false;}
			}

			// distribute surface dofs
			return distribute_surface_dofs();
		}

		bool distribute_surface_dofs()
		{
			return true;

			// update surface distribution
			if(!update_surface_distribution())
				{UG_LOG("Cannot update surface distribution.\n"); return false;}

			// distribute on surface grid
			if(!m_pSurfaceDoFDistribution->distribute_dofs())
				{UG_LOG("Cannot distribute dofs on surface.\n"); return false;}

			return true;
		}

		TDoFDistribution* get_surface_dof_distribution()
		{
			// TODO: SHOULD BE:
			//return m_pSurfaceDoFDistribution;

			if(!level_distribution_required(num_levels()-1)) return NULL;
			return m_vLevelDoFDistribution[num_levels()-1];
		}


		TDoFDistribution* get_level_dof_distribution(size_t level)
		{
			if(!level_distribution_required(level)) return NULL;
			return m_vLevelDoFDistribution[level];
		}

		virtual ~MGDoFManager()
		{
			delete_distributions();
		}

	protected:
		virtual SurfaceView* create_surface_view()
		{
			// serial version
			SurfaceView* pSurfaceView = new SurfaceView(*m_pMultiGrid);
			if(pSurfaceView == NULL)
				{UG_LOG("Allocation of Surface View failed.\n"); return false;}

			// Create surface view for all elements
			CreateTopView(*m_pMultiGrid, *pSurfaceView, m_pMultiGrid->vertices_begin(), m_pMultiGrid->vertices_end());
			CreateTopView(*m_pMultiGrid, *pSurfaceView, m_pMultiGrid->edges_begin(), m_pMultiGrid->edges_end());
			CreateTopView(*m_pMultiGrid, *pSurfaceView, m_pMultiGrid->faces_begin(), m_pMultiGrid->faces_end());
			CreateTopView(*m_pMultiGrid, *pSurfaceView, m_pMultiGrid->volumes_begin(), m_pMultiGrid->volumes_end());

			// set storage manager
			m_surfaceStorageManager.set_subset_handler(*pSurfaceView);

			return pSurfaceView;
		}

		inline bool update_surface_distribution()
		{
			// Create surface view, if it does not exist
			if(m_pSurfaceView == NULL)
				m_pSurfaceView = create_surface_view();
			if(m_pSurfaceView == NULL)
				{UG_LOG("Cannot create surface view."); return false;}

			// Create surface dof distributions
			if(m_pSurfaceDoFDistribution == NULL)
				m_pSurfaceDoFDistribution = new TDoFDistribution(m_pSurfaceView->get_geometric_object_collection(), *m_pSurfaceView,
																m_surfaceStorageManager, *m_pFunctionPattern);
			if(m_pSurfaceDoFDistribution == NULL)
				{UG_LOG("Cannot allocate Surface DoF Distribution.\n"); return false;}
			return true;
		}

		inline bool level_distribution_required(size_t level)
		{
			if(level >= m_pMGSubsetHandler->num_levels())
				{UG_LOG("Level DoF Distribution required on Level " << level << ", but MGSubsetHandler has only " << m_pMGSubsetHandler->num_levels()<< ".\n"); return false;}

			// Create level dof distributions
			for(size_t l = m_vLevelDoFDistribution.size(); l <= level; ++l)
			{
				m_vLevelDoFDistribution.push_back(new TDoFDistribution(m_pMGSubsetHandler->get_goc_by_level(l), *m_pMGSubsetHandler,
																		m_levelStorageManager, *m_pFunctionPattern));
				if(m_vLevelDoFDistribution[l] == NULL)
					{UG_LOG("Cannot allocate Level DoF Distribution on Level " << l << ".\n"); return false;}
			}
			return true;
		}

		inline void delete_distributions()
		{
			// Delete surface dof distributions
			if(m_pSurfaceDoFDistribution != NULL)
				delete m_pSurfaceDoFDistribution;
			m_pSurfaceDoFDistribution = NULL;

			// delete surface view
			if(m_pSurfaceView != NULL)
				delete m_pSurfaceView;
			m_pSurfaceView = NULL;

			// Delete level dof distributions
			for(size_t l = 0; l < m_vLevelDoFDistribution.size(); ++l)
			{
				delete m_vLevelDoFDistribution[l];
				m_vLevelDoFDistribution[l] = NULL;
			}
			m_vLevelDoFDistribution.clear();
		}

	protected:
		// MultiGridSubsetHandler this DofManager works on
		MultiGridSubsetHandler* m_pMGSubsetHandler;

		// MultiGrid associated to the SubsetHandler
		MultiGrid* m_pMultiGrid;

		// Surface View
		SurfaceView* m_pSurfaceView;

		// DoF Pattern
		FunctionPattern* m_pFunctionPattern;

		// Level DoF Distributors
		std::vector<TDoFDistribution*> m_vLevelDoFDistribution;

		// Surface Grid DoF Distributor
		TDoFDistribution* m_pSurfaceDoFDistribution;

		// Storage manager
		typename TDoFDistribution::StorageManager	m_levelStorageManager;
		typename TDoFDistribution::StorageManager	m_surfaceStorageManager;
};









} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__DOF_MANAGER__MG_DOF_MANAGER__ */
