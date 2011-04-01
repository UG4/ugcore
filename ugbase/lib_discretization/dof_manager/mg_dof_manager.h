/*
 * mg_dof_manager.h
 *
 *  Created on: 12.07.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__DOF_MANAGER__MG_DOF_MANAGER__
#define __H__LIB_DISCRETIZATION__DOF_MANAGER__MG_DOF_MANAGER__

#include <vector>

#include "lib_grid/lg_base.h"
#include "./function_pattern.h"
#include "lib_discretization/dof_manager/dof_distribution.h"

namespace ug{

/**
 * A MultiGridDoFManager handles the distribution of degrees of freedom on a
 * MultiGrid. It distributes the dof on each grid level and for the surface grid.
 * Thus, it creates num_level + 1 DoFDistributions
 *
 * \tparam 	TDoFDistribution	Type of DoF Distribution
 */
template <typename TDoFDistribution>
class MGDoFManager : public GridObserver
{
	public:
	///	DoF Distribution type
		typedef IDoFDistribution<TDoFDistribution> dof_distribution_type;

	public:
	///	Default Constructor
		MGDoFManager()
			: m_pMGSubsetHandler(NULL), m_pMultiGrid(NULL), m_pSurfaceView(NULL),
				m_pFuncPattern(NULL), m_pSurfDD(NULL)
		{
			m_vLevelDD.clear();
		};

	///	Constructor setting Function Pattern and Multi Grid Subset Handler
		MGDoFManager(MultiGridSubsetHandler& mgsh, FunctionPattern& dp)
			: m_pMGSubsetHandler(NULL), m_pMultiGrid(NULL), m_pSurfaceView(NULL),
				m_pFuncPattern(NULL), m_pSurfDD(NULL)
		{
			m_vLevelDD.clear();
			assign_multi_grid_subset_handler(mgsh);
			assign_function_pattern(dp);
		};

	/// set multi grid subset handler
		bool assign_multi_grid_subset_handler(MultiGridSubsetHandler& mgsh);

	/// set function pattern
		bool assign_function_pattern(FunctionPattern& dp);

	/// number of levels
		size_t num_levels() const
		{
		//	without SubsetHandler, we have no level information
			if(m_pMGSubsetHandler == NULL) return 0;

		//	forward request
			return m_pMGSubsetHandler->num_levels();
		}

	/// distribute dofs on all levels + surface level
		bool enable_dofs();

	/// distribute dofs on all levels
		bool enable_level_dofs();

	///	distribute dofs on surface grid
		bool enable_surface_dofs();

	///	returns if level dofs are enabled
		bool level_dofs_enabled() const {return m_vLevelDD.size() != 0;}

	///	returns if surface dofs are enabled
		bool surface_dofs_enabled() const {return m_pSurfDD != NULL;}

	///	returns Surface DoF Distribution
		dof_distribution_type* get_surface_dof_distribution()
		{
		// 	update surface distribution
			if(!surface_distribution_required())
			{
				UG_LOG("Cannot update surface distribution.\n");
				throw(UGFatalError("Surface DoF Distribution missing but requested."));
			}

			return m_pSurfDD;
		}

	///	returns Surface DoF Distribution
		const dof_distribution_type* get_surface_dof_distribution() const
		{
		// 	update surface distribution
			if(m_pSurfDD == NULL)
			{
				UG_LOG("surface distribution missing.\n");
				throw(UGFatalError("Surface DoF Distribution missing but requested."));
			}

			return m_pSurfDD;
		}

	///	returns Level DoF Distribution
		dof_distribution_type* get_level_dof_distribution(size_t level)
		{
			if(level < m_vLevelDD.size())
				return m_vLevelDD[level];
			else
				return NULL;
		}

	///	returns Level DoF Distribution
		const dof_distribution_type* get_level_dof_distribution(size_t level) const
		{
			if(level < m_vLevelDD.size())
				return m_vLevelDD[level];
			else
				return NULL;
		}

	///	returns the Level DoF Distributions in a vector
		std::vector<const dof_distribution_type*> get_level_dof_distributions() const
		{
			std::vector<const dof_distribution_type*> vLevelDD;
			for(size_t i = 0; i < m_vLevelDD.size(); ++i)
				vLevelDD.push_back(m_vLevelDD[i]);
			return vLevelDD;
		}

	///	returns the surface view
		const SurfaceView* get_surface_view() const {return m_pSurfaceView;}

	///	print a statistic on dof distribution
		void print_statistic() const;

	///	print a statistic on layout informations
		void print_layout_statistic() const;

	///	Destructor
		virtual ~MGDoFManager()
		{
			m_levelStorageManager.clear_subset_handler();
			m_surfaceStorageManager.clear_subset_handler();

			disable_level_dofs();
			disable_surface_dofs();
		}

	public:
	///////////////////////////////////////////////////
	//	GridObserver Callbacks
	///////////////////////////////////////////////////

	//	creation callbacks
	/**
	 *
	 */
	/// \{
		virtual void vertex_created(Grid* grid, VertexBase* vrt,
									GeometricObject* pParent = NULL,
									bool replacesParent = false);

		virtual void edge_created(Grid* grid, EdgeBase* e,
									GeometricObject* pParent = NULL,
									bool replacesParent = false);

		virtual void face_created(Grid* grid, Face* f,
									GeometricObject* pParent = NULL,
									bool replacesParent = false);

		virtual void volume_created(Grid* grid, Volume* vol,
									GeometricObject* pParent = NULL,
									bool replacesParent = false);
	///	\}

	//	erase callbacks
	/**
	 *
	 * \{ */
		virtual void vertex_to_be_erased(Grid* grid, VertexBase* vrt,
										 VertexBase* replacedBy = NULL);

		virtual void edge_to_be_erased(Grid* grid, EdgeBase* e,
										 EdgeBase* replacedBy = NULL);

		virtual void face_to_be_erased(Grid* grid, Face* f,
										 Face* replacedBy = NULL);

		virtual void volume_to_be_erased(Grid* grid, Volume* vol,
										 Volume* replacedBy = NULL);

	/**	\}	*/

		void defragment()
		{
		//	we add the shadows to the surface view, that may have been created
		//	due to grid adaption

		//	check grid options
			if(!m_pMultiGrid->option_is_enabled(VOLOPT_AUTOGENERATE_FACES)){
			  UG_LOG("WARNING: Auto-enabling grid option VOLOPT_AUTOGENERATE_FACES");
			  m_pMultiGrid->enable_options(VOLOPT_AUTOGENERATE_FACES);
			}
			if(!m_pMultiGrid->option_is_enabled(VOLOPT_AUTOGENERATE_FACES)){
			  UG_LOG("WARNING: Auto-enabling grid option FACEOPT_AUTOGENERATE_EDGES");
			  m_pMultiGrid->enable_options(FACEOPT_AUTOGENERATE_EDGES);
			}
			if(!m_pMultiGrid->option_is_enabled(VOLOPT_AUTOGENERATE_EDGES)){
			  UG_LOG("WARNING: Auto-enabling grid option VOLOPT_AUTOGENERATE_EDGES");
			  m_pMultiGrid->enable_options(FACEOPT_AUTOGENERATE_EDGES);
			}

		//	add missing shadows to surface view
			if(surface_dofs_enabled())
			{
				add_associated_sides_to_surface_view<Volume>();
				add_associated_sides_to_surface_view<Face>();
				add_associated_sides_to_surface_view<EdgeBase>();
			}

		//	defragment dof distributions
			if(surface_dofs_enabled())
				get_surface_dof_distribution()->defragment();
			if(level_dofs_enabled())
				for(size_t lev = 0; lev < num_levels(); ++lev)
					get_level_dof_distribution(lev)->defragment();

		}

		template <class TElem>
		void add_associated_sides_to_surface_view()
		{
			typedef typename geometry_traits<TElem>::const_iterator iterator;
			typedef typename TElem::lower_dim_base_object Side;
			std::vector<Side*> vSides;
			Grid& grid = *m_pSurfaceView->get_assigned_grid();

			for(size_t l  = 0; l < m_pSurfaceView->num_levels(); ++l){
				for(int si = 0; si < m_pSurfaceView->num_subsets(); ++si){
					for(iterator iter = m_pSurfaceView->begin<TElem>(si, l);
						iter != m_pSurfaceView->end<TElem>(si, l); ++iter)
					{
						TElem* e = *iter;
						CollectAssociated(vSides, grid, e);

						for(size_t i = 0; i < vSides.size(); ++i)
						{
							Side* s = vSides[i];

							if(m_pSurfaceView->get_subset_index(s) == -1)
							{
								add_to_surface_view(s);
								get_surface_dof_distribution()->grid_obj_added(s);
							}
						}
					}
				}
			}
		}

	protected:
	///	creates the surface view
		virtual bool surface_view_required();

	///	adds an element to the surface view
		void add_to_surface_view(GeometricObject* obj);
		void add_to_surface_view(VertexBase* vrt);
		void add_to_surface_view(EdgeBase* edge);
		void add_to_surface_view(Face* face);
		void add_to_surface_view(Volume* vol);

	///	removes an element to the surface view
		void remove_from_surface_view(GeometricObject* obj);
		void remove_from_surface_view(VertexBase* vrt);
		void remove_from_surface_view(EdgeBase* edge);
		void remove_from_surface_view(Face* face);
		void remove_from_surface_view(Volume* vol);

	///	creates level DoF Distributions iff needed
		bool level_distribution_required(size_t numLevel);

	///	deletes all level distributions
		void disable_level_dofs();

	///	adds an element to to level dof distribution
		void add_to_level_dof_distribution(GeometricObject* vrt);
		void add_to_level_dof_distribution(VertexBase* vrt);
		void add_to_level_dof_distribution(EdgeBase* edge);
		void add_to_level_dof_distribution(Face* face);
		void add_to_level_dof_distribution(Volume* vol);

	///	removes an element from level dof distribution
		void remove_from_level_dof_distribution(GeometricObject* vrt);
		void remove_from_level_dof_distribution(VertexBase* vrt);
		void remove_from_level_dof_distribution(EdgeBase* edge);
		void remove_from_level_dof_distribution(Face* face);
		void remove_from_level_dof_distribution(Volume* vol);

	///	creates the surface distribution
		bool surface_distribution_required();

	///	deletes the surface distributions
		void disable_surface_dofs();

	/// print statistic for a DoFDistribution
		void print_statistic(const dof_distribution_type& dd) const;

	/// print statistic on layouts for a DoFDistribution
		void print_layout_statistic(const dof_distribution_type& dd) const;

	protected:
	// 	MultiGridSubsetHandler this DofManager works on
		MultiGridSubsetHandler* m_pMGSubsetHandler;

	// 	MultiGrid associated to the SubsetHandler
		MultiGrid* m_pMultiGrid;

	// 	Surface View
		SurfaceView* m_pSurfaceView;

	// 	DoF Pattern
		FunctionPattern* m_pFuncPattern;

	// 	Level DoF Distributors
		std::vector<dof_distribution_type*> m_vLevelDD;

	// 	Surface Grid DoF Distributor
		dof_distribution_type* m_pSurfDD;

	// 	Storage manager
		typename TDoFDistribution::storage_manager_type	m_levelStorageManager;
		typename TDoFDistribution::storage_manager_type	m_surfaceStorageManager;
};

} // end namespace ug

// include implementation
#include "mg_dof_manager_impl.h"

#endif /* __H__LIB_DISCRETIZATION__DOF_MANAGER__MG_DOF_MANAGER__ */
