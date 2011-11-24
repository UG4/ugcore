/*
 * mg_dof_manager.h
 *
 *  Created on: 12.07.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__DOF_MANAGER__MG_DOF_MANAGER__
#define __H__UG__LIB_DISC__DOF_MANAGER__MG_DOF_MANAGER__

#include <vector>

#include "lib_grid/lg_base.h"
#include "./function_pattern.h"
#include "lib_disc/dof_manager/dof_distribution.h"

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
			: m_bGrouped(false), m_pMGSubsetHandler(NULL), m_pMultiGrid(NULL),
			  m_pSurfaceView(NULL), m_pFuncPattern(NULL), m_pSurfDD(NULL)
		{
			m_vLevelDD.clear();
		};

	///	Constructor setting Function Pattern and Multi Grid Subset Handler
		MGDoFManager(MultiGridSubsetHandler& mgsh, FunctionPattern& dp)
			: m_bGrouped(false), m_pMGSubsetHandler(NULL), m_pMultiGrid(NULL),
			  m_pSurfaceView(NULL), m_pFuncPattern(NULL), m_pSurfDD(NULL)
		{
			m_vLevelDD.clear();
			assign_multi_grid_subset_handler(mgsh);
			assign_function_pattern(dp);
		};

	///	set grouped
		void set_grouping(bool bGrouped) {m_bGrouped = bGrouped;}

	/// set multi grid subset handler
		bool assign_multi_grid_subset_handler(MultiGridSubsetHandler& mgsh);

	/// set function pattern
		bool assign_function_pattern(FunctionPattern& dp);

	/// number of levels
		virtual size_t num_levels() const
		{
		//	without SubsetHandler, we have no level information
			if(m_pMGSubsetHandler == NULL) return 0;

		//	forward request
			return m_pMGSubsetHandler->num_levels();
		}

	/// distribute dofs on all levels + surface level
		bool enable_indices();

	/// distribute dofs on all levels
		bool enable_level_indices();

	///	distribute dofs on surface grid
		bool enable_surface_indices();

	///	returns if level dofs are enabled
		bool level_indices_enabled() const {return m_vLevelDD.size() != 0;}

	///	returns if surface dofs are enabled
		bool surface_indices_enabled() const {return m_pSurfDD != NULL;}

	///	returns Surface DoF Distribution
		dof_distribution_type* surface_dof_distribution()
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
		const dof_distribution_type* surface_dof_distribution() const
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
		dof_distribution_type* level_dof_distribution(size_t level)
		{
			if(!level_distribution_required(level+1))
			{
				UG_LOG("Cannot update level distribution.\n");
				throw(UGFatalError("Cannot update level distribution.\n"));
			}

			return m_vLevelDD[level];
		}

	///	returns Level DoF Distribution
		const dof_distribution_type* level_dof_distribution(size_t level) const
		{
			if(!(level < m_vLevelDD.size()))
				throw(UGFatalError("Level DoF Distribution missing"));

			return m_vLevelDD[level];
		}

	///	returns the Level DoF Distributions in a vector
		std::vector<const dof_distribution_type*> level_dof_distributions() const
		{
			std::vector<const dof_distribution_type*> vLevelDD;
			for(size_t i = 0; i < m_vLevelDD.size(); ++i)
				vLevelDD.push_back(m_vLevelDD[i]);
			return vLevelDD;
		}

	///	returns the surface view
		const SurfaceView* surface_view() const {return m_pSurfaceView;}

	///	print a statistic on dof distribution
		virtual void print_statistic(int verboseLev = 1) const;

	///	print a statistic on layout informations
		virtual void print_layout_statistic(int verboseLev = 1) const;

	///	print a statistic on local dofs
		virtual void print_local_dof_statistic(int verboseLev = 1) const;

	///	Destructor
		virtual ~MGDoFManager()
		{
			m_levelStorageManager.clear();
			m_surfaceStorageManager.clear();

			disable_level_indices();
			disable_surface_indices();
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

	///	defragments the index set
		void defragment();

	protected:

	///	adds lower-dim elements ot surface view
		template <class TElem>
		void add_associated_sides_to_surface_view();

	///	creates the surface view
		virtual bool surface_view_required();

	///	returns if an element belongs to the surface view
		bool is_in_surface_view(GeometricObject* obj);
		bool is_in_surface_view(VertexBase* vrt);
		bool is_in_surface_view(EdgeBase* edge);
		bool is_in_surface_view(Face* face);
		bool is_in_surface_view(Volume* vol);

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
		void disable_level_indices();

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
		void disable_surface_indices();

	/// print statistic for a DoFDistribution
		void print_statistic(const dof_distribution_type& dd, int verboseLev = 1) const;

	/// print statistic on local dof distribution
		void print_local_dof_statistic(const dof_distribution_type& dd, int verboseLev = 1) const;

	protected:
	//	group flag
		bool m_bGrouped;

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

#endif /* __H__UG__LIB_DISC__DOF_MANAGER__MG_DOF_MANAGER__ */
