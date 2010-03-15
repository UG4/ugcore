/*
 * dofpattern.h
 *
 *  Created on: 05.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__DOF_MANAGER__DOF_MANAGER__
#define __H__LIBDISCRETIZATION__DOF_MANAGER__DOF_MANAGER__

#include "common/common.h"
#include "lib_grid/lib_grid.h"
#include "lib_algebra/lib_algebra.h"

// trial spaces
#include "lib_discretization/local_shape_function_set/local_dof_pattern.h"
#include "lib_discretization/local_shape_function_set/local_shape_function_set_factory.h"

// multi indices
#include "lib_algebra/multi_index/multi_indices.h"

// dof manager strategies
#include "lib_discretization/dof_manager/dof_manager.h"

namespace ug{

////////////////////////////////////////////////
// P1ConformDoFManager
/// manages the distribution of degrees of freedom for lagrange p1 functions
/**
 * Manages the distribution of DoFs on a Domain, that has been separated into subsets.
 * Given a SubsetHandler, only Lagrange P1 solutions can be added to selected subsets. The
 * LocalShapeFunctionSetID must be LSFS_LAGRANGEP1.
 * All solutions are grouped automatically.
 * Invoking the finalize command let the P1ConformDoFManager distribute the DoFs.
 */
class P1ConformDoFManager : public GridObserver{
	public:
		// type of multiindex used
		typedef MultiIndex<1> multi_index_type;

		// type of single index (e.g. uint or int)
		typedef multi_index_type::single_index_type single_index_type;

		// value container for element local indices
		typedef std::vector<multi_index_type> local_index_type;

		// type of grid
		typedef MultiGrid grid_type;

		// element container type (i.e. grid or subset_handler)
		typedef grid_type element_container_type;

		// type of subset handler
		typedef MultiGridSubsetHandler subset_handler_type;

	public:
		/// constructor
		/**
		 * \param[in] name	Name of this P1ConformDoFManager
		 * \param[in] sh	SubsetHandler
		 */
		P1ConformDoFManager(grid_type& grid);

		/// add a single solution of LocalShapeFunctionSetID to the entire domain (i.e. all elements of the (Multi-)grid)
		/**
		 * \param[in] name			Name of this Single Solution
		 * \param[in] TrialSpace	Trial Space for this function
		 */
		bool add_discrete_function(std::string name, LocalShapeFunctionSetID id, int dim);

		/// add a single solution of LocalShapeFunctionSetID to selected subsets
		/**
		 * \param[in] name			Name of this Single Solution
		 * \param[in] TrialSpace	Trial Space for this function
		 * \param[in] SubsetIndecex	Std::Vector of subset indeces, where this solution lives
		 */
		bool add_discrete_function(std::string name, LocalShapeFunctionSetID id, std::vector<int>& SubsetIndices, int dim);

		/// group single solutions
		/**
		 * By this function a user can group single solutions to a new one. The single solutions will be
		 * removed from the pattern and a new group solution containing those will be added. It is
		 * also possible to group 'Group Solutions'.
		 * The Grouping has an effect on the dof numbering. While single solutions get an own index with only one component for
		 * each dof, a Grouped Solution has on index with several components.
		 *
		 */
		bool group_discrete_functions(std::string name, std::vector<uint>& GroupSolutions);

		/// groups selected dof-groups associated with one Geom Object
		template <typename TGeomObj>
		bool group_dof_groups(std::vector<uint>& dofgroups);

		/// groups all dofgroups on a given geometric object type, where a component of selected functions appear
		template <typename TGeomObj>
		bool group_discrete_functions(std::vector<uint>& selected_functions);

		/// groups all dof on every geometric object for selected functions
		bool group_discrete_functions(std::vector<uint>& selected_functions);

		/// resets grouping for a geometric object
		template <typename TGeomObj>
		bool reset_grouping();

		/// clears grouping for all geometric objects
		bool reset_grouping();

		/// gives informations about the current status
		bool print_info();

		/// performs a finalizing step. The Pattern can not be altered after finishing
		bool finalize();

		/// returns the number of dofs on level 'level'
		inline uint num_dofs(uint level = 0);

		/// returns the number of levels
		inline uint num_levels();

		/// returns the trial space of the discrete function nr_fct
		LocalShapeFunctionSetID get_local_shape_function_set_id(uint nr_fct);

		// returns the number of multi_indices on the Element for the discrete function 'nr_fct'
		template<typename TElem>
		std::size_t num_multi_indices(TElem* elem, uint nr_fct);

		/// returns the indices of the dofs on the Element elem for the discrete function 'nr_fct' and returns num_indices
		template<typename TElem>
		std::size_t get_multi_indices(TElem* elem, uint nr_fct, local_index_type& ind);

		/// returns the index of the dofs on the Geom Obj for the discrete function 'nr_fct'
		template<typename TElem>
		std::size_t get_multi_indices_of_geom_obj(TElem* vrt, uint nr_fct, local_index_type& ind);

		/// returns the name of the discrete function nr_fct
		std::string get_name(uint nr_fct);

		/// returns the number of discrete functions in this dof pattern
		inline uint num_fct();

		/// returns true if the discrete function nr_fct is defined on subset s
		bool fct_def_in_subset(uint nr_fct, int s);

		/// returns the assigned SubsetHandler
		element_container_type& get_assigned_element_container();

		const IndexInfo& get_index_info(uint level)
		{
			return m_IndexInfo[level];
		}

		/// destructor
		~P1ConformDoFManager();

	protected:
		// if locked, pattern can not be changed anymore
		bool m_bLocked;

		// names of functions
		std::vector<std::string> m_SingleDiscreteFunctionNames;
		single_index_type m_num_single_discrete_functions;

		// attachment type
		typedef ug::Attachment<single_index_type> AIndex;
		AIndex m_aIndex;

		// Accessor type
		Grid::AttachmentAccessor<VertexBase, AIndex> m_aaIndex;
		std::vector<single_index_type> m_num_dof_index;

		// number of Levels
		uint m_num_levels;
		std::vector<IndexInfo> m_IndexInfo;

		// associated Grid
		grid_type& m_grid;
};




}

#include "p1conform_dof_manager_impl.h"


#endif /* __H__LIBDISCRETIZATION__DOF_MANAGER__DOF_MANAGER__ */
