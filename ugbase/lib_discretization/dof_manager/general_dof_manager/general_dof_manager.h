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
// GeneralDoFManager
/// manages the distribution of degrees of freedom
/**
 * Manages the distribution of DoFs on a Domain, that has been separated into subsets.
 * Given a SubsetHandler, single solutions can be added to selected subsets. The chosen
 * LocalShapeFunctionSetID determines the number of DoFs on the Geometric Objects in the Subset.
 * The user can group several solutions to a multivalued solution.
 * Invoking the finalize command let the GeneralDoFManager distribute the DoFs.
 */
class GeneralDoFManager : public GridObserver{
	public:
		typedef MultiIndex<4> multi_index_type;

		typedef multi_index_type::single_index_type single_index_type;

	protected:
		typedef multi_index_type::size_type size_type;


	protected:
		typedef MultiGrid grid_type;
		typedef MultiGridSubsetHandler subset_handler_type;

	public:
		/// constructor
		/**
		 * \param[in] name	Name of this GeneralDoFManager
		 * \param[in] sh	SubsetHandler
		 */
		GeneralDoFManager(subset_handler_type& sh);

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
		size_type num_multi_indices(TElem* elem, uint nr_fct);

		/// returns the indices of the dofs on the Element elem for the discrete function 'nr_fct' and returns num_indices
		template<typename TElem>
		size_type get_multi_indices(TElem* elem, uint nr_fct, multi_index_type ind[]);

		/// returns the index of the dofs on the Geom Obj for the discrete function 'nr_fct'
		template<typename TElem>
		size_type get_multi_indices_of_geom_obj(TElem* vrt, uint nr_fct, multi_index_type ind[]);

		/// returns the name of the discrete function nr_fct
		std::string get_name(uint nr_fct);

		/// returns the number of discrete functions in this dof pattern
		inline uint num_fct();

		/// returns true if the discrete function nr_fct is defined on subset s
		bool fct_def_in_subset(uint nr_fct, int s);

		/// returns the assigned SubsetHandler
		subset_handler_type& get_assigned_subset();

		const IndexInfo& get_index_info(uint level)
		{
			return m_IndexInfo[level];
		}

		/// destructor
		~GeneralDoFManager();

	private:
		template <typename TGeomObj>
		bool add_dof_group(uint nr_fct, int subsetIndex, uint num_dofs);

		template <typename TGeomObj>
		bool finalize();

		template <typename TGeomObj>
		uint set_dof_group_ids(int subsetIndex, uint start_id);

		bool get_dof_group_index_info(int subsetIndex, uint group_id, uint level, uint& num_indices, uint& num_comp);

		template <typename TGeomObj>
		bool get_dof_group_index_info(int subsetIndex, uint group_id, uint level, uint& num_indices, uint& num_comp);

	protected:
		class DoFGroupMap{
		public:
			DoFGroupMap();
			DoFGroupMap(const DoFGroupMap& map);
			DoFGroupMap& operator= (const DoFGroupMap& map);

			bool set_num_dofs(uint num_dofs);
			bool set_nr_dofgroup(uint dof, uint group);
			bool set_comp_in_dofgroup(uint dof, uint comp);

			~DoFGroupMap();

			inline uint num_dofs() const;
			inline uint nr_dofgroup(uint i) const;
			inline uint comp_in_dofgroup(uint i) const;

		private:
			uint m_num_dofs;
			uint* m_nr_dofgroup;
			uint* m_comp_in_dofgroup;
		};

		/**
		 * A DoFGroup is a grouping of degrees of freedom associated with a TElem,
		 * that are numbered with the same 'i'-index, but different \alpha index.
		 * Thus all indices have the form (i,0), (i,1), (i,2), ... , (i, num_dof)
		 *
		 * Currently, the num_dofs must be the same for all DoFs
		 *
		 * A DOFGROUP lives on ONE SUBSET
		 *
		 */
		template <typename TElem>
		class DoFGroup
		{
			typedef ug::Attachment<single_index_type> AIndex;

			public:
				DoFGroup();
				DoFGroup(subset_handler_type& sh, int s) ;
				// default copy constructor and default operator= are sufficent

				~DoFGroup();
				bool activate(uint level);

				bool deactivate();

				single_index_type add_dof(uint nr_fct, uint fct_comp);

				inline single_index_type index(TElem* elem) const;
				inline single_index_type num_dofs() const;
				inline single_index_type num_indices(uint level) const;

				uint id() const {return m_dof_group_id;};
				void set_id(uint id) {m_dof_group_id = id;};

				subset_handler_type* get_subset_handler() const;
				int get_subset_index() const;

				uint get_nr_fct(uint comp) const;
				uint get_comp_fct(uint comp) const;

			private:
				AIndex m_aIndex;
				typename ISubsetHandler::AttachmentAccessor<TElem, AIndex> m_aaIndex;

			private:
				subset_handler_type* m_sh;		// subset handler
				int m_subsetIndex;				// subset index, where dof group lives

			private:
				uint m_dof_group_id;

			private:
				bool m_pActivated;
				std::vector<single_index_type> m_num_dof_index;

			private:
				single_index_type m_num_dof;		// number of dof in this dof group
				std::vector<uint> m_nr_fct;		// discrete function associated with each comp (DO WE NEED THIS ???)
				std::vector<uint> m_nr_comp;    // componenent of discrete fct associated with comp
		};

		template <typename TGeomObj>
		inline std::vector<DoFGroup<TGeomObj> >& get_dof_group_container();

		template <typename TGeomObj> class dof_group_traits{};

		/// Single discrete function
		/**
		 * This class represents a single discrete function
		 */
		class SingleDiscreteFunction {
		public:
			void set_name(std::string name);
			std::string get_name() const;

			void set_subset_handler(subset_handler_type& sh);
			subset_handler_type* get_subset_handler() const;

			void set_continuous_dof_pattern(ContinuousDoFPattern& p);
			const ContinuousDoFPattern& get_continuous_dof_pattern() const;

			void set_local_shape_function_set_id(LocalShapeFunctionSetID tsp);
			LocalShapeFunctionSetID get_local_shape_function_set_id() const;

			void set_dim(int dim);
			int get_dim() const;


			void set_subset_indices(std::vector<int> subsetIndices, int num_subsets);

			bool is_def_in_subset(int s) const;

			template <typename TGeomObj>
			inline size_type get_multi_indices_of_geom_obj(TGeomObj* elem, multi_index_type ind[]);

			template <typename TElem>
			inline bool set_num_dofs(int subsetIndex, uint num);
			template <typename TElem>
			inline uint num_dofs(int subsetIndex) const;

			template <typename TElem>
			inline bool set_nr_dofgroup(int subsetIndex, uint dof, uint group);
			template <typename TElem>
			inline uint nr_dofgroup(int subsetIndex, uint dof) const;

			template <typename TElem>
			inline bool set_comp_in_dofgroup(int subsetIndex, uint dof, uint comp);
			template <typename TElem>
			inline uint comp_in_dofgroup(int subsetIndex, uint dof) const;


		private:
			subset_handler_type* m_sh;
			std::string m_name; 									// name of solution
			LocalShapeFunctionSetID m_fct_set_id; 							// associated trial space
			uint m_dim;												// space dimension of trial space
			ContinuousDoFPattern m_continuousDoFManager;					// cached elementdofpattern
			std::vector<int> m_subsetIndices; 						// subset indices where solution lives
			std::vector<DoFGroupMap> m_dofGroupMap[NUM_GEOMETRIC_BASE_OBJECTS]; 	// for each GeomObject: number of dofs, which DoFGroup and which component
		};
		friend class SingleDiscreteFunction;

	protected:
		// if locked, pattern can not be changed anymore
		bool m_bLocked;

		// true if GeomObject has DoFs somewhere at grid
		bool m_bGeomObjectHasDoF[NUM_GEOMETRIC_BASE_OBJECTS];

		// vector containing all Discrete Functions of this pattern (nr_fct refers to these)
		std::vector<SingleDiscreteFunction> m_SingleDiscreteFunctionVec;

		// DoFGroups for each GeomObject
		std::vector<DoFGroup<VertexBase> >  m_VertexBaseDoFGroups;
		std::vector<DoFGroup<EdgeBase> >    m_EdgeBaseDoFGroups;
		std::vector<DoFGroup<Face> >        m_FaceDoFGroups;
		std::vector<DoFGroup<Volume> >      m_VolumeDoFGroups;

		// number of Levels
		uint m_num_levels;
		std::vector<IndexInfo> m_IndexInfo;

		// associated SubsetHandler
		subset_handler_type& m_sh;
		int m_num_subsets;
};




}

#include "general_dof_manager_impl.h"


#endif /* __H__LIBDISCRETIZATION__DOF_MANAGER__DOF_MANAGER__ */
