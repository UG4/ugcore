/*
 * grid_function_space.h
 *
 *  Created on: 19.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__FUNCTION_SPACE__GRID_FUNCTION_SPACE__
#define __H__LIBDISCRETIZATION__FUNCTION_SPACE__GRID_FUNCTION_SPACE__

#include "lib_discretization/dof_manager/p1conform_dof_manager/p1conform_dof_manager.h"

#ifdef UG_PARALLEL
	#include "lib_discretization/parallelization/parallelization.h"
#endif

namespace ug{

// A grid function brings approximation space and algebra together. For a given DoFManager and level, a grid function
// represents the solutions on the level 'level'
template <typename TApproximationSpace, typename TDoFManager, typename TAlgebra>
class GridFunction{
	public:
		// DOMAIN
		// domain type
		typedef typename TApproximationSpace::domain_type domain_type;

		// global coordinate type
		typedef typename domain_type::position_type position_type;

		// subset handler, where DoF Manager is defined
		typedef typename domain_type::grid_type grid_type;

		// subset handler, where DoF Manager is defined
		typedef typename domain_type::subset_handler_type subset_handler_type;


		// ALGEBRA
		// algebra type
		typedef TAlgebra algebra_type;

		// vector type used to store dof values
		typedef typename algebra_type::vector_type vector_type;

		// local vector type
		typedef typename vector_type::local_vector_type local_vector_type;

		// local index type
		typedef typename vector_type::local_index_type local_index_type;


		// APPROXIMATIONSPACE
		// Approximation Space
		typedef TApproximationSpace approximation_space_type;

		// DOFMANAGER
		// dof manager used for this approximation space
		typedef TDoFManager dof_manager_type;

	public:
		// Constructor
		GridFunction(std::string name, approximation_space_type& approximationSpace, dof_manager_type& dof_manager, uint level, bool allocate = true) :
			m_approximationSpace(approximationSpace), m_dof_manager(dof_manager), m_level(level), m_name(name)
		{
			UG_ASSERT(level < m_dof_manager.num_levels(), "Accessing level that does not exist");

			if(allocate)
			{
				// create storage for dofs
				m_dof_storage_vector = new vector_type;
				uint num_dofs = m_dof_manager.num_dofs(level);
				if(m_dof_storage_vector->create(num_dofs) != true)
					UG_ASSERT(0, "Cannot create vector memory for " << num_dofs << " dofs.\n");
			}
		};

		// creates a GridFunction with same pattern as given GridFunction, BUT values are not copied
		GridFunction(const GridFunction& v) :
			m_approximationSpace(v.m_approximationSpace), m_dof_manager(v.m_dof_manager), m_level(v.m_level), m_name(v.m_name)
		{
			// create storage for dofs
			m_dof_storage_vector = new vector_type;
			uint num_dofs = v.m_dof_storage_vector->size();
			if(m_dof_storage_vector->create(num_dofs) != true)
				UG_ASSERT(0, "Cannot create vector memory for " << num_dofs << " dofs.\n");
		};

		bool assign(const GridFunction& v)
		{
			UG_ASSERT(v.m_dof_storage_vector->size() == m_dof_storage_vector->size(), "Size of discrete function does not match. Cannot assign grid function.");
			*m_dof_storage_vector = *v.m_dof_storage_vector;
			return true;
		}

		template <typename TSurfaceDoFManager>
		bool project_surface(GridFunction<TApproximationSpace, TSurfaceDoFManager, TAlgebra>& v)
		{
			UG_DLOG(LIB_DISC_DISCRETE_FUNCTION, 1, "Projecting to surface : " << v.get_name() << " -> " << this->get_name() << ".\n");
			if(&m_approximationSpace != &v.m_approximationSpace)
			{
				UG_DLOG(LIB_DISC_DISCRETE_FUNCTION, 2, "Approximation Spaces are different. Aborting.\n");
				return false;
			}

			if(&m_dof_manager == &v.m_dof_manager)
			{
				UG_DLOG(LIB_DISC_DISCRETE_FUNCTION, 2, "Same DoF Manager used.\n");

				if(m_level == v.m_level)
				{
					UG_DLOG(LIB_DISC_DISCRETE_FUNCTION, 2, "Same level '" << m_level << "' of DoF Manager.\n");
					// Surface Level == finest grid level (full refinement)
					m_dof_storage_vector = v.m_dof_storage_vector;
				}
				else
				{
					UG_DLOG(LIB_DISC_DISCRETE_FUNCTION, 2, "Different levels of DoF Manager. Do nothing.\n");
				}
			}
			else
			{
				UG_DLOG(LIB_DISC_DISCRETE_FUNCTION, 2, "Different DoF Manager used\n");
				UG_ASSERT(0, "Not implemented.");
			}

			return true;
		}

		template <typename TSurfaceDoFManager>
		bool release_surface(GridFunction<TApproximationSpace, TSurfaceDoFManager, TAlgebra>& v)
		{
			UG_DLOG(LIB_DISC_DISCRETE_FUNCTION, 1, "Release link to surface : " << v.get_name() << " -> " << this->get_name() << ".\n");
			if(&m_approximationSpace != &v.m_approximationSpace)
			{
				UG_DLOG(LIB_DISC_DISCRETE_FUNCTION, 2, "Approximation Spaces are different. Aborting.\n");
				return false;
			}

			if(&m_dof_manager == &v.m_dof_manager)
			{
				UG_DLOG(LIB_DISC_DISCRETE_FUNCTION, 2, "Same DoF Manager used.\n");

				if(m_level == v.m_level)
				{
					UG_DLOG(LIB_DISC_DISCRETE_FUNCTION, 2, "Same level '" << m_level << "' of DoF Manager.\n");
					// forget about memory without deleting it.
					m_dof_storage_vector = NULL;
				}
			}
			else
			{
				UG_DLOG(LIB_DISC_DISCRETE_FUNCTION, 2, "Different DoF Manager used\n");
				UG_ASSERT(0, "Not implemented.");
			}

			return true;
		}


		// TODO: use smart pointers / or implement copy constructors etc.
		~GridFunction()
		{
			if(m_dof_storage_vector != NULL)
			{
				m_dof_storage_vector->destroy();
				delete m_dof_storage_vector;
			}
		}

		// norm used to measure convergence in solvers
		inline number norm() const
		{
			return m_dof_storage_vector->two_norm();
		}

		// number of dofs
		inline uint num_dofs() const
		{
			return m_dof_manager.num_dofs(m_level);
		}

		// number of functions
		inline uint num_fct() const
		{
			return m_dof_manager.num_fct();
		}

		// name of function
		inline std::string get_name(uint fct) const
		{
			return m_dof_manager.get_name(fct);
		}

		// dim of functions
		inline uint dim_fct(uint fct) const
		{
			return m_dof_manager.dim_fct(fct);
		}

		// where is function defined
		inline bool fct_is_def_in_subset(uint fct, int subsetIndex) const
		{
			return m_dof_manager.fct_is_def_in_subset(fct, subsetIndex);
		}

		// number of subsets
		inline int num_subsets() const
		{
			return m_dof_manager.num_subsets();
		}

		// iterator for elements where this grid function is defined
		template <typename TElem>
		inline typename geometry_traits<TElem>::iterator begin(int subsetIndex) const
		{
			return m_dof_manager.template begin<TElem>(subsetIndex, m_level);
		}

		// iterator for elements where this grid function is defined
		template <typename TElem>
		inline typename geometry_traits<TElem>::iterator end(int subsetIndex) const
		{
			return m_dof_manager.template end<TElem>(subsetIndex, m_level);
		}

		// number of elements of this type for a subset
		template <typename TElem>
		inline uint num(int subsetIndex) const
		{
			return m_dof_manager.template num<TElem>(subsetIndex, m_level);
		}

		// set all dofs on level 'level' to value 'w'
		bool set(number w)
		{
			return m_dof_storage_vector->set(w);
		}

		// evaluate Grid function on an element
		template <typename TElem>
		number evaluate(TElem* elem, unsigned int fct, unsigned int comp, const MathVector<reference_element_traits<TElem>::dim>& loc_pos) const;

		// evaluate Grid function for a global position
		number evaluate(unsigned int nr_fct, unsigned int comp, const position_type& glob_pos) const;

		// get dof values
		inline bool get_dof_values(local_vector_type& val, local_index_type& ind) const
		{
			m_dof_storage_vector->get(val, ind);
			return true;
		}

		// get multiindices on an finite element in canonical order
		template <typename TElem>
		inline size_t get_multi_indices(TElem* elem, unsigned int fct, local_index_type& ind, std::size_t offset = 0) const
		{
			return m_dof_manager.get_multi_indices(elem, fct, ind, offset);
		}

		// get multiindices on an geometric object in canonical order
		template <typename TGeomObj>
		inline size_t get_multi_indices_of_geom_obj(TGeomObj* elem, unsigned int fct, local_index_type& ind, std::size_t offset = 0) const
		{
			return m_dof_manager.get_multi_indices_of_geom_obj(elem, fct, ind, offset);
		}

		// get local dof values
		template <typename TElem>
		inline int get_dof_values(TElem* elem, unsigned int fct, local_vector_type& val) const
		{
			const int num_ind = m_dof_manager.num_multi_indices(elem, fct);
			local_index_type ind(num_ind);

			m_dof_manager.get_multi_indices(elem, fct, ind);

			m_dof_storage_vector->get(val, ind);

			return num_ind;
		}

		template <typename TGeomObj>
		inline int get_dof_values_of_geom_obj(TGeomObj* elem, unsigned int fct, local_vector_type& val) const
		{
			int num_ind = m_dof_manager.num_multi_indices(elem, fct);
			local_index_type ind(num_ind);

			// TODO: This is a quick hack. Needed: function num_multi_indices_of_geom_obj
			num_ind = m_dof_manager.get_multi_indices_of_geom_obj(elem, fct, ind);

			ind.resize(num_ind);

			m_dof_storage_vector->get(val, ind);

			UG_DEBUG_BEGIN(LIB_DISC_DISCRETE_FUNCTION, 2);
				UG_DLOG(LIB_DISC_DISCRETE_FUNCTION, 2, "\n " << num_ind << " Values read: ");
				for(uint i = 0; i < ind.size(); ++i)
					UG_DLOG(LIB_DISC_DISCRETE_FUNCTION, 2, "(" << ind[i] << ", " << val[i] << ") ");
			UG_DEBUG_END(LIB_DISC_DISCRETE_FUNCTION, 2);

			return num_ind;
		}

		template <typename TElem>
		inline bool set_dof_values(TElem* elem, unsigned int fct, local_vector_type& val)
		{
			const int num_ind = m_dof_manager.num_multi_indices(elem, fct);
			local_index_type ind(num_ind);

			m_dof_manager.get_multi_indices(elem, fct, ind);

			m_dof_storage_vector->set(val, ind);

			UG_DEBUG_BEGIN(LIB_DISC_DISCRETE_FUNCTION, 2);
				UG_DLOG(LIB_DISC_DISCRETE_FUNCTION, 2, "\n " << num_ind << " Values set: ");
				for(uint i = 0; i < ind.size(); ++i)
					UG_DLOG(LIB_DISC_DISCRETE_FUNCTION, 2, "(" << ind[i] << ", " << val[i] << ") ");
			UG_DEBUG_END(LIB_DISC_DISCRETE_FUNCTION, 2);

			return true;
		}

		template <typename TGeomObj>
		inline bool set_dof_values_of_geom_obj(TGeomObj* elem, unsigned int fct, local_vector_type& val)
		{
			int num_ind = m_dof_manager.num_multi_indices(elem, fct);
			local_index_type ind(num_ind);

			// TODO: This is a quick hack. Needed: function num_multi_indices_of_geom_obj
			num_ind = m_dof_manager.get_multi_indices_of_geom_obj(elem, fct, ind);

			ind.resize(num_ind);

			m_dof_storage_vector->set(val, ind);

			return true;
		}

		// approximation space of this gridfunction
		approximation_space_type& get_approximation_space()
		{
			return m_approximationSpace;
		}

		// returns the dofmanager used for this discrete function
		dof_manager_type& get_dof_manager()
		{
			return m_dof_manager;
		}

		// type of fundamental function
		LocalShapeFunctionSetID get_local_shape_function_set_id(uint fct) const
		{
			return m_dof_manager.get_local_shape_function_set_id(fct);
		}

		// returns the level of the dofmanager, that is used for this grid function
		uint get_level() const
		{
			return m_level;
		}

		// export the dof storage of this vector
		vector_type& get_vector()
		{
			return *m_dof_storage_vector;
		}

		// export the dof storage of this vector
		const vector_type& get_vector() const
		{
			return *m_dof_storage_vector;
		}

		// return the domain
		const domain_type& get_domain() const
		{
			return m_approximationSpace.get_domain();
		}

		// return the domain
		domain_type& get_domain()
		{
			return m_approximationSpace.get_domain();
		}

		// this function finalizes the dof pattern. Afterwards the pattern can only be changed by destroying the vector and creating a new one.
		bool finalize()
		{
			return m_dof_storage_vector->finalize();
		}

		GridFunction& operator+=(const GridFunction& v)
		{
			*m_dof_storage_vector += *(v.m_dof_storage_vector);
			return *this;
		}
		GridFunction& operator-=(const GridFunction& v)
		{
			*m_dof_storage_vector -= *(v.m_dof_storage_vector);
			return *this;
		}
		GridFunction& operator*=(number w)
		{
			*m_dof_storage_vector *= w;
			return *this;
		}

		std::string get_name()
		{
			return m_name;
		}

	///	testversion
	/**	If no communicator is passed, then the method will create its
	 *	own temporarily. This may expensive if the method is called
	 *	repeatedly.*/
		void parallel_additive_to_consistent(pcl::ParallelCommunicator<IndexLayout>* pCom = NULL)
		{
		#ifdef UG_PARALLEL
		//	create a new communicator if required.
			pcl::ParallelCommunicator<IndexLayout> tCom;
			if(!pCom)
				pCom = &tCom;
			pcl::ParallelCommunicator<IndexLayout>& com = *pCom;

		//	step 1: add slave values to master
		//	create the required communication policies
			ComPol_VecAdd<vector_type> cpVecAdd(m_dof_storage_vector);

		//	perform communication on the level
			com.send_data(m_dof_manager.get_slave_layout(m_level), cpVecAdd);
			com.receive_data(m_dof_manager.get_master_layout(m_level), cpVecAdd);
			com.communicate();

		//	step 2: copy master values to slaves
		//	create the required communication policies
			ComPol_VecCopy<vector_type> cpVecCopy(m_dof_storage_vector);

		//	perform communication on the level
			com.send_data(m_dof_manager.get_master_layout(m_level), cpVecCopy);
			com.receive_data(m_dof_manager.get_slave_layout(m_level), cpVecCopy);
			com.communicate();
		#endif
		}
		void parallel_additive_to_unique(pcl::ParallelCommunicator<IndexLayout>* pCom = NULL)
		{
		#ifdef UG_PARALLEL
		//	create a new communicator if required.
			pcl::ParallelCommunicator<IndexLayout> tCom;
			if(!pCom)
				pCom = &tCom;
			pcl::ParallelCommunicator<IndexLayout>& com = *pCom;

		//	step 1: add slave values to master and set slave values to zero
		//	create the required communication policies
			ComPol_VecAddSetZero<vector_type> cpVecAddSetZero(m_dof_storage_vector);

		//	perform communication on the level
			com.send_data(m_dof_manager.get_slave_layout(m_level), cpVecAddSetZero);
			com.receive_data(m_dof_manager.get_master_layout(m_level), cpVecAddSetZero);
			com.communicate();
		#endif
		}

	protected:
		// Approximation Space
		approximation_space_type& m_approximationSpace;

		// dof manager of this discrete function
		dof_manager_type& m_dof_manager;

		// level of Approximation Space, where this Grid Function lives
		uint m_level;

		// name
		std::string m_name;

		// vector storage, to store values of local degrees of freedom
		vector_type* m_dof_storage_vector;
};

template <typename TApproximationSpace, typename TDoFManager, typename TAlgebra>
inline std::ostream& operator<< (std::ostream& outStream, const GridFunction<TApproximationSpace, TDoFManager, TAlgebra>& v)
{
	outStream << v.get_vector();
	return outStream;
}


// Surface Grid function lives on the surface grid
template <typename TApproximationSpace, typename TAlgebra>
class SurfaceGridFunction : public GridFunction<TApproximationSpace, typename TApproximationSpace::surface_dof_manager_type, TAlgebra>{
	public:
		// Approximation Space
		typedef TApproximationSpace approximation_space_type;

		// DOFMANAGER
		// dof manager used for this approximation space
		typedef typename TApproximationSpace::surface_dof_manager_type dof_manager_type;

	public:
		// Constructor
		SurfaceGridFunction(approximation_space_type& approximationSpace, dof_manager_type& dof_manager) :
			GridFunction<TApproximationSpace, dof_manager_type, TAlgebra>(approximationSpace, dof_manager, dof_manager.num_levels() - 1)
		{};

		SurfaceGridFunction(const SurfaceGridFunction& v) :
			GridFunction<TApproximationSpace, dof_manager_type, TAlgebra>(v)
		{}
};

// Level Grid function lives on a level grid
template <typename TApproximationSpace, typename TAlgebra>
class LevelGridFunction : public GridFunction<TApproximationSpace, typename TApproximationSpace::level_dof_manager_type, TAlgebra>{
public:
	// Approximation Space
	typedef TApproximationSpace approximation_space_type;

	// DOFMANAGER
	// dof manager used for this approximation space
	typedef typename TApproximationSpace::level_dof_manager_type dof_manager_type;

	public:
		// Constructor
		LevelGridFunction(approximation_space_type& approximationSpace, dof_manager_type& dof_manager, uint level) :
			GridFunction<TApproximationSpace, dof_manager_type, TAlgebra>(approximationSpace, dof_manager, level)
		{};

		LevelGridFunction(const LevelGridFunction& v) :
			GridFunction<TApproximationSpace, dof_manager_type, TAlgebra>(v)
		{}


};

/**
 * This class describes a grid function space on a domain.
 * The Domain gives the partition of the Grid/Multigrid in terms of subsets.
 * Now the user can add fundamental Discrete functions on this subsets or unions of them.
 *
 * Once finalized, this function pattern is fixed. Internally dof indices are created.
 * Using this Approximation Space the user can create GridFunctions of the following types:
 *
 * - surface grid function = a grid function representing the space on the surface grid
 * - level grid function = a grid function representing the space on a level grid
 *  (NOTE: For a fully refined Multigrid a level grid covers the whole domain. However
 *         for a locally/adaptivly refined MultiGrid the level grid solution is only
 *         living on a part of the domain)
 */
struct nonadaptive{};
struct adaptive{};

template <typename TDomain, typename TAlgebra, template <class TDomain> class TDoFManager, typename TCategory = typename TDomain::categroy_type>
class ApproximationSpace;


// non adaptive Approximation space
template <typename TDomain, typename TAlgebra, template <class TDomain> class TDoFManager>
class ApproximationSpace<TDomain, TAlgebra, TDoFManager, nonadaptive>{
	private:
		// to make it more readable
		typedef ApproximationSpace<TDomain, TAlgebra, TDoFManager, nonadaptive> approximation_space_type;

	public:
		// DOMAIN
		// global coordinate type
		typedef TDomain domain_type;

		// global coordinate type
		typedef typename domain_type::position_type position_type;

		// subset handler, where DoF Manager is defined
		typedef typename domain_type::grid_type grid_type;

		// subset handler, where DoF Manager is defined
		typedef typename domain_type::subset_handler_type level_subset_handler_type;
		typedef typename domain_type::subset_handler_type surface_subset_handler_type;


		// GRID FUNCTION of this factory
		typedef TAlgebra algebra_type;

		// dof manager used
		#ifdef UG_PARALLEL
			typedef ParallelDoFManager<TDoFManager<level_subset_handler_type> >
				level_dof_manager_type;
			typedef ParallelDoFManager<TDoFManager<surface_subset_handler_type> >
				surface_dof_manager_type;
		#else
			typedef TDoFManager<level_subset_handler_type> level_dof_manager_type;
			typedef TDoFManager<surface_subset_handler_type> surface_dof_manager_type;
		#endif

		typedef GridFunction<approximation_space_type, level_dof_manager_type, algebra_type> level_function_type;
		typedef GridFunction<approximation_space_type, surface_dof_manager_type, algebra_type> surface_function_type;
		typedef surface_function_type function_type;

	public:
		ApproximationSpace(std::string name, domain_type& domain) :
			m_name(name), m_domain(domain),
			m_subset_handler(domain.get_subset_handler()), m_level_dof_manager(m_subset_handler)
		{
			#ifdef UG_PARALLEL
				m_level_dof_manager.set_grid_layout_map(
					domain.get_distributed_grid_manager()->
					grid_layout_map());
			#endif
		};


		/// add a single solution of LocalShapeFunctionSetID to selected subsets in dimension 'dim'
		/**
		 * \param[in] name			Name of this Single Solution
		 * \param[in] TrialSpace	Trial Space for this function
		 * \param[in] SubsetIndecex	Std::Vector of subset indeces, where this solution lives
		 */
		bool add_fundamental_discrete_function(std::string name, LocalShapeFunctionSetID id, int dim, std::vector<int>& SubsetIndices)
		{
			if(m_level_dof_manager.add_discrete_function(name, id, SubsetIndices, dim) != true)
			{
				UG_LOG("ERROR in add_fundamental_discrete_function: cannot add discret function to level dof manager.\n");
				return false;
			}
			return true;
		}

		/// add a single solution of LocalShapeFunctionSetID to selected subsets for the world dimension
		bool add_fundamental_discrete_function(std::string name, LocalShapeFunctionSetID id, std::vector<int>& SubsetIndices)
		{
			const int dim = TDomain::dim;
			return this->add_fundamental_discrete_function(name, id, dim, SubsetIndices);
		}

		// add function on all domains for the world dimension
		bool add_fundamental_discrete_function(std::string name, LocalShapeFunctionSetID id)
		{
			const int dim = TDomain::dim;
			if(m_level_dof_manager.add_discrete_function(name, id, dim) != true)
			{
				UG_LOG("ERROR in add_fundamental_discrete_function: cannot add discret function to level dof manager.\n");
				return false;
			}
			return true;
		}

		// finish function pattern
		bool finalize()
		{
			if(m_level_dof_manager.finalize() != true)
			{
				UG_LOG("ERROR in finalize(): cannot finalize level dof manager.\n");
				return false;
			}
			return true;
		}

		// adaptive domain used or not ?
		inline bool is_adaptive()
		{
			return true;
		}

		// number of fundamental discrete functions
		inline uint num_fct() const
		{
			return m_level_dof_manager.num_fct();
		}

		// number of subsets
		inline int num_subsets() const
		{
			return m_subset_handler.num_subsets();
		}

		// type of fundamental function
		LocalShapeFunctionSetID get_local_shape_function_set_id(uint fct) const
		{
			return m_level_dof_manager.get_local_shape_function_set_id(fct);
		}

		// where is function defined
		inline bool fct_is_def_in_subset(uint fct, int subsetIndex) const
		{
			return m_level_dof_manager.fct_is_def_in_subset(fct, subsetIndex);
		}

		// number of levels
		uint num_levels()
		{
			return m_level_dof_manager.num_levels();
		}

		// return the domain
		const domain_type& get_domain() const
		{
			return m_domain;
		}

		// return the domain
		domain_type& get_domain()
		{
			return m_domain;
		}

		// create a new grid function of this approximation space
		level_function_type* create_level_grid_function(std::string name, uint level, bool allocate = true)
		{
			if(m_level_dof_manager.is_locked() != true)
			{
				std::cout << "Approximation Space has not been locked. Can not create a grid function of an unlocked space.\n";
				return NULL;
			}

			level_function_type* gridFct = new level_function_type(name, *this, m_level_dof_manager, level, allocate);
			return gridFct;
		}

		// create a new grid function of this approximation space
		surface_function_type* create_surface_grid_function(std::string name, bool allocate = true)
		{
			if(m_level_dof_manager.is_locked() != true)
			{
				std::cout << "Approximation Space has not been locked. Can not create a grid function of an unlocked space.\n";
				return NULL;
			}

			surface_function_type* gridFct = new surface_function_type(name, *this, m_level_dof_manager, m_level_dof_manager.num_levels() - 1, allocate);
			return gridFct;
		}

	protected:
		// name of this Approximation Space
		std::string m_name;

		// Domain, where solution lives
		domain_type& m_domain;

		// grid or multigrid or subsethandler, where elements are stored
		typename domain_type::subset_handler_type& m_subset_handler;

		// dof manager used for this Approximation Space
		level_dof_manager_type m_level_dof_manager;
};


}


#endif /* __H__LIBDISCRETIZATION__FUNCTION_SPACE__GRID_FUNCTION_SPACE__ */
