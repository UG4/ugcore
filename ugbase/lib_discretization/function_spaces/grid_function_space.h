/*
 * grid_function_space.h
 *
 *  Created on: 19.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__FUNCTION_SPACE__GRID_FUNCTION_SPACE__
#define __H__LIBDISCRETIZATION__FUNCTION_SPACE__GRID_FUNCTION_SPACE__

#include "lib_discretization/dof_manager/p1conform_dof_manager/p1conform_dof_manager.h"

namespace ug{
// predeclaration
template <typename TDomain, typename TAlgebra, typename TDoFManager>
class ApproximationSpace;


template <typename TDomain, typename TAlgebra, typename TDoFManager>
class DiscreteMultiGridFunction{
	public:
		// domain type
		typedef TDomain domain_type;

		// global coordinate type
		typedef typename TDomain::position_type position_type;

		// subset handler, where DoF Manager is defined
		typedef typename TDomain::grid_type grid_type;

		// subset handler, where DoF Manager is defined
		typedef typename TDomain::subset_handler_type subset_handler_type;

		// algebra type
		typedef TAlgebra algebra_type;

		// vector type used to store dof values
		typedef typename TAlgebra::vector_type vector_type;

		// local vector type
		typedef typename vector_type::local_vector_type local_vector_type;

		// local index type
		typedef typename vector_type::local_index_type local_index_type;


		// dof manager used for this approximation space
		typedef TDoFManager dof_manager_type;

		// element_container
		typedef typename TDoFManager::element_container_type element_container_type;

		// index type
		typedef typename dof_manager_type::multi_index_type multi_index_type;

		// Approximation Space
		typedef ApproximationSpace<TDomain, TAlgebra, TDoFManager> approximation_space_type;

	public:
		// Constructor
		DiscreteMultiGridFunction(domain_type& domain, approximation_space_type& approximationSpace, dof_manager_type& dof_manager) :
			m_domain(domain), m_dof_manager(dof_manager), m_approximationSpace(approximationSpace)
		{
			int num_levels = m_dof_manager.num_levels();

			m_dof_storage_vector.resize(num_levels);
			for(std::size_t l = 0; l < m_dof_storage_vector.size(); ++l)
			{
				m_dof_storage_vector[l] = new vector_type();
				UG_ASSERT(m_dof_storage_vector[l] != NULL, "Cannot allocate dof vector.");
				uint num_dofs = m_dof_manager.num_dofs(l);
				bool b = m_dof_storage_vector[l]->create(num_dofs);
				UG_ASSERT(b == true, "Cannot create vector memory for " << num_dofs << " dofs.\n");
			}
		};

		DiscreteMultiGridFunction(const DiscreteMultiGridFunction& v) :
			m_domain(v.m_domain), m_dof_manager(v.m_dof_manager), m_approximationSpace(v.m_approximationSpace)
		{
			int num_levels = v.m_dof_storage_vector.size();

			m_dof_storage_vector.resize(num_levels);
			for(std::size_t l = 0; l < m_dof_storage_vector.size(); ++l)
			{
				m_dof_storage_vector[l] = new vector_type();
				UG_ASSERT(m_dof_storage_vector[l] != NULL, "Cannot allocate dof vector.");
				uint num_dofs = v.m_dof_storage_vector[l].size();
				bool b = m_dof_storage_vector[l]->create(num_dofs);
				UG_ASSERT(b == true, "Cannot create vector memory for " << num_dofs << " dofs.\n");
			}
		};

		// TODO: use smart pointers / or implement copy constructors etc.
		~DiscreteMultiGridFunction()
		{
			for(std::size_t i = 0; i < m_dof_storage_vector.size(); ++i)
			{
				if(m_dof_storage_vector[i] != NULL)
				{
					m_dof_storage_vector[i].destroy();
					delete m_dof_storage_vector[i];
				}
			}
		}

		inline uint dim_fct(uint fct)
		{
			return m_dof_manager.dim_fct(fct);
		}

		inline uint num_fct()
		{
			return m_dof_manager.num_fct();
		}

		inline std::string get_name(uint fct)
		{
			return m_dof_manager.get_name(fct);
		}

		inline bool fct_is_def_in_subset(uint fct, int subsetIndex)
		{
			return m_dof_manager.fct_is_def_in_subset(fct, subsetIndex);
		}

		inline uint num_levels()
		{
			return m_dof_manager.num_levels();
		}

		inline int num_subsets()
		{
			return m_dof_manager.num_subsets();
		}

		template <typename TElem>
		inline typename geometry_traits<TElem>::iterator begin(uint level, int subsetIndex) const
		{
			return m_dof_manager.template begin<TElem>(level, subsetIndex);
		}

		template <typename TElem>
		inline typename geometry_traits<TElem>::iterator end(uint level, int subsetIndex) const
		{
			return m_dof_manager.template end<TElem>(level, subsetIndex);
		}

		// set all dofs on level 'level' to value 'w'
		bool set(number w, uint level)
		{
			return m_dof_storage_vector[level]->set(w);
		}

		// set all dofs to value 'w'
		bool set(number w)
		{
			for(std::size_t i = 0; i < m_dof_storage_vector.size(); ++i)
			{
				if(this->set(w, i)!=true) return false;
			}
			return true;
		}

		template <typename TElem>
		number evaluate(TElem* elem, unsigned int fct, unsigned int comp, const MathVector<reference_element_traits<TElem>::dim>& loc_pos) const;

		number evaluate(unsigned int nr_fct, unsigned int comp, const position_type& glob_pos) const;

		inline bool get_dof_values(uint level, local_vector_type& val, local_index_type& ind)
		{
			m_dof_storage_vector[level]->get(val, ind);
			return true;
		}

		template <typename TElem>
		inline uint get_multi_indices(TElem* elem, unsigned int fct, local_index_type& ind)
		{
			return m_dof_manager.get_multi_indices(elem, fct, ind);
		}

		template <typename TElem>
		inline uint get_multi_indices_of_geom_obj(TElem* elem, unsigned int fct, local_index_type& ind)
		{
			return m_dof_manager.get_multi_indices_of_geom_obj(elem, fct, ind);
		}

		template <typename TElem>
		inline int get_dof_values(TElem* elem, unsigned int fct, local_vector_type& val)
		{
			const int num_ind = m_dof_manager.num_multi_indices(elem, fct);
			local_index_type ind(num_ind);

			// TODO: avoid double request of level: here and in get_multi_index (also below ... )
			const uint level = m_dof_manager.get_grid().get_level(elem);
			m_dof_manager.get_multi_indices(elem, fct, ind);

			m_dof_storage_vector[level].get(val, ind);

			return num_ind;
		}

		template <typename TElem>
		inline int get_dof_values_of_geom_obj(TElem* elem, unsigned int fct, local_vector_type& val)
		{
			int num_ind = m_dof_manager.num_multi_indices(elem, fct);
			local_index_type ind(num_ind);

			// TODO: This is a quick hack. Needed: function num_multi_indices_of_geom_obj
			num_ind = m_dof_manager.get_multi_indices_of_geom_obj(elem, fct, ind);

			ind.resize(num_ind);

			const uint level = m_dof_manager.get_grid().get_level(elem);
			m_dof_storage_vector[level]->get(val, ind);

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

			const uint level = m_dof_manager.get_grid().get_level(elem);
			m_dof_storage_vector[level]->set(val, ind);

			UG_DEBUG_BEGIN(LIB_DISC_DISCRETE_FUNCTION, 2);
				UG_DLOG(LIB_DISC_DISCRETE_FUNCTION, 2, "\n " << num_ind << " Values set: ");
				for(uint i = 0; i < ind.size(); ++i)
					UG_DLOG(LIB_DISC_DISCRETE_FUNCTION, 2, "(" << ind[i] << ", " << val[i] << ") ");
			UG_DEBUG_END(LIB_DISC_DISCRETE_FUNCTION, 2);

			return true;
		}

		template <typename TElem>
		inline bool set_dof_values_of_geom_obj(TElem* elem, unsigned int fct, local_vector_type& val)
		{
			int num_ind = m_dof_manager.num_multi_indices(elem, fct);
			local_index_type ind(num_ind);

			// TODO: This is a quick hack. Needed: function num_multi_indices_of_geom_obj
			num_ind = m_dof_manager.get_multi_indices_of_geom_obj(elem, fct, ind);

			ind.resize(num_ind);

			const uint level = m_dof_manager.get_grid().get_level(elem);
			m_dof_storage_vector[level]->set(val, ind);

			return true;
		}

		domain_type& get_domain()
		{
			return m_domain;
		}

		approximation_space_type& get_approximation_space()
		{
			return m_approximationSpace;
		}

		vector_type& get_vector(uint level)
		{
			UG_ASSERT(level < m_dof_storage_vector.size(), "Requested level does not exist in DiscreteMultiGridFunction.");
			return *m_dof_storage_vector[level];
		}

		const vector_type& get_vector(uint level) const
		{
			UG_ASSERT(level < m_dof_storage_vector.size(), "Requested level does not exist in DiscreteMultiGridFunction.");
			return *m_dof_storage_vector[level];
		}

		// this function finalizes the dof pattern. Afterwards the pattern can only be changed by destroying the vector and creating a new one.
		bool finalize()
		{
			bool b = true;
			for(std::size_t i = 0; i < m_dof_storage_vector.size(); ++i)
			{
				b = b && m_dof_storage_vector[i]->finalize();
			}
			return b;
		}

	protected:
		// domain
		domain_type& m_domain;

		// dof manager of this discrete function
		dof_manager_type& m_dof_manager;

		// Approximation Space
		approximation_space_type& m_approximationSpace;

		// vector storage, to store values of local degrees of freedom
		std::vector<vector_type*> m_dof_storage_vector;
};


/**
 * This class describes a function space on a grid.
 * Still, the dimension of the function space may be infinite.
 * Discrete subspaces of a GridFunctionSpace can be implemented
 * by a ApproximationSpace
 */
template <typename TDomain, typename TAlgebra, typename TDoFManager>
class ApproximationSpace{
	public:
		// global coordinate type
		typedef TDomain domain_type;

		// global coordinate type
		typedef typename TDomain::position_type position_type;

		// subset handler, where DoF Manager is defined
		typedef typename TDomain::grid_type grid_type;

		// subset handler, where DoF Manager is defined
		typedef typename TDomain::subset_handler_type subset_handler_type;

		// vector type used to store dof values
		typedef typename TAlgebra::vector_type vector_type;

		// vector type used to store dof values
		typedef typename vector_type::local_index_type local_index_type;

		// dof manager used for this approximation space
		typedef TDoFManager dof_manager_type;

		// element_container
		typedef typename TDoFManager::element_container_type element_container_type;

		// index type
		typedef typename dof_manager_type::multi_index_type multi_index_type;

		// type of Grid functions, that are contained in this Approximation Space
		typedef DiscreteMultiGridFunction<TDomain, TAlgebra, TDoFManager> function_type;

	public:
		ApproximationSpace(std::string name, domain_type& domain, element_container_type& sh) :
			m_lock(false), m_name(name), m_domain(domain), m_elem_container(sh), m_dof_manager(sh)
		{
			// TODO: Check that correct sh == domain.sh
		};


		/// add a single solution of LocalShapeFunctionSetID to selected subsets
		/**
		 * \param[in] name			Name of this Single Solution
		 * \param[in] TrialSpace	Trial Space for this function
		 * \param[in] SubsetIndecex	Std::Vector of subset indeces, where this solution lives
		 */
		bool add_fundamental_discrete_function(std::string name, LocalShapeFunctionSetID id, int dim, std::vector<int>& SubsetIndices)
		{
			return m_dof_manager.add_discrete_function(name, id, SubsetIndices, dim);
		}

		bool add_fundamental_discrete_function(std::string name, LocalShapeFunctionSetID id)
		{
			const int dim = TDomain::dim;
			return m_dof_manager.add_discrete_function(name, id, dim);
		}

		bool finalize()
		{
			if(m_dof_manager.finalize() != true) return false;
			else m_lock = true;
			return true;
		}

		uint num_dofs(uint level)
		{
			return m_dof_manager.num_dofs(level);
		}

		uint num_levels()
		{
			return m_dof_manager.num_levels();
		}

		// create a new grid function of this approximation space
		function_type* create_grid_function(std::string name)
		{
			if(m_lock != true)
			{
				std::cout << "Approximation Space has not been locked. Can not create a grid function of an unlocked space.\n";
				return NULL;
			}
			function_type* gridFct = new function_type(m_domain, *this, m_dof_manager);
			return gridFct;
		}

		LocalShapeFunctionSetID get_local_shape_function_set_id(uint fct)
		{
			return m_dof_manager.get_local_shape_function_set_id(fct);
		}

		template <typename TElem>
		inline uint get_multi_indices_of_geom_obj(TElem* elem, unsigned int fct, local_index_type& ind)
		{
			return m_dof_manager.get_multi_indices_of_geom_obj(elem, fct, ind);
		}

		inline uint num_fct()
		{
			return m_dof_manager.num_fct();
		}

		inline int num_subsets()
		{
			return m_dof_manager.num_subsets();
		}

		inline bool fct_is_def_in_subset(uint fct, int subsetIndex)
		{
			return m_dof_manager.fct_is_def_in_subset(fct, subsetIndex);
		}

		dof_manager_type& get_dof_manager()
		{
			return m_dof_manager;
		}

		domain_type& get_domain()
		{
			return m_domain;
		}

		element_container_type& get_element_container()
		{
			return m_elem_container;
		}

		template <typename TElem>
		inline typename geometry_traits<TElem>::iterator begin(uint level, int subsetIndex) const
		{
			return m_dof_manager.template begin<TElem>(level, subsetIndex);
		}

		template <typename TElem>
		inline typename geometry_traits<TElem>::iterator end(uint level, int subsetIndex) const
		{
			return m_dof_manager.template end<TElem>(level, subsetIndex);
		}


	protected:
		// after lock is set, Approximation space can not be changed
		bool m_lock;

		// name of this Approximation Space
		std::string m_name;

		// Domain, where solution lives
		domain_type& m_domain;

		// grid or multigrid or subsethandler, where elements are stored
		element_container_type& m_elem_container;

		// dof manager used for this Approximation Space
		dof_manager_type m_dof_manager;
};


}


#endif /* __H__LIBDISCRETIZATION__FUNCTION_SPACE__GRID_FUNCTION_SPACE__ */
