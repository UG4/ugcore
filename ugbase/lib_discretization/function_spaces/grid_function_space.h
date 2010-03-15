/*
 * grid_function_space.h
 *
 *  Created on: 19.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__FUNCTION_SPACE__FUNCTION_SPACE__
#define __H__LIBDISCRETIZATION__FUNCTION_SPACE__FUNCTION_SPACE__

#include "lib_discretization/dof_manager/p1conform_dof_manager/p1conform_dof_manager.h"

namespace ug{

template <typename TDomain, typename TAlgebra, typename TDoFManager>
class DiscreteGridFunction{
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

	public:
		// Constructor
		DiscreteGridFunction(dof_manager_type& dof_manager) : m_dof_manager(dof_manager) {};

		template <typename TElem>
		number evaluate(TElem* elem, unsigned int fct, unsigned int comp, const MathVector<reference_element_traits<TElem>::dim>& loc_pos) const;

		number evaluate(unsigned int nr_fct, unsigned int comp, const position_type& glob_pos) const;

		inline bool get_dof_values(local_vector_type& val, local_index_type& ind)
		{
			m_dof_storage_vector.get(val, ind);
			return true;
		}

		template <typename TElem>
		inline uint get_multi_indices(TElem* elem, unsigned int fct, unsigned int comp, local_index_type& ind)
		{
			return m_dof_manager.get_multi_indices(elem, fct, comp, ind);
		}

		template <typename TElem>
		inline int get_dof_values(TElem* elem, unsigned int fct, unsigned int comp, local_vector_type& val)
		{
			const int num_ind = get_num_multi_indices(elem, fct, comp);
			local_index_type ind(num_ind);

			m_dof_manager.get_multi_indices(elem, fct, comp, ind);

			m_dof_storage_vector.get(val, ind);

			return num_ind;
		}

		template <typename TElem>
		inline bool set_dof_values(TElem* elem, unsigned int fct, unsigned int comp, local_vector_type& val)
		{
			const int num_ind = get_num_multi_indices(elem, fct, comp);
			local_index_type ind(num_ind);

			m_dof_manager.get_multi_indices(elem, fct, comp, ind);

			m_dof_storage_vector.get(val, ind);

			return num_ind;
		}

		vector_type& get_vector(uint level = 0)
		{
			return m_dof_storage_vector;
		}

	protected:
		// dof manager of this discrete function
		dof_manager_type& m_dof_manager;

		// vector storage, to store values of local degrees of freedom
		vector_type m_dof_storage_vector;
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
		typedef typename TDomain::position_type position_type;

		// subset handler, where DoF Manager is defined
		typedef typename TDomain::grid_type grid_type;

		// subset handler, where DoF Manager is defined
		typedef typename TDomain::subset_handler_type subset_handler_type;

		// vector type used to store dof values
		typedef typename TAlgebra::vector_type vector_type;

		// dof manager used for this approximation space
		typedef TDoFManager dof_manager_type;

		// element_container
		typedef typename TDoFManager::element_container_type element_container_type;

		// index type
		typedef typename dof_manager_type::multi_index_type multi_index_type;

		// type of Grid functions, that are contained in this Approximation Space
		typedef DiscreteGridFunction<TDomain, TAlgebra, TDoFManager> grid_function_type;

	public:
		ApproximationSpace(std::string name, element_container_type& sh) : m_name(name), m_dof_manager(sh), m_lock(false) {};


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

		// create a new grid function of this approximation space
		grid_function_type* create_grid_function(std::string name)
		{
			if(m_lock != true)
			{
				std::cout << "Approximation Space has not been locked. Can not create a grid function of an unlocked space.\n";
				return NULL;
			}
			grid_function_type* gridFct = new grid_function_type(m_dof_manager);
			return gridFct;
		}

		dof_manager_type& get_dof_manager()
		{
			return m_dof_manager;
		}

	protected:
		// after lock is set, Approximation space can not be changed
		bool m_lock;

		// name of this Approximation Space
		std::string m_name;

		// dof manager used for this Approximation Space
		dof_manager_type m_dof_manager;
};


}


#endif /* __H__LIBDISCRETIZATION__FUNCTION_SPACE__FUNCTION_SPACE__ */
