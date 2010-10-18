/*
 * dirichlet_bnd.h
 *
 *  Created on: 04.08.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__POST_PROCESS__DIRICHLET_BOUNDARY__DIRICHLET_BND__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__POST_PROCESS__DIRICHLET_BOUNDARY__DIRICHLET_BND__

#include "../post_process_interface.h"

namespace ug{

template<typename TDomain, typename TAlgebra>
class DirichletBoundary : public IDirichletPostProcess<TAlgebra> {
	public:
		// type of domain
		typedef TDomain domain_type;

		// type of position coordinates (e.g. position_type)
		typedef typename domain_type::position_type position_type;

		// type of algebra
		typedef TAlgebra algebra_type;

		// type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

		// type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

		// local matrix type
		typedef LocalMatrix<typename TAlgebra::matrix_type::entry_type> local_matrix_type;

		// local vector type
		typedef LocalVector<typename TAlgebra::vector_type::entry_type> local_vector_type;

		// local index type
		typedef LocalIndices local_index_type;

	public:
		typedef bool (*Boundary_fct)(number&, const position_type&, number);

		DirichletBoundary(domain_type& domain) :
		  m_domain(domain)
		  {
			register_assemble_functions();
			m_vFctId.clear(); m_vBNDFct.clear();
		  }

		void add(size_t fct_id, Boundary_fct BNDFunc)
		{
			m_vFctId.push_back(fct_id);
			m_vBNDFct.push_back(BNDFunc);
		}

	public:
		template <typename TElem>
		inline bool prepare_element_loop()
		{
			// remember position attachement
			m_aaPos = m_domain.get_position_accessor();
			m_vCorner.resize(1);
			return true;
		}

		template <typename TElem>
		inline bool prepare_element(TElem* elem)
		{
			m_vCorner[0] = m_aaPos[elem];
			return true;
		}

		template <typename TElem>
		inline bool finish_element_loop()
		{
			return true;
		}

		template <typename TElem>
		inline bool post_process_J(matrix_type& J, const local_index_type& ind, number time=0.0)
		{
			number val;
			for(size_t i = 0; i < m_vBNDFct.size(); ++i)
			{
				const size_t fct = m_vFctId[i];
				for(size_t dof = 0; dof < ind.num_dofs(fct); ++dof)
				{
					if(m_vBNDFct[i](val, m_vCorner[0], time))
					{
						SetDirichletRow(J, ind.global_index(fct, dof), ind.comp(fct, dof));
					}
				}
			}
			return true;
		}

		// TODO: This is currently only for vertices
		template <typename TElem>
		inline bool post_process_d(vector_type& d, const local_index_type& ind, number time=0.0)
		{
			number val;
			for(size_t i = 0; i < m_vBNDFct.size(); ++i)
			{
				const size_t fct = m_vFctId[i];
				for(size_t dof = 0; dof < ind.num_dofs(fct); ++dof)
				{
					if(m_vBNDFct[i](val, m_vCorner[0], time))
					{
						BlockRef(d[ind.global_index(fct, dof)], ind.comp(fct, dof)) = 0.0;
					}
				}
			}
			return true;
		}

		template <typename TElem>
		inline bool post_process_f(local_vector_type& rhs, const local_index_type& ind, number time=0.0)
		{
			number val;
			for(size_t i = 0; i < m_vBNDFct.size(); ++i)
			{
				const size_t fct = m_vFctId[i];
				for(size_t dof = 0; dof < ind.num_dofs(fct); ++dof)
				{
					if(m_vBNDFct[i](val, m_vCorner[0], time))
					{
						BlockRef(rhs[ind.global_index(fct, dof)], ind.comp(fct, dof)) = val;
					}
				}
			}
			return true;
		}

		template <typename TElem>
		inline bool set_solution(vector_type& x, const local_index_type& ind, number time=0.0)
		{
			number val;
			for(size_t i = 0; i < m_vBNDFct.size(); ++i)
			{
				const size_t fct = m_vFctId[i];
				for(size_t dof = 0; dof < ind.num_dofs(fct); ++dof)
				{
					if(m_vBNDFct[i](val, m_vCorner[0], time))
					{
						BlockRef(x[ind.global_index(fct, dof)], ind.comp(fct, dof)) = val;
					}
				}
			}
			return true;
		}

	private:
		///////////////////////////////////////
		// registering for reference elements
		///////////////////////////////////////
		template <int dim> class numType{};

		void register_assemble_functions()
		{
			numType<TDomain::dim> dummy;
			register_assemble_functions(dummy);
		}

		// register for 1D
		void register_assemble_functions(numType<1> dummy)
		{
			register_all_assemble_functions<VertexBase>(ROID_VERTEX);
		}

		// register for 2D
		void register_assemble_functions(numType<2> dummy)
		{
			register_all_assemble_functions<VertexBase>(ROID_VERTEX);
			// TODO: register more elements
		}

		// register for 3D
		void register_assemble_functions(numType<3> dummy)
		{
			register_all_assemble_functions<VertexBase>(ROID_VERTEX);
			// TODO: register more elements
		}

		// help function
		template <typename TElem>
		void register_all_assemble_functions(int id)
		{
			register_prepare_element_loop_function(	id, &DirichletBoundary::template prepare_element_loop<TElem>);
			register_prepare_element_function(		id, &DirichletBoundary::template prepare_element<TElem>);
			register_finish_element_loop_function(	id, &DirichletBoundary::template finish_element_loop<TElem>);
			register_post_process_J_function(		id, &DirichletBoundary::template post_process_J<TElem>);
			register_post_process_d_function(		id, &DirichletBoundary::template post_process_d<TElem>);
			register_post_process_f_function(		id, &DirichletBoundary::template post_process_f<TElem>);
			register_set_solution_function(			id, &DirichletBoundary::template set_solution<TElem>);
		}

	protected:
		std::vector<MultiIndex<2> > m_multInd;
		std::vector<position_type> m_vCorner;

		std::vector<Boundary_fct> m_vBNDFct;
		std::vector<size_t> m_vFctId;

		typename TDomain::position_accessor_type m_aaPos;

		size_t m_fct;
		domain_type& m_domain;
};

} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__POST_PROCESS__DIRICHLET_BOUNDARY__DIRICHLET_BND__ */
