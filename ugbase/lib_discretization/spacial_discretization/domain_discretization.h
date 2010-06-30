/*
 * domain_discretization.h
 *
 *  Created on: 29.06.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__DOMAIN_DISCRETIZATION__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__DOMAIN_DISCRETIZATION__

// other ug4 modules
#include "common/common.h"
#include "lib_algebra/lib_algebra.h"

// library intern headers
#include "lib_discretization/spacial_discretization/domain_discretization_interface.h"

namespace ug {

template <	typename TDiscreteFunction,
			typename TAlgebra = typename TDiscreteFunction::algebra_type >
class DomainDiscretization :
	public IDomainDiscretization<TDiscreteFunction, TAlgebra>
{
	protected:
		// type of algebra
		typedef TDiscreteFunction discrete_function_type;

		// type of algebra
		typedef TAlgebra algebra_type;

		// type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

		// type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

	public:
		DomainDiscretization()
		{};

		// Assemble routines for time independent problems
		IAssembleReturn assemble_jacobian_defect(matrix_type& J, vector_type& d, const discrete_function_type& u)
		{
			for(size_t i = 0; i < m_disc.size(); ++i)
				if((m_disc[i].disc)->assemble_jacobian_defect(J, d, u, m_disc[i].si) != IAssemble_OK)
					return IAssemble_ERROR;

			for(size_t i = 0; i < m_dirichletDisc.size(); ++i)
				if((m_dirichletDisc[i].disc)->clear_dirichlet_jacobian_defect(J, d, u, m_dirichletDisc[i].si) != IAssemble_OK)
					return IAssemble_ERROR;
			return IAssemble_OK;
		}
		IAssembleReturn assemble_jacobian(matrix_type& J, const discrete_function_type& u)
		{
			for(size_t i = 0; i < m_disc.size(); ++i)
				if((m_disc[i].disc)->assemble_jacobian(J, u, m_disc[i].si) != IAssemble_OK)
					return IAssemble_ERROR;

			for(size_t i = 0; i < m_dirichletDisc.size(); ++i)
				if((m_dirichletDisc[i].disc)->clear_dirichlet_jacobian(J, u, m_dirichletDisc[i].si) != IAssemble_OK)
					return IAssemble_ERROR;
			return IAssemble_OK;
		}
		IAssembleReturn assemble_defect(vector_type& d, const discrete_function_type& u)
		{
			for(size_t i = 0; i < m_disc.size(); ++i)
				if((m_disc[i].disc)->assemble_defect(d, u, m_disc[i].si) != IAssemble_OK)
					return IAssemble_ERROR;

			for(size_t i = 0; i < m_dirichletDisc.size(); ++i)
				if((m_dirichletDisc[i].disc)->clear_dirichlet_defect(d, u, m_dirichletDisc[i].si) != IAssemble_OK)
					return IAssemble_ERROR;
			return IAssemble_OK;
		}
		IAssembleReturn assemble_linear(matrix_type& mat, vector_type& rhs, const discrete_function_type& u)
		{
			for(size_t i = 0; i < m_disc.size(); ++i)
				if((m_disc[i].disc)->assemble_linear(mat, rhs, u, m_disc[i].si) != IAssemble_OK)
					return IAssemble_ERROR;

			for(size_t i = 0; i < m_dirichletDisc.size(); ++i)
				if((m_dirichletDisc[i].disc)->set_dirichlet_linear(mat, rhs, u, m_dirichletDisc[i].si) != IAssemble_OK)
					return IAssemble_ERROR;
			return IAssemble_OK;
		}
		IAssembleReturn assemble_solution(discrete_function_type& u)
		{
			for(size_t i = 0; i < m_dirichletDisc.size(); ++i)
				if((m_dirichletDisc[i].disc)->set_dirichlet_solution(u, m_dirichletDisc[i].si) != IAssemble_OK)
					return IAssemble_ERROR;
			return IAssemble_OK;
		}

		// Assemble routines for time dependent problems
		IAssembleReturn assemble_jacobian_defect(matrix_type& J, vector_type& d, const discrete_function_type& u, number time, number s_m, number s_a)
		{
			for(size_t i = 0; i < m_disc.size(); ++i)
				if((m_disc[i].disc)->assemble_jacobian_defect(J, d, u, m_disc[i].si, time, s_m, s_a) != IAssemble_OK)
					return IAssemble_ERROR;

			for(size_t i = 0; i < m_dirichletDisc.size(); ++i)
				if((m_dirichletDisc[i].disc)->clear_dirichlet_jacobian_defect(J, d, u, m_dirichletDisc[i].si, time) != IAssemble_OK)
					return IAssemble_ERROR;
			return IAssemble_OK;
		}

		IAssembleReturn assemble_jacobian(matrix_type& J, const discrete_function_type& u, number time, number s_m, number s_a)
		{
			for(size_t i = 0; i < m_disc.size(); ++i)
				if((m_disc[i].disc)->assemble_jacobian(J, u, m_disc[i].si, time, s_m, s_a) != IAssemble_OK)
					return IAssemble_ERROR;

			for(size_t i = 0; i < m_dirichletDisc.size(); ++i)
				if((m_dirichletDisc[i].disc)->clear_dirichlet_jacobian(J, u, m_dirichletDisc[i].si, time) != IAssemble_OK)
					return IAssemble_ERROR;
			return IAssemble_OK;
		}
		IAssembleReturn assemble_defect(vector_type& d, const discrete_function_type& u, number time, number s_m, number s_a)
		{
			for(size_t i = 0; i < m_disc.size(); ++i)
				if((m_disc[i].disc)->assemble_defect(d, u, m_disc[i].si, time, s_m, s_a) != IAssemble_OK)
					return IAssemble_ERROR;

			for(size_t i = 0; i < m_dirichletDisc.size(); ++i)
				if((m_dirichletDisc[i].disc)->clear_dirichlet_defect(d, u, m_dirichletDisc[i].si, time) != IAssemble_OK)
					return IAssemble_ERROR;
			return IAssemble_OK;
		}
		IAssembleReturn assemble_linear(matrix_type& mat, vector_type& rhs, const discrete_function_type& u, number time, number s_m, number s_a)
		{
			for(size_t i = 0; i < m_disc.size(); ++i)
				if((m_disc[i].disc)->assemble_linear(mat, rhs, u, m_disc[i].si, time, s_m, s_a) != IAssemble_OK)
					return IAssemble_ERROR;

			for(size_t i = 0; i < m_dirichletDisc.size(); ++i)
				if((m_dirichletDisc[i].disc)->set_dirichlet_linear(mat, rhs, u, m_dirichletDisc[i].si, time) != IAssemble_OK)
					return IAssemble_ERROR;
			return IAssemble_OK;
		}
		IAssembleReturn assemble_solution(discrete_function_type& u, number time)
		{
			for(size_t i = 0; i < m_dirichletDisc.size(); ++i)
				if((m_dirichletDisc[i].disc)->set_dirichlet_solution(u, m_dirichletDisc[i].si, time) != IAssemble_OK)
					return IAssemble_ERROR;
			return IAssemble_OK;
		}


	public:
		bool add_disc(IDimensionDomainDiscretization<discrete_function_type>& disc, int si)
		{
			m_disc.push_back(ElemDisc(si, &disc));
			return true;
		}

	protected:
		struct ElemDisc
		{
			ElemDisc(int s, IDimensionDomainDiscretization<discrete_function_type>* di) :
				si(s), disc(di) {};

			int si;
			IDimensionDomainDiscretization<discrete_function_type>* disc;
		};

		std::vector<ElemDisc> m_disc;

	public:
		bool add_dirichlet_bnd(IDirichletBoundaryValues<discrete_function_type>& bnd_disc, int si)
		{
			m_dirichletDisc.push_back(DirichletDisc(si, &bnd_disc));
			return true;
		}

	protected:
		struct DirichletDisc
		{
			DirichletDisc(int s, IDirichletBoundaryValues<discrete_function_type>* di) :
				si(s), disc(di) {};

			int si;
			IDirichletBoundaryValues<discrete_function_type>* disc;
		};

		std::vector<DirichletDisc> m_dirichletDisc;

	protected:
		// TODO: What is this function used for????
		virtual size_t num_fct() const
		{return m_disc[0].disc->num_fct();}

		virtual bool is_dirichlet(int si, size_t fct)
		{
			for(size_t i = 0; i < m_dirichletDisc.size(); ++i)
			{
				if(m_dirichletDisc[i].si != si) continue;
				if((m_dirichletDisc[i].disc)->is_dirichlet(fct) == true) return true;
			}
			return false;
		}
};

} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__DOMAIN_DISCRETIZATION__ */
