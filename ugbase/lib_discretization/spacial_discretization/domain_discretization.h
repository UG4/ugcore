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
#include "./domain_discretization_interface.h"
#include "./elem_disc/elem_disc_assemble_util.h"
#include "./coupled_elem_disc/coupled_elem_disc_assemble_util.h"
#include "./post_process/dirichlet_boundary/subset_dirichlet_post_process_util.h"
#include "./subset_assemble_util.h"
#include "lib_discretization/common/function_group.h"

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

		IAssembleReturn assemble_jacobian(matrix_type& J, const discrete_function_type& u)
		{
			if(!check_solution(u)) return IAssemble_ERROR;

			for(size_t i = 0; i < m_vElemDisc.size(); ++i)
				if(!AssembleJacobian<IElemDisc<TAlgebra>,TDiscreteFunction,TAlgebra>
					(*m_vElemDisc[i].disc, J, u, m_vElemDisc[i].fcts, m_vElemDisc[i].si, m_vElemDisc[i].dim))
						return IAssemble_ERROR;

			for(size_t i = 0; i < m_vCoupledElemDisc.size(); ++i)
				if(!AssembleJacobian<CoupledSystem<TDiscreteFunction, TAlgebra>,TDiscreteFunction,TAlgebra>
					(*m_vCoupledElemDisc[i].disc, J, u, m_vCoupledElemDisc[i].fcts, m_vCoupledElemDisc[i].si, m_vCoupledElemDisc[i].dim))
						return IAssemble_ERROR;

			for(size_t i = 0; i < m_vDirichletDisc.size(); ++i)
				if((m_vDirichletDisc[i].disc)->clear_dirichlet_jacobian(J, u, m_vDirichletDisc[i].si) != IAssemble_OK)
					return IAssemble_ERROR;

			return IAssemble_OK;
		}

		IAssembleReturn assemble_defect(vector_type& d, const discrete_function_type& u)
		{
			if(!check_solution(u)) return IAssemble_ERROR;

			for(size_t i = 0; i < m_vElemDisc.size(); ++i)
				if(!AssembleDefect<IElemDisc<TAlgebra>,TDiscreteFunction,TAlgebra>
					(*m_vElemDisc[i].disc, d, u, m_vElemDisc[i].fcts, m_vElemDisc[i].si, m_vElemDisc[i].dim))
						return IAssemble_ERROR;

			for(size_t i = 0; i < m_vCoupledElemDisc.size(); ++i)
				if(!AssembleDefect<CoupledSystem<TDiscreteFunction, TAlgebra>,TDiscreteFunction,TAlgebra>
					(*m_vCoupledElemDisc[i].disc, d, u, m_vCoupledElemDisc[i].fcts, m_vCoupledElemDisc[i].si,  m_vCoupledElemDisc[i].dim))
						return IAssemble_ERROR;

			for(size_t i = 0; i < m_vDirichletDisc.size(); ++i)
				if((m_vDirichletDisc[i].disc)->clear_dirichlet_defect(d, u, m_vDirichletDisc[i].si) != IAssemble_OK)
					return IAssemble_ERROR;
			return IAssemble_OK;
		}

		IAssembleReturn assemble_linear(matrix_type& mat, vector_type& rhs, const discrete_function_type& u)
		{
			if(!check_solution(u)) return IAssemble_ERROR;

			for(size_t i = 0; i < m_vElemDisc.size(); ++i)
				if(!AssembleLinear<IElemDisc<TAlgebra>,TDiscreteFunction,TAlgebra>
					(*m_vElemDisc[i].disc, mat, rhs, u, m_vElemDisc[i].fcts, m_vElemDisc[i].si, m_vElemDisc[i].dim))
						return IAssemble_ERROR;

			for(size_t i = 0; i < m_vCoupledElemDisc.size(); ++i)
				if(!AssembleLinear<CoupledSystem<TDiscreteFunction, TAlgebra>,TDiscreteFunction,TAlgebra>
					(*m_vCoupledElemDisc[i].disc, mat, rhs, u, m_vCoupledElemDisc[i].fcts, m_vCoupledElemDisc[i].si, m_vCoupledElemDisc[i].dim))
						return IAssemble_ERROR;

			for(size_t i = 0; i < m_vConstraintsPostProcess.size(); ++i)
				if(m_vConstraintsPostProcess[i]->post_process_linear(mat, rhs, u) != IAssemble_OK)
					return IAssemble_ERROR;

			for(size_t i = 0; i < m_vDirichletDisc.size(); ++i)
				if((m_vDirichletDisc[i].disc)->set_dirichlet_linear(mat, rhs, u, m_vDirichletDisc[i].si) != IAssemble_OK)
					return IAssemble_ERROR;
			return IAssemble_OK;
		}

		IAssembleReturn assemble_solution(discrete_function_type& u)
		{
			for(size_t i = 0; i < m_vDirichletDisc.size(); ++i)
				if((m_vDirichletDisc[i].disc)->set_dirichlet_solution(u, m_vDirichletDisc[i].si) != IAssemble_OK)
					return IAssemble_ERROR;

			return IAssemble_OK;
		}

		///////////////////////
		// Time dependent part
		///////////////////////
		IAssembleReturn assemble_jacobian(matrix_type& J, const discrete_function_type& u, number time, number s_m, number s_a)
		{
			if(!check_solution(u)) return IAssemble_ERROR;

			for(size_t i = 0; i < m_vElemDisc.size(); ++i)
				if(!AssembleJacobian<IElemDisc<TAlgebra>,TDiscreteFunction,TAlgebra>
					(*m_vElemDisc[i].disc, J, u, m_vElemDisc[i].fcts, m_vElemDisc[i].si, m_vElemDisc[i].dim, time, s_m, s_a))
						return IAssemble_ERROR;

			for(size_t i = 0; i < m_vCoupledElemDisc.size(); ++i)
				if(!AssembleJacobian<CoupledSystem<TDiscreteFunction, TAlgebra>,TDiscreteFunction,TAlgebra>
					(*m_vCoupledElemDisc[i].disc, J, u, m_vCoupledElemDisc[i].fcts, m_vCoupledElemDisc[i].si, m_vCoupledElemDisc[i].dim, time, s_m, s_a))
						return IAssemble_ERROR;

			for(size_t i = 0; i < m_vDirichletDisc.size(); ++i)
				if((m_vDirichletDisc[i].disc)->clear_dirichlet_jacobian(J, u, m_vDirichletDisc[i].si, time) != IAssemble_OK)
					return IAssemble_ERROR;
			return IAssemble_OK;
		}

		IAssembleReturn assemble_defect(vector_type& d, const discrete_function_type& u, number time, number s_m, number s_a)
		{
			if(!check_solution(u)) return IAssemble_ERROR;

			for(size_t i = 0; i < m_vElemDisc.size(); ++i)
				if(!AssembleDefect<IElemDisc<TAlgebra>,TDiscreteFunction,TAlgebra>
					(*m_vElemDisc[i].disc, d, u, m_vElemDisc[i].fcts, m_vElemDisc[i].si, m_vElemDisc[i].dim, time, s_m, s_a))
						return IAssemble_ERROR;

			for(size_t i = 0; i < m_vCoupledElemDisc.size(); ++i)
				if(!AssembleDefect<CoupledSystem<TDiscreteFunction, TAlgebra>,TDiscreteFunction,TAlgebra>
					(*m_vCoupledElemDisc[i].disc, d, u, m_vCoupledElemDisc[i].fcts, m_vCoupledElemDisc[i].si,  m_vCoupledElemDisc[i].dim, time, s_m, s_a))
						return IAssemble_ERROR;

			for(size_t i = 0; i < m_vDirichletDisc.size(); ++i)
				if((m_vDirichletDisc[i].disc)->clear_dirichlet_defect(d, u, m_vDirichletDisc[i].si, time) != IAssemble_OK)
					return IAssemble_ERROR;
			return IAssemble_OK;
		}

		IAssembleReturn assemble_linear(matrix_type& mat, vector_type& rhs, const discrete_function_type& u, number time, number s_m, number s_a)
		{
			if(!check_solution(u)) return IAssemble_ERROR;

			for(size_t i = 0; i < m_vElemDisc.size(); ++i)
				if(!AssembleLinear<IElemDisc<TAlgebra>,TDiscreteFunction,TAlgebra>
					(*m_vElemDisc[i].disc, mat, rhs, u, m_vElemDisc[i].fcts, m_vElemDisc[i].si, m_vElemDisc[i].dim, time, s_m, s_a))
						return IAssemble_ERROR;

			for(size_t i = 0; i < m_vCoupledElemDisc.size(); ++i)
				if(!AssembleLinear<CoupledSystem<TDiscreteFunction, TAlgebra>,TDiscreteFunction,TAlgebra>
					(*m_vCoupledElemDisc[i].disc, mat, rhs, u, m_vCoupledElemDisc[i].fcts, m_vCoupledElemDisc[i].si, m_vCoupledElemDisc[i].dim, time, s_m, s_a))
						return IAssemble_ERROR;

			for(size_t i = 0; i < m_vDirichletDisc.size(); ++i)
				if((m_vDirichletDisc[i].disc)->set_dirichlet_linear(mat, rhs, u, m_vDirichletDisc[i].si, time) != IAssemble_OK)
					return IAssemble_ERROR;
			return IAssemble_OK;
		}

		IAssembleReturn assemble_solution(discrete_function_type& u, number time)
		{
			if(!check_solution(u)) return IAssemble_ERROR;

			for(size_t i = 0; i < m_vDirichletDisc.size(); ++i)
				if((m_vDirichletDisc[i].disc)->set_dirichlet_solution(u, m_vDirichletDisc[i].si, time) != IAssemble_OK)
					return IAssemble_ERROR;
			return IAssemble_OK;
		}

	public:
		bool add(IElemDisc<TAlgebra>& elemDisc, const FunctionGroup& fcts, int si, int dim)
		{
			// check if number of functions match
			if(elemDisc.num_fct() != fcts.num_fct())
			{
				UG_LOG("Wrong number of local functions given for Elemet Discretization.\n");
				UG_LOG("Needed: " << elemDisc.num_fct() << ", given: " << fcts.num_fct() << ".\n");
				UG_LOG("Cannot add element discretization.\n");
				return false;
			}
			m_vElemDisc.push_back(ElemDisc(si, dim, elemDisc, fcts));
			return true;
		}

	protected:
		struct ElemDisc
		{
			ElemDisc(int si_, int dim_, IElemDisc<TAlgebra>& disc_, const FunctionGroup& fcts_) :
				si(si_), dim(dim_), disc(&disc_), fcts(fcts_) {};

			int si;
			int dim;
			IElemDisc<TAlgebra>* disc;
			FunctionGroup fcts;
		};

		std::vector<ElemDisc> m_vElemDisc;

	public:
		bool add(CoupledSystem<TDiscreteFunction, TAlgebra>& coupledElemDisc, const FunctionGroup& fcts, int si, int dim)
		{
			// check if number of functions match
			if(coupledElemDisc.num_fct() != fcts.num_fct())
			{
				UG_LOG("Wrong number of local functions given for Elemet Discretization.\n");
				UG_LOG("Needed: " << coupledElemDisc.num_fct() << ", given: " << fcts.num_fct() << ".\n");
				UG_LOG("Cannot add element discretization.\n");
				return false;
			}
			m_vCoupledElemDisc.push_back(CoupledDisc(si, dim, coupledElemDisc, fcts));
			return true;
		}

	protected:
		struct CoupledDisc
		{
			CoupledDisc(int si_, int dim_, CoupledSystem<TDiscreteFunction, TAlgebra>& disc_, const FunctionGroup& fcts_) :
				si(si_), dim(dim_), disc(&disc_), fcts(fcts_) {};

			int si;
			int dim;
			CoupledSystem<TDiscreteFunction, TAlgebra>* disc;
			FunctionGroup fcts;
		};

		std::vector<CoupledDisc> m_vCoupledElemDisc;

	public:
		bool add(IConstraintsPostProcess<TDiscreteFunction, TAlgebra>& constraintsPP)
		{
			m_vConstraintsPostProcess.push_back(&constraintsPP);
			return true;
		}

	protected:
		std::vector<IConstraintsPostProcess<TDiscreteFunction, TAlgebra>*> m_vConstraintsPostProcess;


	protected:
		bool check_solution(const discrete_function_type& u)
		{
			for(size_t i = 0; i < m_vElemDisc.size(); ++i)
			{
				for(size_t fct = 0; fct < (m_vElemDisc[i].disc)->num_fct(); ++fct)
				{
					if((u.local_shape_function_set_id((m_vElemDisc[i].fcts)[fct]) !=	(m_vElemDisc[i].disc)->local_shape_function_set_id(fct)))
					{
						UG_LOG("fct " << fct << " of elem disc is fct " << (m_vElemDisc[i].fcts)[fct] << "of solution\n");
						UG_LOG("id " << u.local_shape_function_set_id((m_vElemDisc[i].fcts)[fct]) << ", " << (m_vElemDisc[i].disc)->local_shape_function_set_id(fct) << "\n");
						UG_LOG("Solution does not match requirements of discretization.\n");
						return false;
					}
				}
			}
			return true;
		}

	public:
		bool add_dirichlet_bnd(IDirichletBoundaryValues<discrete_function_type>& bnd_disc, int si)
		{
			m_vDirichletDisc.push_back(DirichletDisc(si, &bnd_disc));
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

		std::vector<DirichletDisc> m_vDirichletDisc;

	protected:
		// TODO: What is this function used for???? Do we have to include it
		virtual size_t num_fct() const
		{
			size_t sum = 0;
			for(size_t i = 0; i < m_vElemDisc.size(); ++i)
				sum += m_vElemDisc[i].fcts.num_fct();

			return sum;
		}

		virtual bool is_dirichlet(int si, size_t fct)
		{
			for(size_t i = 0; i < m_vDirichletDisc.size(); ++i)
			{
				if(m_vDirichletDisc[i].si != si) continue;
				if((m_vDirichletDisc[i].disc)->is_dirichlet(fct) == true) return true;
			}
			return false;
		}
};

} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__DOMAIN_DISCRETIZATION__ */
