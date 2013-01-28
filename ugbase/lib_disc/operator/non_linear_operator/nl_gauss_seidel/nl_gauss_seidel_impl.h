/*
 * nl_gauss_seidel_impl.h
 *
 *  Created on: 07.01.2013
 *  (main parts are based on the structure of
 *  	newton_impl.h by Andreas Vogel)
 *
 *  Author: raphaelprohl
 *
 *  Nonlinear GaussSeidel-method: (c.f. "Iterative Solution of nonlinear
 *  				equations in several variables" by Ortega/Rheinboldt)
 *
 * 	Let L(u) denote a nonlinear functional of n components (l_1,...,l_n).
 * 	Then the basic step of the nonlinear GaussSeidel method is to solve the
 * 	i-th equation
 *
 * 	l_i(u_1^{k+1},...,u_{i-1}^{k+1},u_i,u_{i+1}^{k},...,u_{n}^{k}) = 0
 *
 * 	for u_i and to set u_i^{k+1} = u_i. Here k denotes the iteration-index.
 * 	Note, that the already computed, updated values (.)^{k+1} are used in this
 * 	method.
 * 	Thus, in order to obtain u^{k+1} from u^k, we solve successively the n
 * 	dimensional nonlinear equations for i = 1,...,n. Here this is done
 * 	by a scalar newton step for every i. But every other scalar nonlinear method
 * 	could be applied as well.
 *
 * 	Using a damped version of the nonlinear GaussSeidel method (= nonlinear
 * 	SOR-method) results in the following update of the variables
 *
 * 	u_i^{k+1} = u_i^k + damp * (u_i -u_i^k).
 */

#ifndef NL_GAUSS_SEIDEL_IMPL_H_
#define NL_GAUSS_SEIDEL_IMPL_H_

// extern includes
#include <iostream>
//#include <list>

#include "lib_disc/function_spaces/grid_function_util.h"
#include "lib_disc/common/local_algebra.h"
#include "nl_gauss_seidel.h"

#define PROFILE_NL_GAUSSSEIDEL
#ifdef PROFILE_NL_GAUSSSEIDEL
	#define NL_GAUSSSEIDEL_PROFILE_FUNC()		PROFILE_FUNC_GROUP("NL GaussSeidel")
	#define NL_GAUSSSEIDEL_PROFILE_BEGIN(name)	PROFILE_BEGIN_GROUP(name, "NL GaussSeidel")
	#define NL_GAUSSSEIDEL_PROFILE_END()		PROFILE_END()
#else
	#define NL_GAUSSSEIDEL_PROFILE_FUNC()
	#define NL_GAUSSSEIDEL_PROFILE_BEGIN(name)
	#define NL_GAUSSSEIDEL_PROFILE_END()
#endif

namespace ug{

template <typename TDomain, typename TAlgebra>
NLGaussSeidelSolver<TDomain, TAlgebra>::
NLGaussSeidelSolver(SmartPtr<approx_space_type> spApproxSpace,
			SmartPtr<IConvergenceCheck<vector_type> > spConvCheck) :
			m_spApproxSpace(spApproxSpace),
			m_spConvCheck(spConvCheck),
			m_damp(1.0),
			m_dgbCall(0)
{};

template <typename TDomain, typename TAlgebra>
NLGaussSeidelSolver<TDomain, TAlgebra>::
NLGaussSeidelSolver() :
	m_spApproxSpace(NULL),
	m_spConvCheck(new StdConvCheck<vector_type>(10, 1e-8, 1e-10, true)),
	m_damp(1.0),
	m_dgbCall(0)
{};


template <typename TDomain, typename TAlgebra>
void NLGaussSeidelSolver<TDomain, TAlgebra>::
set_convergence_check(SmartPtr<IConvergenceCheck<vector_type> > spConvCheck)
{
	m_spConvCheck = spConvCheck;
	m_spConvCheck->set_offset(3);
	m_spConvCheck->set_symbol('#');
	m_spConvCheck->set_name("Nonlinear Gauss Seidel Solver");
}

template <typename TDomain, typename TAlgebra>
bool
NLGaussSeidelSolver<TDomain, TAlgebra>::
init(SmartPtr<IOperator<vector_type> > N)
{
	NL_GAUSSSEIDEL_PROFILE_BEGIN(NL_GAUSSSEIDELSolver_init);
	m_N = N.template cast_dynamic<AssembledOperator<TAlgebra> >();
	if(m_N.invalid())
		UG_THROW("NLGaussSeidelSolver: currently only works for AssembledDiscreteOperator.");

	m_pAss = m_N->get_assemble();

	m_gridLevel = m_N->level();

	//	Check for approxSpace
	if(m_spApproxSpace.invalid())
		UG_THROW("NLGaussSeidelSolver::prepare: Approximation Space not set.");

	//	set DoF distribution type
	if(m_gridLevel.type() == GridLevel::LEVEL)
		m_spLevDD = m_spApproxSpace->level_dof_distribution(m_gridLevel.level());
	else if (m_gridLevel.type() == GridLevel::SURFACE)
		m_spSurfDD = m_spApproxSpace->surface_dof_distribution(m_gridLevel.level());
	else
		UG_THROW("Grid Level not recognized.");

	return true;
}

template <typename TDomain, typename TAlgebra>
bool NLGaussSeidelSolver<TDomain, TAlgebra>::prepare(vector_type& u)
{
	return true;
}


/*template <typename TDomain, typename TAlgebra>
bool NLGaussSeidelSolver<TDomain, TDD, TAlgebra>::preprocess(gridfunc_type& u)
{
	//	fill m_DiagMarker by selecting all elements
	//	which have contributions to the diag

	domain_type& dom = *u.domain();
	typename domain_type::grid_type& grid = *dom.grid();

	static const int dim = gridfunc_type::dim;
	typedef typename gridfunc_type::template dim_traits<dim>::geometric_base_object geometric_base_object;
	typedef typename gridfunc_type::template dim_traits<dim>::const_iterator const_iterator;

	const_iterator iter = u.template begin<geometric_base_object>();
	const_iterator end = u.template end<geometric_base_object>();

	LocalIndices ind; LocalVector locU;

	//std::vector<BoolMarker>& vDiagMarker = m_vDiagMarker;

	m_vDiagMarker.resize(u.size());

	UG_LOG("u_size:" << u.size() << "\n");

	//	loop all DoFs
	for (size_t i = 0; i < u.size(); i++)
	{
		m_vDiagMarker[i].assign_grid(grid);
		m_vDiagMarker[i].clear(); // clear necessary? default = unmark?!
	}

	size_t count_elem = 0;

	//	loop over all elements on subset si
	for(;iter != end; ++iter)
	{
		//	get element
		geometric_base_object* elem = *iter;

		size_t count_i = 0;
		size_t count_noti = 0;
		for (size_t i = 0; i < u.size(); i++)
		{
			// unmark
			m_vDiagMarker[i].unmark(*iter);
		}

		if(m_vDiagMarker[1].is_marked(*iter))
		{
			UG_LOG("1-ter index ist markiert!\n");
		}

		for (size_t i = 0; i < u.size(); i++)
		{
			if(m_vDiagMarker[i].is_marked(*iter))
			{
				count_i++;
			}
			else
			{
				count_noti++;
			}

		}
		UG_LOG("Elem has " << count_i << " marks \n");
		UG_LOG(count_noti << " unmarks \n");

		count_i = 0;
		count_noti = 0;

		// 	get global indices
		u.indices(elem, ind);

		// 	adapt local algebra
		locU.resize(ind);

		//	local vector extract -> locU
		GetLocalVector(locU, u);

		for(size_t fct=0; fct < locU.num_all_fct(); ++fct)
			for(size_t dof=0; dof < locU.num_all_dof(fct); ++dof)
			{
				size_t globIndex = ind.index(fct,dof);
				//const size_t comp = ind.comp(fct,dof);
				UG_LOG("index:" << globIndex << "\n");

				if(m_vDiagMarker.at(1).is_marked(*iter))
				{
					UG_LOG("1-ter index ist markiert!\n");
				}
				//	mark elem in order to show that it
				//	has got an effect on globIndex
				//	(these elems wont be skipped in assembling!)
				//BoolMarker* p_vDiagMarker = m_vDiagMarker.data();

				BoolMarker& vDiagMarkerGlobInd = m_vDiagMarker.at(globIndex);

				vDiagMarkerGlobInd.mark(*iter);

				if(m_vDiagMarker.at(1).is_marked(*iter))
				{
					UG_LOG("1-ter index ist markiert!\n");
				}

				if(m_vDiagMarker.at(2).is_marked(*iter))
				{
					UG_LOG("2-ter index ist markiert!\n");
				}

			}

		for (size_t i = 0; i < u.size(); i++)
		{
			if(m_vDiagMarker.at(i).is_marked(*iter))
			{
				count_i++;
				UG_LOG("elem has influence on: " << i << "\n");
			}
			else
			{
				count_noti++;
			}

		}
		UG_LOG("Elem has " << count_i << " influences \n");
		UG_LOG(count_noti << " elems should be skipped \n");
		count_elem++;
		UG_LOG("\n");

	}

	UG_LOG("Loop over " << count_elem << "elems in preprocess \n");

	return true;
}*/


template <typename TDomain, typename TAlgebra>
bool NLGaussSeidelSolver<TDomain, TAlgebra>::apply(vector_type& u)
{
	NL_GAUSSSEIDEL_PROFILE_BEGIN(NL_GAUSSSEIDELSolver_apply);
	//	increase call count
	m_dgbCall++;

	//	Check for approxSpace
	if(m_spApproxSpace.invalid())
		UG_THROW("NLGaussSeidelSolver::apply: Approximation Space not set.");

	//	Check for DoF distribution
	if(m_spLevDD.invalid() && m_spSurfDD.invalid())
		UG_THROW("NLGaussSeidelSolver::apply: DoFDistribution not set."
				" 'NLGaussSeidelSolver::init'-call is necessary!");

	//	Jacobian
	if(m_J.invalid() || m_J->discretization() != m_pAss) {
		m_J = CreateSmartPtr(new AssembledLinearOperator<TAlgebra>(*m_pAss));
		m_J->set_level(m_gridLevel);
	}

	//	get #indices of gridFunction u
	const size_t n_indices = u.size();

	//	resize
	try{
		m_d.resize(n_indices); m_d = u;
		m_c.resize(n_indices); m_c = u;
	}UG_CATCH_THROW("NLGaussSeidelSolver::apply: Resize of Defect/Correction failed.");

	//	Set dirichlet values
	try{
		m_N->prepare(m_d, u);
	}
	UG_CATCH_THROW("NLGaussSeidelSolver::apply: Prepare of Operator failed.");

	// 	Compute first Defect d = L(u)
	try{
		NL_GAUSSSEIDEL_PROFILE_BEGIN(NL_GAUSSSEIDELComputeDefect1);
		m_N->apply(m_d, u);
		NL_GAUSSSEIDEL_PROFILE_END();
	}UG_CATCH_THROW("NLGaussSeidelSolver::apply: "
			"Computation of Start-Defect failed.");

	//	write start defect for debug
	int loopCnt = 0;
	char ext[20]; sprintf(ext, "_iter%03d", loopCnt);
	std::string name("NLGaussSeidel_Defect");
	name.append(ext);
	write_debug(m_d, name.c_str());
	write_debug(u, "NLGaussSeidel_StartSolution");

	// 	start convergence check
	m_spConvCheck->start(m_d);

	matrix_type& J = m_J->get_matrix();
	number damp = m_damp;

	// some vars for BoolMarker
	TDomain& dom = *m_spApproxSpace->domain();
	typename TDomain::grid_type& grid = *dom.grid();

	typedef typename domain_traits<TDomain::dim>::geometric_base_object geometric_base_object;
	typedef typename TDomain::grid_type::template traits<geometric_base_object>::iterator
			ElemIter;

	ElemIter iter;
	ElemIter iterBegin = grid.template begin<geometric_base_object>();
	ElemIter iterEnd = grid.template end<geometric_base_object>();

	LocalIndices ind; LocalVector locU;

	m_vDiagMarker.assign_grid(grid);

	//	creating an ElemList
	//	TODO: use ElemList as a member-variable of NL_GaussSeidel-class?
	std::list<geometric_base_object> ElemList;

	//	loop iteration
	while(!m_spConvCheck->iteration_ended())
	{
		// 	set correction c = 0
		NL_GAUSSSEIDEL_PROFILE_BEGIN(NL_GAUSSSEIDELSetCorretionZero);
		if(!m_c.set(0.0))
		{
			UG_LOG("ERROR in 'NLGaussSeidelSolver::apply':"
					" Cannot reset correction to zero.\n");
			return false;
		}
		NL_GAUSSSEIDEL_PROFILE_END();

		//	loop all DoFs
		for (size_t i = 0; i < n_indices; i++)
		{
			// 	Compute Jacobian J(u) using the updated u-components

			//	we only need J(i,i) and d(i) here!
			//	TODO: for DoF i create an ElemList!
			//	This ElemList should be passed to the assemble_funcs.
			//	iterBegin and iterEnd must be passed to elem_disc_assemble_util.h!

			// 	some debug-vars
			size_t count_i = 0;
			size_t count_noti = 0;

			//	loop over all elements on grid
			for(iter = iterBegin; iter != iterEnd; ++iter)
			{
				//	get element
				geometric_base_object* elem = *iter;

				if(m_gridLevel.type() == GridLevel::LEVEL)
					m_spLevDD->indices(elem, ind);
				else if (m_gridLevel.type() == GridLevel::SURFACE)
					m_spSurfDD->indices(elem, ind);
				else
					UG_THROW("Grid Level not recognized.");

				// 	adapt local algebra
				locU.resize(ind);

				//	local vector extract -> locU
				GetLocalVector(locU, u);

				bool skip_loop = false;
				m_vDiagMarker.unmark(*iter);

				for(size_t fct=0; fct < locU.num_all_fct() && !skip_loop; ++fct)
					for(size_t dof=0; dof < locU.num_all_dof(fct); ++dof)
					{
						size_t globIndex = ind.index(fct,dof);
						//UG_LOG("index:" << globIndex << "\n");

						//	mark elem in order to show that it
						//	has got an effect on globIndex
						//	(these elems won't be skipped in assembling!)

						if (globIndex == i)
						{
							m_vDiagMarker.mark(*iter);
							count_i++;
							//TODO: push_front oder push_back?
							//ElemList.push_front(*iter);
							skip_loop = true;
							break;
						}
						else{count_noti++;}
					}
			} //end(elem)

			//TODO: ElemList needs to be sorted by elem-types (edge,triangle,quad,...)

			//UG_LOG("DoF " << i << " is influenced by " << count_i << " elems \n");
			//UG_LOG(count_noti << " elems should be skipped \n");
			m_pAss->set_selector(&m_vDiagMarker);

			try{
				NL_GAUSSSEIDEL_PROFILE_BEGIN(NL_GAUSSSEIDELComputeJacobian);
				m_J->init(u); //TODO: pass ElemList here!
				NL_GAUSSSEIDEL_PROFILE_END();
			}UG_CATCH_THROW("NLGaussSeidelSolver::apply: "
					"Initialization of Jacobian failed.");

			//	Write Jacobian for debug
			std::string matname("NLGaussSeidel_Jacobian");
			matname.append(ext);
			write_debug(m_J->get_matrix(), matname.c_str());

			//	get i,i-th block of J: J(i,i)
			//	depending on the AlgebraType J(i,i) is a 1x1, 2x2, 3x3 Matrix
			//	m_c_i = m_damp * d_i /J_ii
			InverseMatMult(m_c[i], damp, J(i,i), m_d[i]);

			// 	update i-th block of solution
			u[i] -= m_c[i];

			m_pAss->set_selector(NULL);
			// 	Compute d = L(u) using the updated u-blocks
			//	TODO: replace prepare(m_d,u) & apply(m_d,u)!
			//	We only need the new m_d due
			//  to the updated block u_i! (not due to u!)
			NL_GAUSSSEIDEL_PROFILE_BEGIN(NL_GAUSSSEIDELComputeLastCompDefect);
			m_N->prepare(m_d, u);
			m_N->apply(m_d, u);
			NL_GAUSSSEIDEL_PROFILE_END();
		}

		// 	check convergence
		m_spConvCheck->update(m_d);
	}

	return m_spConvCheck->post();
}

template <typename TDomain, typename TAlgebra>
void NLGaussSeidelSolver<TDomain, TAlgebra>::write_debug(const vector_type& vec, const char* filename)
{
//	add iter count to name
	std::string name(filename);
	char ext[20]; sprintf(ext, "_call%03d", m_dgbCall);
	name.append(ext).append(".vec");

//	write
	base_writer_type::write_debug(vec, name.c_str());
}

template <typename TDomain, typename TAlgebra>
void NLGaussSeidelSolver<TDomain, TAlgebra>::write_debug(const matrix_type& mat, const char* filename)
{
//	add iter count to name
	std::string name(filename);
	char ext[20]; sprintf(ext, "_call%03d", m_dgbCall);
	name.append(ext).append(".mat");

//	write
	base_writer_type::write_debug(mat, name.c_str());
}

}

#endif /* NL_GAUSS_SEIDEL_IMPL_H_ */
