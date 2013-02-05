/*
 * nl_gauss_seidel_impl.h
 *
 *  Created on: 07.01.2013
 *  (main parts are based on the structure of
 *  	newton_impl.h and some ideas of Sebastian Reiter & Andreas Vogel)
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
	if(m_pAss == NULL)
		UG_THROW("AssembledLinearOperator: Assembling routine not set.");

	//	Check for approxSpace
	if(m_spApproxSpace.invalid())
		UG_THROW("NLGaussSeidelSolver::prepare: Approximation Space not set.");

	m_gridLevel = m_N->level();

	//	set DoF distribution type
	if(m_gridLevel.type() == GridLevel::LEVEL)
	{
		m_spLevDD = m_spApproxSpace->level_dof_distribution(m_gridLevel.level());
		//	note: #ftcs needs to be constant over all subsets!!!
		//m_num_fct = m_spLevDD
	}
	else if (m_gridLevel.type() == GridLevel::SURFACE)
		m_spSurfDD = m_spApproxSpace->surface_dof_distribution(m_gridLevel.level());
	else
		UG_THROW("Grid Level not recognized.");

	/*TDomain& dom = *m_spApproxSpace->domain();
	typename TDomain::grid_type& grid = *dom.grid();

	grid.attach_to_vertices(m_aElemList);
	m_aaElemList.access(grid, m_aElemList);*/

	return true;
}


template <typename TDomain, typename TAlgebra>
bool NLGaussSeidelSolver<TDomain, TAlgebra>::prepare(vector_type& u)
{
	//	In this method a elemList is created for every DoF i. These element-list
	//	incorporate by selecting all elements
	//	which have contributions to the diag

	//	Check for approxSpace
	if(m_spApproxSpace.invalid())
		UG_THROW("NLGaussSeidelSolver::apply: Approximation Space not set.");

	//	Check for DoF distribution
	if(m_spLevDD.invalid() && m_spSurfDD.invalid())
		UG_THROW("NLGaussSeidelSolver::apply: DoFDistribution not set."
				" 'NLGaussSeidelSolver::init'-call is necessary!");

	// some vars for BoolMarker
	TDomain& dom = *m_spApproxSpace->domain();
	typename TDomain::grid_type& grid = *dom.grid();

	typedef typename domain_traits<TDomain::dim>::geometric_base_object geometric_base_object;
	typedef typename TDomain::grid_type::template traits<geometric_base_object>::iterator
			ElemIter;

	ElemIter iterBegin = grid.template begin<geometric_base_object>();
	ElemIter iterEnd = grid.template end<geometric_base_object>();

	typedef typename elemList::iterator ListIter;

	/* VERSION †BER ATTACHMENT:
	typename Grid::traits<geometric_base_object>::secure_container elems;

	int vtr_count = 0;

	//	loop over all vertices on grid/multigrid
	for(VertexBaseIterator iter = grid.vertices_begin(); iter != grid.vertices_end(); ++iter)
	{
		grid.associated_elements(elems, *iter);
		for(size_t i = 0; i < elems.size(); ++i) m_aaElemList[*iter].push_back(elems[i]);
	}*/

	LocalVector locU; LocalIndices ind;

	m_vElemList.resize(u.size());

	//	loop over all elements on grid/multigrid
	for(ElemIter elemIter = iterBegin; elemIter != iterEnd; ++elemIter)
	{
		//	get element
		geometric_base_object* elem = *elemIter;

		if(m_gridLevel.type() == GridLevel::LEVEL)
			m_spLevDD->indices(elem, ind);
		else if (m_gridLevel.type() == GridLevel::SURFACE)
			m_spSurfDD->indices(elem, ind);
		else
			UG_THROW("Grid Level not recognized.");

		// 	adapt local algebra
		locU.resize(ind);

		//TODO: bool skip_loop = false; ?
		//for(size_t fct=0; fct < locU.num_all_fct() && !skip_loop; ++fct)

		bool elem_in_list = false;

		for(size_t fct=0; fct < locU.num_all_fct(); ++fct)
			for(size_t dof=0; dof < locU.num_all_dof(fct); ++dof)
			{
				size_t globIndex = ind.index(fct,dof);
				size_t globComp = ind.comp(fct,dof);
				//UG_LOG("index:" << globIndex << "\n");
				//UG_LOG("comp:" << globComp << "\n");

				//	check if the globIndex-th elemList incorporates the element *iter
				for(ListIter listIter = m_vElemList[globIndex].begin();
						listIter != m_vElemList[globIndex].end(); ++listIter)
					if (*listIter == *elemIter)
						elem_in_list = true;

				//	only add the element *elemIter to the globIndex-th elemList,
				//	if the list does not incorporate the element already
				if (!elem_in_list) m_vElemList[globIndex].push_back(*elemIter);

			} //end (dof)
	} //end (elem)

	/*//if (m_vElemList[i].size() > 8 )
	for (size_t i = 0; i < u.size(); i++)
		UG_LOG(" FŸr DoF " << i << " haben wir " << m_vElemList[i].size() << " Elemente in der Liste \n ");*/

	return true;
}


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

	//	resize
	try{
		m_d.resize(u.size()); m_d.copy_layouts(u);
		m_c_comp.resize(1); m_c_comp.copy_layouts(u);
		m_u_comp.resize(1); m_u_comp.copy_layouts(u);
		m_d_comp.resize(1); m_d_comp.copy_layouts(u);
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

	matrix_type J_block;

	// some vars for BoolMarker
	TDomain& dom = *m_spApproxSpace->domain();
	typename TDomain::grid_type& grid = *dom.grid();

	typedef typename domain_traits<TDomain::dim>::geometric_base_object geometric_base_object;
	typedef typename TDomain::grid_type::template traits<geometric_base_object>::iterator
			ElemIter;

	ElemIter iter;
	ElemIter iterBegin = grid.template begin<geometric_base_object>();
	ElemIter iterEnd = grid.template end<geometric_base_object>();

	m_sel.assign_grid(grid);

	//	loop iteration
	while(!m_spConvCheck->iteration_ended())
	{
		// 	set correction c = 0
		NL_GAUSSSEIDEL_PROFILE_BEGIN(NL_GAUSSSEIDELSetCorretionZero);
		if(!m_c_comp.set(0.0))
		{
			UG_LOG("ERROR in 'NLGaussSeidelSolver::apply':"
					" Cannot reset correction to zero.\n");
			return false;
		}
		NL_GAUSSSEIDEL_PROFILE_END();

		//	loop all DoFs
		for (size_t i = 0; i < u.size(); i++)
		{
			// 	Compute Jacobian J(u) using the updated u-components

			//	we only need J(i,i) and d(i) here!
			//	use the i-th ElemList to select those elements
			//	which are associated to (in the neighborhood of) DoF i!

			m_sel.clear();
			m_sel.select(m_vElemList[i].begin(), m_vElemList[i].end());

			//int nSelElem = m_sel.num<geometric_base_object>();
			//UG_LOG("for DoF: " << i << " " << nSelElem << " elems are selected\n");

			m_u_comp[0] = u[i];
			m_d_comp[0] = m_d[i];
			//m_J_comp[0][0] = J(i,i);

			m_pAss->set_selector(&m_sel);
			//m_pAss->ass_index(i);

			try{
				NL_GAUSSSEIDEL_PROFILE_BEGIN(NL_GAUSSSEIDELComputeJacobian);
				//m_J->init(m_u_comp);
				m_J->init(u);
				//m_pAss->assemble_jacobian(J, u, m_gridLevel.level());
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
			//m_J_comp = J(i,i);
			InverseMatMult(m_c_comp[0], damp, J(i,i) , m_d_comp[0]);

			//UG_LOG("J(" << i << "," << i << "): " << J(i,i) << "\n");

			// 	update i-th block of solution
			u[i] -= m_c_comp[0];

			m_pAss->set_selector(NULL);
			// 	Compute d = L(u) using the updated u-blocks
			//	TODO: replace prepare(m_d,u) & apply(m_d,u)!
			//	We only need the new m_d due
			//  to the updated block u_i! (not due to u!)
			//m_d[i] = m_d_comp;


			//if (i == (n_indices - 1))
			//{
				NL_GAUSSSEIDEL_PROFILE_BEGIN(NL_GAUSSSEIDELComputeLastCompDefect);
				m_N->prepare(m_d, u);
				m_N->apply(m_d, u);
				NL_GAUSSSEIDEL_PROFILE_END();
				//m_d[i] = d_comp;
			//}
		}

		/*NL_GAUSSSEIDEL_PROFILE_BEGIN(NL_GAUSSSEIDELComputeLastCompDefect);
		m_N->prepare(m_d, u);
		m_N->apply(m_d, u);
		NL_GAUSSSEIDEL_PROFILE_END();*/

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
