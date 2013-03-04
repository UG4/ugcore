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


template <typename TAlgebra>
void LocalToGlobalMapper_NL_GS<TAlgebra>::AddLocalVec(vector_type& vec, const LocalVector& lvec, ConstSmartPtr<DoFDistribution> dd)
{
	const LocalIndices& ind = lvec.get_indices();

	for(size_t fct=0; fct < lvec.num_all_fct(); ++fct)
		for(size_t dof=0; dof < lvec.num_all_dof(fct); ++dof)
		{
			const size_t index = ind.index(fct,dof);
			if (index == m_assIndex)
			{
				const size_t comp = ind.comp(fct,dof);
				BlockRef(vec[0], comp) += lvec.value(fct,dof);
			}

		}
}

template <typename TAlgebra>
void LocalToGlobalMapper_NL_GS<TAlgebra>::AddLocalMatToGlobal(matrix_type& mat, const LocalMatrix& lmat, ConstSmartPtr<DoFDistribution> dd)
{
	const LocalIndices& rowInd = lmat.get_row_indices();
	const LocalIndices& colInd = lmat.get_col_indices();

	for(size_t fct1=0; fct1 < lmat.num_all_row_fct(); ++fct1)
		for(size_t dof1=0; dof1 < lmat.num_all_row_dof(fct1); ++dof1)
		{
			const size_t rowIndex = rowInd.index(fct1,dof1);
			if (rowIndex != m_assIndex)
				continue;

			const size_t rowComp = rowInd.comp(fct1,dof1);
			for(size_t fct2=0; fct2 < lmat.num_all_col_fct(); ++fct2)
				for(size_t dof2=0; dof2 < lmat.num_all_col_dof(fct2); ++dof2)
				{
					const size_t colIndex = colInd.index(fct2,dof2);
					if (colIndex == m_assIndex)
					{
						const size_t colComp = colInd.comp(fct2,dof2);
						BlockRef(mat(0, 0), rowComp, colComp)
									+= lmat.value(fct1,dof1,fct2,dof2);
					}
				}
		}
}

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

	m_spAss = m_N->discretization();
	if(m_spAss.invalid())
		UG_THROW("AssembledLinearOperator: Assembling routine not set.");

	//	Check for approxSpace
	if(m_spApproxSpace.invalid())
		UG_THROW("NLGaussSeidelSolver::prepare: Approximation Space not set.");

	m_gridLevel = m_N->level();

	//	set DoF distribution type
	if(m_gridLevel.type() == GridLevel::LEVEL)
		m_spLevDD = m_spApproxSpace->level_dof_distribution(m_gridLevel.level());
	else if (m_gridLevel.type() == GridLevel::SURFACE)
		m_spSurfDD = m_spApproxSpace->surface_dof_distribution(m_gridLevel.level());
	else
		UG_THROW("Grid Level not recognized.");

	//	set specific LocalToGlobalMapping
	//m_spMapping =

	/*TDomain& dom = *m_spApproxSpace->domain();
	typename TDomain::grid_type& grid = *dom.grid();

	grid.attach_to_vertices(m_aElemList);
	m_aaElemList.access(grid, m_aElemList);*/

	return true;
}


template <typename TDomain, typename TAlgebra>
bool NLGaussSeidelSolver<TDomain, TAlgebra>::prepare(vector_type& u)
{
	//	In this method a elemList is created for every DoF 'globIndex'.
	//	This element-list incorporates all elements which have contributions
	//	to 'globIndex'.

	//	Check for approxSpace
	if(m_spApproxSpace.invalid())
		UG_THROW("NLGaussSeidelSolver::apply: Approximation Space not set.");

	//	Check for DoF distribution
	if(m_spLevDD.invalid() && m_spSurfDD.invalid())
		UG_THROW("NLGaussSeidelSolver::apply: DoFDistribution not set."
				" 'NLGaussSeidelSolver::init'-call is necessary!");

	// some vars for Element iterators
	TDomain& dom = *m_spApproxSpace->domain();
	typename TDomain::grid_type& grid = *dom.grid();

	typedef typename domain_traits<TDomain::dim>::geometric_base_object geometric_base_object;
	typedef typename TDomain::grid_type::template traits<geometric_base_object>::iterator
			ElemIter;

	ElemIter iterBegin = grid.template begin<geometric_base_object>(grid.top_level());
	ElemIter iterEnd = grid.template end<geometric_base_object>(grid.top_level());

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

	//	loop over all elements in toplevel of the grid/multigrid
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

		bool elem_in_list = false;

		for(size_t fct=0; fct < locU.num_all_fct(); ++fct)
			for(size_t dof=0; dof < locU.num_all_dof(fct); ++dof)
			{
				size_t globIndex = ind.index(fct,dof);

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

	return true;
}


template <typename TDomain, typename TAlgebra>
bool NLGaussSeidelSolver<TDomain, TAlgebra>::apply(vector_type& u)
{
	NL_GAUSSSEIDEL_PROFILE_BEGIN(NL_GAUSSSEIDELSolver_apply);
	//	increase call count
	m_dgbCall++;

	//	Check for approxSpac
	if(m_spApproxSpace.invalid())
		UG_THROW("NLGaussSeidelSolver::apply: Approximation Space not set.");

	//	Check for DoF distribution
	if(m_spLevDD.invalid() && m_spSurfDD.invalid())
		UG_THROW("NLGaussSeidelSolver::apply: DoFDistribution not set."
				" 'NLGaussSeidelSolver::init'-call is necessary!");

	//	Jacobian
	if(m_J_block.invalid() || m_J_block->discretization() != m_spAss) {
		m_J_block = CreateSmartPtr(new AssembledLinearOperator<TAlgebra>(m_spAss));
		m_J_block->set_level(m_gridLevel);
	}

	//	resize
	try{
		m_d.resize(u.size()); m_c_block.resize(1); m_d_block.resize(1);
		#ifdef UG_PARALLEL
			m_d.copy_layouts(u); m_c_block.copy_layouts(u); m_d_block.copy_layouts(u);
		#endif
	}UG_CATCH_THROW("NLGaussSeidelSolver::apply: Resize of Defect/Correction failed.");

	//	Set dirichlet values
	try{
		m_N->prepare(u);
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

	// 	assign selector to grid
	TDomain& dom = *m_spApproxSpace->domain();
	typename TDomain::grid_type& grid = *dom.grid();
	m_sel.assign_grid(grid);

	matrix_type& J_block = m_J_block->get_matrix();

	AssAdapter<TAlgebra>& assAdapt = m_spAss->get_ass_adapter();

	//	loop iteration
	while(!m_spConvCheck->iteration_ended())
	{
		//bool active_set_changed = false;
		assAdapt.set_mapping(&m_map);

		//	loop all DoFs
		for (size_t i = 0; i < u.size(); i++)
		{
			// 	Compute Jacobian J(u) using the updated u-components

			//	since we only need J(i,i) and d(i) in every DoF loop,
			//	we access the i-th Element list indicating the neighboorhood of DoF i!
			m_sel.clear();
			m_sel.select(m_vElemList[i].begin(), m_vElemList[i].end());

			//	by passing the selector to the assembling the assemble operators
			//	are build up only by looping over the elements which has been selected
			assAdapt.set_selector(&m_sel);

			//	assemble only with respect to DoF i (causes resizing of matrices/vectors)
			assAdapt.set_ass_index(i);
			m_map.set_ass_index(i);

			try{
				NL_GAUSSSEIDEL_PROFILE_BEGIN(NL_GAUSSSEIDELComputeJacobian);
				m_J_block->init(u);
				NL_GAUSSSEIDEL_PROFILE_END();
			}UG_CATCH_THROW("NLGaussSeidelSolver::apply: "
					"Initialization of Jacobian failed.");

			//	Write Jacobian for debug
			/*std::string matname("NLGaussSeidel_Jacobian");
			sprintf(ext, "_iter%03d_DoF%03lu", loopCnt, (unsigned long) i);
			matname.append(ext);
			write_debug(m_J_block->get_matrix(), matname.c_str());*/

			//	get i-th block of defect d: d(i) =: m_d_block
			NL_GAUSSSEIDEL_PROFILE_BEGIN(NL_GAUSSSEIDELComputeLastCompDefect);
			m_N->apply(m_d_block, u);
			NL_GAUSSSEIDEL_PROFILE_END();

			//	get i,i-th block of J: J(i,i)
			//	depending on the AlgebraType J(i,i) is a 1x1, 2x2, 3x3 Matrix
			//	m_c_i = m_damp * d_i /J_ii
			NL_GAUSSSEIDEL_PROFILE_BEGIN(NL_GAUSSSEIDELInvertJ);
				InverseMatMult(m_c_block[0], m_damp, J_block(0,0) , m_d_block[0]);
			NL_GAUSSSEIDEL_PROFILE_END();

			/*if (m_bObs)
			{
				std::vector<size_t> active_set;
				vector_type diff;
				diff.resize(0);
				diff.resize(1);
				diff[0] = u[i] - m_obsVec[i];

				bool penetrate = false;
				//	loop fcts
				if (diff[0] > 0)
				{
					bool DoF_is_in_activeSet = false;
					penetrate = true;

					std::vector<size_t>::iterator iter;

					//	adds DoF i to active_set
					for (iter = active_set.begin(); iter != active_set.end(); ++iter)
						if (*iter == i)
							DoF_is_in_activeSet = true;

					if (!DoF_is_in_activeSet)
					{
						active_set.push_back(i);
						active_set_changed = true;
					}
				}

				if (penetrate)
					m_c_block[0] = 0.0;
			}*/

			// 	update i-th block of solution
			u[i] -= m_c_block[0];

		}

		//	if active set not changed, iteration converged!
	//	if (!active_set_changed)
		//	break;

		//	set mapping, selector and ass_index to NULL
		assAdapt.set_mapping();
		assAdapt.set_selector();
		assAdapt.set_ass_index();

		NL_GAUSSSEIDEL_PROFILE_BEGIN(NL_GAUSSSEIDELComputeLastCompDefect);
		m_N->prepare(u);
		m_N->apply(m_d, u);
		NL_GAUSSSEIDEL_PROFILE_END();

		//	update counter
		loopCnt++;
		sprintf(ext, "_iter%03d", loopCnt);

		// 	check convergence
		m_spConvCheck->update(m_d);

		//	write defect for debug
		std::string name("NLGaussSeidel_Defect"); name.append(ext);
		write_debug(m_d, name.c_str());
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
