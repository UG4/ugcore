/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Raphael Prohl
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

/*
 *  (main parts are based on the structure of
 *  	newton_impl.h and some ideas of Sebastian Reiter & Andreas Vogel)
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__NL_GAUSS_SEIDEL__NL_GAUSS_SEIDEL_IMPL_H_
#define __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__NL_GAUSS_SEIDEL__NL_GAUSS_SEIDEL_IMPL_H_

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
void LocalToGlobalMapperNLGS<TAlgebra>::add_local_vec_to_global(vector_type& vec,
		const LocalVector& lvec, ConstSmartPtr<DoFDistribution> dd)
{
	const LocalIndices& ind = lvec.get_indices();

	for(size_t fct=0; fct < lvec.num_all_fct(); ++fct)
		for(size_t dof=0; dof < lvec.num_all_dof(fct); ++dof)
		{
			const size_t index = ind.index(fct,dof);
			if (index == m_assemblingIndex)
			{
				const size_t comp = ind.comp(fct,dof);
				BlockRef(vec[0], comp) += lvec.value(fct,dof);
			}
		}

}

template <typename TAlgebra>
void LocalToGlobalMapperNLGS<TAlgebra>::add_local_mat_to_global(matrix_type& mat,
		const LocalMatrix& lmat, ConstSmartPtr<DoFDistribution> dd)
{
	const LocalIndices& rowInd = lmat.get_row_indices();
	const LocalIndices& colInd = lmat.get_col_indices();

	for(size_t fct1=0; fct1 < lmat.num_all_row_fct(); ++fct1)
		for(size_t dof1=0; dof1 < lmat.num_all_row_dof(fct1); ++dof1)
		{
			const size_t rowIndex = rowInd.index(fct1,dof1);

			if (rowIndex == m_assemblingIndex)
			{
				const size_t rowComp = rowInd.comp(fct1,dof1);
				for(size_t fct2=0; fct2 < lmat.num_all_col_fct(); ++fct2)
					for(size_t dof2=0; dof2 < lmat.num_all_col_dof(fct2); ++dof2)
					{
						const size_t colIndex = colInd.index(fct2,dof2);
						if (colIndex == m_assemblingIndex)
						{
							const size_t colComp = colInd.comp(fct2,dof2);
							BlockRef(mat(0, 0), rowComp, colComp)
									+= lmat.value(fct1,dof1,fct2,dof2);
						}
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
			m_bProjectedGS(false),
			m_dgbCall(0)
{};

template <typename TDomain, typename TAlgebra>
NLGaussSeidelSolver<TDomain, TAlgebra>::
NLGaussSeidelSolver() :
	m_spApproxSpace(NULL),
	m_spConvCheck(new StdConvCheck<vector_type>(10, 1e-8, 1e-10, true)),
	m_damp(1.0),
	m_bProjectedGS(false),
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
init(SmartPtr<IOperator<vector_type> > op)
{
	NL_GAUSSSEIDEL_PROFILE_BEGIN(NL_GAUSSSEIDELSolver_init);
	m_spAssOp = op.template cast_dynamic<AssembledOperator<TAlgebra> >();
	if(m_spAssOp.invalid())
		UG_THROW("NLGaussSeidelSolver: currently only works for AssembledDiscreteOperator.");

	m_spAss = m_spAssOp->discretization();
	if(m_spAss.invalid())
		UG_THROW("AssembledLinearOperator: Assembling routine not set.");

	//	Check for approxSpace
	if(m_spApproxSpace.invalid())
		UG_THROW("NLGaussSeidelSolver::prepare: Approximation Space not set.");

	m_gridLevel = m_spAssOp->level();

	//	set DoF distribution type
	if(m_gridLevel.type() == GridLevel::LEVEL)
		m_spLevDD = m_spApproxSpace->dof_distribution(m_gridLevel);
	else if (m_gridLevel.type() == GridLevel::SURFACE)
		m_spSurfDD = m_spApproxSpace->dof_distribution(m_gridLevel);
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

	typedef typename domain_traits<TDomain::dim>::grid_base_object grid_base_object;
	typedef typename TDomain::grid_type::template traits<grid_base_object>::iterator
			ElemIter;

	ElemIter iterBegin = grid.template begin<grid_base_object>(grid.top_level());
	ElemIter iterEnd = grid.template end<grid_base_object>(grid.top_level());

	typedef typename elemList::iterator ListIter;

	/* VERSION �BER ATTACHMENT:
	typename Grid::traits<grid_base_object>::secure_container elems;

	int vtr_count = 0;

	//	loop over all vertices on grid/multigrid
	for(VertexIterator iter = grid.vertices_begin(); iter != grid.vertices_end(); ++iter)
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
		grid_base_object* elem = *elemIter;

		if(m_gridLevel.type() == GridLevel::LEVEL)
			m_spLevDD->indices(elem, ind);
		else if (m_gridLevel.type() == GridLevel::SURFACE)
			m_spSurfDD->indices(elem, ind);
		else
			UG_THROW("Grid Level not recognized.");

		// 	adapt local algebra
		locU.resize(ind);

		bool elemInList = false;

		for(size_t fct=0; fct < locU.num_all_fct(); ++fct)
			for(size_t dof=0; dof < locU.num_all_dof(fct); ++dof)
			{
				size_t globIndex = ind.index(fct,dof);

				//	check if the globIndex-th elemList incorporates the element *iter
				for(ListIter listIter = m_vElemList[globIndex].begin();
						listIter != m_vElemList[globIndex].end(); ++listIter)
					if (*listIter == *elemIter)
						elemInList = true;

				//	only add the element *elemIter to the globIndex-th elemList,
				//	if the list does not incorporate the element already
				if (!elemInList) m_vElemList[globIndex].push_back(*elemIter);

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

	//	Check for approxSpace
	if(m_spApproxSpace.invalid())
		UG_THROW("NLGaussSeidelSolver::apply: Approximation Space not set.");

	//	Check for DoF distribution
	if(m_spLevDD.invalid() && m_spSurfDD.invalid())
		UG_THROW("NLGaussSeidelSolver::apply: DoFDistribution not set."
				" 'NLGaussSeidelSolver::init'-call is necessary!");

	//	Jacobian
	if(m_spJBlock.invalid() || m_spJBlock->discretization() != m_spAss) {
		m_spJBlock = make_sp(new AssembledLinearOperator<TAlgebra>(m_spAss));
		m_spJBlock->set_level(m_gridLevel);
	}

	//	create tmp vectors
	SmartPtr<vector_type> spD = u.clone_without_values();
	SmartPtr<vector_type> spC = u.clone_without_values();
	SmartPtr<vector_type> spDBlock = u.clone_without_values();

	//	resize vectors, because they are only used as block-vectors with one entry
	(*spC).resize(1); (*spDBlock).resize(1);

	//	Set dirichlet values
	try{
		m_spAssOp->prepare(u);
	}
	UG_CATCH_THROW("NLGaussSeidelSolver::apply: Prepare of Operator failed.");

	// 	Compute first Defect d = L(u)
	try{
		NL_GAUSSSEIDEL_PROFILE_BEGIN(NL_GAUSSSEIDELComputeDefect1);
		m_spAssOp->apply(*spD, u);
		NL_GAUSSSEIDEL_PROFILE_END();
	}UG_CATCH_THROW("NLGaussSeidelSolver::apply: "
			"Computation of Start-Defect failed.");

	//	write start defect for debug
	int loopCnt = 0;
	char ext[20]; sprintf(ext, "_iter%03d", loopCnt);
	std::string name("NLGaussSeidel_Defect");
	name.append(ext);
	write_debug(*spD, name.c_str());
	write_debug(u, "NLGaussSeidel_StartSolution");

	// 	start convergence check
	m_spConvCheck->start(*spD);

	// 	assign selector to grid
	TDomain& dom = *m_spApproxSpace->domain();
	typename TDomain::grid_type& grid = *dom.grid();
	m_sel.assign_grid(grid);

	matrix_type& JBlock = m_spJBlock->get_matrix();

	//	loop iteration
	while(!m_spConvCheck->iteration_ended())
	{
		//bool activeSet_changed = false;
		m_spAss->ass_tuner()->set_mapping(&m_map);

		//	loop all indizes
		for (size_t i = 0; i < u.size(); i++)
		{
			// 	Compute Jacobian J(u) using the updated u-components

			//	since we only need J(i,i) and d(i) in every DoF loop,
			//	we access the i-th Element list indicating the neighborhood of DoF i!
			m_sel.clear();
			m_sel.select(m_vElemList[i].begin(), m_vElemList[i].end());

			//	by passing the selector to the assembling the assemble operators
			//	are build up only by looping over the elements which has been selected
			m_spAss->ass_tuner()->set_selector(&m_sel);

			//	assemble only with respect to index i (causes resizing of matrices/vectors)
			m_spAss->ass_tuner()->set_single_index_assembling(i);
			m_map.set_assembling_index(i);

			try{
				NL_GAUSSSEIDEL_PROFILE_BEGIN(NL_GAUSSSEIDELComputeJacobian);
				m_spJBlock->init(u);
				NL_GAUSSSEIDEL_PROFILE_END();
			}UG_CATCH_THROW("NLGaussSeidelSolver::apply: "
					"Initialization of Jacobian failed.");

			//	get i-th block of defect d: d(i) =: m_d_block
			NL_GAUSSSEIDEL_PROFILE_BEGIN(NL_GAUSSSEIDELComputeLastCompDefect);
			m_spAssOp->apply(*spDBlock, u);
			NL_GAUSSSEIDEL_PROFILE_END();

			//	get i,i-th block of J: J(i,i)
			//	depending on the AlgebraType J(i,i) is a 1x1, 2x2, 3x3 Matrix
			//	m_c_i = m_damp * d_i /J_ii
			NL_GAUSSSEIDEL_PROFILE_BEGIN(NL_GAUSSSEIDELInvertJ);
				InverseMatMult((*spC)[0], m_damp, JBlock(0,0) , (*spDBlock)[0]);
			NL_GAUSSSEIDEL_PROFILE_END();

			// 	update i-th block of solution
			u[i] -= (*spC)[0];

			//	should projected GS be performed?
			if (m_bProjectedGS)
			{
				//	call ProjectVectorCorrection

				/*std::vector<DoFIndex> vActiveSet;
				value_type diff;
				diff = u[i] - m_ConsVec[i];
				//	get number of unknowns per value_type
				//	(e.g. if CPU == 3 -> nrFcts = 3!)
				size_t nrFcts = GetSize(diff);

				//	loop fcts
				for (size_t fct = 0; fct < nrFcts; fct++)
				{
					bool penetrate = false;

					if (BlockRef(diff,fct) > 0) //	i.e.: u > m_ConsVec
					{
						bool MultiIndex_is_in_activeSet = false;
						penetrate = true;

						std::vector<DoFIndex>::iterator iter;

						//	adds MultiIndex-pair (i,fct) to vActiveSet
						for (iter = vActiveSet.begin(); iter != vActiveSet.end(); ++iter)
						{
							DoFIndex activeMultiIndex = *iter;

							if (activeMultiIndex[0] == i && activeMultiIndex[1] == fct)
								MultiIndex_is_in_activeSet = true;
							//	TODO: anstatt 'MultiIndex_is_in_activeSet' hier zu verwenden,
							//	mit continue in loop weiterspringen
						}

						if (!MultiIndex_is_in_activeSet)
						{
							DoFIndex activeMultiIndex(i,fct);
							vActiveSet.push_back(activeMultiIndex);
							//activeSet_changed = true;
						}
					}

					if (penetrate)
						BlockRef(u[i],fct) = BlockRef(m_ConsVec[i],fct);

				} //end (fcts)*/
			} //end(m_bProjectedGS)

		} //end(DoFs)

		//	TODO: if active set not changed, iteration converged!
		//if (!activeSet_changed)
		//	break;

		//	set mapping, selector and ass_index to NULL
		m_spAss->ass_tuner()->set_mapping();
		m_spAss->ass_tuner()->set_selector();
		m_spAss->ass_tuner()->disable_single_index_assembling();

		NL_GAUSSSEIDEL_PROFILE_BEGIN(NL_GAUSSSEIDELComputeLastCompDefect);
		m_spAssOp->prepare(u);
		m_spAssOp->apply(*spD, u);
		NL_GAUSSSEIDEL_PROFILE_END();

		//	update counter
		loopCnt++;
		sprintf(ext, "_iter%03d", loopCnt);

		// 	check convergence
		m_spConvCheck->update(*spD);

		//	write defect for debug
		std::string name("NLGaussSeidel_Defect"); name.append(ext);
		write_debug(*spD, name.c_str());
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

#endif /* __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__NL_GAUSS_SEIDEL__NL_GAUSS_SEIDEL_IMPL_H_ */
