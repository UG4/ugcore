/*
 * sqp_impl.h
 *
 *  Created on: 11.01.2012
 *      Author: Raphael Prohl
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__SQP_METHOD__SQP_IMPL__
#define __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__SQP_METHOD__SQP_IMPL__

#include <iostream>
#include <sstream>

#include "sqp.h"
#include "lib_disc/function_spaces/grid_function_util.h"

#include "lib_disc/common/groups_util.h"
#include "sqp_elem_util.h"

#define PROFILE_SQP
#ifdef PROFILE_SQP
	#define SQP_PROFILE_FUNC()		PROFILE_FUNC()
	#define SQP_PROFILE_BEGIN(name)	PROFILE_BEGIN(name)
	#define SQP_PROFILE_END()		PROFILE_END()
#else
	#define SQP_PROFILE_FUNC()
	#define SQP_PROFILE_BEGIN(name)
	#define SQP_PROFILE_END()
#endif

namespace ug{

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
bool
SQPMethod<TDomain,TDoFDistribution, TAlgebra>::
init()
{
	return true;
}


template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
bool
SQPMethod<TDomain,TDoFDistribution, TAlgebra>::
prepare()
{
//	Check if Tolerance Check has been set
	if(m_tolerance_check == 0.0)
	{
		UG_LOG("ERROR in 'SQPMethod::prepare': Tolerance Check not set.\n");
		return false;
	}

	return true;
}

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
bool
SQPMethod<TDomain,TDoFDistribution, TAlgebra>::
check_tolerance(const vector_type& u, const dof_distribution_type& dd)
{

	//soll auf ElemDisc (FE1NonlinearElasticity) - Methode zugreifen,
	//die elementweise über die attached Data läuft

	//	update the elem discs
	/*	if(!update_disc_items())
			UG_THROW_FATAL("Cannot update disc items.");*/

	//	Union of Subsets
		SubsetGroup unionSubsets;
		std::vector<SubsetGroup> vSSGrp;

	//	create list of all subsets
		const ISubsetHandler& sh = dd.get_function_pattern().subset_handler();
		if(!CreateSubsetGroups(vSSGrp, unionSubsets, m_vElemDisc, sh))
			UG_THROW_FATAL("ERROR in 'SQPMethod':"
					" Can not Subset Groups and union.\n");

	//	loop subsets
		for(size_t i = 0; i < unionSubsets.num_subsets(); ++i)
		{
		//	get subset
			const int si = unionSubsets[i];

		//	get dimension of the subset
			const int dim = unionSubsets.dim(i);

		//	request if subset is regular grid
			bool bNonRegularGrid = !unionSubsets.regular_grid(i);

		//	overrule by regular grid if required
			if(m_bForceRegGrid) bNonRegularGrid = false;

		//	Elem Disc on the subset
			std::vector<IElemDisc*> vSubsetElemDisc;

		//	get all element discretizations that work on the subset
			GetElemDiscOnSubset(vSubsetElemDisc, m_vElemDisc, vSSGrp, si);

		//	success flag
			bool bSuc = true;

		//	assemble on suitable elements
			switch(dim)
			{
			case 1:
				bSuc &= CheckTolerance<Edge,TDoFDistribution,TAlgebra>
					(vSubsetElemDisc, dd, si, bNonRegularGrid, u);
				break;
			case 2:
				bSuc &= CheckTolerance<Triangle,TDoFDistribution,TAlgebra>
					(vSubsetElemDisc, dd, si, bNonRegularGrid, u);
				bSuc &= CheckTolerance<Quadrilateral,TDoFDistribution,TAlgebra>
					(vSubsetElemDisc, dd, si, bNonRegularGrid, u);
				break;
			case 3:
				bSuc &= CheckTolerance<Tetrahedron,TDoFDistribution,TAlgebra>
					(vSubsetElemDisc, dd, si, bNonRegularGrid, u);
				bSuc &= CheckTolerance<Pyramid,TDoFDistribution,TAlgebra>
					(vSubsetElemDisc, dd, si, bNonRegularGrid, u);
				bSuc &= CheckTolerance<Prism,TDoFDistribution,TAlgebra>
					(vSubsetElemDisc, dd, si, bNonRegularGrid, u);
				bSuc &= CheckTolerance<Hexahedron,TDoFDistribution,TAlgebra>
					(vSubsetElemDisc, dd, si, bNonRegularGrid, u);
				break;
			default:
				UG_THROW_FATAL("ERROR in 'SQPMethod::check_tolerance':"
							"Dimension "<<dim<<" (subset="<<si<<") not supported.\n");
			}

		//	check success
			if(!bSuc)
				UG_THROW_FATAL("ERROR in 'SQPMethod::check_tolerance':"
								" Assembling of elements of Dimension " << dim << " in "
								" subset "<<si<< " failed.\n");
		}

	return true;
}

/*template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
bool
SQPMethod<TDomain,TDoFDistribution, TAlgebra>::
update(vector_type& u)
{
//	increase call count
	m_dgbCall++;

// 	Compute first Defect
	NEWTON_PROFILE_BEGIN(NewtonComputeDefect1);
	if(m_N->apply(m_d, u) != true)
	{
		UG_LOG("ERROR in 'SQPMethod::apply': Cannot apply Non-linear"
				" Operator to compute start defect.\n");
		return false;
	}
	NEWTON_PROFILE_END();

//	write start defect for debug
	int loopCnt = 0;
	char ext[20]; sprintf(ext, "_iter%03d", loopCnt);
	std::string name("NEWTON_Defect");
	name.append(ext);
	write_debug(m_d, name.c_str());
	write_debug(u, "NEWTON_StartSolution");

// 	increase offset of output for linear solver
	IConvergenceCheck* pLinConvCheck = m_pLinearSolver->get_convergence_check();
	int iLinSolverOffset = 0;
	if(pLinConvCheck != NULL)
	{
		iLinSolverOffset = pLinConvCheck->get_offset();
		pLinConvCheck->set_offset( m_pConvCheck->get_offset() + 3);
	}

// 	set info string indicating the used linear solver
	std::stringstream ss; ss << "(Linear Solver: " << m_pLinearSolver->name() << ")";
	m_pConvCheck->set_info(ss.str());

// 	copy pattern
	vector_type s; s.resize(u.size()); s = u;

// 	start convergence check
	m_pConvCheck->start(m_d);

//	loop iteration
	while(!m_pConvCheck->iteration_ended())
	{
	// 	set c = 0
		NEWTON_PROFILE_BEGIN(NewtonSetCorretionZero);
		if(!m_c.set(0.0))
		{
			UG_LOG("ERROR in 'SQPMethod::apply':"
					" Cannot reset correction to zero.\n");
			return false;
		}
		NEWTON_PROFILE_END();

	// 	Compute Jacobian
		NEWTON_PROFILE_BEGIN(NewtonComputeJacobian);
		if(!m_J->init(u))
		{
			UG_LOG("ERROR in 'SQPMethod::apply':"
					" Cannot prepare Jacobi Operator.\n");
			return false;
		}
		NEWTON_PROFILE_END();

	//	Write Jacobian for debug
		std::string matname("NEWTON_Jacobian");
		matname.append(ext);
		write_debug(m_J->get_matrix(), matname.c_str());

	// 	Init Jacobi Inverse
		NEWTON_PROFILE_BEGIN(NewtonPrepareLinSolver);
		if(!m_pLinearSolver->init(*m_J, u))
		{
			UG_LOG("ERROR in 'SQPMethod::apply': Cannot init Inverse Linear "
					"Operator for Jacobi-Operator.\n");
			return false;
		}
		NEWTON_PROFILE_END();

	// 	Solve Linearized System
		NEWTON_PROFILE_BEGIN(NewtonApplyLinSolver);
		if(!m_pLinearSolver->apply(m_c, m_d))
		{
			UG_LOG("ERROR in 'SQPMethod::apply': Cannot apply Inverse Linear "
					"Operator for Jacobi-Operator.\n");
			return false;
		}
		NEWTON_PROFILE_END();

	// 	Line Search
		if(m_pLineSearch != NULL)
		{
			m_pLineSearch->set_offset("   #  ");
			NEWTON_PROFILE_BEGIN(NewtonLineSearch);
			if(!m_pLineSearch->search(*m_N, u, m_c, m_d, m_pConvCheck->defect()))
			{
				UG_LOG("ERROR in 'SQPMethod::apply': "
						"Newton Solver did not converge.\n");
				return false;
			}
			NEWTON_PROFILE_END();
		}
	// 	No line search: Compute new defect
		else
		{
		// 	update solution
			u -= m_c;

		// 	compute new Defect
			NEWTON_PROFILE_BEGIN(NewtonComputeDefect);
			if(!m_N->prepare(m_d, u))
			{
				UG_LOG("ERROR in 'SQPMethod::apply': Cannot prepare Non-linear"
						" Operator for defect computation.\n");
				return false;
			}
			if(!m_N->apply(m_d, u))
			{
				UG_LOG("ERROR in 'SQPMethod::apply': Cannot apply Non-linear "
						"Operator to compute defect.\n");
				return false;
			}
			NEWTON_PROFILE_END();
		}

	//	update counter
		loopCnt++;
		sprintf(ext, "_iter%03d", loopCnt);

	//	write defect for debug
		std::string name("NEWTON_Defect"); name.append(ext);
		write_debug(m_d, name.c_str());
		std::string name2("NEWTON_Correction"); name2.append(ext);
		write_debug(m_c, name2.c_str());

	// 	check convergence
		m_pConvCheck->update(m_d);
	}

	// reset offset of output for linear solver to previous value
	if(pLinConvCheck != NULL)
	{
		pLinConvCheck->set_offset(iLinSolverOffset);
	}

	return m_pConvCheck->post();
}*/

}

#endif /* __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__SQP_METHOD__SQP_IMPL__ */

