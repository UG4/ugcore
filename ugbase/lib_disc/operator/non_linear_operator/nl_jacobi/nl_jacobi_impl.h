#ifndef __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__NL_JACOBI__NL_JACOBIL_IMPL_H_
#define __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__NL_JACOBI__NL_JACOBIL_IMPL_H_

// extern includes
#include <iostream>

#include "lib_disc/function_spaces/grid_function_util.h"
#include "lib_disc/common/local_algebra.h"
#include "nl_jacobi.h"

#define PROFILE_NL_JACOBI
#ifdef PROFILE_NL_JACOBI
	#define NL_JACOBI_PROFILE_FUNC()		PROFILE_FUNC_GROUP("NL Jacobi")
	#define NL_JACOBI_PROFILE_BEGIN(name)	PROFILE_BEGIN_GROUP(name, "NL Jacobi")
	#define NL_JACOBI_PROFILE_END()		PROFILE_END()
#else
	#define NL_JACOBI_PROFILE_FUNC()
	#define NL_JACOBI_PROFILE_BEGIN(name)
	#define NL_JACOBI_PROFILE_END()
#endif

namespace ug{

template <typename TAlgebra>
NLJacobiSolver<TAlgebra>::
NLJacobiSolver(SmartPtr<IConvergenceCheck<vector_type> > spConvCheck) :
			m_spConvCheck(spConvCheck),
			m_damp(1.0),
			m_dgbCall(0)
{};

template <typename TAlgebra>
NLJacobiSolver<TAlgebra>::
NLJacobiSolver() :
	m_spConvCheck(new StdConvCheck<vector_type>(10, 1e-8, 1e-10, true)),
	m_damp(1.0),
	m_dgbCall(0)
{};

template <typename TAlgebra>
void NLJacobiSolver<TAlgebra>::
set_convergence_check(SmartPtr<IConvergenceCheck<vector_type> > spConvCheck)
{
	m_spConvCheck = spConvCheck;
	m_spConvCheck->set_offset(3);
	m_spConvCheck->set_symbol('#');
	m_spConvCheck->set_name("Nonlinear Jacobi Solver");
}

template <typename TAlgebra>
bool
NLJacobiSolver<TAlgebra>::
init(SmartPtr<IOperator<vector_type> > op)
{
	NL_JACOBI_PROFILE_BEGIN(NL_JACOBISolver_init);
	m_spAssOp = op.template cast_dynamic<AssembledOperator<TAlgebra> >();
	if(m_spAssOp.invalid())
		UG_THROW("NLJacobiSolver: currently only works for AssembledDiscreteOperator.");

	m_spAss = m_spAssOp->discretization();
	if(m_spAss.invalid())
		UG_THROW("AssembledLinearOperator: Assembling routine not set.");

	return true;
}

template <typename TAlgebra>
bool NLJacobiSolver<TAlgebra>::prepare(vector_type& u)
{
	return true;
}


template <typename TAlgebra>
bool NLJacobiSolver<TAlgebra>::apply(vector_type& u)
{
	NL_JACOBI_PROFILE_BEGIN(NL_JACOBISolver_apply);
	//	increase call count
	m_dgbCall++;

//	Jacobian
	if(m_spJ.invalid() || m_spJ->discretization() != m_spAss) {
		m_spJ = make_sp(new AssembledLinearOperator<TAlgebra>(m_spAss));
		m_spJ->set_level(m_spAssOp->level());
	}

//	create tmp vectors
	SmartPtr<vector_type> spD = u.clone_without_values();
	SmartPtr<vector_type> spC = u.clone_without_values();

//	Set dirichlet values
	try{
		m_spAssOp->prepare(u);
	}
	UG_CATCH_THROW("NLJacobiSolver::apply: Prepare of Operator failed.");

// 	Compute first Defect d = L(u)
	try{
		NL_JACOBI_PROFILE_BEGIN(NL_JACOBIComputeDefect1);
		m_spAssOp->apply(*spD, u);
		NL_JACOBI_PROFILE_END();
	}UG_CATCH_THROW("NLJacobiSolver::apply: "
			"Computation of Start-Defect failed.");

//	write start defect for debug
	int loopCnt = 0;
	char ext[20]; sprintf(ext, "_iter%03d", loopCnt);
	std::string name("NLJacobi_Defect");
	name.append(ext);
	write_debug(*spD, name.c_str());
	write_debug(u, "NLJacobi_StartSolution");

// 	start convergence check
	m_spConvCheck->start(*spD);

	matrix_type& J = m_spJ->get_matrix();

//	loop iteration
	while(!m_spConvCheck->iteration_ended())
	{
		// 	set correction c = 0
		NL_JACOBI_PROFILE_BEGIN(NL_JACOBISetCorretionZero);
		spC->set(0.0);
		NL_JACOBI_PROFILE_END();

		// 	Compute Jacobian J(u)
		try{
			NL_JACOBI_PROFILE_BEGIN(NL_JACOBIComputeJacobian);
			m_spJ->init(u);
			NL_JACOBI_PROFILE_END();
		}UG_CATCH_THROW("NLJacobiSolver::apply: "
				"Initialization of Jacobian failed.");

		//	Write Jacobian for debug
		std::string matname("NLJacobi_Jacobian");
		matname.append(ext);
		write_debug(m_spJ->get_matrix(), matname.c_str());

		NL_JACOBI_PROFILE_BEGIN(NL_JACOBIInvertBlocks);

		//	loop all indizes
		for (size_t i = 0; i < u.size(); i++)
		{
			//	get i,i-th block of J: J(i,i)
			//	m_c_i = m_damp * d_i /J_ii
			InverseMatMult((*spC)[i], m_damp, J(i,i), (*spD)[i]);

			// 	update solution
			u[i] -= (*spC)[i];
		}
		NL_JACOBI_PROFILE_END();

	// 	compute new Defect
		NL_JACOBI_PROFILE_BEGIN(NL_JACOBIComputeDefect);
		m_spAssOp->prepare(u);
		m_spAssOp->apply(*spD, u);
		NL_JACOBI_PROFILE_END();

	//	update counter
		loopCnt++;
		sprintf(ext, "_iter%03d", loopCnt);

	// 	check convergence
		m_spConvCheck->update(*spD);

	//	write defect for debug
		std::string name("NLJacobi_Defect"); name.append(ext);
		write_debug(*spD, name.c_str());
	}

	return m_spConvCheck->post();
}

template <typename TAlgebra>
void NLJacobiSolver<TAlgebra>::write_debug(const vector_type& vec, const char* filename)
{
//	add iter count to name
	std::string name(filename);
	char ext[20]; sprintf(ext, "_call%03d", m_dgbCall);
	name.append(ext).append(".vec");

//	write
	base_writer_type::write_debug(vec, name.c_str());
}

template <typename TAlgebra>
void NLJacobiSolver<TAlgebra>::write_debug(const matrix_type& mat, const char* filename)
{
//	add iter count to name
	std::string name(filename);
	char ext[20]; sprintf(ext, "_call%03d", m_dgbCall);
	name.append(ext).append(".mat");

//	write
	base_writer_type::write_debug(mat, name.c_str());
}

}

#endif /* __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__NL_JACOBI__NL_JACOBIL_IMPL_H_ */
