/**
 * \file amg_impl.h
 *
 * \author Martin Rupp
 *
 * \date 16.11.2010
 *
 * implementation file for amg
 *
 * Goethe-Center for Scientific Computing 2009-2010.
 */

#ifndef __H__LIB_ALGEBRA__AMG__FAMG_IMPL_H__
#define __H__LIB_ALGEBRA__AMG__FAMG_IMPL_H__

//#include "sparsematrix_util.h"

//#include "famg_nodeinfo.h"
#include "common/stopwatch.h"
#include "common/assert.h"
//#include "maxheap.h"

namespace ug{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// createFAMGLevel:
//-------------------------
/**
 * create FAMG matrix R, P, and AH = R A P
 * \param AH
 * \param R
 * \param A
 * \param P
 * \param level
 */
template<typename TAlgebra>
void FAMG<TAlgebra>::create_AMG_level(matrix_type &AH, prolongation_matrix_type &R, const matrix_type &A,
		prolongation_matrix_type &P, size_t level)
{
	c_create_AMG_level(AH, R, A, P, level);
}

template<typename TAlgebra>
FAMG<TAlgebra>::FAMG() : AMGBase<TAlgebra>()
{
	m_theta = 0.95;
	m_delta = 0.5;

	m_dDampingForSmootherInInterpolationCalculation = 0.8;
	m_bAggressiveCoarsening = false;
	m_writeTestvectors = false;
	m_bTestvectorsFromMatrixRows = false;
	m_bWriteFValues = false;

	m_dProlongationTruncation = 0.0;
	m_dHReduceInterpolationNodesParameter = 0.0;
	m_dPrereduceAToStrongParameter = 0.0;
	m_dGalerkinTruncation = 1e-12;

	m_dStrongConnectionExternal = 0.1;


	m_bExternalCoarsening = false;
	m_bUsePrecalculate = true;

	iDebugLevelOverlapAMG = iDebugLevelOverlapMatrix = iDebugLevelTestvectorCalc = iDebugLevelPhase3 =
	iDebugLevelCalculateParentPairs = iDebugLevelColoring = iDebugLevelRecvCoarsening = iDebugLevelGetRatings = iDebugLevelPrecalculateCoarsening =
				iDebugLevelAggressiveCoarsening = iDebugLevelSendCoarsening = iDebugLevelCommunicateProlongation = iDebugLevelAfterCommunicateProlongation = 0;

	m_bWriteFValues = false;

	m_bProjectedEVP = false;
#ifdef UG_PARALLEL
	m_pParallelCoarsening = NULL;
#endif
}

template<typename TAlgebra>
FAMG<TAlgebra>::~FAMG()
{

}

template<typename TAlgebra>
void FAMG<TAlgebra>::tostring() const
{
	AMGBase<TAlgebra>::tostring();
	UG_LOG("FAMG Preconditioner:\n");

	UG_LOG(" Delta: " << m_delta << " (forces interpolation quality measure F < delta.) " << std::endl);
	UG_LOG(" Theta: " << m_theta << " (with multiple parents pairs, discard pairs with m_theta * F > min F.) " << std::endl);
	UG_LOG(" Damping for Smoother in interpolation calculation: " << m_dDampingForSmootherInInterpolationCalculation << std::endl);
	UG_LOG(" Aggressive Coarsening is " << (m_bAggressiveCoarsening ? "[ON]" : "OFF"));
	UG_LOG(", external Coarsening is " << (m_bExternalCoarsening ? "[ON]" : "OFF"));
	UG_LOG(", precalculate Coarsening is " << (m_bUsePrecalculate ? "[ON]" : "OFF"));
	UG_LOG(", truncation of interpolation = " << m_dProlongationTruncation << "\n");
	if(m_bExternalCoarsening)
		UG_LOG("Strong Connection External Coarsening: " << m_dStrongConnectionExternal << "\n");
	UG_LOG("H-Reduce Interpolation Nodes Parameter:" << m_dHReduceInterpolationNodesParameter << "\n");
	UG_LOG("Galerkin Truncation: " << m_dGalerkinTruncation << "\n");
	UG_LOG("prereduce A parameter: " << m_dPrereduceAToStrongParameter << "\n");

	UG_LOG(" \n");
	UG_LOG(m_testvectors.size() + m_vVectorWriters.size() << " test vectors. Nr. of testvector damps: " << m_iTestvectorDamps << std::endl);
	if(m_writeMatrices && m_writeTestvectors)
		UG_LOG(" Write Testvectors is on.\n")
	UG_LOG(" \n");
}


template<typename TAlgebra>
void FAMG<TAlgebra>::get_testvectors(stdvector<vector_type> &testvectors, stdvector<double> &omega)
{
	AMG_PROFILE_FUNC();

	testvectors.resize(m_vVectorWriters.size() + m_testvectors.size());
	omega.resize(m_vVectorWriters.size() + m_testvectors.size());

	for(size_t i=0; i<m_testvectors.size(); i++)
	{
		testvectors[i] = m_testvectors[i];
		omega[i] = m_omegaVectors[i];
	}

	for(size_t i=0; i<m_vVectorWriters.size(); i++)
	{
		size_t index = i+m_testvectors.size();
		m_vVectorWriters[i]->update(testvectors[index]);
		omega[index] = m_omegaVectorWriters[i];

	}
}

template<typename TAlgebra>
void FAMG<TAlgebra>::get_testvectors_from_matrix_rows
	(const matrix_type &A, stdvector<vector_type> &testvectors, stdvector<double> &omega)
{
	AMG_PROFILE_FUNC();

	testvectors.resize(1);
	omega.resize(1);
	vector_type &v = testvectors[0];
	v.resize(A.num_rows());
	for(size_t i=0; i<A.num_rows(); i++)
	{
		if(A.is_isolated(i))
			v[i] = 0.0;
		else
			v[i] = 1.0;
	}
}


template<typename TAlgebra>
void FAMG<TAlgebra>::precalc_level(size_t level)
{

	AMG_PROFILE_FUNC();
	//UG_LOG("\n\n\nprecalc!\n\n\n");
	UG_ASSERT(m_testvectorsmoother != NULL, "please provide a testvector smoother.");


	if(level == 0)
	{
		testvectors.clear();
		if(m_bTestvectorsFromMatrixRows)
			get_testvectors_from_matrix_rows(*levels[level]->pAgglomeratedA, testvectors, omega);
		else
			get_testvectors(testvectors, omega);
	}

#ifdef UG_PARALLEL
	if((level != 0 || m_bTestvectorsFromMatrixRows == false)
			&& levels[level]->bHasBeenMerged)
	{
		vector_type t;
		if(super::isMergingMaster(level))
			t.resize(levels[level]->pAgglomeratedA->num_rows());
		for(size_t i=0; i<testvectors.size(); i++)
		{
			t.set_storage_type(PST_CONSISTENT);
			testvectors[i].set_storage_type(PST_CONSISTENT);
			super::gather_vertical(testvectors[i], t, level, PST_CONSISTENT);
			testvectors[i] = t;
		}
	}
#endif
}

template<typename TAlgebra>
void FAMG<TAlgebra>::c_create_AMG_level(matrix_type &AH, prolongation_matrix_type &R, const matrix_type &A,
		prolongation_matrix_type &P, size_t level)
{
	AMG_PROFILE_FUNC();

	UG_ASSERT(m_testvectorsmoother != NULL, "please provide a testvector smoother.");

	if(m_bProjectedEVP && level == 0)
	{

	}
	UG_ASSERT(testvectors.size() != 0, "testvectors?");
	UG_ASSERT(testvectors[0].size() == A.num_rows(), "testvectorsize = " << testvectors[0].size() << " != A.num_rows() = " << A.num_rows() << " ?");

	if(m_writeMatrices && m_writeTestvectors)
	{
		for(size_t i=0; i<testvectors.size(); i++)
			WriteVectorToConnectionViewer((m_writeMatrixPath + ToString("testvector") + ToString(i) + ToString("_L") + ToString(level) + ".vec").c_str(),
					testvectors[i], &m_amghelper.positions[level][0], m_dbgDimension);
	}

	//UG_ASSERT(testvectors.size() > 0, "we need at least one testvector.");

	// testvectors will be altered by FAMGLevelCalculator

	// UG_SET_DEBUG_LEVEL(LIB_ALG_AMG, 4);

	if(m_dPrereduceAToStrongParameter > 0.0)
	{
		matrix_type Aeps;
		ReduceToStrongConnections(Aeps, A, m_dPrereduceAToStrongParameter);

		FAMGLevelCalculator<matrix_type, prolongation_matrix_type, vector_type> dummy(*this, AH, R, Aeps, A, P, level, testvectors, omega);
		dummy.do_calculation();
	}
	else
	{
		FAMGLevelCalculator<matrix_type, prolongation_matrix_type, vector_type> dummy(*this, AH, R, A, A, P, level, testvectors, omega);
		dummy.do_calculation();
	}

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//       check_testvector
//--------------------------------
/**
 */
template<typename TAlgebra>
bool FAMG<TAlgebra>::check_testvector(size_t fromlevel, size_t tolevel)
{
	AMG_PROFILE_FUNC();

	UG_LOG("\n");
	UG_LOG("            check_testvector\n");
	UG_LOG("==========================================\n")
	UG_LOG("check if testvectors are interpolated exactly\n")
	if (0){
		int level = super::m_usedLevels-1;
		matrix_type &A = *super::levels[level]->pA;
		for(size_t i=0; i<A.num_rows(); i++)
		{
			typename matrix_type::value_type sum=0;
			for(typename matrix_type::row_iterator it = A.begin_row(i); it != A.end_row(i); ++it)
				sum += it.value();
			if(BlockNorm(sum) > 0.1)
			{
				UG_LOG("Row " << i << " has sum " << sum << ": ");
				A.pr(i);
				UG_LOG("\n");
			}
		}
	}

	size_t preSmooth = super::get_num_presmooth();
	size_t postSmooth = super::get_num_postsmooth();
	super::set_num_presmooth(get_testvector_smooths());
	super::set_num_postsmooth(get_testvector_smooths());

	for(size_t level=0; ; level++)
	{
		if(level > tolevel) break;
		// get testvectors, get agglomerated tv
		precalc_level(level);
#ifdef UG_PARALLEL
		if(super::isMergingSlave(level))
		{
			UG_LOG("merged.\n");
			break;
		}
#endif
		if(level >=fromlevel)
		{
			UG_LOG("Check Level " << level << "\n-----------------------------------\n")
		}

		SmartPtr<matrix_type> pA, pAH;
		pA = super::levels[level]->pAgglomeratedA;
		pAH = super::levels[level+1]->pAgglomeratedA;

		matrix_type &A = *pA;
		matrix_type &AH = *pAH;
		matrix_type &R = super::levels[level]->R;
		matrix_type &P = super::levels[level]->P;
		vector_type c, d, tv;

	#ifdef UG_PARALLEL
		SetParallelVectorAsMatrix(c, A, PST_CONSISTENT);
		SetParallelVectorAsMatrix(d, A, PST_CONSISTENT);
		SetParallelVectorAsMatrix(tv, A, PST_CONSISTENT);
	#endif
		tv.resize(A.num_rows());
		d.resize(A.num_rows());
		c.resize(A.num_rows());
		for(size_t i=0; i<d.size(); i++)
			tv[i] = testvectors[0][i];
		A.apply(d, tv);
		c.set(0.0);
		typename super::checkResult res;

		/*PRINTLAYOUT(A.get_process_communicator(), A.get_communicator(), d.get_master_layout(), d.get_slave_layout());
		d.change_storage_type(PST_CONSISTENT);

		UG_LOG("\n\ntv level " << level << "\n");
		tv.print();
		UG_LOG("\n\nA level " << level << "\n");
		A.print();
		UG_LOG("\n\ntestvectors level " << level << "\n");
		testvectors[0].print();
		UG_LOG("\n\nd level " << level << "\n");
		d.print();

		d.change_storage_type(PST_ADDITIVE);*/

		if(level >=fromlevel)
			super::check_level(c, d, A, level, res, &tv);

		if(level+1 < super::m_usedLevels-1)
		{
			matrix_type Aeps;
			ReduceToStrongConnections(Aeps, A, m_dPrereduceAToStrongParameter);
			FAMGLevelCalculator<matrix_type, prolongation_matrix_type, vector_type> dummy(*this, AH, R, Aeps, A, P, level, testvectors, omega);
			dummy.onlyTV();
		}
		else break;

	}
	super::set_num_presmooth(preSmooth);
	super::set_num_postsmooth(postSmooth);
	return true;
}


template<typename TAlgebra>
void FAMG<TAlgebra>::init_for_EVP(const vector_type &u, const vector_type &b)
{
	CloneVector(evpU, u);
	CloneVector(evpB, b);
	m_bProjectedEVP=true;
}


} // namespace ug

#include "famg_level_calculator.h"
#endif //  __H__LIB_ALGEBRA__AMG__FAMG_IMPL_H__

