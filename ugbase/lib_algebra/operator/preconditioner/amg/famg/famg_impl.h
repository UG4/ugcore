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
#include "../stopwatch.h"
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



} // namespace ug

#endif //  __H__LIB_ALGEBRA__AMG__FAMG_IMPL_H__
