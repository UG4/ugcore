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
#include "stopwatch.h"
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
void famg<TAlgebra>::create_AMG_level(matrix_type &AH, prolongation_matrix_type &R, const matrix_type &A,
		prolongation_matrix_type &P, size_t level)
{
	c_create_AMG_level(AH, R, A, P, level);
}

template<typename TAlgebra>
famg<TAlgebra>::famg() : amg_base<TAlgebra>()
{
	m_theta = 0.95;
	m_delta = 0.5;
	m_dEpsilonTr = 0.3;
	m_dDampingForSmootherInInterpolationCalculation = 0.8;
	m_bAggressiveCoarsening = false;
	m_writeTestvectors = false;
}

template<typename TAlgebra>
famg<TAlgebra>::~famg()
{

}

template<typename TAlgebra>
void famg<TAlgebra>::tostring() const
{
	amg_base<TAlgebra>::tostring();
	UG_LOG("FAMG Preconditioner:\n");

	UG_LOG(" Delta: " << m_delta << " (forces interpolation quality measure F < delta.) " << std::endl);
	UG_LOG(" Theta: " << m_theta << " (forces theta * F < minimimum F in this node) " << std::endl);
	UG_LOG(" Damping for Smoother in interpolation calculation: " << m_dDampingForSmootherInInterpolationCalculation << std::endl);
	UG_LOG(" Aggressive Coarsening is " << (m_bAggressiveCoarsening ? "[ON]\n" : "OFF\n"));
	UG_LOG(" epsilon_tr (truncation of interpolation) = " << m_dEpsilonTr << std::endl);
	UG_LOG(" \n");
	UG_LOG(" testvector is " << (m_bTestvectorZeroAtDirichlet ? "0" : "1") << " at dirichlet nodes" << std::endl);
	UG_LOG(" Nr. of testvector damps: " << m_iTestvectorDamps << std::endl);
	if(m_writeMatrices && m_writeTestvectors)
		UG_LOG(" Write Testvectors is on.\n")
	UG_LOG(" \n");
}


} // namespace ug

#endif //  __H__LIB_ALGEBRA__AMG__FAMG_IMPL_H__
