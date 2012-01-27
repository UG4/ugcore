/**
 * \file amg.h
 *
 * \author Martin Rupp
 *
 * \date 06.08.2010
 *
 * class declaration for amg_base
 *
 * Goethe-Center for Scientific Computing 2009-2011.
 *
 *
 */


#ifndef __H__UG__LIB_ALGEBRA__AMG_SOLVER__AMG_H__
#define __H__UG__LIB_ALGEBRA__AMG_SOLVER__AMG_H__

#include <vector>
#include <iostream>

#include "../amg_base.h"
#include "rsamg_nodeinfo.h"

#ifdef UG_PARALLEL
#include "rsamg_parallel_coarsening.h"
#endif

/**
 * \brief Algebraic Multigrid Functions.
 *
 *
 * \defgroup lib_algebra_AMG AMG
 * \ingroup lib_algebra
 */

namespace ug{

/// \addtogroup lib_algebra_AMG
///	@{

// AMG
//---------------------------------
/**
 * algebraic multigrid (AMG) class.
 *
 * scientific references:
 * [AMGKS99]: Algebraic Multigrid (AMG): An Introduction with Applications, K Stüben. GMD Report 70, November 1999
 *
 */
template <typename TAlgebra>
class RSAMG:
	public AMGBase< TAlgebra >
{
public:
	typedef AMGBase<TAlgebra> super;
	using super::m_amghelper;
	using super::m_parentIndex;
	using super::m_writeMatrices;
	using super::m_writeMatrixPath;
	using super::levels;

public:
//	Algebra type
	typedef TAlgebra algebra_type;

//	Vector type
	typedef typename TAlgebra::vector_type vector_type;

//	Matrix type
	typedef typename TAlgebra::matrix_type matrix_type;

	typedef typename TAlgebra::matrix_type prolongation_matrix_type;

	typedef typename matrix_type::value_type value_type;
	
//  functions
	RSAMG() ;
	virtual ILinearIterator<vector_type,vector_type>* clone()
	{
		RSAMG<algebra_type>* clone = new RSAMG<algebra_type>();
		return dynamic_cast<ILinearIterator<vector_type,vector_type>* >(clone);
	}
	virtual ~RSAMG() { }

	virtual const char* name() const {return "AMGPreconditioner";}


	//! 	sets prolongation_strong, a measure for strong connectivity
	void 	set_epsilon_strong(double new_theta) 		{ m_dTheta = new_theta; }
	double 	get_epsilon_strong() const					{ return m_dTheta; }


	//!		sets prolongation_trunction, used in truncation of the interpolation [AMGKS99] 7.2.4
	void 	set_prolongation_truncation(double prolongationTr) 	{ m_dProlongationTr = prolongationTr; }
	double 	get_prolongation_truncation() const				{ return m_dProlongationTr; }


	/**
	 * enables Aggressive Coarsening on the first level
	 * two coarse nodes are strong connected if there exist nrOfPaths of length 2 on the original graph
	 * \param nrOfPaths
	 * \note only 1 and 2 supported at the moment (A1 and A2).
	 * \sa CreateAggressiveCoarseningGraph
	 * [AMGKS99] 7.1.2
	 */
	void enable_aggressive_coarsening_A(int nrOfPaths)
	{
		m_bAggressiveCoarsening = true; m_iAggressiveCoarseningNrOfPaths = nrOfPaths;
		UG_ASSERT(m_iAggressiveCoarseningNrOfPaths == 1 || m_iAggressiveCoarseningNrOfPaths == 2, "only A1 and A2 supported, not A" << m_iAggressiveCoarseningNrOfPaths);
	}
	void disable_aggressive_coarsening() 	{ m_bAggressiveCoarsening = false; }
	bool is_aggressive_coarsening() const 	{ return m_bAggressiveCoarsening; }
	bool is_aggressive_coarsening_A(int nrOfPaths) const { return m_bAggressiveCoarsening && m_iAggressiveCoarseningNrOfPaths == nrOfPaths; }

	void tostring() const;

protected:
//  functions
	virtual void create_AMG_level(matrix_type &AH, prolongation_matrix_type &R, const matrix_type &A,
					prolongation_matrix_type &P, size_t level);

#ifdef UG_PARALLEL
public:
	void set_parallel_coarsening(IParallelCoarsening *pPC)
	{
		m_pParallelCoarsening = pPC;
	}
#endif

// data
	double m_dProlongationTr;	///< parameter used for truncation of interpolation
	double m_dTheta; 		///< measure for strong connectivity

	bool m_bAggressiveCoarsening;				///< true if aggressive coarsening is used on first level
	int m_iAggressiveCoarseningNrOfPaths;  		///<
#ifdef UG_PARALLEL
	IParallelCoarsening *m_pParallelCoarsening; ///<
#endif
};
	
	
///	@}

} // namespace ug


#include "rsamg_impl.h"

#endif // __H__UG__LIB_ALGEBRA__AMG_SOLVER__AMG_H__
