/**
 * \file amg.h
 *
 * \author Martin Rupp
 *
 * \date 06.08.2010
 *
 * class declaration for amg_base
 *
 * Goethe-Center for Scientific Computing 2009-2010.
 */


#ifndef __H__LIB_DISCRETIZATION__AMG_SOLVER__AMG_H__
#define __H__LIB_DISCRETIZATION__AMG_SOLVER__AMG_H__

#include <vector>
#include <iostream>

#include "amg_base.h"


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
//! algebraic multigrid class.
//!

template <typename TAlgebra>
class amg:
	public amg_base< TAlgebra >
{
public:
	typedef amg_base<TAlgebra> super;
	using super::m_amghelper;
	using super::m_parentIndex;
	using super::m_writeMatrices;
	using super::m_writeMatrixPath;

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
	amg() ;
	virtual ILinearIterator<vector_type,vector_type>* clone()
	{
		amg<algebra_type>* clone = new amg<algebra_type>();
		return dynamic_cast<ILinearIterator<vector_type,vector_type>* >(clone);
	}
	virtual ~amg() { }

	virtual const char* name() const {return "AMGPreconditioner";}

	void set_theta(double new_theta) 		{ m_dTheta = new_theta; }
	double get_theta() const				{ return m_dTheta; }
	void set_sigma(double new_sigma) 		{ m_dSigma = new_sigma; }
	double get_sigma() const				{ return m_dSigma; }

	void set_epsilon(double new_epsilon) 	{ m_dEpsilon = new_epsilon; }
	double get_epsilon() const				{ return m_dEpsilon; }

	void enable_aggressive_coarsening_A(int nrOfPaths)
	{
		m_bAggressiveCoarsening = true; m_iAggressiveCoarseningNrOfPaths = nrOfPaths;
		UG_ASSERT(m_iAggressiveCoarseningNrOfPaths >= 0 && m_iAggressiveCoarseningNrOfPaths < 2, "only A1 and A2 supported, not A" << m_iAggressiveCoarseningNrOfPaths);
	}
	void disable_aggressive_coarsening() 	{ m_bAggressiveCoarsening = false; }
	bool is_aggressive_coarsening() const 	{ return m_bAggressiveCoarsening; }
	bool is_aggressive_coarsening_A(int nrOfPaths) const { return m_bAggressiveCoarsening && m_iAggressiveCoarseningNrOfPaths == nrOfPaths; }

	void tostring() const;

protected:
//  functions
	virtual void create_AMG_level(matrix_type &AH, prolongation_matrix_type &R, const matrix_type &A,
					prolongation_matrix_type &P, size_t level);

// data
	double m_dEpsilon;	///< parameter used for truncation of interpolation
	double m_dTheta; 	///< measure for strong connectivity
	double m_dSigma;

	bool m_bAggressiveCoarsening;				///< true if aggressive coarsening is used on first level
	int m_iAggressiveCoarseningNrOfPaths;  	///<
};
	
	
///	@}

} // namespace ug


#include "amg_impl.h"

#endif // __H__LIB_DISCRETIZATION__AMG_SOLVER__AMG_H__
