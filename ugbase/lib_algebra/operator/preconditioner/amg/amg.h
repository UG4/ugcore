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
	using super::amghelper;
	using super::parentIndex;
	using super::m_writeMatrices;
	using super::m_writeMatrixPath;

public:
//	Algebra type
	typedef TAlgebra algebra_type;

//	Vector type
	typedef typename TAlgebra::vector_type vector_type;

//	Matrix type
	typedef typename TAlgebra::matrix_type matrix_type;

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

	void set_theta(double new_theta) { theta = new_theta; }
	void set_sigma(double new_sigma) { sigma = new_sigma; }

	void enable_aggressive_coarsening_A_2() { aggressiveCoarsening = true; aggressiveCoarseningNrOfPaths = 2;}
	void enable_aggressive_coarsening_A_1() { aggressiveCoarsening = true; aggressiveCoarseningNrOfPaths = 1;}
	void disable_aggressive_coarsening() { aggressiveCoarsening = false; }

	void tostring() const;

protected:
//  functions
	virtual void create_AMG_level(matrix_type &AH, SparseMatrix<double> &R, const matrix_type &A,
							SparseMatrix<double> &P, int level);

// data
	double eps_truncation_of_interpolation;	///< parameter used for truncation of interpolation
	double theta; 							///< measure for strong connectivity
	double sigma;

	bool aggressiveCoarsening;				///< true if aggressive coarsening is used on first level
	int aggressiveCoarseningNrOfPaths;  	///<
};
	
	
///	@}

} // namespace ug


#include "amg_impl.h"

#endif // __H__LIB_DISCRETIZATION__AMG_SOLVER__AMG_H__
