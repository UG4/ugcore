/**
 * \file famg.h
 *
 * \author Martin Rupp
 *
 * \date 16.11.2010
 *
 * class declaration for famg
 *
 * Goethe-Center for Scientific Computing 2010.
 */


#ifndef __H__LIB_DISCRETIZATION__AMG_SOLVER__FAMG_H__
#define __H__LIB_DISCRETIZATION__AMG_SOLVER__FAMG_H__

#include <vector>
#include <iostream>

#include "amg_base.h"



/**
 * \brief Filtering Algebraic Multigrid Functions.
 *
 *
 * \defgroup lib_algebra_FAMG FAMG
 * \ingroup lib_algebra
 */



namespace ug{

/// \addtogroup lib_algebra_FAMG
///	@{

#define FAMG_MAX_LEVELS 32


// AMG
//---------------------------------
//! algebraic multigrid class.
//!

template <typename TAlgebra>
class famg:
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
	famg() ;
	virtual ILinearIterator<vector_type,vector_type>* clone()
	{
		famg<algebra_type>* clone = new famg<algebra_type>();
		return dynamic_cast<ILinearIterator<vector_type,vector_type>* >(clone);
	}
	//	Name of preconditioner
	virtual ~famg();
	void cleanup();

	virtual const char* name() const {return "FAMGPreconditioner";}
	void tostring() const;

	void set_aggressive_coarsening(bool bAggressiveCoarsening) { m_bAggressiveCoarsening = bAggressiveCoarsening; }
	void set_damping_for_smoother_in_interpolation_calculation(double dDampingForSmootherInInterpolationCalculation)
	{
		m_dDampingForSmootherInInterpolationCalculation = dDampingForSmootherInInterpolationCalculation;
	}

	void set_delta(double delta) { m_delta = delta; }
	void set_theta(double theta) { m_theta = theta; }

	void set_testvector_zero_at_dirichlet(bool bZeroAtDirichlet) { m_bTestvectorZeroAtDirichlet = bZeroAtDirichlet; }
	void set_testvector_damps(int testvectordamps) { m_iTestvectorDamps = testvectordamps; }

private:
//  functions
	virtual void create_AMG_level(matrix_type &AH, SparseMatrix<double> &R, const matrix_type &A,
							SparseMatrix<double> &P, int level);

	void c_create_AMG_level(SparseMatrix<double> &AH, SparseMatrix<double> &R, const SparseMatrix<double> &A,
			SparseMatrix<double> &P, int level);

private:
// data
	bool m_bAggressiveCoarsening;
	double m_dDampingForSmootherInInterpolationCalculation;
	double m_delta;								///< "Interpolation quality" F may not be worse than this (F < m_delta)
	double m_theta;								///< clip all interpolations with m_theta * F > min F.

	int m_iTestvectorDamps;
	bool m_bTestvectorZeroAtDirichlet;

};


///	@}

} // namespace ug


#include "famg_impl.h"

#endif // __H__LIB_DISCRETIZATION__AMG_SOLVER__FAMG_H__
