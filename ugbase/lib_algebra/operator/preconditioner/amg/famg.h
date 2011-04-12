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
#include "amg_profiling.h"
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

template<typename matrix_type, typename prolongation_matrix_type, typename vector_type>
class FAMGLevelCalculator;

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

	typedef typename matrix_type::value_type value_type;

	typedef typename TAlgebra::matrix_type prolongation_matrix_type;

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
	bool get_aggressive_coarsening() const { return m_bAggressiveCoarsening; }

	void set_damping_for_smoother_in_interpolation_calculation(double dDampingForSmootherInInterpolationCalculation)
	{
		m_dDampingForSmootherInInterpolationCalculation = dDampingForSmootherInInterpolationCalculation;
	}
	double get_damping_for_smoother_in_interpolation_calculation() const { return m_dDampingForSmootherInInterpolationCalculation; }

	void set_delta(double delta) { m_delta = delta; }
	void set_theta(double theta) { m_theta = theta; }
	double get_delta() const { return m_delta;}
	double get_theta() const { return m_theta;}


	void set_testvector_damps(size_t testvectordamps) { m_iTestvectorDamps = testvectordamps; }
	size_t get_testvector_damps() const { return m_iTestvectorDamps; }


	void reset_testvectors()
	{
		m_testvectors.clear();
		m_vVectorWriters.clear();
		m_omegaVectors.clear();
		m_omegaVectorWriters.clear();
	}

	// void add_std_testvector
	// with callback calculate(vector_type &c, double &weight, stdvector<positions>, stdvector<bool> isinner) o€.

	void add_testvector(vector_type& c, double weight)
	{
		UG_ASSERT(m_testvectors.size() == m_omegaVectors.size(), "???");
		m_testvectors.push_back(c);
		m_omegaVectors.push_back(weight);
	}

	void add_vector_writer(IVectorWriter<vector_type> *vw, double weight)
	{
		m_vVectorWriters.push_back(vw);
		m_omegaVectorWriters.push_back(weight);
	}

	void write_testvectors(bool wt)
	{
		m_writeTestvectors = wt;
	}

private:
//  functions
	virtual void create_AMG_level(matrix_type &AH, prolongation_matrix_type &R, const matrix_type &A,
				prolongation_matrix_type &P, size_t level);

	void c_create_AMG_level(matrix_type &AH, prolongation_matrix_type &R, const matrix_type &A,
			prolongation_matrix_type &P, size_t level);

	void get_testvectors(stdvector<vector_type> &testvectors, stdvector<double> &omega);

private:
// data
	bool m_bAggressiveCoarsening;
	double m_dDampingForSmootherInInterpolationCalculation;
	double m_delta;								///< "Interpolation quality" F may not be worse than this (F < m_delta)
	double m_theta;								///< clip all interpolations with m_theta * F > min F.

	size_t m_iTestvectorDamps;
	bool m_bTestvectorZeroAtDirichlet;
	bool m_writeTestvectors;

	friend class FAMGLevelCalculator<matrix_type, matrix_type, vector_type >;

	stdvector<vector_type> m_testvectors;
	stdvector<double> m_omegaVectors; // testvectorWeights
	stdvector< IVectorWriter<vector_type> * > m_vVectorWriters;
	stdvector<double> m_omegaVectorWriters; // testvectorWeights
};


///	@}

} // namespace ug


#include "famg_impl.h"

#endif // __H__LIB_DISCRETIZATION__AMG_SOLVER__FAMG_H__
