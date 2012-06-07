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


#ifndef __H__UG__LIB_DISC__AMG_SOLVER__FAMG_H__
#define __H__UG__LIB_DISC__AMG_SOLVER__FAMG_H__

#include <vector>
#include <iostream>

#include "../amg_base.h"
#include "../amg_profiling.h"
#ifdef UG_PARALLEL
#include "../rsamg/rsamg_parallel_coarsening.h"
#endif
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
class FAMG:
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

	typedef typename matrix_type::value_type value_type;

	typedef typename TAlgebra::matrix_type prolongation_matrix_type;

//  functions
	FAMG() ;
	virtual SmartPtr<ILinearIterator<vector_type> > clone()
	{
		SmartPtr<FAMG<algebra_type> > clone(new FAMG<algebra_type>);
		return clone;
	}
	//	Name of preconditioner
	virtual ~FAMG();
	void cleanup();

	virtual const char* name() const {return "FAMGPreconditioner";}
	virtual void tostring() const;

	void set_aggressive_coarsening(bool bAggressiveCoarsening) { m_bAggressiveCoarsening = bAggressiveCoarsening; }
	bool get_aggressive_coarsening() const { return m_bAggressiveCoarsening; }

	void set_damping_for_smoother_in_interpolation_calculation(double dDampingForSmootherInInterpolationCalculation)
	{
		m_dDampingForSmootherInInterpolationCalculation = dDampingForSmootherInInterpolationCalculation;
	}


	void set_external_coarsening(bool bUse) { m_bExternalCoarsening = bUse; }
	void set_use_precalculate(bool bUse) { m_bUsePrecalculate = bUse; }

#ifdef UG_PARALLEL
public:
	void set_parallel_coarsening(IParallelCoarsening *pPC)
	{
		m_pParallelCoarsening = pPC;
	}
#endif


	size_t iDebugLevelOverlapAMG, iDebugLevelOverlapMatrix, iDebugLevelTestvectorCalc, iDebugLevelPhase3,
			iDebugLevelCalculateParentPairs, iDebugLevelColoring, iDebugLevelRecvCoarsening, iDebugLevelGetRatings, iDebugLevelPrecalculateCoarsening,
			iDebugLevelAggressiveCoarsening, iDebugLevelSendCoarsening, iDebugLevelCommunicateProlongation, iDebugLevelAfterCommunicateProlongation;
	void set_debug_level_overlap(size_t amg, size_t matrix) { iDebugLevelOverlapAMG = amg; iDebugLevelOverlapMatrix = matrix; }
	void set_debug_level_testvector_calc(size_t lvl) { iDebugLevelTestvectorCalc = lvl; }
	void set_debug_level_phase_3(size_t lvl) { iDebugLevelPhase3 = lvl; }
	void set_debug_level_calculate_parent_pairs(size_t lvl) { iDebugLevelCalculateParentPairs = lvl; }
	void set_debug_level_coloring(size_t lvl) { iDebugLevelColoring = lvl; }
	void set_debug_level_recv_coarsening(size_t lvl) { iDebugLevelRecvCoarsening = lvl; }
	void set_debug_level_get_ratings(size_t lvl) { iDebugLevelGetRatings = lvl; }
	void set_debug_level_precalculate_coarsening(size_t lvl) { iDebugLevelPrecalculateCoarsening = lvl; }
	void set_debug_level_aggressive_coarsening(size_t lvl) { iDebugLevelAggressiveCoarsening = lvl; }
	void set_debug_level_send_coarsening(size_t lvl) { iDebugLevelSendCoarsening = lvl; }
	void set_debug_level_communicate_prolongation(size_t lvl) { iDebugLevelCommunicateProlongation = lvl; }
	void set_debug_level_after_communciate_prolongation(size_t lvl) { iDebugLevelAfterCommunicateProlongation = lvl; }


	double get_damping_for_smoother_in_interpolation_calculation() const { return m_dDampingForSmootherInInterpolationCalculation; }

	void set_delta(double delta) { m_delta = delta; }
	void set_theta(double theta) { m_theta = theta; }
	double get_delta() const { return m_delta;}
	double get_theta() const { return m_theta;}



	void set_testvector_smooths(size_t testvectordamps) { m_iTestvectorDamps = testvectordamps; }
	size_t get_testvector_smooths() const { return m_iTestvectorDamps; }


	void reset_testvectors()
	{
		m_testvectors.clear();
		m_vVectorWriters.clear();
		m_omegaVectors.clear();
		m_omegaVectorWriters.clear();
	}

	// void add_std_testvector
	// with callback calculate(vector_type &c, double &weight, stdvector<positions>, stdvector<bool> isinner) o�.

	void add_testvector(vector_type& c, double weight)
	{
		UG_ASSERT(m_testvectors.size() == m_omegaVectors.size(), "???");
		m_testvectors.push_back(c);
		m_omegaVectors.push_back(weight);
	}

	void add_testvector(IVectorWriter<vector_type> *vw, double weight)
	{
		m_vVectorWriters.push_back(vw);
		m_omegaVectorWriters.push_back(weight);
	}

	void set_write_testvectors(bool wt)
	{
		m_writeTestvectors = wt;
	}

	void set_testvector_from_matrix_rows(bool bEnable)
	{
		m_bTestvectorsFromMatrixRows = bEnable;
	}

	void init_for_EVP(const vector_type &u, const vector_type &b);


	void set_testvector_smoother(SmartPtr<ILinearIterator<vector_type> > testvectorsmoother) { m_testvectorsmoother = testvectorsmoother; }

	//!		sets epsilon_trunction, used in truncation of the interpolation [AMGKS99] 7.2.4
	//!		prolongation truncation: set all P(i,j) = 0 if they are too small compared to greatest P(i,.)
	void 	set_prolongation_truncation(double prolongationTr) 	{ m_dProlongationTruncation = prolongationTr; }
	double 	get_prolongation_truncation() const				{ return m_dProlongationTruncation; }

	void 	set_prereduce_A_parameter(double d) 	{	m_dPrereduceAToStrongParameter = d; }
	double 	get_prereduce_A_parameter() const		{	return m_dPrereduceAToStrongParameter; }

	void 	set_H_reduce_interpolation_nodes_parameter(double d)	{	m_dHReduceInterpolationNodesParameter = d; 	}
	double	get_H_reduce_interpolation_nodes_parameter() const		{	return m_dHReduceInterpolationNodesParameter; }

	void	set_galerkin_truncation(double d)	{ m_dGalerkinTruncation = d; }
	double	get_galerkin_truncation() const		{ return m_dGalerkinTruncation; }

	void	set_strong_connection_external(double d) 	{ m_dStrongConnectionExternal = d; }
	double	get_strong_connection_external() 			{ return m_dStrongConnectionExternal; }

	void	set_write_f_values(bool b)	{ m_bWriteFValues = b; }

	bool check_testvector(size_t fromlevel, size_t tolevel);
	bool check_testvector(size_t fromlevel)
	{
		return check_testvector(fromlevel, AMGBase<TAlgebra>::m_usedLevels);
	}
	bool check_testvector()
	{
		return check_testvector(0);
	}

private:
//  functions
	virtual void create_AMG_level(matrix_type &AH, prolongation_matrix_type &R, const matrix_type &A,
				prolongation_matrix_type &P, size_t level);

	void c_create_AMG_level(matrix_type &AH, prolongation_matrix_type &R, const matrix_type &A,
			prolongation_matrix_type &P, size_t level);

	virtual void precalc_level(size_t level);

	void get_testvectors(const matrix_type &A, stdvector<vector_type> &testvectors, stdvector<double> &omega);
	void get_testvectors(stdvector<vector_type> &testvectors, stdvector<double> &omega);
	void get_testvectors_from_matrix_rows
			(const matrix_type &A, stdvector<vector_type> &testvectors, stdvector<double> &omega);

private:
// data
	double m_dDampingForSmootherInInterpolationCalculation;
	double m_delta;				///< "Interpolation quality" F may not be worse than this (F < m_delta)
	double m_theta;				///< with multiple parents, discard pairs with m_theta * F > min F.

	double m_dProlongationTruncation;		///< parameter used for truncation of interpolation

	double m_dHReduceInterpolationNodesParameter; ///< parameter used to reduce the number of interpolation nodes by looking at the H(i,j) value.

	double m_dPrereduceAToStrongParameter;	///< parameter used to reduce the matrix A to
	double m_dGalerkinTruncation;

	double m_dStrongConnectionExternal;


	size_t m_iTestvectorDamps;
	bool m_writeTestvectors;


	bool m_bAggressiveCoarsening;
	bool m_bExternalCoarsening;
	bool m_bUsePrecalculate;
	bool m_bTestvectorsFromMatrixRows;
	bool m_bWriteFValues;

	SmartPtr<ILinearIterator<vector_type> > m_testvectorsmoother;

	friend class FAMGLevelCalculator<matrix_type, matrix_type, vector_type >;

	stdvector<vector_type> m_testvectors;
	stdvector<double> m_omegaVectors; // testvectorWeights
	stdvector< IVectorWriter<vector_type> * > m_vVectorWriters;
	stdvector<double> m_omegaVectorWriters; // testvectorWeights

	stdvector< vector_type > testvectors;
	stdvector<double> omega;

	vector_type evpU, evpB;
	bool m_bProjectedEVP;

#ifdef UG_PARALLEL
	IParallelCoarsening *m_pParallelCoarsening;
#endif
};


///	@}

} // namespace ug


#include "famg_impl.h"

#endif // __H__UG__LIB_DISC__AMG_SOLVER__FAMG_H__
