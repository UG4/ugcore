/**
 * \file amg_debug.h
 *
 * \author Martin Rupp
 *
 * \date 24.11.10
 *
 * Goethe-Center for Scientific Computing 2010.
 */

#include "common/common.h"
#ifdef UG_PARALLEL
#include "lib_algebra/parallelization/parallelization_util.h"
#endif

#ifndef __H__LIB_ALGEBRA__FAMG_SOLVER__FAMG_TESTVECTORS_H__
#define __H__LIB_ALGEBRA__FAMG_SOLVER__FAMG_TESTVECTORS_H__

namespace ug {

#if 0
static void get_testvector_constant(stdvector<double> &testvector, stdvector<size_t> &N, double &i_value, size_t i)
{
	for(size_t j=0; j<N.size(); j++)
		testvector[j] = 1.0;
	i_value = 1.0;
}
#endif


/*
template<typename value_type, typename pos_type>
void get_testvector_value_special(value_type &value, const pos_type &myPos, const pos_type &referencePos)
{
	BlockRef(value, TcomponentOut) = myPos[TcomponentIn] - referencePos[TcomponentIn];
}


inline void get_testvector_xx(DenseVector<VariableArray1<double> > &testvector, stdvector<size_t> &N2, int i_index)
{
	return get_testvector_special<0, 0>(testvector, N2, i_index);
}
inline void get_testvector_xy(DenseVector<VariableArray1<double> > &testvector, stdvector<size_t> &N2, int i_index)
{
	return get_testvector_special<0, 1>(testvector, N2, i_index);
}
...
2d :
(x 0), (y 0), (0 x), (0 y)
3d :
(x 0 0) (0 y 0) (0 0 z)  ..?

*/



template<typename matrix_type, typename vector_type>
double ScalProd(const vector_type &v1, const matrix_type &M, const vector_type &v2)
{
	double sum=0.0;
	UG_ASSERT(v1.size() == M.num_rows() && M.num_cols() == v2.size(), "size mismatch");
	for(size_t i=0; i < v1.size(); ++i)
	{
		for(typename matrix_type::const_row_iterator it = M.begin_row(i); it != M.end_row(i); ++it)
			sum += v1[i] * it.value() * v2[i];
	}
	return sum;
}

template<typename matrix_type, typename vector_type>
double EnergyProd(const vector_type &v1, const matrix_type &M)
{
	double sum=0.0;
	UG_ASSERT(v1.size() == M.num_rows() && M.num_cols() == v1.size(), "size mismatch");
	for(size_t i=0; i < v1.size(); ++i)
	{
		double s=0.0;
		for(typename matrix_type::const_row_iterator it = M.begin_row(i); it != M.end_row(i); ++it)
			s += v1[it.index()] * it.value();
		sum += v1[i]*s;
	}
	return sum;
}



template<typename matrix_type, typename vector_type>
void CalculateTestvector(const matrix_type &A_OL2, vector_type &big_testvector,
		size_t iTestvectorDamps)
{
	AMG_PROFILE_FUNC();
	vector_type d; d.resize(A_OL2.num_rows());
	for(size_t jj=0; jj < iTestvectorDamps; jj++)
	{
		MatMult(d, 1.0, A_OL2, big_testvector);
		for(size_t i=0; i<A_OL2.num_rows(); i++)
			big_testvector[i] = big_testvector[i] - 0.6*d[i]/A_OL2(i,i);
	}

	//big_testvector.set_storage_type(PST_CONSISTENT);
	//big_testvector.change_storage_type(PST_UNIQUE);
	//UG_LOG("Norm of testvector " << big_testvector.norm() << "\n");
	//big_testvector.change_storage_type(PST_CONSISTENT);
	//VecScaleAssign(big_testvector, 1/sqrt(EnergyProd(big_testvector, A_OL2)), big_testvector);
}


template<typename algebra_type>
void FAMGLevelCalculator<algebra_type>::calculate_testvectors()
{
	UG_SET_DEBUG_LEVEL(LIB_ALG_AMG, m_famg.iDebugLevelTestvectorCalc);

	// todo: all global?
	UG_DLOG(LIB_ALG_AMG, 1, "\ncalculating testvector... ");
	Stopwatch SW;
	if(bTiming) SW.start();

	vector_type d; d.resize(A.num_rows());
	vector_type c; c.resize(A.num_rows());
#ifdef UG_PARALLEL
	c.set_layouts(A.master_layout(), A.slave_layout());
	d.set_layouts(A.master_layout(), A.slave_layout());
#endif

	for(size_t i=0; i<m_testvectors.size(); i++)
	{

		UG_ASSERT(m_testvectors[i].size() == A.num_rows(), m_testvectors[i].size() << " != " << A.num_rows());

#ifdef UG_PARALLEL
		m_testvectors[i].set_storage_type(PST_CONSISTENT);
#endif

		m_famg.m_testvectorsmoother->init(m_famg.levels[level]->pAgglomeratedA);

		if(m_famg.m_writeMatrices && m_famg.m_writeTestvectors)
					for(size_t i=0; i<m_testvectors.size(); i++)
						WriteVectorToConnectionViewer((m_famg.m_writeMatrixPath + ToString("testvector_pS") + ToString(i) + ToString("_L") + ToString(level) + ".vec").c_str(),
							m_testvectors[i], &m_famg.m_amghelper.positions[level][0], m_famg.m_dbgDimension);

		for(size_t jj=0; jj < m_famg.get_testvector_smooths(); jj++)
		{
			d.set(0.0);
			c.set(0.0);
			A.apply(d, m_testvectors[i]);
			m_famg.m_testvectorsmoother->apply(c, d);
			m_testvectors[i] -= c;
		}

		if(m_famg.m_writeMatrices && m_famg.m_writeTestvectors)
			for(size_t i=0; i<m_testvectors.size(); i++)
				WriteVectorToConnectionViewer((m_famg.m_writeMatrixPath + ToString("testvector_S") + ToString(i) + ToString("_L") + ToString(level) + ".vec").c_str(),
					m_testvectors[i], &m_famg.m_amghelper.positions[level][0], m_famg.m_dbgDimension);
		//VecScaleAdd(m_testvectors[i], 1.0, m_testvectors[i], -1.0, d);
	}


	if(bTiming) UG_DLOG(LIB_ALG_AMG, 1, "took " << SW.ms() << " ms");
}
template<typename algebra_type>
void FAMGLevelCalculator<algebra_type>::get_testvectors_on_OL2()
{
#ifdef UG_PARALLEL
	for(size_t i=0; i<m_testvectors.size(); i++)
	{
		SetLayoutValues(&m_testvectors[i], A.slave_layout(), 0.0);
		m_testvectors[i].resize(A_OL2.num_rows());
		m_testvectors[i].set_master_layout(A_OL2.master_layout());
		m_testvectors[i].set_slave_layout(A_OL2.slave_layout());
		for(size_t j=A.num_rows(); j<A_OL2.num_rows(); j++)
			m_testvectors[i][j]=0.0;

		m_testvectors[i].set_storage_type(PST_ADDITIVE);
		m_testvectors[i].change_storage_type(PST_CONSISTENT);
	}
#endif
}


template<typename algebra_type>
void FAMGLevelCalculator<algebra_type>::calculate_next_testvectors()
{
	vector_type t;
	t.resize(AH.num_rows());
	// todo: remove dynamic cast, change big_testvector to parallel
	for(size_t i=0; i<m_testvectors.size(); i++)
	{
		m_famg.injection(t, m_testvectors[i], level);
		m_testvectors[i] = t;
	}
}

}
#endif // __H__LIB_ALGEBRA__FAMG_SOLVER__FAMG_TESTVECTORS_H__
