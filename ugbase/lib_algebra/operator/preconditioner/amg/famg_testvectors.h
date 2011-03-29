/**
 * \file amg_debug.h
 *
 * \author Martin Rupp
 *
 * \date 24.11.10
 *
 * Goethe-Center for Scientific Computing 2010.
 */


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
		for(typename matrix_type::const_row_iterator it = M.begin_row(i); it != M.end_row(i); ++it)
			sum += (v1[i]*v1[i]) * it.value();
	}
	return sum;
}



template<typename matrix_type, typename vector_type>
void CalculateTestvector(const matrix_type &A_OL2, vector_type &big_testvector,
		size_t iTestvectorDamps)
{
	FAMG_PROFILE_FUNC();
	vector_type d; d.resize(A_OL2.num_rows());
	for(size_t jj=0; jj < iTestvectorDamps; jj++)
	{
		MatMult(d, 1.0, A_OL2, big_testvector);
		for(size_t i=0; i<A_OL2.num_rows(); i++)
			big_testvector[i] = big_testvector[i] - 0.6*d[i]/A_OL2(i,i);
	}

	VecScaleAssign(big_testvector, 1/sqrt(EnergyProd(big_testvector, A_OL2)), big_testvector);
}

template<typename matrix_type, typename vector_type>
void CalculateNextTestvector(const matrix_type &R, vector_type &big_testvector)
{
	FAMG_PROFILE_FUNC();
	vector_type t;
	t.resize(R.num_rows());

	MatMult(t, 1.0, R, big_testvector);
	big_testvector.resize(R.num_rows());
	big_testvector = t;
}

}

#endif // __H__LIB_ALGEBRA__FAMG_SOLVER__FAMG_TESTVECTORS_H__
