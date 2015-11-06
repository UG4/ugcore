/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef FAMG_ADDITIONAL_MATH_H_
#define FAMG_ADDITIONAL_MATH_H_
#include "lib_algebra/small_algebra/small_matrix/densematrix.h"
namespace ug{

inline void vecSum(double &erg, double alpha, double vec)
{
	erg = alpha * vec;
}

template<typename T>
inline void vecSum(typename T::value_type &erg, double alpha, const T &vec)
{
	erg = alpha * vec[0];
	for (size_t i = 1; i < vec.size(); i++)
		erg += vec[i];
	erg *= alpha;
}

inline double vecSum(double alpha, double vec)
{
	return alpha * vec;
}
template<typename T>
inline typename T::value_type vecSum(double alpha, const T &vec)
{
	typename T::value_type erg;
	vecSum(erg, alpha, vec);
	return erg;
}

inline void matSum(double &erg, double alpha, double vec)
{
	erg = alpha * vec;
}

template<typename T1, typename T2>
inline void matSum(T1 &erg, double alpha, T2 &mat)
{
	for (size_t r = 0; mat.num_rows(); r++)
	{
		erg[r] = mat(r, 0);
		for (size_t c = 1; c < mat.num_cols(); c++)
			erg[r] += mat(r, c);
		erg[r] *= alpha;
	}
}

template<typename T1, typename T2>
inline T1 matSum(double alpha, T2 &mat)
{
	T1 erg;
	matSum(erg, alpha, mat);
	return erg;
}

template<typename T1>
inline typename DenseMatrix<T1>::value_type Sum1Mat1(const DenseMatrix<T1> &mat)
{
	typename DenseMatrix<T1>::value_type ret = 0.0;
	for (size_t r = 0; mat.num_rows(); r++)
		for (size_t c = 1; c < mat.num_cols(); c++)
			ret += mat(r, c);
	return ret;
}

inline double matTrace(const double d)
{
	return d;
}

template<typename T1>
inline double matTrace(const DenseMatrix<T1> &mat)
{
	double tr=0.0;
	const size_t rk = (mat.num_rows() < mat.num_cols()) ? mat.num_rows() : mat.num_cols();
	for (size_t k = 0; k<rk; k++)
	{
		tr += mat(k, k);
	}
	return tr;
}

inline double matDiagMax(const double d)
{
	return d;
}

template<typename T1>
inline double matDiagMax(const DenseMatrix<T1> &mat)
{
	double val=0.0;
	double max=0.0;
	const size_t rk = (mat.num_rows() < mat.num_cols()) ? mat.num_rows() : mat.num_cols();
	for (size_t k = 0; k<rk; k++)
	{
		double abs=fabs(mat(k, k));
		if (abs>max) {
			val = mat(k, k);
			max = abs;
		}
	}
	return val;
}


inline double Sum1Mat1(double d)
{
	return d;
}

inline void GetDiag(double &a, double b)
{
	a = b;
}

template<typename T1, typename T2>
inline void GetDiag(T1 &m1, const T2 &m)
{
	UG_ASSERT(m.num_rows()==m.num_cols(), "");
	m1.resize(m.num_rows(), m.num_rows());
	m1=0.0;
	for (size_t i = 0; i < m.num_rows(); i++)
		m1(i, i) = m(i, i);
}

template<typename T1, typename T2>
inline void GetDiagSqrt(T1 &v, const T2 &m)
{
	UG_ASSERT(m.num_rows()==m.num_cols(), "");
	v.resize(m.num_rows());
	for (size_t i = 0; i < m.num_rows(); i++)
		v[i] = sqrt(m(i, i));
}

inline void GetDiagSqrt(double &a, double b)
{
	a = sqrt(b);
}

inline double EnergyProd(double v1, double M, double v2)
{
	return v1 * M * v2;
}

template<typename T1, typename T2>
inline double EnergyProd(const T1 &v1, const DenseMatrix<T2> &M, const T1 &v2)
{
	double sum = 0;
	for (size_t r = 0; r < M.num_rows(); r++)
	{
		double t = 0;
		for (size_t c = 0; c < M.num_cols(); c++)
			t += M(r, c) * v2[c];
		sum += t * v1[r];
	}
	return sum;
}

template<typename TMatrix>
void BlockMatrixToDoubleMatrix(DenseMatrix<VariableArray2<double> > &Ad, TMatrix &Ab)
{
	const size_t blockSize = block_traits<typename TMatrix::value_type>::static_num_rows;
	Ad.resize(blockSize*Ab.num_rows(), blockSize*Ab.num_cols());

	for(size_t i=0; i<Ab.num_rows(); i++)
	{
		for(size_t j=0; j<Ab.num_cols(); j++)
			for(size_t r=0; r<blockSize; r++)
				for(size_t c=0; c<blockSize; c++)
					Ad(i*blockSize+r, j*blockSize+c) = BlockRef(Ab(i, j), r, c);
	}
}

}
#endif /* ADDITIONAL_MATH_H_ */
