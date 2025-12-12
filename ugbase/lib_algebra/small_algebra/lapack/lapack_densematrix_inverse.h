/*
 * Copyright (c) 2010-2014:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__SMALL_ALGEBRA__DENSEMATRIX_INVERT_H__
#define __H__UG__SMALL_ALGEBRA__DENSEMATRIX_INVERT_H__

#include "../small_matrix/densematrix.h"
#include "../small_matrix/densevector.h"
#include "../../common/operations.h"

namespace ug {
///////////////////////////////////////////////////////////////////////////////////////
/**
 * smallInverse<size_t n>
 * A class to hold an inverse of a smallMatrix<n>
 * implemented with LAPACKs LU-Decomposition dgetrf
 * (uses double[n*n] for LU and interchange[n] for pivoting
 */
template<typename TStorage>
class DenseMatrixInverse
{
private: // storage
	DenseMatrix<TStorage> densemat;
	std::vector<lapack_int> interchange;

public:
	inline size_t num_cols() const
	{
		return densemat.num_cols();
	}

	inline size_t num_rows() const
	{
		return densemat.num_rows();
	}

	inline void resize(size_t k)
	{
		densemat.resize(k,k);
		densemat = 0.0;
	}
	inline void resize(size_t k, size_t k2)
	{
		UG_COND_THROW(k!=k2, "only square matrices supported");
		resize(k);
	}

	double &operator () (int r, int c)
	{
		return densemat(r,c);
	}
	const double &operator () (int r, int c) const
	{
		return densemat(r,c);
	}

public:
	//! initializes this object as inverse of mat
	bool set_as_inverse_of(const DenseMatrix<TStorage> &mat)
	{
		UG_COND_THROW(mat.num_rows() != mat.num_cols(), "only for square matrices");

		densemat = mat;
		return invert();
	}

	bool invert()
	{
		if(densemat.num_rows() == 0) return false;

		interchange.resize(densemat.num_rows());
		int info = getrf(densemat.num_rows(), densemat.num_cols(), &densemat(0,0),
				densemat.num_rows(), &interchange[0]);
		if(info != 0)
		{
			if(info > 0)
			{UG_THROW("ERROR in 'DenseMatrixInverse::invert': Matrix singular in U(i,i), with i="<<info<<"\n");}
			if(info < 0)
			{UG_THROW("ERROR in 'DenseMatrixInverse::invert':  i-th argument had had illegal value, with i="<<info<<"\n");}
		}
		return info == 0;
	}

	template<typename vector_t>
	void apply(DenseVector<vector_t> &vec) const
	{
		if(vec.size() == 0) return;
		int info = getrs(ModeNoTrans, num_rows(), 1, &densemat(0,0), num_rows(), &interchange[0], &vec[0], num_rows());
		(void) info;
		UG_COND_THROW(info != 0, "DenseMatrixInverse::mat_mult: getrs failed.");
	}

	// todo: implement operator *=

	template<typename T> friend std::ostream &operator << (std::ostream &out, const DenseMatrixInverse<T> &mat);
};


template<typename T>
std::ostream &operator << (std::ostream &out, const DenseMatrixInverse<T> &mat)
{
	out << "[ ";
	for(size_t r=0; r<mat.num_rows(); ++r)
	{
		for(size_t c=0; c<mat.num_cols(); ++c)
			out << mat.densemat(r, c) << " ";
		if(r != mat.num_rows()-1) out << "| ";
	}
	out << "]";
	out << " (DenseMatrixInverse " << mat.num_rows() << "x" << mat.num_cols() << ", " << ((T::ordering == ColMajor) ? "ColMajor)" : "RowMajor)");

	return out;
}


template<typename T>
struct matrix_algebra_type_traits<DenseMatrixInverse<T> >
{
	static constexpr int type = MATRIX_USE_GLOBAL_FUNCTIONS;
};

}
#endif
