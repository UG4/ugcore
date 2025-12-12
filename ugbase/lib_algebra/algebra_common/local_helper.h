/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
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

/*
 *	Diese Funktionen sollen das Benutzen von SparseMatrix::add(M&m) erleichtern.
 */
#ifndef __H__UG__MARTIN_ALGEBRA__LOCAL_HELPER__
#define __H__UG__MARTIN_ALGEBRA__LOCAL_HELPER__

namespace ug {


template<typename M>
class localMatrix_from_mat_and_array
{
public:
	localMatrix_from_mat_and_array(M &m_, const size_t *rows_, const size_t *cols_) : rows(rows_), cols(cols_)
	{
		m = &m_;
	}

	~localMatrix_from_mat_and_array() = default;

	[[nodiscard]] size_t num_rows() const { return m->num_rows(); }
	[[nodiscard]] size_t num_cols() const { return m->num_cols(); }
	[[nodiscard]] size_t row_index(size_t i) const { return rows[i]; }
	[[nodiscard]] size_t col_index(size_t i) const { return cols[i]; }
	typename M::value_type &operator () (size_t i, size_t j) { return (*m)(i,j); }
	const typename M::value_type &operator () (size_t i, size_t j) const { return (*m)(i,j); }

private:
	M *m;
	const size_t *rows;
	const size_t *cols;
};

template<typename M>
class const_localMatrix_from_mat_and_array
{
public:
	const_localMatrix_from_mat_and_array(const M &m_, const size_t *rows_, const size_t *cols_) : m(m_), rows(rows_), cols(cols_)
	{	}

	[[nodiscard]] size_t num_rows() const { return m.num_rows(); }
	[[nodiscard]] size_t num_cols() const { return m.num_cols(); }
	[[nodiscard]] size_t row_index(size_t i) const { return rows[i]; }
	[[nodiscard]] size_t col_index(size_t i) const { return cols[i]; }
	const typename M::value_type &operator () (size_t i, size_t j) const { return m(i,j); }

private:
	const M &m;
	const size_t *rows;
	const size_t *cols;
};

template<typename T>
class localMatrix_from_col_major_and_array
{
public:
	localMatrix_from_col_major_and_array(size_t numrows_, size_t numcols_, T *m_, const size_t *rows_, const size_t *cols_)
			: m(m_), rows(rows_), cols(cols_)
	{
		numrows = numrows_;
		numcols = numcols_;
	}

	[[nodiscard]] size_t num_rows() const { return numrows; }
	[[nodiscard]] size_t num_cols() const { return numcols; }
	[[nodiscard]] size_t row_index(size_t i) const { return rows[i]; }
	[[nodiscard]] size_t col_index(size_t i) const { return cols[i]; }
	T &operator () (size_t i, size_t j) { return m[i + j*numcols]; }
	const T &operator () (size_t i, size_t j) const { return m[i + j*numcols]; }

private:
	T *m;
	size_t numrows;
	size_t numcols;
	const size_t *rows;
	const size_t *cols;
};

template<typename T>
class localMatrix_from_row_major_and_array
{
public:
	localMatrix_from_row_major_and_array(size_t numrows_, size_t numcols_, T *m_, const size_t *rows_, const size_t *cols_)
			: m(m_), rows(rows_), cols(cols_)
	{
		numrows = numrows_;
		numcols = numcols_;
	}

	[[nodiscard]] size_t num_rows() const { return numrows; }
	[[nodiscard]] size_t num_cols() const { return numcols; }
	[[nodiscard]] size_t row_index(size_t i) const { return rows[i]; }
	[[nodiscard]] size_t col_index(size_t i) const { return cols[i]; }
	T &operator () (size_t i, size_t j) { return m[i*numrows + j]; }
	const T &operator () (size_t i, size_t j) const { return m[i*numrows + j]; }

private:
	T *m;
	const size_t *rows;
	const size_t *cols;
	size_t numrows;
	size_t numcols;
};

template<typename T>
class localVector_from_array
{
public:
	localVector_from_array(size_t N_, T *v_, const size_t *indices_)
			: N(N_), v(v_), indices(indices_)
	{
	}

	[[nodiscard]] size_t size() const { return N; }
	[[nodiscard]] size_t index(size_t i) const { return indices[i]; }

	T &operator [] (size_t i) { return v[i]; }
	const T &operator [] (size_t i) const { return v[i]; }

private:
	size_t N;
	T *v;
	const size_t *indices;
};



////////////////////////////////////////////////////////////////////////////////////////////////////


/** Add a local matrix
 *
 * The local matrix type must declare the following members:
 * - num_rows()
 * - num_cols()
 * - row_index(size_t i)
 * - col_index(size_t j)
 * - operator () (size_t i, size_t j)
 * so that mat(i,j) will go to SparseMat(mat.row_index(i), mat.col_index(j))
 * \param mat the whole local matrix type
 */
template<typename TGlobalMatrix, typename TLocalMatrix>
inline bool AddLocalMatrix(TGlobalMatrix &mat, const TLocalMatrix &localMat)
{
	mat.add(localMat);
	return true;
}
template<typename TGlobalMatrix, typename TLocalMatrix>
inline bool SetLocalMatrix(TGlobalMatrix &mat, const TLocalMatrix &localMat)
{
	mat.set(localMat);
	return true;
}
template<typename TGlobalMatrix, typename TLocalMatrix>
inline bool GetLocalMatrix(const TGlobalMatrix &mat, TLocalMatrix &localMat)
{
	mat.get(localMat);
	return true;
}

template<typename TGlobalMatrix, typename TLocalMatrix>
inline bool AddLocalMatrix(TGlobalMatrix &mat, const TLocalMatrix &localMat, const size_t *rowIndices, const size_t *colIndices)
{
	const_localMatrix_from_mat_and_array<TLocalMatrix> loc(localMat, rowIndices, colIndices);
	return AddLocalMatrix(mat, loc);
}

template<typename TGlobalMatrix, typename TLocalMatrix>
inline bool SetLocalMatrix(TGlobalMatrix &mat, TLocalMatrix &localMat, const size_t *rowIndices, const size_t *colIndices)
{
	const_localMatrix_from_mat_and_array<TLocalMatrix> loc(localMat, rowIndices, colIndices);
	return AddLocalMatrix(mat, loc);
}

template<typename TGlobalMatrix, typename TLocalMatrix>
inline bool GetLocalMatrix(const TGlobalMatrix &mat, TLocalMatrix &localMat, const size_t *rowIndices, const size_t *colIndices)
{
	localMatrix_from_mat_and_array<TLocalMatrix> loc(localMat, rowIndices, colIndices);
	return GetLocalMatrix(mat, loc);
}


}

#endif
