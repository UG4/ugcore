/*
 *  local_helper.h
 *
 *  Created by Martin Rupp on 27.7.2010.
 *  Copyright 2010 G-CSC, University of Frankfurt. All rights reserved.
 *
 *	Diese Funktionen sollen das Benutzen von SparseMatrix::add(M&m) erleichtern.
 */
#ifndef __H__UG__MARTIN_ALGEBRA__LOCAL_HELPER__
#define __H__UG__MARTIN_ALGEBRA__LOCAL_HELPER__

template<typename M>
class localMatrix_from_mat_and_array
{
public:
	localMatrix_from_mat_and_array(M &m_, size_t *rows_, size_t *cols_) : rows(rows_), cols(cols_)
	{
		m = &m_;
	}

	~localMatrix_from_mat_and_array()
	{

	}

	size_t num_rows() const { return m->num_rows(); }
	size_t num_cols() const { return m->num_cols(); }
	size_t row_index(size_t i) const { return rows[i]; }
	size_t col_index(size_t i) const { return cols[i]; }
	typename M::value_type &operator()(size_t i, size_t j) { return (*m)(i,j); }
	const typename M::value_type &operator()(size_t i, size_t j) const { return (*m)(i,j); }

private:
	M *m;
	size_t *rows;
	size_t *cols;
};

template<typename M>
localMatrix_from_mat_and_array<M> LocalMatrix(M &m, size_t *rows, size_t *cols)
{
	return localMatrix_from_mat_and_array<M> (m, rows, cols);
}

template<typename M>
class const_localMatrix_from_mat_and_array
{
public:
	const_localMatrix_from_mat_and_array(const M &m_, size_t *rows_, size_t *cols_) : m(m_), rows(rows_), cols(cols_)
	{	}

	size_t num_rows() const { return m.num_rows(); }
	size_t num_cols() const { return m.num_cols(); }
	size_t row_index(size_t i) const { return rows[i]; }
	size_t col_index(size_t i) const { return cols[i]; }
	const typename M::value_type &operator()(size_t i, size_t j) const { return m(i,j); }

private:
	const M &m;
	size_t *rows;
	size_t *cols;
};

template<typename M>
const_localMatrix_from_mat_and_array<M> LocalMatrix(const M &m, size_t *rows, size_t *cols)
{
	return const_localMatrix_from_mat_and_array<M> (m, rows, cols);
}

template<typename T>
class localMatrix_from_col_major_and_array
{
public:
	localMatrix_from_col_major_and_array(size_t numrows_, size_t numcols_, T *m_, size_t *rows_, size_t *cols_)
			: m(m_), rows(rows_), cols(cols_)
	{
		numrows = numrows_;
		numcols = numcols_;
	}

	size_t num_rows() const { return numrows; }
	size_t num_cols() const { return numcols; }
	size_t row_index(size_t i) const { return rows[i]; }
	size_t col_index(size_t i) const { return cols[i]; }
	T &operator()(size_t i, size_t j) { return m[i + j*numcols]; }
	const T &operator()(size_t i, size_t j) const { return m[i + j*numcols]; }

private:
	T *m;
	size_t numrows;
	size_t numcols;
	size_t *rows;
	size_t *cols;
};

template<typename T>
class localMatrix_from_row_major_and_array
{
public:
	localMatrix_from_row_major_and_array(size_t numrows_, size_t numcols_, T *m_, size_t *rows_, size_t *cols_)
			: m(m_), rows(rows_), cols(cols_)
	{
		numrows = numrows_;
		numcols = numcols_;
	}

	size_t num_rows() const { return numrows; }
	size_t num_cols() const { return numcols; }
	size_t row_index(size_t i) const { return rows[i]; }
	size_t col_index(size_t i) const { return cols[i]; }
	T &operator()(size_t i, size_t j) { return m[i*numrows + j]; }
	const T &operator()(size_t i, size_t j) const { return m[i*numrows + j]; }

private:
	T *m;
	size_t *rows;
	size_t *cols;
	size_t numrows;
	size_t numcols;
};

template<typename T>
class localVector_from_array
{
public:
	localVector_from_array(size_t N_, T *v_, size_t *indices_)
			: N(N_), v(v_), indices(indices_)
	{
	}

	size_t size() const { return N; }
	size_t index(size_t i) const { return indices[i]; }

	T &operator[](size_t i) { return v[i]; }
	const T &operator[](size_t i) const { return v[i]; }

private:
	size_t N;
	T *v;
	size_t *indices;
};


#endif
