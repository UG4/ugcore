/*
 * sparsematrix_interface.h
 *
 *  Created on: 22.03.2011
 *      Author: mrupp
 */

#ifndef SPARSEMATRIX_INTERFACE_H_
#define SPARSEMATRIX_INTERFACE_H_

template<typename T> class SparseMatrix
{
public:
	// export value_type
	typedef T value_type;

public:
	SparseMatrix();
	virtual ~SparseMatrix();

	bool resize(size_t newRows, size_t newCols);

	// blas2-funcionts

	// element access functions
public:

	/** operator() (size_t r, size_t c) const
	 * access connection (r, c)
	 * \param r row
	 * \param c column
	 * \note it is assert'ed that connection (r,c) is there
	 * use operator()(r,c,bConnectionFound) to check.
	 * \return SparseMat(r, c)
	 */
	value_type &operator() (size_t r, size_t c);

	/** operator() (size_t r, size_t c) const
	 * access or create connection (r, c)
	 * \param r row
	 * \param c column
	 * \note (r,c) is added to sparsity pattern if not already there
	 * use operator()(r,c,bConnectionFound) to prevent
	 * \return SparseMat(r, c)=0.0 if connection created, otherwise SparseMat(r, c)
	 */
	const value_type &operator() (size_t r, size_t c) const;

	// information functions
public:
	inline size_t num_connections(size_t row) const;
	size_t num_rows() const { return rows; }
	size_t num_cols() const { return cols; }
	size_t total_num_connections() const { return iTotalNrOfConnections; }


	// iterators
	const_row_iterator begin_row(size_t row) const;
	row_iterator begin_row(size_t row);
	const_row_iterator end_row(size_t row) const;
	row_iterator end_row(size_t row);

	// connectivity functions
	//-------------------------
	const_row_iterator get_connection(size_t r, size_t c, bool &bFound) const;
	row_iterator get_connection(size_t r, size_t c, bool &bFound);

	const_row_iterator get_connection(size_t r, size_t c) const;
	row_iterator get_connection(size_t r, size_t c);

	// finalizing functions
	//----------------------
	void defragment();

}


// external functions which use the interface

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CreateAsMultiplyOf:
//-------------------------
/**
 * \brief Calculates M = A*B*C.
 * \param M (out) Matrix M, M = A*B*C$
 * \param A (in) Matrix A
 * \param B (in) Matrix B
 * \param C (in) Matrix C
 */
template<typename ABC_type, typename A_type, typename B_type, typename C_type>
void CreateAsProductOf(ABC_type &M, const A_type &A, const B_type &B, const C_type &C);


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CreateAsMultiplyOf:
//-------------------------
/**
 * \brief Calculates M = A*B.
 * \param M (out) Matrix M, M = A*B$
 * \param A (in) Matrix A
 * \param B (in) Matrix B
 */
template<typename AB_type, typename A_type, typename B_type>
void CreateAsMultiplyOf(AB_type &M, const A_type &A, const B_type &B);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// MatAdd:
//-------------------------
/**
 * \brief Calculates M = A + B
 * \param M (out) Matrix M, M = A + B
 * \param A (in) Matrix A
 * \param B (in) Matrix B
 * note: A and/or B may be equal to M.
 */
template<typename matrix_type>
void MatAdd(matrix_type &M, number &alpha1, const matrix_type &A, number &alpha2, const matrix_type &B);


// various SetDirichletRow and GetNeighborhood methods (see sparsematrix_util.h)


#endif /* SPARSEMATRIX_INTERFACE_H_ */
