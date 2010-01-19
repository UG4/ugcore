/*
 *  matrixRow.h
 *  flexamg
 *
 *  Created by Martin Rupp on 18.01.10.
 *  Copyright 2010 . All rights reserved.
 *
 */

///////////////////////////////////////////////////////////////////

//!
//! class matrixrow
//! templated class, parameter is a blockmatrix type like double or blockDenseMatrix
//! you get a matrixrow variable when A is a SparseMatrix and you do A[i] or A.getrow(i)
//! you can do A[i].getDiag() or A[i].addMatrixRow(c, nr);
//! iterators: matrixrow<entry_type>::citerator it(A[i]), ++it, (*it).iIndex, dValue.
//! also possible: double a = A[i]*x. (thats why i like this concept)
//! its not 100% atm, because that const matrix_type *A would rather be a const matrix_type &A,
//! but this gives all problems with const and non-const objects (like in iterator and citerator). Ideas?
template<typename entry_type>
class matrixrow
{
	typedef SparseMatrix<entry_type> matrix_type;
	//typedef typename SparseMatrix<entry_type>::vec_type vec_type;
	typedef typename SparseMatrix<entry_type>::connection connection;
	//typedef Vector<vec_type> Vector_type;
	
public:
	typedef typename matrix_type::cRowIterator cRowIterator;
	typedef typename matrix_type::cLowerLeftIterator cLowerLeftIterator;
	typedef typename matrix_type::cUpperRightIterator cUpperRightIterator;
	
	// functions
public:
	matrixrow(const matrix_type &A_, const int row_) : row(row_), A(A_)
	{
		ASSERT2(row < A.rows, *this);
	}
	
	// get the sum of the entries of the SparseMatrix row
	//entry_type sum() const;
	
	
	//!
	//! operator []:
	//! get connection nr i.
	inline const connection &operator [] (int i) const;	
	//inline connection &operator [] (int i);
	
	inline int getConNr(int index) const;
	
	//!
	//! operator * with Vector. This Vector is templated, since then one can do
	//! blockDenseMatrix<3> = matrixrow<double> * Vector<blockDenseMatrix<3> >, or
	//! Vector<blockDenseMatrix<3> > = SparseMatrix<double> * Vector<blockDenseMatrix<3> >
	//! like the prolongation/restriction in MG.
	template<typename vec_type>
	inline vec_type operator * (const Vector<vec_type> &x) const;
	
	template<typename vec_type>
	inline void copyToMult(vec_type &d, const Vector<vec_type> &x) const;
	template<typename vec_type>
	inline void addToMult(vec_type &d, const Vector<vec_type> &x) const;
	template<typename vec_type>
	inline void substractFromMult(vec_type &d, const Vector<vec_type> &x) const;
	
	
	friend ostream &operator<<(ostream &output, const matrixrow &r)
	{
		output << "matrixrow[row = " << r.row << "] of " << r.A << ". ";
		return output;
	}
	
	// wrapped functions of SparseMatrix
	inline bool indexWithinBounds(int i) const	{return i >=0 && i < getNrOfConnections(); }
	inline int getRow() const { return row; }	
	inline int getNrOfConnections() const {	return A.getNrOfConnections(row);	}
	inline bool isUnconnected() const {	return A.isUnconnected(row); }
	
	void printtype() const { cout << *this; }	
	void print() const { A.printrow(row); }
	void p() const { A.printrow(row); } //gdb
	
	cRowIterator beginRow() const {	return A.beginRow(row); }
	cLowerLeftIterator beginLowerLeftRow(int row)  const { return cLowerLeftIterator(*this, row); }
	cUpperRightIterator beginUpperRightRow(int row)  const { return cUpperRightIterator(*this, row); }		
	//     data
	//----------------
	
private:
	const matrix_type &A;
	const int row;
};


#include "matrixRow.hpp"