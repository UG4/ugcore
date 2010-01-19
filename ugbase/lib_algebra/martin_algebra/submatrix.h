/*
 *  submatrix.h
 *  flexamg
 *
 *  Created by Martin Rupp on 18.01.10.
 *  Copyright 2010 . All rights reserved.
 *
 */

//!
//! submatrix
//! provide indices and nr of unknowns, and get submatrix of a bigger matrix
//! ATTENTION: indices and unknowns are NOT copied and have to be valid pointers
//! for the whole lifecycle of submatrix (performance reasons).
template <typename mat_type>
class submatrix
{
public:
	submatrix(int *indices_Rows, int *indices_Cols, int rows, int cols, bool zero = true)
	{
		nrRows = rows;
		nrCols = cols;
		values = new mat_type[nrRows*nrCols];
		
		if(zero)
		{
			for(int a=0; a<nrRows; a++) 
				for(int b=0; b<nrCols; b++)
					values[a*nrCols + b] = 0.0;
		}
		
		indRows = indices_Rows;
		indCols = indices_Cols;
		unknownsRows = unknownsCols = NULL;
	}
	
	submatrix(int *indices_Rows, int *indices_Cols, int *unknowns_Rows, int *unknowns_Cols, int rows, int cols, bool zero = true)
	{
		submatrix(indices_Rows, indices_Cols, rows, cols, zero);
		unknownsRows = unknowns_Rows;
		unknownsCols = unknowns_Cols;
	}

	submatrix(int *indices, int *unknowns_, int nrindices, bool zero = true)
	{ 	
		submatrix(indices, indices, nrindices, nrindices, zero)	;
		unknownsRows = unknowns_;
		unknownsCols = unknowns_;
	}

	
	submatrix(int *indices, int nrindices, bool zero = true) : submatrix(indices, indices, nrindices, nrindices, zero)
	{ 	}
	
	
	~submatrix()
	{
		delete[] values;
	}
	
	mat_type &operator () (int from, int to)
	{
		ASSERT1(from < nrRows && to >= 0 && to < nrCols && from >= 0);
		
		if(unknownsRows) setSize(values[from*nrCols + to], unknownsRows[from], unknownsCols[to]);
		return values[from*nrCols + to];
	}
	
	const mat_type &operator () (int from, int to) const 
	{
		ASSERT1(from < nrRows && to >= 0 && to < nrCols && from >= 0);
		return values[from*nrCols + to];
	}
	
	int getCols() const { return nrCols; }
	int getRows() const { return nrRows; }
	
	int getRowIndex(int localIndex) const { return indRows[localIndex]; }
	int getColIndex(int localIndex) const { return indCols[localIndex]; }
	
	
private:	
	int *indRows;
	int *indCols;

	int *unknownsRows;
	int *unknownsCols;
	int nrCols, nrRows;
	mat_type *values;
};


////////////////////////////////////////////////////


template <typename vec_type>
class subvector
{
public:
	subvector(int *indices_, int nr_)
	{
		nr = nr_;
		values = new vec_type[nr];
		
		for(int i=0; i<nr; i++) values[i] = 0.0;
		indices = indices_;
	}
	
	subvector(int *indices_, int *unknowns_, int nr_)
	{
		nr = nr_;
		values = new vec_type[nr];
		
		unknowns = unknowns_;
		indices = indices_;
	}
	subvector()
	{
		delete[] values;
	}
	
	vec_type &operator () (int i)
	{
		ASSERT1(i < nr && i >= 0);
		
		if(unknowns) setSize(values[i], unknowns[i]);
		return values[i];
	}
	
	int getNr() { return nr; }
	int getIndex(int localIndex) { return indices[localIndex]; }
	int getLocalFromIndex(int index)
	{
		for(int i=0; i<nr; i++)
			if(indices[i] == index) 
				return i;
		
		return -1;
	}
	
	
private:	
	int *unknowns;
	int *indices;
	int nr;
	vec_type *values;
};

