/*
 *  submatrix.h
 *  flexamg
 *
 *  Created by Martin Rupp on 18.01.10.
 *  Copyright 2010 . All rights reserved.
 *
 */
#include <iomanip> 
//!
//! submatrix
//! provide indices and nr of unknowns, and get submatrix of a bigger matrix
//! ATTENTION: indices and unknowns are NOT copied and have to be valid pointers
//! for the whole lifecycle of submatrix (performance reasons) (ideas?)
//! nr unknowns have to be known in some cases to init variable blockmatrices
//! (ideas?)
template <typename entry_type>
class submatrix
{
public:
	//! create the submatrix
	//! @param indices_Rows pointer to array of integers of indices of the rows of the submatrix
	//! @param indices_Cols pointer to array of integers of indices of the columns of the submatrix
	//! @param rows number of rows
	//! @param cosl number of columns
	//! @param zero if true, matrix is initialized with 0.0.	
	void create(int *indices_Rows, int *indices_Cols, int rows, int cols, bool zero = true)
	{
		nrRows = rows;
		nrCols = cols;
		values = new entry_type[nrRows*nrCols];
		
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
	
public:
	//! constructor for different rows and columns
	submatrix(int *indices_Rows, int *indices_Cols, int rows, int cols, bool zero = true)
	{
		create(indices_Rows, indices_Cols, rows, cols, zero);
	}
	
	//! constructor for different rows and columns, and nr of unknowns specified
	submatrix(int *indices_Rows, int *indices_Cols, int *unknowns_Rows, int *unknowns_Cols, int rows, int cols, bool zero = true)
	{
		create(indices_Rows, indices_Cols, rows, cols, zero);
		unknownsRows = unknowns_Rows;
		unknownsCols = unknowns_Cols;
	}

	//! constructor for rows=columns and unknowns specified
	submatrix(int *indices, int *unknowns_, int nrindices, bool zero = true)
	{ 	
		create(indices, indices, nrindices, nrindices, zero);
		unknownsRows = unknowns_;
		unknownsCols = unknowns_;
	}

	//! constructor for rows=columns
	submatrix(int *indices, int nrindices, bool zero = true)
	{
		create(indices, indices, nrindices, nrindices, zero);
	}
	
	//! empty constructor
	submatrix()
	{
		indRows = indCols = NULL;
		unknownsCols = unknownsRows = NULL;
		nrCols = nrRows = 0;
		values = NULL;
	}
	
	~submatrix()
	{
		delete[] values;
	}
	

	
	//! access element with LOCAL indices
	entry_type &operator () (int from, int to)
	{
		ASSERT1(from < nrRows && to >= 0 && to < nrCols && from >= 0);
		
		if(unknownsRows) setSize(values[from*nrCols + to], unknownsRows[from], unknownsCols[to]);
		return values[from*nrCols + to];
	}
	
	//! access element with LOCAL indices (const)
	const entry_type &operator () (int from, int to) const 
	{
		ASSERT1(from < nrRows && to >= 0 && to < nrCols && from >= 0);
		return values[from*nrCols + to];
	}
	
	int getCols() const { return nrCols; }
	int getRows() const { return nrRows; }
	
	int getRowIndex(int localIndex) const { return indRows[localIndex]; }
	int getColIndex(int localIndex) const { return indCols[localIndex]; }
	
	int getLocalRowIndex(int globalIndex) const 
	{ 
		for(int i=0; i<nrRows; i++)
			if(indRows[i] == globalIndex) return i;
		return -1;
	}
	int getLocalColIndex(int globalIndex) const 
	{ 
		for(int i=0; i<nrCols; i++)
			if(indCols[i] == globalIndex) return i;
		return -1;
	}
	
	
	//ostream &operator << (ostream out, const submatrix &SM)
	void print()
	{
		int w = 12;
		cout << setw(w) << "|";		
		for(int i=0; i<nrCols; i++)
			cout << setw(w-1) << indCols[i] << "|";
		cout << endl;
		for(int i=0; i<nrCols+1; i++)
			cout << setw(w) << setfill('-') << "|";
		cout << endl;
		for(int i=0; i<nrRows; i++)
		{
			cout << setfill(i % 2 ? '.' : ' ') << setw(w-1) << indRows[i] << "|";
			for(int j=0; j < nrCols; j++)
				cout << right << setw(w-1) << values[i*nrCols + j] << "|";
			cout << endl;
		}
//		return out;
	}
	
	
	void amg_inv()
	{
		static vector<__CLPK_integer> interchange(40);
		interchange.resize(max(getCols(), getRows()));
		__CLPK_integer info = 0;
		__CLPK_integer rows = getRows();
		__CLPK_integer cols = getCols();
		
		dgetrf_(&rows, &cols, values, &rows, &interchange[0], &info);
		ASSERT2(info == 0, "info is " << info << ( info > 0 ? ": SparseMatrix singular in U(i,i)" : ": i-th argument had had illegal value"));
		
		/*double worksize; __CLPK_integer iWorksize = -1;
		dgetri_(&rows, &SM(0,0), &rows, &interchange[0], &worksize, &iWorksize, &info);
		ASSERT1(info == 0);
		iWorksize = worksize;
		double *work = new double[iWorksize];*/
		
		__CLPK_integer iWorksize = 1024;
		vector<double> work(1024);
		dgetri_(&rows, values, &rows, &interchange[0], &work[0], &iWorksize, &info);
		ASSERT1(info == 0);
	}
	
private:	
	int *indRows;		///< stores indices of the rows
	int *indCols;		///< stores indices of the columns

	int *unknownsRows;	///< stores nr of unknowns of the rows
	int *unknownsCols;
	int nrCols, nrRows;	///< stores nr of unknowns of the rows
	entry_type *values;
	bool copiedIndices;
};


////////////////////////////////////////////////////

//!
//! subvector
//! provide indices and nr of unknowns, and get subvector of a bigger vector
//! ATTENTION: indices and unknowns are NOT copied and have to be valid pointers
//! for the whole lifecycle of subvector (performance reasons).
template <typename vec_type>
class subvector
{
public:
	void create(int *indices_, int nr_)
	{
		nr = nr_;
		values = new vec_type[nr];
		
		for(int i=0; i<nr; i++) values[i] = 0.0;
		indices = indices_;
	}
	
	subvector()
	{
		nr = 0; 
		values = NULL;
		indices = NULL;
	}
	
	subvector(int *indices_, int nr_)
	{
		create(indices_, nr);
	}
		
	
	subvector(int *indices_, int *unknowns_, int nr_)
	{
		create(indices_, nr);		
		unknowns = unknowns_;
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

