/**
 * \file sparsematrix_print.h
 *
 * \author Martin Rupp
 *
 * \date 04.11.2009
 *
 * Goethe-Center for Scientific Computing 2009-2010.
 */

#ifndef __H__UG__CPU_ALGEBRA__SPARSEMATRIX_PRINT__
#define  __H__UG__CPU_ALGEBRA__SPARSEMATRIX_PRINT__

namespace ug {

//!
//! print to console whole SparseMatrix
template<typename T>
void SparseMatrix<T>::print(const char * const text) const
{
	cout << "================= SparseMatrix " << rows << "x" << cols << " =================" << endl;
	for(size_t i=0; i < rows; i++)
		getrow(i).print();
}


//!
//! print the row row to the console
template<typename T>
void SparseMatrix<T>::printrow(size_t row) const
{
#ifdef FLEXAMG
	cout << "row " << row << " [" << GetOriginalIndex(tolevel, row) << "] : ";
#else
	cout << "row " << row << ": ";
#endif
	for(cRowIterator it=beginRow(row); !it.isEnd(); ++it)
	{
		if((*it).dValue == 0.0) continue;
		cout << " ";
#ifdef FLEXAMG
		cout << "(" << (*it).iIndex << " [" << GetOriginalIndex(fromlevel, (*it).iIndex) << "] -> " << (*it).dValue << ")";
#else
		cout << "(" << (*it).iIndex << " -> " << (*it).dValue << ")";
#endif
	}
	cout << endl;
}

template<typename T>
void SparseMatrix<T>::printtype() const
{
	cout << *this;
}



#define CONNECTION_VIEWER_VERSION 1

// WriteMatrixToConnectionViewer
//--------------------------------------------------
/**
 * \brief writes to a file in somewhat SparseMatrix-market format (for connection viewer)
 * \param filename Filename to write matrix to
 * \param A SparseMatrix A.
 * \param positions Positions, there has to be one position for each i in (0, ..., max(A.num_rows(), A.num_cols())).
 */
template<typename Matrix_type, typename postype>
void WriteMatrixToConnectionViewer(const char *filename, const Matrix_type &A, postype *positions, int dimensions)
{
	fstream file(filename, ios::out);
	file << CONNECTION_VIEWER_VERSION << endl;
	file << dimensions << endl;

	int rows = A.num_rows();
	// write positions
	file << rows << endl;
	if(dimensions == 2)
		for(int i=0; i < rows; i++)
			file << positions[i][0] << " " << positions[i][1] << endl;
	else
		for(int i=0; i < rows; i++)
		  file << positions[i][0] << " " << positions[i][1] << " " << positions[i][2] << endl;

	file << 1 << endl; // show all cons
	// write connections
	for(int i=0; i < rows; i++)
	{
		for(typename Matrix_type::cRowIterator conn = A.beginRow(i); !conn.isEnd(); ++conn)
			if((*conn).dValue != 0.0)
				file << i << " " << (*conn).iIndex << " " << (*conn).dValue <<		endl;
	}
}

// WriteMatrixToConnectionViewer
//--------------------------------------------------
/**
 * \brief writes to a file in somewhat SparseMatrix-market format (for connection viewer)
 * \param filename Filename to write matrix to
 * \param A SparseMatrix A.
 * \param positions Positions, there has to be one position for each i in (0, ..., max(A.num_rows(), A.num_cols())).
 */
template<typename Vector_type, typename postype>
void WriteVectorToConnectionViewer(const char *filename, const Vector_type &b, postype *positions, int dimensions)
{
	fstream file(filename, ios::out);
	file << CONNECTION_VIEWER_VERSION << endl;
	file << dimensions << endl;

	int rows = b.size();
	// write positions
	file << rows << endl;
	if(dimensions == 2)
		for(int i=0; i < rows; i++)
			file << positions[i][0] << " " << positions[i][1] << endl;
	else
		for(int i=0; i < rows; i++)
		  file << positions[i][0] << " " << positions[i][1] << " " << positions[i][2] << endl;

	file << 1 << endl; // show all cons
	// write connections
	for(int i=0; i < rows; i++)
	{
		file << i << " " << i << " " << b[i] <<		endl;
	}
}

}
#endif // __H__UG__CPU_ALGEBRA__SPARSEMATRIX_PRINT__
