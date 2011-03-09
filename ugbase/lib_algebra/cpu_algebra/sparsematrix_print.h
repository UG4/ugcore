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

#include "sparsematrix.h"

namespace ug {

//!
//! print to console whole SparseMatrix
template<typename T>
void SparseMatrix<T>::print(const char * const text) const
{
	UG_LOG("================= SparseMatrix " << rows << "x" << cols << " =================\n");
	for(size_t i=0; i < rows; i++)
		printrow(i);
}


//!
//! print the row row to the console
template<typename T>
void SparseMatrix<T>::printrow(size_t row) const
{
	UG_LOG("row " << row << ": ");
	for(const_row_iterator it=begin_row(row); it != end_row(row); ++it)
	{
		if(it.value() == 0.0) continue;
		UG_LOG(" ");
		UG_LOG("(" << it.index() << " -> " << it.value() << ")");
	}

	UG_LOG("\n");
}

template<typename T>
void SparseMatrix<T>::printtype() const
{
	std::cout << *this;
}



#define CONNECTION_VIEWER_VERSION 1

// WriteMatrixToConnectionViewer
//--------------------------------------------------
/**
 * \brief writes to a file in somewhat SparseMatrix-market format (for connection viewer)
 * \param filename Filename to write matrix to
 * \param A SparseMatrix A.
 * \param positions Positions, there has to be one position for each i in (0, ..., max(A.num_rows(), A.num_cols())).
 * \param dimensions Dimension of positions
 */
template<typename Matrix_type, typename postype>
void WriteMatrixToConnectionViewer(const char *filename, const Matrix_type &A, postype *positions, int dimensions)
{
	std::fstream file(filename, std::ios::out);
	file << CONNECTION_VIEWER_VERSION << std::endl;
	file << dimensions << std::endl;

	size_t rows = A.num_rows();

	// write positions
	file << rows << std::endl;
	if(dimensions == 1)
			for(size_t i=0; i < rows; i++)
				file << positions[i][0] << " 0.0"  << std::endl;
	else if(dimensions == 2)
		for(size_t i=0; i < rows; i++)
			file << positions[i][0] << " " << positions[i][1] << std::endl;
	else
		for(size_t i=0; i < rows; i++)
		  file << positions[i][0] << " " << positions[i][1] << " " << positions[i][2] << std::endl;

	file << 1 << std::endl; // show all cons
	// write connections
	for(size_t i=0; i < rows; i++)
	{
		for(typename Matrix_type::const_row_iterator conn = A.begin_row(i); conn != A.end_row(i); ++conn)
			if(conn.value() != 0.0)
				file << i << " " << conn.index() << " " << conn.value() <<		std::endl;
	}
}

// WriteMatrixToConnectionViewer
//--------------------------------------------------
/**
 * this version can handle different from and to spaces
 */
template <typename Matrix_type, typename postype>
bool WriteMatrixToConnectionViewer(	const char *filename,
									const Matrix_type &A,
									std::vector<postype> &positionsFrom, std::vector<postype> &positionsTo, size_t dimensions)
{
	const char * p = strstr(filename, ".mat");
	if(p == NULL)
	{
		UG_LOG("Currently only '.mat' format supported for domains.\n");
		return false;
	}

	if(positionsFrom.size() != A.num_cols())
	{
		UG_LOG("uFrom.size() != A.num_cols() !\n");
		return false;
	}
	if(positionsTo.size() != A.num_rows())
	{
		UG_LOG("uTo.size() != A.num_rows() !\n");
		return false;
	}

	std::vector<postype> positions;
	std::vector<size_t> mapFrom, mapTo;
	mapFrom.resize(positionsFrom.size());
	mapTo.resize(positionsTo.size());

	if(positionsFrom.size() > positionsTo.size())
	{
		positions.resize(positionsFrom.size());
		for(size_t i=0; i<positionsFrom.size(); i++)
		{
			positions[i] = positionsFrom[i];
			mapFrom[i] = i;
		}


		for(size_t i=0; i<positionsTo.size(); i++)
		{
			size_t j;
			for(j=0; j<positionsFrom.size(); j++)
			{
				if(positionsTo[i] == positionsFrom[j])
					break;
			}
			mapTo[i] = j;
			if(j == positionsFrom.size())
				positions.push_back(positionsTo[i]);
		}
	}
	else
	{
		positions.resize(positionsTo.size());
		for(size_t i=0; i<positionsTo.size(); i++)
		{
			positions[i] = positionsTo[i];
			mapTo[i] = i;
		}

		for(size_t  i=0; i<positionsFrom.size(); i++)
		{
			size_t j;
			for(j=0; j<positionsTo.size(); j++)
			{
				if(positionsFrom[i] == positionsTo[j])
					break;
			}
			mapFrom[i] = j;
			if(j == positionsTo.size())
				positions.push_back(positionsFrom[i]);
		}
	}


	std::fstream file(filename, std::ios::out);
	file << CONNECTION_VIEWER_VERSION << std::endl;
	file << dimensions << std::endl;

	// write positions
	file << positions.size() << std::endl;

	if(dimensions == 1)
		for(size_t i=0; i < positions.size(); i++)
			file << positions[i][0] << " 0.0"  << std::endl;
	else if(dimensions == 2)

		for(size_t i=0; i < positions.size(); i++)
			file << positions[i][0] << " " << positions[i][1] << std::endl;
	else
		for(size_t i=0; i < positions.size(); i++)
		  file << positions[i][0] << " " << positions[i][1] << " " << positions[i][2] << std::endl;

	file << 1 << std::endl; // show all cons
	// write connections
	for(size_t i=0; i < A.num_rows(); i++)
	{
		for(typename Matrix_type::const_row_iterator conn = A.begin_row(i); conn != A.end_row(i); ++conn)
			if(conn.value() != 0.0)
				file << mapTo[i] << " " << mapFrom[conn.index()] << " " << conn.value() <<		std::endl;
	}
	return true;
}


// WriteMatrixToConnectionViewer
//--------------------------------------------------
/**
 * \brief writes to a file in somewhat SparseMatrix-market format (for connection viewer)
 * \param filename Filename to write matrix to
 * \param b Vector
 * \param positions Positions, there has to be one position for each i in (0, ..., max(A.num_rows(), A.num_cols())).
 * \param dimensions	Dimensions of Positions
 */
template<typename Vector_type, typename postype>
void WriteVectorToConnectionViewer(const char *filename, const Vector_type &b, postype *positions, int dimensions)
{
	std::fstream file(filename, std::ios::out);
	file << CONNECTION_VIEWER_VERSION << std::endl;
	file << dimensions << std::endl;

	size_t rows = b.size();
	// write positions
	file << rows << std::endl;

	if(dimensions == 1)
		for(size_t i=0; i < rows; i++)
			file << positions[i][0] << " 0.0"  << std::endl;
	else if(dimensions == 2)
		for(size_t i=0; i < rows; i++)
			file << positions[i][0] << " " << positions[i][1] << std::endl;
	else
		for(size_t i=0; i < rows; i++)
		  file << positions[i][0] << " " << positions[i][1] << " " << positions[i][2] << std::endl;

	file << 1 << std::endl; // show all cons
	// write connections
	for(size_t i=0; i < rows; i++)
	{
		file << i << " " << i << " " << b[i] <<		std::endl;
	}
}


#if 0
template<typename Vector_type, typename postype>
void WriteVectorToConnectionViewerNEW(const char *filename, const Vector_type &b, postype *positions, int dimensions)
{
	std::fstream file(filename, std::ios::out);
	file << CONNECTION_VIEWER_VERSION << std::endl;
	file << 3 << std::endl;

	double nmax=0;
	for(size_t i=0; i<b.size(); i++)
		if(nmax < BlockNorm(b[i])) nmax = BlockNorm(b[i]);

	nmax*=4;
	double scale = 1.0/nmax;
	size_t rows = b.size();
	// write positions
	file << rows << std::endl;
	if(dimensions == 1)
		for(size_t i=0; i < rows; i++)
			file << positions[i][0] << " 0.0" << std::endl;
	else if(dimensions == 2)
		for(size_t i=0; i < rows; i++)
			file << positions[i][0] << " " << positions[i][1] << " " << b[i] * scale << std::endl;
	else
		for(size_t i=0; i < rows; i++)
		  file << positions[i][0] << " " << positions[i][1] << " " << positions[i][2] << std::endl;

	file << 1 << std::endl; // show all cons
	// write connections
	for(size_t i=0; i < rows; i++)
	{
		file << i << " " << i << " " << b[i] <<		std::endl;
	}
}
#endif

}
#endif // __H__UG__CPU_ALGEBRA__SPARSEMATRIX_PRINT__
