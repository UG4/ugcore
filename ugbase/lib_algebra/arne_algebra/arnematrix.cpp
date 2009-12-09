#include "arnematrix.h"
#include <iostream>

namespace ug{

bool ArneMatrix::create_matrix(int nrow, int ncol)
{
	_Matrix = new ScalarMatrix(nrow, ncol);

	if(_Matrix == NULL) return false;
	else return true;
}

bool ArneMatrix::delete_matrix()
{
	if(_Matrix != NULL)
		delete _Matrix;
	return true;
}

bool ArneMatrix::set_values(int nrows, int* ncols, int* rows, int* cols, double* values)
{
	if(_Matrix == NULL) return false;

	int j = 0;
	for(int i = 0; i < nrows; ++i)
	{
		for(int k = 0; k < ncols[i]; ++k)
		{
			(*_Matrix)(rows[i], cols[j]) = values[j];
			++j;
		}
	}
	return true;
}

bool ArneMatrix::add_values(int nrows, int* ncols, int* rows, int* cols, double* values)
{
	if(_Matrix == NULL) return false;

	int j = 0;
	for(int i = 0; i < nrows; ++i)
	{
		for(int k = 0; k < ncols[i]; ++k)
		{
			(*_Matrix)(rows[i], cols[j]) += values[j];
			++j;
		}
	}
	return true;
}

bool ArneMatrix::set_dirichletrows(int nrows, int* rows)
{
	if(_Matrix == NULL) return false;

	typedef ublas::matrix_row< ScalarMatrix > row_type;

	for(int i = 0; i < nrows; ++i)
	{
		row_type row(*_Matrix, rows[i]);
		for(row_type::iterator iter_ij = row.begin(); iter_ij != row.end(); ++iter_ij)
		{
			if(iter_ij.index() == rows[i]) *iter_ij = 1.0;
			else *iter_ij = 0.0;
		}
	}
	return true;
}

bool ArneMatrix::set(number w)
{
	if(_Matrix == NULL) return false;

    typedef ScalarMatrix::value_type value_type;

    typedef ScalarMatrix::iterator1 row_iter_type;
    typedef ScalarMatrix::iterator2 entry_iter_type;

    for (row_iter_type rowi = _Matrix->begin1(); rowi != _Matrix->end1(); ++rowi)
    {
		for(entry_iter_type iter_ij = rowi.begin(); iter_ij != rowi.end(); ++iter_ij)
		{
			*iter_ij = w;
		}
	}
	return true;
}


bool ArneMatrix::finalize()
{
	return true;
}

bool ArneMatrix::printToFile(const char* filename)
{
	if(_Matrix == NULL) return false;
	FILE* file;

	file = fopen(filename, "w");

    typedef ScalarMatrix::value_type value_type;

    typedef ScalarMatrix::const_iterator1 row_iter_type;
    typedef ScalarMatrix::const_iterator2 entry_iter_type;

    for (row_iter_type rowi = _Matrix->begin1(); rowi != _Matrix->end1(); ++rowi)
    {
		for(entry_iter_type iter_ij = rowi.begin(); iter_ij != rowi.end(); ++iter_ij)
		{
			fprintf(file, "%i, %i: %e\n", iter_ij.index1(), iter_ij.index2(), *iter_ij);
		}
	}

	return true;
}

ArneMatrix::ScalarMatrix* ArneMatrix::getStorage()
{
	return _Matrix;
}

ArneMatrix::~ArneMatrix()
{
	if(_Matrix != NULL)
	delete _Matrix;
}


};
