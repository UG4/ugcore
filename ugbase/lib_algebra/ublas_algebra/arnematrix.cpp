#include "arnematrix.h"
#include <iostream>

namespace ug{

bool
ArneMatrix::
create(uint nrow, uint ncol)
{
	_Matrix = new ScalarMatrix(nrow, ncol, 8);

	if(_Matrix == NULL) return false;
	else return true;
}

bool
ArneMatrix::
create(const ArneMatrix& v)
{
	if(_Matrix == NULL) return false;

	uint nrow = v.row_size();
	uint ncol = v.col_size();
	_Matrix = new ScalarMatrix(nrow, ncol, 8);

	if(_Matrix == NULL) return false;
	else return true;
}


bool
ArneMatrix::
destroy()
{
	if(_Matrix != NULL)
		delete _Matrix;
	_Matrix = NULL;
	return true;
}

bool
ArneMatrix::
set(const local_matrix_type& mat, const local_index_type& I, const local_index_type& J)
{
	if(_Matrix == NULL) return false;

	for(uint i = 0; i < I.size(); ++i)
	{
		for(uint j = 0; j < J.size(); ++j)
		{
			(*_Matrix)(I[i][0], J[j][0]) = mat(i,j);
		}
	}
	return true;
}

bool
ArneMatrix::
add(const local_matrix_type& mat, const local_index_type& I, const local_index_type& J)
{
	if(_Matrix == NULL) return false;

	for(uint i = 0; i < I.size(); ++i)
	{
		for(uint j = 0; j < J.size(); ++j)
		{
			(*_Matrix)(I[i][0], J[j][0]) += mat(i,j);
		}
	}
	return true;
}

bool
ArneMatrix::
get(local_matrix_type& mat, const local_index_type& I, const local_index_type& J) const
{
	if(_Matrix == NULL) return false;

	for(uint i = 0; i < I.size(); ++i)
	{
		for(uint j = 0; j < J.size(); ++j)
		{
			mat(i,j) = (*_Matrix)(I[i][0], J[j][0]);
		}
	}
	return true;
}

bool ArneMatrix::set_dirichlet_rows(const local_index_type& I)
{
	if(_Matrix == NULL) return false;

	typedef ublas::matrix_row< ScalarMatrix > row_type;

	for(uint i = 0; i < I.size(); ++i)
	{
		row_type row(*_Matrix, I[i][0]);
		for(row_type::iterator iter_ij = row.begin(); iter_ij != row.end(); ++iter_ij)
		{
			if(iter_ij.index() == I[i][0]) *iter_ij = 1.0;
			else *iter_ij = 0.0;
		}
	}
	return true;
}

bool
ArneMatrix::
set(number w)
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


bool
ArneMatrix::
operator=(number w)
{
	return this->set(w);
}

bool
ArneMatrix::
finalize()
{
	return true;
}

bool
ArneMatrix::
apply(ArneVector&b, const ArneVector& x)
{
	ublas::axpy_prod(*_Matrix, *x.getStorage(), *b.getStorage(), true);
	return true;
}

bool
ArneMatrix::
apply_transposed(ArneVector&b, const ArneVector& x)
{
	ublas::axpy_prod(*x.getStorage(), *_Matrix, *b.getStorage(), true);
	return true;
}

bool
ArneMatrix::
matmul_minus(ArneVector&b, const ArneVector& x)
{
	ArneVector& x2 = const_cast<ArneVector&>(x);

	*x2.getStorage() *= -1.0;
	ublas::axpy_prod(*_Matrix, *x2.getStorage(), *b.getStorage(), false);
	*x2.getStorage() *= -1.0;
	return true;
}

uint
ArneMatrix::
row_size() const
{
	if(_Matrix == NULL) return 0;
	return _Matrix->size1();
}

uint
ArneMatrix::
col_size() const
{
	if(_Matrix == NULL) return 0;
	return _Matrix->size2();
}

bool
ArneMatrix::
printToFile(const char* filename)
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
			fprintf(file, "%i, %i: %e\n", (int) iter_ij.index1(), (int) iter_ij.index2(), *iter_ij);
		}
	}

    fclose(file);
	return true;
}

ArneMatrix::ScalarMatrix*
ArneMatrix::
getStorage()
{
	return _Matrix;
}

const ArneMatrix::ScalarMatrix*
ArneMatrix::
getStorage() const
{
	return _Matrix;
}


ArneMatrix::
~ArneMatrix()
{
	destroy();
}

std::ostream& operator<< (std::ostream& outStream, const ug::ArneMatrix& m)
{
	if(m._Matrix == NULL) return outStream;

	typedef ublas::compressed_matrix<double, ublas::row_major> ScalarMatrix;
    typedef ScalarMatrix::value_type value_type;
    typedef ScalarMatrix::const_iterator1 row_iter_type;
    typedef ScalarMatrix::const_iterator2 entry_iter_type;

    for (row_iter_type rowi = m._Matrix->begin1(); rowi != m._Matrix->end1(); ++rowi)
    {
		for(entry_iter_type iter_ij = rowi.begin(); iter_ij != rowi.end(); ++iter_ij)
		{
			outStream << "[" << iter_ij.index1() << ", " << iter_ij.index2() << "]: " <<  *iter_ij << std::endl;
		}
	}

 	return outStream;
}
};
