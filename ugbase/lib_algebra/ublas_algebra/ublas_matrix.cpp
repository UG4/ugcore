#include "ublas_matrix.h"
#include <iostream>

namespace ug{

bool
UblasMatrix::
create(size_t nrow, size_t ncol)
{
	m_pMatrix = new ScalarMatrix(nrow, ncol, 8);

	if(m_pMatrix == NULL) return false;
	else return true;
}

bool
UblasMatrix::
create(const UblasMatrix& v)
{
	if(m_pMatrix == NULL) return false;

	size_t nrow = v.row_size();
	size_t ncol = v.col_size();
	m_pMatrix = new ScalarMatrix(nrow, ncol, 8);

	if(m_pMatrix == NULL) return false;
	else return true;
}


bool
UblasMatrix::
destroy()
{
	if(m_pMatrix != NULL)
		delete m_pMatrix;
	m_pMatrix = NULL;
	return true;
}

bool
UblasMatrix::
set(const local_matrix_type& mat, const local_index_type& I, const local_index_type& J)
{
	if(m_pMatrix == NULL) return false;

	for(size_t i = 0; i < I.size(); ++i)
	{
		for(size_t j = 0; j < J.size(); ++j)
		{
			(*m_pMatrix)(I[i][0], J[j][0]) = mat(i,j);
		}
	}
	return true;
}

bool
UblasMatrix::
add(const local_matrix_type& mat, const local_index_type& I, const local_index_type& J)
{
	if(m_pMatrix == NULL) return false;

	for(size_t i = 0; i < I.size(); ++i)
	{
		for(size_t j = 0; j < J.size(); ++j)
		{
			(*m_pMatrix)(I[i][0], J[j][0]) += mat(i,j);
		}
	}
	return true;
}

bool
UblasMatrix::
get(local_matrix_type& mat, const local_index_type& I, const local_index_type& J) const
{
	if(m_pMatrix == NULL) return false;

	for(size_t i = 0; i < I.size(); ++i)
	{
		for(size_t j = 0; j < J.size(); ++j)
		{
			mat(i,j) = (*m_pMatrix)(I[i][0], J[j][0]);
		}
	}
	return true;
}

bool UblasMatrix::set_dirichlet_rows(const local_index_type& I)
{
	if(m_pMatrix == NULL) return false;

	typedef ublas::matrix_row< ScalarMatrix > row_type;

	for(size_t i = 0; i < I.size(); ++i)
	{
		row_type row(*m_pMatrix, I[i][0]);
		for(row_type::iterator iter_ij = row.begin(); iter_ij != row.end(); ++iter_ij)
		{
			if(iter_ij.index() == I[i][0]) *iter_ij = 1.0;
			else *iter_ij = 0.0;
		}
	}
	return true;
}

bool
UblasMatrix::
set(number w)
{
	if(m_pMatrix == NULL) return false;

    typedef ScalarMatrix::value_type value_type;

    typedef ScalarMatrix::iterator1 row_iter_type;
    typedef ScalarMatrix::iterator2 entry_iter_type;

    for (row_iter_type rowi = m_pMatrix->begin1(); rowi != m_pMatrix->end1(); ++rowi)
    {
		for(entry_iter_type iter_ij = rowi.begin(); iter_ij != rowi.end(); ++iter_ij)
		{
			*iter_ij = w;
		}
	}
	return true;
}


bool
UblasMatrix::
operator=(number w)
{
	return this->set(w);
}

bool
UblasMatrix::
finalize()
{
	return true;
}

bool
UblasMatrix::
apply(UblasVector&b, const UblasVector& x)
{
	ublas::axpy_prod(*m_pMatrix, *x.getStorage(), *b.getStorage(), true);
	return true;
}

bool
UblasMatrix::
apply_transposed(UblasVector&b, const UblasVector& x)
{
	ublas::axpy_prod(*x.getStorage(), *m_pMatrix, *b.getStorage(), true);
	return true;
}

bool
UblasMatrix::
matmul_minus(UblasVector&b, const UblasVector& x)
{
	UblasVector& x2 = const_cast<UblasVector&>(x);

	*x2.getStorage() *= -1.0;
	ublas::axpy_prod(*m_pMatrix, *x2.getStorage(), *b.getStorage(), false);
	*x2.getStorage() *= -1.0;
	return true;
}

size_t
UblasMatrix::
row_size() const
{
	if(m_pMatrix == NULL) return 0;
	return m_pMatrix->size1();
}

size_t
UblasMatrix::
col_size() const
{
	if(m_pMatrix == NULL) return 0;
	return m_pMatrix->size2();
}

bool
UblasMatrix::
printToFile(const char* filename)
{
	if(m_pMatrix == NULL) return false;
	FILE* file;

	file = fopen(filename, "w");

    typedef ScalarMatrix::value_type value_type;

    typedef ScalarMatrix::const_iterator1 row_iter_type;
    typedef ScalarMatrix::const_iterator2 entry_iter_type;

    for (row_iter_type rowi = m_pMatrix->begin1(); rowi != m_pMatrix->end1(); ++rowi)
    {
		for(entry_iter_type iter_ij = rowi.begin(); iter_ij != rowi.end(); ++iter_ij)
		{
			fprintf(file, "%i, %i: %e\n", (int) iter_ij.index1(), (int) iter_ij.index2(), *iter_ij);
		}
	}

    fclose(file);
	return true;
}

UblasMatrix::ScalarMatrix*
UblasMatrix::
getStorage()
{
	return m_pMatrix;
}

const UblasMatrix::ScalarMatrix*
UblasMatrix::
getStorage() const
{
	return m_pMatrix;
}


UblasMatrix::
~UblasMatrix()
{
	destroy();
}

std::ostream& operator<< (std::ostream& outStream, const ug::UblasMatrix& m)
{
	if(m.m_pMatrix == NULL) return outStream;

	typedef ublas::compressed_matrix<double, ublas::row_major> ScalarMatrix;
    typedef ScalarMatrix::value_type value_type;
    typedef ScalarMatrix::const_iterator1 row_iter_type;
    typedef ScalarMatrix::const_iterator2 entry_iter_type;

    for (row_iter_type rowi = m.m_pMatrix->begin1(); rowi != m.m_pMatrix->end1(); ++rowi)
    {
		for(entry_iter_type iter_ij = rowi.begin(); iter_ij != rowi.end(); ++iter_ij)
		{
			outStream << "[" << iter_ij.index1() << ", " << iter_ij.index2() << "]: " <<  *iter_ij << std::endl;
		}
	}

 	return outStream;
}
};
