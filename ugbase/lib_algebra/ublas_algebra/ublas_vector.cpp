
#include "ublas_vector.h"
#include "stdio.h"
#include <cmath>

namespace ug{

bool
UblasVector::
create(uint nentries)
{
	int err = 0;

	m_pVector = new ScalarVector((int)nentries);

	if(m_pVector == NULL) err = 1;

	if(err) return false;
	else return true;
}

bool
UblasVector::
create(const UblasVector& v)
{
	uint nentries = v.size();

	m_pVector = new ScalarVector((int)nentries);

	if(m_pVector == NULL) return false;
	return true;
}


bool UblasVector::destroy()
{
	if(m_pVector != NULL)
		delete m_pVector;

	m_pVector = NULL;
	return true;
}

bool
UblasVector::
set(const local_vector_type& u, const local_index_type& ind)
{
	if(m_pVector == NULL) return false;
	for(std::size_t i = 0; i < ind.size(); ++i)
	{
		UG_ASSERT(ind[i][0] < m_pVector->size() && ind[i][0] >= 0, "ind = " << ind[i][0] << ", size = " << m_pVector->size());
		(*m_pVector)(ind[i][0]) = u[i];
	}
	return true;
}

bool
UblasVector::
add(const local_vector_type& u, const local_index_type& ind)
{
	if(m_pVector == NULL) return false;
	for(std::size_t i = 0; i < ind.size(); ++i)
	{
		assert(ind[i][0] < m_pVector->size() && ind[i][0] >= 0);
		(*m_pVector)(ind[i][0]) += u[i];
	}
	return true;
}

bool
UblasVector::
get(local_vector_type& u, const local_index_type& ind) const
{
	if(m_pVector == NULL) return false;
	for(std::size_t i = 0; i < ind.size(); ++i)
	{
		assert(ind[i][0] < m_pVector->size() && ind[i][0] >= 0);
		u[i] = (*m_pVector)(ind[i][0]);
	}
	return true;
}

bool UblasVector::finalize()
{
	return true;
}

bool UblasVector::printToFile(const char* filename) const
{
	if(m_pVector == NULL) return false;

	FILE* file;
	file = fopen(filename, "w");
	if(file == NULL) return false;

	for(uint i = 0; i < m_pVector->size(); ++i)
	{
		fprintf(file, "%i: %e\n", i, (*m_pVector)(i));
	}

	fclose(file);
	return true;
}

UblasVector::
~UblasVector()
{
	if(m_pVector != NULL)
		delete m_pVector;

	m_pVector = NULL;
}

UblasVector&
UblasVector::
operator+= (const UblasVector& v) {

		(*m_pVector) += *(v.m_pVector);

		return *this;
}

UblasVector&
UblasVector::
operator-= (const UblasVector& v) {

		(*m_pVector) -= *(v.m_pVector);

		return *this;
}

UblasVector&
UblasVector::
operator= (const UblasVector& v) {

		for(unsigned int i = 0; i < v.m_pVector->size(); ++i)
			(*m_pVector)(i) = (*(v.m_pVector))(i);

		return *this;
}

number
UblasVector::
operator *(const UblasVector& v)
{
	number val = 0;
	for(unsigned int i = 0; i < v.m_pVector->size(); ++i)
	{
		val += (*m_pVector)(i) * (*(v.m_pVector))(i);
	}
	return val;
}

number UblasVector::one_norm() const
	{
		double norm = 0;
		for(uint i = 0; i < m_pVector->size(); ++i)
		{
			norm += fabs( (*m_pVector)(i) );
		}

		return norm;
}

number
UblasVector::
two_norm() const
{
	double norm = 0;
	for(uint i = 0; i < m_pVector->size(); ++i)
	{
		norm += (*m_pVector)(i) * (*m_pVector)(i);
	}

	return (number) sqrt(norm);
}

bool
UblasVector::
set(number w)
{
	for(uint i = 0; i < m_pVector->size(); ++i)
	{
		(*m_pVector)(i) = w;
	}

	return true;
}

bool
UblasVector::
operator*=(number w)
{
	for(uint i = 0; i < m_pVector->size(); ++i)
	{
		(*m_pVector)(i) *= w;
	}

	return true;
}

bool
UblasVector::
operator=(number w)
{
	return this->set(w);
}

uint
UblasVector::
size() const
{
	if(m_pVector != NULL) return m_pVector->size();
	else return 0;
}

UblasVector::ScalarVector*
UblasVector::
getStorage()
{
	return m_pVector;
}

const UblasVector::ScalarVector*
UblasVector::
getStorage() const
{
	return m_pVector;
}

std::ostream& operator<< (std::ostream& outStream, const ug::UblasVector& v)
{
	if(v.m_pVector == NULL) return outStream;

	for(std::size_t i = 0; i < v.m_pVector->size(); ++i)
		outStream << "[" << i << "]: " << (*(v.m_pVector))(i) << std::endl;
	return outStream;
}

std::ostream& operator<< (std::ostream& outStream, const ug::UblasVector::local_index_type& ind)
{
	for(std::size_t i = 0; i < ind.size(); ++i)
		outStream << "[" << i << "]: " << ind[i] << std::endl;
	return outStream;
}

} // namespace ug


