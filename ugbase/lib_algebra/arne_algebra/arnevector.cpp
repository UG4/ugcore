
#include "arnevector.h"
#include "stdio.h"
#include <cmath>

namespace ug{

bool
ArneVector::
create(uint nentries)
{
	int err = 0;

	_Vector = new ScalarVector((int)nentries);

	if(_Vector == NULL) err = 1;

	if(err) return false;
	else return true;
}

bool
ArneVector::
create(const ArneVector& v)
{
	uint nentries = v.size();

	_Vector = new ScalarVector((int)nentries);

	if(_Vector == NULL) return false;
	return true;
}


bool ArneVector::destroy()
{
	delete _Vector;
	return true;
}

bool
ArneVector::
set(const local_vector_type& u, const local_index_type& ind)
{
	if(_Vector == NULL) return false;
	for(std::size_t i = 0; i < ind.size(); ++i)
	{
		assert(ind[i][0] < _Vector->size() && ind[i][0] >= 0);
		(*_Vector)(ind[i][0]) = u[i];
	}
	return true;
}

bool
ArneVector::
add(const local_vector_type& u, const local_index_type& ind)
{
	if(_Vector == NULL) return false;
	for(std::size_t i = 0; i < ind.size(); ++i)
	{
		assert(ind[i][0] < _Vector->size() && ind[i][0] >= 0);
		(*_Vector)(ind[i][0]) += u[i];
	}
	return true;
}

bool
ArneVector::
get(local_vector_type& u, const local_index_type& ind) const
{
	if(_Vector == NULL) return false;
	for(std::size_t i = 0; i < ind.size(); ++i)
	{
		assert(ind[i][0] < _Vector->size() && ind[i][0] >= 0);
		u[i] = (*_Vector)(ind[i][0]);
	}
	return true;
}

bool ArneVector::finalize()
{
	return true;
}

bool ArneVector::printToFile(const char* filename) const
{
	if(_Vector == NULL) return false;

	FILE* file;
	file = fopen(filename, "w");
	if(file == NULL) return false;

	for(uint i = 0; i < _Vector->size(); ++i)
	{
		fprintf(file, "%i: %e\n", i, (*_Vector)(i));
	}

	fclose(file);
	return true;
}

ArneVector::
~ArneVector()
{
	if(_Vector != NULL)
	delete _Vector;
}

ArneVector&
ArneVector::
operator+= (const ArneVector& v) {

		(*_Vector) += *(v._Vector);

		return *this;
}

ArneVector&
ArneVector::
operator-= (const ArneVector& v) {

		(*_Vector) -= *(v._Vector);

		return *this;
}

ArneVector&
ArneVector::
operator= (const ArneVector& v) {

		for(unsigned int i = 0; i < v._Vector->size(); ++i)
			(*_Vector)(i) = (*(v._Vector))(i);

		return *this;
}

number
ArneVector::
operator *(const ArneVector& v)
{
	number val = 0;
	for(unsigned int i = 0; i < v._Vector->size(); ++i)
	{
		val += (*_Vector)(i) * (*(v._Vector))(i);
	}
	return val;
}

number ArneVector::one_norm() const
	{
		double norm = 0;
		for(uint i = 0; i < _Vector->size(); ++i)
		{
			norm += abs( (*_Vector)(i) );
		}

		return norm;
}

number
ArneVector::
two_norm() const
{
	double norm = 0;
	for(uint i = 0; i < _Vector->size(); ++i)
	{
		norm += (*_Vector)(i) * (*_Vector)(i);
	}

	return sqrt(norm);
}

bool
ArneVector::
set(number w)
{
	for(uint i = 0; i < _Vector->size(); ++i)
	{
		(*_Vector)(i) = w;
	}

	return true;
}

bool
ArneVector::
operator*=(number w)
{
	for(uint i = 0; i < _Vector->size(); ++i)
	{
		(*_Vector)(i) *= w;
	}

	return true;
}

bool
ArneVector::
operator=(number w)
{
	return this->set(w);
}

uint
ArneVector::
size() const
{
	return _Vector->size();
}

ArneVector::ScalarVector*
ArneVector::
getStorage()
{
	return _Vector;
}

const ArneVector::ScalarVector*
ArneVector::
getStorage() const
{
	return _Vector;
}

} // namespace ug


