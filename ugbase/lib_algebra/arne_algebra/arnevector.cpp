
#include "arnevector.h"
#include "stdio.h"

namespace ug{

bool ArneVector::create_vector(int nentries)
{
	int err = 0;

	_Vector = new ScalarVector(nentries);

	if(_Vector == NULL) err = 1;

	if(err) return false;
	else return true;
}

bool ArneVector::delete_vector()
{
	delete _Vector;
	return true;
}

bool ArneVector::set_values(int nvalues, int* indices, double* values)
{
	if(_Vector == NULL) return false;
	for(int i = 0; i < nvalues; ++i)
	{
		assert(indices[i] < _Vector->size() && indices[i] >= 0);
		(*_Vector)(indices[i]) = values[i];
	}
	return true;
}

bool ArneVector::add_values(int nvalues, int* indices, double* values)
{
	if(_Vector == NULL) return false;
	for(int i = 0; i < nvalues; ++i)
	{
		assert(indices[i] < _Vector->size() && indices[i] >= 0);
		(*_Vector)(indices[i]) += values[i];
	}
	return true;
}

bool ArneVector::get_values(int nvalues, int* indices, double* values) const
{
	if(_Vector == NULL) return false;
	for(int i = 0; i < nvalues; ++i)
	{
		assert(indices[i] < _Vector->size() && indices[i] >= 0);
		values[i] = (*_Vector)(indices[i]);
	}
	return true;
}

bool ArneVector::finalize()
{
	return true;
}

bool ArneVector::printToFile(const char* filename)
{
	if(_Vector == NULL) return false;

	FILE* file;
	file = fopen(filename, "w");
	if(file == NULL) return false;
	
	for(int i = 0; i < _Vector->size(); ++i)
	{
		fprintf(file, "%i: %e\n", i, (*_Vector)(i));
	}
	
	fclose(file);
	return true;
}

ArneVector::~ArneVector()
{
	if(_Vector != NULL)
	delete _Vector;
}

ArneVector& ArneVector::operator+= (const ArneVector& v) {
		
		(*_Vector) += *(v._Vector);

		return *this;
}

ArneVector& ArneVector::operator-= (const ArneVector& v) {

		(*_Vector) -= *(v._Vector);
		
		return *this;
}

number ArneVector::norm2()
	{
		double norm = 0;
		for(int i = 0; i < _Vector->size(); ++i)
		{
			norm += (*_Vector)(i) * (*_Vector)(i);
		}

		return sqrt(norm);
}

bool ArneVector::set(number w)
	{
		for(int i = 0; i < _Vector->size(); ++i)
		{
			(*_Vector)(i) = w;
		}

		return true;
}

int ArneVector::length()
{
	return _Vector->size();
}

ArneVector::ScalarVector* ArneVector::getStorage()
{
	return _Vector;
}

}


