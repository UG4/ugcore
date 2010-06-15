
#include "hyprevector.h"
#include <cmath>

namespace ug{

bool HypreVector::create(size_t nentries)
{
	UG_DLOG(LIB_ALG_VECTOR, 1, "Creating Hypre Vector of length " << nentries << ". \n");
	int n = (int) nentries;
	int err = 0;
	err += HYPRE_IJVectorCreate(MPI_COMM_WORLD, 0, n-1, &m_hyprex);
	err += HYPRE_IJVectorSetObjectType(m_hyprex, HYPRE_PARCSR);
	err += HYPRE_IJVectorInitialize(m_hyprex);
	if(err) return false;
	else return true;
}

bool HypreVector::destroy()
{
	UG_DLOG(LIB_ALG_VECTOR, 1, "Destroying Hypre Vector. \n");
	int err = HYPRE_IJVectorDestroy(m_hyprex);
	if(err) return false;
	else
	{
		m_hyprex = NULL;
				return true;
	}
}


bool HypreVector::set(const local_vector_type& u, local_index_type& ind)
{
	bool b;

	const int nvalues = ind.size();
	int *indices = new int[nvalues];
	double *values = new double[nvalues];

	UG_DLOG(LIB_ALG_VECTOR, 4, "\n ---- Setting values in Hypre Vector---- \n");
	for(int i = 0; i < nvalues; ++i)
	{
		indices[i] = ind[i][0];
		values[i] = (double) u[i];
		UG_DLOG(LIB_ALG_VECTOR, 4, "Hypre-Index: " << indices[i] << ", Value: " << values[i] << "\n");
	}
	UG_DLOG(LIB_ALG_VECTOR, 4, " ---- End ----\n");

	b = set_values(nvalues, indices, values);

	delete[] indices;
	delete[] values;

	return b;
}
bool HypreVector::add(const local_vector_type& u, local_index_type& ind)
{
	bool b;

	const int nvalues = ind.size();
	int *indices = new int[nvalues];
	double *values = new double[nvalues];

	UG_DLOG(LIB_ALG_VECTOR, 4, "\n ---- Adding values in Hypre Vector---- \n");
	for(int i = 0; i < nvalues; ++i)
	{
		indices[i] = ind[i][0];
		values[i] = (double) u[i];
		UG_DLOG(LIB_ALG_VECTOR, 4, "Hypre-Index: " << indices[i] << ", Value: " << values[i] << "\n");
	}
	UG_DLOG(LIB_ALG_VECTOR, 4, " ---- End ----\n");

	b = add_values(nvalues, indices, values);

	delete[] indices;
	delete[] values;

	return b;
}

bool HypreVector::get(local_vector_type& u, local_index_type& ind) const
{
	bool b;

	const int nvalues = ind.size();
	int *indices = new int[nvalues];
	double *values = new double[nvalues];

	for(int i = 0; i < nvalues; ++i)
	{
		indices[i] = ind[i][0];
	}

	b = get_values(nvalues, indices, values);

	UG_DLOG(LIB_ALG_VECTOR, 4, "\n ---- Getting values in Hypre Vector---- \n");
	for(int i = 0; i < nvalues; ++i)
	{
		u[i] = (number) values[i];
		UG_DLOG(LIB_ALG_VECTOR, 4, "Hypre-Index: " << indices[i] << ", Value: " << values[i] << " [Converted to value: " << u[i] << "]\n");
	}
	UG_DLOG(LIB_ALG_VECTOR, 4, " ---- End ----\n");

	delete[] indices;
	delete[] values;

	return b;
}


bool HypreVector::set_values(int nvalues, int* indices, double* values)
{
	if(m_hyprex == NULL) return false;
	return !(bool)HYPRE_IJVectorSetValues(m_hyprex, nvalues, indices, values);
}

bool HypreVector::add_values(int nvalues, int* indices, double* values)
{
	if(m_hyprex == NULL) return false;
	return !(bool)HYPRE_IJVectorAddToValues(m_hyprex, nvalues, indices, values);
}

bool HypreVector::get_values(int nvalues, int* indices, double* values) const
{
	if(m_hyprex == NULL) return false;
	return !(bool)HYPRE_IJVectorGetValues(m_hyprex, nvalues, indices, values);
}

bool HypreVector::finalize()
{
	return !(bool)HYPRE_IJVectorAssemble(m_hyprex);
}

bool HypreVector::printToFile(const char* filename)
{
	if(m_hyprex == NULL) return false;
	return !(bool)HYPRE_IJVectorPrint(m_hyprex, filename);
}

HYPRE_IJVector HypreVector::getStorage()
{
	return m_hyprex;
}

HypreVector::~HypreVector()
{
	if(m_hyprex != NULL)
		HYPRE_IJVectorDestroy(m_hyprex);
}


	HypreVector& HypreVector::operator+= (const HypreVector& v) {
		int jlower, jupper, nvalues;

		HYPRE_IJVectorGetLocalRange(m_hyprex, &jlower, &jupper);

		nvalues = jupper - jlower + 1;
		double* values = new double[nvalues];
		int* indices = new int[nvalues];

		for(int i = 0; i < nvalues; ++i)
		{
			indices[i] = i;
		}

		if( v.get_values(nvalues, indices, values) == false)
			assert(0 && "ERROR in operator+=");

		if(add_values(nvalues, indices, values) == false)
			assert(0 && "ERROR in operator+=");

		return *this;
	}

	HypreVector& HypreVector::operator-= (const HypreVector& v) {
		int jlower, jupper, nvalues;

		HYPRE_IJVectorGetLocalRange(m_hyprex, &jlower, &jupper);

		nvalues = jupper - jlower + 1;
		double* values = new double[nvalues];
		int* indices = new int[nvalues];

		for(int i = 0; i < nvalues; ++i)
		{
			indices[i] = i;
		}

		if(v.get_values(nvalues, indices, values) == false)
			assert(0 && "ERROR in operator-=");

		for(int i = 0; i < nvalues; ++i)
		{
			values[i] *= -1;
		}

		if(add_values(nvalues, indices, values) == false)
			assert(0 && "ERROR in operator-=");

		return *this;
	}

number HypreVector::two_norm() const
{
		int jlower, jupper, nvalues;

		HYPRE_IJVectorGetLocalRange(m_hyprex, &jlower, &jupper);

		nvalues = jupper - jlower + 1;
		double* values = new double[nvalues];
		int* indices = new int[nvalues];

		for(int i = 0; i < nvalues; ++i)
		{
			indices[i] = i;
		}

		if(get_values(nvalues, indices, values) == false) return false;

		double norm = 0;
		for(int i = 0; i < nvalues; ++i)
		{
			norm += values[i]*values[i];
		}

		return sqrt(norm);
}

number HypreVector::one_norm() const
{
		int jlower, jupper, nvalues;

		HYPRE_IJVectorGetLocalRange(m_hyprex, &jlower, &jupper);

		nvalues = jupper - jlower + 1;
		double* values = new double[nvalues];
		int* indices = new int[nvalues];

		for(int i = 0; i < nvalues; ++i)
		{
			indices[i] = i;
		}

		if(get_values(nvalues, indices, values) == false) return false;

		double norm = 0;
		for(int i = 0; i < nvalues; ++i)
		{
			norm += fabs((double)values[i]);
		}

		return sqrt(norm);
}

bool HypreVector::set(number w)
{
		int jlower, jupper, nvalues;

		HYPRE_IJVectorGetLocalRange(m_hyprex, &jlower, &jupper);

		nvalues = jupper - jlower + 1;
		double* values = new double[nvalues];
		int* indices = new int[nvalues];

		for(int i = 0; i < nvalues; ++i)
		{
			indices[i] = i;
			values[i] = w;
		}

		if(set_values(nvalues, indices, values) == false) return false;

		return true;
}

size_t HypreVector::size() const
{
		int jlower, jupper;

		HYPRE_IJVectorGetLocalRange(m_hyprex, &jlower, &jupper);

		return jupper - jlower + 1;
}

} // namespace ug
