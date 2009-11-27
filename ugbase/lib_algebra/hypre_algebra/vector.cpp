
#include "vector.h"

namespace ug{

bool Vector::create_vector(int nentries)
{
	int err = 0;
	err += HYPRE_IJVectorCreate(MPI_COMM_WORLD, 0, nentries-1, &m_hyprex);
	err += HYPRE_IJVectorSetObjectType(m_hyprex, HYPRE_PARCSR);
	err += HYPRE_IJVectorInitialize(m_hyprex);
	if(err) return false;
	else return true;
}

bool Vector::delete_vector()
{
	int err = HYPRE_IJVectorDestroy(m_hyprex);
	if(err) return false;
	else
	{
		m_hyprex = NULL;
				return true;
	}
}

bool Vector::set_values(int nvalues, int* indices, double* values)
{
	if(m_hyprex == NULL) return false;
	return !(bool)HYPRE_IJVectorSetValues(m_hyprex, nvalues, indices, values);
}

bool Vector::add_values(int nvalues, int* indices, double* values)
{
	if(m_hyprex == NULL) return false;
	return !(bool)HYPRE_IJVectorAddToValues(m_hyprex, nvalues, indices, values);
}

bool Vector::get_values(int nvalues, int* indices, double* values) const
{
	if(m_hyprex == NULL) return false;
	return !(bool)HYPRE_IJVectorGetValues(m_hyprex, nvalues, indices, values);
}

bool Vector::finalize()
{
	return !(bool)HYPRE_IJVectorAssemble(m_hyprex);
}

bool Vector::printToFile(const char* filename)
{
	if(m_hyprex == NULL) return false;
	return !(bool)HYPRE_IJVectorPrint(m_hyprex, filename);
}

HYPRE_IJVector Vector::getStorage()
{
	return m_hyprex;
}

Vector::~Vector()
{
	if(m_hyprex != NULL)
		HYPRE_IJVectorDestroy(m_hyprex);
}


}


