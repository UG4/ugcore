#include "hyprematrix.h"

namespace ug{

bool HypreMatrix::create_matrix(int nrow, int ncol)
{
	int err = 0;
	err += HYPRE_IJMatrixCreate(MPI_COMM_WORLD, 0, nrow-1,0, ncol-1, &m_hypreA);
	err += HYPRE_IJMatrixSetObjectType(m_hypreA, HYPRE_PARCSR);
	err += HYPRE_IJMatrixInitialize(m_hypreA);
	if(err) return false;
	else return true;
}

bool HypreMatrix::delete_matrix()
{
	int err = HYPRE_IJMatrixDestroy(m_hypreA);
	if(err) return false;
	else
	{
		m_hypreA = NULL;
		return true;
	}
}

bool HypreMatrix::set_values(int nrows, int* ncols, int* rows, int* cols, double* values)
{
	if(m_hypreA == NULL) return false;
	return !(bool)HYPRE_IJMatrixSetValues(m_hypreA, nrows, ncols, rows, cols, values);
}

bool HypreMatrix::add_values(int nrows, int* ncols, int* rows, int* cols, double* values)
{
	if(m_hypreA == NULL) return false;
	return !(bool)HYPRE_IJMatrixAddToValues(m_hypreA, nrows, ncols, rows, cols, values);
}

bool HypreMatrix::set_dirichletrows(int nrows, int* rows)
{
	if(m_hypreA == NULL) return false;
	int one_ncols;
	int* one_cols;
	int* ncols = new int[nrows];

	// Close Matrix (assemble="make ready","finalize") to be able to read columns
	if(HYPRE_IJMatrixAssemble(m_hypreA)){
		delete[] ncols;
		return false;
	}

	void     *object;
	HYPRE_IJMatrixGetObject(m_hypreA, &object);

	// read number of columns for each row
	for(int i=0; i< nrows; i++)
	{
		HYPRE_ParCSRMatrixGetRow((HYPRE_ParCSRMatrix) object,rows[i], &one_ncols, &one_cols, NULL);
		ncols[i] = one_ncols;
		HYPRE_ParCSRMatrixRestoreRow((HYPRE_ParCSRMatrix) object,rows[i], &one_ncols, &one_cols, NULL);
	}

	// count number of matrix entries to be handled
	int nvalues = 0;
	for(int i=0; i<nrows; i++)
	{
		nvalues += ncols[i];
	}

	if(nvalues<=0){
		delete[] ncols;
		return false;
	}

	// create Arrays with length: number of matrix entries to be set
	int* cols = new int[nvalues];
	double* values = new double[nvalues];

	// Get Column number of entries
	int count = 0;
	for(int i=0; i< nrows; i++)
	{
		HYPRE_ParCSRMatrixGetRow((HYPRE_ParCSRMatrix) object,rows[i], &one_ncols, &one_cols, NULL);
		for(int j=0;j < one_ncols; j++)
		{
			cols[count++] = one_cols[j];
		}
		HYPRE_ParCSRMatrixRestoreRow((HYPRE_ParCSRMatrix) object,rows[i], &one_ncols, &one_cols, NULL);
	}

	// Choose new values (i.e. identity row)
	count = 0;
	for(int i=0; i<nrows; i++)
		for(int j=0; j<ncols[i]; j++)
		{
			if(rows[i] == cols[count]) values[count++] = 1.0;
			else values[count++] = 0.0;
		}

	// Reinitialize Matrix
	if(HYPRE_IJMatrixInitialize(m_hypreA)){
		delete[] ncols;
		delete[] cols;
		delete[] values;
		return false;
	}


	// Write new entries
	if(HYPRE_IJMatrixSetValues(m_hypreA, nrows, ncols, rows, cols, values)){
		delete[] ncols;
		delete[] cols;
		delete[] values;
		return false;
	}

	delete[] ncols;
	delete[] cols;
	delete[] values;
	return true;
}

bool HypreMatrix::set(number w)
{
	int ilower, iupper, jlower, jupper, nrows;

	HYPRE_IJMatrixGetLocalRange (m_hypreA, &ilower,&iupper, &jlower, &jupper);

	nrows = iupper - ilower + 1;
	int* rows = new int[nrows];

	if(m_hypreA == NULL) return false;
	int one_ncols;
	int* one_cols;
	int* ncols = new int[nrows];

	// Close Matrix (assemble="make ready","finalize") to be able to read columns
	if(HYPRE_IJMatrixAssemble(m_hypreA)){
		delete[] ncols;
		return false;
	}

	void     *object;
	HYPRE_IJMatrixGetObject(m_hypreA, &object);

	// read number of columns for each row
	for(int i=0; i< nrows; i++)
	{
		HYPRE_ParCSRMatrixGetRow((HYPRE_ParCSRMatrix) object,rows[i], &one_ncols, &one_cols, NULL);
		ncols[i] = one_ncols;
		HYPRE_ParCSRMatrixRestoreRow((HYPRE_ParCSRMatrix) object,rows[i], &one_ncols, &one_cols, NULL);
	}

	// count number of matrix entries to be handled
	int nvalues = 0;
	for(int i=0; i<nrows; i++)
	{
		nvalues += ncols[i];
	}

	if(nvalues<=0){
		delete[] ncols;
		return false;
	}

	// create Arrays with length: number of matrix entries to be set
	int* cols = new int[nvalues];
	double* values = new double[nvalues];

	// Get Column number of entries
	int count = 0;
	for(int i=0; i< nrows; i++)
	{
		HYPRE_ParCSRMatrixGetRow((HYPRE_ParCSRMatrix) object,rows[i], &one_ncols, &one_cols, NULL);
		for(int j=0;j < one_ncols; j++)
		{
			cols[count++] = one_cols[j];
		}
		HYPRE_ParCSRMatrixRestoreRow((HYPRE_ParCSRMatrix) object,rows[i], &one_ncols, &one_cols, NULL);
	}

	// Choose new values (i.e. identity row)
	count = 0;
	for(int i=0; i<nrows; i++)
		for(int j=0; j<ncols[i]; j++)
		{
			values[count++] = w;
		}

	// Reinitialize Matrix
	if(HYPRE_IJMatrixInitialize(m_hypreA)){
		delete[] ncols;
		delete[] cols;
		delete[] values;
		return false;
	}


	// Write new entries
	if(HYPRE_IJMatrixSetValues(m_hypreA, nrows, ncols, rows, cols, values)){
		delete[] ncols;
		delete[] cols;
		delete[] values;
		return false;
	}

	delete[] ncols;
	delete[] cols;
	delete[] values;
	return true;
}


bool HypreMatrix::finalize()
{
	return !(bool)HYPRE_IJMatrixAssemble(m_hypreA);
}

bool HypreMatrix::printToFile(const char* filename)
{
	if(m_hypreA == NULL) return false;
	return !(bool)HYPRE_IJMatrixPrint(m_hypreA, filename);
}

HYPRE_IJMatrix HypreMatrix::getStorage()
{
	return m_hypreA;
}

HypreMatrix::~HypreMatrix()
{
	if(m_hypreA != NULL)
		HYPRE_IJMatrixDestroy(m_hypreA);
}


};
