/*
 * vector.h
 *
 *  Created on: 02.07.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIB_ALGEBRA__VECTOR__
#define __H__LIB_ALGEBRA__VECTOR__

#include <HYPRE.h>
#include <_hypre_utilities.h>
#include <HYPRE_krylov.h>
#include <HYPRE_parcsr_ls.h>
#include "../../common/types.h"
#include <cmath>
#include "assert.h"

namespace ug{

class Vector{

	public:
	bool create_vector(int nentries);

	bool delete_vector();

	bool set_values(int nvalues, int* indices, double* values);

	bool add_values(int nvalues, int* indices, double* values);

	bool get_values(int nvalues, int* indices, double* values) const;

	bool finalize();

	bool printToFile(const char* filename);

	HYPRE_IJVector getStorage();

	~Vector();

	Vector& operator+= (const Vector& v) {
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

	Vector& operator-= (const Vector& v) {
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

	number norm2()
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
	bool set(number w)
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
	int length()
	{
		int jlower, jupper, nvalues;

		HYPRE_IJVectorGetLocalRange(m_hyprex, &jlower, &jupper);

		return jupper - jlower + 1;
	}


	private:
		HYPRE_IJVector m_hyprex;


};

}

#endif /* __H__LIB_ALGEBRA__VECTOR__ */
