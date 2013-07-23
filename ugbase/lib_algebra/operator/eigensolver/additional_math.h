/**
 * \file ugbase/lib_algebra/operator/eigensolver/additional_math.h
 *
 *  \date	 	28.03.2013
 *  \author 	mrupp
 */

#ifndef __UG__ADDITIONAL_MATH_H__
#define __UG__ADDITIONAL_MATH_H__

#include "smart_ptr_vector.h"

#define PINVIT_PROFILE_FUNC() PROFILE_FUNC_GROUP("pinvit algebra")
#define PINVIT_PROFILE_BEGIN(t) PROFILE_BEGIN_GROUP(t, "pinvit algebra")
#define PINVIT_PROFILE_END() PROFILE_END()

namespace ug{

/*template<typename mat_type, typename vec_type, typename densematrix_type>
void MultiEnergyProd(const SparseMatrix<mat_type> &A,
			Vector<vec_type> **x,
			DenseMatrix<densematrix_type> &rA, size_t n)
{
	UG_ASSERT(n == rA.num_rows() && n == rA.num_cols(), "");
	vec_type Ai_xc;
	rA = 0.0;
	for(size_t i=0; i<A.num_rows(); i++)
	{
		for(size_t c=0; c<n; c++)
		{
			// Ai_xc = A[i] * x[c].
			Ai_xc = 0.0;
			A.mat_mult_add_row(i, Ai_xc, 1.0, (*x[c]));
			for(size_t r=0; r<n; r++)
				rA(c, r) += VecDot((*x[r])[i], Ai_xc);
		}
	}
}*/


inline bool absCompare(double a, double b)
{
	return abs(a) < abs(b);
}



template<typename vector_type, typename densematrix_type>
void MultiScalProd(vector_type &px,
			DenseMatrix<densematrix_type> &rA, size_t n)
{
	PINVIT_PROFILE_FUNC();
	UG_ASSERT(0, "");
	UG_ASSERT(n == rA.num_rows() && n == rA.num_cols(), "");
	for(size_t r=0; r<n; r++)
		for(size_t c=r; c<n; c++)
			rA(r, c) = px[c]->dotprod(*px[r]);

	for(size_t r=0; r<n; r++)
		for(size_t c=0; c<r; c++)
			rA(r,c) = rA(c, r);
}

template<typename matrix_type, typename vector_type>
double EnergyProd(vector_type &v1, matrix_type &A, vector_type &v2, vector_type &tmp)
{
	PINVIT_PROFILE_FUNC();

#ifdef UG_PARALLEL
	pcl::ProcessCommunicator pc;
	v2.change_storage_type(PST_CONSISTENT);
#endif
	A.apply(tmp, v2);
	// tmp is additive, v1 is consistent
	double a = v1.dotprod(tmp);
	//UG_LOG("EnergyProd " << a << "\n");

	return a;
}

template<typename matrix_type, typename vector_type>
double EnergyProd(vector_type &v1, matrix_type &A, vector_type &v2)
{
	vector_type t;
	CloneVector(t, v1);
	return EnergyProd(v1, A, v2, t);
}

template<typename matrix_type, typename vector_type>
double EnergyNorm(vector_type &x, matrix_type &A, vector_type &tmp)
{
	return sqrt( EnergyProd(x, A, x, tmp) );
}

template<typename matrix_type, typename vector_type>
double EnergyNorm(vector_type &x, matrix_type &A)
{
	vector_type tmp;
	CloneVector(tmp, x);
	return sqrt( EnergyProd(x, A, x, tmp) );
}


template<typename matrix_type, typename vector_type, typename densematrix_type>
void MultiEnergyProd(matrix_type &A,
			SmartPtrVector<vector_type> &px,
			DenseMatrix<densematrix_type> &rA, size_t n)
{
	PINVIT_PROFILE_FUNC();
#ifdef UG_PARALLEL
	pcl::ProcessCommunicator pc;
#endif
	UG_ASSERT(n == rA.num_rows() && n == rA.num_cols(), "");
	vector_type t;
	CloneVector(t, *px[0]);


	for(size_t r=0; r<n; r++)
	{
		// todo: why is SparseMatrix<T>::apply not const ?!?
#ifdef UG_PARALLEL
		px[r]->change_storage_type(PST_CONSISTENT);
#endif
		A.apply(t, *px[r]);
#ifdef UG_PARALLEL
		t.change_storage_type(PST_CONSISTENT);
#endif
		for(size_t c=r; c<n; c++)
		{
			//rA(r, c) = VecProd((*px[c]), t);
			rA(r, c) = px[c]->dotprod(t);
			//UG_LOG("MultiEnergyProd : (" << r << ", " << c << ") = " << rA(r, c) << "\n");
		}
	}


	for(size_t r=0; r<n; r++)
		for(size_t c=0; c<r; c++)
			rA(r,c) = rA(c, r);
}


template<typename tvector>
void PrintStorageType(const tvector &v)
{
#ifdef UG_PARALLEL
	if(v.has_storage_type(PST_UNDEFINED))
		UG_LOG("PST_UNDEFINED ");
	if(v.has_storage_type(PST_CONSISTENT))
		UG_LOG("PST_CONSISTENT ");
	if(v.has_storage_type(PST_ADDITIVE))
		UG_LOG("PST_ADDITIVE ");
	if(v.has_storage_type(PST_UNIQUE))
		UG_LOG("PST_UNIQUE ");
#else
	UG_LOG("serial ");
#endif
}


template<typename matrix_type>
void PrintMatrix(const matrix_type &mat, const char *name)
{
	UG_LOG(name << ":\n" << name << " := matrix([\n");
	for(size_t r=0; r<mat.num_rows(); r++)
	{
		UG_LOG("[");
		for(size_t c=0; c<mat.num_cols(); c++)
		{
			UG_LOG(mat(r, c));
			if(c < mat.num_cols()-1) UG_LOG(",\t");
		}
		UG_LOG("]\n");
	}
	UG_LOG("]);\n");

}

template<typename matrix_type>
void PrintMaple(const matrix_type &mat, const char *name)
{
	UG_LOG(name << ":\n" << name << " := matrix([");
	for(size_t r=0; r<mat.num_rows(); r++)
	{
		UG_LOG("[");
		for(size_t c=0; c<mat.num_cols(); c++)
		{
			UG_LOG(mat(r, c));
			if(c < mat.num_cols()-1) UG_LOG(", ");
		}
		UG_LOG("]");
		if(r < mat.num_rows()-1) UG_LOG(", ");
	}
	UG_LOG("]);\n");
	UG_LOG("(" << mat.num_rows() << " x " << mat.num_cols() << ")\n");

}

template<typename T>
void MemSwap(T &a, T &b)
{
	char c[sizeof(T)];
	memcpy(c, &a, sizeof(T));
	memcpy(&a, &b, sizeof(T));
	memcpy(&b, c, sizeof(T));
}

} // namespace ug

#endif /* __UG__ADDITIONAL_MATH_H__ */
