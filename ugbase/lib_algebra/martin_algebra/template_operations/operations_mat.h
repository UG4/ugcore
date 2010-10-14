/*
 * operations_mat.h
 *
 *  Created on: 29.09.2010
 *      Author: mrupp
 */

#ifndef __H__UG__MARTIN_ALGEBRA__OPERATIONS_MAT__
#define __H__UG__MARTIN_ALGEBRA__OPERATIONS_MAT__

namespace ug
{

// operations for matrices
//-----------------------------------------------------------------------------
// these functions execute matrix operations by using the operations on the elements (when TE_MAT_RowApplicable)

//template<typename vector_t, typename matrix_t>
//inline void MatMult(vector_t &dest, const number &alpha1, const matrix_t &A1, const vector_t &v1) // to be specialized

//! calculates dest = alpha1 * A1 *v1;
template<typename vector_t, typename matrix_t>
inline void MatMult(vector_t &dest, const number &alpha1, const TE_MAT_RowApplicable<matrix_t> &A1, const vector_t &v1)
{
	if(OPERATIONS_VERBOSE) cout << "dest: " << dest << endl << "\t= " << alpha1 << " * " << A1.cast() << "\t* " << v1 << endl;
	for(size_t i=0; i<dest.size(); i++)
	{
		dest[i] = 0.0;
		A1.cast().mat_mult_add_row(i, dest[i], alpha1, v1);
	}
}

//! calculates dest = alpha1 * A1 *v1 + beta1*w1;
template<typename vector_t, typename matrix_t>
inline void MatMultAdd(vector_t &dest, const number &alpha1, const TE_MAT<matrix_t> &A1, const vector_t &v1, const number &beta1, const vector_t &w1)
{
	if(OPERATIONS_VERBOSE) cout << "dest: " << dest << endl <<
		"\t= " << alpha1 << " * " << A1.cast() << "\t * " << v1 << endl <<
		"\t+ " << beta1 << " * " << w1 << endl;

	VecScaleAssign(dest, beta1, w1);
	MatMultAdd(dest, alpha1, A1.cast(), v1);
}

//! calculates dest = alpha1 * A1 *v1 + beta1*w1;
template<typename vector_t, typename matrix_t>
inline void MatMultAdd(vector_t &dest, const number &alpha1, const TE_MAT_RowApplicable<matrix_t> &A1, const vector_t &v1, const number &beta1, const vector_t &w1)
{
	if(OPERATIONS_VERBOSE) cout << "dest: " << dest << endl <<
		"\t= " << alpha1 << " * " << A1.cast() << "\t * " << v1 << endl <<
		"\t+ " << beta1 << " * " << w1 << endl;

	for(size_t i=0; i<dest.size(); i++)
	{
		VecScaleAssign(dest[i], beta1, w1[i]);
		A1.cast().mat_mult_add_row(i, dest[i], alpha1, v1);
	}
}

//! calculates dest = alpha*A1*v1 + beta1*w1, and norm = norm^2(dest)
template<typename vector_t, typename matrix_t>
inline void MatMultAddNorm(vector_t &dest,
		const number &alpha1, const TE_MAT_RowApplicable<matrix_t> &A1, const vector_t &v1,
		const number &beta1, const vector_t &w1,
		double &norm)
{
	if(OPERATIONS_VERBOSE) cout << "dest: " << dest << endl <<
			"\t= " << alpha1 << " * " << A1.cast() << "\t * " << v1 << endl <<
			"\t+ " << beta1 << " * " << w1 << endl;

	norm=0;
	for(size_t i=0; i<dest.size(); i++)
	{
		VecScaleAssign(dest[i], beta1, w1[i]);
		A1.cast().mat_mult_add_row(i, dest[i], alpha1, v1);
		VecNormSquaredAdd(dest[i], norm);
	}
}


//! calculates dest = alpha1 * A1 *v1 + alpha2 * A2*v2;
template<typename vector_t, typename matrix_t>
inline void MatMultAdd(vector_t &dest, const number &alpha1, const TE_MAT_RowApplicable<matrix_t> &A1, const vector_t &v1, const number &alpha2, const TE_MAT_RowApplicable<matrix_t> &A2, const vector_t &v2)
{
	if(OPERATIONS_VERBOSE) cout << "dest: " << dest << endl <<
		"\t= " << alpha1 << " * " << A1.cast() << "\t * " << v1 << endl <<
		"\t+ " << alpha2 << " * " << A2.cast() << endl << "\t\t * " << v2 << endl;

	for(size_t i=0; i<dest.size(); i++)
	{
		dest[i] = 0.0;
		A1.cast().mat_mult_add_row(i, dest[i], alpha1, v1);
		A2.cast().mat_mult_add_row(i, dest[i], alpha2, v2);
	}
}

//! calculates dest = alpha1 * A1 *v1 + alpha2 * A2*v2;
template<typename vector_t, typename matrix_t>
inline void MatMultAdd(vector_t &dest, const number &alpha1, const TE_MAT<matrix_t> &A1, const vector_t &v1, const number &alpha2, const TE_MAT<matrix_t> &A2, const vector_t &v2)
{
	if(OPERATIONS_VERBOSE) cout << "dest: " << dest << endl <<
		"\t= " << alpha1 << " * " << A1.cast() << "\t * " << v1 << endl <<
		"\t+ " << alpha2 << " * " << A2.cast() << endl << "\t\t * " << v2 << endl;

	MatMult(dest, alpha1, A1.cast(), v1);
	MatMultAdd(dest, alpha2, A2, v2, 1.0, dest);
}


//! calculates dest = alpha1 * A1 *v1 + alpha2 * A2*v2 + beta1*w1;
template<typename vector_t, typename matrix_t>
inline void MatMultAdd(vector_t &dest, const number &alpha1, const TE_MAT_RowApplicable<matrix_t> &A1, const vector_t &v1, const number &alpha2, const TE_MAT_RowApplicable<matrix_t> &A2, const vector_t &v2,
				const number &beta1, const vector_t &w1)
{
	if(OPERATIONS_VERBOSE) cout << "dest: " << dest << endl <<
		"\t= " << alpha1 << " * " << A1.cast() << "\t * " << v1 << endl <<
		"\t+ " << alpha2 << " * " << A2.cast() << "\t * " << v2 << endl <<
		"\t+ " << beta1 << " * " << w1 << endl;

	for(size_t i=0; i<dest.size(); i++)
	{
		VecScaleAssign(dest[i], beta1, w1[i]);
		A1.cast().mat_mult_add_row(i, dest[i], alpha1, v1);
		A2.cast().mat_mult_add_row(i, dest[i], alpha2, v2);
	}
}

//! calculates dest = alpha1 * A1 *v1 + alpha2 * A2*v2 + beta1*w1;
template<typename vector_t, typename matrix_t>
inline void MatMultAdd(vector_t &dest, const number &alpha1, const TE_MAT<matrix_t> &A1, const vector_t &v1, const number &alpha2, const TE_MAT<matrix_t> &A2, const vector_t &v2,
				const number &beta1, const vector_t &w1)
{
	if(OPERATIONS_VERBOSE) cout << "dest: " << dest << endl <<
		"\t= " << alpha1 << " * " << A1.cast() << "\t * " << v1 << endl <<
		"\t+ " << alpha2 << " * " << A2.cast() << "\t * " << v2 << endl <<
		"\t+ " << beta1 << " * " << w1 << endl;;

	MatMultAdd(dest, alpha1, A1, v1, alpha2, A2, v2);
	VecScaleAdd(dest, 1.0, dest, beta1, w1);
}


//! calculates dest = alpha1 * A1 *v1 + beta1*w1 + beta2*w2;
template<typename vector_t, typename matrix_t>
inline void MatMultAdd(vector_t &dest, const number &alpha1, const TE_MAT_RowApplicable<matrix_t> &A1, const vector_t &v1, const number &beta1, const vector_t &w1, const number &beta2, const vector_t &w2)
{
	if(OPERATIONS_VERBOSE) cout << "dest: " << dest << endl <<
		"\t= " << alpha1 << " * " << A1.cast() << "\t * " << v1 << endl <<
		"\t+ " << beta1 << " * " << w1 << endl <<
		"\t+ " << beta2 << " * " << w2 << endl;

	for(size_t i=0; i<dest.size(); i++)
	{
		VecScaleAdd(dest[i], beta1, w1[i], beta2, w2[i]);
		A1.cast().mat_mult_add_row(i, dest[i], alpha1, v1);
	}
}

//! calculates dest = alpha1 * A1 *v1 + beta1*w1 + beta2*w2;
template<typename vector_t, typename matrix_t>
inline void MatMultAdd(vector_t &dest, const number &alpha1, const TE_MAT<matrix_t> &A1, const vector_t &v1, const number &beta1, const vector_t &w1,
				const number &beta2, const vector_t &w2)
{
	if(OPERATIONS_VERBOSE) cout << "dest: " << dest << endl <<
		"\t= " << alpha1 << " * " << A1.cast() << "\t * " << v1 << endl <<
		"\t+ " << beta1 << " * " << w1 << endl <<
		"\t+ " << beta2 << " * " << w2 << endl;

	MatMult(dest, alpha1, A1, v1);
	VecScaleAdd(dest, 1.0, dest, beta1, w1, beta2, w2);
}



} // namespace ug
#endif /* __H__UG__MARTIN_ALGEBRA__OPERATIONS_MAT__ */
