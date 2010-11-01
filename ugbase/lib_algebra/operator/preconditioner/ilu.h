/*
 * ilu.h
 *
 *  Created on: 04.07.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__ILU__
#define __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__ILU__

#include "lib_algebra/lib_algebra.h"
#include "lib_algebra/operator/operator_interface.h"
#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallelization.h"
#endif

namespace ug{

template<typename Matrix_type>
bool FactorizeILU(Matrix_type &A)
{
	typedef typename Matrix_type::value_type block_type;
	// for all rows
	for(size_t i=1; i < A.num_rows(); i++)
	{
		// eliminate all entries A(i, k) with k<i with rows A(k, .) and k<i
		for(typename Matrix_type::rowIterator it_k = A.beginRow(i); !it_k.isEnd() && ((*it_k).iIndex < i); ++it_k)
		{
			const size_t k = (*it_k).iIndex;
			block_type &a_ik = (*it_k).dValue;
			if(BlockNorm(a_ik) < 1e-7)	continue;

			// add row k to row i by A(i, .) -= A(k,.)  A(i,k) / A(k,k)
			// so that A(i,k) is zero.
			// safe A(i,k)/A(k,k) in A(i,k)
			a_ik /= A(k,k);

			typename Matrix_type::rowIterator it_j = it_k;
			for(++it_j; !it_j.isEnd(); ++it_j)
			{
				const size_t j = (*it_j).iIndex;
				block_type& a_ij = (*it_j).dValue;
				bool bFound;
				typename Matrix_type::rowIterator p = A.get_connection(k,j, bFound);
				if(bFound)
				{
					const block_type a_kj = (*p).dValue;
					a_ij -= a_kj * a_ik;
				}
			}
		}
	}

	return true;
}

template<typename Matrix_type>
bool FactorizeILUSorted(Matrix_type &A)
{
	typedef typename Matrix_type::value_type block_type;

	// for all rows
	for(size_t i=1; i < A.num_rows(); i++)
	{

		// eliminate all entries A(i, k) with k<i with rows A(k, .) and k<i
		for(typename Matrix_type::rowIterator it_k = A.beginRow(i); !it_k.isEnd() && ((*it_k).iIndex < i); ++it_k)
		{
			const size_t k = (*it_k).iIndex;
			block_type &a_ik = (*it_k).dValue;
			if(BlockNorm(a_ik) < 1e-7)	continue;
			block_type &a_kk = A(k,k);

			// add row k to row i by A(i, .) -= A(k,.)  A(i,k) / A(k,k)
			// so that A(i,k) is zero.
			// safe A(i,k)/A(k,k) in A(i,k)
			a_ik /= a_kk;

			typename Matrix_type::rowIterator it_ij = it_k; // of row i
			++it_ij; // skip a_ik
			typename Matrix_type::rowIterator it_kj = A.beginRow(k); // of row k

			while(!it_ij.isEnd() && !it_kj.isEnd())
			{
				if((*it_ij).iIndex > (*it_kj).iIndex)
					++it_kj;
				else if((*it_ij).iIndex < (*it_kj).iIndex)
					++it_ij;
				else
				{
					block_type &a_ij = (*it_ij).dValue;
					const block_type &a_kj = (*it_kj).dValue;
					a_ij -= a_kj * a_ik;
					++it_kj; ++it_ij;
				}
			}
		}
	}

	return true;
}


// solve x = L^-1 b
template<typename Matrix_type, typename Vector_type>
bool invert_L(const Matrix_type &A, Vector_type &x, const Vector_type &b)
{
	typename Vector_type::value_type s;
	for(size_t i=0; i < x.size(); i++)
	{
		s = b[i];
		for(typename Matrix_type::cRowIterator it = A.beginRow(i); !it.isEnd(); ++it)
		{
			if((*it).iIndex >= i) continue;
			MatMultAdd(s, 1.0, s, -1.0, (*it).dValue, x[(*it).iIndex]);
		}
		x[i] = s;
	}

	return true;
}

// solve x = U^-1 * b
template<typename Matrix_type, typename Vector_type>
bool invert_U(const Matrix_type &A, Vector_type &x, const Vector_type &b)
{
	typename Vector_type::value_type s;
	for(size_t i = x.size()-1; ; --i)
	{
		s = b[i];
		for(typename Matrix_type::cRowIterator it = A.beginRow(i); !it.isEnd(); ++it)
		{
			if((*it).iIndex <= i) continue;
			// s -= (*it).dValue * x[(*it).iIndex];
			MatMultAdd(s, 1.0, s, -1.0, (*it).dValue, x[(*it).iIndex]);

		}
		// x[i] = s/A(i,i);
		InverseMatMult(x[i], 1.0, A(i,i), s);
		if(i == 0) break;
	}

	return true;
}


template <typename TAlgebra>
class ILUPreconditioner : public IPreconditioner<TAlgebra>
{
	public:
	//	Algebra type
		typedef TAlgebra algebra_type;

	//	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	//	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

	public:
	//	Constructor
		ILUPreconditioner() {};

	// 	Clone
		ILinearIterator<vector_type,vector_type>* clone()
		{
			return new ILUPreconditioner<algebra_type>();
		}

	//	Destructor
		~ILUPreconditioner()
		{
			m_ILU.destroy();
			m_h.destroy();
		}

	protected:
	//	Name of preconditioner
		virtual const char* name() const {return "ILUPreconditioner";}

	//	Preprocess routine
		virtual bool preprocess(matrix_type& mat)
		{
		//	TODO: error handling / memory check

		//  Rename Matrix for convenience

#if 1
			matrix_type& A = *this->m_pMatrix;
			m_ILU = A;
			m_h.resize(A.num_cols());
#else
			matrix_type& A = *this->m_pMatrix;

		//	Resize Matrix
			if(	m_ILU.num_rows() != mat.num_rows() ||
				m_ILU.num_cols() != mat.num_cols())
			{
			//	destroy memory
				m_ILU.destroy();
				m_h.destroy();

			//	create new memory
				m_ILU.create(A.num_rows(), A.num_cols());

			//	Create help vector
				m_h.create(A.num_cols());
			}

		// 	Copy matrix
			for(size_t i=0; i < A.num_rows(); i++)
			{
				for(typename matrix_type::rowIterator it_k = A.beginRow(i); !it_k.isEnd(); ++it_k)
				{
					const size_t k = (*it_k).iIndex;
					m_ILU(i,k) = (*it_k).dValue;
				}
			}
#endif
		// 	Compute ILU Factorization
		if(matrix_type::rows_sorted)
			FactorizeILUSorted(m_ILU);
		else
			FactorizeILU(m_ILU);

			return true;
		}

	//	Stepping routine
		virtual bool step(matrix_type& mat, vector_type& c, const vector_type& d)
		{
			// apply iterator: c = LU^{-1}*d (damp is not used)
			invert_L(m_ILU, m_h, d); // h := L^-1 d
			invert_U(m_ILU, c, m_h); // c := U^-1 h = (LU)^-1 d

#ifdef 	UG_PARALLEL
		//	Correction is always consistent
		//	todo: We set here correction to consistent, but it is not. Think about how to use ilu in parallel.
			c.set_storage_type(PST_CONSISTENT);
#endif

			return true;
		}

	//	Postprocess routine
		virtual bool postprocess() {return true;}

	protected:
		matrix_type m_ILU;
		vector_type m_h;
};


} // end namespace ug

#endif
