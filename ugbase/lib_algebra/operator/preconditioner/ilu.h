/*
 * ilu.h
 *
 *  Created on: 04.07.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__ILU__
#define __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__ILU__

#include "lib_algebra/operator/operator_interface.h"
#ifdef UG_PARALLEL
	#include "pcl/pcl_util.h"
	#include "lib_algebra/cpu_algebra/sparsematrix_util.h"
	#include "lib_algebra/parallelization/parallelization_util.h"
#endif

namespace ug{


// ILU(0) solver, i.e. static pattern ILU w/ P=P(A)
// (cf. Y Saad, Iterative methods for Sparse Linear Systems, p. 270)
template<typename Matrix_type>
bool FactorizeILU(Matrix_type &A)
{
	typedef typename Matrix_type::value_type block_type;
	// for all rows
	for(size_t i=1; i < A.num_rows(); i++)
	{
		// eliminate all entries A(i, k) with k<i with rows A(k, .) and k<i
		const typename Matrix_type::row_iterator rowEnd = A.end_row(i);
		for(typename Matrix_type::row_iterator it_k = A.begin_row(i); it_k != rowEnd && (it_k.index() < i); ++it_k)
		{
			const size_t k = it_k.index();
			block_type &a_ik = it_k.value();
			if(BlockNorm(a_ik) < 1e-7)	continue; //Nae: this may be dangerous

			// add row k to row i by A(i, .) -= A(k,.)  A(i,k) / A(k,k)
			// so that A(i,k) is zero.
			// safe A(i,k)/A(k,k) in A(i,k)
			a_ik /= A(k,k);

			typename Matrix_type::row_iterator it_j = it_k;
			for(++it_j; it_j != rowEnd; ++it_j)
			{
				const size_t j = it_j.index();
				block_type& a_ij = it_j.value();
				bool bFound;
				typename Matrix_type::row_iterator p = A.get_connection(k,j, bFound);
				if(bFound)
				{
					const block_type a_kj = p.value();
					a_ij -= a_ik *a_kj;
				}
			}
		}
	}

	return true;
}

// ILU(0)-beta solver, i.e.
// -> static pattern ILU w/ P=P(A)
// -> Fill-in is computed and lumped onto the diagonal
template<typename Matrix_type>
bool FactorizeILUBeta(Matrix_type &A, number beta)
{
	//assert(0);
	typedef typename Matrix_type::value_type block_type;

	// for all rows i=1:n do
	for(size_t i=1; i < A.num_rows(); i++)
	{
		block_type &Aii = A(i,i);
		block_type Nii(Aii); Nii*=0.0;

		// for k=1:(i-1) do
		// eliminate all entries A(i, k) with k<i with rows A(k, .) and k<i
		const typename Matrix_type::row_iterator it_iEnd = A.end_row(i);
		for(typename Matrix_type::row_iterator it_ik = A.begin_row(i); it_ik != it_iEnd && (it_ik.index() < i); ++it_ik)
		{

			// add row k to row i by A(i, .) -=  [A(i,k) / A(k,k)] A(k,.)
			// such that A(i,k) is zero.

			// 1) Contribution to L part:
			// store A(i,k)/A(k,k) in A(i,k)
			const size_t k = it_ik.index();
			block_type &a_ik = it_ik.value();
			a_ik /= A(k,k);

			if(BlockNorm(a_ik) < 1e-7)	continue; //Nae: this may be dangerous

			// 2) Contribution to U part:
			// compute contributions from row k for j=k:N
			const typename Matrix_type::row_iterator it_kEnd = A.end_row(k);
			for (typename Matrix_type::row_iterator it_kj=A.begin_row(k); it_kj != it_kEnd ;++it_kj)
			{
				const size_t j = it_kj.index();
				if (j<i) continue;  // index j belongs L part?

				// this index j belongs U part
				const block_type& a_kj = it_kj.value();

				bool aijFound;
				typename Matrix_type::row_iterator pij = A.get_connection(i,j, aijFound);
				if(aijFound) {
					// entry belongs to pattern
					// -> proceed with standard elimination
					block_type a_ij = pij.value();
					a_ij -= a_ik *a_kj ;

				} else {
					// entry DOES NOT belong to pattern
					// -> we lump it onto the diagonal
					// TODO : non square matrices!!!
					Nii -=  a_ik * a_kj;
				}

			}
		}

		// add fill-in to diagonal
		AddMult(Aii, beta, Nii);
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
		for(typename Matrix_type::row_iterator it_k = A.begin_row(i); it_k != A.end_row(i) && (it_k.index() < i); ++it_k)
		{
			const size_t k = it_k.index();
			block_type &a_ik = it_k.value();
			if(BlockNorm(a_ik) < 1e-7)	continue;
			block_type &a_kk = A(k,k);

			// add row k to row i by A(i, .) -= A(k,.)  A(i,k) / A(k,k)
			// so that A(i,k) is zero.
			// safe A(i,k)/A(k,k) in A(i,k)
			a_ik /= a_kk;

			typename Matrix_type::row_iterator it_ij = it_k; // of row i
			++it_ij; // skip a_ik
			typename Matrix_type::row_iterator it_kj = A.begin_row(k); // of row k

			while(it_ij != A.end_row(i) && it_kj != A.end_row(k))
			{
				if(it_ij.index() > it_kj.index())
					++it_kj;
				else if(it_ij.index() < it_kj.index())
					++it_ij;
				else
				{
					block_type &a_ij = it_ij.value();
					const block_type &a_kj = it_kj.value();
					a_ij -= a_ik * a_kj;
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
		for(typename Matrix_type::const_row_iterator it = A.begin_row(i); it != A.end_row(i); ++it)
		{
			if(it.index() >= i) continue;
			MatMultAdd(s, 1.0, s, -1.0, it.value(), x[it.index()]);
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
		for(typename Matrix_type::const_row_iterator it = A.begin_row(i); it != A.end_row(i); ++it)
		{
			if(it.index() <= i) continue;
			// s -= it.value() * x[it.index()];
			MatMultAdd(s, 1.0, s, -1.0, it.value(), x[it.index()]);

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

	///	Matrix Operator type
		typedef typename IPreconditioner<TAlgebra>::matrix_operator_type matrix_operator_type;

	public:
	//	Constructor
		ILUPreconditioner(double beta=0.0) : m_pDebugWriter(NULL), m_beta(beta) {};


	// 		Clone
		ILinearIterator<vector_type,vector_type>* clone()
		{
			return new ILUPreconditioner<algebra_type>(m_beta);
		}

		//	Destructor
		virtual ~ILUPreconditioner()
		{
		}

		//	set debug output
		void set_debug(IDebugWriter<algebra_type>* debugWriter)
		{
			m_pDebugWriter = debugWriter;
		}

		// set factor for ILU_{\beta}
		void set_beta(double beta)
		{
			m_beta = beta;
		}

	protected:
	//	Name of preconditioner
		virtual const char* name() const {return "ILUPreconditioner";}

	//	Preprocess routine
		virtual bool preprocess(matrix_operator_type& mat)
		{
		//	TODO: error handling / memory check

		//  Rename Matrix for convenience
			matrix_type& A = mat;

#ifdef 	UG_PARALLEL
		//	copy original matrix
			MakeConsistent(A, m_ILU);

		//	set zero on slaves
			std::vector<IndexLayout::Element> vIndex;
			CollectUniqueElements(vIndex,  m_ILU.get_slave_layout());
			SetDirichletRow(m_ILU, vIndex);

		//	Debug output of matrices
			if(m_pDebugWriter != NULL)
			{
				m_pDebugWriter->write_matrix(m_ILU,
											 "ILUBeforeFactorize");
			}

#else
		//	copy original matrix
			m_ILU = A;
#endif

		//	resize help vector
			m_h.resize(A.num_cols());


		// 	Compute ILU Factorization
		if (m_beta!=0.0)
			FactorizeILUBeta(m_ILU, m_beta);

		else if(matrix_type::rows_sorted)
		{
			FactorizeILUSorted(m_ILU);
		}
		else
		{
			FactorizeILU(m_ILU);
		}

		//	Debug output of matrices
			if(m_pDebugWriter != NULL)
			{
				m_pDebugWriter->write_matrix(m_ILU,
											 "ILUAfterFactorize");
			}

		//	we're done
			return true;
		}

	//	Stepping routine
		virtual bool step(matrix_operator_type& mat, vector_type& c, const vector_type& d)
		{
		//	todo: Think about how to use ilu in parallel.

#ifdef UG_PARALLEL
			static bool first = true;

		//	make defect unique
			vector_type dhelp;
			dhelp.resize(d.size()); dhelp = d;
			dhelp.change_storage_type(PST_UNIQUE);

		// 	apply iterator: c = LU^{-1}*d (damp is not used)
			invert_L(m_ILU, m_h, dhelp); // h := L^-1 d

			// TODO: m_h does not have a storage type???
			//if(first) write_debug(m_h, "ILU_mh");

			invert_U(m_ILU, c, m_h); // c := U^-1 h = (LU)^-1 d

			if(first) write_debug(c, "ILU_c");

		//	Correction is always consistent
			c.set_storage_type(PST_ADDITIVE);
			c.change_storage_type(PST_CONSISTENT);

			if(first)
			{
				write_debug(c, "ILU_cConsistent");
				first = false;
			}

#else
		// 	apply iterator: c = LU^{-1}*d (damp is not used)
			invert_L(m_ILU, m_h, d); // h := L^-1 d
			invert_U(m_ILU, c, m_h); // c := U^-1 h = (LU)^-1 d
#endif

		//	we're done
			return true;
		}

	//	Postprocess routine
		virtual bool postprocess() {return true;}

	protected:
		bool write_debug(const vector_type& vec, const char* filename)
		{
		//	if no debug writer set, we're done
			if(m_pDebugWriter == NULL) return true;

		//	write
			return m_pDebugWriter->write_vector(vec, filename);
		}

	protected:
		matrix_type m_ILU;
		vector_type m_h;

		
	//	Debug Writer
		IDebugWriter<algebra_type>* m_pDebugWriter;
		
		// Factor for ILU-beta
		number m_beta;

};


} // end namespace ug

#endif
