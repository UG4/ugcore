/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Authors: Andreas Vogel, Arne Nägel
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */


#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__ILU__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__ILU__

#include "common/util/smart_pointer.h"
#include "lib_algebra/operator/interface/preconditioner.h"

#ifdef UG_PARALLEL
	#include "pcl/pcl_util.h"
	#include "lib_algebra/parallelization/parallelization_util.h"
	#include "lib_algebra/parallelization/parallel_matrix_overlap_impl.h"
#endif
#include "lib_algebra/algebra_common/permutation_util.h"

namespace ug{


// ILU(0) solver, i.e. static pattern ILU w/ P=P(A)
// (cf. Y Saad, Iterative methods for Sparse Linear Systems, p. 270)
template<typename Matrix_type>
bool FactorizeILU(Matrix_type &A)
{
	PROFILE_FUNC_GROUP("algebra ILU");
	typedef typename Matrix_type::row_iterator row_iterator;
	typedef typename Matrix_type::value_type block_type;

	// for all rows
	for(size_t i=1; i < A.num_rows(); i++)
	{
		// eliminate all entries A(i, k) with k<i with rows A(k, .) and k<i
		const row_iterator rowEnd = A.end_row(i);
		for(row_iterator it_k = A.begin_row(i);
								it_k != rowEnd && (it_k.index() < i); ++it_k)
		{
			const size_t k = it_k.index();
			block_type &a_ik = it_k.value();
		
			// add row k to row i by A(i, .) -= A(k,.)  A(i,k) / A(k,k)
			// so that A(i,k) is zero.
			// save A(i,k)/A(k,k) in A(i,k)
			if(fabs(BlockNorm(A(k,k))) < 1e-15*BlockNorm(A(i,k)))
				UG_THROW("Diag is Zero for k="<<k<<", cannot factorize ILU.");

			a_ik /= A(k,k);

			row_iterator it_j = it_k;
			for(++it_j; it_j != rowEnd; ++it_j)
			{
				const size_t j = it_j.index();
				block_type& a_ij = it_j.value();
				bool bFound;
				row_iterator p = A.get_connection(k,j, bFound);
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
	PROFILE_FUNC_GROUP("algebra ILUBeta");
	typedef typename Matrix_type::row_iterator row_iterator;
	typedef typename Matrix_type::value_type block_type;

	// for all rows i=1:n do
	for(size_t i=1; i < A.num_rows(); i++)
	{
		block_type &Aii = A(i,i);
		block_type Nii(Aii); Nii*=0.0;

		// for k=1:(i-1) do
		// eliminate all entries A(i, k) with k<i with rows A(k, .) and k<i
		const row_iterator it_iEnd = A.end_row(i);
		for(row_iterator it_ik = A.begin_row(i);
							it_ik != it_iEnd && (it_ik.index() < i); ++it_ik)
		{

			// add row k to row i by A(i, .) -=  [A(i,k) / A(k,k)] A(k,.)
			// such that A(i,k) is zero.

			// 1) Contribution to L part:
			// store A(i,k)/A(k,k) in A(i,k)
			const size_t k = it_ik.index();
			block_type &a_ik = it_ik.value();
			a_ik /= A(k,k);

			// 2) Contribution to U part:
			// compute contributions from row k for j=k:N
			const row_iterator it_kEnd = A.end_row(k);
			for (row_iterator it_kj=A.begin_row(k); it_kj != it_kEnd ;++it_kj)
			{
				const size_t j = it_kj.index();
				if (j<=k) continue;  // index j belongs L part?

				// this index j belongs U part
				const block_type& a_kj = it_kj.value();

				bool aijFound;
				row_iterator pij = A.get_connection(i,j, aijFound);
				if(aijFound) {
					// entry belongs to pattern
					// -> proceed with standard elimination
					block_type& a_ij = pij.value();
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
bool FactorizeILUSorted(Matrix_type &A, const number eps = 1e-50)
{
	PROFILE_FUNC_GROUP("algebra ILU");
	typedef typename Matrix_type::row_iterator row_iterator;
	typedef typename Matrix_type::value_type block_type;

	// for all rows
	for(size_t i=1; i < A.num_rows(); i++)
	{

		// eliminate all entries A(i, k) with k<i with rows A(k, .) and k<i
		for(row_iterator it_k = A.begin_row(i);
							it_k != A.end_row(i) && (it_k.index() < i); ++it_k)
		{
			const size_t k = it_k.index();
			block_type &a_ik = it_k.value();
			block_type &a_kk = A(k,k);

			// add row k to row i by A(i, .) -= A(k,.)  A(i,k) / A(k,k)
			// so that A(i,k) is zero.
			// safe A(i,k)/A(k,k) in A(i,k)
			if(fabs(BlockNorm(A(k,k))) < eps * BlockNorm(A(i,k)))
				UG_THROW("ILU: Blocknorm of diagonal is near-zero for k="<<k<<
				         " with eps: "<< eps <<", ||A_kk||="<<fabs(BlockNorm(A(k,k)))
				         <<", ||A_ik||="<<BlockNorm(A(i,k)));

			try {a_ik /= a_kk;}
			UG_CATCH_THROW("Failed to calculate A_ik /= A_kk "
				"with i = " << i << " and k = " << k << ".");


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
	PROFILE_FUNC_GROUP("algebra ILU");
	typedef typename Matrix_type::const_row_iterator const_row_iterator;

	typename Vector_type::value_type s;
	for(size_t i=0; i < x.size(); i++)
	{
		s = b[i];
		for(const_row_iterator it = A.begin_row(i); it != A.end_row(i); ++it)
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
bool invert_U(const Matrix_type &A, Vector_type &x, const Vector_type &b,
			  const number eps = 1e-8)
{
	PROFILE_FUNC_GROUP("algebra ILU");
	typedef typename Matrix_type::const_row_iterator const_row_iterator;

	typename Vector_type::value_type s;
	
	size_t numNearZero=0;

	// last row diagonal U entry might be close to zero with corresponding close to zero rhs
	// when solving Navier Stokes system, therefore handle separately
	if(x.size() > 0)
	{
		size_t i=x.size()-1;
		s = b[i];

		// check if diag part is significantly smaller than rhs
		// This may happen when matrix is indefinite with one eigenvalue
		// zero. In that case, the factorization on the last row is
		// nearly zero due to round-off errors. In order to allow ill-
		// scaled matrices (i.e. small matrix entries row-wise) this
		// is compared to the rhs, that is small in this case as well.
		if (BlockNorm(A(i,i)) <= eps * BlockNorm(s))
		{
			if(numNearZero++<5)
			{	UG_LOG("ILU Warning: Near-zero diagonal entry "
					"with norm "<<BlockNorm(A(i,i))<<" in last row of U "
					" with corresponding non-near-zero rhs with norm "
					<< BlockNorm(s) << ". Setting rhs to zero.\n");
				UG_LOG("NOTE: Call this method with a smaller 'eps' parameter "
					   "to avoid this warning. (current eps: " << eps <<
					   "). If this method is called from the "
					   "ILU preconditioner class, you may want to call "
					   "ILU::set_inversion_eps(...) with a smaller threshold.\n")
			}
			// set correction to zero
			x[i] = 0;
		} else {
			// c[i] = s/uii;
			InverseMatMult(x[i], 1.0, A(i,i), s);
		}
	}
	if(x.size() <= 1) return true;

	// handle all other rows
	for(size_t i = x.size()-2; ; --i)
	{
		s = b[i];
		for(const_row_iterator it = A.begin_row(i); it != A.end_row(i); ++it)
		{
			if(it.index() <= i) continue;
			// s -= it.value() * x[it.index()];
			MatMultAdd(s, 1.0, s, -1.0, it.value(), x[it.index()]);

		}
		// x[i] = s/A(i,i);
		InverseMatMult(x[i], 1.0, A(i,i), s);
		if(i == 0) break;
	}

	if(numNearZero>=5)
	{	UG_LOG("...\nILU Warning: " << numNearZero << " ( out of " << x.size() << ") near-zero diagonal entries in last row of U.\n");	}

	return true;
}

///	ILU / ILU(beta) preconditioner
template <typename TAlgebra>
class ILU : public IPreconditioner<TAlgebra>
{
	public:
	///	Algebra type
		typedef TAlgebra algebra_type;

	///	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	///	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

	///	Matrix Operator type
		typedef typename IPreconditioner<TAlgebra>::matrix_operator_type matrix_operator_type;

	///	Base type
		typedef IPreconditioner<TAlgebra> base_type;

	protected:
		using base_type::set_debug;
		using base_type::debug_writer;
		using base_type::write_debug;

	public:
	//	Constructor
		ILU(double beta=0.0) :
			m_beta(beta),
			m_sortEps(1.e-50),
			m_invEps(1.e-8),
			m_bSort(false),
			m_bDisablePreprocessing(false) {};

	/// clone constructor
		ILU( const ILU<TAlgebra> &parent )
			: base_type(parent),
			  m_beta(parent.m_beta),
			  m_sortEps(parent.m_sortEps),
			  m_invEps(parent.m_invEps),
			  m_bSort(parent.m_bSort),
			  m_bDisablePreprocessing(parent.m_bDisablePreprocessing)
		{	}

	///	Clone
		virtual SmartPtr<ILinearIterator<vector_type> > clone()
		{
			return make_sp(new ILU<algebra_type>(*this));
		}

	///	Destructor
		virtual ~ILU(){}

	///	returns if parallel solving is supported
		virtual bool supports_parallel() const {return true;}

	///	set factor for \f$ ILU_{\beta} \f$
		void set_beta(double beta) {m_beta = beta;}

	/// set cuthill-mckee sort on/off
		void set_sort(bool b)
		{
			m_bSort = b;
		}

	/// disable preprocessing (if underlying matrix has not changed)
		void set_disable_preprocessing(bool bDisable)	{m_bDisablePreprocessing = bDisable;}

	///	sets the smallest allowed value for sorted factorization
		void set_sort_eps(number eps)					{m_sortEps = eps;}

	///	sets the smallest allowed value for the Aii/Bi quotient
		void set_inversion_eps(number eps)				{m_invEps = eps;}

	protected:
	//	Name of preconditioner
		virtual const char* name() const {return "ILU";}

	protected:
		// cuthill-mckee sorting
		void calc_cuthill_mckee()
		{
			PROFILE_BEGIN_GROUP(ILU_ReorderCuthillMcKey, "ilu algebra");
			GetCuthillMcKeeOrder(m_ILU, m_newIndex);
			m_bSortIsIdentity = GetInversePermutation(m_newIndex, m_oldIndex);

			if(!m_bSortIsIdentity)
			{
				matrix_type mat;
				mat = m_ILU;
				SetMatrixAsPermutation(m_ILU, mat, m_newIndex);
			}
		}

	protected:

	//	Preprocess routine
		virtual bool preprocess(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp)
		{
			// do not do a thing if preprocessing disabled
			if (m_bDisablePreprocessing) return true;

			matrix_type &mat = *pOp;
			PROFILE_BEGIN_GROUP(ILU_preprocess, "algebra ILU");
		//	Debug output of matrices
			write_debug(mat, "ILU_prep_01_BeforeMakeUnique");

			m_ILU = mat;
#ifdef 	UG_PARALLEL
			MatAddSlaveRowsToMasterRowOverlap0(m_ILU);

		//	set zero on slaves
			std::vector<IndexLayout::Element> vIndex;
			CollectUniqueElements(vIndex,  m_ILU.layouts()->slave());
			SetDirichletRow(m_ILU, vIndex);

		//DEBUG: find dirichlet master rows
			// vIndex.clear();
			// CollectUniqueElements(vIndex,  m_ILU.layouts()->master());
			// typedef typename matrix_type::row_iterator row_iterator;
			// typedef typename matrix_type::value_type block_type;

			// // for all rows
			// UG_LOG("Dirichlet Masters: ");
			// for(size_t i_index=0; i_index < vIndex.size(); i_index++)
			// {
			// 	const size_t i = vIndex[i_index];

			// 	size_t numNonzero = 0;
			// 	size_t lastNonzero = 0;
			// 	for(row_iterator it_k = m_ILU.begin_row(i); it_k != m_ILU.end_row(i); ++it_k)
			// 	{
			// 		const size_t k = it_k.index();
			// 		block_type &a = it_k.value();

			// 	//todo:	Use a block-compatible comparison
			// 		if(a != 0){
			// 			++numNonzero;
			// 			lastNonzero = k;
			// 		}
			// 	}

			// 	if(numNonzero == 1){
			// 		m_ILU(i, lastNonzero) = 1;
			// 	}
			// }
			// UG_LOG("\n");
#endif

			write_debug(m_ILU, "ILU_prep_02_AfterMakeUnique");

			if(m_bSort)
				calc_cuthill_mckee();

		//	Debug output of matrices
			write_debug(m_ILU, "ILU_prep_03_BeforeFactorize");

		//	resize help vector
			m_h.resize(mat.num_cols());

		// 	Compute ILU Factorization
			if (m_beta!=0.0) FactorizeILUBeta(m_ILU, m_beta);
			else if(matrix_type::rows_sorted) FactorizeILUSorted(m_ILU, m_sortEps);
			else FactorizeILU(m_ILU);
			m_ILU.defragment();

		//	Debug output of matrices
			write_debug(m_ILU, "ILU_prep_04_AfterFactorize");

		//	we're done
			return true;
		}


		void applyLU(vector_type &c, const vector_type &d, vector_type &tmp)
		{
			if(!m_bSort || m_bSortIsIdentity)
			{
				// 	apply iterator: c = LU^{-1}*d
				invert_L(m_ILU, tmp, d); // h := L^-1 d
				invert_U(m_ILU, c, tmp, m_invEps); // c := U^-1 h = (LU)^-1 d
			}
			else
			{
				// we save one vector here by renaming
				SetVectorAsPermutation(tmp, d, m_newIndex);
				invert_L(m_ILU, c, tmp); // c = L^{-1} d
				invert_U(m_ILU, tmp, c, m_invEps); // tmp = (LU)^{-1} d
				SetVectorAsPermutation(c, tmp, m_oldIndex);
			}
		}

	//	Stepping routine
		virtual bool step(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp, vector_type& c, const vector_type& d)
		{
			PROFILE_BEGIN_GROUP(ILU_step, "algebra ILU");
		//	\todo: introduce damping
#ifdef UG_PARALLEL
		//	for debug output (only for application is written)
			static bool first = true;

		//	make defect unique
			SmartPtr<vector_type> spDtmp = d.clone();
			spDtmp->change_storage_type(PST_UNIQUE);

			applyLU(c, *spDtmp, m_h);

		//	Correction is always consistent
			c.set_storage_type(PST_ADDITIVE);

		//	write debug
			if(first) write_debug(c, "ILU_c");

			c.change_storage_type(PST_CONSISTENT);

		//	write debug
			if(first) {write_debug(c, "ILU_cConsistent"); first = false;}

#else
			applyLU(c, d, m_h);
#endif

		//	we're done
			return true;
		}

	///	Postprocess routine
		virtual bool postprocess() {return true;}

	protected:
	///	storage for factorization
		matrix_type m_ILU;

	///	help vector
		vector_type m_h;
		
	/// Factor for ILU-beta
		number m_beta;

	///	smallest allowed value for sorted factorization
		number m_sortEps;

	///	smallest allowed value for the Aii/Bi quotient
		number m_invEps;

	/// for cuthill-mckee reordering
		std::vector<size_t> m_newIndex, m_oldIndex;
		bool m_bSortIsIdentity;
		bool m_bSort;

	/// whether or not to disable preprocessing
		bool m_bDisablePreprocessing;
};

} // end namespace ug

#endif
