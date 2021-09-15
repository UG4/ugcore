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

#include <limits>
#include "common/error.h"
#include "common/util/smart_pointer.h"
#include "lib_algebra/operator/interface/preconditioner.h"

#ifdef UG_PARALLEL
	#include "pcl/pcl_util.h"
	#include "lib_algebra/parallelization/parallelization_util.h"
	#include "lib_algebra/parallelization/matrix_overlap.h"
	#include "lib_algebra/parallelization/overlap_writer.h"
#endif

#include "lib_algebra/ordering_strategies/algorithms/IOrderingAlgorithm.h"
#include "lib_algebra/ordering_strategies/algorithms/native_cuthill_mckee.h" // for backward compatibility

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
// Returns true on success, or false on issues that lead to some changes in the solution
// (the solution is computed unless no exceptions are thrown)
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
// Returns true on success, or false on issues that lead to some changes in the solution
// (the solution is computed unless no exceptions are thrown)
template<typename Matrix_type, typename Vector_type>
bool invert_U(const Matrix_type &A, Vector_type &x, const Vector_type &b,
			  const number eps = 1e-8)
{
	PROFILE_FUNC_GROUP("algebra ILU");
	typedef typename Matrix_type::const_row_iterator const_row_iterator;

	typename Vector_type::value_type s;

	bool result = true;

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
		//TODO: Note that this may happen for problems with naturally
		// non-zero kernels, e.g. for the Stokes equation. One should
		// probably suppress this message in those cases but set the
		// rhs to 0.
		if (BlockNorm(A(i,i)) <= eps * BlockNorm(s))
		{
			UG_LOG("ILU Warning: Near-zero last diagonal entry "
					"with norm "<<BlockNorm(A(i,i))<<" in U "
					"for non-near-zero rhs entry with norm "
					<< BlockNorm(s) << ". Setting rhs to zero.\n"
					"NOTE: Reduce 'eps' using e.g. ILU::set_inversion_eps(...) "
					"to avoid this warning. Current eps: " << eps << ".\n")
			// set correction to zero
			x[i] = 0;
			result = false;
		} else {
			// c[i] = s/uii;
			InverseMatMult(x[i], 1.0, A(i,i), s);
		}
	}
	if(x.size() <= 1) return result;

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

	return result;
}


#ifdef UG_PARALLEL
inline void
LayoutEntriesToEndPermutation(std::vector<size_t>& newIndexOut,
                              const IndexLayout& layout,
                       		  const size_t numRows)
{
	using namespace std;
	vector<size_t> inds;
	CollectElements(inds, layout);

	newIndexOut.clear();
	newIndexOut.resize(numRows, 0);

	size_t numUniqueInds = 0;
	for(size_t i = 0; i < inds.size(); ++i){
		if(newIndexOut[inds[i]] != 1){
			newIndexOut[inds[i]] = 1;
			++numUniqueInds;
		}
	}

	size_t curNormalInd = 0;
	size_t curSlaveInd = numRows - numUniqueInds;

//	0 marks normal index, 1 marks layout index
	for(size_t irow = 0; irow < numRows; ++irow){
		if(newIndexOut[irow] == 0)
			newIndexOut[irow] = curNormalInd++;
		else
			newIndexOut[irow] = curSlaveInd++;
	}
}
#endif


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

	///	Ordering container type
		typedef std::vector<size_t> ordering_container_type;

	///	Ordering algorithm type
		typedef IOrderingAlgorithm<matrix_type, ordering_container_type> ordering_algo_type;

	protected:
		using base_type::set_debug;
		using base_type::debug_writer;
		using base_type::write_debug;
		using base_type::print_debugger_message;

	public:
	//	Constructor
		ILU(double beta=0.0) :
			m_beta(beta),
			m_sortEps(1.e-50),
			m_invEps(1.e-8),
			m_bSort(false),
			m_bDisablePreprocessing(false),
			m_useConsistentInterfaces(false),
			m_useOverlap(false) {};

	/// clone constructor
		ILU( const ILU<TAlgebra> &parent )
			: base_type(parent),
			  m_beta(parent.m_beta),
			  m_sortEps(parent.m_sortEps),
			  m_invEps(parent.m_invEps),
			  m_bSort(parent.m_bSort),
			  m_bDisablePreprocessing(parent.m_bDisablePreprocessing),
			  m_useConsistentInterfaces(parent.m_useConsistentInterfaces),
			  m_useOverlap(parent.m_useOverlap)
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

	/// 	sets an ordering algorithm
		void set_ordering_algorithm(SmartPtr<ordering_algo_type> ordering_algo){
			m_spOrderingAlgo = ordering_algo;
		}

	/// set cuthill-mckee sort on/off
		void set_sort(bool b)
		{
			m_bSort = b;

/* obsolete?
			if(m_spOrderingAlgo.valid()){
				delete m_spOrderingAlgo;
			}
*/

			if(m_bSort){
				m_spOrderingAlgo = make_sp(new NativeCuthillMcKeeOrdering<matrix_type, ordering_container_type>());
			}

			std::cerr << "please use 'set_ordering_algorithm(..)' in the future" << std::endl;
		}

	/// disable preprocessing (if underlying matrix has not changed)
		void set_disable_preprocessing(bool bDisable)	{m_bDisablePreprocessing = bDisable;}

	///	sets the smallest allowed value for sorted factorization
		void set_sort_eps(number eps)					{m_sortEps = eps;}

	///	sets the smallest allowed value for the Aii/Bi quotient
		void set_inversion_eps(number eps)				{m_invEps = eps;}

	///	enables consistent interfaces.
	/**	Connections between coefficients which lie in the same parallel interface
	 * are made consistent between processes.*/
		void enable_consistent_interfaces (bool enable)	{m_useConsistentInterfaces = enable;}

		void enable_overlap (bool enable)				{m_useOverlap = enable;}

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
			#ifdef UG_PARALLEL
			write_overlap_debug(mat, "ILU_prep_01_A_BeforeMakeUnique");
			#else
			write_debug(mat, "ILU_PreProcess_orig_A");
			#endif

			m_ILU = mat;

		//	this is an experimental trigger. Best to leave it false for now.
			const bool sortSlaveToEnd = false;

			#ifdef 	UG_PARALLEL
				if(m_useOverlap){
					CreateOverlap(m_ILU);
					m_oD.set_layouts(m_ILU.layouts());
					m_oC.set_layouts(m_ILU.layouts());
					m_oD.resize(m_ILU.num_rows(), false);
					m_oC.resize(m_ILU.num_rows(), false);

					if(debug_writer().valid()){
						m_overlapWriter = make_sp(new OverlapWriter<TAlgebra>());
						m_overlapWriter->init (*m_ILU.layouts(),
					                     	   *debug_writer(),
				                      		   m_ILU.num_rows());
					}

				//	slave-overlap rows shall be at the end of the matrix
					if(sortSlaveToEnd){
						m_bSort = true;
						LayoutEntriesToEndPermutation(	m_newIndex,
						                              	m_ILU.layouts()->slave(),
						                       			m_ILU.num_rows());
						m_bSortIsIdentity = GetInversePermutation(m_newIndex, m_oldIndex);

						if(!m_bSortIsIdentity)
						{
							matrix_type mat;
							mat = m_ILU;
							SetMatrixAsPermutation(m_ILU, mat, m_newIndex);
						}
					}
				}
				else if(m_useConsistentInterfaces){
					MatMakeConsistentOverlap0(m_ILU);
				}
				else {
					MatAddSlaveRowsToMasterRowOverlap0(m_ILU);
				//	set dirichlet rows on slaves
					std::vector<IndexLayout::Element> vIndex;
					CollectUniqueElements(vIndex,  m_ILU.layouts()->slave());
					SetDirichletRow(m_ILU, vIndex);
				}
			#endif

			m_h.resize(m_ILU.num_cols());

			#ifdef UG_PARALLEL
			write_overlap_debug(m_ILU, "ILU_prep_02_A_AfterMakeUnique");
			#endif

		//	if using overlap we already sort in a different way
			if(m_spOrderingAlgo.valid() && !(m_useOverlap && sortSlaveToEnd)){
				m_spOrderingAlgo->set_matrix(&mat);
				m_spOrderingAlgo->compute();
				//m_spOrderingAlgo->check();
				m_ordering = m_spOrderingAlgo->ordering();
				SetMatrixAsPermutation(m_ILU, mat, m_ordering); //see SetMatrixAsPermutation

				m_bSortIsIdentity = GetInversePermutation(m_ordering, m_old_ordering);
			}

		//	Debug output of matrices
			#ifdef UG_PARALLEL
			write_overlap_debug(m_ILU, "ILU_prep_03_A_BeforeFactorize");
			#else
			write_debug(m_ILU, "ILU_PreProcess_U_BeforeFactor");
			#endif


		// 	Compute ILU Factorization
			if (m_beta!=0.0) FactorizeILUBeta(m_ILU, m_beta);
			else if(matrix_type::rows_sorted) FactorizeILUSorted(m_ILU, m_sortEps);
			else FactorizeILU(m_ILU);
			m_ILU.defragment();

		//	Debug output of matrices
			#ifdef UG_PARALLEL
			write_overlap_debug(m_ILU, "ILU_prep_04_A_AfterFactorize");
			#else
			write_debug(m_ILU, "ILU_PreProcess_U_AfterFactor");
			#endif

		//	we're done
			return true;
		}


		void applyLU(vector_type &c, const vector_type &d, vector_type &tmp)
		{	
			if(m_spOrderingAlgo.invalid() || m_bSortIsIdentity)
			{
				// 	apply iterator: c = LU^{-1}*d
				if(! invert_L(m_ILU, tmp, d)) // h := L^-1 d
					print_debugger_message("ILU: There were issues at inverting L\n");
				if(! invert_U(m_ILU, c, tmp, m_invEps)) // c := U^-1 h = (LU)^-1 d
					print_debugger_message("ILU: There were issues at inverting U\n");
			}
			else
			{
				// we save one vector here by renaming
				SetVectorAsPermutation(tmp, d, m_ordering);
				if(! invert_L(m_ILU, c, tmp)) // c = L^{-1} d
					print_debugger_message("ILU: There were issues at inverting L (after permutation)\n");
				if(! invert_U(m_ILU, tmp, c, m_invEps)) // tmp = (LU)^{-1} d
					print_debugger_message("ILU: There were issues at inverting U (after permutation)\n");
				SetVectorAsPermutation(c, tmp, m_old_ordering);
			}
		}

	//	Stepping routine
		virtual bool step(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp,
		                  vector_type& c,
		                  const vector_type& d)
		{
			PROFILE_BEGIN_GROUP(ILU_step, "algebra ILU");

			
		//	\todo: introduce damping
			#ifdef UG_PARALLEL
			//	for debug output (only for application is written)
				static bool first = true;

				if(first) write_overlap_debug(d, "ILU_step_1_d");
				if(m_useOverlap){
					for(size_t i = 0; i < d.size(); ++i)
						m_oD[i] = d[i];
					for(size_t i = d.size(); i < m_oD.size(); ++i)
						m_oD[i] = 0;
					m_oD.set_storage_type(PST_ADDITIVE);
					m_oD.change_storage_type(PST_CONSISTENT);

					if(first) write_overlap_debug(m_oD, "ILU_step_2_oD_consistent");

					applyLU(m_oC, m_oD, m_h);

					for(size_t i = 0; i < c.size(); ++i)
						c[i] = m_oC[i];
					SetLayoutValues(&c, c.layouts()->slave(), 0);
					c.set_storage_type(PST_UNIQUE);
				}
				else if(m_useConsistentInterfaces){
					// make defect consistent
					SmartPtr<vector_type> spDtmp = d.clone();
					spDtmp->change_storage_type(PST_CONSISTENT);
					applyLU(c, *spDtmp, m_h);

					// declare c unique to enforce that only master correction is used
					// when it is made consistent below
					c.set_storage_type(PST_UNIQUE);
				}
				else{
				//	make defect unique
					SmartPtr<vector_type> spDtmp = d.clone();
					spDtmp->change_storage_type(PST_UNIQUE);
					if(first) write_debug(*spDtmp, "ILU_step_2_d_unique");
					applyLU(c, *spDtmp, m_h);
					c.set_storage_type(PST_ADDITIVE);
				}

			//	write debug
				if(first) write_overlap_debug(c, "ILU_step_3_c");

				c.change_storage_type(PST_CONSISTENT);

			//	write debug
				if(first) {write_overlap_debug(c, "ILU_step_4_c_consistent"); first = false;}

			#else
				write_debug(d, "ILU_step_d");
				applyLU(c, d, m_h);
				write_debug(c, "ILU_step_c");
			#endif

		//	we're done
			return true;
		}

	///	Postprocess routine
		virtual bool postprocess() {return true;}

	private:
	#ifdef UG_PARALLEL
		template <class T> void write_overlap_debug(const T& t, std::string name)
		{
			if(debug_writer().valid()){
				if(m_useOverlap && m_overlapWriter.valid() && t.layouts()->overlap_enabled())
					m_overlapWriter->write(t, name);
				else
					write_debug(t, name.c_str());
			}
		}
	#endif

	protected:
	///	storage for factorization
		matrix_type m_ILU;

	///	help vector
		vector_type m_h;

	///	for overlaps only
		vector_type m_oD;
		vector_type m_oC;
		#ifdef UG_PARALLEL
			SmartPtr<OverlapWriter<TAlgebra> >	m_overlapWriter;
		#endif

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

	/// for ordering algorithms
		SmartPtr<ordering_algo_type> m_spOrderingAlgo;
		ordering_container_type m_ordering, m_old_ordering;

	/// whether or not to disable preprocessing
		bool m_bDisablePreprocessing;

		bool m_useConsistentInterfaces;
		bool m_useOverlap;
};

} // end namespace ug

#endif
