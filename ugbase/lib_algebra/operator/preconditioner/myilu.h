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


#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__myILU__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__myILU__

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


//#include "MatrixView.h"


// ordering stuff
#include "lib_algebra/ordering_strategies/algorithms/IOrderingAlgorithm.h"
#include "lib_algebra/ordering_strategies/execution/matrix_ordering.h"
#include "lib_algebra/ordering_strategies/typedefs.h"

#include "permutation_util.h"

#include "common/code_marker.h"

namespace ug{


// (cf. Y Saad, Iterative methods for Sparse Linear Systems, p. 270)
// fill-in is computed and lumped onto the diagonal if beta != 0
template<typename Matrix_type>
void Factorize(Matrix_type &A, number beta, number eps = 1e-50){
	PROFILE_FUNC_GROUP("algebra myILU");
	typedef typename Matrix_type::row_iterator row_iterator;
	typedef typename Matrix_type::value_type block_type;

	STATIC_ASSERT(Matrix_type::rows_sorted, Matrix_has_to_have_sorted_rows);

	bool bSorted = Matrix_type::rows_sorted;
	if(bSorted){
		eps = 1e-15;
		std::cout << "[Factorize] sorted" << std::endl;
	}
	else{
		std::cout << "[Factorize] unsorted" << std::endl;
	}

	bool bFill = beta != 0;

	if(bFill){
		eps = 1e-50;
		std::cout << "[Factorize] fill-in" << std::endl;
	}

	// for all rows
	for(size_t i=1; i < A.num_rows(); i++){
		block_type Nii(A(i,i)); Nii*=0.0;

		// eliminate all entries A(i, k) with k<i with rows A(k, .) and k<i
		const row_iterator rowEnd = A.end_row(i);

		for(row_iterator it_k = A.begin_row(i); it_k != rowEnd && (it_k.index() < i); ++it_k){
			const size_t k = it_k.index();
			block_type &a_ik = it_k.value();
			block_type &a_kk = A(k,k);

			if(fabs(BlockNorm(A(k,k))) < eps * BlockNorm(A(i,k))){
				UG_THROW("myILU: Blocknorm of diagonal is near-zero for k="<<k<<
				         " with eps: "<< eps <<", ||A_kk||="<<fabs(BlockNorm(A(k,k)))
				         <<", ||A_ik||="<<BlockNorm(A(i,k)));
			}

			// 1) Contribution to L part:
			// store A(i,k)/A(k,k) in A(i,k)
			try{ a_ik /= a_kk; }
			UG_CATCH_THROW("Failed to calculate A_ik /= A_kk "
				"with i = " << i << " and k = " << k << ".");

			row_iterator it_ij = it_k; // of row i
			++it_ij; // skip a_ik

			// 2) Contribution to U part:
			// compute contributions from row k for j=k:N
			// add row k to row i by A(i, .) -= A(k,.)  A(i,k) / A(k,k)
			// so that A(i,k) is zero.
			if(bSorted){
				row_iterator it_kj = A.begin_row(k); // of row k

				while(it_ij != rowEnd && it_kj != A.end_row(k)){
					if(it_ij.index() > it_kj.index()){
						++it_kj;
					}
					else if(it_ij.index() < it_kj.index()){
						++it_ij;
					}
					else{
						block_type &a_ij = it_ij.value();
						const block_type &a_kj = it_kj.value();
						a_ij -= a_ik * a_kj;
						++it_kj; ++it_ij;
					}
				}
			}
			else{
				error();
			}
		}
		// add fill-in to diagonal
		if(bFill){
			AddMult(A(i,i), beta, Nii);
		}
	}
}


// solve x = L^-1 b
// Returns true on success, or false on issues that lead to some changes in the solution
// (the solution is computed unless no exceptions are thrown)
template<typename Matrix_type, typename Vector_type>
bool invert_L_(const Matrix_type &A, Vector_type &x, const Vector_type &b){
	PROFILE_FUNC_GROUP("algebra myILU");
	typedef typename Matrix_type::const_row_iterator const_row_iterator;

	typename Vector_type::value_type s;
	for(size_t i=0; i < x.size(); i++){
		s = b[i];
		for(const_row_iterator it = A.begin_row(i); it != A.end_row(i); ++it){
			if(it.index() >= i){
				continue;
			}
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
bool invert_U_(const Matrix_type &A, Vector_type &x, const Vector_type &b,
			  const number eps = 1e-8)
{
	PROFILE_FUNC_GROUP("algebra myILU");
	typedef typename Matrix_type::const_row_iterator const_row_iterator;

	typename Vector_type::value_type s;
	
	bool result = true;
	
	// last row diagonal U entry might be close to zero with corresponding close to zero rhs
	// when solving Navier Stokes system, therefore handle separately
	if(x.size() > 0){
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
		if(BlockNorm(A(i,i)) <= eps * BlockNorm(s)){
			UG_LOG("myILU Warning: Near-zero last diagonal entry "
					"with norm "<<BlockNorm(A(i,i))<<" in U "
					"for non-near-zero rhs entry with norm "
					<< BlockNorm(s) << ". Setting rhs to zero.\n"
					"NOTE: Reduce 'eps' using e.g. myILU::set_inversion_eps(...) "
					"to avoid this warning. Current eps: " << eps << ".\n")
			// set correction to zero
			x[i] = 0;
			result = false;
		}
		else{
			// c[i] = s/uii;
			InverseMatMult(x[i], 1.0, A(i,i), s);
		}
	}
	if(x.size() <= 1){
		return result;
	}

	// handle all other rows
	for(size_t i = x.size()-2; i > 0; --i){
		s = b[i];
		for(const_row_iterator it = A.begin_row(i); it != A.end_row(i); ++it){
			if(it.index() <= i){
				continue;
			}
			// s -= it.value() * x[it.index()];
			MatMultAdd(s, 1.0, s, -1.0, it.value(), x[it.index()]);

		}
		// x[i] = s/A(i,i);
		InverseMatMult(x[i], 1.0, A(i,i), s);
	}

	return result;
}


///	myILU / myILU(beta) preconditioner
template <typename TAlgebra>
class myILU : public IPreconditioner<TAlgebra>{
public:
	typedef typename TAlgebra::vector_type vector_type;
	typedef typename TAlgebra::matrix_type matrix_type;
	typedef IPreconditioner<TAlgebra> base_type;
	typedef typename base_type::matrix_operator_type matrix_operator_type;

	typedef std::vector<size_t> ordering_container_type;
	typedef IOrderingAlgorithm<matrix_type, weighted_base_graph_type, ordering_container_type> weighted_ordering_algo_type;
	typedef IOrderingAlgorithm<matrix_type, base_graph_type, ordering_container_type> unweighted_ordering_algo_type;

protected:
	using base_type::set_debug;
	using base_type::debug_writer;
	using base_type::write_debug;
	using base_type::print_debugger_message;

public:
//	Constructor
	myILU(double beta=0.0) :
		m_beta(beta),
		m_invEps(1.e-8),
		m_bDisablePreprocessing(false)/*,
		m_useConsistentInterfaces(false),
		m_useOverlap(false)*/{};

/// clone constructor
	myILU( const myILU<TAlgebra> &parent )
		: base_type(parent),
		  m_beta(parent.m_beta),
		  m_invEps(parent.m_invEps),
		  m_bDisablePreprocessing(parent.m_bDisablePreprocessing)/*,
		  m_useConsistentInterfaces(parent.m_useConsistentInterfaces),
		  m_useOverlap(parent.m_useOverlap)*/{}

///	Clone
	virtual SmartPtr<ILinearIterator<vector_type> > clone(){
		return make_sp(new myILU<TAlgebra>(*this));
	}

///	Destructor
	virtual ~myILU(){}

///	returns if parallel solving is supported
	virtual bool supports_parallel() const{ return false; }

///	set factor for Factorize
	void set_beta(double beta){ m_beta = beta; }

/// 	set ordering on/off
	void set_ordering_algorithm_weighted(SmartPtr<weighted_ordering_algo_type> ordering_algo){
		m_spOrderingAlgo_weighted = ordering_algo;
	}

/// 	set ordering on/off
	void set_ordering_algorithm_unweighted(SmartPtr<unweighted_ordering_algo_type> ordering_algo){
		m_spOrderingAlgo_unweighted = ordering_algo;
	}

/// 	disable preprocessing (if underlying matrix has not changed)
	void set_disable_preprocessing(bool bDisable){ m_bDisablePreprocessing = bDisable; }

///	sets the smallest allowed value for the Aii/Bi quotient
	void set_inversion_eps(number eps){ m_invEps = eps; }


/*

///	enables consistent interfaces.
///	Connections between coefficients which are located in the same parallel interface
///	 are made consistent between processes.
	void enable_consistent_interfaces(bool enable){ m_useConsistentInterfaces = enable; }

	void enable_overlap(bool enable){ m_useOverlap = enable; }
*/

protected:
//	Name of preconditioner
	virtual const char* name() const {return "myILU";}

protected:
//	Preprocess routine
	virtual bool preprocess(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp){
		// do not do a thing if preprocessing disabled
		if(m_bDisablePreprocessing){
			return true;
		}

		PROFILE_BEGIN_GROUP(myILU_preprocess, "algebra myILU");

		matrix_type &orig_matrix = *pOp;

		m_h.resize(orig_matrix.num_cols());
		
		if(m_spOrderingAlgo_weighted.invalid() && m_spOrderingAlgo_unweighted.invalid())
		{
			std::cout << "no ordering" << std::endl;
			m_matrix = orig_matrix;
		}
		else
		{
			if(!m_spOrderingAlgo_weighted.invalid()){
				std::cout << "weighted ordering" << std::endl;

				std::cout << "	m_spOrderingAlgo_weighted init;" << std::endl;

				if(m_spOrderingAlgo_weighted->type() == weighted_ordering_algo_type::Type::GRAPH_BASED){
					std::cout << "graph based" << std::endl;
					weighted_graph_type<TAlgebra> meta_graph(orig_matrix);
					g = meta_graph.graph();
					m_spOrderingAlgo_weighted->set_graph(g);
				}
				else if(m_spOrderingAlgo_weighted->type() == weighted_ordering_algo_type::Type::MATRIX_BASED){
					std::cout << "matrix based" << std::endl;
					m_spOrderingAlgo_weighted->set_matrix(&orig_matrix);
				}

				std::cout << "	m_spOrderingAlgo_weighted->compute();" << std::endl;
				m_spOrderingAlgo_weighted->compute();
				m_spOrderingAlgo_weighted->check();

				m_ordering = m_spOrderingAlgo_weighted->ordering();

				std::cout << "	reorder.execute();" << std::endl;
				reorder_matrix(m_matrix, orig_matrix, m_ordering);

				bool is_identity = inverse_permutation(m_ordering, m_old_ordering);

				if(is_identity){
					std::cout << "identity permutation" << std::endl;
				}
			}

			else if(!m_spOrderingAlgo_unweighted.invalid()){
				std::cout << "unweighted ordering" << std::endl;

				std::cout << "	m_spOrderingAlgo_unweighted init;" << std::endl;

				if(m_spOrderingAlgo_unweighted->type() == unweighted_ordering_algo_type::Type::GRAPH_BASED){
					std::cout << "graph based" << std::endl;
					unweighted_graph_type<TAlgebra> meta_graph(orig_matrix);
					g_unweighted = meta_graph.graph();
					m_spOrderingAlgo_unweighted->set_graph(g_unweighted);
				}
				else if(m_spOrderingAlgo_unweighted->type() == unweighted_ordering_algo_type::Type::MATRIX_BASED){
					std::cout << "matrix based" << std::endl;
					m_spOrderingAlgo_unweighted->set_matrix(&orig_matrix);
				}

				std::cout << "	m_spOrderingAlgo_weighted->compute();" << std::endl;
				m_spOrderingAlgo_unweighted->compute();
				m_spOrderingAlgo_unweighted->check();

				m_ordering = m_spOrderingAlgo_unweighted->ordering();

				std::cout << "	reorder.execute();" << std::endl;
				reorder_matrix(m_matrix, orig_matrix, m_ordering);

				bool is_identity = inverse_permutation(m_ordering, m_old_ordering);

				if(is_identity){
					std::cout << "identity permutation" << std::endl;
				}
			}
		}


	// 	Compute Factorization

		Factorize(m_matrix, m_beta);

		m_matrix.defragment(); //what is this?!

		return true;
	}


	void applyLU(vector_type &c, const vector_type &d, vector_type &tmp)
	{
		if(m_spOrderingAlgo_weighted.invalid() && m_spOrderingAlgo_unweighted.invalid())			////////////
		{
			// 	apply iterator: c = LU^{-1}*d
			if(! invert_L_(m_matrix, tmp, d)) // h := L^-1 d
				print_debugger_message("myILU: There were issues at inverting L\n");
			if(! invert_U_(m_matrix, c, tmp, m_invEps)) // c := U^-1 h = (LU)^-1 d
				print_debugger_message("myILU: There were issues at inverting U\n");
		}
		else
		{
			// we save one vector here by renaming
			permute_vector(tmp, d, m_ordering);

			if(! invert_L_(m_matrix, c, tmp)) // c = L^{-1} d
				print_debugger_message("myILU: There were issues at inverting L (after permutation)\n");
			if(! invert_U_(m_matrix, tmp, c, m_invEps)) // tmp = (LU)^{-1} d
				print_debugger_message("myILU: There were issues at inverting U (after permutation)\n");

			permute_vector(c, tmp, m_old_ordering);
		}
	}

//	Stepping routine
	//TODO: introduce damping
	virtual bool step(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp, vector_type& c, const vector_type& d)
	{
		PROFILE_BEGIN_GROUP(myILU_step, "algebra myILU");

		write_debug(d, "myILU_step_d");
		applyLU(c, d, m_h);
		write_debug(c, "myILU_step_c");

		return true;
	}

///	Postprocess routine
	virtual bool postprocess(){ return true; }

protected:
	matrix_type m_matrix;

///	help vector
	vector_type m_h;

/// Parameter for fill-In
	number m_beta;

/// smallest allowed value for the Aii/Bi quotient
	number m_invEps;

/// whether or not to disable preprocessing
	bool m_bDisablePreprocessing;


/// ordering stuff
	weighted_base_graph_type* g;
	unweighted_base_graph_type* g_unweighted;							//////////////////
	SmartPtr<weighted_ordering_algo_type> m_spOrderingAlgo_weighted;
	SmartPtr<unweighted_ordering_algo_type> m_spOrderingAlgo_unweighted;
	
	ordering_container_type m_ordering, m_old_ordering;

#if 0
	bool m_useConsistentInterfaces;
	bool m_useOverlap;

///	for overlaps only
	vector_type m_oD;
	vector_type m_oC;
#endif

#if 0
#ifdef UG_PARALLEL
	SmartPtr<OverlapWriter<TAlgebra> > m_overlapWriter;
#endif
#endif
};

} //namespace

#endif //guard
