/*
 * schur_complement_operator.cpp
 *
 *  Created on: 18.12.2013
 *      Author: anaegel
 */

//	THIS FILE IS ONLY TEMPORARY!
//  When everything works as
//	expected one should move the code in this file to
//	schur_impl.h, so that it will work with all the different algebras.
//
//	In the moment, template instantiations are invoked at the end of this file.

#ifdef UG_PARALLEL


// extern headers
#include <cmath>
#include <sstream>  // added for 'stringstream'

// algebra types
#include "lib_algebra/cpu_algebra_types.h"
#include "lib_algebra/operator/algebra_debug_writer.h"
#include "lib_algebra/operator/preconditioner/preconditioners.h"
#include "lib_algebra/parallelization/parallel_index_layout.h"

#include "lib_algebra/algebra_common/sparsematrix_util.h"  // DenseMatrixFromSparseMatrix
#include "lib_algebra/small_algebra/small_algebra.h"   // DEnseMAtrix...

#include "pcl/pcl_layout_tests.h"

#include "common/profiler/profiler.h"        // additions for profiling

// own header
#include "schur.h"
#include "common/progress.h"

namespace ug{


template <class VT>
void UG_LOG_Vector(const VT &vec)
{
	for (size_t i=0; i<vec.size(); ++i)
	{ UG_LOG( vec[i] << " "); }
	UG_LOG(std::endl);
	/*{ std::cerr << vec[i] << " "; }
		std::cerr << std::endl;*/
}

template <class MT>
void UG_LOG_Matrix(const MT &A)
{
	UG_DEBUG_BEGIN(SchurDebug, 8)
	for (size_t i = 0; i<A.num_rows(); ++i)
	{
		UG_LOG ("{ "<<i<<": ");
		for(typename MT::const_row_iterator it = A.begin_row(i); it != A.end_row(i); ++it)
		{ UG_LOG ("("<< i << "," << it.index() <<" ) =" << A(i, it.index()) << std::endl); }
		UG_LOG ("}\n");
	}
	UG_DEBUG_END(SchurDebug, 8)
};

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	LocalSchurComplement implementation
template <typename TAlgebra>
void SchurComplementOperator<TAlgebra>::
init()
{
	try{
	SCHUR_PROFILE_BEGIN(SCHUR_Op_init);
	PROGRESS_START(prog, 0, "SchurComplementOperator::init");

	UG_DLOG(SchurDebug, 5, "\n% 'SchurComplementOperator::init()':");
//	check that operator has been set
	if(m_spOperator.invalid())
		UG_THROW("SchurComplementOperator::init: No Operator A set.");


//	save matrix from which we build the Schur complement
	typename TAlgebra::matrix_type &mat = m_spOperator->get_matrix();

	//	get matrix from dirichlet operator

	UG_DLOG(SchurDebug, 2, "SchurComplementOperator::m_spOperator.layouts():\n" << *mat.layouts() << "\n")
	// init sub matrices
	m_slicing.get_matrix(mat, SD_INNER, SD_INNER, sub_matrix(SD_INNER, SD_INNER));
	m_slicing.get_matrix(mat, SD_INNER, SD_SKELETON, sub_matrix(SD_INNER, SD_SKELETON));
	m_slicing.get_matrix(mat, SD_SKELETON, SD_INNER, sub_matrix(SD_SKELETON, SD_INNER));
	m_slicing.get_matrix(mat, SD_SKELETON, SD_SKELETON, sub_matrix(SD_SKELETON, SD_SKELETON));

	sub_matrix(SD_SKELETON, SD_SKELETON).set_layouts(m_slicing.get_slice_layouts(mat.layouts(), SD_SKELETON));
	sub_matrix(SD_SKELETON, SD_SKELETON).set_storage_type(PST_ADDITIVE);

	IF_DEBUG(SchurDebug, 5)
	{
		// debug screen
		UG_LOG("A             ="); UG_LOG_Matrix(mat); UG_LOG(std::endl);
		UG_LOG("A_I,I         ="); UG_LOG_Matrix(sub_matrix(SD_INNER, SD_INNER)); UG_LOG(std::endl);
		UG_LOG("A_I,Gamma     ="); UG_LOG_Matrix(sub_matrix(SD_INNER, SD_SKELETON)); UG_LOG(std::endl);
		UG_LOG("A_Gamma,I     ="); UG_LOG_Matrix(sub_matrix(SD_SKELETON, SD_INNER)); UG_LOG(std::endl);
		UG_LOG("A_Gamma,Gamma ="); UG_LOG_Matrix(sub_matrix(SD_SKELETON, SD_SKELETON)); UG_LOG(std::endl);

		//	debug output of matrices
		write_debug(sub_matrix(SD_INNER, SD_INNER), "Schur_II.mat");
		write_debug(sub_matrix(SD_INNER, SD_SKELETON), "Schur_IB.mat");
		write_debug(sub_matrix(SD_SKELETON, SD_INNER), "Schur_BI.mat");
		write_debug(sub_matrix(SD_SKELETON, SD_SKELETON), "Schur_BB.mat");
	}

	sub_matrix(SD_INNER, SD_INNER).set_storage_type(PST_ADDITIVE);
	sub_matrix(SD_INNER, SD_SKELETON).set_storage_type(PST_ADDITIVE);
	sub_matrix(SD_SKELETON, SD_INNER).set_storage_type(PST_ADDITIVE);
	sub_matrix(SD_SKELETON, SD_SKELETON).set_storage_type(PST_ADDITIVE);

	SCHUR_PROFILE_BEGIN(SchurComplementOperator_init_DirSolver);
	//	init solver for Dirichlet problem
		if(m_spDirichletSolver.valid())
			if(!m_spDirichletSolver->init(sub_operator(SD_INNER, SD_INNER)))
				UG_THROW("LocalSchurComplement::init: Cannot init "
						"Dirichlet solver for operator A.");
	SCHUR_PROFILE_END_(SchurComplementOperator_init_DirSolver);

//	reset apply counter
	m_applyCnt = 0;
	}UG_CATCH_THROW("SchurComplementOperator::" << __FUNCTION__ << " failed")

} /* end 'SchurComplementOperator::init()' */


/// apply schur complement f_{\Gamma} = (A_{\Gamma, \Gamma} -  A_{\Gamma, I}  A_{I, I}^{-1}  A_{I, \Gamma} )u_{\Gamma}
template <typename TAlgebra>
void SchurComplementOperator<TAlgebra>::
apply(vector_type& fskeleton, const vector_type& uskeleton)
{
	try{
	SCHUR_PROFILE_BEGIN(SCHUR_Op_apply);

	UG_DLOG(SchurDebug, 5, "\n% 'SchurComplementOperator::apply()':"<< std::endl);

	SCHUR_PROFILE_BEGIN(SCHUR_apply);
//	check that matrix has been set
	if(m_spOperator.invalid())
		UG_THROW("SchurComplementOperator::apply: Matrix A not set.");

//	check Dirichlet solver
	if(m_spDirichletSolver.invalid())
		UG_THROW("SchurComplementOperator::apply: No Dirichlet Solver set.");

//	Check parallel storage type of matrix
	if(!sub_matrix(SD_INNER, SD_INNER).has_storage_type(PST_ADDITIVE))
		UG_THROW("SchurComplementOperator::apply: Inadequate storage format of matrix.");

// DEBUG ONLY!!!
/*	 vector_type uskeleton(uskeleton2);
	 uskeleton.set(0.0);
	 if (pcl::GetProcRank()==0) uskeleton[2] = 1.0;
*/

//	Check parallel storage type of vectors
	if (!uskeleton.has_storage_type(PST_CONSISTENT))
		UG_THROW("SchurComplementOperator::apply: Inadequate storage format of vec 'uskeleton' (should be consistent).");


	//	aux vectors
	const int n_inner = sub_size(SD_INNER);
	vector_type uinner; uinner.create(n_inner);
	vector_type finner; finner.create(n_inner);



	// A. compute first contribution
	//    $A_{\Gamma, \Gamma} u_{Gamma}$
	fskeleton.set_storage_type(PST_ADDITIVE);
	IF_DEBUG(SchurDebug, 5) { UG_LOG("\nSchurU_Gamma"); UG_LOG_Vector(uskeleton);}

	sub_operator(SD_SKELETON, SD_SKELETON)->apply(fskeleton, uskeleton);
	IF_DEBUG(SchurDebug, 5) { UG_LOG("\nSchurF_Gamma1="); UG_LOG_Vector(fskeleton);}

	/*
	if(debug_writer().valid())
	{
		//	Debug output of vector
		//	add iter count to name
			std::string name("Schur ");
			char ext[20]; sprintf(ext, "_a_rhs_skeleton%03d.vec", m_applyCnt);
			name.append(ext);
			//UG_LOG("fskelton="<<fskeleton);
			//debug_writer()->write_vector(fskeleton, name.c_str());
	}
*/

	// B. compute second contribution
	// f_{\Gamma} -=  A_{\Gamma, I}  A_{I, I}^{-1}  A_{I, \Gamma} )u_{\Gamma}
    // Note: requires no communication in the interior!

	sub_operator(SD_INNER, SD_SKELETON)->apply(finner, uskeleton);
	// UG_LOG("\nSchurF_Inner="); UG_LOG_Vector<vector_type>(finner);


	// Storage types do not matter,
    // but are required for solver!
	finner.set_storage_type(PST_ADDITIVE);
	uinner.set_storage_type(PST_CONSISTENT);

	if(!m_spDirichletSolver->apply_return_defect(uinner, finner))
	{
		UG_LOG_ALL_PROCS("ERROR in 'SchurComplementOperator::apply': "
						 "Could not solve Dirichlet problem (step 3.b) on Proc "
							<< pcl::GetProcRank() << /*" (m_statType = '" << m_statType <<*/ "').\n");
		UG_LOG_ALL_PROCS("ERROR in 'SchurComplementOperator::apply':"
						" Last defect was " << m_spDirichletSolver->defect() <<
						" after " << m_spDirichletSolver->step() << " steps.\n");

		UG_THROW("Cannot solve Local Schur Complement.");
	}
	// UG_LOG("\nSchurU_Inner="); UG_LOG_Vector<vector_type>(uinner);

	// modify fskeleton
	sub_operator(SD_SKELETON, SD_INNER)->apply_sub(fskeleton, uinner);
	IF_DEBUG(SchurDebug, 5){UG_LOG("\nSchurF_Gamma2="); UG_LOG_Vector(fskeleton); UG_LOG(std::endl);}

	fskeleton.set_storage_type(PST_ADDITIVE);
	if(!fskeleton.has_storage_type(PST_ADDITIVE))
		UG_THROW("SchurComplementOperator::apply: Inadequate storage format of vec 'fskeleton' (should be additive).");

//	UG_THROW("STOPP!");
	m_applyCnt++;

	}UG_CATCH_THROW("SchurComplementOperator::" << __FUNCTION__ << " failed")
} /* end 'SchurComplementOperator::apply()' */

template <typename TAlgebra>
void SchurComplementOperator<TAlgebra>::
apply_sub(vector_type& fskeleton, const vector_type& uskeleton)
{
	SCHUR_PROFILE_BEGIN(SCHUR_Op_apply_sub);
	UG_DLOG(SchurDebug, 5, "\n% 'SchurComplementOperator::apply_sub()':");

//	create new rhs
	vector_type dskeleton(fskeleton);

//	solve
	apply(dskeleton, uskeleton);

//	subtract from vector
	fskeleton -= dskeleton;
}

template <typename TAlgebra>
void SchurComplementOperator<TAlgebra>::
debug_compute_matrix()
{
	SCHUR_PROFILE_BEGIN(SCHUR_Op_debug_compute_matrix);
	typedef typename TAlgebra::matrix_type sparse_matrix_type;
	typedef typename DenseMatrixFromSparseMatrix<sparse_matrix_type>::type dense_matrix_type;

	//sparse_matrix_type schur_matrix1(n_skeleton, n_skeleton);

	const int n_skeleton = sub_size(SD_SKELETON);

	dense_matrix_type schur_matrix;
	schur_matrix.resize(n_skeleton, n_skeleton);

	// create temporary vectors
	vector_type sol; sol.create(n_skeleton);
	vector_type rhs; rhs.create(n_skeleton);

	// resize matrix
//	schur_matrix.resize_and_clear(n_skeleton, n_skeleton);



	// compute columns s_k = S e_k
	for (int i=0; i<n_skeleton; ++i)
	{
		sol.set(0.0); sol[i] = 1.0;
		apply(rhs, sol);

		// copy to matrix
		for (int j=0; j<n_skeleton; ++j)
			schur_matrix(j, i) = rhs[j];
	}

	// print matrix
	UG_LOG(JuliaString(schur_matrix, "Schur"));
}

template <typename TAlgebra>
void SchurComplementOperator<TAlgebra>::
compute_matrix(matrix_type &schur_matrix)
{
	try{
	SCHUR_PROFILE_BEGIN(SCHUR_Op_compute_matrix);

	const int n_skeleton = sub_size(SD_SKELETON);

	schur_matrix.resize_and_clear(n_skeleton, n_skeleton);

	// create temporary vectors
	vector_type sol; sol.create(n_skeleton);
	vector_type rhs; rhs.create(n_skeleton);

	// resize matrix
//	schur_matrix.resize_and_clear(n_skeleton, n_skeleton);

	typename TAlgebra::matrix_type &mat = m_spOperator->get_matrix();
	schur_matrix.set_layouts(m_slicing.get_slice_layouts(mat.layouts(), SD_SKELETON));
	schur_matrix.set_storage_type(PST_ADDITIVE);

	UG_DLOG(SchurDebug, 2, "SchurMatrix Layouts: " << *schur_matrix.layouts());

	PROGRESS_START(prog, n_skeleton, "computing explicit Schur Matrix ( " << n_skeleton << " )");
	// compute columns s_k = S e_k
	for (int i=0; i<n_skeleton; ++i)
	{
		PROGRESS_UPDATE(prog, i);
		sol.set(0.0); sol[i] = 1.0;
		apply(rhs, sol);

		// copy to matrix
		for (int j=0; j<n_skeleton; ++j)
		{
			//schur_matrix(j, i) = rhs[j];
			if(rhs[j] != 0.0)
				schur_matrix(j, i) = rhs[j];
		}
	}

	PROGRESS_FINISH(prog);
	IF_DEBUG(SchurDebug, 2)
	{ schur_matrix.print("Schur"); }


	{
		SCHUR_PROFILE_BEGIN(SCHUR_Op_compute_matrix_wait);
		PROGRESS_START(prog2, 1, "Computing explicit Schur Matrix: Waiting for other processes..");
		mat.layouts()->proc_comm().barrier();

	}

	}UG_CATCH_THROW("SchurComplementOperator::" << __FUNCTION__ << " failed")
}




////////////////////////////////////////////////////////////////////////
//	template instantiations for all current algebra types.




#ifdef UG_CPU_1
template class SchurComplementOperator<CPUAlgebra>;
#endif
#ifdef UG_CPU_2
template class SchurComplementOperator<CPUBlockAlgebra<2> >;
#endif
#ifdef UG_CPU_3
template class SchurComplementOperator<CPUBlockAlgebra<3> >;
#endif

};  // end of namespace

#endif
