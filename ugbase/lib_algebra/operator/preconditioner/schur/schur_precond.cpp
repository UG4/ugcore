/*
 * schur_precond.cpp
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

namespace ug{

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	SchurSolver implementation
template <typename TAlgebra>
SchurPrecond<TAlgebra>::SchurPrecond() :
	//m_spOperator(NULL),
	m_spSchurComplementOp(NULL),
	m_spDirichletSolver(NULL),
	m_spSkeletonSolver(NULL)
{
	// clear aux vector smart ptrs
	// (will be initialized in first step)
	m_aux_rhs[0] = m_aux_rhs[1] =NULL;
	m_aux_sol[0] = m_aux_sol[1] = NULL;
}


template <typename TAlgebra>
bool SchurPrecond<TAlgebra>::
postprocess()
{
	// clear aux vector smart ptrs
	// (were initialized in first step)
	m_aux_rhs[0] = m_aux_rhs[1] =NULL;
	m_aux_sol[0] = m_aux_sol[1] = NULL;
	return true;
}

template <typename TAlgebra>
bool SchurPrecond<TAlgebra>::
create_and_init_local_schur_complement(SmartPtr<MatrixOperator<matrix_type, vector_type> > A,
		std::vector<slice_desc_type> &skeletonMark)
{
	try{
	SCHUR_PROFILE_BEGIN(SchurPrecondInit_CreateInitLocalSchurComplement);

	m_spSchurComplementOp = new SchurComplementOperator<TAlgebra>(A, skeletonMark);

//	set dirichlet solver for local Schur complement
	m_spSchurComplementOp->set_dirichlet_solver(m_spDirichletSolver);

	if(debug_writer().valid())
		m_spSchurComplementOp->set_debug(debug_writer());

//	init
	UG_DLOG(SchurDebug, 1, "\n%   - Init local Schur complement ... ");

	m_spSchurComplementOp->init();
	UG_DLOG(SchurDebug, 1, "done.\n");


//	1.4 check all procs
	/*bool bSuccess = true;
	if(!pcl::AllProcsTrue(bSuccess))
	{
		UG_LOG("ERROR in SchurPrecond::init: Some processes could not init"
				" local Schur complement.\n");
		return false;
	}*/
	return true;

	}UG_CATCH_THROW("SchurPrecond::" << __FUNCTION__ << " failed")
	return false;
}


template <typename TAlgebra>
void SchurPrecond<TAlgebra>::
init_skeleton_solver()
{
	try{
	SCHUR_PROFILE_BEGIN(SchurPrecond_InitSkeletonSolver);

	if(!m_spSkeletonSolver->init(m_spSchurComplementOp))
		UG_THROW("SchurPrecond::init: Failed to init skeleton solver.");

//	2.5 check all procs
/*	if(!pcl::AllProcsTrue(bSuccess))
	{
		UG_LOG("ERROR in SchurPrecond::init: Some processes could not init"
				" Schur complement inverse.\n");
		return false;
	}
	*/

	}UG_CATCH_THROW("SchurPrecond::" << __FUNCTION__ << " failed")

}

template <typename TAlgebra>
bool SchurPrecond<TAlgebra>::
check_requirements()
{
	if(m_spDirichletSolver.invalid())
	{
		UG_LOG("ERROR in SchurSolver: No dirichlet solver set "
				" for inversion of A_{II} in Local Schur complement.\n");
		return false;
	}

	if(m_spSkeletonSolver.invalid())
	{
		UG_LOG("ERROR in SchurPrecond: No skeleton solver set.\n");
		return false;
	}

	return true;
}

template <typename TAlgebra>
void SchurPrecond<TAlgebra>::
get_skeleton_slicing(SmartPtr<MatrixOperator<matrix_type, vector_type> > A,
		std::vector<slice_desc_type> &skeletonMark)
{
	matrix_type &Amat = A->get_matrix();
	const int N = Amat.num_rows();
	ConstSmartPtr<AlgebraLayouts> layouts = Amat.layouts();
	skeletonMark.clear();
	skeletonMark.resize(N, SD_INNER);
	MarkAllFromLayout<slice_desc_type> (skeletonMark, layouts->master(), SD_SKELETON);
	MarkAllFromLayout(skeletonMark, layouts->slave(), SD_SKELETON);
}

template <typename TAlgebra>
bool SchurPrecond<TAlgebra>::
preprocess(SmartPtr<MatrixOperator<matrix_type, vector_type> > A)
{
	try{
//	status
	UG_DLOG(SchurDebug, 2, "\n% Initializing SCHUR precond: \n");

	if(check_requirements() == false)
		return false;

//	Determine slicing for SchurComplementOperator
	std::vector<slice_desc_type> skeletonMark;
	get_skeleton_slicing(A, skeletonMark);

//	create & init local Schur complement object

	if(create_and_init_local_schur_complement(A, skeletonMark) == false)
		return false;

//  configure schur complement solver
	init_skeleton_solver();

//	status
	UG_DLOG(SchurDebug, 1, "\n% 'SchurPrecond::init()' done!\n");

//	we're done
	return true;

	}UG_CATCH_THROW("SchurPrecond::" << __FUNCTION__ << " failed");
	return false;
} /* end 'SchurPrecond::preprocess()' */


template <typename TAlgebra>
void SchurPrecond<TAlgebra>::
create_aux_vectors(const vector_type& d)
{
	const SlicingData sd =m_spSchurComplementOp->slicing();
	const size_t n_inner=m_spSchurComplementOp->sub_size(SD_INNER);
	const size_t n_skeleton=m_spSchurComplementOp->sub_size(SD_SKELETON);
	(void)n_skeleton; // warning fix
	if (m_aux_rhs[SD_SKELETON].invalid())
	{
		UG_DLOG(SchurDebug, 1, "% Creating skeleton defect vector of size " << n_skeleton << std::endl);
		//m_aux_rhs[SD_SKELETON] = new vector_type(n_skeleton);
		m_aux_rhs[SD_SKELETON] = sd.slice_clone_without_values(d, SD_SKELETON);
		m_aux_rhs[SD_SKELETON]->set_storage_type(PST_ADDITIVE);

		//std::cerr<< "Skeleton f:\n" <<*m_aux_rhs[SD_SKELETON]->layouts();
	}

	if (m_aux_sol[SD_SKELETON].invalid())
	{
		UG_DLOG(SchurDebug, 1, "% Creating skeleton corr vector of size " << n_skeleton << std::endl);
		//m_aux_sol[SD_SKELETON] = new vector_type(n_skeleton);
		m_aux_sol[SD_SKELETON] = sd.slice_clone_without_values(d, SD_SKELETON);
		m_aux_sol[SD_SKELETON]->set_storage_type(PST_CONSISTENT);

		//std::cerr<< "Skeleton u:\n" << *m_aux_sol[SD_SKELETON]->layouts();
	}

	if (m_aux_rhs[SD_INNER].invalid())
	{
		UG_DLOG(SchurDebug, 1, "% Creating inner defect vector of size " << n_inner << std::endl);
		m_aux_rhs[SD_INNER] = new vector_type(n_inner);
		m_aux_rhs[SD_INNER]->set_storage_type(PST_ADDITIVE);
	}

	if (m_aux_sol[SD_INNER].invalid())
	{
		UG_DLOG(SchurDebug, 1, "% Creating inner corr vector of size " << n_inner << std::endl);
		m_aux_sol[SD_INNER] = new vector_type(n_inner);
		m_aux_sol[SD_INNER]->set_storage_type(PST_CONSISTENT);
	}
}


//	Stepping routine
template <typename TAlgebra>
void SchurPrecond<TAlgebra>::
schur_solver_forward(vector_type &u_inner, vector_type &f_inner)
{
	UG_DLOG(SchurDebug, 3, "\n% 'SchurPrecond::step() - forward':");
	SCHUR_PROFILE_BEGIN(SchurSolverStep_Forward);

	// solve
	//UG_LOG("\nf_inner1="); UG_LOG_Vector<vector_type>(f_inner);

	m_spDirichletSolver->apply_return_defect(u_inner, f_inner);
	// store first correction -> will be used again
	//UG_LOG("\nu_inner1="); UG_LOG_Vector<vector_type>(u_inner);

	//UG_LOG("\nf_skeleton="); UG_LOG_Vector<vector_type>(f_skeleton);

	// slicing.subtract_vector_slice(d, SD_SKELETON, f_skeleton);
	/// f_skeleton *= -1.0;
}

template <typename TAlgebra>
void SchurPrecond<TAlgebra>::
schur_solve_skeleton(vector_type &u_skeleton, const vector_type &f_skeleton)
{
	SCHUR_PROFILE_BEGIN(SchurSolverStep_SchurSolve);

	UG_DLOG(SchurDebug, 3, "\n% 'SchurPrecond::step() - skeleton solve':");

	if(!f_skeleton.has_storage_type(PST_ADDITIVE))
	{ UG_THROW("ERROR: In 'SchurPrecond::step':Inadequate storage format of 'f_skeleton'.\n"); }

	if (!m_spSkeletonSolver->apply(u_skeleton, f_skeleton))
	{ UG_LOG("SchurPrecond: Failed to solve skeleton system!\n"); }

	if(!u_skeleton.has_storage_type(PST_CONSISTENT))
	{ UG_THROW("ERROR: In 'SchurPrecond::step':Inadequate storage format of 'u_skeleton'.\n"); }

	//UG_LOG("\nu_skeleton="); UG_LOG_Vector<vector_type>(u_skeleton);

}

template <typename TAlgebra>
void SchurPrecond<TAlgebra>::
schur_solver_backward(vector_type &u_inner, vector_type &f_inner, vector_type &u_skeleton)
{
	SCHUR_PROFILE_BEGIN(SchurSolverStep_Backward);

	UG_DLOG(SchurDebug, 3, "\n% 'SchurPrecond::step() - backward':\n");
	m_spSchurComplementOp->sub_operator(SD_INNER, SD_SKELETON)->apply(f_inner, u_skeleton);
	if(!m_spDirichletSolver->apply_return_defect(u_inner, f_inner) )
	{ UG_LOG("SchurPrecond: Failed to solve inner system!\n"); }
}

template <typename TAlgebra>
bool SchurPrecond<TAlgebra>::
step(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp, vector_type& c, const vector_type& d)
{
	try{
	bool bSuccess = true;	//	status

	c.set_storage_type(PST_UNIQUE);

	UG_DLOG(SchurDebug, 2, "\n% 'SchurPrecond::step()':");

	const size_t n_skeleton=m_spSchurComplementOp->sub_size(SD_SKELETON);

	//	check storage type
	if(!d.has_storage_type(PST_ADDITIVE))
	{
		UG_LOG("ERROR: In 'SchurPrecond::step':Inadequate storage format of 'd'.\n");
		return false;
	}

	// no skeleton => direct solve is sufficient
	if (n_skeleton == 0)
	{
		UG_LOG("\n% 'SchurPrecond::step() - direct solve':'\n");
		c.set_storage_type(PST_CONSISTENT);
		bSuccess = m_spDirichletSolver->apply(c, d);
		return bSuccess;
	}


	UG_ASSERT(n_skeleton > 0, "HUHH: #skeleton dof should be positive ");

	// now we have a non-trivial skeleton
	c.set(0.0);

	create_aux_vectors(d);


	// create short cuts
	const SlicingData &slicing = m_spSchurComplementOp->slicing();

	vector_type &f_skeleton=*m_aux_rhs[SD_SKELETON];
	vector_type &u_skeleton=*m_aux_sol[SD_SKELETON];

	vector_type &f_inner=*m_aux_rhs[SD_INNER];
	vector_type &u_inner=*m_aux_sol[SD_INNER];


	// get defect
	slicing.get_vector_slice(d, SD_INNER, f_inner);
	m_aux_rhs[SD_INNER]->set_storage_type(PST_ADDITIVE);

	slicing.get_vector_slice(d, SD_SKELETON, f_skeleton);
	m_aux_rhs[SD_SKELETON]->set_storage_type(PST_ADDITIVE);

	// get initial guess
	slicing.get_vector_slice(c, SD_INNER, u_inner);
	slicing.get_vector_slice(c, SD_SKELETON, u_skeleton);

	// A.  Reduce defect to skeleton (forward solve)
	schur_solver_forward(u_inner, f_inner);

	slicing.set_vector_slice(u_inner, c, SD_INNER);


	// update defect on skeleton
	m_spSchurComplementOp->sub_operator(SD_SKELETON, SD_INNER)->apply_sub(f_skeleton, u_inner);

	// B. solve system on skeleton
	schur_solve_skeleton(u_skeleton, f_skeleton);

	// set correction on skeleton slice of c
	slicing.set_vector_slice(u_skeleton, c, SD_SKELETON);

	// C. Prolongate correction back (backward solve)
	schur_solver_backward(u_inner, f_inner, u_skeleton);

	slicing.subtract_vector_slice(u_inner, c, SD_INNER);

	// c is consistent, since both slices are mutex and consistent
	c.set_storage_type(PST_CONSISTENT);

	return true;

	}UG_CATCH_THROW("SchurPrecond::" << __FUNCTION__ << " failed");
	return false;

} /* end 'SchurPrecond::step()' */





////////////////////////////////////////////////////////////////////////
//	template instantiations for all current algebra types.

#ifdef UG_CPU_1
template class SchurPrecond<CPUAlgebra>;
#endif
#ifdef UG_CPU_2
template class SchurPrecond<CPUBlockAlgebra<2> >;
#endif
#ifdef UG_CPU_3
template class SchurPrecond<CPUBlockAlgebra<3> >;
#endif

};  // end of namespace

#endif
