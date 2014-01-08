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
	m_spSkeletonMatrix(NULL),
	m_spDirichletSolver(NULL),
	m_spSkeletonSolver(NULL)
{
	m_bExactSchurComplement = false;
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
preprocess(SmartPtr<MatrixOperator<matrix_type, vector_type> > A)
{
//	status
	UG_DLOG(SchurDebug, 2, "\n% Initializing SCHUR precond: \n");

	const SlicingData::slice_desc_type SD_INNER=SlicingData::SD_INNER;
	const SlicingData::slice_desc_type SD_SKELETON=SlicingData::SD_SKELETON;


//	bool flag
	bool bSuccess = true;

//	remember A
	//m_spOperator = A;

//	0. get matrix
	//matrix_type *m_pMatrix = &m_spOperator->get_matrix();
	matrix_type &Amat = A->get_matrix();
	const int N = Amat.num_rows();



//	check that DDInfo has been set
/*	if(m_pDDInfo == NULL)
	{
		UG_LOG("ERROR in SchurSolver::init: DDInfo not set.\n");
		return false;
	}*/

	bool debugLayouts = (debug_writer()==NULL) ? false : true;
// Determine splitting


//  ----- 1. CONFIGURE LOCAL SCHUR COMPLEMENT  ----- //

//	1.1 Determine slicing for SchurComplementOperator
	ConstSmartPtr<AlgebraLayouts> layouts = Amat.layouts();
	std::vector<SlicingData::slice_desc_type> skeletonMark(N, SlicingData::SD_INNER);
	MarkAllFromLayout<SlicingData::slice_desc_type> (skeletonMark, layouts->master(), SlicingData::SD_SKELETON);
	MarkAllFromLayout(skeletonMark, layouts->slave(), SlicingData::SD_SKELETON);

//	1.2 init Dirichlet system solver
	if(m_spDirichletSolver.invalid())
	{
		UG_LOG("ERROR in SchurSolver::init: No dirichlet solver set "
				" for inversion of A_{II} in Local Schur complement.\n");
		return false;
	}

	if(m_spSkeletonSolver.invalid())
	{
		UG_LOG("ERROR in SchurPrecond::init: No skeleton solver set.\n");
		return false;
	}

//	1.3 create & init local Schur complement object
	m_spSchurComplementOp = new SchurComplementOperator<TAlgebra>(A, skeletonMark);

//	set dirichlet solver for local Schur complement
	m_spSchurComplementOp->set_dirichlet_solver(m_spDirichletSolver);

//	m_spSchurComplementOp->set_exact_schur_complement(m_bExactSchurComplement);


//	init
	UG_DLOG(SchurDebug, 1, "\n%   - Init local Schur complement ... ");
	SCHUR_PROFILE_BEGIN(SchurPrecondInit_InitLocalSchurComplement);
	m_spSchurComplementOp->init();
	UG_DLOG(SchurDebug, 1, "done.\n");

	//if (debug_writer().valid())
	//		m_spSchurComplementOp->set_debug(debug_writer());

//	1.4 check all procs
	if(!pcl::AllProcsTrue(bSuccess))
	{
		UG_LOG("ERROR in SchurPrecond::init: Some processes could not init"
				" local Schur complement.\n");
		return false;
	}
	SCHUR_PROFILE_END_(SchurPrecondInit_InitLocalSchurComplement);

//  ----- 2. CONFIGURE SCHUR COMPLEMENT SOLVER  ----- //

// Extension for preconditioned Schur complement solve
	//SmartPtr<IPreconditioner<TAlgebra> > precond = new Jacobi<TAlgebra>();
/*
	SmartPtr<IPreconditioner<TAlgebra> > precond = new Jacobi<TAlgebra>(0.0);
	precond->set_approximation(m_spSkeletonMatrix);
	precond->set_damp(0.5);
	m_spSkeletonSolver->set_preconditioner(precond);

*/
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

//	status
	UG_DLOG(SchurDebug, 1, "\n% 'SchurPrecond::init()' done!\n");

//	we're done
	return true;
} /* end 'SchurPrecond::preprocess()' */




//	Stepping routine


template <typename TAlgebra>
bool SchurPrecond<TAlgebra>::
step(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp, vector_type& c, const vector_type& d)
{

	bool bSuccess = true;	//	status

	c.set_storage_type(PST_UNIQUE);

	UG_DLOG(SchurDebug, 2, "\n% 'SchurPrecond::step()':");

	const SlicingData::slice_desc_type SD_INNER=SlicingData::SD_INNER;
	const SlicingData::slice_desc_type SD_SKELETON=SlicingData::SD_SKELETON;

	const size_t n_inner=m_spSchurComplementOp->sub_size(SD_INNER);
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
		UG_LOG("\n% 'SchurPrecond::step() - direct solve':");
		c.set_storage_type(PST_CONSISTENT);
		m_spDirichletSolver->apply(c, d);
		return bSuccess;
	}

	const SlicingData sd =m_spSchurComplementOp->slicing();
	UG_ASSERT(n_skeleton > 0, "HUHH: #skeleton dof should be positive ");

	// now we have a non-trivial skeleton
	c.set(0.0);
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



	// create short cuts
	const SlicingData&slicing = m_spSchurComplementOp->slicing();

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
	UG_DLOG(SchurDebug, 3, "\n% 'SchurPrecond::step() - forward':");
	SCHUR_PROFILE_BEGIN(SchurSolverStep_Forward);

	// solve
	//UG_LOG("\nf_inner1="); UG_LOG_Vector<vector_type>(f_inner);

	m_spDirichletSolver->apply_return_defect(u_inner, f_inner);
	// store first correction -> will be used again
	slicing.set_vector_slice(u_inner, c, SD_INNER);
	//UG_LOG("\nu_inner1="); UG_LOG_Vector<vector_type>(u_inner);

	// update defect on skeleton
	m_spSchurComplementOp->sub_operator(SD_SKELETON, SD_INNER)->apply_sub(f_skeleton, u_inner);
	//UG_LOG("\nf_skeleton="); UG_LOG_Vector<vector_type>(f_skeleton);

	// slicing.subtract_vector_slice(d, SD_SKELETON, f_skeleton);
	/// f_skeleton *= -1.0;
	SCHUR_PROFILE_END_(SchurSolverStep_Forward);

	// B. solve system on skeleton
	SCHUR_PROFILE_BEGIN(SchurSolverStep_SchurSolve);
	UG_DLOG(SchurDebug, 3, "\n% 'SchurPrecond::step() - skeleton solve':");

	if(!f_skeleton.has_storage_type(PST_ADDITIVE))
	{ UG_THROW("ERROR: In 'SchurPrecond::step':Inadequate storage format of 'f_skeleton'.\n"); }

	if (!m_spSkeletonSolver->apply(u_skeleton, f_skeleton))
	{ UG_LOG("Failed to solve system!"); }

	if(!u_skeleton.has_storage_type(PST_CONSISTENT))
	{ UG_THROW("ERROR: In 'SchurPrecond::step':Inadequate storage format of 'u_skeleton'.\n"); }

	slicing.set_vector_slice(u_skeleton, c, SD_SKELETON);

	//UG_LOG("\nu_skeleton="); UG_LOG_Vector<vector_type>(u_skeleton);
	SCHUR_PROFILE_END_(SchurSolverStep_SchurSolve);

	// C. Prolongate correction back (backward solve)
	SCHUR_PROFILE_BEGIN(SchurSolverStep_Backward);
	UG_DLOG(SchurDebug, 3, "\n% 'SchurPrecond::step() - backward':\n");
	m_spSchurComplementOp->sub_operator(SD_INNER, SD_SKELETON)->apply(f_inner, u_skeleton);
	m_spDirichletSolver->apply_return_defect(u_inner, f_inner);
	slicing.subtract_vector_slice(u_inner, c, SD_INNER);

	// c is consistent, since both slices are mutex and consistent
	c.set_storage_type(PST_CONSISTENT);

	//	check all procs
	if(!pcl::AllProcsTrue(bSuccess))
	{
		UG_LOG("ERROR in SchurSolver::apply: Some processes could not back solve.\n");
		return false;
	}
	SCHUR_PROFILE_END_(SchurSolverStep_Backward);

	//m_spSchurComplementOp->debug_compute_matrix();
	return bSuccess;

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
