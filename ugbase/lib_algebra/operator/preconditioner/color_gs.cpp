/*
 * Copyright (c) 2014-2015:  G-CSC, Goethe University Frankfurt
 * Author: Arne Nägel
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

//	THIS FILE IS ONLY TEMPORARY!
//  When everything works as
//	expected one should move the code in this file to
//	schur_impl.h, so that it will work with all the different algebras.
//
//	In the moment, template instantiations are invoked at the end of this file.




// extern headers
#include <cmath>
#include <sstream>  // added for 'stringstream'

// algebra types
#include "lib_algebra/cpu_algebra_types.h"
#include "lib_algebra/operator/algebra_debug_writer.h"
#include "lib_algebra/operator/preconditioner/preconditioners.h"

#include "lib_algebra/algebra_common/sparsematrix_util.h"  // DenseMatrixFromSparseMatrix
#include "lib_algebra/small_algebra/small_algebra.h"   // DEnseMAtrix...
#include "lib_algebra/algebra_template_define_helper.h"

#ifdef UG_PARALLEL
#include "lib_algebra/parallelization/parallel_index_layout.h"
#include "pcl/pcl_layout_tests.h"
#endif

#include "common/profiler/profiler.h"        // additions for profiling

// own header
//#include "schur/schur.h"
//#include "schur/schur_precond.h"
//#include "schur/schur_complement_inverse_interface.h"
#include "color_gs.h"

namespace ug{

namespace ColorGS{


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	SchurSolver implementation
template <typename TAlgebra>
ColorPrecond<TAlgebra>::ColorPrecond()
{
	clear_aux_vectors();
}



/*

template <typename TAlgebra>
bool ColorPrecond<TAlgebra>::
create_and_init_local_schur_complement(SmartPtr<MatrixOperator<matrix_type, vector_type> > A,
		std::vector<schur_slice_desc_type> &skeletonMark)
{
	SCHUR_PROFILE_BEGIN(ColorPrecondInit_CreateInitLocalSchurComplement);

	try{


	/*m_spSchurComplementOp = make_sp(new SchurComplementOperator<TAlgebra>(A, skeletonMark));
	if (m_spSchurComplementOp.invalid())
	{
		UG_ASSERT(m_spSchurComplementOp.invalid(), "Failed creating operator!")
	}

//	set dirichlet solver for local Schur complement
	m_spSchurComplementOp->set_dirichlet_solver(m_spDirichletSolver);

	if(debug_writer().valid())
		m_spSchurComplementOp->set_debug(debug_writer());

//	init
	UG_DLOG(SchurDebug, 1, "\n%   - Init local Schur complement ... ");

	m_spSchurComplementOp->init();
	UG_DLOG(SchurDebug, 1, "done.\n");



//	1.4 check all procs
	bool bSuccess = true;
	if(!pcl::AllProcsTrue(bSuccess))
	{
		UG_LOG("ERROR in ColorPrecond::init: Some processes could not init"
				" local Schur complement.\n");
		return false;
	}
	return true;

	}UG_CATCH_THROW("ColorPrecond::" << __FUNCTION__ << " failed")
	return false;
}*/



template <typename TAlgebra>
bool ColorPrecond<TAlgebra>::
check_requirements()
{
	/*
	if(m_spDirichletSolver.invalid())
	{
		UG_LOG("ERROR in SchurSolver: No dirichlet solver set "
				" for inversion of A_{II} in Local Schur complement.\n");
		return false;
	}

	if(m_spSkeletonSolver.invalid())
	{
		UG_LOG("ERROR in ColorPrecond: No skeleton solver set.\n");
		return false;
	}
*/
	return true;
}


template <typename TAlgebra>
bool ColorPrecond<TAlgebra>::
preprocess(SmartPtr<MatrixOperator<matrix_type, vector_type> > A)
{

	try{
//	status
	UG_DLOG(SchurDebug, 2, "\n% Initializing SCHUR precond: \n");

	m_pA = A;
	if(check_requirements() == false)
		return false;

	/*
//	Determine slicing for SchurComplementOperator
	std::vector<schur_slice_desc_type> skeletonMark;
	get_skeleton_slicing(A, skeletonMark);

//	create & init local Schur complement object

	if(create_and_init_local_schur_complement(A, skeletonMark) == false)
		return false;
*/
//  configure schur complement solver
// 	init_skeleton_solver();

	create_diag_vectors(*A); // may happen on device

//	status
	UG_DLOG(SchurDebug, 1, "\n% 'ColorPrecond::preprocess(A)' done!\n");

//	we're done
	return true;

	}UG_CATCH_THROW("ColorPrecond::" << __FUNCTION__ << " failed");
	return false;
} /* end 'ColorPrecond::preprocess()' */



template <typename TAlgebra>
bool ColorPrecond<TAlgebra>::
postprocess()
{
	clear_aux_vectors();
	return true;
}


template <typename TAlgebra>
void ColorPrecond<TAlgebra>::
create_aux_vectors(const vector_type& d)
{
	const SlicingData &sd = slicing();

	for (size_t i = 0; i<ColorGS::MAX_COLORS; i++)
	{

		const size_t nelems = slicing().get_num_elems(i); // elements for this color


		// create aux vectors defect
		if (m_aux_rhs[i].invalid())
		{
			UG_DLOG(SchurDebug, 1, "% Creating rhs vector of size " << n_skeleton << std::endl);
			//m_aux_rhs[SD_SKELETON] = new vector_type(n_skeleton);
			m_aux_rhs[i] = sd.slice_clone_without_values(d, i);
#ifdef UG_PARALLEL
			m_aux_rhs[i]->set_storage_type(PST_ADDITIVE);
#endif
		}

		// create aux vectors solution
		if (m_aux_sol[i].invalid())
		{
			UG_DLOG(SchurDebug, 1, "% Creating corr vector of size " << n_skeleton << std::endl);
			m_aux_sol[i] = sd.slice_clone_without_values(d, i);
#ifdef UG_PARALLEL
			m_aux_sol[i]->set_storage_type(PST_CONSISTENT);
#endif
		}



	} // for all colors

}

template <typename TAlgebra>
void ColorPrecond<TAlgebra>::
clear_aux_vectors()
{
	for (size_t i = 0; i<ColorGS::MAX_COLORS; i++)
	{
		m_aux_rhs[i] = m_aux_sol[i] = SPNULL;
	}

}

template <typename TAlgebra>
void ColorPrecond<TAlgebra>::
create_diag_vectors(const matrix_type& mat)
{
	const SlicingData &sd = slicing();
	size_t count = 0;

	// Create storage.
	for (size_t c = 0; c<ColorGS::MAX_COLORS; c++)
	{

		const size_t nelems = slicing().get_num_elems(c); // elements for this color
		m_diag[c].reserve(nelems);
		count += nelems;

	} // for all colors


	// Extract matrix diagonal.
#ifdef UG_PARALLEL
					//	temporary vector for the diagonal
	ParallelVector<Vector< typename matrix_type::value_type > > diag;
	diag.resize(size);

	//	copy the layouts+communicator into the vector
	diag.set_layouts(mat.layouts());

	// 	copy diagonal
	for(size_t i = 0; i < diag.size(); ++i){
		diag[i] = mat(i, i);
	}

	//	make diagonal consistent
	diag.set_storage_type(PST_ADDITIVE);
	diag.change_storage_type(PST_CONSISTENT);


	if ((diag.size() > 0) && (CheckVectorInvertible(diag) == false))
		return false;
//			UG_ASSERT(CheckVectorInvertible(diag), "Jacobi: A has noninvertible diagonal");

	for(size_t i = 0; i < mat.size(); ++i){
		size_t color=sd.get_type(i);
		m_diag[color].push_back(diag[i]);
	}

	diag.clear(); // will be executed automatically anyway...
#else
	for(size_t i = 0; i < mat.num_rows(); ++i){
		size_t color=sd.get_type(i);
		m_diag[color].push_back(mat(i, i));
	}
	UG_ASSERT(count == mat.num_rows(), "Sizes do not match!")
#endif











}
template <typename TAlgebra>
void ColorPrecond<TAlgebra>::
clear_diag_vectors()
{
	for (size_t i = 0; i<ColorGS::MAX_COLORS; i++)
	{
		m_diag[i].clear();
	}

}

/*
//	Stepping routine
template <typename TAlgebra>
void ColorPrecond<TAlgebra>::
schur_solver_forward(vector_type &u_inner, vector_type &f_inner)
{
	UG_DLOG(SchurDebug, 3, "\n% 'ColorPrecond::step() - forward':");
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
*/



// solve a*c = beta*d
void InverseMatMult(number & ci, double beta, const number &Aii, const number & di)
{
	ci = beta*di/Aii;
}




template <typename TAlgebra>
bool ColorPrecond<TAlgebra>::
step(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp, vector_type& c, const vector_type& d)
{
	PROFILE_BEGIN(ColorPrecond_step)
	try{
	bool bSuccess = true;	//	status
#ifdef UG_PARALLEL
	c.set_storage_type(PST_UNIQUE);

	UG_DLOG(SchurDebug, 2, "\n% 'ColorPrecond::step()':");

	// const size_t n_skeleton=m_spSchurComplementOp->sub_size(SD_SKELETON);

	//	check storage type

	if(!d.has_storage_type(PST_ADDITIVE))
	{
		UG_LOG("ERROR: In 'ColorPrecond::step':Inadequate storage format of 'd'.\n");
		return false;
	}
#endif

	// now we have a non-trivial skeleton
	c.set(0.0);

	// create all auxiliary vectors (if neccessary)
	create_aux_vectors(d);


	// create short cuts
	const SlicingData &slicing = this->slicing();


	for (size_t i = 0; i<ColorGS::MAX_COLORS; i++)
	{
		const size_t nelems = slicing.get_num_elems(i);

		vector_type &di=*m_aux_rhs[i];
		vector_type &ci=*m_aux_sol[i];
		const auto DiagAi = m_diag[i];


		// Some checks.
		UG_ASSERT(nelems == m_diag[i].size(),  "Vector sizes do not match...")
		UG_ASSERT(nelems == di.size(), "Vector sizes do not match...")
		UG_ASSERT(nelems == ci.size(), "Vector sizes do not match...")

		// UG_DLOG("Iter: "<< i << "(" << nelems << ")" << std::endl);


		// Get (original) defect.
		slicing.get_vector_slice(d, i, di);


#ifdef UG_PARALLEL
		m_aux_rhs[i]->set_storage_type(PST_ADDITIVE);
#endif

		// Get initial guess (should be 0).
		slicing.get_vector_slice(c, i, ci);



		// Consider all previous updates (TODO: this is slow...)
		for (size_t j = 0; i<j; ++j)
		{

			// di = di - Aij*cj
			vector_type &cj=*m_aux_sol[j];

			matrix_type Aij;
			slicing.get_matrix(pOp->get_matrix(), i,j, Aij);

			Aij.axpy(di, 1.0, di, -1.0, cj);

		}


		//  Jacobi iteration (in PARALLEL!).
		#pragma omp parallel for
		for (size_t ii=0; ii<nelems; ++ii)
		{
			InverseMatMult(ci[ii], 1.0, DiagAi[ii], di[ii]);
		}

		// Set correction.
		slicing.set_vector_slice(ci, c, i);


		// Update defect (for j>i).




	}  // end of loop over colors.

#ifdef UG_PARALLEL
	// c is consistent, as qll slices are consistent
	c.set_storage_type(PST_CONSISTENT);
#endif
	return true;

	}UG_CATCH_THROW("ColorPrecond::" << __FUNCTION__ << " failed");
	return false;

} /* end 'ColorPrecond::step()' */






}  // end of namespace ColorGS


////////////////////////////////////////////////////////////////////////
//	template instantiations for all current algebra types.


// define ColorPrecond< ALGEBRA > for all Algebra Types (see algebra_template_define_helper.h)
UG_ALGEBRA_CPP_TEMPLATE_DEFINE_ALL(ColorGS::ColorPrecond)

}  // end of namespace ug


