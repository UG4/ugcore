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

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__SCHUR_PRECOND__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__SCHUR_PRECOND__


#ifdef UG_PARALLEL

#include <iostream>
#include <sstream>
#include <string>
#include <set>

#include "lib_algebra/operator/interface/preconditioner.h"
#include "lib_algebra/operator/interface/linear_operator.h"
#include "lib_algebra/operator/interface/linear_operator_inverse.h"
#include "lib_algebra/operator/interface/preconditioned_linear_operator_inverse.h"
#include "lib_algebra/operator/interface/matrix_operator.h"
#include "lib_algebra/operator/interface/matrix_operator_inverse.h"
#include "lib_algebra/parallelization/parallelization.h"
#include "lib_algebra/operator/debug_writer.h"
#include "schur_complement_inverse_interface.h"
#include "pcl/pcl.h"

#include "common/log.h"

#include "schur.h"


#define PROFILE_SCHUR
#ifdef PROFILE_SCHUR
	#define SCHUR_PROFILE_FUNC()			PROFILE_FUNC_GROUP("algebra schur")
	#define SCHUR_PROFILE_BEGIN(name)	PROFILE_BEGIN_GROUP(name, "algebra schur")
	#define SCHUR_PROFILE_END_(name)			PROFILE_END_(name)
#else
	#define SCHUR_PROFILE_FUNC()
	#define SCHUR_PROFILE_BEGIN(name)
	#define SCHUR_PROFILE_END_(name)
#endif

namespace ug{

extern DebugID SchurDebug;


/// operator implementation of the DD Schur complement solver
/**
 * This operator implements a Schur complement solver */
template <typename TAlgebra>
class SchurPrecond: public IPreconditioner<TAlgebra>
{
	public:
	// 	Algebra type
		using algebra_type = TAlgebra;

	// 	Vector type
		using vector_type = typename TAlgebra::vector_type;

	// 	Matrix type
		using matrix_type = typename TAlgebra::matrix_type;

	///	Base type
		using base_type = IPreconditioner<TAlgebra>;

	protected:
		//using base_type::set_debug;
		using base_type::write_debug;
		using base_type::debug_writer;

	public:
	///	constructor
		SchurPrecond();

		SchurPrecond(const SchurPrecond &parent) : base_type(parent)
		{
			m_spDirichletSolver = parent.m_spDirichletSolver;
			m_spSkeletonSolver = parent.m_spSkeletonSolver;
			m_spSchurComplementOp = parent.m_spSchurComplementOp;
		}

		SmartPtr<ILinearIterator<vector_type> > clone() override {
			return make_sp(new SchurPrecond(*this));
		}
		~SchurPrecond() override = default;

	protected:
	///	name of solver
		[[nodiscard]] const char* name() const override {return "Schur complement";}

		//	Preprocess routine
		bool preprocess(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp) override;

		//	Stepping routine
		bool step(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp, vector_type& c, const vector_type& d) override;

		//	Postprocess routine
		bool postprocess() override;//  {return true;}

	///	returns if parallel solving is supported
		bool supports_parallel() const override {
			if(m_spDirichletSolver.valid()
				&& (!m_spDirichletSolver->supports_parallel()))
					return false;

			if(m_spSkeletonSolver.valid()
					&& (!m_spSkeletonSolver->supports_parallel()))
					return false;

			return true;
		}

public:

	///	sets the Dirichlet solver (forward to Schur complement)
		void set_dirichlet_solver(SmartPtr<ILinearOperatorInverse<vector_type> > dirichletSolver)
		{ m_spDirichletSolver = dirichletSolver; }

	///	sets the coarse problem solver
		void set_skeleton_solver(SmartPtr<ISchurComplementInverse<algebra_type> > skeletonSolver)
		{ m_spSkeletonSolver = skeletonSolver; }

		void set_schur_complement_operator(SmartPtr<SchurComplementOperator<algebra_type> > scop)
		{ m_spSchurComplementOp = scop; }


	//	set debug output
		void set_debug(SmartPtr<IDebugWriter<algebra_type> > spDebugWriter) override {
			base_type::set_debug(spDebugWriter);
			//m_spSchurComplementOp->set_debug(spDebugWriter);
		}

		std::string config_string() const override {
			std::stringstream ss; ss << name() << "\n";
			ss << " Dirichlet Solver: ";
			if(m_spDirichletSolver.valid()) ss << ConfigShift(m_spDirichletSolver->config_string()) << "\n";
			else ss << "  NOT SET!\n";
			ss << " Skeleton Solver: ";
			if(m_spSkeletonSolver.valid()) ss << ConfigShift(m_spSkeletonSolver->config_string()) << "\n";
			else ss << "  NOT SET!\n";

			return ss.str();
		}

public:
		/// returns the size of the skeleton problem
		/// note that this is the global size, i.e. without counting slaves
		size_t num_global_skeleton()
		{
			matrix_type &Amat = m_pA->get_matrix();
			ConstSmartPtr<AlgebraLayouts> layouts = Amat.layouts();
			return layouts->proc_comm().allreduce(m_myMasterSkeleton, PCL_RO_SUM);
		}
		/// returns the local skeleton size.
		/// note that the sum over these is bigger than num_global_skeleton,
		/// since this function includes also slaves
		size_t num_local_skeleton() const {
			return m_myTotalSkeleton;
		}


private:
		bool create_and_init_local_schur_complement(SmartPtr<MatrixOperator<matrix_type, vector_type> > A,
				std::vector<schur_slice_desc_type> &skeletonMark);

		void init_skeleton_solver();
		bool check_requirements();
		void get_skeleton_slicing(SmartPtr<MatrixOperator<matrix_type, vector_type> > A,
				std::vector<schur_slice_desc_type> &skeletonMark);

		void create_aux_vectors(const vector_type& d);
		void schur_solver_forward(vector_type &u_inner, vector_type &f_inner);
		void schur_solve_skeleton(vector_type &u_skeleton, const vector_type &f_skeleton);
		void schur_solver_backward(vector_type &u_inner, vector_type &f_inner, vector_type &u_skeleton);


	protected:
	// 	Reference to operator that is inverted by this Inverse Operator
	//	SmartPtr<MatrixOperator<matrix_type,vector_type> > m_spOperator;

	/// Solver Dirichlet problems \f$A_{II}\f$ (also used in Schur complement)
	SmartPtr<ILinearOperatorInverse<vector_type> > m_spDirichletSolver;

	///	Solver for coarse (skeleton) problem
	SmartPtr< ISchurComplementInverse<TAlgebra> > m_spSkeletonSolver;

	///	Local Schur complement for each subdomain
	SmartPtr<SchurComplementOperator<algebra_type> > m_spSchurComplementOp;

	// temporary vectors for correction/defect
	SmartPtr<vector_type> m_aux_rhs[2];
	SmartPtr<vector_type> m_aux_sol[2];

	size_t m_myMasterSkeleton, m_myTotalSkeleton;
	SmartPtr<MatrixOperator<matrix_type, vector_type> > m_pA;
	//	pointer to Domain decomposition info object
	//	pcl::IDomainDecompositionInfo* m_pDDInfo;
};


}
#endif
#endif
