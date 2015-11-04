
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

#include "slicing.h"

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
//template <typename TAlgebra>
//class SchurSolver : public IMatrixOperatorInverse<	typename TAlgebra::matrix_type,
//													typename TAlgebra::vector_type>,
//	public DebugWritingObject<TAlgebra>

template <typename TAlgebra>
class SchurPrecond: public IPreconditioner<TAlgebra>
{
	public:
	// 	Algebra type
		typedef TAlgebra algebra_type;

	// 	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	// 	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

	///	Base type
		typedef IPreconditioner<TAlgebra> base_type;

	protected:
		//using base_type::set_debug;
		using base_type::write_debug;
		using base_type::debug_writer;

	public:
	///	constructor
		SchurPrecond();

		SchurPrecond(const SchurPrecond<TAlgebra> &parent) : base_type(parent)
		{
			m_spDirichletSolver = parent.m_spDirichletSolver;
			m_spSkeletonSolver = parent.m_spSkeletonSolver;
			m_spSchurComplementOp = parent.m_spSchurComplementOp;
		}

		virtual SmartPtr<ILinearIterator<vector_type> > clone()
		{
			return make_sp(new SchurPrecond<algebra_type>(*this));
		}

	protected:
	///	name of solver
		virtual const char* name() const {return "Schur complement";}

		//	Preprocess routine
		virtual bool preprocess(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp);

		//	Stepping routine
		virtual bool step(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp, vector_type& c, const vector_type& d);

		//	Postprocess routine
		virtual bool postprocess();//  {return true;}

	///	returns if parallel solving is supported
		virtual bool supports_parallel() const
		{
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
		void set_debug(SmartPtr<IDebugWriter<algebra_type> > spDebugWriter)
		{
			base_type::set_debug(spDebugWriter);
			//m_spSchurComplementOp->set_debug(spDebugWriter);
		}

		virtual std::string config_string() const
		{
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
		size_t num_local_skeleton()
		{
			return m_myTotalSkeleton;
		}


private:
		bool create_and_init_local_schur_complement(SmartPtr<MatrixOperator<matrix_type, vector_type> > A,
				std::vector<slice_desc_type> &skeletonMark);

		void init_skeleton_solver();
		bool check_requirements();
		void get_skeleton_slicing(SmartPtr<MatrixOperator<matrix_type, vector_type> > A,
				std::vector<slice_desc_type> &skeletonMark);

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
#endif /* UG_PARALLEL */
#endif /* __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__SCHUR_PRECOND__ */
