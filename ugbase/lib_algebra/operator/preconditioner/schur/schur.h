/*
 * schur.h
 *
 *  Created on: 18.12.2013
 *      Author: anaegel
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__SCHUR__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__SCHUR__


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
#include "pcl/pcl.h"

#include "common/log.h"

#include "slicing.h"


namespace ug{

extern DebugID SchurDebug;

#if 0
template <typename TAlgebra>
class SchurComplementOperator;

template<typename TAlgebra>
class ISchurComplementApproximateInverse : public ILinearOperatorInverse<typename TAlgebra::vector_type>
{
public:
	virtual void init(SchurComplementOperator<TAlgebra> op) = 0;
};

template<typename TAlgebra>
class ExactSchurInverse : public ISchurComplementApproximateInverse<TAlgebra>
{
public:
	typedef TAlgebra algebra_type;
	typedef typename TAlgebra::vector_type vector_type;
	typedef typename TAlgebra::matrix_type matrix_type;

	ExactSchurInverse(SmartPtr<ILinearOperatorInverse<vector_type> > linOpInv )
	{
		m_linOpInv = linOpInv;
	}

	virtual void init(SchurComplementOperator<TAlgebra> op)
	{
		m_exactSchurOp = new MatrixOperator<matrix_type, vector_type>;
		op.compute_matrix(m_exactSchurOp->get_matrix());
		m_linOpInv.init(m_exactSchurOp);
	}

	virtual bool apply(vector_type& u, const vector_type& f)
	{
		return m_linOpInv->apply(u, f);
	}

	virtual bool apply_return_defect(vector_type& u, vector_type& f)
	{
		return m_linOpInv->apply_return_defect(u, f);
	}

protected:
	SmartPtr<MatrixOperator<matrix_type, vector_type> > m_exactSchurOp;
	SmartPtr<ILinearOperatorInverse<vector_type> > m_linOpInv;

};
#endif

/*
template<typename TAlgebra>
class ApproximateSchurInverse : public ISchurComplementApproximateInverse<TAlgebra>
{
public:
	typedef TAlgebra algebra_type;
	typedef typename TAlgebra::vector_type vector_type;
	typedef typename TAlgebra::matrix_type matrix_type;

	ApproximateSchurInverse(SmartPtr<IPreconditionedLinearOperatorInverse<vector_type> > linSolver)
	{
		m_linSolver = linSolver;
	}

	virtual void init(SchurComplementOperator<TAlgebra> op)
	{
		SmartPtr<IPreconditioner<TAlgebra> > precond = new Jacobi<TAlgebra>;
		precond->set_approximation(op.sub_operator(SD_SKELETON, SD_SKELETON));
		precond->set_damp(0.5);
		m_linSolver->set_preconditioner(precond);
	}
	virtual bool apply(vector_type& u, const vector_type& f)
	{
		return m_linSolver->apply(u, f);
	}

	virtual bool apply_return_defect(vector_type& u, vector_type& f)
	{
		return m_linSolver->apply_return_defect(u, f);
	}
protected:
	SmartPtr<IPreconditionedLinearOperatorInverse<vector_type> > m_linSolver;
}
*/


template <typename TAlgebra>
class SchurComplementOperator
	: public ILinearOperator<	typename TAlgebra::vector_type,
	  	  	  	  	  	  	  	typename TAlgebra::vector_type>,
	  public DebugWritingObject<TAlgebra>
{
	public:

	// 	Algebra type
		typedef TAlgebra algebra_type;

	// 	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	// 	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

	protected:
		using DebugWritingObject<TAlgebra>::write_debug;
		using DebugWritingObject<TAlgebra>::debug_writer;

	public:
	///	constructor
	SchurComplementOperator(SmartPtr<MatrixOperator<matrix_type, vector_type> > Alocal,
							SlicingData::slice_desc_type_vector &sdv)
	: m_spOperator(Alocal),
	  m_slicing(sdv)
	{
		m_op[0][0] = new MatrixOperator<matrix_type, vector_type>();
		m_op[0][1] = new MatrixOperator<matrix_type, vector_type>();
		m_op[1][0] = new MatrixOperator<matrix_type, vector_type>();
		m_op[1][1] = new MatrixOperator<matrix_type, vector_type>();
	}

	// destructor
	virtual ~SchurComplementOperator() {};

	///	name of solver
	virtual const char* name() const {return "My local Schur complement Solver";}


	/// implementation of the operator for the solution dependent initialization.
	void init(const vector_type& u) {init();}

	///	initializes the solver for operator A
	virtual void init();

	///	applies the Schur complement built from matrix operator set via 'set_matrix()'
	/// to 'u' and returns the result 'f := S times u'
	virtual void apply(vector_type& f, const vector_type& u);

	///	applies the Schur complement built from matrix operator set via 'set_matrix()'
	/// to 'u' and returns the result 'f := f - S times u'
	virtual void apply_sub(vector_type& f, const vector_type& u);

	//	save current operator
	void set_matrix(SmartPtr<MatrixOperator<matrix_type, vector_type> > A)
	{ m_spOperator = A; }

	///	sets a Dirichlet solver
	void set_dirichlet_solver(SmartPtr<ILinearOperatorInverse<vector_type> > dirichletSolver)
	{ m_spDirichletSolver = dirichletSolver; }


	matrix_type &sub_matrix(int r, int c)
	{return sub_operator(r,c)->get_matrix();}

	SmartPtr<MatrixOperator<matrix_type, vector_type> > sub_operator(int r, int c)
	{return m_op[r][c];}

	size_t sub_size(SlicingData::slice_desc_type type)
	{return m_slicing.get_num_elems(type);}

	const SlicingData &slicing() const
	{return m_slicing;}


	// for debugging: computes schur operator
	void debug_compute_matrix();

	void compute_matrix(matrix_type &schur_matrix);

protected:
	// 	Operator that is inverted by this Inverse Operator
	SmartPtr<MatrixOperator<matrix_type,vector_type> > m_spOperator;

	// slices from matrix
	const SlicingData m_slicing;

	// 	Linear Solver to invert the local Dirichlet problems
	SmartPtr<ILinearOperatorInverse<vector_type> > m_spDirichletSolver;



	// sub matrices/operator (will be set by init)
	SmartPtr<MatrixOperator<matrix_type,vector_type> > m_op[2][2];

	int m_applyCnt;

}; /* end class 'LocalSchurComplement' */



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

		virtual SmartPtr<ILinearIterator<vector_type> > clone()
		{
			SmartPtr<SchurPrecond<algebra_type> > newInst(new SchurPrecond<algebra_type>());
			UG_THROW("Implement SchurPrecond::clone()!")
			return newInst;
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
		void set_exact_schur_complement(bool b)
		{
			m_bExactSchurComplement = b;
		}
		void set_schur_complement_operator(SmartPtr<SchurComplementOperator<algebra_type> > scop)
		{ m_spSchurComplementOp = scop; }

		//	define an approximation for schur complement
		void set_schur_complement_approx(SmartPtr<MatrixOperator<matrix_type, vector_type> > S)
		{ m_spSkeletonMatrix = S; }

	///	sets the Dirichlet solver (forward to Schur complement)
		void set_dirichlet_solver(SmartPtr<ILinearOperatorInverse<vector_type> > dirichletSolver)
		{ m_spDirichletSolver = dirichletSolver; }

	///	sets the coarse problem solver
		void set_skeleton_solver(SmartPtr<IPreconditionedLinearOperatorInverse<vector_type> > skeletonSolver)
		{ m_spSkeletonSolver = skeletonSolver; }

	//	set debug output
		void set_debug(SmartPtr<IDebugWriter<algebra_type> > spDebugWriter)
		{
			base_type::set_debug(spDebugWriter);
			//m_spSchurComplementOp->set_debug(spDebugWriter);
		}

		virtual std::string config_string() const
		{
			std::stringstream ss; ss << name() << "\n";
			if(m_bExactSchurComplement)
				ss << " Using Exact Schur Complement for Inverse\n";
			else
				ss << " Using approximated Schur Complement for Inverse\n";
			ss << " Dirichlet Solver: ";
			if(m_spDirichletSolver.valid()) ss << ConfigShift(m_spDirichletSolver->config_string()) << "\n";
			else ss << "  NOT SET!\n";
			ss << " Skeleton Solver: ";
			if(m_spSkeletonSolver.valid()) ss << ConfigShift(m_spSkeletonSolver->config_string()) << "\n";
			else ss << "  NOT SET!\n";

			return ss.str();
		}

	/*///	initializes the solver for operator A
		virtual bool init(SmartPtr<MatrixOperator<matrix_type, vector_type> > A);

	///	solves the reduced system \f$F \lambda = d\f$ with preconditioned cg method
		virtual bool apply_return_defect(vector_type& x, vector_type& d);

	///	solves the system
		virtual bool apply(vector_type& x, const vector_type& b)
		{
		//	copy defect
			vector_type d; d.resize(b.size());
			d = b;

		//	solve on copy of defect
			return apply_return_defect(x, d);
		}*/

		// destructor
//		virtual ~SchurPrecond() {};


			/*void set_domain_decomp_info(pcl::IDomainDecompositionInfo& ddInfo)
			{
				m_pDDInfo = &ddInfo;
			}*/




	protected:
	// 	Reference to operator that is inverted by this Inverse Operator
	//	SmartPtr<MatrixOperator<matrix_type,vector_type> > m_spOperator;

	///	Local Schur complement for each subdomain
	SmartPtr<SchurComplementOperator<algebra_type> > m_spSchurComplementOp;

	/// Approximation of the Schur complement for preconditioner
	SmartPtr<MatrixOperator<matrix_type,vector_type> > m_spSkeletonMatrix;

	/// Solver Dirichlet problems \f$A_{II}\f$ (also used in Schur complement)
	SmartPtr<ILinearOperatorInverse<vector_type> > m_spDirichletSolver;

	///	Solver for coarse (skeleton) problem
	SmartPtr<ILinearOperatorInverse<vector_type> > m_spSkeletonSolver;

	// temporary vectors for correction/defect
	SmartPtr<vector_type> m_aux_rhs[2];
	SmartPtr<vector_type> m_aux_sol[2];

	//	pointer to Domain decomposition info object
	//	pcl::IDomainDecompositionInfo* m_pDDInfo;


	int m_iterCnt;
	bool m_bExactSchurComplement;

}; /* end class 'SchurPrecond' */

} // end namespace ug

#endif /* UG_PARALLEL */

#endif /* __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__FETI__ */
