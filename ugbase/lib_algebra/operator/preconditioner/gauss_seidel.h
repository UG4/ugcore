/*
 * gauss_seidel.h
 *
 *  Created on: 14.07.2010
 *      Author: Martin Rupp
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__GAUSS_SEIDEL__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__GAUSS_SEIDEL__

#include "lib_algebra/operator/interface/preconditioner.h"
#include "lib_algebra/algebra_common/core_smoothers.h"
#include "lib_algebra/algebra_common/sparsematrix_util.h"
#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallelization.h"
	#include "lib_algebra/parallelization/parallel_matrix_overlap_impl.h"
#endif

namespace ug{

/// Gauss-Seidel preconditioner for the 'forward' ordering of the dofs
/**
 * This class implements the Gauss-Seidel preconditioner (and smoother) for the
 * 'forward' ordering of the dofs. When a relaxation parameter is set by the method
 * 'set_sor_relax', the resulting preconditioner is better known as (forward) 'SOR'-method.
 * References:
 * <ul>
 * <li> W. Hackbusch. Iterative solution of large sparse systems of equations. New York: Springer, 1994
 * </ul>
 *
 * \tparam	TAlgebra	Algebra type
 */
template <typename TAlgebra>
class GaussSeidel : public IPreconditioner<TAlgebra>
{
	public:
	//	Algebra type
		typedef TAlgebra algebra_type;

	//	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	//	Matrix type
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
		GaussSeidel() { m_relax = 1.0; };

	// 	Clone
		virtual SmartPtr<ILinearIterator<vector_type> > clone()
		{
			SmartPtr<GaussSeidel<algebra_type> > newInst(new GaussSeidel<algebra_type>());
			newInst->set_debug(debug_writer());
			newInst->set_damp(this->damping());
			return newInst;
		}

	//	set relaxation parameter to define a SOR-method
		void set_sor_relax(number relaxFactor){ m_relax = relaxFactor;}

	protected:
	//	Name of preconditioner
		virtual const char* name() const {return "Gauss-Seidel";}

	//	Preprocess routine
		virtual bool preprocess(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp)
		{
			PROFILE_BEGIN_GROUP(GaussSeidel_preprocess, "algebra gaussseidel");
#ifdef UG_PARALLEL
			if(pcl::GetNumProcesses() > 1)
			{
				//	copy original matrix
				MakeConsistent(*pOp, m_A);
				//	set zero on slaves
				std::vector<IndexLayout::Element> vIndex;
				CollectUniqueElements(vIndex,  m_A.layouts()->slave());
				SetDirichletRow(m_A, vIndex);
			}
			matrix_type &A = m_A;
#else
			matrix_type &A = *pOp;
#endif
//			UG_ASSERT(CheckDiagonalInvertible(A), "GS: A has noninvertible diagonal");
			CheckDiagonalInvertible(A);
			return true;
		}

	//	Postprocess routine
		virtual bool postprocess() {return true;}

	//	Stepping routine
		virtual bool step(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp, vector_type& c, const vector_type& d)
		{
			PROFILE_BEGIN_GROUP(GaussSeidel_step, "algebra gaussseidel");
#ifdef UG_PARALLEL
			if(pcl::GetNumProcesses() > 1)
			{
			//	make defect unique
				SmartPtr<vector_type> spDtmp = d.clone();
				spDtmp->change_storage_type(PST_UNIQUE);

				gs_step_LL(m_A, c, *spDtmp, m_relax);
				c.set_storage_type(PST_UNIQUE);
				return true;
			}
			else
#endif
			{
				gs_step_LL(*pOp, c, d, m_relax);
#ifdef UG_PARALLEL
				c.set_storage_type(PST_UNIQUE);
#endif
				return true;
			}
		}

	protected:
#ifdef UG_PARALLEL
		matrix_type m_A;
#endif

	private:
		//	relaxation parameter
		number m_relax;
};

/// Gauss-Seidel preconditioner for the 'backward' ordering of the dofs
/**
 * This class implements the Gauss-Seidel preconditioner (and smoother) for the
 * 'backward' ordering of the dofs. When a relaxation parameter is set by the method
 * 'set_sor_relax', the resulting preconditioner is better known as backward 'SOR'-method.
 * References:
 * <ul>
 * <li> W. Hackbusch. Iterative solution of large sparse systems of equations. New York: Springer, 1994
 * </ul>
 *
 * \tparam	TAlgebra	Algebra type
 */
template <typename TAlgebra>
class BackwardGaussSeidel : public IPreconditioner<TAlgebra>
{
	public:
	//	Algebra type
		typedef TAlgebra algebra_type;

	//	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	//	Matrix type
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
		BackwardGaussSeidel() { m_relax = 1.0; };

	// 	Clone
		virtual SmartPtr<ILinearIterator<vector_type> > clone()
		{
			SmartPtr<BackwardGaussSeidel<algebra_type> > newInst(new BackwardGaussSeidel<algebra_type>());
			newInst->set_debug(debug_writer());
			newInst->set_damp(this->damping());
			return newInst;
		}


	//	set relaxation parameter to define a SOR-method
		void set_sor_relax(number relaxFactor){ m_relax = relaxFactor;}

	protected:
	//	Name of preconditioner
		virtual const char* name() const {return "Backward Gauss-Seidel";}

	//	Preprocess routine
		virtual bool preprocess(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp)
		{
			PROFILE_BEGIN_GROUP(BackwardGaussSeidel_preprocess, "algebra BackwardGaussSeidel");
#ifdef UG_PARALLEL
			if(pcl::GetNumProcesses() > 1)
			{
				//	copy original matrix
				MakeConsistent(*pOp, m_A);
				//	set zero on slaves
				std::vector<IndexLayout::Element> vIndex;
				CollectUniqueElements(vIndex,  m_A.layouts()->slave());
				SetDirichletRow(m_A, vIndex);
			}
			matrix_type &A = m_A;
#else
			matrix_type &A = *pOp;
#endif
//			UG_ASSERT(CheckDiagonalInvertible(A), "BackwardGS: A has noninvertible diagonal");
			CheckDiagonalInvertible(A);
			return true;
		}

	//	Postprocess routine
		virtual bool postprocess() {return true;}

	//	Stepping routine
		virtual bool step(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp, vector_type& c, const vector_type& d)
		{
			PROFILE_BEGIN_GROUP(GaussSeidel_step, "algebra BackwardGaussSeidel");
#ifdef UG_PARALLEL
			if(pcl::GetNumProcesses() > 1)
			{
				//	make defect unique
				// todo: change that copying
				vector_type dhelp;
				dhelp.resize(d.size()); dhelp = d;
				dhelp.change_storage_type(PST_UNIQUE);

				gs_step_UR(m_A, c, dhelp, m_relax);
				c.set_storage_type(PST_UNIQUE);
				return true;
			}
			else
#endif
			{
				gs_step_UR(*pOp, c, d, m_relax);
#ifdef UG_PARALLEL
				c.set_storage_type(PST_UNIQUE);
#endif
				return true;
			}
		}

	protected:
#ifdef UG_PARALLEL
		matrix_type m_A;
#endif

	private:
		//	relaxation parameter
		number m_relax;
};

/// Symmetric Gauss-Seidel preconditioner
/**
 * This class implements the symmetric Gauss-Seidel preconditioner (and smoother).
 * When a relaxation parameter is set by the method 'set_sor_relax', the resulting
 * preconditioner is better known as 'SSOR'-method.
 *
 * References:
 * <ul>
 * <li> W. Hackbusch. Iterative solution of large sparse systems of equations. New York: Springer, 1994
 * </ul>
 *
 * \tparam	TAlgebra	Algebra type
 */
template <typename TAlgebra>
class SymmetricGaussSeidel : public IPreconditioner<TAlgebra>
{
	public:
	//	Algebra type
		typedef TAlgebra algebra_type;

	//	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	//	Matrix type
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
		SymmetricGaussSeidel() { m_relax = 1.0; };

	// 	Clone
		virtual SmartPtr<ILinearIterator<vector_type> > clone()
		{
			SmartPtr<SymmetricGaussSeidel<algebra_type> > newInst(new SymmetricGaussSeidel<algebra_type>());
			newInst->set_debug(debug_writer());
			newInst->set_damp(this->damping());
			return newInst;
		}

	//	set relaxation parameter to define a SSOR-method
		void set_sor_relax(number relaxFactor){ m_relax = relaxFactor;}

	protected:
	//	Name of preconditioner
		virtual const char* name() const {return "Symmetric Gauss-Seidel";}

	//	Preprocess routine
		virtual bool preprocess(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp)
		{
			PROFILE_BEGIN_GROUP(SymmetricGaussSeidel_preprocess, "algebra SymmetricGaussSeidel");
#ifdef UG_PARALLEL
			if(pcl::GetNumProcesses() > 1)
			{
				//	copy original matrix
				MakeConsistent(*pOp, m_A);
				//	set zero on slaves
				std::vector<IndexLayout::Element> vIndex;
				CollectUniqueElements(vIndex,  m_A.layouts()->slave());
				SetDirichletRow(m_A, vIndex);
			}
			matrix_type &A = m_A;
#else
			matrix_type &A = *pOp;
#endif
			//UG_ASSERT(CheckDiagonalInvertible(A), "SGS: A has noninvertible diagonal");
			CheckDiagonalInvertible(A);
			return true;
		}

	//	Postprocess routine
		virtual bool postprocess() {return true;}

	//	Stepping routine
		virtual bool step(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp, vector_type& c, const vector_type& d)
		{
			PROFILE_BEGIN_GROUP(SymmetricGaussSeidel_step, "algebra SymmetricGaussSeidel");
#ifdef UG_PARALLEL
			if(pcl::GetNumProcesses() > 1)
			{
				//	make defect unique
				// todo: change that copying
				vector_type dhelp;
				dhelp.resize(d.size()); dhelp = d;
				dhelp.change_storage_type(PST_UNIQUE);

				sgs_step(m_A, c, dhelp, m_relax);
				c.set_storage_type(PST_UNIQUE);
				return true;
			}
			else
#endif
			{
				sgs_step(*pOp, c, d, m_relax);
#ifdef UG_PARALLEL
				c.set_storage_type(PST_UNIQUE);
#endif
				return true;
			}
		}

	protected:
#ifdef UG_PARALLEL
		matrix_type m_A;
#endif

	private:
		//	relaxation parameter
		number m_relax;
};


} // end namespace ug

#endif // __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__GAUSS_SEIDEL__
