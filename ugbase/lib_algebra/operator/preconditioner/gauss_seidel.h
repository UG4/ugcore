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

template<typename TAlgebra>
class GaussSeidelBase : public IPreconditioner<TAlgebra>
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
		GaussSeidelBase() { m_relax = 1.0; };

	//	set relaxation parameter to define a SOR-method
		void set_sor_relax(number relaxFactor){ m_relax = relaxFactor;}

		virtual const char* name() const = 0;
	protected:
		void copy_config(GaussSeidelBase &newInst)
		{
			newInst.set_debug(base_type::debug_writer());
			newInst.set_damp(base_type::damping());
			newInst.set_sor_relax(m_relax);
		}

		virtual bool supports_parallel() const {return true;}

	//	Preprocess routine
		virtual bool preprocess(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp)
		{
			PROFILE_BEGIN_GROUP(GaussSeidel_preprocess, "algebra gaussseidel");
			matrix_type *pA;
#ifdef UG_PARALLEL
			if(pcl::NumProcs() > 1)
			{
				//	copy original matrix
				MakeConsistent(*pOp, m_A);
				//	set zero on slaves
				std::vector<IndexLayout::Element> vIndex;
				CollectUniqueElements(vIndex,  m_A.layouts()->slave());
				SetDirichletRow(m_A, vIndex);
				pA = &m_A;
			}
			else
				pA = &(*pOp);
#else
			pA = &(*pOp);
#endif
			THROW_IF_NOT_EQUAL(pA->num_rows(), pA->num_cols());
//			UG_ASSERT(CheckDiagonalInvertible(A), "GS: A has noninvertible diagonal");
			UG_COND_THROW(CheckDiagonalInvertible(*pA) == false, name() << ": A has noninvertible diagonal");
			return true;
		}

	//	Postprocess routine
		virtual bool postprocess() {return true;}

		virtual void step(const matrix_type &A, vector_type &c, const vector_type &d, const number relax) = 0;

	//	Stepping routine
		virtual bool step(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp, vector_type& c, const vector_type& d)
		{
			PROFILE_BEGIN_GROUP(GaussSeidel_step, "algebra gaussseidel");

#ifdef UG_PARALLEL
			if(pcl::NumProcs() > 1)
			{
				// todo: do not clone every time
			//	make defect unique
				SmartPtr<vector_type> spDtmp = d.clone();
				spDtmp->change_storage_type(PST_UNIQUE);

				THROW_IF_NOT_EQUAL_3(c.size(), spDtmp->size(), m_A.num_rows());
				step(m_A, c, *spDtmp, m_relax);
				c.set_storage_type(PST_UNIQUE);
				return true;
			}
			else
#endif
			{
				matrix_type &A = *pOp;
				THROW_IF_NOT_EQUAL_4(c.size(), d.size(), A.num_rows(), A.num_cols());
				step(A, c, d, m_relax);
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
class GaussSeidel : public GaussSeidelBase<TAlgebra>
{
	typedef TAlgebra algebra_type;
	typedef typename TAlgebra::vector_type vector_type;
	typedef typename TAlgebra::matrix_type matrix_type;
	typedef GaussSeidelBase<TAlgebra> base_type;

	protected:
	//	Name of preconditioner
		virtual const char* name() const {return "Gauss-Seidel";}

	// 	Clone
		virtual SmartPtr<ILinearIterator<vector_type> > clone()
		{
			SmartPtr<GaussSeidel<algebra_type> > newInst(new GaussSeidel<algebra_type>());
			base_type::copy_config(*newInst);
			return newInst;
		}


	//	Stepping routine
		virtual void step(const matrix_type &A, vector_type &c, const vector_type &d, const number relax)
		{
			gs_step_LL(A, c, d, relax);
		}
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
class BackwardGaussSeidel : public GaussSeidelBase<TAlgebra>
{
	typedef TAlgebra algebra_type;
	typedef typename TAlgebra::vector_type vector_type;
	typedef typename TAlgebra::matrix_type matrix_type;
	typedef GaussSeidelBase<TAlgebra> base_type;

	protected:
	//	Name of preconditioner
		virtual const char* name() const {return "Backward Gauss-Seidel";}

	// 	Clone
		virtual SmartPtr<ILinearIterator<vector_type> > clone()
		{
			SmartPtr<BackwardGaussSeidel<algebra_type> > newInst(new BackwardGaussSeidel<algebra_type>());
			base_type::copy_config(*newInst);
			return newInst;
		}

	//	Stepping routine
		virtual void step(const matrix_type &A, vector_type &c, const vector_type &d, const number relax)
		{
			gs_step_UR(A, c, d, relax);
		}
};


template <typename TAlgebra>
class SymmetricGaussSeidel : public GaussSeidelBase<TAlgebra>
{
	typedef TAlgebra algebra_type;
	typedef typename TAlgebra::vector_type vector_type;
	typedef typename TAlgebra::matrix_type matrix_type;
	typedef GaussSeidelBase<TAlgebra> base_type;

	protected:
	//	Name of preconditioner
		virtual const char* name() const {return "Symmetric Gauss-Seidel";}

	// 	Clone
		virtual SmartPtr<ILinearIterator<vector_type> > clone()
		{
			SmartPtr<SymmetricGaussSeidel<algebra_type> > newInst(new SymmetricGaussSeidel<algebra_type>());
			base_type::copy_config(*newInst);
			return newInst;
		}

	//	Stepping routine
		virtual void step(const matrix_type &A, vector_type &c, const vector_type &d, const number relax)
		{
			sgs_step(A, c, d, relax);
		}
};

} // end namespace ug

#endif // __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__GAUSS_SEIDEL__
