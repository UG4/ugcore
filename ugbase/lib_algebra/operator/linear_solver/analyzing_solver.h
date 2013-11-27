/*
 * analyzing solver
 *
 *  Created on: 19.11.2013
 *      Author: mrupp
 */

#ifndef __H__LIB_ALGEBRA__LAPACK_ANALYZING_SOLVER__
#define __H__LIB_ALGEBRA__LAPACK_ANALYZING_SOLVER__

#include "../interface/linear_operator_inverse.h"
#include "../interface/matrix_operator.h"
#include "common/error.h"
#include "common/util/smart_pointer.h"

namespace ug{

template <typename M, typename X, typename Y = X>
class AnalyzingSolver
	: public virtual ILinearOperatorInverse<X,Y>
{
	public:
	///	Domain space
		typedef X domain_function_type;

	///	Range space
		typedef Y codomain_function_type;

	///	Matrix type
		typedef M matrix_type;

	public:
		virtual bool apply(Y& u, const X& f)
		{
			return m_pLinearOperatorInverse->apply(u, f);
		}

		virtual bool apply_return_defect(Y& u, X& f)
		{
			return m_pLinearOperatorInverse->apply_return_defect(u, f);
		}

		AnalyzingSolver(SmartPtr<ILinearOperatorInverse<X,Y> > pLinearOperatorInverse)
		{
			m_pLinearOperatorInverse = pLinearOperatorInverse;
		}

	/// virtual destructor
		virtual ~AnalyzingSolver() {};

	///	returns if parallel solving is supported
		virtual bool supports_parallel() const
		{
			return m_pLinearOperatorInverse->supports_parallel();
		}

	public:
		void check(const matrix_type &A)
		{
			UG_LOG("ANALYZING SOLVER:\n");
			UG_LOG(" Matrix is of dimension " << A.num_rows() << "\n");
			if(A.num_rows() != A.num_cols())
			{	UG_LOG(" Matrix is not quadratic???\n");	}
			if(GetRows(A(0,0)) == 1)
			{ UG_LOG(" Submatrices are DOUBLE.\n")}
			else
			{UG_LOG(" Submatrices are Matrices of size " << GetRows(A(0,0)) << " x " << GetCols(A(0,0)) << "\n");}
		}
		virtual bool init(SmartPtr<ILinearOperator<Y,X> > A, const Y&u)
		{
			check(A);
			return m_pLinearOperatorInverse->init(A, u);
		}

		virtual bool init(SmartPtr<ILinearOperator<Y,X> > A)
		{
			check(A);
			return m_pLinearOperatorInverse->init(A);
		}

		void check(SmartPtr<ILinearOperator<Y,X> > A)
		{
		//	cast operator
			SmartPtr<MatrixOperator<M,Y,X> > op =
									A.template cast_dynamic<MatrixOperator<M,Y,X> >();

		//	check if correct types are present
			if(op.invalid())
				UG_THROW("IMatrixOperatorInverse::init:"
						" Passed operator is not matrix-based.");

		//	forward request
			check(*op);
		}
		virtual const char *name() const { return m_pLinearOperatorInverse->name(); }
	private:
		SmartPtr<ILinearOperatorInverse<X,Y> > m_pLinearOperatorInverse;
};


} // end namespace ug

#endif /* __H__LIB_ALGEBRA__LAPACK_LU_OPERATOR__ */
