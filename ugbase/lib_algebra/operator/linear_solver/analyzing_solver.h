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
#include "common/util/histogramm.h"

namespace ug{

void checksub(const CPUAlgebra::matrix_type &A);

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


			////
			// check symmetry


			/////


			const size_t nrOfRows = block_traits<typename matrix_type::value_type>::static_num_rows;
			size_t m_size = A.num_rows() * nrOfRows;

			CPUAlgebra::matrix_type mat;
			mat.resize_and_clear(m_size, m_size);

			for(size_t r=0; r<A.num_rows(); r++)
				for(typename matrix_type::const_row_iterator it = A.begin_row(r); it != A.end_row(r); ++it)
				{
					size_t rr = r*nrOfRows;
					size_t cc = it.index()*nrOfRows;
					for(size_t r2=0; r2<nrOfRows; r2++)
						for(size_t c2=0; c2<nrOfRows; c2++)
						{
							if(BlockRef(it.value(), r2, c2) != 0.0)
								mat(rr + r2, cc + c2) = BlockRef(it.value(), r2, c2);
						}
				}
			mat.defragment();
			checksub(mat);
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
