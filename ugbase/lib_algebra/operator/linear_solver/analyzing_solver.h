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

		void checksub(const CPUAlgebra::matrix_type &A)
		{
			UG_LOG(reset_floats);
			size_t N = A.num_rows();
			// check isolated
			size_t iIsolated = 0;
			for(size_t r=0; r<A.num_rows(); r++)
				if(A.is_isolated(r)) iIsolated++;

			typedef typename CPUAlgebra::matrix_type::const_row_iterator row_it;
			typedef typename CPUAlgebra::matrix_type::value_type value_type;

			UG_LOG("Nr of dirichlet nodes: " << iIsolated << " (" << iIsolated*100.0/N << "% )\n");

			// check symmetric

			std::vector<double> alpha(N, 0);

			size_t iUnsymmetric=0;
			const double unsymmetricEps = 1e-8;
			for(size_t r=0; r<A.num_rows(); r++)
			{
				double dUnsymmetric = 0;
				if(A.is_isolated(r) ) continue;

				for(row_it it = A.begin_row(r); it != A.end_row(r); ++it)
				{
					size_t c = it.index();
					if(A.is_isolated(c) ) continue;
					double T = A(c, r);
					if(A(c, r) != MatrixTranspose(it.value()))
					{
						dUnsymmetric += dabs(T-it.value());
					}

				}
				double diag = A(r, r);
				if(diag == 0.0) continue;

				dUnsymmetric /= diag;

				alpha[r] = dUnsymmetric;
				if(dUnsymmetric > unsymmetricEps)
					iUnsymmetric++;

			}
			std::sort(alpha.begin(), alpha.end());

			// check sign condition

			if(iUnsymmetric==0)
			{	UG_LOG("Matrix is symmetric! (maximum alpha = " << alpha[N-1] << ")\n"); }
			else
			{
				UG_LOG("Matrix is unsymmetric in " << (iUnsymmetric*100)/(N-iIsolated) << "% of the non-dirchlet rows (" << iUnsymmetric << " total)\n");
				UG_LOG(" row i alpha-unsymmetric means: sum_{A_{ij} != 0} |A_{ij}-A{ji}| / A_{ii} >= alpha\n")

				UG_LOG("> alpha distribution:\n");
				UG_LOG(DistributionPercentage(alpha));
			}
			size_t signConditionMet=0;
			size_t zeroDiagonal=0;

			double minEW = 1e20;
			double maxEW = 1e-20;

			for(size_t r=0; r<A.num_rows(); r++)
			{
				if(A(r, r)==0.0)
				{
					zeroDiagonal++;
					continue;
				}
				bool bPos = A(r, r) > 0;
				double s=0.0;
				bool bSignCondMet = true;
				for(row_it it = A.begin_row(r); it != A.end_row(r); ++it)
				{
					if(it.index() == r) continue;
					if(bPos)
					{
						if(it.value() > 0)
							bSignCondMet = false;
					}
					else
					{
						if(it.value() < 0)
							bSignCondMet = false;
					}
					s += dabs(it.value());
				}
				if(bSignCondMet && A.is_isolated(r) == false)
					signConditionMet++;

				minEW = std::min(minEW, A(r,r)-s);
				maxEW = std::max(maxEW, A(r,r)+s);
			}

			if(signConditionMet == N-iIsolated)
			{	UG_LOG("Sign condition met in all nodes\n"); }
			else
			{
				UG_LOG("Sign condition met in " << (signConditionMet*100.0)/(N-iIsolated) << "% of the non-dirichlet rows (" << signConditionMet << " total)\n");
			}
			UG_LOG("Gershgorin Eigenvalues are within [" << minEW << ", " << maxEW << "]\n");
		}


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
