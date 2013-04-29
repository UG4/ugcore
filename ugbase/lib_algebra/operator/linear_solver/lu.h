/*
 * lu.h
 *
 *  Created on: 16.06.2010
 *      Author: mrupp
 */

#ifndef __H__LIB_ALGEBRA__LAPACK_LU_OPERATOR__
#define __H__LIB_ALGEBRA__LAPACK_LU_OPERATOR__
#include <iostream>
#include <sstream>

#include "common/common.h"
#include "lib_algebra/operator/interface/matrix_operator_inverse.h"

#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallelization.h"
#endif
#include "../preconditioner/ilut.h"
#include "../interface/preconditioned_linear_operator_inverse.h"
#include "linear_solver.h"

#include "../operator_util.h"
namespace ug{

template <typename TAlgebra>
class LU
	: public IMatrixOperatorInverse<typename TAlgebra::matrix_type,
	  	  	  	  	  	  	  	    typename TAlgebra::vector_type>
{
	public:
	///	Algebra type
		typedef TAlgebra algebra_type;

	///	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	///	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

	///	Base type
		typedef IMatrixOperatorInverse<matrix_type,vector_type> base_type;

	protected:
		using base_type::convergence_check;

	public:
	///	constructor
		LU() : m_spOperator(NULL), m_mat() {};

	///	returns name of solver
		virtual const char* name() const {return "LU";}
		bool init_dense(const matrix_type &A)
		{
			m_bDense = true;
			const size_t nrOfRows = block_traits<typename matrix_type::value_type>::static_num_rows;
			m_mat.resize(m_size);
			for(size_t r=0; r<A.num_rows(); r++)
				for(typename matrix_type::const_row_iterator it = A.begin_row(r); it != A.end_row(r); ++it)
				{
					size_t rr = r*nrOfRows;
					size_t cc = it.index()*nrOfRows;
					for(size_t r2=0; r2<nrOfRows; r2++)
							for(size_t c2=0; c2<nrOfRows; c2++)
							  m_mat(rr + r2, cc + c2) = BlockRef(it.value(), r2, c2);
				}

			return m_mat.invert();
		}

		void init_var(const matrix_type &A)
		{
			UG_ASSERT(0, "not yet tested");
			/*m_size = 0;
			std::vector<size_t> blockbegin(A.num_rows()+1);

			for(size_t i=0; i<A.num_rows(); i++)
			{
				bool bFound;
				typename matrix_type::const_row_iterator it = A.get_connection(i,i, bFound);
				UG_ASSERT(bFound, "Matrix has to have entry A(" << i << ", " << i << ")");
				size_t s = GetRows(it.value());
				UG_ASSERT(s == GetCols(it.value()), "diagonal elements have to be square");
				if(i == 0)
				blockbegin[i] = m_size;
				m_size += s;
			}
			blockbegin[A.num_rows()] = m_size;

			m_mat.resize(m_size);

			for(size_t r=0; r<A.num_rows(); r++)
				for(typename matrix_type::const_row_iterator it = A.begin_row(r); it != A.end_row(r); ++it)
				{
					size_t c = it.index();
					const typename matrix_type::value_type &val = it.value();
					UG_ASSERT(blockbegin[r]+GetRows(val) == blockbegin[r+1], "blocksizes in matrix inconsistent");
					UG_ASSERT(blockbegin[c]+GetCols(val) == blockbegin[c+1], "blocksizes in matrix inconsistent");
					for(size_t r2=0; r2 < GetRows(val); r2++)
						for(size_t c2=0; c2 < GetCols(val); c2++)
							m_mat(blockbegin[r] + r2, blockbegin[c]+c2) = BlockRef(val, r2, c2);
				}
*/
		}


		bool m_bDense;
		SmartPtr<ILUTPreconditioner<CRSAlgebra> > ilut;
		SmartPtr<MatrixOperator<CRSAlgebra::matrix_type, CRSAlgebra::vector_type> > mo;
		SmartPtr<IPreconditionedLinearOperatorInverse<CRSAlgebra::vector_type> > linearSolver;
		bool init_sparse(const matrix_type &A)
		{
			m_bDense = false;
			ilut = new ILUTPreconditioner<CRSAlgebra>(0);

			mo = new MatrixOperator<CRSAlgebra::matrix_type, CRSAlgebra::vector_type>;
			CRSAlgebra::matrix_type &mat = mo->get_matrix();
			mat.resize(m_size, m_size);
#ifdef UG_PARALLEL
			mat.set_storage_type(PST_ADDITIVE);
#endif

			const size_t nrOfRows = block_traits<typename matrix_type::value_type>::static_num_rows;
			for(size_t r=0; r<A.num_rows(); r++)
				for(typename matrix_type::const_row_iterator it = A.begin_row(r); it != A.end_row(r); ++it)
				{
					size_t rr = r*nrOfRows;
					size_t cc = it.index()*nrOfRows;
					for(size_t r2=0; r2<nrOfRows; r2++)
						for(size_t c2=0; c2<nrOfRows; c2++)
							mat(rr + r2, cc + c2) = BlockRef(it.value(), r2, c2);
				}
			mat.defragment();

			SmartPtr<StdConvCheck<CRSAlgebra::vector_type> > convCheck = new StdConvCheck<CRSAlgebra::vector_type>;
			convCheck->set_maximum_steps(100);
			convCheck->set_minimum_defect(1e-12);
			convCheck->set_reduction(1e-16);
			convCheck->set_verbose(false);



			linearSolver = new LinearSolver<CRSAlgebra::vector_type>;
			linearSolver->set_preconditioner(ilut);
			linearSolver->set_convergence_check(convCheck);
			linearSolver->init(mo);
			return true;
		}

		bool solve_dense(vector_type &x, const vector_type &b)
		{
			x = b;
			m_tmp.resize(m_size);
			for(size_t i=0, k=0; i<b.size(); i++)
			{
				for(size_t j=0; j<GetSize(b[i]); j++)
					m_tmp[k++] = BlockRef(b[i],j);
			}
			m_mat.apply(m_tmp);


			for(size_t i=0, k=0; i<b.size(); i++)
			{
				for(size_t j=0; j<GetSize(b[i]); j++)
					BlockRef(x[i],j) = m_tmp[k++];
			}

			return true;
		}

		bool solve_sparse(vector_type &x, const vector_type &b)
		{
			m_u.resize(m_size);
			m_b.resize(m_size);

#ifdef UG_PARALLEL
			m_b.set_storage_type(PST_ADDITIVE);
			m_u.set_storage_type(PST_CONSISTENT);
#endif
			for(size_t i=0, k=0; i<b.size(); i++)
			{
				for(size_t j=0; j<GetSize(b[i]); j++)
					m_b[k++] = BlockRef(b[i],j);
			}
			for(size_t i=0; i<m_u.size(); i++) m_u[i] = 0.0;




			//ApplyLinearSolver(mo, m_u, m_b, linearSolver);
			linearSolver->apply_return_defect(m_u,m_b);

			for(size_t i=0, k=0; i<b.size(); i++)
			{
				for(size_t j=0; j<GetSize(b[i]); j++)
					BlockRef(x[i],j) = m_u[k++];
			}
			return true;
		}

	///	initializes the solver for a matrix A
		bool init_lu(const matrix_type &A)
		{
			PROFILE_BEGIN_GROUP(LU_init_lu, "algebra lu");
			if(block_traits<typename vector_type::value_type>::is_static)
			{
				const size_t nrOfRows = block_traits<typename matrix_type::value_type>::static_num_rows;
				UG_ASSERT(nrOfRows == block_traits<typename matrix_type::value_type>::static_num_cols, "only square matrices supported");
				m_size = A.num_rows() * nrOfRows;

				if(m_size > 4000)
					init_sparse(A);
				else
					init_dense(A);
			}
			else
				init_var(A);
			return true;
		}

		bool apply_lu(vector_type &x, const vector_type &b)
		{
			PROFILE_BEGIN_GROUP(LU_apply_lu, "algebra lu");
#ifndef NDEBUG
			if(block_traits<typename vector_type::value_type>::is_static)
			{
				const size_t static_size = block_traits<typename vector_type::value_type>::static_size;
				UG_ASSERT(m_size == b.size() * static_size && m_size == x.size() * static_size,
						" wrong size! has to be " << m_size << ", but is " << b << " and " << x);
			}
			else
			{
				size_t b_size = 0;
				for(size_t i=0; i<b.size(); i++)
				{
					UG_ASSERT(GetSize(b[i]) == GetSize(x[i]), "wrong size! Sizes of b and x must be the same, but is "
							<< GetSize(b[i]) << " and " << GetSize(x[i]) << "!");
					b_size += GetSize(b[i]);
				}
				UG_ASSERT(m_size == b_size, " wrong size! has to be " << m_size << ", but is " << b_size << "!");
			}
#endif

			if(m_bDense)
				return solve_dense(x, b);
			else
				return solve_sparse(x, b);
		}

	///	set operator L, that will be inverted
		virtual bool init(SmartPtr<MatrixOperator<matrix_type, vector_type> > Op)
		{
		// 	remember operator
			m_spOperator = Op;

		//	get matrix of Operator
			m_pMatrix = &m_spOperator->get_matrix();

		//	check that matrix exist
			if(m_pMatrix == NULL)
			{
				UG_LOG("ERROR in 'LU::init': No Matrix given.\n");
				return false;
			}

		//	init LU operator
			if(!init_lu(*m_pMatrix))
			{
				UG_LOG("ERROR in 'LU::init': Cannot init LU Decomposition.\n");
				return false;
			}

		//	we're done
			return true;
		}

	///	Compute u = L^{-1} * f
		virtual bool apply(vector_type& u, const vector_type& f)
		{
			convergence_check()->set_symbol('%');
			convergence_check()->set_name("LU Solver");

#ifdef UG_PARALLEL
			if(!f.has_storage_type(PST_ADDITIVE))
			{
				UG_LOG("ERROR: In 'LU::apply': "
						"Inadequate storage format of Vector f.\n");
				return false;
			}
			if(!u.has_storage_type(PST_CONSISTENT))
			{
				UG_LOG("ERROR: In 'LU::apply': "
						"Inadequate storage format of Vector u.\n");
				return false;
			}
#endif
			UG_ASSERT(f.size() == m_pMatrix->num_rows(), "Vector ["<<f.size()<<"] and Row "<<m_pMatrix->num_rows()<<" size mismatch");
			UG_ASSERT(u.size() == m_pMatrix->num_cols(), "Vector ["<<u.size()<<"] and Col "<<m_pMatrix->num_cols()<<" size mismatch");
			UG_ASSERT(f.size() == u.size(), "Vector sizes have to match!");

			if(!apply_lu(u, f))
			{
				UG_LOG("ERROR in 'LU::apply': "
						"Cannot apply LU decomposition.\n");
				return false;
			}

#ifdef UG_PARALLEL
			// todo: we set solution to consistent here, but that is only true for
			//			serial case. Handle parallel case.
			u.set_storage_type(PST_CONSISTENT);
#endif

		//	we're done
			return true;
		}

	/// Compute u = L^{-1} * f AND return defect f := f - L*u
		virtual bool apply_return_defect(vector_type& u, vector_type& f)
		{
			PROFILE_BEGIN_GROUP(LU_apply_return_defect, "algebra lu");
		//	solve u
			if(!apply(u, f)) return false;

		//	update defect
			if(!m_pMatrix->matmul_minus(f, u))
			{
				UG_LOG("ERROR in 'LU::apply_return_defect': "
						"Cannot apply matmul_minus.\n");
				return false;
			}

		//	we're done
			return true;
		}

	///	Destructor
		virtual ~LU() {};

	protected:
	/// Operator to invert
		SmartPtr<MatrixOperator<matrix_type, vector_type> > m_spOperator;

	/// matrix to invert
		matrix_type* m_pMatrix;

	/// inverse
		DenseMatrixInverse<DenseMatrix<VariableArray2<double> > > m_mat;
		DenseVector<VariableArray1<double> > m_tmp;
		CRSAlgebra::vector_type m_u;
		CRSAlgebra::vector_type m_b;
		size_t m_size;
};

} // end namespace ug

#endif /* __H__LIB_ALGEBRA__LAPACK_LU_OPERATOR__ */
