/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
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

#ifndef __H__LIB_ALGEBRA__LAPACK_LU_OPERATOR__
#define __H__LIB_ALGEBRA__LAPACK_LU_OPERATOR__
#include <iostream>
#include <sstream>

#include "common/common.h"
#include "lib_algebra/operator/interface/matrix_operator_inverse.h"

#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallelization.h"
#endif
#include "../preconditioner/ilut_scalar.h"
#include "../interface/preconditioned_linear_operator_inverse.h"
#include "linear_solver.h"

#include "lib_algebra/cpu_algebra_types.h"

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
		LU() : m_spOperator(NULL), m_mat(), m_bSortSparse(true), m_bInfo(false)
		{
#ifdef LAPACK_AVAILABLE
			m_iMinimumForSparse = 4000;
#else
			m_iMinimumForSparse = 1000;
#endif
		};

	///	returns if parallel solving is supported
		virtual bool supports_parallel() const {return false;}

	///
		void set_minimum_for_sparse(size_t N)
		{
			m_iMinimumForSparse=N;
		}

		void set_sort_sparse(bool b)
		{
			m_bSortSparse = b;
		}

		void set_info(bool b)
		{
			m_bInfo = b;
		}

	private:

		void print_info(const matrix_type &A)
		{
			size_t blockSize =  block_traits<typename matrix_type::value_type>::static_num_rows;
			UG_LOG("Matrix " << A.num_rows() << " x " << A.num_rows() << " with " << blockSize << " x " << blockSize << " blocks");
		}

	///	returns name of solver
		virtual const char* name() const {return "LU";}
		bool init_dense(const matrix_type &A)
		{
			PROFILE_FUNC();
			m_bDense = true;

			if(m_bInfo)
			{
				UG_LOG("LU, using DenseLU on ");
				print_info(A);
				UG_LOG("\n	DenseLU needs " << GetBytesSizeString(m_size*m_size*sizeof(double)) << " of memory.\n");
			}

			m_size = GetDenseDoubleFromSparse(m_mat, A);

			if(m_mat.invert() == false)
			{
				UG_THROW("ERROR in Matrix is singular");
				return false;
			}
			else
				return true;
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


		bool init_sparse(const matrix_type &A)
		{
			try{
			PROFILE_FUNC();
			m_bDense = false;

			if(m_bInfo)
			{
				UG_LOG("LU using Sparse LU on ");
				print_info(A);
				UG_LOG("\n");
			}
			ilut_scalar = make_sp(new ILUTScalarPreconditioner<algebra_type>(0.0));
			ilut_scalar->set_sort(m_bSortSparse);
			ilut_scalar->set_info(m_bInfo);
			ilut_scalar->preprocess(A);

			}UG_CATCH_THROW("LU::" << __FUNCTION__ << " failed")
			return true;
		}

		bool solve_dense(vector_type &x, const vector_type &b)
		{
			try{
			PROFILE_FUNC();
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

			}UG_CATCH_THROW("LU::" << __FUNCTION__ << " failed")
			return true;
		}

		bool solve_sparse(vector_type &x, const vector_type &b)
		{
			PROFILE_FUNC();
			ilut_scalar->solve(x, b);
			return true;
		}

	public:
	///	initializes the solver for a matrix A
		bool init_lu(const matrix_type *pA)
		{
			try{
		//	get matrix of Operator
			m_pMatrix = pA;
			if(m_pMatrix->num_rows() == 0) return true;

		//	check that matrix exist
			if(m_pMatrix == NULL)
			{
				UG_LOG("ERROR in 'LU::init': No Matrix given.\n");
				return false;
			}

			const matrix_type &A = *pA;

			PROFILE_FUNC();
			PROFILE_BEGIN_GROUP(LU_init_lu, "algebra lu");
			if(block_traits<typename vector_type::value_type>::is_static)
			{
				const size_t nrOfRows = block_traits<typename matrix_type::value_type>::static_num_rows;
				UG_ASSERT(nrOfRows == block_traits<typename matrix_type::value_type>::static_num_cols, "only square matrices supported");
				m_size = A.num_rows() * nrOfRows;

				if(m_size > m_iMinimumForSparse)
					init_sparse(A);
				else
					init_dense(A);
			}
			else
				init_var(A);
			}UG_CATCH_THROW("LU::" << __FUNCTION__ << " failed")
			return true;
		}

		bool apply_lu(vector_type &x, const vector_type &b)
		{
			if(m_pMatrix->num_rows() == 0) return true;
			try{
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

			}UG_CATCH_THROW("LU::" << __FUNCTION__ << " failed")
		}

	///	set operator L, that will be inverted
		virtual bool init(SmartPtr<MatrixOperator<matrix_type, vector_type> > Op)
		{
		// 	remember operator
			m_spOperator = Op;

		//	init LU operator
			if(!init_lu(&m_spOperator->get_matrix()))
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
			PROFILE_FUNC();
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

		virtual std::string config_string() const
		{
			std::stringstream ss;
			ss << "LU Decomposition: Direct Solver for Linear Equation Systems.\n";
			ss << " Minimum Entries for Sparse LU: " << m_iMinimumForSparse;
			if(m_iMinimumForSparse==0)
				ss << " (= always Sparse LU)";
			return ss.str();
		}


	///	Destructor
		virtual ~LU() {};

	protected:
	/// Operator to invert
		SmartPtr<MatrixOperator<matrix_type, vector_type> > m_spOperator;

	/// matrix to invert
		const matrix_type* m_pMatrix;

	/// inverse
		DenseMatrixInverse<DenseMatrix<VariableArray2<double> > > m_mat;
		DenseVector<VariableArray1<double> > m_tmp;
		CPUAlgebra::vector_type m_u;
		CPUAlgebra::vector_type m_b;
		size_t m_size;

		bool m_bDense;
		SmartPtr<ILUTScalarPreconditioner<algebra_type> > ilut_scalar;
		size_t m_iMinimumForSparse;
		bool m_bSortSparse, m_bInfo;
};

} // end namespace ug

#endif /* __H__LIB_ALGEBRA__LAPACK_LU_OPERATOR__ */
