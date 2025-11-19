/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
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
		using domain_function_type = X;

	///	Range space
		using codomain_function_type = Y;

	///	Matrix type
		using matrix_type = M;

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

		virtual std::string config_string() const
		{
			return "AnalyzingSolver " + m_pLinearOperatorInverse->config_string();
		}
	private:
		SmartPtr<ILinearOperatorInverse<X,Y> > m_pLinearOperatorInverse;
};


} // end namespace ug

#endif