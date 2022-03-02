/*
 * Copyright (c) 2010-2022:  G-CSC, Goethe University Frankfurt
 * Author: Lukas Larisch
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

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_SOLVER__FAST_MARCHING__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_SOLVER__FAST_MARCHING__
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


/*

	///	initializes this inverse operator for a matrix-based operator

	 * This method passes the operator A that is inverted by this operator. In
	 * addition some preparation step can be made.
	 *
	 * \param[in]	A		linear matrix-basewd operator to invert
	 * \returns		bool	success flag


		virtual bool init(SmartPtr<MatrixOperator<M,Y,X> > A) = 0;



	/// applies the inverse operator, i.e. returns u = A^{-1} * f
	 * This method applies the inverse operator.
	 *
	 * \param[out]	u		solution
	 * \param[in]	f		right-hand side
	 * \returns		bool	success flag


		virtual bool apply(Y& u, const X& f) = 0;



	/// applies the inverse operator and updates the defect
	 * This method applies the inverse operator and updates the defect, i.e.
	 * returns u = A^{-1} * f and in f the last defect d:= f - A*u is returned.
	 *
	 * \param[out]	u		solution
	 * \param[in]	f		right-hand side on entry, defect on exit
	 * \returns		bool	success flag


		virtual bool apply_return_defect(Y& u, X& f) = 0;



*/

template <typename TAlgebra>
class FastMarching
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

		using base_type::init;

	protected:
		using base_type::convergence_check;

	public:
	///	constructor
		FastMarching() : m_spOperator(NULL), m_mat(){};

	///	returns if parallel solving is supported
		virtual bool supports_parallel() const {return false;}

		const char* name() const {return "FastMarching";}

	public:
	///	initializes the solver for a matrix A
		bool init_fast_marching(const matrix_type *pA)
		{
			try{
		//	get matrix of Operator
			m_pMatrix = pA;
			if(m_pMatrix->num_rows() == 0) return true;

		//	check that matrix exist
			if(m_pMatrix == NULL)
			{
				UG_LOG("ERROR in '" << name() << "::init': No Matrix given.\n");
				return false;
			}

			const matrix_type &A = *pA;

			}UG_CATCH_THROW(name() << "::" << __FUNCTION__ << " failed")
			return true;
		}

		bool apply_fast_marching(vector_type &x, const vector_type &b)
		{
			try{

			if(m_pMatrix->num_rows() == 0) return true;

			//todo

			}UG_CATCH_THROW(name() << "::" << __FUNCTION__ << " failed")
			return true;
		}

	///	set operator L, that will be inverted
		virtual bool init(SmartPtr<MatrixOperator<matrix_type, vector_type> > Op)
		{
		// 	remember operator
			m_spOperator = Op;

		//	init LU operator
			if(!init_fast_marching(&m_spOperator->get_matrix()))
			{
				UG_LOG("ERROR in '" << name() << "::init': Cannot init LU Decomposition.\n");
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
			convergence_check()->set_name(name());

#ifdef UG_PARALLEL
			if(!f.has_storage_type(PST_ADDITIVE))
			{
				UG_LOG("ERROR: In '" << name() << "::apply': "
						"Inadequate storage format of Vector f.\n");
				return false;
			}
			if(!u.has_storage_type(PST_CONSISTENT))
			{
				UG_LOG("ERROR: In '" << name() << "::apply': "
						"Inadequate storage format of Vector u.\n");
				return false;
			}
#endif
			UG_ASSERT(f.size() == m_pMatrix->num_rows(), "Vector ["<<f.size()<<"] and Row "<<m_pMatrix->num_rows()<<" size mismatch");
			UG_ASSERT(u.size() == m_pMatrix->num_cols(), "Vector ["<<u.size()<<"] and Col "<<m_pMatrix->num_cols()<<" size mismatch");
			UG_ASSERT(f.size() == u.size(), "Vector sizes have to match!");

			if(!apply_fast_marching(u, f))
			{
				UG_LOG("ERROR in '" << name() << "::apply': "
						"Cannot apply FastMarching.\n");
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
		//	solve u
			if(!apply(u, f)) return false;

		//	update defect
			if(!m_pMatrix->matmul_minus(f, u))
			{
				UG_LOG("ERROR in '" << name() << "::apply_return_defect': "
						"Cannot apply matmul_minus.\n");
				return false;
			}

		//	we're done
			return true;
		}

		virtual std::string config_string() const
		{
			std::stringstream ss;
			ss << name() << ": Direct Solver for eikonal equations.\n";
			return ss.str();
		}

	///	Destructor
		virtual ~FastMarching() {};

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
};

} // end namespace ug

#endif
