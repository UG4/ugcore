/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
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

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__JACOBI__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__JACOBI__

#include "lib_algebra/operator/interface/preconditioner.h"
#include "lib_algebra/small_algebra/additional_math.h"
#include "lib_algebra/cpu_algebra/vector.h"

#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallelization.h"
#endif

namespace ug{

/////////////////////////////////////////////////////////////////////////////////////////////
///		Jacobi-Iteration
/**
 * Here, the Jacobi-iteration is described for solving the linear equation
 *
 * 		\f$ A * x = b.			A \in R^{nxn}, x \in R^n, b \in R^n \f$.
 *
 * 	Most of the common linear iteration-methods base on the decomposition of A into
 * 	its diagonal (D) and strict-upper(-U) and strict-lower part (-L),
 *
 * 		\f$ A = D - L - U \f$.
 *
 * 	Among others, W. Hackbusch ('Iterative Loesung grosser Gleichungssysteme'),
 * 	distinguishes three different forms for describing a linear iteration scheme.
 * 	The general 'first normal-form' of a linear iteration scheme takes the form
 *
 * 		\f$ x^{m+1} = M * x^m + N * b \f$,
 *
 * 	with some Matrices \f$ M \f$ and \f$ N \in R^{nxn} \f$. m denotes the iteration index.
 * 	The general 'second normal-form' of a linear iteration scheme takes the form
 *
 * 		\f$ x^{m+1} = x^m - N * (A * x^m - b) \f$.
 *
 * 	Those linear iteration schemes, which can be represented by the second normal-form
 * 	are the linear, consistent iteration schemes.
 *
 * 	Introducing the correction \f$ c{m+1} := x^{m+1} - x^m \f$ and the defect
 * 	\f$ d^m := b - A * x^m \f$ the second normal-form can be rewritten as
 *
 * 		\f$ c = N * d \f$.
 *
 *	The matrix of the second normal-form for the Jacobi-method takes the simple form
 *
 *		\f$ N = D^{-1} \f$. 			.
 *
 *	References:
 * <ul>
 * <li> W. Hackbusch. Iterative Loesung grosser Gleichungssysteme
 * </ul>
 */


///	Jacobi Preconditioner
template <typename TAlgebra>
class Jacobi : public IPreconditioner<TAlgebra>
{
	public:
	///	Algebra type
		using algebra_type = TAlgebra;

	///	Vector type
		using vector_type = typename TAlgebra::vector_type;

	///	Matrix type
		using matrix_type = typename TAlgebra::matrix_type;

	///	Matrix Operator type
		using matrix_operator_type = typename IPreconditioner<TAlgebra>::matrix_operator_type;

	///	Base type
		using base_type = IPreconditioner<TAlgebra>;

	protected:
		using base_type::set_debug;
		using base_type::debug_writer;
		using base_type::write_debug;
		using base_type::damping;
		using base_type::approx_operator;

	public:
	///	default constructor
		Jacobi() {this->set_damp(1.0); m_bBlock = true;};

	///	constructor setting the damping parameter
		Jacobi(number damp) {this->set_damp(damp); m_bBlock = true;};

	/// clone constructor
		Jacobi( const Jacobi<TAlgebra> &parent )
			: base_type(parent)
		{
			set_block(parent.m_bBlock);
		}

	///	Clone
		virtual SmartPtr<ILinearIterator<vector_type> > clone()
		{
			return make_sp(new Jacobi<algebra_type>(*this));
		}


	///	returns if parallel solving is supported
		virtual bool supports_parallel() const {return true;}


	///	Destructor
		virtual ~Jacobi()
		{};

	/// sets if blocked jacobi is used (inverting block-diagonal), or plain (scalar) diagonal if false
		void set_block(bool b)
		{
			m_bBlock = b;
		}

	protected:
	///	Name of preconditioner
		virtual const char* name() const {return "Jacobi";}

	///	Preprocess routine
		virtual bool preprocess(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp)
		{

			PROFILE_BEGIN_GROUP(Jacobi_preprocess, "algebra Jacobi");

			matrix_type &mat = *pOp;
		// 	create help vector to apply diagonal
			size_t size = mat.num_rows();
			if(size != mat.num_cols())
			{
				UG_LOG("Square Matrix needed for Jacobi Iteration.\n");
				return false;
			}

			//	resize
			m_diagInv.resize(size);
#ifdef UG_PARALLEL
					//	temporary vector for the diagonal
			ParallelVector<Vector< typename matrix_type::value_type > > diag;
			diag.resize(size);

		//	copy the layouts+communicator into the vector
			diag.set_layouts(mat.layouts());

		// 	copy diagonal
			for(size_t i = 0; i < diag.size(); ++i){
				diag[i] = mat(i, i);
			}

		//	make diagonal consistent
			diag.set_storage_type(PST_ADDITIVE);
			diag.change_storage_type(PST_CONSISTENT);


			if(diag.size() > 0)
				if(CheckVectorInvertible(diag) == false)
					return false;
//			UG_ASSERT(CheckVectorInvertible(diag), "Jacobi: A has noninvertible diagonal");

#endif

//	get damping in constant case to damp at once
			number damp = 1.0;
			if(damping()->constant_damping())
				damp = damping()->damping();

			typename matrix_type::value_type m;
		// 	invert diagonal and multiply by damping
			for(size_t i = 0; i < mat.num_rows(); ++i)
			{
#ifdef UG_PARALLEL
				typename matrix_type::value_type &d = diag[i];
#else
				typename matrix_type::value_type &d = mat(i,i);
#endif
				if(!m_bBlock)
					GetDiag(m, d);
				else
					m = d;
				m *= 1./damp;
				GetInverse(m_diagInv[i], m);
			}

		//	done
			return true;
		}

		virtual bool step(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp, vector_type& c, const vector_type& d)
		{
			PROFILE_BEGIN_GROUP(Jacobi_step, "algebra Jacobi");

		// 	multiply defect with diagonal, c = damp * D^{-1} * d
		//	note, that the damping is already included in the inverse diagonal
			for(size_t i = 0; i < m_diagInv.size(); ++i)
			{
			// 	c[i] = m_diagInv[i] * d[i];
				MatMult(c[i], 1.0, m_diagInv[i], d[i]);
			}

#ifdef UG_PARALLEL

		// 	the computed correction is additive
			c.set_storage_type(PST_ADDITIVE);

		//	we make it consistent
			if(!c.change_storage_type(PST_CONSISTENT))
			{
				UG_LOG("ERROR in 'JacobiPreconditioner::apply': "
						"Cannot change parallel status of correction to consistent.\n");
				return false;
			}
#endif
		//	done
			return true;
		}

	///	Postprocess routine
		virtual bool postprocess() {return true;}


	//	overwrite function in order to specially treat constant damping
		virtual bool apply(vector_type& c, const vector_type& d)
		{
			PROFILE_BEGIN_GROUP(Jacobi_apply, "algebra Jacobi");
		//	Check that operator is initialized
			if(!this->m_bInit)
			{
				UG_LOG("ERROR in '"<<name()<<"::apply': Iterator not initialized.\n");
				return false;
			}

		//	Check parallel status
			#ifdef UG_PARALLEL
			if(!d.has_storage_type(PST_ADDITIVE))
				UG_THROW(name() << "::apply: Wrong parallel "
				               "storage format. Defect must be additive.");
			#endif

		//	Check sizes
			THROW_IF_NOT_EQUAL_4(c.size(), d.size(), approx_operator()->num_rows(), approx_operator()->num_cols());

		// 	apply iterator: c = B*d
			if(!step(approx_operator(), c, d))
			{
				UG_LOG("ERROR in '"<<name()<<"::apply': Step Routine failed.\n");
				return false;
			}

		//	apply scaling
			if(!damping()->constant_damping()){
				const number kappa = damping()->damping(c, d, approx_operator());
				if(kappa != 1.0){
					c *= kappa;
				}
			}

		//	Correction is always consistent
			#ifdef 	UG_PARALLEL
			if(!c.change_storage_type(PST_CONSISTENT))
				UG_THROW(name() << "::apply': Cannot change "
						"parallel storage type of correction to consistent.");
			#endif

		//	we're done
			return true;
		}

	protected:
	///	type of block-inverse
		using inverse_type = typename block_traits<typename matrix_type::value_type>::inverse_type;

	///	storage of the inverse diagonal in parallel
		std::vector<inverse_type> m_diagInv;
		bool m_bBlock;


};

} // end namespace ug

//#include "gpujacobi.h"

#endif
