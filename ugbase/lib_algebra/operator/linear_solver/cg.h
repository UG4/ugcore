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

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_SOLVER__CG__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_SOLVER__CG__

#include <iostream>
#include <string>

#include "lib_algebra/operator/interface/operator.h"
#include "common/profiler/profiler.h"
#include "lib_algebra/operator/interface/pprocess.h"
#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallelization.h"
#endif

namespace ug{

///	the CG method as a solver for linear operators
/**
 * This class implements the CG - method for the solution of linear operator
 * problems like A*x = b, where the solution x = A^{-1} b is computed.
 *
 * For detailed description of the algorithm, please refer to:
 *
 * - Barrett, Berry, Chan, Demmel, Donatom Dongarra, Eijkhout, Pozo, Romine,
 * 	 Van der Vorst, "Templates for the Solution of Linear Systems: Building
 * 	 Blocks for Iterative Methods", p.13, Fig, 2.5
 *
 * - Saad, "Iterative Methods For Sparse Linear Systems", p277, Alg. 9.1
 *
 * \tparam 	TVector		vector type
 */
template <typename TVector>
class CG
	: public IPreconditionedLinearOperatorInverse<TVector>
{
	public:
	///	Vector type
		typedef TVector vector_type;

	///	Base type
		typedef IPreconditionedLinearOperatorInverse<vector_type> base_type;

	protected:
		using base_type::convergence_check;
		using base_type::linear_operator;
		using base_type::preconditioner;
		using base_type::write_debug;

	public:
	///	constructors
		CG() : base_type() {}

		CG(SmartPtr<ILinearIterator<vector_type,vector_type> > spPrecond)
			: base_type ( spPrecond )  {}

		CG(SmartPtr<ILinearIterator<vector_type,vector_type> > spPrecond, SmartPtr<IConvergenceCheck<vector_type> > spConvCheck)
			: base_type ( spPrecond, spConvCheck)  {}

	///	name of solver
		virtual const char* name() const {return "CG";}

	///	returns if parallel solving is supported
		virtual bool supports_parallel() const
		{
			if(preconditioner().valid())
				return preconditioner()->supports_parallel();
			return true;
		}

	///	Solve J(u)*x = b, such that x = J(u)^{-1} b
		virtual bool apply_return_defect(vector_type& x, vector_type& b)
		{
			PROFILE_BEGIN_GROUP(CG_apply_return_defect, "CG algebra");
		//	check parallel storage types
			#ifdef UG_PARALLEL
			if(!b.has_storage_type(PST_ADDITIVE) || !x.has_storage_type(PST_CONSISTENT))
				UG_THROW("CG::apply_return_defect:"
								"Inadequate storage format of Vectors.");
			#endif

		// 	rename r as b (for convenience)
			vector_type& r = b;

		// 	Build defect:  r := b - J(u)*x
			linear_operator()->apply_sub(r, x);

		// 	create help vector (h will be consistent r)
			SmartPtr<vector_type> spQ = r.clone_without_values(); vector_type& q = *spQ;
			SmartPtr<vector_type> spZ = x.clone_without_values(); vector_type& z = *spZ;
			SmartPtr<vector_type> spP = x.clone_without_values(); vector_type& p = *spP;

		// 	Preconditioning
			if(preconditioner().valid())
			{
				// apply z = M^-1 * s
				if(!preconditioner()->apply(z, r))
				{
					UG_LOG("ERROR in 'CG::apply_return_defect': "
							"Cannot apply preconditioner. Aborting.\n");
					return false;
				}
			}
			else z = r;
			
		// 	make z consistent
			#ifdef UG_PARALLEL
			if(!z.change_storage_type(PST_CONSISTENT))
				UG_THROW("CG::apply_return_defect: "
								"Cannot convert z to consistent vector.");
			#endif

		//	post-process the correction
			m_corr_post_process.apply (z);

		//	compute start defect
			prepare_conv_check();
			convergence_check()->start(r);

		// 	start search direction
			p = z;

		// 	start rho
			number rhoOld = VecProd(z, r), rho;

		// 	Iteration loop
			while(!convergence_check()->iteration_ended())
			{
			// 	Build q = A*p (q is additive afterwards)
				linear_operator()->apply(q, p);

			// 	lambda = (q,p)
				const number lambda = VecProd(q, p);

			//	check lambda
				if(lambda == 0.0)
				{
					UG_LOG("ERROR in 'CG::apply_return_defect': lambda=" <<
					       lambda<< " is not admitted. Aborting solver.\n");
					return false;
				}

			//	alpha = rho / (q,p)
				const number alpha = rhoOld/lambda;

			// 	Update x := x + alpha*p
				VecScaleAdd(x, 1.0, x, alpha, p);

			// 	Update r := r - alpha*t
				VecScaleAdd(r, 1.0, r, -alpha, q);

			// 	Check convergence
				convergence_check()->update(r);
				if(convergence_check()->iteration_ended()) break;

			// 	Preconditioning
				if(preconditioner().valid())
				{
				// 	apply z = M^-1 * r
					if(!preconditioner()->apply(z, r))
					{
						UG_LOG("ERROR in 'CG::apply_return_defect': "
								"Cannot apply preconditioner. Aborting.\n");
						return false;
					}
				}
				else z = r;

				#ifdef UG_PARALLEL
			// 	make z consistent
				if(!z.change_storage_type(PST_CONSISTENT))
					UG_THROW("CG::apply_return_defect': "
									"Cannot convert z to consistent vector.");
				#endif

			//	post-process the correction
				m_corr_post_process.apply (z);

			// 	new rho = (z,r)
				rho = VecProd(z, r);

			// 	new beta = rho / rhoOld
				const number beta = rho/rhoOld;

			// 	new direction p := beta * p + z
				VecScaleAdd(p, beta, p, 1.0, z);

			// 	remember old rho
				rhoOld = rho;
			}

		//	post output
			return convergence_check()->post();
		}
		
	///	adds a post-process for the iterates
		void add_postprocess_corr (SmartPtr<IPProcessVector<vector_type> > p)
		{
			m_corr_post_process.add (p);
		}

	///	removes a post-process for the iterates
		void remove_postprocess_corr (SmartPtr<IPProcessVector<vector_type> > p)
		{
			m_corr_post_process.remove (p);
		}

	protected:
	///	adjust output of convergence check
		void prepare_conv_check()
		{
		//	set iteration symbol and name
			convergence_check()->set_name(name());
			convergence_check()->set_symbol('%');

		//	set preconditioner string
			std::string s;
			if(preconditioner().valid())
			  s = std::string(" (Precond: ") + preconditioner()->name() + ")";
			else
				s = " (No Preconditioner) ";
			convergence_check()->set_info(s);
		}

	protected:
		number VecProd(vector_type& a, vector_type& b)
		{
			return a.dotprod(b);
		}
	
	protected:
	///	postprocessor for the correction in the iterations
		/**
		 * These postprocess operations are applied to the preconditioned
		 * defect before the orthogonalization. The goal is to prevent the
		 * useless kernel parts to prevail in the (floating point) arithmetics.
		 */
		PProcessChain<vector_type> m_corr_post_process;
};

} // end namespace ug

#endif /* __H__UG__LIB_DISC__OPERATOR__LINEAR_SOLVER__CG__ */
