/*
 * cg.h
 *
 *  Created on: 22.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_SOLVER__CG__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_SOLVER__CG__

#include <iostream>
#include <sstream>

#include "lib_algebra/operator/interface/operator.h"
#include "common/profiler/profiler.h"
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
	///	default constructor
		CG() {}

	///	constructor setting preconditioner and convergence check
		CG( SmartPtr<ILinearIterator<vector_type> > spPrecond,
		    SmartPtr<IConvergenceCheck<vector_type> > spConvCheck) :
		    base_type(spPrecond, spConvCheck)
		{};

	///	name of solver
		virtual const char* name() const {return "CG";}

	///	Solve J(u)*x = b, such that x = J(u)^{-1} b
		virtual bool apply_return_defect(vector_type& x, vector_type& b)
		{
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
		//	todo: 	It would be sufficient to copy only the pattern and
		//			without initializing, but in parallel we have to copy communicators
			vector_type q; q.create(r.size()); q = r;
			vector_type z; z.create(x.size()); z = x;
			vector_type p; p.create(x.size()); p = x;

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

	protected:
	///	adjust output of convergence check
		void prepare_conv_check()
		{
		//	set iteration symbol and name
			convergence_check()->set_name(name());
			convergence_check()->set_symbol('%');

		//	set preconditioner string
			std::stringstream ss;
			if(preconditioner().valid())
				ss <<" (Precond: "<<preconditioner()->name()<<")";
			else
				ss << " (No Preconditioner) ";
			convergence_check()->set_info(ss.str());
		}

	protected:
		number VecProd(vector_type& a, vector_type& b)
		{
			return a.dotprod(b);
		}
};

} // end namespace ug

#endif /* __H__UG__LIB_DISC__OPERATOR__LINEAR_SOLVER__CG__ */
