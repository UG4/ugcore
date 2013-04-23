/*
 * linear_solver.h
 *
 *  Created on: 22.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__LINEAR_SOLVER__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__LINEAR_SOLVER__
#include <iostream>
#include <string>

#include "lib_algebra/operator/interface/preconditioned_linear_operator_inverse.h"
#include "lib_algebra/operator/interface/linear_solver_profiling.h"
#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallelization.h"
#endif

namespace ug{

/// linear solver using abstract preconditioner interface
/**
 * This class is a linear iterating scheme, that uses any implementation
 * of the ILinearIterator interface to precondition the iteration.
 *
 * \tparam 		TAlgebra		algebra type
 */
template <typename TVector>
class LinearSolver
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
	///	returns the name of the solver
		virtual const char* name() const {return "Iterative Linear Solver";}

	///	solves the system and returns the last defect
		virtual bool apply_return_defect(vector_type& x, vector_type& b)
		{
			LS_PROFILE_BEGIN(LS_ApplyReturnDefect);

			#ifdef UG_PARALLEL
			if(!b.has_storage_type(PST_ADDITIVE) || !x.has_storage_type(PST_CONSISTENT))
				UG_THROW("LinearSolver::apply: Inadequate storage format of Vectors.");
			#endif

		// 	rename b as d (for convenience)
			vector_type& d = b;

		// 	build defect:  d := b - J(u)*x
			LS_PROFILE_BEGIN(LS_BuildDefect);
			linear_operator()->apply_sub(d, x);
			LS_PROFILE_END(); //LS_BuildDefect

		// 	create correction
			LS_PROFILE_BEGIN(LS_CreateCorrection);
			SmartPtr<vector_type> spC = x.clone_without_values();
			vector_type& c = *spC;
			LS_PROFILE_END();

			LS_PROFILE_BEGIN(LS_ComputeStartDefect);
			prepare_conv_check();
			convergence_check()->start(d);
			LS_PROFILE_END();

			int loopCnt = 0;
			char ext[20]; sprintf(ext, "_iter%03d", loopCnt);
			std::string name("LS_Defect_"); name.append(ext).append(".vec");
			write_debug(d, name.c_str());
			name = std::string("LS_Solution_"); name.append(ext).append(".vec");
			write_debug(x, name.c_str());

			// 	Iteration loop
			while(!convergence_check()->iteration_ended())
			{
				char ext[20]; sprintf(ext, "_iter%03d", ++loopCnt);

			// 	Compute a correction c := B*d using one iterative step
			// 	Internally the defect is updated d := d - A*c = d - A*(x+c)
				if(preconditioner().valid()) {
					LS_PROFILE_BEGIN(LS_ApplyPrecond);
					if(!preconditioner()->apply_update_defect(c, d))
					{
						UG_LOG("ERROR in 'LinearSolver::apply': Iterator "
								"Operator applied incorrectly. Aborting.\n");
						return false;
					}
					LS_PROFILE_END(); //LS_ApplyPrecond
				}

			// 	add correction to solution: x += c
				LS_PROFILE_BEGIN(LS_AddCorrection);
				x += c;
				LS_PROFILE_END(); //LS_AddCorrection

				name = std::string("LS_Defect_"); name.append(ext).append(".vec");
				write_debug(d, name.c_str());
				name = std::string("LS_Correction_"); name.append(ext).append(".vec");
				write_debug(c, name.c_str());
				name = std::string("LS_Solution_"); name.append(ext).append(".vec");
				write_debug(x, name.c_str());

			// 	compute new defect (in parallel) d := d - A*c
				LS_PROFILE_BEGIN(LS_ComputeNewDefect);
				convergence_check()->update(d);
				LS_PROFILE_END(); //LS_ComputeNewDefect
			}

		//	write some information when ending the iteration
			if(!convergence_check()->post())
			{
				UG_LOG("ERROR in 'LinearSolver::apply': post-convergence-check "
						"signaled failure. Aborting.\n");
				return false;
			}

		//	end profiling of whole function
			LS_PROFILE_END(); //LS_ApplyReturnDefect

		//	we're done
			return true;
		}

	protected:
	///	prepares the convergence check output
		void prepare_conv_check()
		{
			convergence_check()->set_name(name());
			convergence_check()->set_symbol('%');
			if(preconditioner().valid())
            {
                std::string s;
                if(preconditioner().valid())
                    s = std::string(" (Precond: ") + preconditioner()->name() + ")";
                else
                    s = " (No Preconditioner) ";
                convergence_check()->set_info(s);
            }
		}
};

} // end namespace ug

#endif /* __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__LINEAR_SOLVER__ */
