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

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__LINEAR_SOLVER__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__LINEAR_SOLVER__
#include <iostream>
#include <string>

#include "lib_algebra/operator/interface/preconditioned_linear_operator_inverse.h"
#include "lib_algebra/operator/interface/linear_solver_profiling.h"
#include "lib_algebra/operator/interface/pprocess.h"
#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallelization.h"
#endif

namespace ug{

/// linear solver using abstract preconditioner interface
/**
 * This class is a linear iterating scheme, that uses any implementation
 * of the ILinearIterator interface to precondition the iteration.
 *
 * \tparam 		TVector		vector type
 */
template <typename TVector>
class LinearSolver
	: public IPreconditionedLinearOperatorInverse<TVector>
{
	public:
	///	Vector type
		using vector_type = TVector;

	///	Base type
		using base_type = IPreconditionedLinearOperatorInverse<vector_type>;

	///	constructors
		LinearSolver() : base_type() {}

		LinearSolver(SmartPtr<ILinearIterator<vector_type> > spPrecond)
			: base_type ( spPrecond )  {}

		LinearSolver(SmartPtr<ILinearIterator<vector_type> > spPrecond, SmartPtr<IConvergenceCheck<vector_type> > spConvCheck)
			: base_type ( spPrecond, spConvCheck)  {}

		~LinearSolver() override = default;

	protected:
		using base_type::convergence_check;
		using base_type::linear_operator;
		using base_type::preconditioner;
		using base_type::write_debug;

	public:
	///	returns the name of the solver
		[[nodiscard]] const char* name() const override {return "Iterative Linear Solver";}

	///	returns if parallel solving is supported
		[[nodiscard]] bool supports_parallel() const override {
			if(preconditioner().valid())
				return preconditioner()->supports_parallel();
			else return true;
		}

		/**
		 * Compute a correction c := B*d using one iterative step
		 * Internally the defect is updated d := d - A*c = b - A*(x+c)
		 * @param c
		 * @param d
		 */
		bool compute_correction(vector_type &c, vector_type &d)
		{
			if(preconditioner().valid()) {
				LS_PROFILE_BEGIN(LS_ApplyPrecond);

				if(!preconditioner()->apply_update_defect(c, d))
				{
					UG_LOG("ERROR in 'LinearSolver': Could not apply preconditioner. Aborting.\n");
					return false;
				}
				LS_PROFILE_END(LS_ApplyPrecond);
			}
			return true;
		}

	///	solves the system and returns the last defect
		bool apply_return_defect(vector_type& x, vector_type& b) override {

			LS_PROFILE_BEGIN(LS_ApplyReturnDefect);

			#ifdef UG_PARALLEL
			if(!b.has_storage_type(PST_ADDITIVE) || !x.has_storage_type(PST_CONSISTENT))
				UG_THROW("LinearSolver::apply: Inadequate parallel storage format of Vectors: "
							<< b.get_storage_type() << " for b (expected " << PST_ADDITIVE << "), "
							<< x.get_storage_type() << " for x (expected " << PST_CONSISTENT << ")");
			#endif

		//	debug output
			if(this->vector_debug_writer_valid())
				write_debug(b, std::string("LS_RHS") + ".vec");
			
		// 	rename b as d (for convenience)
			vector_type& d = b;

		// 	build defect:  d := b - J*x
			LS_PROFILE_BEGIN(LS_BuildDefect);
			linear_operator()->apply_sub(d, x);
			LS_PROFILE_END(LS_BuildDefect);

		// 	create correction
			LS_PROFILE_BEGIN(LS_CreateCorrection);
			SmartPtr<vector_type> spC = x.clone_without_values();
			vector_type& c = *spC;
			#ifdef UG_PARALLEL
				// this is ok if clone_without_values() inits with zeros
				c.set_storage_type(PST_CONSISTENT);
			#endif
			LS_PROFILE_END(LS_CreateCorrection);

			LS_PROFILE_BEGIN(LS_ComputeStartDefect);
			prepare_conv_check();
			convergence_check()->start(d);
			LS_PROFILE_END(LS_ComputeStartDefect);

			int loopCnt = 0;
			write_debugXCD(x, c, d, loopCnt, false);

		// 	Iteration loop
			while(!convergence_check()->iteration_ended())
			{
				enter_precond_debug_section(loopCnt);
				if( !compute_correction(c, d) )
				{
					this->leave_vector_debug_writer_section();
					return false;
				}
				this->leave_vector_debug_writer_section();

			//	post-process the correction
				m_corr_post_process.apply (c);

			// 	add correction to solution: x += c
				LS_PROFILE_BEGIN(LS_AddCorrection);
				x += c;
				LS_PROFILE_END(LS_AddCorrection);
				
				write_debugXCD(x, c, d, ++loopCnt, true);

			// 	compute norm of new defect (in parallel)
				LS_PROFILE_BEGIN(LS_ComputeNewDefectNorm);
				convergence_check()->update(d);
				LS_PROFILE_END(LS_ComputeNewDefectNorm);
			}

		//	write some information when ending the iteration
			if(!convergence_check()->post())
			{
				UG_LOG("ERROR in 'LinearSolver::apply': post-convergence-check "
						"signaled failure. Aborting.\n");
				return false;
			}

		//	end profiling of whole function
			LS_PROFILE_END(LS_ApplyReturnDefect);

		//	we're done
			return true;
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
	
	/// debugger output: solution, correction, defect
		void write_debugXCD(vector_type &x, vector_type &c, vector_type &d, int loopCnt, bool bWriteC)
		{
			if(!this->vector_debug_writer_valid()) return;
			char ext[20]; snprintf(ext, sizeof(ext),"_iter%03d", loopCnt);
			write_debug(d, std::string("LS_Defect") + ext + ".vec");
			if(bWriteC) write_debug(c, std::string("LS_Correction") + ext + ".vec");
			write_debug(x, std::string("LS_Solution") + ext + ".vec");
		}
		
	/// debugger section for the preconditioner
		void enter_precond_debug_section(int loopCnt)
		{
			if(!this->vector_debug_writer_valid()) return;
			char ext[20]; snprintf(ext, sizeof(ext),"_iter%03d", loopCnt);
			this->enter_vector_debug_writer_section(std::string("LS_Precond_") + ext);
		}

	protected:
	///	postprocessor for the correction in the iterations
		PProcessChain<vector_type> m_corr_post_process;
};

} // end namespace ug

#endif