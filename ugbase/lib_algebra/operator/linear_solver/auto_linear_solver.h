/*
 * linear_solver.h
 *
 *  Created on: 30.10.2013
 *      Author: martin rupp
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__AUTO_LINEAR_OPERATOR__LINEAR_SOLVER__
#define __H__UG__LIB_DISC__OPERATOR__AUTO_LINEAR_OPERATOR__LINEAR_SOLVER__
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
class AutoLinearSolver
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
	AutoLinearSolver(double reductionAlwaysAccept, double worseThenAverage)
	{
		m_bInited=-1;
		m_N = 0;
		m_reductionAlwaysAccept = reductionAlwaysAccept;
		m_worseThenAverage = worseThenAverage;
	}
	AutoLinearSolver()
	{
		m_bInited=-1;
		m_N = 0;
		m_reductionAlwaysAccept = 0.001;
		m_worseThenAverage = 2.0;
	}

///	returns if parallel solving is supported
	virtual bool supports_parallel() const
	{
		if(preconditioner().valid())
			return preconditioner()->supports_parallel();
		else return true;
	}

private:
	double m_reductionAlwaysAccept;
	double m_worseThenAverage;
	int m_bInited;
	size_t m_N;
	SmartPtr<ILinearOperator<vector_type,vector_type> > pJ;
	vector_type m_u;
	double m_avgReduction;

	public:

	void set_reduction_always_accept(double d)
	{
		m_reductionAlwaysAccept = d;
	}

	void set_reinit_when_worse_then_average(double d)
	{
		m_worseThenAverage = d;
	}

	///	returns the name of the solver
		virtual const char* name() const {return "Auto Iterative Linear Solver";}


		virtual bool init(SmartPtr<ILinearOperator<vector_type,vector_type> > J, const vector_type& u)
		{
			if(m_bInited == -1)
			{
				m_bInited=false;
				pJ = J;
				reinit(u);
			}
			else
			{
				m_bInited=false;
				pJ = J;
				if(u.size() != m_N)
					reinit(u);
				else
					m_u = u;
			}

			return true;
		}

		bool reinit(const vector_type &u)
		{
			if(m_bInited) return true;
			m_bInited = true;

			if(!base_type::init(pJ, u)) return false;

			m_N = u.size();

			return true;
//			LS_PROFILE_END();
		}

		bool compute_correction(vector_type &c, vector_type &d)
		{
						// 	Compute a correction c := B*d using one iterative step
			// 	Internally the defect is updated d := d - A*c = b - A*(x+c)
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

			return true;
		}

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


			if(m_bInited == false)
			{
				vector_type x2 = x;
				vector_type d2 = d;
				try{
					while(!convergence_check()->iteration_ended())
					{
						if( !compute_correction(c, d) ) return false;
						x += c;
						convergence_check()->update(d);
						double r = convergence_check()->rate();
						if(r > m_reductionAlwaysAccept
							&& (r >= 1 || r * m_worseThenAverage > m_avgReduction))
						{
							UG_LOG("REINIT because of " << r << " and avg is " << m_avgReduction << " (worseThenAverage = " << m_worseThenAverage << "\n");
							reinit(m_u);
							break;
						}
						else
						{
							x2 = x;
							d2 = d;
						}
					}
				}
				catch(...)
				{

				}
				x = x2;
				d = d2;
			}

			if(!convergence_check()->iteration_ended())
			{
				convergence_check()->start(d);
				while(!convergence_check()->iteration_ended())
				{
					if( !compute_correction(c, d) ) return false;
					x += c;
					convergence_check()->update(d);
				}
				m_avgReduction = convergence_check()->avg_rate();
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
