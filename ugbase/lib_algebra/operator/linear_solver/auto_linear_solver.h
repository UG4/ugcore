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
#include "common/util/ostream_util.h"
namespace ug{

// this is for matrix operators which recalculate their "real" operator
// if necessary
class UpdateableMatrixOperator
{
public:
	virtual ~UpdateableMatrixOperator() {}
	virtual void calculate_matrix() = 0;
};

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
	AutoLinearSolver //(double reductionAlwaysAccept, double worseThenAverage, double desiredDefect)
	(double desiredDefect, double desiredReduction)
	{
		m_bInited=-1;
//		m_N = 0;
//		m_reductionAlwaysAccept = reductionAlwaysAccept;
//		m_worseThenAverage = worseThenAverage;
		m_desiredDefect = desiredDefect;
		m_desiredReduction = desiredReduction;
		m_initCalled = m_initsDone = 0;
		m_savedTime = 0.0;
	}
	AutoLinearSolver()
	{
		m_bInited=-1;
//		m_N = 0;
		m_reductionAlwaysAccept = 0.001;
		m_worseThenAverage = 2.0;
		m_initCalled = m_initsDone = 0;
		m_savedTime = 0.0;
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
//	size_t m_N;
	SmartPtr<ILinearOperator<vector_type,vector_type> > pJ;
	vector_type m_u;
	double m_avgReduction;
	size_t m_initsDone;
	size_t m_initCalled;

	double m_lastInitTime;
	double m_reductionPerTime, m_lastCallTime, m_lastCallReduction;
	double m_desiredReduction, m_desiredDefect;
	double m_savedTime;

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
			m_u = u;
			return init_op(J);
		}
		virtual bool init(SmartPtr<ILinearOperator<vector_type,vector_type> > J)
		{
			m_u.resize(0);
			return init_op(J);
		}

		bool init_op(SmartPtr<ILinearOperator<vector_type,vector_type> > J)
		{
			ILinearOperatorInverse<vector_type, vector_type>::init(J);
			//UG_LOG("ALS:init\n");
			m_initCalled++;
			if(m_bInited == -1)
			{
				m_bInited=false;
				pJ = J;
				reinit();
			}
			else
			{
				m_bInited=false;
				pJ = J;
//				if(u.size() != m_N)
//					reinit(u);
//				else
//					m_u = u;
			}

			return true;
		}

		bool reinit()
		{
			if(m_bInited) return true;
			m_bInited = true;
			//UG_LOG("ALS:reinit\n");

			double tStart = get_clock_s();
			// this does not work with stuff.
			SmartPtr<UpdateableMatrixOperator> uo = pJ.template cast_dynamic<UpdateableMatrixOperator> ();
			if(uo.valid()) uo->calculate_matrix();
			if(m_u.size() != 0.0)
			{ 	if(!base_type::init(pJ, m_u)) return false; }
			else if(!base_type::init(pJ)) return false;
			m_lastInitTime = get_clock_s()-tStart;

			m_initsDone++;

//			m_N = u.size();

			return true;
		}

		bool compute_correction(vector_type &c, vector_type &d)
		{
						// 	Compute a correction c := B*d using one iterative step
			// 	Internally the defect is updated d := d - A*c = b - A*(x+c)
			if(preconditioner().valid()) {
				LS_PROFILE_BEGIN(LS_ApplyPrecond);

				if(!preconditioner()->apply(c, d))
				{
					UG_LOG("ERROR in 'LinearSolver::apply': Iterator "
							"Operator applied incorrectly. Aborting.\n");
					return false;
				}
				linear_operator()->apply_sub(d, c);
				LS_PROFILE_END(LS_ApplyPrecond);
			}

			return true;
		}

		virtual bool apply(vector_type& x, const vector_type& b)
		{
			//UG_LOG("ALS:apply\n");
			SmartPtr<vector_type> spB = b.clone_without_values();
			*spB = b;
			return apply_return_defect(x, *spB);
		}
	///	solves the system and returns the last defect
		virtual bool apply_return_defect(vector_type& x, vector_type& b)
		{
			//UG_LOG("ALS:return_defect\n");
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
			LS_PROFILE_END(LS_BuildDefect);

		// 	create correction
			LS_PROFILE_BEGIN(LS_CreateCorrection);
			SmartPtr<vector_type> spC = x.clone_without_values();
			vector_type& c = *spC;
			LS_PROFILE_END(LS_CreateCorrection);

			LS_PROFILE_BEGIN(LS_ComputeStartDefect);
			prepare_conv_check();
			convergence_check()->start(d);
			LS_PROFILE_END(LS_ComputeStartDefect);

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
					double tStartIterationTime = get_clock_s();
					while(!convergence_check()->iteration_ended())
					{
						//double tComputeTime = get_clock_s();
						if( !compute_correction(c, d) ) return false;
						x += c;
						convergence_check()->update(d);


						double r = convergence_check()->rate();
						double reduction = convergence_check()->reduction();
						double defect = convergence_check()->defect();
						double steps = convergence_check()->step();
//						if(r > m_reductionAlwaysAccept
//							&& (r >= 1 || r * m_worseThenAverage > m_avgReduction))
						double spentTime = get_clock_s() - tStartIterationTime;

						double tComputeTime = spentTime/steps; //get_clock_s()-tComputeTime;

//						double timeForOriginalAlgorithm = reduction/m_reductionPerTime + m_lastInitTime;
						double approxSteps =
								std::min(log(m_desiredDefect/defect)/log(r),
										log(m_desiredReduction/reduction)/log(r) );
						double approxRemainingTimeForSolution = approxSteps*tComputeTime;
						/*UG_LOG("\n");
						UG_LOG("m_lastCallReduction = " << reset_floats << m_lastCallReduction << "\n");
						UG_LOG("reduction = " << reset_floats << reduction << "\n");
						UG_LOG("reduction remain = " << reset_floats << m_lastCallReduction/reduction << "\n");
						UG_LOG("approxSteps = " << reset_floats << approxSteps << "\n");
						UG_LOG("approxTimeForSolution = " << approxTimeForSolution << "\n");
						UG_LOG("spentTime = " << spentTime << "\n");
						UG_LOG("tComputeTime = " << tComputeTime << "\n");
						UG_LOG("rate = " << r << "\n");
						UG_LOG("reduction = " << reduction << "\n");*/
						if(r > 1)
						{
							m_savedTime -= spentTime;
							UG_LOG("AutoLinearSolver: REINIT because reduction rate >= 1.")
							UG_LOG(" [ Inits called: " << m_initCalled << ", inits done: " << reset_floats << m_initsDone+1
									<< " (" << (100.0*(m_initsDone+1))/m_initCalled << " %), Total Time saved: " << m_savedTime << " s]\n");


							reinit();
							break;
						}
						double approxResolveTime = m_lastCallTime + m_lastInitTime;
						if(approxResolveTime < approxRemainingTimeForSolution)
						{
							m_savedTime -= spentTime;
							UG_LOG("AutoLinearSolver: REINIT because\n" << reset_floats )
							UG_LOG("  approximated remaining Time for Solution with old preconditioner = " << approxRemainingTimeForSolution << " s\n");
							UG_LOG("  > approximated Time for Reinit preconditioner and solve          = " << approxResolveTime << " s\n");
							UG_LOG(" with old preconditioner:\n");
							UG_LOG("  reduction rate = " << r << " avg reduction = " << convergence_check()->avg_rate() << "\n")
							UG_LOG("  steps = " << steps << ", reduction = " << reduction << "\n")
							UG_LOG("  spentTime = " << spentTime << ", time per Step = " << tComputeTime << "\n");
							UG_LOG(" approximation of solution with old:\n")
							UG_LOG("  reduction remain = " << reset_floats << m_lastCallReduction/reduction << "\n");
							UG_LOG("  approxSteps with last reduction rate = " << approxSteps << "\n");
							UG_LOG(" [ Inits called: " << m_initCalled << ", inits done: " << reset_floats << m_initsDone+1
									<< " (" << (100.0*(m_initsDone+1))/m_initCalled << " %), Total Time saved: " << m_savedTime << " s]\n");
							reinit();

							break;
						}
						else
						{
							x2 = x;
							d2 = d;
						}
					}
					if(!m_bInited)
					{
						double spentTime = get_clock_s() - tStartIterationTime;
						m_savedTime += m_lastCallTime+m_lastInitTime-spentTime;
						UG_LOG("AutoLinearSolver solved with old preconditioner [")
						UG_LOG(" Inits called: " << m_initCalled << ", inits done: " << reset_floats << m_initsDone+1
								<< " (" << (100.0*(m_initsDone+1))/m_initCalled << " %), Total Time saved: " << m_savedTime << " s]\n");
						UG_LOG(" time spent with old    = " << reset_floats << spentTime << " s, ");
						UG_LOG(" last time init + solve = " << m_lastCallTime+m_lastInitTime << " s, ");

						if(m_lastCallTime+m_lastInitTime > spentTime)
						{	UG_LOG("saved " << m_lastCallTime+m_lastInitTime - spentTime << " s!.\n"); }
						else
						{	UG_LOG("spent too much time with old! " <<  spentTime  - m_lastCallTime+m_lastInitTime<< " s!\n"); }

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
				double T = get_clock_s();
				convergence_check()->start(d);
				while(!convergence_check()->iteration_ended())
				{
					if( !compute_correction(c, d) ) return false;
					x += c;
					convergence_check()->update(d);
				}
				m_avgReduction = convergence_check()->avg_rate();
				T = get_clock_s() - T;
				m_reductionPerTime = convergence_check()->reduction()/T;
				m_lastCallTime = T;
				m_lastCallReduction = convergence_check()->reduction();
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

		void print_information()
		{
			UG_LOG("AutoLinearSolver:\n");
			UG_LOG(" avg reduction is " << m_avgReduction << "\n");
			UG_LOG(" Inits called: " << m_initCalled << ", inits done: " << reset_floats << m_initsDone
					<< " (" << (100.0*(m_initsDone))/m_initCalled << " %)\n");
			UG_LOG(" m_reductionPerTime              = " << m_reductionPerTime << " s\n");
			UG_LOG(" reductionTime for 0.1 reduction = " << 0.1/m_reductionPerTime << " s\n");
			UG_LOG(" m_lastInitTime                  = " << m_lastInitTime << " s\n");
			UG_LOG(" SAVED TIME: " << m_savedTime << " s\n");
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
