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

#ifndef __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__LINE_SEARCH__
#define __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__LINE_SEARCH__

#include <ostream>
#include <string>
#include <vector>
#include <cmath>
#include <sstream>
#include <limits>

#include "common/common.h"

//#include "lib_disc/operator/non_linear_operator/newton_solver/nestedNewtonRFSwitch.h"

//#if ENABLE_NESTED_NEWTON_RESOLFUNC_UPDATE

#include "lib_disc/operator/non_linear_operator/newton_solver/newtonUpdaterGeneric.h"

//#endif



namespace ug{

///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
// Line Search
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////

/** LineSearch
 *
 * This is the base class for a line search object. An instance is
 * passed to a newton solver to control the line search.
 *
 */
template <typename TVector>
class ILineSearch
{
	public:
		typedef TVector vector_type;

	public:
	/// set string to be printed before each output of line search
		virtual void set_offset(std::string offset) = 0;

	/**
	 *	Performs a line search to a given direction.
	 *
	 * \param[in]		Op		Non-linear operator
	 * \param[in]  		u 		current solution
	 * \param[in]		p		search direction
	 * \param[in,out]	d		defect
	 * \param[in]		defect	norm of current defect
	 *
	 * \return 	true 		if line search successful
	 * 			false 		if line search failed
	 */
		virtual bool search(SmartPtr<IOperator<vector_type> > spOp,
		                    vector_type& u, vector_type& p,
		                    vector_type& d, number defect) = 0;

	///	returns information about configuration parameters
		/**
		 * this should return necessary information about parameters and possibly
		 * calling config_string of subcomponents.
		 *
		 * \returns std::string	necessary information about configuration parameters
		 */

		virtual std::string config_string() const = 0;

	/// virtual destructor
		virtual ~ILineSearch() {}

//#if ENABLE_NESTED_NEWTON_RESOLFUNC_UPDATE
		virtual void setNewtonUpdater( SmartPtr<NewtonUpdaterGeneric<TVector> > nU ) = 0;
//#endif
		virtual bool createNewtonUpdater() = 0;

};

/// standard implementation of the line search based on the "sufficient descent"
template <typename TVector>
class StandardLineSearch : public ILineSearch<TVector>
{
	public:
	//	type of
		typedef TVector vector_type;

	public:
	///	default constructor (setting default values)
		StandardLineSearch()
		 :	 m_maxSteps(10), m_lambdaStart(1.0), m_lambdaReduce(0.5), m_alpha(0.25),
		  	 m_maxDefect(1e+10), m_verbose(true), m_bAcceptBest(false), m_bCheckAll(false), m_offset(""),
			 m_newtonUpdater(SPNULL)
//#if ENABLE_NESTED_NEWTON_RESOLFUNC_UPDATE
//			 ,
//			m_newtonUpdater(new NewtonUpdaterGeneric<vector_type>{})
//#endif
			 {};

	///	constructor
		StandardLineSearch(int maxSteps, number lambdaStart, number lambdaReduce, bool bAcceptBest)
		 :	 m_maxSteps(maxSteps), m_lambdaStart(lambdaStart), m_lambdaReduce(lambdaReduce), m_alpha(0.25),
			 m_maxDefect(1e+10), m_verbose(true), m_bAcceptBest(bAcceptBest), m_bCheckAll(false), m_offset(""),
			 m_newtonUpdater(SPNULL)
//#if ENABLE_NESTED_NEWTON_RESOLFUNC_UPDATE
//			 ,
//			 m_newtonUpdater(new NewtonUpdaterGeneric<vector_type>{})
//#endif
			 {};

	///	constructor
		StandardLineSearch(int maxSteps, number lambdaStart, number lambdaReduce, bool bAcceptBest, bool bCheckAll)
		 :	 m_maxSteps(maxSteps), m_lambdaStart(lambdaStart), m_lambdaReduce(lambdaReduce), m_alpha(0.25),
			 m_maxDefect(1e+10), m_verbose(true), m_bAcceptBest(bAcceptBest), m_bCheckAll(bCheckAll), m_offset(""),
			 m_newtonUpdater(SPNULL)
//#if ENABLE_NESTED_NEWTON_RESOLFUNC_UPDATE
//			 ,
//			 m_newtonUpdater(new NewtonUpdaterGeneric<vector_type>{})
//#endif
			 {};

	///	returns information about configuration parameters
		virtual std::string config_string() const
		{
			std::stringstream ss;
			ss << "StandardLineSearch( maxSteps = " << m_maxSteps << ", lambdaStart = " << m_lambdaStart << ", lambdaReduce = " << m_lambdaReduce << ", accept best = " <<
					(m_bAcceptBest ? "true" : "false") << " check all = " << (m_bCheckAll ? "true" : "false");
			return ss.str();

		}

	///	sets maximum number of line search steps
		void set_maximum_steps(int steps) {m_maxSteps = steps;}

	///	sets start factor
		void set_lambda_start(number start) {m_lambdaStart = start;}

	///	sets factor by which line search factor is multiplied in each step
		void set_reduce_factor(number factor) {m_lambdaReduce = factor;}
		
	///	sets the factor controlling the sufficient descent
		void set_suff_descent_factor(number factor) {m_alpha = factor;}

	///	sets iff after max_steps the best try is used
		void set_accept_best(bool bAcceptBest) {m_bAcceptBest = bAcceptBest;}

	///	sets iff all the max_steps line search steps must be tested even if the sufficient descent is achieved
		void set_check_all(bool bCheckAll) {m_bCheckAll = bCheckAll;}

	///	sets maximum allowed norm of the defect (an exception is thrown if this value if exceeded)
		void set_maximum_defect(number maxDef) {m_maxDefect = maxDef;}

	///	sets if info should be printed
		void set_verbose(bool level) {m_verbose = level;}

	///	\copydoc ILineSearch::set_offset
		virtual void set_offset(std::string offset) {m_offset = offset;};

	///	\copydoc ILineSearch::search
		virtual bool search(SmartPtr<IOperator<vector_type> > spOp,
		                    vector_type& u, vector_type& p,
		                    vector_type& d, number defect)
		{
			PROFILE_BEGIN_GROUP(StandardLineSearch_search, ""); // group?
		// 	clone pattern for s
			s.resize(u.size());

		//	start factor
			number lambda = m_lambdaStart;

		//	some values
			number norm, norm_old = defect;
			bool converged = false;
			std::vector<number> vRho;

		// remember u
			s = u;

		//	print heading line
		if(m_verbose)
			UG_LOG(m_offset << "   ++++ Line Search:\n"
							<< "   +  Iter       lambda        Defect          Rate \n");


		//	loop line search steps
			for(int k = 1; k <= m_maxSteps; ++k)
			{
//#if ENABLE_NESTED_NEWTON_RESOLFUNC_UPDATE
				if( m_newtonUpdater != SPNULL )
				{
					// 	try on line u := u - lambda*p

					bool acceptedNewtonUpdate = m_newtonUpdater->updateSolution(u, 1.0, u, (-1)*lambda, p);

					if( ! acceptedNewtonUpdate )
					{
						UG_LOG("Update in Line Search did not work.\n");
						norm = std::numeric_limits<number>::max();
							//VecScaleAdd(u, 1.0, u, (-1)*lambda, p);
						//return false;
					}
					else
					{
						// 	compute new Defect
						spOp->prepare(u);
						spOp->apply(d, u);

						//	compute new Residuum
						norm = d.norm();
					}
				}
				else
				{
//#else
				// 	try on line u := u - lambda*p
					VecScaleAdd(u, 1.0, u, (-1)*lambda, p);

				// 	compute new Defect
					spOp->prepare(u);
					spOp->apply(d, u);

				//	compute new Residuum
					norm = d.norm();
				}
//#endif


			// 	compute reduction
				vRho.push_back(norm/norm_old);

			//	print rate
				if(m_verbose)
					UG_LOG(m_offset << "   + " << std::setw(4)
							<< k << ":   " << std::setw(11)
							<< std::resetiosflags( ::std::ios::scientific )<<
							lambda << "     "
							<< std::scientific << norm << "   " << vRho.back() <<"\n");

			// 	check if reduction fits
				if(vRho.back() <= 1 - m_alpha * std::fabs(lambda))
				{
					converged = true;
					if(!m_bCheckAll) break;
				}

				lambda *= m_lambdaReduce;

			//	check if maximum number of steps reached
				if(k == m_maxSteps)
				{
				//	if not accept best, line search failed
					if(!m_bAcceptBest)
					{
						UG_LOG(m_offset << "   ++++ Line Search did not converge.\n");
						return false;
					}

				//	search minimum
					size_t best = 0;
					number rho_min = vRho.front();
					for(size_t i = 1; i < vRho.size(); ++i)
					{
						if(rho_min > vRho[i])
						{
							rho_min = vRho[i];
							best = i;
						}
					}

				/*	check if best is converging (i.e. rho < 1)
					if(vRho[best] >= 1)
					{
						UG_LOG(m_offset << "   ++++ Accept Best: No try with "
								"Rate < 1, cannot accept any line search step.\n");
						UG_LOG(m_offset << "   ++++ Line Search did not converge.\n");
						return false;
					}*/

				//	accept best
					if(m_verbose)
						UG_LOG(m_offset << "   ++++ Accept Best: Accepting step " <<
							best+1 << ", Rate = "<< vRho[best] <<".\n");

//#if ENABLE_NESTED_NEWTON_RESOLFUNC_UPDATE
					if( m_newtonUpdater != SPNULL )
					{
						// 	try on line u := u - lambda*p

						if( ! m_newtonUpdater->updateSolution(u, 1.0, s, (-1)*m_lambdaStart*std::pow(m_lambdaReduce, (number)best), p) )
						{
							UG_LOG("Update in Line Search kmax did not work.\n");

							norm = std::numeric_limits<number>::max();
							//return false;
						}
						else
						{
							spOp->prepare(u);
							spOp->apply(d, u);

							// compute new Residuum
							norm = d.norm();
						}
					}
					else
					{
//#else
					// 	try on line u := u - lambda*p
						VecScaleAdd(u, 1.0, s, (-1)*m_lambdaStart*std::pow(m_lambdaReduce, (number)best), p);
					// 	compute new Defect
						spOp->prepare(u);
						spOp->apply(d, u);

						// compute new Residuum
						norm = d.norm();
					}
//#endif


					// check if defect-norm is smaller than maximum allowed defect value
					if (norm > m_maxDefect)
					{
						UG_LOG("ERROR in 'StandardLineSearch::search':"
								" maximum defect-limit is reached.\n");
						return false;
					}

				//	break to finish
					break;
				}

//#if ENABLE_NESTED_NEWTON_RESOLFUNC_UPDATE

				if( m_newtonUpdater != SPNULL )
				{
					// 	reset u and eventual local variables
					m_newtonUpdater->resetSolution(u,s);
				}
				else
				{
//#else
				// 	reset u
					u = s;
				}
//#endif

			}

		//	print end line
			if(m_verbose)
			{
				//only for rate < 1, we call it "Line Search converged"
				if(converged)
					UG_LOG(m_offset << "   ++++ Line Search converged.\n");
			}


//#if ENABLE_NESTED_NEWTON_RESOLFUNC_UPDATE
			if( m_newtonUpdater != SPNULL )
			{
				if( ! m_newtonUpdater->tellAndFixUpdateEvents(u) )
				{
					UG_LOG("unable to fix local Newton updates" << std::endl );
					return false;
				}
			}
//#endif
		//	we're done
			return true;
		}


//#if ENABLE_NESTED_NEWTON_RESOLFUNC_UPDATE

		virtual void setNewtonUpdater( SmartPtr<NewtonUpdaterGeneric<TVector> > nU )
		{
			m_newtonUpdater = nU;
		}

		virtual bool createNewtonUpdater()
		{
			if( m_newtonUpdater != SPNULL )
			{
				m_newtonUpdater = SmartPtr<NewtonUpdaterGeneric<TVector> >
										  (new NewtonUpdaterGeneric<TVector>{});

				return true;
			}

			return false;

		}


//#endif

	protected:
	/// solution in line direction
		vector_type s;

	protected:
	/// maximum number of steps to be performed
		int m_maxSteps;

	/// initial step length scaling
		number m_lambdaStart;

	/// reduction factor for the step length
		number m_lambdaReduce;
		
	///	sufficient descent factor
		number m_alpha;

	/// maximum allowed defect
		number m_maxDefect;

	/// verbosity level
		bool m_verbose;

	///	accept best
		bool m_bAcceptBest;

	///	check all
		bool m_bCheckAll;

	/// number of spaces inserted before output
		std::string m_offset;

//#if ENABLE_NESTED_NEWTON_RESOLFUNC_UPDATE
private:
		SmartPtr<NewtonUpdaterGeneric<TVector> > m_newtonUpdater;
//#endif
};

} // end namespace ug

#endif /* __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__LINE_SEARCH__ */
