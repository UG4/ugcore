/*
 * line_search.h
 *
 *  Created on: 18.07.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__LINE_SEARCH__
#define __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__LINE_SEARCH__

#include <ostream>
#include <string>
#include <vector>
#include <cmath>

#include "common/common.h"

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

	/// virtual destructor
		virtual ~ILineSearch() {}
};

/// standard implementation for line search
template <typename TVector>
class StandardLineSearch : public ILineSearch<TVector>
{
	public:
	//	type of
		typedef TVector vector_type;

	public:
	///	default constructor (setting default values)
		StandardLineSearch()
		 :	 m_maxSteps(10), m_lambdaStart(1.0), m_lambdaReduce(0.5),
		  	 m_maxDefect(1e+10), m_verbose(true), m_bAcceptBest(false), m_bCheckAll(false), m_offset("")
			 {};

	///	constructor
		StandardLineSearch(int maxSteps, number lambdaStart, number lambdaReduce, bool bAcceptBest)
		 :	 m_maxSteps(maxSteps), m_lambdaStart(lambdaStart), m_lambdaReduce(lambdaReduce),
			 m_maxDefect(1e+10), m_verbose(true), m_bAcceptBest(bAcceptBest), m_bCheckAll(false), m_offset("")
			 {};

	///	constructor
		StandardLineSearch(int maxSteps, number lambdaStart, number lambdaReduce, bool bAcceptBest, bool bCheckAll)
		 :	 m_maxSteps(maxSteps), m_lambdaStart(lambdaStart), m_lambdaReduce(lambdaReduce),
			 m_maxDefect(1e+10), m_verbose(true), m_bAcceptBest(bAcceptBest), m_bCheckAll(bCheckAll), m_offset("")
			 {};

	///	sets maximum number of line search steps
		void set_maximum_steps(int steps) {m_maxSteps = steps;}

	///	sets start factor
		void set_lambda_start(number start) {m_lambdaStart = start;}

	///	sets factor by which line search factor is multiplied in each step
		void set_reduce_factor(number factor) {m_lambdaReduce = factor;}

	///	sets iff after max_steps the best try is used
		void set_accept_best(bool bAcceptBest) {m_bAcceptBest = bAcceptBest;}

	///	sets iff after max_steps the best try is used
		void set_check_all(bool bCheckAll) {m_bCheckAll = bCheckAll;}

	///	sets maximum allowed defect
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
			number alpha = 0.25;

		//	some values
			number norm, norm_old = defect;
			bool converged = false;
			std::vector<number> vRho;

		// remember u
			s = u;



		// check if defect-norm is already smaller than maximum allowed defect value
		/*if (norm_old < m_maxDefect)
		{
			UG_LOG("ERROR in 'StandardLineSearch::search':"
					" no computation required.\n");
									return true;
		}
*/

		//	print heading line
		if(m_verbose)
			UG_LOG(m_offset << "   ++++ Line Search:  Iter       lambda        Defect          Rate \n");


		//	loop line search steps
			for(int k = 1; k <= m_maxSteps; ++k)
			{
			// 	try on line u := u - lambda*p
				VecScaleAdd(u, 1.0, u, (-1)*lambda, p);

			// 	compute new Defect
				spOp->prepare(u);
				spOp->apply(d, u);

			//	compute new Residuum
				norm = d.norm();

			// 	compute reduction
				vRho.push_back(norm/norm_old);

			//	print rate
				if(m_verbose)
					UG_LOG(m_offset << "   +                 " << std::setw(4)
							<< k << ":   " << std::setw(11)
							<< std::resetiosflags( ::std::ios::scientific )<<
							lambda << "     "
							<< std::scientific << norm << "   " << vRho.back() <<"\n");

			// 	check if reduction fits
				if(vRho.back() <= 1 - alpha * fabs(lambda))
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
					UG_LOG(m_offset << "   ++++ Accept Best: Accepting step " <<
					       best+1 << ", Rate = "<< vRho[best] <<".\n");

				// 	try on line u := u - lambda*p
					VecScaleAdd(u, 1.0, s, (-1)*m_lambdaStart*std::pow(m_lambdaReduce, (number)best), p);

				// 	compute new Defect
					spOp->prepare(u);
					spOp->apply(d, u);

					// compute new Residuum
					norm = d.norm();

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

			// 	reset u
				u = s;
			}

		//	print end line
			if(m_verbose)
			{
				//only for rate < 1, we call it "Line Search converged"
				if(converged)
					UG_LOG(m_offset << "   ++++ Line Search converged.\n");
			}

		//	we're done
			return true;
		}

	protected:
	/// solution in line direction
		vector_type s;

	protected:
	/// maximum number of steps to be performed
		int m_maxSteps;

	/// lambda start
		number m_lambdaStart;

	/// lambda reduce
		number m_lambdaReduce;

	/// maximum allowed defect
		number m_maxDefect;

	/// verbose level
		bool m_verbose;

	///	accept best
		bool m_bAcceptBest;

	///	check all
		bool m_bCheckAll;

	/// number of spaces inserted before output
		std::string m_offset;
};

} // end namespace ug

#endif /* __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__LINE_SEARCH__ */
