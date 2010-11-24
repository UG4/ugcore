/*
 * convergence_check.h
 *
 *  Created on: 18.07.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__OPERATOR__NON_LINEAR_OPERATOR__LINE_SEARCH__
#define __H__LIB_DISCRETIZATION__OPERATOR__NON_LINEAR_OPERATOR__LINE_SEARCH__

#include <ostream>
#include <string>

#include "common/common.h"
#include "lib_discretization/operator/operator.h"

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
template <typename TFunction>
class ILineSearch
{
	public:
		typedef TFunction function_type;

	public:
		/// set string to be printed before each output of line search
		virtual void set_offset(std::string offset) = 0;

		/** search
		 *
		 *	Performs a line search to a given direction.
		 *
		 * \param[in]		Op		Non-linear operator
		 * \param[in]  		u 		current solution
		 * \param[in]		p		search direction
		 * \param[in|out]	d		defect
		 * \param[in]		defect	norm of current defect
		 *
		 * \return 	true 		if line search successful
		 * 			false 		if line search failed
		 */
		virtual bool search(IOperator<function_type, function_type>& Op, function_type& u, function_type& p, function_type& d, number defect) = 0;

		/// virtual destructor
		virtual ~ILineSearch() {}
};

template <typename TFunction>
class StandardLineSearch : public ILineSearch<TFunction>
{
	public:
		typedef TFunction function_type;

	public:
		StandardLineSearch(int maxSteps, number lambdaStart, number lambdaReduce, bool verbose)
		 :	 m_maxSteps(maxSteps), m_lambdaStart(lambdaStart), m_lambdaReduce(lambdaReduce),
			 m_verbose(verbose), m_offset("")
			 {};

		StandardLineSearch()
		 :	 m_maxSteps(10), m_lambdaStart(1.0), m_lambdaReduce(0.5),
			 m_verbose(true), m_offset("")
			 {};

		void set_maximum_steps(int steps) {m_maxSteps = steps;}
		void set_lambda_start(number start) {m_lambdaStart = start;}
		void set_reduce_factor(number factor) {m_lambdaReduce = factor;}
		void set_verbose_level(bool level) {m_verbose = level;}
		void set_offset(std::string offset) {m_offset = offset;};

		bool search(IOperator<function_type, function_type>& Op, function_type& u, function_type& p, function_type& d, number defect)
		{
			// clone pattern for s
			s.resize(u.size());

			number lambda = m_lambdaStart;
			number alpha = 0.25;

			number norm, norm_old;
			norm_old = defect;

			// remember u
			s = u;

			if(m_verbose)
				UG_LOG(m_offset << "   ++++ Line Search:  Iter       lambda        Defect          Rate \n");

			for(int k = 1; k <= m_maxSteps; ++k)
			{
				// try on line u := u - lambda*p
				VecScaleAdd(u, 1.0, u, (-1)*lambda, p);

				// compute new Defect
				if(!Op.prepare(d, u))
					{UG_LOG("StandardLineSearch: Cannot prepare Non-linear"
							" Operator for defect computation.\n"); return false;}
				if(!Op.apply(d, u))
					{UG_LOG("StandardLineSearch: Cannot apply Non-linear Operator "
							"to compute defect.\n"); return false;}

				//compute new Residuum
				norm = d.two_norm();

				// compute reduction
				number rho = norm/norm_old;

				// print rate
				if(m_verbose)
					UG_LOG(m_offset << "   +                 " << std::setw(4) << k << ":   " << std::setw(10)
									<< std::resetiosflags( ::std::ios::scientific )<< lambda << "     "
									<< std::scientific << norm << "   " << rho <<"\n");

				// check if reduction fits
				if(rho <= 1 - alpha * fabs(lambda)) break;
				else lambda *= m_lambdaReduce;

				if(k == m_maxSteps)
					{UG_LOG(m_offset << "   ++++ Line Search did not converge.\n"); return false;}

				// reset u
				u = s;
			}
			if(m_verbose)
				UG_LOG(m_offset << "   ++++ Line Search converged.\n");
			return true;
		}

	protected:
		// solution in line direction
		function_type s;

	protected:
		// maximum number of steps to be performed
		int m_maxSteps;

		// lambda start
		number m_lambdaStart;

		// lambda reduce
		number m_lambdaReduce;

		// verbose level
		bool m_verbose;

		// number of spaces inserted before output
		std::string m_offset;
};

} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__OPERATOR__NON_LINEAR_OPERATOR__LINE_SEARCH__ */
