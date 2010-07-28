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
#include "lib_discretization/operator/operator_interface.h"

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
class LineSearch
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
		 * \param[in/out]	d		defect
		 * \param[in]		defect	norm of current defect
		 *
		 * \return 	true 		if line search successful
		 * 			false 		if line search failed
		 */
		virtual bool search(IOperator<function_type, function_type>& Op, function_type& u, function_type& p, function_type& d, number defect) = 0;

		/// virtual destructor
		virtual ~LineSearch() {}
};

template <typename TFunction>
class StandardLineSearch : public LineSearch<TFunction>
{
	public:
		typedef TFunction function_type;

	public:
		StandardLineSearch(int maxSteps, number lambdaStart, number lambdaReduce, bool verbose)
		 :	 m_maxSteps(maxSteps), m_lambdaStart(lambdaStart), m_lambdaReduce(lambdaReduce),
			 m_verbose(verbose), m_offset("")
			 {};

		void set_offset(std::string offset) {m_offset = offset;};

		bool search(IOperator<function_type, function_type>& Op, function_type& u, function_type& p, function_type& d, number defect)
		{
			// clone pattern for s
			s.clone_pattern(u);

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
				VecScaleAppend(u, p, (-1)*lambda);

				// compute new Defect
				if(!Op.prepare(u, d))
					{UG_LOG("StandardLineSearch: Cannot prepare Non-linear"
							" Operator for defect computation.\n"); return false;}
				if(!Op.apply(u, d))
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

		bool VecScaleAppend(TFunction& a_func, TFunction& b_func, number s)
		{
			#ifdef UG_PARALLEL
			if(a_func.has_storage_type(PST_UNIQUE) && b_func.has_storage_type(PST_UNIQUE));
			else if(a_func.has_storage_type(PST_CONSISTENT) && b_func.has_storage_type(PST_CONSISTENT));
			else if (a_func.has_storage_type(PST_ADDITIVE) && b_func.has_storage_type(PST_ADDITIVE))
			{
				a_func.set_storage_type(PST_ADDITIVE);
				b_func.set_storage_type(PST_ADDITIVE);
			}
			#endif
			typename TFunction::vector_type& a = a_func.get_vector();
			typename TFunction::vector_type& b = b_func.get_vector();

            for(size_t i = 0; i < a.size(); ++i)
				a[i] += s*b[i];
			return true;
		}

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
