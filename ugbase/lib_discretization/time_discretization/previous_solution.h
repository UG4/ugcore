/*
 * previous_solution.h
 *
 *  Created on: 23.10.2009
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISCRETIZATION__TIME_DISCRETIZATION__PREVIOUS_SOLUTION__
#define __H__UG__LIB_DISCRETIZATION__TIME_DISCRETIZATION__PREVIOUS_SOLUTION__

// extern libraries
#include <deque>

// other ug libraries
#include "common/common.h"

namespace ug{

/// \ingroup lib_disc_time_assemble
/// @{

/// old solutions and time steps
/**
 * This class holds solutions and corresponding points in time. It is
 * intended to group previous computed solutions for a time stepping scheme, such
 * that previous steps can be passed to a time stepping scheme at once. Internally,
 * this object is basically a deque of old solutions, and adding a newly
 * computed solution lets the object pop the oldest stored solution.
 */
template <typename TVector>
class PreviousSolutions
{
	public:
	///	vector type of solutions
		typedef TVector vector_type;

	public:

	///	returns number of previous time steps handled
		size_t size() const {return m_vPreviousSolution.size();}

	///	adds new time point, oldest solution is discarded and returned
		vector_type* push_discard_oldest(vector_type& vec, number time)
		{
			vector_type* discardVec = m_vPreviousSolution.back().solution();
			m_vPreviousSolution.pop_back();
			m_vPreviousSolution.push_front(PreviousSolution(vec, time));
			return discardVec;
		}

	///	adds new time point, not discarding the oldest
		void push(vector_type& vec, number time)
		{
			m_vPreviousSolution.push_front(PreviousSolution(vec, time));
		}

	///	returns previous solution
		const vector_type& solution(size_t i) const {return *(m_vPreviousSolution.at(i).solution());}

	///	returns previous solution
		vector_type& solution(size_t i) {return *(m_vPreviousSolution.at(i).solution());}

	///	returns point in time for previous solution
		number time(size_t i) const {return m_vPreviousSolution.at(i).time();}

	///	returns oldest solution
		vector_type& oldest_solution() {return *(m_vPreviousSolution.back().solution());}

	protected:
		class PreviousSolution
		{
			public:
				PreviousSolution() : vec(NULL), t(0.0) {}

				PreviousSolution(vector_type& vec_, number t_)
					: vec(&vec_), t(t_) {}

			///	access solution
				vector_type* solution() {return vec;}

			///	const access solution
				const vector_type* solution() const {return vec;}

			///	access time
				number& time() {return t;}

			///	const access time
				const number& time() const {return t;}

			protected:
			//	solution vector at time point
				vector_type* vec;

			//	point in time
				number t;
		};

	//	deque of previous solutions
		std::deque<PreviousSolution> m_vPreviousSolution;
};

/// @}

} // end namespace ug

#endif /* __H__UG__LIB_DISCRETIZATION__TIME_DISCRETIZATION__PREVIOUS_SOLUTION__ */
