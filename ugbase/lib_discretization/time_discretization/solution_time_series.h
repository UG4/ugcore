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

#include "lib_discretization/common/local_algebra.h"

namespace ug{

/// \ingroup lib_disc_time_assemble
/// @{

/// time series of solutions and corresponding time point
/**
 * This class holds solutions and corresponding points in time. It is
 * intended to group previous computed solutions for a time stepping scheme, such
 * that previous steps can be passed to a time stepping scheme at once. Internally,
 * this object is basically a deque of old solutions, and adding a newly
 * computed solution lets the object pop the oldest stored solution.
 */
template <typename TVector>
class SolutionTimeSeries
{
	public:
	///	vector type of solutions
		typedef TVector vector_type;

	public:

	///	returns number of time steps handled
		size_t size() const {return m_vTimeSolution.size();}

	///	returns point in time for solution
		number time(size_t i) const {return m_vTimeSolution.at(i).time();}

	///	returns solution
		vector_type& solution(size_t i) {return *(m_vTimeSolution.at(i).solution());}

	///	returns solution
		const vector_type& solution(size_t i) const {return *(m_vTimeSolution.at(i).solution());}

	///	returns oldest solution
		vector_type& oldest() {return *(m_vTimeSolution.back().solution());}

	/// const access to oldest solution
		const vector_type& oldest() const {return *(m_vTimeSolution.back().solution());}

	///	returns latest solution
		vector_type& latest() {return *(m_vTimeSolution.front().solution());}

	///	const access to latest solution
		const vector_type& latest() const {return *(m_vTimeSolution.front().solution());}

	///	adds new time point, not discarding the oldest
		void push(vector_type& vec, number time) {m_vTimeSolution.push_front(TimeSol(vec, time));}

	///	adds new time point, oldest solution is discarded and returned
		vector_type* push_discard_oldest(vector_type& vec, number time)
		{
			vector_type* discardVec = m_vTimeSolution.back().solution();
			remove_oldest(); push(vec, time);
			return discardVec;
		}

	///	removes latest time point
		void remove_latest() {m_vTimeSolution.pop_front();}

	///	removes oldest time point
		void remove_oldest() {m_vTimeSolution.pop_back();}

	protected:
	///	grouping of solution and time point
		class TimeSol
		{
			public:
				TimeSol() : vec(NULL), t(0.0) {}

				TimeSol(vector_type& vec_, number t_)
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
		std::deque<TimeSol> m_vTimeSolution;
};


template <typename TVector>
class LocalVectorTimeSeries
{
	public:
	///	Vector type
		typedef TVector vector_type;

	///	Entry type
		typedef typename vector_type::value_type value_type;

	public:
	///	constructor
		LocalVectorTimeSeries(const SolutionTimeSeries<TVector>& solTimeSeries)
			: m_rSolTimeSeries(solTimeSeries)
		{}

	///	returns number of time points
		size_t size() const {return m_rSolTimeSeries.size();}

	///	returns time point i
		number time(size_t i) const {return m_rSolTimeSeries.time(i);}

	///	returns the local vector for the i'th time point
		const LocalVector<value_type>& solution(size_t i) const
			{return m_vLocalVector.at(i);}

	///	update for indices
		void read_values(LocalIndices& ind)
		{
			for(size_t i = 0; i < m_vLocalVector.size(); ++i)
			{
				m_vLocalVector[i].set_indices(ind);
				m_vLocalVector[i].read_values(m_rSolTimeSeries.solution(i));
			}
		}

	protected:
	///	time series
		const SolutionTimeSeries<TVector>& m_rSolTimeSeries;

	///	vector of local vectors (one for each time point)
		std::vector<LocalVector<value_type> > m_vLocalVector;
};

/// @}

} // end namespace ug

#endif /* __H__UG__LIB_DISCRETIZATION__TIME_DISCRETIZATION__PREVIOUS_SOLUTION__ */
