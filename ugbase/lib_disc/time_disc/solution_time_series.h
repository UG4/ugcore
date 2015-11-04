
#ifndef __H__UG__LIB_DISC__TIME_DISC__PREVIOUS_SOLUTION__
#define __H__UG__LIB_DISC__TIME_DISC__PREVIOUS_SOLUTION__

// extern libraries
#include <deque>

// other ug libraries
#include "common/common.h"

#include "lib_disc/common/local_algebra.h"

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
class VectorTimeSeries
{
	public:
	///	vector type of solutions
		typedef TVector vector_type;

	public:
		virtual ~VectorTimeSeries() {}

	//! clones the object (deep-copy) including values
		SmartPtr<VectorTimeSeries<vector_type> > clone() const
		{
			SmartPtr<VectorTimeSeries<TVector> > cloneTimeSol
				= SmartPtr<VectorTimeSeries<TVector> >(new VectorTimeSeries<TVector>);

			for(int i = m_vTimeSol.size()-1; i >= 0 ; --i)
				cloneTimeSol->push(solution(i)->clone(), time(i));

			return cloneTimeSol;
		}

	/// clears the content of the member m_vTimeSol
		void clear() {m_vTimeSol.clear();}

	///	returns number of time steps handled
		size_t size() const {return m_vTimeSol.size();}

	///	returns point in time for solution
		number time(size_t i) const {return m_vTimeSol.at(i).time();}

	///	returns solution
		SmartPtr<vector_type> solution(size_t i) {return m_vTimeSol.at(i).solution();}

	///	returns solution
		ConstSmartPtr<vector_type> solution(size_t i) const {return m_vTimeSol.at(i).solution();}

	///	returns oldest solution
		SmartPtr<vector_type> oldest() {return m_vTimeSol.back().solution();}

	/// const access to oldest solution
		ConstSmartPtr<vector_type> oldest() const {return m_vTimeSol.back().solution();}

	///	returns latest solution
		SmartPtr<vector_type> latest() {return m_vTimeSol.front().solution();}

	///	const access to latest solution
		ConstSmartPtr<vector_type> latest() const {return m_vTimeSol.front().solution();}

	///	adds new time point, not discarding the oldest
		void push(SmartPtr<vector_type> vec, number time) {m_vTimeSol.push_front(TimeSol(vec, time));}

	///	adds new time point, oldest solution is discarded and returned
		SmartPtr<vector_type> push_discard_oldest(SmartPtr<vector_type> vec, number time)
		{
			SmartPtr<vector_type> discardVec = m_vTimeSol.back().solution();
			remove_oldest(); push(vec, time);
			return discardVec;
		}

	///	removes latest time point
		void remove_latest() {m_vTimeSol.pop_front();}

	///	removes oldest time point
		void remove_oldest() {m_vTimeSol.pop_back();}

	protected:
	///	grouping of solution and time point
		class TimeSol
		{
			public:
				TimeSol() : vec(NULL), t(0.0) {}

				TimeSol(SmartPtr<vector_type> vec_, number t_)
					: vec(vec_), t(t_) {}

			///	access solution
				SmartPtr<vector_type> solution() {return vec;}

			///	const access solution
				ConstSmartPtr<vector_type> solution() const {return vec;}

			///	access time
				number& time() {return t;}

			///	const access time
				const number& time() const {return t;}

			protected:
			//	solution vector at time point
				SmartPtr<vector_type> vec;

			//	point in time
				number t;
		};

	protected:

	//	deque of previous solutions
		std::deque<TimeSol> m_vTimeSol;
};

/// time series of local vectors
class LocalVectorTimeSeries
{
	public:
	///	constructor
		LocalVectorTimeSeries() {}

	///	returns number of time points
		size_t size() const {return m_vLocalVector.size();}

	///	returns time point i
		number time(size_t i) const {return m_vTime.at(i);}

	///	returns time points
		const std::vector<number>& times() const {return m_vTime;}

	///	returns the local vector for the i'th time point
		const LocalVector& solution(size_t i) const {return m_vLocalVector.at(i);}

	///	returns the local vector for the i'th time point
		LocalVector& solution(size_t i) {return m_vLocalVector.at(i);}

	///	access dofs by map
		void access_by_map(const FunctionIndexMapping& funcMap)
		{
			for(size_t t=0; t < size(); ++t)
				solution(t).access_by_map(funcMap);
		}

	/// extract local values from global vectors
		template <typename TVector>
		void read_values(ConstSmartPtr<VectorTimeSeries<TVector> > vecTimeSeries, LocalIndices& ind)
		{
			m_vLocalVector.resize(vecTimeSeries->size());
			for(size_t i = 0; i < m_vLocalVector.size(); ++i)
			{
				m_vLocalVector[i].resize(ind);
				GetLocalVector(m_vLocalVector[i], *vecTimeSeries->solution(i));
			}
		}

	///	extract time points
		template <typename TVector>
		void read_times(ConstSmartPtr<VectorTimeSeries<TVector> > vecTimeSeries)
		{
			m_vTime.resize(vecTimeSeries->size());
			for(size_t i = 0; i < m_vTime.size(); ++i)
				m_vTime[i] = vecTimeSeries->time(i);
		}

	protected:
	///	time series
		std::vector<number> m_vTime;

	///	vector of local vectors (one for each time point)
		std::vector<LocalVector> m_vLocalVector;
};

/// @}

} // end namespace ug

#endif /* __H__UG__LIB_DISC__TIME_DISC__PREVIOUS_SOLUTION__ */
