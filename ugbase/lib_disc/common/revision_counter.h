/*
 * revision_counter.h
 *
 *  Created on: 18.11.2013
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__COMMON__REVISION_COUNTER__
#define __H__UG__LIB_DISC__COMMON__REVISION_COUNTER__

#include "common/common.h"

namespace ug{

/// Class used to identify a state of adaption of a grid, approx-space, ...
/**
 * This class is used as a state counter of an object. E.g., this is very
 * useful to track the state of an adaptive multigrid. Based on this counter
 * reinitialization of dependent structures can be triggered, that should only
 * be reinitialized when actually used again (and not on every state change
 *  - use observers or msg listeners in that case).
 *
 * NOTE: the current implementation allows only std::numeric_limits<unit64>::max()
 * 		 states (approx 10^20 states), which should be pretty enough for most
 * 		 considered uses. If this is not enough a different implementation,
 * 		 handling counter overflow, should be used.
 * 		 (However, if the state is altered every nano-second, the overflow will
 * 		  appear after more than 3000 years)
 */
class RevisionCounter
{
	public:
	///	constructor (with invalid state initialization)
		RevisionCounter() : m_pObj(0), m_cnt(0) {};

	///	constructor (with valid state initialization)
		RevisionCounter(const void* pObj) : m_pObj(pObj), m_cnt(1) {}

	///	constructor (with valid state initialization)
		template <typename T>
		RevisionCounter(const T* pObj) : m_pObj(static_cast<const void*>(pObj)), m_cnt(1) {}

	///	increase state (prefix)
		RevisionCounter& operator++() {
			if(invalid())
				UG_THROW("AdaptState: increasing invalid state not admissible.")

			++m_cnt;

			if(invalid())
				UG_THROW("AdaptState: counter overflow. Alter implementation.")
			return *this;
		}

	///	increase state (postfix)
		RevisionCounter operator++(int) {
			RevisionCounter tmp(*this);
			++(*this);
			return tmp;
		}

	///	compare two states
		bool operator==(const RevisionCounter& rhs) const{
			if(invalid() || rhs.invalid()) return false;
			if(m_pObj != rhs.m_pObj) return false;
			return (m_cnt == rhs.m_cnt);
		}

	///	compare two states
		bool operator!=(const RevisionCounter& rhs) const{
			return !((*this) == rhs);
		}

	///	compare two states
		bool operator<(const RevisionCounter& rhs) const{
			if(m_pObj != rhs.m_pObj) return m_pObj < rhs.m_pObj;
			if(m_cnt != rhs.m_cnt) return m_cnt < rhs.m_cnt;
			return false;
		}

	///	compare two states
		bool operator>(const RevisionCounter& rhs) const{
			if(m_pObj != rhs.m_pObj) return m_pObj > rhs.m_pObj;
			if(m_cnt != rhs.m_cnt) return m_cnt > rhs.m_cnt;
			return false;
		}

	///	returns if state is valid
		bool valid() const {return m_pObj != 0 && m_cnt != 0;}

	///	returns if state is invalid
		bool invalid() const {return !valid();}

	///	invalidates state
		void invalidate() {m_pObj = 0; m_cnt = 0;}

	///	returns the associated object
		const void* obj() const {return m_pObj;}

	protected:
		const void* m_pObj; ///< associated object
		uint64 m_cnt; ///< state counter (0 = invalid)
};

} // end namespace ug

#endif /* __H__UG__LIB_DISC__COMMON__REVISION_COUNTER__ */
