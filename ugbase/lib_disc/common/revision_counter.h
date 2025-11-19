/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
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
		RevisionCounter() : m_pObj(nullptr), m_cnt(0) {};

	///	constructor (with valid state initialization)
		RevisionCounter(const void* pObj) : m_pObj(pObj), m_cnt(1) {}

	///	constructor (with valid state initialization)
		template <typename T>
		RevisionCounter(const T* pObj) : m_pObj(static_cast<const void*>(pObj)), m_cnt(1) {}

	///	increase state (prefix)
		RevisionCounter& operator ++ () {
			if(invalid())
				UG_THROW("AdaptState: increasing invalid state not admissible.")

			++m_cnt;

			if(invalid())
				UG_THROW("AdaptState: counter overflow. Alter implementation.")
			return *this;
		}

	///	increase state (postfix)
		RevisionCounter operator ++ (int) {
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
		bool valid() const {return m_pObj != nullptr && m_cnt != 0;}

	///	returns if state is invalid
		bool invalid() const {return !valid();}

	///	invalidates state
		void invalidate() {m_pObj = nullptr; m_cnt = 0;}

	///	returns the associated object
		const void* obj() const {return m_pObj;}

	protected:
		const void* m_pObj; ///< associated object
		uint64 m_cnt; ///< state counter (0 = invalid)
};

} // end namespace ug

#endif