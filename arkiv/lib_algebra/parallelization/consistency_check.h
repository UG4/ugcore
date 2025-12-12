/*
 * Copyright (c) 2012-2014:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
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



#ifndef CONSISTENCY_CHECK_H
#define CONSISTENCY_CHECK_H

#ifdef UG_PARALLEL

#include <string>
#include "pcl/pcl.h"
#include "common/log.h"
#include "common/assert.h"
#include "communication_scheme.h"
#include "lib_algebra/parallelization/parallel_index_layout.h" // for IndexLayout

namespace ug {

/**
 * \defgroup lib_algebra_parallel_consistencycheck Parallel Algebra Consistency Check
 * \ingroup lib_algebra_parallelization
 * \brief with ConsistencyCheck, you can check the consistency of any array over layouts
 * \{
 */

template<typename TVec, typename TValue>
class ConsistencyCheckClassSend
{
public:
	ConsistencyCheckClassSend(const TVec &_vec) : vec(_vec) {}
	const TValue &send(int pid, int index) const
	{
		return vec[index];
	}

	const TVec &vec;
};

// for std::vector<bool> and others
template<typename TVec>
class ConsistencyCheckClassSend<TVec, bool>
{
public:
	ConsistencyCheckClassSend(const TVec &_vec) : vec(_vec) {}
	bool send(int pid, int index) const
	{
		return vec[index];
	}

	const TVec &vec;
};

template<typename TVec, typename TValue>
class ConsistencyCheckClass
: public CommunicationScheme<ConsistencyCheckClass<TVec, TValue>, TValue>,
  public ConsistencyCheckClassSend<TVec, TValue>
{
public:
	using ConsistencyCheckClassSend<TVec, TValue>::vec;
	ConsistencyCheckClass(const TVec &_vec, std::string _name = "") :
		ConsistencyCheckClassSend<TVec, TValue>(_vec), name(_name)
	{
		bOK = true;
	}

	void receive(int pid, int index, TValue &v)
	{
		if(vec[index] != v)
		{
			if(bOK)
			{ UG_LOG("\n\n----------------------\n" << name << " not consistent:\n"); bOK = false; }
			UG_LOG("index " << index << " is " << vec[index] << " on this proc (" << pcl::ProcRank() <<
					"), but " << v << " on master (proc " << pid << ").\n");
		}
	}

	bool isOK()
	{
		return bOK;
	}

	int get_element_size()
	{
		if(block_traits<TValue>::is_static) return sizeof(TValue);
		else return -1;
	}
private:
	bool bOK;

	std::string name;
};

/** ConsistencyCheck
 * \brief receives data over a interface based on a CommunicationScheme on a subgroup of processes
 * \tparam 	TVec			vector type to check
 * \param 	vec				vec to check
 * \param	com				InterfaceCommunicator used to send data
 * \param	pc				ProcessCommunicator used for AllProcsTrue
 * \param	masterLayout	layout to send data
 * \param	slaveLayout		layout to receive data and check if equal
 * \param	name			name to be given out if not consistent. default ""
 */
template<typename TVec>
void ConsistencyCheck(const TVec &vec, pcl::InterfaceCommunicator<IndexLayout> &com,
		const pcl::ProcessCommunicator &pc, const IndexLayout &masterLayout,
		const IndexLayout &slaveLayout, std::string name="")
{
	PROFILE_FUNC_GROUP("algebra parallelization debug");
	ConsistencyCheckClass<TVec, typename TVec::value_type> scheme(vec, name);
	CommunicateOnInterfaces(com, masterLayout, slaveLayout, scheme);

	UG_COND_THROW(!AllProcsTrue(scheme.isOK(), pc), name << " not consistent!");
}

template<typename TVec>
void ConsistencyCheck(const TVec &vec, const HorizontalAlgebraLayouts &layout, std::string name="")
{
	PROFILE_FUNC_GROUP("algebra parallelization debug");
	ConsistencyCheckClass<TVec, typename TVec::value_type> scheme(vec, name);
	CommunicateOnInterfaces(layout.comm(), layout.master(), layout.slave(), scheme);

	UG_COND_THROW(!AllProcsTrue(scheme.isOK(), layout.proc_comm()), name << " not consistent!");
}

// end group lib_algebra_parallel_consistencycheck
/// \}

}
#endif
#endif
