#ifdef UG_PARALLEL

#ifndef CONSISTENCY_CHECK_H
#define CONSISTENCY_CHECK_H

#include <string>
#include "pcl/pcl.h"
#include "common/log.h"
#include "common/assert.h"
#include "communication_scheme.h"
#include "lib_algebra/parallelization/parallel_index_layout.h" // for IndexLayout

/**
 * with ConsistencyCheck, you can check the consistency of any array over layouts
*/
namespace ug{

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
			{ UG_LOG(name << " not consistent:\n"); bOK = false; }
			UG_LOG("index " << index << " is " << vec[index] << " on this proc (" << pcl::GetProcRank() <<
					", but " << v << " on master (proc " << pid << ".\n");
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

	UG_ASSERT(AllProcsTrue(scheme.isOK(), pc), name << " not consistent!");
}

template<typename TVec>
void ConsistencyCheck(const TVec &vec, const HorizontalAlgebraLayouts &layout, std::string name="")
{
	PROFILE_FUNC_GROUP("algebra parallelization debug");
	ConsistencyCheckClass<TVec, typename TVec::value_type> scheme(vec, name);
	CommunicateOnInterfaces(layout.comm(), layout.master(), layout.slave(), scheme);

	UG_ASSERT(AllProcsTrue(scheme.isOK(), layout.proc_comm()), name << " not consistent!");
}

}
#endif // CONSISTENCY_CHECK_H
#endif // UG_PARALLEL
