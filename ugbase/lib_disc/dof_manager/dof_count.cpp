/*
 * dof_count.cpp
 *
 *  Created on: 06.11.2013
 *      Author: andreasvogel
 */

#include "dof_count.h"
#ifdef UG_PARALLEL
	#include "pcl/pcl_process_communicator.h"
#endif

namespace ug{

////////////////////////////////////////////////////////////////////////////////
// DoFCount
////////////////////////////////////////////////////////////////////////////////

DoFCount::DoFCount(const GridLevel& gl, ConstSmartPtr<DoFDistributionInfo> spDDInfo)
	:	DoFDistributionInfoProvider(spDDInfo), m_gridLevel(gl)
{
	vvCmpSubset.resize(num_fct());
	for(size_t fct = 0; fct < vvCmpSubset.size(); ++fct)
		vvCmpSubset[fct].resize(num_subsets());
}

void DoFCount::add(int fct, int si, SurfaceView::SurfaceState ss, byte is, uint64 numDoF)
{
	if(!(fct < (int)vvCmpSubset.size()))
		UG_THROW("DoFCount: fct index "<<fct<<" invalid. NumFct: "<<vvCmpSubset.size());

	if(!(si < (int)vvCmpSubset[fct].size()))
		UG_THROW("DoFCount: subset index "<<si<<" invalid. NumSubset: "<<vvCmpSubset[fct].size());

	vvCmpSubset[fct][si].add(numDoF, ss, is);
}

void DoFCount::sum_values_over_procs(int proc)
{
	PROFILE_FUNC();
#ifdef UG_PARALLEL
	pcl::ProcessCommunicator commWorld;

	// collect all serial values
	std::vector<uint64> vNumLocal;
	collect_values(vNumLocal);

	// sum
	if(proc == -1 || true){
		std::vector<uint64> vNumGlobal(vNumLocal.size());
		commWorld.allreduce(&vNumLocal[0], &vNumGlobal[0], vNumLocal.size(),
							PCL_DT_UNSIGNED_LONG_LONG, PCL_RO_SUM);
		set_values(vNumGlobal);
	}
	else{
		UG_THROW("Implement reduce")
	}
#endif
}

void DoFCount::collect_values(std::vector<uint64>& vNum) const
{
	PROFILE_FUNC();
	for(size_t fct = 0; fct < vvCmpSubset.size(); ++fct)
		for(size_t si = 0; si < vvCmpSubset[fct].size(); ++si)
			vvCmpSubset[fct][si].collect_values(vNum);
}

void DoFCount::set_values(const std::vector<uint64>& vNum)
{
	PROFILE_FUNC();
	size_t cnt = 0;
	for(size_t fct = 0; fct < vvCmpSubset.size(); ++fct)
		for(size_t si = 0; si < vvCmpSubset[fct].size(); ++si)
			vvCmpSubset[fct][si].set_values(vNum, cnt);
}

uint64 DoFCount::num(int fct, int si, SurfaceView::SurfaceState ss, byte is) const
{
	if(fct == ALL_FCT){
		uint64 cnt = 0;
		for(size_t fct = 0; fct < vvCmpSubset.size(); ++fct)
			cnt += num(fct, si, ss, is);
		return cnt;
	}

	if(si == ALL_SUBSET){
		uint64 cnt = 0;
		for(size_t si = 0; si < vvCmpSubset[fct].size(); ++si)
			cnt += num(fct, si, ss, is);
		return cnt;
	}

	return vvCmpSubset[fct][si].num(ss, is);
}

uint64 DoFCount::num_contains(int fct, int si, SurfaceView::SurfaceState ss, byte is) const
{
	if(fct == ALL_FCT){
		uint64 cnt = 0;
		for(size_t fct = 0; fct < vvCmpSubset.size(); ++fct)
			cnt += num_contains(fct, si, ss, is);
		return cnt;
	}

	if(si == ALL_SUBSET){
		uint64 cnt = 0;
		for(size_t si = 0; si < vvCmpSubset[fct].size(); ++si)
			cnt += num_contains(fct, si, ss, is);
		return cnt;
	}

	return vvCmpSubset[fct][si].num_contains(ss, is);
}


////////////////////////////////////////////////////////////////////////////////
// DoFCount::Cnt
////////////////////////////////////////////////////////////////////////////////

DoFCount::Cnt::Cnt()
{
	vNumSS.resize(DoFCount::SS_MAX + 1);
}

void DoFCount::Cnt::collect_values(std::vector<uint64>& vNum) const
{
	PROFILE_FUNC();
	for(size_t i = 0; i < vNumSS.size(); ++i)
		vNumSS[i].collect_values(vNum);
}

void DoFCount::Cnt::set_values(const std::vector<uint64>& vNum, size_t& cnt)
{
	PROFILE_FUNC();
	for(size_t i = 0; i < vNumSS.size(); ++i)
		vNumSS[i].set_values(vNum, cnt);
}


void DoFCount::Cnt::add(uint64 num, SurfaceView::SurfaceState ss, byte is)
{
	// restrict to only considered flags:
	size_t ss_index = (size_t)(ss & SS_MAX)();
	size_t is_index = (size_t)(is & ES_MAX);

	if(!(ss_index < vNumSS.size()))
		UG_THROW("Something wrong with surface state storage: is: "<<
		         ss_index<<", max: "<<vNumSS.size());

	if(!(is_index < vNumSS[ss_index].vNumIS.size()))
		UG_THROW("Something wrong with interface state storage: is: "<<
		         is_index<<", max: "<<vNumSS[ss_index].vNumIS.size());

	vNumSS[ss_index].vNumIS[is_index] += num;
}


uint64 DoFCount::Cnt::num(SurfaceView::SurfaceState ss, byte is) const
{
	if(ss == ALL_SS){
		uint64 cnt = 0;
		for(byte l = 0; l <= SS_MAX; ++l)
			cnt += num(l, is);
		return cnt;
	}

	if(ss == UNIQUE_SS){
		return num(SurfaceView::UNASSIGNED, is)
				+ num(SurfaceView::PURE_SURFACE, is)
				+ num(SurfaceView::SHADOWING, is)
	//			+ num(SurfaceView::SHADOW_COPY, is) 	// copies not counted
				+ num(SurfaceView::SHADOW_NONCOPY, is);
	}

	return vNumSS[ss()].num(is);
}

uint64 DoFCount::Cnt::num_contains(SurfaceView::SurfaceState ss, byte is) const
{
	if(ss == ALL_SS){
		uint64 cnt = 0;
		for(byte l = 0; l <= SS_MAX; ++l)
			cnt += vNumSS[l].num_contains(is);
		return cnt;
	}

	if(ss == UNIQUE_SS){
		return vNumSS[SurfaceView::UNASSIGNED].num_contains(is)
				+ vNumSS[SurfaceView::PURE_SURFACE].num_contains(is)
				+ vNumSS[SurfaceView::SHADOWING].num_contains(is)
	//			+ vNumSS[SurfaceView::SHADOW_COPY].num_contains(is) // copies not counted
				+ vNumSS[SurfaceView::SHADOW_NONCOPY].num_contains(is);
	}

	uint64 cnt = 0;
	for(byte l = 0; l <= SS_MAX; ++l){
		if(l & is)
			cnt += vNumSS[l].num_contains(is);
	}
	return cnt;
}

////////////////////////////////////////////////////////////////////////////////
// DoFCount::Cnt::PCnt
////////////////////////////////////////////////////////////////////////////////

uint64 DoFCount::Cnt::PCnt::num(byte is) const
{
	if(is == ALL_ES){
		uint64 cnt = 0;
		for(byte l = 0; l <= ES_MAX; ++l)
			cnt += num(l);
		return cnt;
	}

	if(is == UNIQUE_ES){
		return 	num(ES_NONE) 					// one-proc-only dofs
				+ num(ES_V_SLAVE)				// pure vert. slaves
				+ num_contains(ES_H_MASTER); 	// all horiz. masters
	}

	return vNumIS[is];
}

uint64 DoFCount::Cnt::PCnt::num_contains(byte is) const
{
	if(is == ALL_ES){
		uint64 cnt = 0;
		for(byte l = 0; l <= ES_MAX; ++l)
			cnt += vNumIS[l];
		return cnt;
	}

	if(is == UNIQUE_ES){
		return 	num(ES_NONE) 					// one-proc-only dofs
				+ num(ES_V_SLAVE)				// pure vert. slaves
				+ num_contains(ES_H_MASTER); 	// all horiz. masters
	}

	uint64 cnt = 0;
	for(byte l = 0; l <= ES_MAX; ++l){
		if(l & is)
			cnt += vNumIS[l];
	}
	return cnt;
}

DoFCount::Cnt::PCnt::PCnt()
{
	vNumIS.resize(DoFCount::ES_MAX + 1, 0);
}

void DoFCount::Cnt::PCnt::collect_values(std::vector<uint64>& vNum) const
{
	PROFILE_FUNC();
	for(byte l = 0; l <= ES_MAX; ++l){
//		if(l & (ES_V_MASTER | ES_V_SLAVE)) continue;
//		if(l & (ES_H_MASTER | ES_H_SLAVE)) continue;

		vNum.push_back(vNumIS[l]);
	}
}

void DoFCount::Cnt::PCnt::set_values(const std::vector<uint64>& vNum, size_t& cnt)
{
	PROFILE_FUNC();
	size_t used = 0;
	for(byte l = 0; l <= ES_MAX; ++l){
//		if(l & (ES_V_MASTER | ES_V_SLAVE)) continue;
//		if(l & (ES_H_MASTER | ES_H_SLAVE)) continue;

		vNumIS[l] = vNum[cnt + l]; ++used;
	}

	cnt += used;
}

} // end namespace ug
