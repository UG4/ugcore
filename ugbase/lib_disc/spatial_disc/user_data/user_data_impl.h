/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__USER_DATA__USER_DATA_IMPL__
#define __H__UG__LIB_DISC__SPATIAL_DISC__USER_DATA__USER_DATA_IMPL__

#include "user_data.h"
#include "lib_disc/common/groups_util.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
//	ICplUserData
////////////////////////////////////////////////////////////////////////////////

template <int dim>
ICplUserData<dim>::ICplUserData()
:	m_locPosDim(-1), m_timePoint(0), m_defaultTimePoint(-1), m_si(-1)
{
	m_vNumIP.clear();
	m_vMayChange.clear();
	m_locPosDim = -1;
	m_pvLocIP1d.clear(); m_pvLocIP2d.clear(); m_pvLocIP3d.clear();
	m_vTime.clear(); m_vTime.push_back(0.0);
}

template <int dim>
void ICplUserData<dim>::clear()
{
	local_ip_series_to_be_cleared();
	m_vNumIP.clear();
	m_vMayChange.clear();
	m_locPosDim = -1;
	m_pvLocIP1d.clear(); m_pvLocIP2d.clear(); m_pvLocIP3d.clear();
	m_timePoint = 0;
	m_vTime.clear(); m_vTime.push_back(0.0);
	m_si = -1;
}

template <int dim>
template <int ldim>
size_t ICplUserData<dim>::register_local_ip_series(const MathVector<ldim>* vPos,
                                         const size_t numIP,
                                         const int timePointSpec,
                                         bool bMayChange)
{
//	check, that dimension is ok.
	if(m_locPosDim == -1) m_locPosDim = ldim;
	else if(m_locPosDim != ldim)
		UG_THROW("Local IP dimension conflict");
	
//	get the "right" time point specification
	int theTimePoint = (m_defaultTimePoint >= 0)? m_defaultTimePoint : timePointSpec;

//	get local positions
	std::vector<const MathVector<ldim>*>& vvIP = get_local_ips(Int2Type<ldim>());

//	search for ips
//	we only identify ip series if the local ip positions will not change
	if(!bMayChange && numIP != 0)
		for(size_t s = 0; s < vvIP.size(); ++s)
		{
		//	return series number iff exists and local ips remain constant
			if(!m_vMayChange[s])
				if(vvIP[s] == vPos && m_vNumIP[s] == numIP && m_vTimePoint[s] == theTimePoint)
					return s;
		}

//	if series not yet registered, add it
	vvIP.push_back(vPos);
	m_vNumIP.push_back(numIP);
	m_vTimePoint.push_back(theTimePoint);
	m_vMayChange.push_back(bMayChange);

//	invoke callback:
//	This callback is called, whenever the local_ip_series have changed. It
//	allows derived classes to react on this changes. For example, the data
//	linker must himself request local_ip_series from the data inputs of
//	the linker. In addition value fields and derivative fields must be adjusted
//	in UserData<TData, dim> etc.
	local_ip_series_added(m_vNumIP.size() - 1);

//	return new series id
	return m_vNumIP.size() - 1;
}


template <int dim>
template <int ldim>
void ICplUserData<dim>::set_local_ips(const size_t seriesID,
                            const MathVector<ldim>* vPos,
                            const size_t numIP)
{
//	check series id
	if(seriesID >= num_series())
		UG_THROW("Trying to set new ips for invalid seriesID "<<seriesID);

//	check that series is changeable
	if(!m_vMayChange[seriesID])
		UG_THROW("Local IP is not changable, but trying to set new ips.");

//	check, that dimension is ok.
	if(m_locPosDim == -1) m_locPosDim = ldim;
	else if(m_locPosDim != ldim)
		UG_THROW("Local IP dimension conflict");

//	get local positions
	std::vector<const MathVector<ldim>*>& vvIP = get_local_ips(Int2Type<ldim>());

//	check if still at same position and with same numIP. In that case the
//	positions have not changed. We have nothing to do
	if(vvIP[seriesID] == vPos && m_vNumIP[seriesID] == numIP) return;

//	remember new positions and numIP
	vvIP[seriesID] = vPos;
	m_vNumIP[seriesID] = numIP;

//	invoke callback:
//	This callback is called, whenever the local_ip_series have changed. It
//	allows derived classes to react on this changes. For example, the data
//	linker must himself request local_ip_series from the data inputs of
//	the linker. In addition value fields and derivative fields must be adjusted
//	in UserData<TData, dim> etc.
	local_ips_changed(seriesID, numIP);
}

template <int dim>
void ICplUserData<dim>::set_time_point(const size_t seriesID,
                            		const int timePointSpec)
{
//	check series id
	if(seriesID >= num_series())
		UG_THROW("Trying to set new ips for invalid seriesID "<<seriesID);

//	check that series is changeable
	if(!m_vMayChange[seriesID])
		UG_THROW("Time point specification is not changable, but trying to set a new one.");
	
//	set the new time point specification (if it is not prescribed by the object)
	m_vTimePoint[seriesID] = (m_defaultTimePoint >= 0)? m_defaultTimePoint : timePointSpec;
	
//TODO: Should we call the callback here? (No data sizes are changed!)
}

template <int dim>
template <int ldim>
const MathVector<ldim>* ICplUserData<dim>::local_ips(size_t s) const
{
//	check, that dimension is ok.
	if(m_locPosDim != ldim) UG_THROW("Local IP dimension conflict");

	UG_ASSERT(s < num_series(), "Wrong series id");

//	NOTE: local ips may be nullptr, if no ip position given, i.e. num_ip(s) == 0
	return get_local_ips(Int2Type<ldim>())[s];
}

template <int dim>
template <int ldim>
const MathVector<ldim>& ICplUserData<dim>::local_ip(size_t s, size_t ip) const
{
//	check, that dimension is ok.
	if(m_locPosDim != ldim) UG_THROW("Local IP dimension conflict");

	UG_ASSERT(s < num_series(), "Wrong series id");
	UG_ASSERT(ip < num_ip(s), "Invalid index.");

	return get_local_ips(Int2Type<ldim>())[s][ip];
}

template <int dim>
int ICplUserData<dim>::time_point_specification(size_t s) const
{
	UG_ASSERT(s < num_series(), "Wrong series id");

	return m_vTimePoint[s];
}

template <int dim>
size_t ICplUserData<dim>::time_point(size_t s) const
{
	UG_ASSERT(s < num_series(), "Wrong series id:" << s << ">=" << num_series());

//	size_t time_spec;
//	if ((time_spec = m_vTimePoint[s]) >= 0)
	if (m_vTimePoint[s] >= 0)
		return m_vTimePoint[s];
	return m_timePoint;
}

template <int dim>
bool ICplUserData<dim>::at_current_time(size_t s) const
{
	UG_ASSERT(s < num_series(), "Wrong series id:" << s << ">=" << num_series());
	
	int time_spec;
	if ((time_spec = m_vTimePoint[s]) >= 0)
		return ((size_t) time_spec) == m_timePoint;
	return true;
}

template <int dim>
void ICplUserData<dim>::set_global_ips(size_t s, const MathVector<dim>* vPos, size_t numIP)
{
	UG_ASSERT(s < num_series(), "Wrong series id: "<<s<<" (numSeries: "<<num_series()<<")");

//	check number of ips (must match local ip number)
	if(numIP != num_ip(s))
		UG_THROW("UserData::set_global_ips: Num Local IPs is " << num_ip(s)
		               << ", but trying to set Num Global IPs: " << numIP <<
		               " for series "<< s);

//	remember global positions
	m_vvGlobPos[s] = vPos;

//	invoke callback:
//	this callback is called every time the global position changes. It gives
//	derived classes the possibility to react on this fact. E.g. the data
//	linker must forward the global positions to its own imports.
	global_ips_changed(s, vPos, numIP);
}

template <int dim>
inline void ICplUserData<dim>::check_s(size_t s) const
{
	UG_ASSERT(s < num_series(), "Wrong series id");
	UG_ASSERT(s < m_vvGlobPos.size(), "Invalid index.");
}

template <int dim>
inline void ICplUserData<dim>::check_s_ip(size_t s, size_t ip) const
{
	check_s(s);
	UG_ASSERT(ip < num_ip(s), "Invalid index.");
	UG_ASSERT(m_vvGlobPos[s] != nullptr, "Global IP not set.");
}

////////////////////////////////////////////////////////////////////////////////
//	UserData
////////////////////////////////////////////////////////////////////////////////

template <typename TData, int dim, typename TRet>
void CplUserData<TData,dim,TRet>::
register_storage_callback(DataImport<TData,dim>* obj, void (DataImport<TData,dim>::*func)())
{
	using Pair = std::pair<DataImport<TData,dim>*, CallbackFct>;
	m_vCallback.push_back(Pair(obj, std::bind(func, obj)));
}

template <typename TData, int dim, typename TRet>
void CplUserData<TData,dim,TRet>::
unregister_storage_callback(DataImport<TData,dim>* obj)
{
	using VecType = std::vector<std::pair<DataImport<TData,dim>*, CallbackFct> >;
	using iterator = typename VecType::iterator;
	iterator iter = m_vCallback.begin();
	while(iter != m_vCallback.end())
	{
		if((*iter).first == obj) iter = m_vCallback.erase(iter);
		else ++iter;
	}
}

template <typename TData, int dim, typename TRet>
void CplUserData<TData,dim,TRet>::
call_storage_callback() const
{
	using VecType = std::vector<std::pair<DataImport<TData,dim>*, CallbackFct> >;
	using iterator = typename VecType::const_iterator;
	for(iterator iter = m_vCallback.begin(); iter != m_vCallback.end(); ++iter)
	{
		//		(((*iter).first)->*((*iter).second))();
		((*iter).second)();
	}
}

template <typename TData, int dim, typename TRet>
inline void CplUserData<TData,dim,TRet>::check_series(size_t s) const
{
	UG_ASSERT(s < num_series(), "Wrong series id"<<s);
	UG_ASSERT(s < m_vvValue.size(), "Invalid index "<<s);
}

template <typename TData, int dim, typename TRet>
inline void CplUserData<TData,dim,TRet>::check_series_ip(size_t s, size_t ip) const
{
	check_series(s);
	UG_ASSERT(ip < num_ip(s), "Invalid index "<<ip);
	UG_ASSERT(ip < m_vvValue[s].size(), "Invalid index "<<ip);
}

template <typename TData, int dim, typename TRet>
void CplUserData<TData,dim,TRet>::local_ip_series_added(const size_t seriesID)
{
	const size_t s = seriesID;

//	check, that only increasing the data, this is important to guarantee,
//	that the allocated memory pointer remain valid. They are used outside of
//	the class as well to allow fast access to the data.
	if(s < m_vvValue.size())
		UG_THROW("Decrease is not implemented. Series: "<<s<<
		         	 	 ", currNumSeries: "<<m_vvValue.size());

//	increase number of series if needed
	m_vvValue.resize(s+1);
	m_vvBoolFlag.resize(s+1);

//	allocate new storage
	m_vvValue[s].resize(num_ip(s));
	m_vvBoolFlag[s].resize(num_ip(s), true);
	value_storage_changed(s);
	call_storage_callback();

//	call base class callback
	base_type::local_ip_series_added(seriesID);
}

template <typename TData, int dim, typename TRet>
void CplUserData<TData,dim,TRet>::local_ip_series_to_be_cleared()
{
//	free the memory
//	clear all series
	m_vvValue.clear();
	m_vvBoolFlag.clear();

//	call base class callback (if implementation given)
//	base_type::local_ip_series_to_be_cleared();
}

template <typename TData, int dim, typename TRet>
void CplUserData<TData,dim,TRet>::local_ips_changed(const size_t seriesID, const size_t newNumIP)
{
//	resize only when more data is needed than actually allocated
	if(newNumIP >= m_vvValue[seriesID].size())
	{
	//	resize
		m_vvValue[seriesID].resize(newNumIP);
		m_vvBoolFlag[seriesID].resize(newNumIP, true);

	//	invoke callback
		value_storage_changed(seriesID);
		call_storage_callback();
	}

//	call base class callback (if implementation given)
//	base_type::local_ips_changed(seriesID);
}

////////////////////////////////////////////////////////////////////////////////
//	DependentUserData
////////////////////////////////////////////////////////////////////////////////

template <typename TData, int dim>
void DependentUserData<TData,dim>::set_function_pattern(ConstSmartPtr<FunctionPattern> fctPatt)
{
	this->m_fctGrp.set_function_pattern(fctPatt);
	extract_fct_grp();
}

template <typename TData, int dim>
void DependentUserData<TData,dim>::set_functions(const char* symbFct)
{
	set_functions(std::string(symbFct));
}

template <typename TData, int dim>
void DependentUserData<TData,dim>::set_functions(const std::string& symbFct)
{
	set_functions(TokenizeTrimString(symbFct));
}

template <typename TData, int dim>
void DependentUserData<TData,dim>::set_functions(const std::vector<std::string>& symbFct)
{
	m_SymbFct = symbFct;
	extract_fct_grp();
}

template <typename TData, int dim>
void DependentUserData<TData,dim>::extract_fct_grp()
{
	//	if associated infos missing return
	ConstSmartPtr<FunctionPattern> spFctPatt = this->m_fctGrp.function_pattern();
	if(spFctPatt.invalid()) return;

	//	if no function passed, clear functions
	if(m_SymbFct.size() == 1 && m_SymbFct[0].empty()) m_SymbFct.clear();

	//	if functions passed with separator, but not all tokens filled, throw error
	for(size_t i = 0; i < m_SymbFct.size(); ++i)
	{
		if(m_SymbFct.empty())
			UG_THROW("Error while setting functions in a DependentUserData: passed "
					"function string lacks a "
					"function specification at position "<<i<<"(of "
					<<m_SymbFct.size()-1<<")");
	}

	if(m_SymbFct.empty()){
		this->m_fctGrp.clear();
		return;
	}

	//	create function group of this elem disc
	try{
		this->m_fctGrp.clear();
		this->m_fctGrp.add(m_SymbFct);
	}UG_CATCH_THROW("DependentUserData: Cannot find some symbolic function "
			"name.");

	//	create a mapping between all functions and the function group of this
	//	element disc.
	try{
		CreateFunctionIndexMapping(this->m_map, this->m_fctGrp, spFctPatt);
	}UG_CATCH_THROW("DependentUserData: Cannot create Function Index Mapping.");

	this->check_setup();
}

template <typename TData, int dim>
void DependentUserData<TData,dim>::update_dof_sizes(const LocalIndices& ind)
{
//	check size
	const FunctionIndexMapping& map = this->map();
	UG_ASSERT(map.num_fct() == this->num_fct(), "Number function mismatch.");

//	cache numFct and their numDoFs
	m_vvNumDoFPerFct.resize(map.num_fct());
	for(size_t fct = 0; fct < m_vvNumDoFPerFct.size(); ++fct)
		m_vvNumDoFPerFct[fct] = ind.num_dof(map[fct]);

	resize_deriv_array();
}

template <typename TData, int dim>
void DependentUserData<TData,dim>::resize_deriv_array()
{
//	resize num fct
	for(size_t s = 0; s < m_vvvvDeriv.size(); ++s)
		resize_deriv_array(s);
}

template <typename TData, int dim>
void DependentUserData<TData,dim>::resize_deriv_array(const size_t s)
{
//	resize ips
	m_vvvvDeriv[s].resize(num_ip(s));

	for(size_t ip = 0; ip < m_vvvvDeriv[s].size(); ++ip)
	{
	//	resize num fct
		m_vvvvDeriv[s][ip].resize(m_vvNumDoFPerFct.size());

	//	resize dofs
		for(size_t fct = 0; fct < m_vvNumDoFPerFct.size(); ++fct)
			m_vvvvDeriv[s][ip][fct].resize(m_vvNumDoFPerFct[fct]);
	}
}

template <typename TData, int dim>
void DependentUserData<TData,dim>::set_zero(std::vector<std::vector<TData> > vvvDeriv[], const size_t nip)
{
	for(size_t ip = 0; ip < nip; ++ip)
		for(size_t fct = 0; fct <  vvvDeriv[ip].size(); ++fct)
			for(size_t sh = 0; sh <  vvvDeriv[ip][fct].size(); ++sh)
			{
				vvvDeriv[ip][fct][sh] = 0.0;
			}
}

template <typename TData, int dim>
inline void DependentUserData<TData,dim>::check_s_ip(size_t s, size_t ip) const
{
	UG_ASSERT(s < this->num_series(), "Wrong series id"<<s);
	UG_ASSERT(s < m_vvvvDeriv.size(), "Invalid index "<<s);
	UG_ASSERT(ip < this->num_ip(s), "Invalid index "<<ip);
	UG_ASSERT(ip < m_vvvvDeriv[s].size(), "Invalid index "<<ip);
}

template <typename TData, int dim>
inline void DependentUserData<TData,dim>::check_s_ip_fct(size_t s, size_t ip, size_t fct) const
{
	check_s_ip(s,ip);
	UG_ASSERT(fct < m_vvvvDeriv[s][ip].size(), "Invalid index.");
}

template <typename TData, int dim>
inline void DependentUserData<TData,dim>::check_s_ip_fct_dof(size_t s, size_t ip, size_t fct, size_t dof) const
{
	check_s_ip_fct(s,ip,fct);
	UG_ASSERT(dof < m_vvvvDeriv[s][ip][fct].size(), "Invalid index.");
}

template <typename TData, int dim>
void DependentUserData<TData,dim>::local_ip_series_added(const size_t seriesID)
{
//	adjust data arrays
	m_vvvvDeriv.resize(seriesID+1);

//	forward change signal to base class
	base_type::local_ip_series_added(seriesID);
}

template <typename TData, int dim>
void DependentUserData<TData,dim>::local_ip_series_to_be_cleared()
{
//	adjust data arrays
	m_vvvvDeriv.clear();

//	forward change signal to base class
	base_type::local_ip_series_to_be_cleared();
}

template <typename TData, int dim>
void DependentUserData<TData,dim>::local_ips_changed(const size_t seriesID, const size_t newNumIP)
{
	UG_ASSERT(seriesID < m_vvvvDeriv.size(), "wrong series id.");

//	resize only when more data is needed than actually allocated
	if(newNumIP >= m_vvvvDeriv[seriesID].size())
		resize_deriv_array(seriesID);

//	call base class callback (if implementation given)
	base_type::local_ips_changed(seriesID, newNumIP);
}

} // end namespace ug

#endif