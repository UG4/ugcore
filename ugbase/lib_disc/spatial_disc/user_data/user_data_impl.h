/*
 * user_data_impl.h
 *
 *  Created on: 04.07.2011
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__USER_DATA__USER_DATA_IMPL__
#define __H__UG__LIB_DISC__SPATIAL_DISC__USER_DATA__USER_DATA_IMPL__

#include "user_data.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
//	IUserData
////////////////////////////////////////////////////////////////////////////////

inline IUserData::IUserData() : m_locPosDim(-1), m_time(0.0), m_si(-1)
{
	m_vNumIP.clear();
	m_vMayChange.clear();
	m_locPosDim = -1;
	m_pvLocIP1d.clear(); m_pvLocIP2d.clear(); m_pvLocIP3d.clear();
}

inline void IUserData::clear()
{
	local_ip_series_to_be_cleared();
	m_vNumIP.clear();
	m_vMayChange.clear();
	m_locPosDim = -1;
	m_pvLocIP1d.clear(); m_pvLocIP2d.clear(); m_pvLocIP3d.clear();
	m_time = 0.0;
	m_si = -1;
}

template <int ldim>
size_t IUserData::register_local_ip_series(const MathVector<ldim>* vPos,
                                         const size_t numIP,
                                         bool bMayChange)
{
//	check, that dimension is ok.
	if(m_locPosDim == -1) m_locPosDim = ldim;
	else if(m_locPosDim != ldim)
		UG_THROW("Local IP dimension conflict");

//	get local positions
	std::vector<const MathVector<ldim>*>& vvIP = get_local_ips(Int2Type<ldim>());

//	search for ips
//	we only identify ip series if the local ip positions will not change
	if(!bMayChange && numIP != 0)
		for(size_t s = 0; s < vvIP.size(); ++s)
		{
		//	return series number iff exists and local ips remain constant
			if(!m_vMayChange[s])
				if(vvIP[s] == vPos && m_vNumIP[s] == numIP) return s;
		}

//	if series not yet registered, add it
	vvIP.push_back(vPos);
	m_vNumIP.push_back(numIP);
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


template <int ldim>
void IUserData::set_local_ips(const size_t seriesID,
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

template <int ldim>
const MathVector<ldim>* IUserData::local_ips(size_t s) const
{
//	check, that dimension is ok.
	if(m_locPosDim != ldim) UG_THROW("Local IP dimension conflict");

	UG_ASSERT(s < num_series(), "Wrong series id");

//	NOTE: local ips may be NULL, if no ip position given, i.e. num_ip(s) == 0
	return get_local_ips(Int2Type<ldim>())[s];
}

template <int ldim>
const MathVector<ldim>& IUserData::local_ip(size_t s, size_t ip) const
{
//	check, that dimension is ok.
	if(m_locPosDim != ldim) UG_THROW("Local IP dimension conflict");

	UG_ASSERT(s < num_series(), "Wrong series id");
	UG_ASSERT(ip < num_ip(s), "Invalid index.");

	return get_local_ips(Int2Type<ldim>())[s][ip];
}

////////////////////////////////////////////////////////////////////////////////
//	IIPDimData
////////////////////////////////////////////////////////////////////////////////

template <int dim>
void IDimUserData<dim>::set_global_ips(size_t s, const MathVector<dim>* vPos, size_t numIP)
{
	numIP = numIP +1;
	numIP = numIP -1;

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
inline void IDimUserData<dim>::check_s(size_t s) const
{
	UG_ASSERT(s < num_series(), "Wrong series id");
	UG_ASSERT(s < m_vvGlobPos.size(), "Invalid index.");
}

template <int dim>
inline void IDimUserData<dim>::check_s_ip(size_t s, size_t ip) const
{
	check_s(s);
	UG_ASSERT(ip < num_ip(s), "Invalid index.");
	UG_ASSERT(m_vvGlobPos[s] != NULL, "Local IP not set.");
}

////////////////////////////////////////////////////////////////////////////////
//	UserData
////////////////////////////////////////////////////////////////////////////////

template <typename TData, int dim, typename TRet>
TRet UserData<TData,dim,TRet>::
operator() (TData& value,
            const MathVector<dim>& globIP,
            number time, int si) const
{
	UG_THROW("UserData: operator()(TData, MathVector<dim>, "
			"time, si) not implemented.");
}

template <typename TData, int dim, typename TRet>
void UserData<TData,dim,TRet>::
operator() (TData vValue[],
            const MathVector<dim> vGlobIP[],
            number time, int si, const size_t nip) const
{
	UG_THROW("UserData: operator()(TData[], MathVector<dim>[], "
			"time, si) not implemented.");
}

template <typename TData, int dim, typename TRet>
TRet UserData<TData,dim,TRet>::
operator() (TData& value,
            const MathVector<dim>& globIP,
            number time, int si,
            LocalVector& u,
            GeometricObject* elem,
            const MathVector<dim> vCornerCoords[],
            const MathVector<1>& locIP) const
{
	UG_THROW("UserData: operator()(TData, MathVector<dim>, "
			"time, si, LocalVector, GeometricObject*, MathVector<1>) not implemented.");
}

template <typename TData, int dim, typename TRet>
void UserData<TData,dim,TRet>::
operator()(TData vValue[],
           const MathVector<dim> vGlobIP[],
           number time, int si,
           LocalVector& u,
           GeometricObject* elem,
           const MathVector<dim> vCornerCoords[],
           const MathVector<1> vLocIP[],
           const size_t nip,
           const MathMatrix<1, dim>* vJT) const
{
	for(size_t ip = 0; ip < nip; ++ip)
	{
		operator()(vValue[ip], vGlobIP[ip], time, si,
		                  u, elem, vCornerCoords, vLocIP[ip]);
	}
}

template <typename TData, int dim, typename TRet>
TRet UserData<TData,dim,TRet>::
operator() (TData& value,
            const MathVector<dim>& globIP,
            number time, int si,
            LocalVector& u,
            GeometricObject* elem,
            const MathVector<dim> vCornerCoords[],
            const MathVector<2>& locIP) const
{
	UG_THROW("UserData: operator()(TData, MathVector<dim>, "
			"time, si, LocalVector, GeometricObject*, MathVector<2>) not implemented.");
}

template <typename TData, int dim, typename TRet>
void UserData<TData,dim,TRet>::
operator()(TData vValue[],
           const MathVector<dim> vGlobIP[],
           number time, int si,
           LocalVector& u,
           GeometricObject* elem,
           const MathVector<dim> vCornerCoords[],
           const MathVector<2> vLocIP[],
           const size_t nip,
           const MathMatrix<2, dim>* vJT) const
{
	for(size_t ip = 0; ip < nip; ++ip)
	{
		operator()(vValue[ip], vGlobIP[ip], time, si,
		                  u, elem, vCornerCoords, vLocIP[ip]);
	}
}

template <typename TData, int dim, typename TRet>
TRet UserData<TData,dim,TRet>::
operator() (TData& value,
            const MathVector<dim>& globIP,
            number time, int si,
            LocalVector& u,
            GeometricObject* elem,
            const MathVector<dim> vCornerCoords[],
            const MathVector<3>& locIP) const
{
	UG_THROW("UserData: operator()(TData, MathVector<dim>, "
			"time, si, LocalVector, GeometricObject*, MathVector<3>) not implemented.");
}

template <typename TData, int dim, typename TRet>
void UserData<TData,dim,TRet>::
operator()(TData vValue[],
           const MathVector<dim> vGlobIP[],
           number time, int si,
           LocalVector& u,
           GeometricObject* elem,
           const MathVector<dim> vCornerCoords[],
           const MathVector<3> vLocIP[],
           const size_t nip,
           const MathMatrix<3, dim>* vJT) const
{
	for(size_t ip = 0; ip < nip; ++ip)
	{
		operator()(vValue[ip], vGlobIP[ip], time, si,
		                  u, elem, vCornerCoords, vLocIP[ip]);
	}
}

template <typename TData, int dim, typename TRet>
void UserData<TData,dim,TRet>::
register_storage_callback(DataImport<TData,dim>* obj, void (DataImport<TData,dim>::*func)())
{
	typedef std::pair<DataImport<TData,dim>*, CallbackFct> Pair;
	m_vCallback.push_back(Pair(obj,func));
}

template <typename TData, int dim, typename TRet>
void UserData<TData,dim,TRet>::
unregister_storage_callback(DataImport<TData,dim>* obj)
{
	typedef typename std::vector<std::pair<DataImport<TData,dim>*, CallbackFct> > VecType;
	typedef typename VecType::iterator iterator;
	iterator iter = m_vCallback.begin();
	while(iter != m_vCallback.end())
	{
		if((*iter).first == obj) iter = m_vCallback.erase(iter);
		else ++iter;
	}
}

template <typename TData, int dim, typename TRet>
void UserData<TData,dim,TRet>::
call_storage_callback() const
{
	typedef typename std::vector<std::pair<DataImport<TData,dim>*, CallbackFct> > VecType;
	typedef typename VecType::const_iterator iterator;
	for(iterator iter = m_vCallback.begin(); iter != m_vCallback.end(); ++iter)
	{
		(((*iter).first)->*((*iter).second))();
	}
}

template <typename TData, int dim, typename TRet>
inline void UserData<TData,dim,TRet>::check_series(size_t s) const
{
	UG_ASSERT(s < num_series(), "Wrong series id"<<s);
	UG_ASSERT(s < m_vvValue.size(), "Invalid index "<<s);
}

template <typename TData, int dim, typename TRet>
inline void UserData<TData,dim,TRet>::check_series_ip(size_t s, size_t ip) const
{
	check_series(s);
	UG_ASSERT(ip < num_ip(s), "Invalid index "<<ip);
	UG_ASSERT(ip < m_vvValue[s].size(), "Invalid index "<<ip);
}

template <typename TData, int dim, typename TRet>
void UserData<TData,dim,TRet>::local_ip_series_added(const size_t seriesID)
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
void UserData<TData,dim,TRet>::local_ip_series_to_be_cleared()
{
//	free the memory
//	clear all series
	m_vvValue.clear();
	m_vvBoolFlag.clear();

//	call base class callback (if implementation given)
//	base_type::local_ip_series_to_be_cleared();
}

template <typename TData, int dim, typename TRet>
void UserData<TData,dim,TRet>::local_ips_changed(const size_t seriesID, const size_t newNumIP)
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
void DependentUserData<TData,dim>::set_dof_sizes(const LocalIndices& ind,
                                               const FunctionIndexMapping& map)
{
//	check size
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
void DependentUserData<TData,dim>::clear_derivative_values()
{
	for(size_t s = 0; s < m_vvvvDeriv.size(); ++s)
		for(size_t ip = 0; ip < m_vvvvDeriv[s].size(); ++ip)
			for(size_t fct = 0; fct <  m_vvvvDeriv[s][ip].size(); ++fct)
				for(size_t sh = 0; sh <  m_vvvvDeriv[s][ip][fct].size(); ++sh)
				{
					m_vvvvDeriv[s][ip][fct][sh] = 0.0;
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
	UserData<TData, dim>::local_ip_series_added(seriesID);
}

template <typename TData, int dim>
void DependentUserData<TData,dim>::local_ip_series_to_be_cleared()
{
//	adjust data arrays
	m_vvvvDeriv.clear();

//	forward change signal to base class
	UserData<TData, dim>::local_ip_series_to_be_cleared();
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

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__USER_DATA__USER_DATA_IMPL__ */
