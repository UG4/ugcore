/*
 * ip_data_impl.h
 *
 *  Created on: 04.07.2011
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__IP_DATA__IP_DATA_IMPL__
#define __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__IP_DATA__IP_DATA_IMPL__

#include "ip_data.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
//	IIPData
////////////////////////////////////////////////////////////////////////////////

inline IIPData::IIPData() : m_locPosDim(-1), m_time(0.0)
{
	m_vNumIP.clear();
	m_pvLocIP1d.clear(); m_pvLocIP2d.clear(); m_pvLocIP3d.clear();
}

inline void IIPData::clear_ips()
{
	m_vNumIP.clear();
	m_locPosDim = -1;
	m_pvLocIP1d.clear(); m_pvLocIP2d.clear(); m_pvLocIP3d.clear();
	adjust_global_ips_and_data(m_vNumIP);
}

template <int ldim>
size_t IIPData::register_local_ip_series(const MathVector<ldim>* vPos, size_t numIP)
{
//	check, that dimension is ok.
	if(m_locPosDim == -1) m_locPosDim = ldim;
	else if(m_locPosDim != ldim)
		throw(UGFatalError("Local IP dimension conflict"));

//	get local positions
	std::vector<const MathVector<ldim>*>& vvIP = get_local_ips(Int2Type<ldim>());

//	search for ips
	for(size_t s = 0; s < vvIP.size(); ++s)
	{
	//	return series number iff exists
		if(vvIP[s] == vPos && m_vNumIP[s] == numIP) return s;
	}

//	if series not yet registered, add it
	vvIP.push_back(vPos);
	m_vNumIP.push_back(numIP);

//	resize global ips and data
	adjust_global_ips_and_data(m_vNumIP);

//	return new series id
	return m_vNumIP.size() - 1;
}

template <int ldim>
const MathVector<ldim>* IIPData::local_ips(size_t s) const
{
//	check, that dimension is ok.
	if(m_locPosDim != ldim) throw(UGFatalError("Local IP dimension conflict"));

	UG_ASSERT(s < num_series(), "Wrong series id");
	UG_ASSERT(get_local_ips(Int2Type<ldim>())[s] != NULL, "Zero local ip pointer.");

	return get_local_ips(Int2Type<ldim>())[s];
}

template <int ldim>
const MathVector<ldim>& IIPData::local_ip(size_t s, size_t ip) const
{
//	check, that dimension is ok.
	if(m_locPosDim != ldim) throw(UGFatalError("Local IP dimension conflict"));

	UG_ASSERT(s < num_series(), "Wrong series id");
	UG_ASSERT(ip < num_ip(s), "Invalid index.");

	return get_local_ips(Int2Type<ldim>())[s][ip];
}

////////////////////////////////////////////////////////////////////////////////
//	IIPDimData
////////////////////////////////////////////////////////////////////////////////

template <int dim>
void IIPDimData<dim>::set_global_ips(size_t s, const MathVector<dim>* vPos, size_t numIP)
{
	UG_ASSERT(s < num_series(), "Wrong series id");

//	check number of ips (must match local ip number)
	if(numIP != num_ip(s))
	{
		UG_LOG("ERROR in 'IPData::set_global_ips':"
				" Num Local IPs is " << num_ip(s)  << ", but trying to set"
				" Num Global IPs: " << numIP << " for series "<< s<< ".\n");
		throw(UGFatalError("Num ip does not match."));
	}

//	remember global positions
	m_vvGlobPos[s] = vPos;

//	invoke callback
	global_ips_changed(s, vPos, numIP);
}

template <int dim>
inline void IIPDimData<dim>::check_s(size_t s) const
{
	UG_ASSERT(s < num_series(), "Wrong series id");
	UG_ASSERT(s < m_vvGlobPos.size(), "Invalid index.");
}

template <int dim>
inline void IIPDimData<dim>::check_s_ip(size_t s, size_t ip) const
{
	check_s(s);
	UG_ASSERT(ip < num_ip(s), "Invalid index.");
	UG_ASSERT(m_vvGlobPos[s] != NULL, "Local IP not set.");
}

////////////////////////////////////////////////////////////////////////////////
//	IPData
////////////////////////////////////////////////////////////////////////////////

template <typename TData, int dim>
inline void IPData<TData,dim>::check_series(size_t s) const
{
	UG_ASSERT(s < num_series(), "Wrong series id"<<s);
	UG_ASSERT(s < m_vvValue.size(), "Invalid index "<<s);
}

template <typename TData, int dim>
inline void IPData<TData,dim>::check_series_ip(size_t s, size_t ip) const
{
	check_series(s);
	UG_ASSERT(ip < num_ip(s), "Invalid index "<<ip);
	UG_ASSERT(ip < m_vvValue[s].size(), "Invalid index "<<ip);
}

template <typename TData, int dim>
void IPData<TData,dim>::adjust_global_ips_and_data(const std::vector<size_t>& vNumIP)
{
//	adjust data arrays
	m_vvValue.resize(vNumIP.size());
	for(size_t s = 0; s < vNumIP.size(); ++s)
		m_vvValue[s].resize(vNumIP[s]);

	base_type::adjust_global_ips_and_data(vNumIP);
}

////////////////////////////////////////////////////////////////////////////////
//	DependentIPData
////////////////////////////////////////////////////////////////////////////////

template <typename TData, int dim>
void DependentIPData<TData,dim>::resize(const LocalIndices& ind,
                                        const FunctionIndexMapping& map)
{
//	resize num fct
	for(size_t s = 0; s < m_vvvDeriv.size(); ++s)
		for(size_t ip = 0; ip < m_vvvDeriv[s].size(); ++ip)
		{
		//	number of functions
			const size_t numFct = map.num_fct();

		//	resize num fct
			m_vvvDeriv[s][ip].resize(numFct);

		//	resize dofs
			for(size_t fct = 0; fct < numFct; ++fct)
				m_vvvDeriv[s][ip][fct].resize(ind.num_dof(map[fct]));
		}
}

template <typename TData, int dim>
void DependentIPData<TData,dim>::clear_derivative_values()
{
	for(size_t s = 0; s < m_vvvDeriv.size(); ++s)
		for(size_t ip = 0; ip < m_vvvDeriv[s].size(); ++ip)
			for(size_t fct = 0; fct <  m_vvvDeriv[s][ip].size(); ++fct)
				for(size_t sh = 0; sh <  m_vvvDeriv[s][ip][fct].size(); ++sh)
				{
					m_vvvDeriv[s][ip][fct][sh] = 0.0;
				}

}

template <typename TData, int dim>
inline void DependentIPData<TData,dim>::check_s_ip(size_t s, size_t ip) const
{
	UG_ASSERT(s < num_series(), "Wrong series id"<<s);
	UG_ASSERT(s < m_vvvDeriv.size(), "Invalid index "<<s);
	UG_ASSERT(ip < num_ip(s), "Invalid index "<<ip);
	UG_ASSERT(ip < m_vvvDeriv[s].size(), "Invalid index "<<ip);
}

template <typename TData, int dim>
inline void DependentIPData<TData,dim>::check_s_ip_fct(size_t s, size_t ip, size_t fct) const
{
	check_s_ip(s,ip);
	UG_ASSERT(fct < m_vvvDeriv[s][ip].size(), "Invalid index.");
}

template <typename TData, int dim>
inline void DependentIPData<TData,dim>::check_s_ip_fct_dof(size_t s, size_t ip, size_t fct, size_t dof) const
{
	check_s_ip_fct(s,ip,fct);
	UG_ASSERT(dof < m_vvvDeriv[s][ip][fct].size(), "Invalid index.");
}

template <typename TData, int dim>
void DependentIPData<TData,dim>::adjust_global_ips_and_data(const std::vector<size_t>& vNumIP)
{
//	adjust data arrays
	m_vvvDeriv.resize(vNumIP.size());
	for(size_t s = 0; s < vNumIP.size(); ++s)
		m_vvvDeriv[s].resize(vNumIP[s]);

//	resize values
	IPData<TData, dim>::adjust_global_ips_and_data(vNumIP);
}


} // end namespace ug

#endif /* __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__IP_DATA__IP_DATA_IMPL__ */
