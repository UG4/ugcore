/*
 * data_import_export_impl.h
 *
 *  Created on: 04.07.2011
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DATA_IMPORT_EXPORT_IMPL__
#define __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DATA_IMPORT_EXPORT_IMPL__

#include "data_import_export.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
// IDataImport
////////////////////////////////////////////////////////////////////////////////

inline bool IDataImport::set_geometric_object_type(int id)
{
//	if lin defect is not supposed to be computed, we're done
	if(!m_bCompLinDefect) return true;

//	Check for evaluation function and choose it if present
	if(id < (int)m_vLinDefectFunc.size() && m_vLinDefectFunc[id] != NULL)
	{
		m_id = id;
		return true;
	}
//	return error else
	else
	{
		UG_LOG("No or not all lin defect functions registered "
				"for object with reference object id " << id << ".\n");
		m_id = -1; return false;
	}
}

template <typename TFunc>
void IDataImport::reg_lin_defect_fct(int id, IElemDisc* obj, TFunc func)
{
//	make sure that there is enough space
	if((size_t)id >= m_vLinDefectFunc.size())
		m_vLinDefectFunc.resize(id+1, NULL);

	m_vLinDefectFunc[id] = (LinDefectFunc)func;
	m_pObj = obj;
}

////////////////////////////////////////////////////////////////////////////////
// DataImport
////////////////////////////////////////////////////////////////////////////////

template <typename TData, int dim>
void DataImport<TData,dim>::set_data(IPData<TData, dim>& data)
{
//	remember IPData
	m_pIPData = &data;

//	remember iexport
	this->m_pIDependentIPData = dynamic_cast<IDependentIPData*>(&data);

//	remember dependent data (i.e. is NULL iff no dependent data given)
	m_pDependentIPData = dynamic_cast<DependentIPData<TData, dim>*>(&data);
}

template <typename TData, int dim>
template <int ldim>
void DataImport<TData,dim>::set_local_ips(const MathVector<ldim>* vPos, size_t numIP)
{
//	if no data set, skip
	if(m_pIPData == NULL) return;

//	request series
	m_seriesID = m_pIPData->template
				register_local_ip_series<ldim>(vPos,numIP);

//	cache the pointer to the data field. This is possible, since once a
//	local ip series is registered it can not be removed ot altered. In the same
//	way the memory starge is not changed but always only increased. Therefore,
//	we can request the data now and it will remain valid until IIPData::clear()
//	is called.
	m_vValue = m_pIPData->values(m_seriesID);

//	in addition we cache the number of ips
	m_numIP = m_pIPData->num_ip(m_seriesID);

//	check that num ip is correct
	UG_ASSERT(m_numIP == numIP, "Different number of ips than requested.");
}

template <typename TData, int dim>
void DataImport<TData,dim>::set_global_ips(const MathVector<dim>* vPos, size_t numIP)
{
//  if no data set, skip
	if(m_pIPData == NULL) return;

//	set global ips for series ID
	UG_ASSERT(m_seriesID >= 0, "Wrong series id.");
	m_pIPData->set_global_ips(m_seriesID,vPos,numIP);
}

template <typename TData, int dim>
void DataImport<TData,dim>::clear_lin_defect()
{
	for(size_t ip = 0; ip < m_vvvLinDefect.size(); ++ip)
		for(size_t fct = 0; fct < m_vvvLinDefect[ip].size(); ++fct)
			for(size_t sh = 0; sh < m_vvvLinDefect[ip][fct].size(); ++sh)
				m_vvvLinDefect[ip][fct][sh] = 0.0;
}

template <typename TData, int dim>
void DataImport<TData,dim>::assemble_jacobian(local_matrix_type& J)
{
	UG_ASSERT(m_pDependentIPData != NULL, "No Export set.");

//	loop integration points
	for(size_t ip = 0; ip < num_ip(); ++ip)
	{
//	loop all functions
	for(size_t fct1 = 0; fct1 < num_fct(); ++fct1)
		for(size_t fct2 = 0; fct2 < m_pDependentIPData->num_fct(); ++fct2)
		{
//	get array of linearized defect and derivative
	const TData* LinDef = lin_defect(ip, fct1);
	const TData* Deriv = m_pDependentIPData->deriv(m_seriesID, ip, fct2);

//	loop shapes of functions
	for(size_t sh1 = 0; sh1 < num_sh(fct1); ++sh1)
		for(size_t sh2 = 0; sh2 < m_pDependentIPData->num_sh(m_seriesID, fct2); ++sh2)
		{
			J(fct1, sh1, fct2, sh2) += LinDef[sh1]*Deriv[sh2];
		}
		}
	}
}

template <typename TData, int dim>
void DataImport<TData,dim>::resize(const LocalIndices& ind, const FunctionIndexMapping& map)
{
//	resize ips
	//\todo: Move this call to some place, where num_ip is changed.
	m_vvvLinDefect.resize(num_ip());

//	resize num fct
	for(size_t ip = 0; ip < num_ip(); ++ip)
	{
	//	resize num fct
		m_vvvLinDefect[ip].resize(map.num_fct());

	//	resize dofs
		for(size_t fct = 0; fct < map.num_fct(); ++fct)
			m_vvvLinDefect[ip][fct].resize(ind.num_dof(map[fct]));
	}
}

template <typename TData, int dim>
inline void DataImport<TData,dim>::check_ip_fct(size_t ip, size_t fct) const
{
	check_ip(ip);
	UG_ASSERT(ip  < m_vvvLinDefect.size(), "Invalid index.");
	UG_ASSERT(fct < m_vvvLinDefect[ip].size(), "Invalid index.");
}

template <typename TData, int dim>
inline void DataImport<TData,dim>::check_ip_fct_sh(size_t ip, size_t fct, size_t sh) const
{
	check_ip_fct(ip, fct);
	UG_ASSERT(sh < m_vvvLinDefect[ip][fct].size(), "Invalid index.");
}

template <typename TData, int dim>
inline void DataImport<TData,dim>::check_ip(size_t ip) const
{
	UG_ASSERT(ip < m_numIP, "Invalid index.");
}

template <typename TData, int dim>
inline void DataImport<TData,dim>::check_values() const
{
	UG_ASSERT(m_vValue != NULL, "Data Value field not set.");
}

////////////////////////////////////////////////////////////////////////////////
// IDataExport
////////////////////////////////////////////////////////////////////////////////

inline bool IDataExport::set_geometric_object_type(int id)
{
	if(id < 0 || (size_t)id >= m_vExportFunc.size() || m_vExportFunc[id] == NULL)
	{
		UG_LOG("No or not all functions registered "
				"for object with reference object id " << id << ".\n");
		m_id = -1; return false;
	}
	else{m_id = id;	return true;}
}

template <typename TFunc>
void IDataExport::reg_export_fct(int id, IElemDisc* obj, TFunc func)
{
//	make sure that there is enough space
	if((size_t)id >= m_vExportFunc.size())
		m_vExportFunc.resize(id+1, NULL);

	m_vExportFunc[id] = (ExportFunc)func;

	if(m_pObj == NULL) m_pObj = obj;
	else if(m_pObj != obj)
		throw(UGFatalError("Exports assume to be used by on object for all functions."));
}

inline bool IDataExport::is_ready() const
{
	if(m_id < 0) {
		UG_LOG("ERROR in 'IDataExport::is_ready':"
				"ElemType id is not set correctly."); return false;
	}
	if((size_t)m_id >= m_vExportFunc.size()) {
		UG_LOG("ERROR in 'IDataExport::is_ready': There is no evaluation "
				"function registered for export and elem type "<<m_id<<".\n");
		return false;
	}
	if(m_vExportFunc[m_id] == NULL) {
		UG_LOG("ERROR in 'IDataExport::is_ready': There is a evaluation "
				"function registered for export and elem type "<<m_id<<
				", but function pointer is zero. \n");
		return false;
	}

//	everything is ok
	return true;
}

} // end namespace ug

#endif /* __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DATA_IMPORT_EXPORT_IMPL__ */
