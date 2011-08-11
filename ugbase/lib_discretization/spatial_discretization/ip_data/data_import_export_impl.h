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

inline bool IDataImport::set_geometric_object_type(ReferenceObjectID id)
{
//	if lin defect is not supposed to be computed, we're done
	if(!m_bCompLinDefect) return true;

//	Check for evaluation function and choose it if present
	if(m_vLinDefectFunc[id] != NULL)
	{
		m_id = id;
		return true;
	}

//	return error else
	else
	{
		UG_LOG("ERROR in 'IDataImport::set_geometric_object_type':"
				"No lin defect functions registered "
				"for object with reference object id " << id << ".\n");
		m_id = ROID_INVALID;
		return false;
	}
}

template <typename TFunc>
void IDataImport::reg_lin_defect_fct(ReferenceObjectID id, IElemDisc* obj, TFunc func)
{
	m_vLinDefectFunc[id] = static_cast<LinDefectFunc>(func);
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
// DataExport
////////////////////////////////////////////////////////////////////////////////

template <typename TData, int dim>
DataExport<TData, dim>::DataExport() : m_pObj(NULL)
{
	this->m_bCompNeedsSol = true;

	clear_export_fct();
}

template <typename TData, int dim>
void DataExport<TData, dim>::clear_export_fct()
{
	for(size_t roid = 0; roid < NUM_REFERENCE_OBJECTS; ++roid)
		m_vExportFunc[roid] = NULL;
}

template <typename TData, int dim>
template <typename TFunc>
void DataExport<TData, dim>::reg_export_fct(ReferenceObjectID id,
                                            IElemDisc* obj, TFunc func)
{
	m_vExportFunc[id] = static_cast<ExportFunc>(func);

	if(m_pObj == NULL) m_pObj = obj;
	else if(m_pObj != obj)
		throw(UGFatalError("Exports assume to be used by on object for all functions."));
}

template <typename TData, int dim>
bool DataExport<TData, dim>::set_geometric_object_type(ReferenceObjectID id)
{
	if(m_vExportFunc[m_id] == NULL) {
		UG_LOG("ERROR in 'DataExport::is_ready': There is no evaluation "
				"function registered for export and elem type "<<m_id<<".\n");
		return false;
	}

	m_id = id;
	return true;
}

template <typename TData, int dim>
bool DataExport<TData, dim>::is_ready() const
{
//	check base
	if(m_vExportFunc[m_id] == NULL) {
		UG_LOG("ERROR in 'DataExport::is_ready': There is no evaluation "
				"function registered for export and elem type "<<m_id<<".\n");
		return false;
	}

//	everything is ok
	return true;
}

template <typename TData, int dim>
bool DataExport<TData, dim>::compute(bool bDeriv)
{
	UG_LOG("ERROR in 'DataExport::compute()': Computation of Export "
		 	 "without current solution called. Cannot evaluate.\n");
	return false;
}

template <typename TData, int dim>
bool DataExport<TData, dim>::compute(const local_vector_type& u, bool bDeriv)
{
	UG_ASSERT(m_vExportFunc[m_id] != NULL, "Func pointer is NULL");
	return (m_pObj->*(m_vExportFunc[m_id]))(u, bDeriv);
}


} // end namespace ug

#endif /* __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DATA_IMPORT_EXPORT_IMPL__ */
