/*
 * data_import_export_impl.h
 *
 *  Created on: 04.07.2011
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__DATA_IMPORT_IMPL__
#define __H__UG__LIB_DISC__SPATIAL_DISC__DATA_IMPORT_IMPL__

#include "data_import.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
// DataImport
////////////////////////////////////////////////////////////////////////////////

template <typename TData, int dim>
DataImport<TData,dim>::~DataImport()
{
	if(data_given()) m_spUserData->unregister_storage_callback(this);
}

template <typename TData, int dim>
void DataImport<TData,dim>::set_roid(ReferenceObjectID id)
{
//	if lin defect is not supposed to be computed, we're done
	if(!this->m_bCompLinDefect) return;

//	Check for evaluation function and choose it if present
	if(m_vLinDefectFunc[id] != NULL)
	{
		m_id = id;
		return;
	}

	UG_THROW("DataImport::set_roid: No lin defect functions registered for "<<id);
}

template <typename TData, int dim>
template <typename TClass>
void
DataImport<TData,dim>::
set_fct(ReferenceObjectID id, TClass* obj,
        void (TClass::*func)(const LocalVector& u,
        					 std::vector<std::vector<TData> > vvvLinDefect[],
        					 const size_t nip))
{
	if(id >= NUM_REFERENCE_OBJECTS)
		UG_THROW("Reference Object id invalid: "<<id);

	m_vLinDefectFunc[id] = boost::bind(func, obj, _1, _2, _3);
}

template <typename TData, int dim>
void
DataImport<TData,dim>::
set_fct(ReferenceObjectID id,
             void (*func)(const LocalVector& u,
            		 	  std::vector<std::vector<TData> > vvvLinDefect[],
            		 	  const size_t nip))
{
	if(id >= NUM_REFERENCE_OBJECTS)
		UG_THROW("Reference Object id invalid: "<<id);

	m_vLinDefectFunc[id] = func;
}

template <typename TData, int dim>
void DataImport<TData,dim>::clear_fct()
{
	for(size_t i = 0; i < NUM_REFERENCE_OBJECTS; ++i)
		m_vLinDefectFunc[i] = NULL;
}


template <typename TData, int dim>
void DataImport<TData,dim>::set_data(SmartPtr<CplUserData<TData, dim> > spData)
{
//	remember UserData
	m_spUserData = spData;

//	remember iexport
	this->m_spICplUserData = spData;

//	remember dependent data (i.e. is NULL iff no dependent data given)
	m_spDependentUserData = m_spUserData.template cast_dynamic<DependentUserData<TData, dim> >();
}

template <typename TData, int dim>
void DataImport<TData,dim>::cache_data_access()
{
	//	cache the pointer to the data field.
		m_vValue = m_spUserData->values(m_seriesID);

	//	in addition we cache the number of ips
		m_numIP = m_spUserData->num_ip(m_seriesID);
}

template <typename TData, int dim>
template <int ldim>
void DataImport<TData,dim>::set_local_ips(const MathVector<ldim>* vPos, size_t numIP,
                                          bool bMayChange)
{
//	if no data set, skip
	if(!data_given()) return;

//	request series if first time requested
	if(m_seriesID == -1)
	{
		m_seriesID = m_spUserData->template
					register_local_ip_series<ldim>(vPos,numIP,bMayChange);

	//	register callback, invoked when data field is changed
		m_spUserData->register_storage_callback(this, &DataImport<TData,dim>::cache_data_access);

	//	cache access to the data
		cache_data_access();

	//	resize also lin defect array
		resize_defect_array();

	//	check that num ip is correct
		UG_ASSERT(m_numIP == numIP, "Different number of ips than requested.");
	}
	else
	{
		if(!bMayChange)
			UG_THROW("DataImport: Setting different local ips to non-changable ip series.");

	//	set new local ips
		m_spUserData->template set_local_ips<ldim>(m_seriesID, vPos,numIP);

		if(numIP != m_numIP)
		{
		//	cache access to the data
			cache_data_access();

		//	resize also lin defect array
			resize_defect_array();
		}

	//	check that num ip is correct
		UG_ASSERT(m_numIP == numIP, "Different number of ips than requested.");
	}
}

template <typename TData, int dim>
void DataImport<TData,dim>::set_local_ips(const MathVector<dim>* vPos, size_t numIP,
                                          bool bMayChange)
{
	set_local_ips<dim>(vPos, numIP, bMayChange);
}

template <typename TData, int dim>
void DataImport<TData,dim>::set_global_ips(const MathVector<dim>* vPos, size_t numIP)
{
//  if no data set, skip
	if(!data_given()) return;

//	set global ips for series ID
	UG_ASSERT(m_seriesID >= 0, "Wrong series id.");
	m_spUserData->set_global_ips(m_seriesID,vPos,numIP);
}

template <typename TData, int dim>
void DataImport<TData,dim>::clear_ips()
{
	if(data_given()) m_spUserData->unregister_storage_callback(this);
	m_seriesID = -1;
	m_vValue = 0;
	m_numIP = 0;
	m_vvvLinDefect.resize(num_ip());
}

template <typename TData, int dim>
void DataImport<TData,dim>::add_jacobian(LocalMatrix& J, const number scale)
{
	UG_ASSERT(m_spDependentUserData.valid(), "No Export set.");

//	access jacobian by maps
	J.access_by_map(this->map(), this->conn_map());

//	loop integration points
	for(size_t ip = 0; ip < num_ip(); ++ip)
	{
//	loop all functions
	for(size_t fct1 = 0; fct1 < this->num_fct(); ++fct1)
		for(size_t fct2 = 0; fct2 < m_spDependentUserData->num_fct(); ++fct2)
		{
//	get array of linearized defect and derivative
	const TData* LinDef = lin_defect(ip, fct1);
	const TData* Deriv = m_spDependentUserData->deriv(m_seriesID, ip, fct2);

//	loop shapes of functions
	for(size_t sh1 = 0; sh1 < num_sh(fct1); ++sh1)
		for(size_t sh2 = 0; sh2 < m_spDependentUserData->num_sh(fct2); ++sh2)
		{
			J(fct1, sh1, fct2, sh2) += scale*(LinDef[sh1]*Deriv[sh2]);
		}
		}
	}
}

template <typename TData, int dim>
void DataImport<TData,dim>::update_dof_sizes(const LocalIndices& ind)
{
//	check size
	const FunctionIndexMapping& map = this->map();
	UG_ASSERT(map.num_fct() == this->num_fct(), "Number function mismatch.");

//	cache numFct and their numDoFs
	m_vvNumDoFPerFct.resize(map.num_fct());
	for(size_t fct = 0; fct < m_vvNumDoFPerFct.size(); ++fct)
		m_vvNumDoFPerFct[fct] = ind.num_dof(map[fct]);

	m_vvvLinDefect.clear();
	resize_defect_array();
}

template <typename TData, int dim>
void DataImport<TData,dim>::resize_defect_array()
{
//	get old size
//	NOTE: for all ips up to oldSize the arrays are already resized
	const size_t oldSize = m_vvvLinDefect.size();

//	resize ips
	m_vvvLinDefect.resize(num_ip());

//	resize num fct
	for(size_t ip = oldSize; ip < num_ip(); ++ip)
	{
	//	resize num fct
		m_vvvLinDefect[ip].resize(m_vvNumDoFPerFct.size());

	//	resize dofs
		for(size_t fct = 0; fct < m_vvNumDoFPerFct.size(); ++fct)
			m_vvvLinDefect[ip][fct].resize(m_vvNumDoFPerFct[fct]);
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

} // end namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__DATA_IMPORT_IMPL__ */
