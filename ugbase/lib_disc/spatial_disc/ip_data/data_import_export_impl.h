/*
 * data_import_export_impl.h
 *
 *  Created on: 04.07.2011
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__DATA_IMPORT_EXPORT_IMPL__
#define __H__UG__LIB_DISC__SPATIAL_DISC__DATA_IMPORT_EXPORT_IMPL__

#include "data_import_export.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
// DataImport
////////////////////////////////////////////////////////////////////////////////

template <typename TData, int dim>
DataImport<TData,dim>::~DataImport()
{
	if(data_given()) m_spIPData->unregister_storage_callback(this);
}

template <typename TData, int dim>
bool DataImport<TData,dim>::set_roid(ReferenceObjectID id)
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
		UG_LOG("ERROR in 'DataImport::set_roid':"
				"No lin defect functions registered for " << id << ".\n");
		m_id = ROID_UNKNOWN;
		return false;
	}
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
void DataImport<TData,dim>::set_data(SmartPtr<IPData<TData, dim> > spData)
{
//	remember IPData
	m_spIPData = spData;

//	remember iexport
	this->m_spIDependentIPData = spData;

//	remember dependent data (i.e. is NULL iff no dependent data given)
	m_spDependentIPData = m_spIPData.template cast_dynamic<DependentIPData<TData, dim> >();
}

template <typename TData, int dim>
void DataImport<TData,dim>::cache_data_access()
{
	//	cache the pointer to the data field.
		m_vValue = m_spIPData->values(m_seriesID);

	//	in addition we cache the number of ips
		m_numIP = m_spIPData->num_ip(m_seriesID);
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
		m_seriesID = m_spIPData->template
					register_local_ip_series<ldim>(vPos,numIP,bMayChange);

	//	register callback, invoked when data field is changed
		m_spIPData->register_storage_callback(this, &DataImport<TData,dim>::cache_data_access);

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
		m_spIPData->template set_local_ips<ldim>(m_seriesID, vPos,numIP);

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
	m_spIPData->set_global_ips(m_seriesID,vPos,numIP);
}

template <typename TData, int dim>
void DataImport<TData,dim>::clear_ips()
{
	if(data_given()) m_spIPData->unregister_storage_callback(this);
	m_seriesID = -1;
	m_vValue = 0;
	m_numIP = 0;
	m_vvvLinDefect.resize(num_ip());
}

template <typename TData, int dim>
void DataImport<TData,dim>::assemble_jacobian(LocalMatrix& J)
{
	UG_ASSERT(m_spDependentIPData.valid(), "No Export set.");

	if(m_bInRhsPart){
	//	loop integration points
		for(size_t ip = 0; ip < num_ip(); ++ip)
		{
	//	loop all functions
		for(size_t fct1 = 0; fct1 < num_fct(); ++fct1)
			for(size_t fct2 = 0; fct2 < m_spDependentIPData->num_fct(); ++fct2)
			{
	//	get array of linearized defect and derivative
		const TData* LinDef = lin_defect(ip, fct1);
		const TData* Deriv = m_spDependentIPData->deriv(m_seriesID, ip, fct2);

	//	loop shapes of functions
		for(size_t sh1 = 0; sh1 < num_sh(fct1); ++sh1)
			for(size_t sh2 = 0; sh2 < m_spDependentIPData->num_sh(fct2); ++sh2)
			{
				J(fct1, sh1, fct2, sh2) -= LinDef[sh1]*Deriv[sh2];
			}
			}
		}
	}
	else{
	//	loop integration points
		for(size_t ip = 0; ip < num_ip(); ++ip)
		{
	//	loop all functions
		for(size_t fct1 = 0; fct1 < num_fct(); ++fct1)
			for(size_t fct2 = 0; fct2 < m_spDependentIPData->num_fct(); ++fct2)
			{
	//	get array of linearized defect and derivative
		const TData* LinDef = lin_defect(ip, fct1);
		const TData* Deriv = m_spDependentIPData->deriv(m_seriesID, ip, fct2);

	//	loop shapes of functions
		for(size_t sh1 = 0; sh1 < num_sh(fct1); ++sh1)
			for(size_t sh2 = 0; sh2 < m_spDependentIPData->num_sh(fct2); ++sh2)
			{
				J(fct1, sh1, fct2, sh2) += LinDef[sh1]*Deriv[sh2];
			}
			}
		}
	}
}

template <typename TData, int dim>
void DataImport<TData,dim>::set_dof_sizes(const LocalIndices& ind,
                                          const FunctionIndexMapping& map)
{
//	check size
	UG_ASSERT(map.num_fct() == num_fct(), "Number function mismatch.");

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

////////////////////////////////////////////////////////////////////////////////
// DataExport
////////////////////////////////////////////////////////////////////////////////

template <typename TData, int dim>
DataExport<TData, dim>::DataExport() : m_id(ROID_UNKNOWN), m_pObj(NULL)
{
//	reset all evaluation functions
	clear_fct();
}

template <typename TData, int dim>
void DataExport<TData, dim>::clear_fct()
{
	for(size_t roid = 0; roid < NUM_REFERENCE_OBJECTS; ++roid)
		m_vExportFunc[roid] = NULL;
}

template <typename TData, int dim>
template <typename T, int refDim>
void DataExport<TData, dim>::
set_fct(ReferenceObjectID id, IElemDisc* obj,
        void (T::*func)(const LocalVector& u,
        				const MathVector<dim> vGlobIP[],
        				const MathVector<refDim> vLocIP[],
        				const size_t nip,
        				TData vValue[],
        				bool bDeriv,
        				std::vector<std::vector<TData> > vvvDeriv[]))
{
//	store the method pointer casted to some generic (incompatible) type
	m_vExportFunc[id] = reinterpret_cast<DummyMethod>(func);

//	store the evaluation forwarder
	m_vCompFct[id] = &DataExport<TData, dim>::template comp<T,refDim>;

//	store the base object needed for invocation
	if(m_pObj == NULL) m_pObj = obj;
	else if(m_pObj != obj)
		throw(UGError("Exports assume to be used by on object for all functions."));
}


template <typename TData, int dim>
template <typename T, int refDim>
inline void DataExport<TData, dim>::
comp(const LocalVector& u, bool bDeriv)
{
	typedef void (T::*ExpFunc)(	const LocalVector& u,
								const MathVector<dim> vGlobIP[],
								const MathVector<refDim> vLocIP[],
								const size_t nip,
								TData vValue[],
								bool bDeriv,
								std::vector<std::vector<TData> > vvvDeriv[]);


	typedef void (IElemDisc::*ElemDiscFunc)(
								const LocalVector& u,
								const MathVector<dim> vGlobIP[],
								const MathVector<refDim> vLocIP[],
								const size_t nip,
								TData vValue[],
								bool bDeriv,
								std::vector<std::vector<TData> > vvvDeriv[]);

//	cast the method pointer back to correct type
	ExpFunc func = reinterpret_cast<ExpFunc>(m_vExportFunc[m_id]);

//	cast if for evaluation using base class
	ElemDiscFunc elemDiscfunc = static_cast<ElemDiscFunc>(func);

	std::vector<std::vector<TData> >* vvvDeriv = NULL;

//	evaluate for each ip series
	for(size_t s = 0; s < this->num_series(); ++s)
	{
		if(bDeriv && this->m_vvvvDeriv[s].size() > 0)
			vvvDeriv = &this->m_vvvvDeriv[s][0];
		else
			vvvDeriv = NULL;

		(m_pObj->*(elemDiscfunc))
				(u,
				 this->ips(s),
				 this->template local_ips<refDim>(s),
				 this->num_ip(s),
				 this->values(s),
				 bDeriv,
				 vvvDeriv);
	}
}

template <typename TData, int dim>
bool DataExport<TData, dim>::set_roid(ReferenceObjectID id)
{
	if(m_vExportFunc[id] == NULL) {
		UG_LOG("ERROR in 'DataExport::set_roid': There is no evaluation "
				"function registered for export and elem type "<<id<<".\n");
		return false;
	}

	if(m_vCompFct[id] == NULL) {
		UG_LOG("ERROR in 'DataExport::set_roid': There is no evaluation forward"
				"function registered for export and elem type "<<m_id<<".\n");
		return false;
	}

	m_id = id;
	return true;
}

template <typename TData, int dim>
bool DataExport<TData, dim>::is_ready() const
{
	if(m_id == ROID_UNKNOWN) {
		UG_LOG("ERROR in 'DataExport::is_ready': The reference element "
				"type has not been set for evaluation.\n");
		return false;
	}

	if(m_vCompFct[m_id] == NULL) {
		UG_LOG("ERROR in 'DataExport::is_ready': There is no evaluation forward"
				"function registered for export and elem type "<<m_id<<".\n");
		return false;
	}

	if(m_vExportFunc[m_id] == NULL) {
		UG_LOG("ERROR in 'DataExport::is_ready': There is no evaluation "
				"function registered for export and elem type "<<m_id<<".\n");
		return false;
	}

//	everything is ok
	return true;
}

template <typename TData, int dim>
void DataExport<TData, dim>::compute(LocalVector* u, GeometricObject* elem, bool bDeriv)
{
	UG_ASSERT(m_vExportFunc[m_id] != NULL, "Func pointer is NULL");
	UG_ASSERT(m_vCompFct[m_id] != NULL, "Func pointer is NULL");
	UG_ASSERT(u != NULL, "LocalVector pointer is NULL");
	(this->*m_vCompFct[m_id])(*u, bDeriv);
}

template <typename TData, int dim>
void DataExport<TData, dim>::
add_needed_data(SmartPtr<IIPData> data)
{
	m_vDependData.push_back(data);
}

template <typename TData, int dim>
void DataExport<TData, dim>::
remove_needed_data(SmartPtr<IIPData> data)
{
	m_vDependData.erase(remove(m_vDependData.begin(),
	                           m_vDependData.end(),
	                           data),
	                           m_vDependData.end());
}


} // end namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__DATA_IMPORT_EXPORT_IMPL__ */
