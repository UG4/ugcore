/*
 * data_export_impl.h
 *
 *  Created on: 04.07.2012
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__DATA_EXPORT_IMPL__
#define __H__UG__LIB_DISC__SPATIAL_DISC__DATA_EXPORT_IMPL__

#include "data_export.h"

namespace ug{
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
		UG_THROW("Exports assume to be used by on object for all functions.");
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
void DataExport<TData, dim>::set_roid(ReferenceObjectID id)
{
	if(m_vExportFunc[id] == NULL)
		UG_THROW("DataExport::set_roid: There is no evaluation "
				"function registered for export and elem type "<<id);

	if(m_vCompFct[id] == NULL)
		UG_THROW("DataExport::set_roid: There is no evaluation forward"
				"function registered for export and elem type "<<m_id);

	m_id = id;
}

template <typename TData, int dim>
void DataExport<TData, dim>::check_setup() const
{
	if(m_id == ROID_UNKNOWN)
		UG_THROW("DataExport::check_setup: The reference element "
				"type has not been set for evaluation.");

	if(m_vCompFct[m_id] == NULL)
		UG_THROW("DataExport::check_setup: There is no evaluation forward"
				"function registered for export and elem type "<<m_id);

	if(m_vExportFunc[m_id] == NULL)
		UG_THROW("DataExport::check_setup: There is no evaluation "
				"function registered for export and elem type "<<m_id);
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
add_needed_data(SmartPtr<IUserData> data)
{
	m_vDependData.push_back(data);
}

template <typename TData, int dim>
void DataExport<TData, dim>::
remove_needed_data(SmartPtr<IUserData> data)
{
	m_vDependData.erase(remove(m_vDependData.begin(),
	                           m_vDependData.end(),
	                           data),
	                           m_vDependData.end());
}

} // end namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__DATA_EXPORT_IMPL__ */
