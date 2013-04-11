/*
 * data_export_impl.h
 *
 *  Created on: 04.07.2012
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__DATA_EXPORT_IMPL__
#define __H__UG__LIB_DISC__SPATIAL_DISC__DATA_EXPORT_IMPL__

#include "data_export.h"
#include "lib_disc/reference_element/reference_element_util.h"

namespace ug{
////////////////////////////////////////////////////////////////////////////////
// DataExport
////////////////////////////////////////////////////////////////////////////////

template <typename TData, int dim>
DataExport<TData, dim>::DataExport() : m_id(ROID_UNKNOWN)
{
//	reset all evaluation functions
	clear_fct();
}

template <typename TData, int dim>
void DataExport<TData, dim>::clear_fct()
{
	for(size_t roid = 0; roid < NUM_REFERENCE_OBJECTS; ++roid){
		eval_fct<1>((ReferenceObjectID)roid).invalidate();
		eval_fct<2>((ReferenceObjectID)roid).invalidate();
		eval_fct<3>((ReferenceObjectID)roid).invalidate();
	}
}

template <typename TData, int dim>
template <typename TClass, int refDim>
void DataExport<TData, dim>::
set_fct(ReferenceObjectID id, TClass* obj,
        void (TClass::*func)(	TData vValue[],
								const MathVector<dim> vGlobIP[],
								number time, int si,
								const LocalVector& u,
								GeometricObject* elem,
								const MathVector<dim> vCornerCoords[],
								const MathVector<refDim> vLocIP[],
								const size_t nip,
								bool bDeriv,
								std::vector<std::vector<TData> > vvvDeriv[]))
{
	if(id >= NUM_REFERENCE_OBJECTS)
		UG_THROW("Reference Object id invalid: "<<id);

	eval_fct<refDim>(id) = Functor<refDim>(obj, func);
}

template <typename TData, int dim>
template <int refDim>
void DataExport<TData, dim>::
set_fct(ReferenceObjectID id,
        void (*func)(	TData vValue[],
        				const MathVector<dim> vGlobIP[],
						number time, int si,
						const LocalVector& u,
						GeometricObject* elem,
						const MathVector<dim> vCornerCoords[],
						const MathVector<refDim> vLocIP[],
						const size_t nip,
        				bool bDeriv,
        				std::vector<std::vector<TData> > vvvDeriv[]))
{
	if(id >= NUM_REFERENCE_OBJECTS)
		UG_THROW("Reference Object id invalid: "<<id);

	eval_fct<refDim>(id) = Functor<refDim>(func);
}


template <typename TData, int dim>
template <int refDim>
inline void DataExport<TData, dim>::
comp(const LocalVector& u, GeometricObject* elem,
     const MathVector<dim> vCornerCoords[], bool bDeriv)
{
	Functor<refDim>& func = eval_fct<refDim>(m_id);

	std::vector<std::vector<TData> >* vvvDeriv = NULL;

//	evaluate for each ip series
	for(size_t s = 0; s < this->num_series(); ++s)
	{
		if(bDeriv && this->m_vvvvDeriv[s].size() > 0)
			vvvDeriv = &this->m_vvvvDeriv[s][0];
		else
			vvvDeriv = NULL;

		(func)(	this->values(s),
				this->ips(s),
				this->time(),
				this->subset(),
				u,
				elem,
				vCornerCoords,
				this->template local_ips<refDim>(s),
				this->num_ip(s),
				bDeriv,
				vvvDeriv);
	}
}

template <typename TData, int dim>
void DataExport<TData, dim>::set_roid(ReferenceObjectID id)
{
	if(id == ROID_UNKNOWN)
		UG_THROW("DataExport::set_roid: Setting unknown ReferenceObjectId.");

	m_id = id;
}

template <typename TData, int dim>
void DataExport<TData, dim>::check_setup() const
{
	if(m_id == ROID_UNKNOWN)
		UG_THROW("DataExport::check_setup: The reference element "
				"type has not been set for evaluation.");

	if(!eval_fct_set(m_id))
		UG_THROW("DataExport::check_setup: There is no evaluation "
				"function registered for data export and element type "<<m_id<<
				", but required. (world dim: "<<dim<<", ref dim: "<<
				ReferenceElementDimension(m_id)<<")");
}

template <typename TData, int dim>
bool DataExport<TData, dim>::eval_fct_set(ReferenceObjectID id) const
{
	const int d = ReferenceElementDimension(id);
	bool bRes = false;
	switch(d){
		case 1: if(eval_fct<1>(id).valid()) bRes = true; break;
		case 2: if(eval_fct<2>(id).valid()) bRes = true; break;
		case 3: if(eval_fct<3>(id).valid()) bRes = true; break;
		default: UG_THROW("DataExport: Dimension "<<d<<" not supported.");
	}
	return bRes;
}


template <typename TData, int dim>
void DataExport<TData, dim>::compute(LocalVector* u, GeometricObject* elem,
                                     const MathVector<dim> vCornerCoords[], bool bDeriv)
{
	UG_ASSERT(m_id != ROID_UNKNOWN, "Invalid RefElem");
	UG_ASSERT(m_id == elem->reference_object_id(), "Wrong element type");
	UG_ASSERT(u != NULL, "LocalVector pointer is NULL");

	switch(this->dim_local_ips()){
		case 1:	comp<1>(*u, elem, vCornerCoords, bDeriv); return;
		case 2:	comp<2>(*u, elem, vCornerCoords, bDeriv); return;
		case 3:	comp<3>(*u, elem, vCornerCoords, bDeriv); return;
		default: UG_THROW("DataExport: Dimension not supported.")
	}
}

template <typename TData, int dim>
void DataExport<TData, dim>::
add_needed_data(SmartPtr<ICplUserData<dim> > data)
{
	m_vDependData.push_back(data);
}

template <typename TData, int dim>
void DataExport<TData, dim>::
remove_needed_data(SmartPtr<ICplUserData<dim> > data)
{
	m_vDependData.erase(remove(m_vDependData.begin(),
	                           m_vDependData.end(),
	                           data),
	                           m_vDependData.end());
}

} // end namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__DATA_EXPORT_IMPL__ */
