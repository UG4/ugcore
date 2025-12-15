/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__DATA_EXPORT_IMPL__
#define __H__UG__LIB_DISC__SPATIAL_DISC__DATA_EXPORT_IMPL__

#include "data_export.h"

#include "lib_disc/reference_element/reference_element_util.h"

namespace ug {

////////////////////////////////////////////////////////////////////////////////
// DataExport
////////////////////////////////////////////////////////////////////////////////

template <typename TData, int dim>
DataExport<TData, dim>::DataExport(const char* functions)
: StdDependentUserData<DataExport<TData,dim>, TData, dim>(functions)
{
//	reset all evaluation functions
	clear_fct();
}

template <typename TData, int dim>
void DataExport<TData, dim>::clear_fct()
{
	for(size_t roid = 0; roid < NUM_REFERENCE_OBJECTS; ++roid){
		eval_fct<1>(static_cast<ReferenceObjectID>(roid)).invalidate();
		eval_fct<2>(static_cast<ReferenceObjectID>(roid)).invalidate();
		eval_fct<3>(static_cast<ReferenceObjectID>(roid)).invalidate();
	}
}

template <typename TData, int dim>
template <int refDim>
void DataExport<TData, dim>::eval_and_deriv(TData vValue[],
                    const MathVector<dim> vGlobIP[],
                    number time, int si,
                    GridObject* elem,
                    const MathVector<dim> vCornerCoords[],
                    const MathVector<refDim> vLocIP[],
                    const size_t nip,
                    LocalVector* u,
                    bool bDeriv,
                    int s,
                    std::vector<std::vector<TData> > vvvDeriv[],
                    const MathMatrix<refDim, dim>* vJT) const
{
	u->access_by_map(this->map());
	const ReferenceObjectID roid = elem->reference_object_id();
	const Functor<refDim>& func = eval_fct<refDim>(roid);
	if(func.invalid())
		UG_THROW("DataExport: no evaluation function set for "<<
		         roid<<" (world dim: "<<dim<<", ref dim: "<<refDim<<").");

	(func)(vValue,vGlobIP,time,si,*u,elem,
			vCornerCoords,vLocIP,nip, bDeriv, vvvDeriv);
}

template <typename TData, int dim>
template <typename TClass, int refDim>
void DataExport<TData, dim>::
set_fct(ReferenceObjectID id, TClass* obj,
        void (TClass::*func)(	TData vValue[],
								const MathVector<dim> vGlobIP[],
								number time, int si,
								const LocalVector& u,
								GridObject* elem,
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
						GridObject* elem,
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

#endif